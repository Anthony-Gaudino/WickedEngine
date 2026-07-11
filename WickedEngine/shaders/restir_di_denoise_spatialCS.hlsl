#include "globals.hlsli"
#include "restir_diHF.hlsli"

/**
 * ReSTIR DI - spatial visibility denoise pass (edge-aware a-trous).
 *
 * Runs after the temporal denoise, which leaves the temporally accumulated
 * shadow visibility (still noisy, because a single light is resampled per
 * frame) in reservoir.visibility, and its temporal variance in the moments
 * texture.
 * This filter smooths that visibility across the surface with a 5-tap B3-spline
 * a-trous kernel, stopping at geometry edges (depth + normal) and at real
 * shadow edges (a variance-guided visibility stop): it blurs freely where the
 * temporal variance is high, which is exactly the sampling noise we want to
 * remove, while keeping crisp penumbra gradients where the variance is low.
 * Filtering in space is what lets the temporal stage stay short (low lag /
 * little ghosting).
 *
 * The caller dispatches RESTIR_DI_DENOISE_ATROUS_PASSES times with a doubling
 * tap stride. To keep the reservoir's temporal history clean, this pass never
 * writes the filtered value into the history: the first pass reads the temporal
 * mean straight from the reservoir (and variance from the moments),
 * intermediate
 * passes ping-pong two (visibility, variance) scratch textures, and the last
 * pass writes the filtered visibility back into reservoir.visibility, which is
 * the value forward shading reads.
 *
 * References:
 * Schied et al. 2017, "Spatiotemporal Variance-Guided Filtering" (SVGF).
 * https://research.nvidia.com/publication/2017-07_spatiotemporal-variance-guided-filtering-real-time-reconstruction-path-traced
 */

PUSHCONSTANT(push, RESTIRDIDenoisePushConstants);

// (moment2, historyLength, unused, slowMean); variance = moment2 - slowMean^2.
Texture2D<float4> moments : register(t0);

RWTexture2D<float2> atrousOutput : register(u0); // (visibility, variance)
RWByteAddressBuffer reservoir : register(u1);    // reservoir_final[cur]
RWTexture2D<float2> atrousInput : register(u2);  // (visibility, variance)

// B3-spline a-trous kernel, indexed by |offset|: {3/8, 1/4, 1/16}.
static const float ATROUS_KERNEL[3] = { 0.375, 0.25, 0.0625 };

/**
 * Fetches the (visibility, variance) pair for a pixel from this pass's input.
 *
 * On the first pass the input is the moments texture: the temporally denoised
 * visibility (slow mean, .w) plus its variance (second moment .x - mean^2). The
 * temporal pass leaves it there rather than in the reservoir so its
 * neighborhood prefilter stays race-free. Afterwards the input is the ping-pong
 * texture.
 *
 * @param[in] p - Pixel coordinate.
 *
 * @return float2(visibility, variance) at that pixel.
 */
float2 FetchVisVar(uint2 p)
{
	if (push.inputFromReservoir)
	{
		const float4 m = moments[p];
		const float vis = m.w;
		return float2(vis, max(0.0, m.x - vis * vis));
	}
	return atrousInput[p];
}

/**
 * Estimates the visibility variance from a 5x5 spatial neighborhood.
 *
 * With only a couple of temporally accumulated frames the temporal variance is
 * a poor estimate, and the binary per-ray visibility (one shadow ray, so 0 or
 * 1) reads as spurious edges everywhere: a variance near zero would collapse
 * the luminance edge-stop and the filter would preserve the noise instead of
 * blurring it. Estimating the variance spatially (SVGF's low-history fallback)
 * recovers a usable value so the stop widens over noisy regions while still
 * tightening on real, low-variance shadow edges.
 *
 * @param[in] center - Pixel coordinate.
 *
 * @return Estimated variance of the visibility around the pixel (>= 0).
 */
float EstimateSpatialVariance(uint2 center)
{
	float m1 = 0;
	float m2 = 0;
	float count = 0;

	[unroll]
	for (int dy = -2; dy <= 2; ++dy)
	{
		[unroll]
		for (int dx = -2; dx <= 2; ++dx)
		{
			const int2 tap = int2(center) + int2(dx, dy);
			if (tap.x < 0 || tap.y < 0 ||
				tap.x >= (int)push.resolution.x ||
				tap.y >= (int)push.resolution.y)
				continue;
			if (texture_depth[(uint2)tap] <= 0)
				continue;

			const float v = FetchVisVar((uint2)tap).x;
			m1 += v;
			m2 += v * v;
			count += 1;
		}
	}

	if (count <= 0)
		return 0;

	const float mean = m1 / count;
	return max(0.0, m2 / count - mean * mean);
}

[numthreads(8, 8, 1)]
void main(uint3 DTid : SV_DispatchThreadID)
{
	if (DTid.x >= push.resolution.x || DTid.y >= push.resolution.y)
		return;

	const uint2 pixel = DTid.xy;

	const float depth = texture_depth[pixel];
	const float2 centerVV = FetchVisVar(pixel);

	[branch]
	if (depth <= 0)
	{
		// Sky / background: nothing to filter, pass the input through.
		if (!push.outputToReservoir)
			atrousOutput[pixel] = centerVV;
		return;
	}

	const float centerZ = compute_lineardepth(depth);
	const float3 centerN = decode_normal(texture_normal_roughness[pixel]);

	// While the temporal history is short the temporal variance is unreliable,
	// so estimate it spatially (SVGF low-history fallback); otherwise use the
	// variance carried in the ping-pong texture.
	const float historyLength =
		push.inputFromReservoir ? moments[pixel].y : 0.0;
	const float centerVar = (push.inputFromReservoir && historyLength < 4.0)
		? max(centerVV.y, EstimateSpatialVariance(pixel))
		: centerVV.y;
	const float sigmaLuma =
		RESTIR_DI_DENOISE_LUMINANCE_SCALE * sqrt(centerVar) + 1e-3;

	float sumVis = 0;
	float sumVar = 0;
	float sumWeight = 0;

	const int step = (int)push.stepSize;

	[unroll]
	for (int dy = -2; dy <= 2; ++dy)
	{
		[unroll]
		for (int dx = -2; dx <= 2; ++dx)
		{
			const int2 tap = int2(pixel) + int2(dx, dy) * step;
			if (tap.x < 0 || tap.y < 0 ||
				tap.x >= (int)push.resolution.x ||
				tap.y >= (int)push.resolution.y)
				continue;

			const uint2 tapPixel = (uint2)tap;
			const float tapDepth = texture_depth[tapPixel];
			if (tapDepth <= 0)
				continue; // Sky tap.

			const float2 tapVV = FetchVisVar(tapPixel);
			const float tapZ = compute_lineardepth(tapDepth);
			const float3 tapN = decode_normal(texture_normal_roughness[tapPixel]);

			const float wKernel = ATROUS_KERNEL[abs(dx)] * ATROUS_KERNEL[abs(dy)];
			const float wNormal = pow(
				saturate(dot(centerN, tapN)), RESTIR_DI_DENOISE_NORMAL_POWER);
			const float wDepth = exp(
				-abs(centerZ - tapZ) /
				(RESTIR_DI_DENOISE_DEPTH_SCALE * centerZ + 1e-3));
			const float wLuma =
				exp(-abs(centerVV.x - tapVV.x) / sigmaLuma);

			const float w = wKernel * wNormal * wDepth * wLuma;

			sumVis += tapVV.x * w;
			// Variance transforms with the square of the weights.
			sumVar += tapVV.y * w * w;
			sumWeight += w;
		}
	}

	const float outVis = (sumWeight > 0) ? sumVis / sumWeight : centerVV.x;
	// Seed the ping-pong variance with the (spatial) center estimate on the
	// first pass; afterwards carry the filtered variance so later passes stay
	// guided.
	const float outVar = push.inputFromReservoir
		? centerVar
		: ((sumWeight > 0) ? sumVar / (sumWeight * sumWeight) : centerVV.y);

	[branch]
	if (push.outputToReservoir)
	{
		const uint index = pixel.y * push.resolution.x + pixel.x;
		RESTIRDIReservoir r = RESTIRDIReservoirLoad(reservoir, index);
		r.visibility = outVis;
		RESTIRDIReservoirStore(reservoir, index, r);
	}
	else
	{
		atrousOutput[pixel] = float2(outVis, outVar);
	}
}
