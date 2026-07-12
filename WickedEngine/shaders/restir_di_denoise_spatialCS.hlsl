#include "globals.hlsli"
#include "restir_diHF.hlsli"

/**
 * ReSTIR DI - diffuse irradiance spatial denoise pass (edge-aware a-trous).
 *
 * Runs after the temporal denoise, which leaves the temporally accumulated
 * diffuse irradiance (still noisy, because a single light is resampled per
 * frame) in the accumulation texture. This filter smooths that irradiance
 * across the surface with a 5-tap B3-spline a-trous kernel, stopping at
 * geometry edges (depth + normal) and at real shadow edges (a variance-guided
 * luminance stop): it blurs freely where the temporal variance is high, which
 * is exactly the sampling noise we want to remove, while keeping crisp penumbra
 * gradients where the variance is low. Filtering in space is what lets the
 * temporal stage stay short (low lag / little ghosting).
 *
 * The signal is the full diffuse contribution (radiance * W * visibility *
 * NdotL, no albedo); forward shading multiplies the filtered result by albedo
 * and adds specular from the reservoir.
 *
 * The caller dispatches RESTIR_DI_DENOISE_ATROUS_PASSES times with a doubling
 * tap stride. The first pass reads the temporal mean from the accumulation
 * texture, intermediate passes ping-pong two (irradiance, variance) scratch
 * textures, and the last pass writes the filtered irradiance to a dedicated
 * output texture that forward shading samples.
 *
 * References:
 * Schied et al. 2017, "Spatiotemporal Variance-Guided Filtering" (SVGF).
 * https://research.nvidia.com/publication/2017-07_spatiotemporal-variance-guided-filtering-real-time-reconstruction-path-traced
 */

PUSHCONSTANT(push, RESTIRDIDenoisePushConstants);

// (rgb = temporally accumulated irradiance mean, a = history length).
Texture2D<float4> accum : register(t0);
Texture2D<float> moment2 : register(t1); // accumulated luma 2nd moment

RWTexture2D<float4> atrousOutput : register(u0); // (rgb irradiance, variance)
RWTexture2D<float3> irradianceFinal : register(u1); // read by forward shading
RWTexture2D<float4> atrousInput : register(u2);  // (rgb irradiance, variance)

// B3-spline a-trous kernel, indexed by |offset|: {3/8, 1/4, 1/16}.
static const float ATROUS_KERNEL[3] = { 0.375, 0.25, 0.0625 };

/**
 * Perceptual luminance of a linear RGB irradiance value.
 *
 * @param[in] c - Linear RGB value.
 *
 * @return Rec. 709 luminance.
 */
float Luma(float3 c)
{
	return dot(c, float3(0.2126, 0.7152, 0.0722));
}

/**
 * Fetches the (irradiance, variance) pair for a pixel from this pass's input.
 *
 * On the first pass the input is the accumulation texture (temporally denoised
 * irradiance mean, .rgb); afterwards the input is the ping-pong texture, which
 * carries the filtered variance in .a.
 *
 * @param[in] p - Pixel coordinate.
 * @param[out] irradiance - Irradiance at that pixel.
 * @param[out] variance - Variance of the irradiance luminance at that pixel.
 */
void FetchIrrVar(uint2 p, out float3 irradiance, out float variance)
{
	if (push.inputFromReservoir)
	{
		irradiance = accum[p].rgb;
		variance = 0; // Estimated spatially by the caller on the first pass.
	}
	else
	{
		const float4 v = atrousInput[p];
		irradiance = v.rgb;
		variance = v.a;
	}
}

/**
 * Estimates the irradiance-luminance variance from a 5x5 spatial neighborhood.
 *
 * With only a couple of temporally accumulated frames the temporal variance is
 * a poor estimate and the noisy single-sample irradiance reads as spurious
 * edges everywhere: a variance near zero would collapse the luminance
 * edge-stop and the filter would preserve the noise instead of blurring it.
 * Estimating the variance spatially (SVGF's low-history fallback) recovers a
 * usable value so the stop widens over noisy regions while still tightening on
 * real, low-variance shadow edges.
 *
 * @param[in] center - Pixel coordinate.
 *
 * @return Estimated variance of the irradiance luminance around it (>= 0).
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

			float3 irr;
			float var;
			FetchIrrVar((uint2)tap, irr, var);
			const float l = Luma(irr);
			m1 += l;
			m2 += l * l;
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
	float3 centerIrr;
	float centerVarIn;
	FetchIrrVar(pixel, centerIrr, centerVarIn);

	[branch]
	if (depth <= 0)
	{
		// Sky / background: nothing to filter, pass the input through.
		if (push.outputToReservoir)
			irradianceFinal[pixel] = centerIrr;
		else
			atrousOutput[pixel] = float4(centerIrr, centerVarIn);
		return;
	}

	const float centerZ = compute_lineardepth(depth);
	const float3 centerN = decode_normal(texture_normal_roughness[pixel]);
	const float centerLuma = Luma(centerIrr);

	// On the first pass the variance is the temporal luminance variance
	// (moment2 - mean^2); this is what keeps the luminance stop open over
	// residual noise even once history is long (static), so the filter never
	// freezes the noise in place - while staying tight on a real, static shadow
	// edge, whose temporal variance is low. While the history is still short
	// the temporal estimate is unreliable, so fall back to a spatial estimate
	// (SVGF low-history fallback). Later passes carry the filtered variance.
	const float historyLength =
		push.inputFromReservoir ? accum[pixel].a : 0.0;
	float centerVar;
	if (push.inputFromReservoir)
	{
		const float temporalVar =
			max(0.0, moment2[pixel] - centerLuma * centerLuma);
		centerVar = (historyLength < 4.0)
			? max(temporalVar, EstimateSpatialVariance(pixel))
			: temporalVar;
	}
	else
	{
		centerVar = centerVarIn;
	}
	const float sigmaLuma =
		RESTIR_DI_DENOISE_LUMINANCE_SCALE * sqrt(centerVar) + 1e-3;

	float3 sumIrr = 0;
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

			float3 tapIrr;
			float tapVar;
			FetchIrrVar(tapPixel, tapIrr, tapVar);
			const float tapZ = compute_lineardepth(tapDepth);
			const float3 tapN = decode_normal(texture_normal_roughness[tapPixel]);

			const float wKernel = ATROUS_KERNEL[abs(dx)] * ATROUS_KERNEL[abs(dy)];
			const float wNormal = pow(
				saturate(dot(centerN, tapN)), RESTIR_DI_DENOISE_NORMAL_POWER);
			const float wDepth = exp(
				-abs(centerZ - tapZ) /
				(RESTIR_DI_DENOISE_DEPTH_SCALE * centerZ + 1e-3));
			const float wLuma =
				exp(-abs(centerLuma - Luma(tapIrr)) / sigmaLuma);

			const float w = wKernel * wNormal * wDepth * wLuma;

			sumIrr += tapIrr * w;
			// Variance transforms with the square of the weights.
			sumVar += tapVar * w * w;
			sumWeight += w;
		}
	}

	const float3 outIrr = (sumWeight > 0) ? sumIrr / sumWeight : centerIrr;
	// Seed the ping-pong variance with the (spatial) center estimate on the
	// first pass; afterwards carry the filtered variance so later passes stay
	// guided.
	const float outVar = push.inputFromReservoir
		? centerVar
		: ((sumWeight > 0) ? sumVar / (sumWeight * sumWeight) : centerVarIn);

	[branch]
	if (push.outputToReservoir)
		irradianceFinal[pixel] = outIrr;
	else
		atrousOutput[pixel] = float4(outIrr, outVar);
}
