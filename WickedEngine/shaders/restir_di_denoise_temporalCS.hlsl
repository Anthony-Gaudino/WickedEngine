#include "globals.hlsli"
#include "restir_diHF.hlsli"

/**
 * ReSTIR DI - visibility temporal denoise pass.
 *
 * Reprojects the previous frame's accumulated visibility and blends the
 * fresh 1-sample shadow toward a running mean (weight 1 / historyLength),
 * giving a temporally denoised shadow. Depth/normal disocclusion resets
 * the history on geometry changes.
 *
 * A shadow sweeping across static geometry is invisible to depth/normal
 * disocclusion and would smear into a trail under a long history. The
 * spatial reuse pass measures that change exactly (A-SVGF: it re-traces
 * the previous frame's sample under the current scene) and writes it to
 * the gradient texture; here that gradient is dilated and drives a
 * per-pixel history reset, so a changed shadow refreshes at once while
 * static shadows keep their full, clean history.
 *
 * The pass is read-only on the reservoir: it writes the accumulated mean
 * only into the moments texture, and the spatial pass reads it back and
 * writes the final filtered value into reservoir.visibility (so the
 * spatial blur cannot feed back into temporal history).
 */

PUSHCONSTANT(push, RESTIRDIPushConstants);

// (moment2, historyLength, unused, slowMean). The slow mean is both the
// temporal history the next frame reprojects and the denoised visibility the
// spatial pass filters; keeping it here (not in the reservoir) also stops the
// spatial blur feeding back into temporal history.
Texture2D<float4> momentsHistory : register(t0);
ByteAddressBuffer reservoirCurrent : register(t1); // reservoir_final[cur]
Texture2D<float> gradientTexture : register(t2);   // A-SVGF temporal gradient

RWTexture2D<float4> momentsOutput : register(u0);

/**
 * Maximum temporal gradient over a small neighborhood.
 *
 * The exact gradient fires only on the thin band of pixels whose occlusion
 * flipped; dilating with a max spreads the reset across the moving edge's
 * neighborhood so the whole region refreshes together.
 *
 * @param[in] center - Pixel coordinate.
 *
 * @return The maximum gradient in the (2r+1)x(2r+1) window (in [0, 1]).
 */
float DilateGradient(uint2 center)
{
	float g = 0;

	const int r = RESTIR_DI_DENOISE_GRADIENT_DILATE;
	[unroll]
	for (int dy = -r; dy <= r; ++dy)
	{
		[unroll]
		for (int dx = -r; dx <= r; ++dx)
		{
			const int2 tap = int2(center) + int2(dx, dy);
			if (tap.x < 0 || tap.y < 0 ||
				tap.x >= (int)push.resolution.x ||
				tap.y >= (int)push.resolution.y)
				continue;
			g = max(g, gradientTexture[(uint2)tap]);
		}
	}

	return g;
}

[numthreads(8, 8, 1)]
void main(uint3 DTid : SV_DispatchThreadID)
{
	if (DTid.x >= push.resolution.x || DTid.y >= push.resolution.y)
		return;

	const uint2 pixel = DTid.xy;
	const uint flatIndex = pixel.y * push.resolution.x + pixel.x;

	RESTIRDIReservoir reservoir =
		RESTIRDIReservoirLoad(reservoirCurrent, flatIndex);

	const float depth = texture_depth[pixel];
	const float W = RESTIRDIReservoirGetInvPdf(reservoir);
	[branch]
	if (depth <= 0 || !(W > 0 && reservoir.targetPdf > 0))
	{
		// Sky or no valid sample: nothing to denoise, reset history.
		momentsOutput[pixel] = 0;
		return;
	}

	const float freshVis = reservoir.visibility;

	// Exact shadow change here (dilated), measured by the spatial reuse pass.
	const float reset = saturate(
		DilateGradient(pixel) * RESTIR_DI_DENOISE_GRADIENT_SENSITIVITY);

	// Start from this frame's single sample (history length 1).
	float mean = freshVis;
	float moment2 = freshVis * freshVis;
	float historyLength = 1;

	const float2 uv = (pixel + 0.5) * push.resolutionRcp;
	const float3 N = decode_normal(texture_normal_roughness[pixel]);
	const float2 prevUV = uv + texture_velocity[pixel];

	[branch]
	if (push.frameIndex > 0 && all(prevUV >= 0) && all(prevUV <= 1))
	{
		const uint2 prevPixel =
			min(uint2(prevUV * push.resolution), push.resolution - 1);

		// Disocclusion: reset history on depth or orientation mismatch.
		const float prevDepth = texture_depth_history.SampleLevel(
			sampler_point_clamp, prevUV, 0);
		const float linearCur = compute_lineardepth(depth);
		const float linearPrev = compute_lineardepth(prevDepth);
		const float3 prevN = decode_normal(texture_normal_roughness[prevPixel]);

		[branch]
		if (abs(linearCur - linearPrev) <= 0.05 * linearCur &&
			dot(prevN, N) >= 0.9)
		{
			const float4 prevMoments = momentsHistory[prevPixel];
			const float prevMeanSlow = prevMoments.w;
			const float prevHistory = prevMoments.y;

			// The exact temporal gradient resets the history where the shadow
			// actually changed; everywhere else it accumulates up to the cap.
			historyLength = lerp(
				min(prevHistory + 1, RESTIR_DI_DENOISE_MAX_HISTORY),
				1.0, reset);

			const float alpha = 1.0 / historyLength;
			mean = lerp(prevMeanSlow, freshVis, alpha);
			moment2 = lerp(prevMoments.x, freshVis * freshVis, alpha);
		}
	}

	// The moments carry the second moment (variance for the spatial filter),
	// the history length, and the slow mean. The slow mean is both this frame's
	// denoised visibility (the spatial pass filters it) and the temporal
	// history the next frame reprojects.
	momentsOutput[pixel] = float4(moment2, historyLength, 0, mean);
}
