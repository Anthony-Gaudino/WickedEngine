#include "globals.hlsli"
#include "restir_diHF.hlsli"

/**
 * ReSTIR DI - visibility temporal denoise pass.
 *
 * Runs after the spatial pass, which leaves a fresh but noisy 1-sample shadow
 * visibility in reservoir_final[cur].visibility. This pass reprojects the
 * previous frame's accumulated visibility and blends the fresh value toward a
 * running mean (weight 1 / historyLength, capped so it stays responsive),
 * giving a temporally denoised shadow. Disocclusion (depth + normal) resets the
 * history.
 *
 * Moving shadows are the hard case: when a shadow sweeps across a surface the
 * geometry is unchanged, so depth/normal disocclusion cannot catch it and a
 * long history smears the shadow into a trail. To keep a long (clean) history
 * everywhere it is safe while staying ghost-free where it is not, the pass runs
 * a per-pixel temporal-gradient antilag: it prefilters the current visibility
 * over a small neighborhood (the raw per-pixel value is a single binary shadow
 * ray) and compares that estimate against the reprojected history. When the
 * gradient exceeds what the estimation noise explains the shadow truly changed,
 * so the effective history is collapsed toward the fresh sample (the spatial
 * a-trous then cleans that low-history region using its spatial variance
 * fallback); otherwise the full history accumulates and the shadow stays clean.
 *
 * This pass is read-only on the reservoir (it reads the fresh visibility and a
 * small neighborhood of it): the accumulated mean is written only into the
 * moments texture, and the spatial pass reads it from there and writes the
 * final filtered value into reservoir.visibility. Writing the mean back into
 * the reservoir here would create a read-after-write race with the neighborhood
 * prefilter (neighbor threads overwriting the fresh value mid-dispatch).
 */

PUSHCONSTANT(push, RESTIRDIPushConstants);

// (moment2, historyLength, unused, slowMean). The slow mean is both the
// temporal history the next frame reprojects and the denoised visibility the
// spatial pass filters; keeping it here (not in the reservoir) also stops the
// spatial blur feeding back into temporal history.
Texture2D<float4> momentsHistory : register(t0);
ByteAddressBuffer reservoirCurrent : register(t1); // reservoir_final[cur], read-only

RWTexture2D<float4> momentsOutput : register(u0);

/**
 * Spatially prefiltered estimate of the current visibility for change
 * detection.
 *
 * The raw per-pixel visibility is a single binary shadow ray, far too noisy to
 * compare against the smooth reprojected history directly. Averaging a small
 * neighborhood yields a lower-noise estimate of the true mean plus the spatial
 * variance, which together bound how much of a history/current gradient is just
 * sampling noise.
 *
 * @param[in] center - Pixel coordinate.
 * @param[out] mean - Neighborhood mean visibility.
 * @param[out] variance - Neighborhood variance of the visibility.
 * @param[out] count - Number of valid (non-sky) taps averaged.
 */
void PrefilterCurrent(
	uint2 center, out float mean, out float variance, out float count)
{
	float m1 = 0;
	float m2 = 0;
	count = 0;

	const int r = RESTIR_DI_DENOISE_PREFILTER_RADIUS;
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
			if (texture_depth[(uint2)tap] <= 0)
				continue;

			const uint index = tap.y * push.resolution.x + tap.x;
			const float v = RESTIRDIReservoirLoad(reservoirCurrent, index).visibility;
			m1 += v;
			m2 += v * v;
			count += 1;
		}
	}

	mean = (count > 0) ? m1 / count : 0;
	variance = (count > 0) ? max(0.0, m2 / count - mean * mean) : 0;
}

[numthreads(8, 8, 1)]
void main(uint3 DTid : SV_DispatchThreadID)
{
	if (DTid.x >= push.resolution.x || DTid.y >= push.resolution.y)
		return;

	const uint2 pixel = DTid.xy;
	const uint flatIndex = pixel.y * push.resolution.x + pixel.x;

	RESTIRDIReservoir reservoir = RESTIRDIReservoirLoad(reservoirCurrent, flatIndex);

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

	// Spatially prefiltered current estimate + its noise, for change detection.
	float curMean;
	float curVariance;
	float curCount;
	PrefilterCurrent(pixel, curMean, curVariance, curCount);

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
		const uint2 prevPixel = min(uint2(prevUV * push.resolution), push.resolution - 1);

		// Disocclusion: reject history on depth or orientation mismatch.
		const float prevDepth = texture_depth_history.SampleLevel(sampler_point_clamp, prevUV, 0);
		const float linearCur = compute_lineardepth(depth);
		const float linearPrev = compute_lineardepth(prevDepth);
		const float3 prevN = decode_normal(texture_normal_roughness[prevPixel]);

		[branch]
		if (abs(linearCur - linearPrev) <= 0.05 * linearCur && dot(prevN, N) >= 0.9)
		{
			const float4 prevMoments = momentsHistory[prevPixel];

			const float prevMeanSlow = prevMoments.w;
			const float prevHistory = prevMoments.y;

			// Temporal-gradient antilag. Compare the prefiltered current
			// estimate against the reprojected history; the standard error of
			// each mean bounds how much of the gap is just sampling noise. When
			// the gap exceeds that band the shadow truly moved, so collapse the
			// history toward the fresh sample; otherwise keep accumulating up
			// to the cap.
			const float prevVar =
				max(0.0, prevMoments.x - prevMeanSlow * prevMeanSlow);
			const float stderrCur = sqrt(curVariance / max(curCount, 1.0));
			const float stderrHist = sqrt(prevVar / max(prevHistory, 1.0));
			const float noise =
				sqrt(stderrCur * stderrCur + stderrHist * stderrHist) + 1e-3;

			const float gradient = abs(curMean - prevMeanSlow);
			const float rejection = saturate(
				(gradient - noise) / (RESTIR_DI_DENOISE_ANTILAG_SCALE * noise));

			historyLength =
				lerp(min(prevHistory + 1, RESTIR_DI_DENOISE_MAX_HISTORY),
					1.0, rejection);

			// Clamp the reprojected history to the current estimate's
			// confidence interval before blending, so the residual stale value
			// the weight rejection leaves behind at a moving edge is pulled
			// into range.
			const float clampHalf = RESTIR_DI_DENOISE_CLAMP_SCALE * noise;
			const float prevMeanClamped =
				clamp(prevMeanSlow, curMean - clampHalf, curMean + clampHalf);

			const float alpha = 1.0 / historyLength;
			mean = lerp(prevMeanClamped, freshVis, alpha);
			moment2 = lerp(prevMoments.x, freshVis * freshVis, alpha);
		}
	}

	// The moments carry the second moment (variance for the spatial filter),
	// the history length, and the slow mean. The slow mean is both this frame's
	// denoised visibility (the spatial pass filters it) and the temporal
	// history the next frame reprojects.
	momentsOutput[pixel] = float4(moment2, historyLength, 0, mean);
}
