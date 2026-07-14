#include "globals.hlsli"
#include "restir_diHF.hlsli"

/**
 * ReSTIR GI - indirect irradiance temporal denoise pass.
 *
 * Reprojects the previous frame's accumulated indirect irradiance and blends
 * the fresh per-frame estimate toward a running mean (weight 1 /
 * historyLength),  with a depth/normal disocclusion reset. The accumulated
 * luminance second moment feeds the spatial a-trous variance estimate.
 *
 * Unlike the shared ReSTIR DI temporal denoise, this one has no A-SVGF gradient
 * (GI has no exact per-pixel change signal) and instead clamps the reprojected
 * history to the current frame's local irradiance color box before blending
 * (see RESTIR_GI_HISTORY_CLAMP_K). That anti-ghost clamp is what stops a stale
 * value - a dark trail a camera move drags across an area, or a lingering
 * bright one - from ghosting, since the history is pulled back into the range
 * the current frame actually observes.
 *
 * References: EmbarkStudios/kajiya - rtdgi/temporal_filter.hlsl (color bbox).
 */

PUSHCONSTANT(push, RESTIRDIPushConstants);

// (rgb = accumulated indirect irradiance mean, a = history length).
Texture2D<float4> accumHistory : register(t0);
Texture2D<float3> rawIrradiance : register(t1);  // per-frame E (spatial resolve)
Texture2D<float> moment2History : register(t2);  // accumulated luma 2nd moment

RWTexture2D<float4> accumOutput : register(u0);
RWTexture2D<float> moment2Output : register(u1); // luma 2nd moment for variance

/**
 * Rec. 709 luminance of a linear RGB value.
 *
 * @param[in] c - Linear RGB value.
 *
 * @return Scalar luminance.
 */
float RESTIRLuma(float3 c)
{
	return dot(c, float3(0.2126, 0.7152, 0.0722));
}

[numthreads(8, 8, 1)]
void main(uint3 DTid : SV_DispatchThreadID)
{
	if (DTid.x >= push.resolution.x || DTid.y >= push.resolution.y)
		return;

	const uint2 pixel = DTid.xy;

	const float depth = texture_depth[pixel];
	[branch]
	if (depth <= 0)
	{
		// Sky / background: nothing to denoise, reset history.
		accumOutput[pixel] = 0;
		moment2Output[pixel] = 0;
		return;
	}

	const float3 freshE = rawIrradiance[pixel];
	const float freshLuma = RESTIRLuma(freshE);

	// Local irradiance color box (mean +/- k*stddev over a 3x3 neighborhood of
	// the resolved irradiance) used to clamp the reprojected history below.
	float3 nSum = 0;
	float3 nSum2 = 0;
	float nCount = 0;
	[unroll]
	for (int dy = -1; dy <= 1; ++dy)
	{
		[unroll]
		for (int dx = -1; dx <= 1; ++dx)
		{
			const int2 tap = int2(pixel) + int2(dx, dy);
			if (tap.x < 0 || tap.y < 0 ||
				tap.x >= (int)push.resolution.x ||
				tap.y >= (int)push.resolution.y)
				continue;
			if (texture_depth[(uint2)tap] <= 0)
				continue;
			const float3 e = rawIrradiance[(uint2)tap];
			nSum += e;
			nSum2 += e * e;
			nCount += 1;
		}
	}
	const float3 nMean = nSum / max(nCount, 1);
	const float3 nDev =
		sqrt(max(0.0, nSum2 / max(nCount, 1) - nMean * nMean));

	// Start from this frame's estimate (history length 1). The second luminance
	// moment tracks the temporal variance the spatial filter needs to tell
	// residual noise (high variance -> blur) from real detail (low variance ->
	// keep sharp).
	float3 mean = freshE;
	float moment2 = freshLuma * freshLuma;
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
			const float4 prevAccum = accumHistory[prevPixel];

			// Anti-ghost: clamp the reprojected history into the current local
			// color box so a stale value cannot trail across motion.
			const float3 clampedHistory = clamp(prevAccum.rgb,
				nMean - nDev * RESTIR_DENOISE_HISTORY_CLAMP_K,
				nMean + nDev * RESTIR_DENOISE_HISTORY_CLAMP_K);

			// How far the history lay outside the color box (a soft,
			// color-based disocclusion). A clearly stale value resets the
			// history length toward 1 so it fades in about a frame instead of
			// over the whole history; a value already inside the box keeps its
			// full confidence.
			const float clampDelta = RESTIRLuma(abs(clampedHistory - prevAccum.rgb));
			const float ref =
				RESTIRLuma(nMean) + RESTIRLuma(prevAccum.rgb) + 1e-3;
			const float staleness =
				saturate(clampDelta / ref * RESTIR_DENOISE_STALE_SENSITIVITY);

			historyLength = lerp(
				min(prevAccum.a + 1, RESTIR_DI_DENOISE_MAX_HISTORY),
				1.0, staleness);

			const float alpha = 1.0 / historyLength;
			mean = lerp(clampedHistory, freshE, alpha);
			moment2 = lerp(
				moment2History[prevPixel], freshLuma * freshLuma, alpha);
		}
	}

	accumOutput[pixel] = float4(mean, historyLength);
	moment2Output[pixel] = moment2;
}
