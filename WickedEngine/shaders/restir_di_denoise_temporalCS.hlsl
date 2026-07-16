#include "globals.hlsli"
#include "restir_diHF.hlsli"

/**
 * ReSTIR DI - diffuse irradiance temporal denoise pass.
 *
 * Reprojects the previous frame's accumulated diffuse irradiance and blends
 * the fresh per-frame sample toward a running mean (weight 1 / historyLength),
 * giving a temporally denoised irradiance. Depth/normal disocclusion resets
 * the history on geometry changes.
 *
 * The signal denoised here is the full diffuse contribution
 * (radiance * W * visibility * NdotL, no albedo) produced by the spatial reuse
 * pass, not the bare shadow visibility. Because standard (unbiased) reuse
 * reselects a light every frame - bright when the visible light wins, black
 * when the occluded one does - only the full contribution averages to the
 * correct color; a visibility channel cannot. Forward shading multiplies the
 * denoised irradiance by albedo and adds specular from the reservoir.
 *
 * A shadow sweeping across static geometry is invisible to depth/normal
 * disocclusion and would smear into a trail under a long history. The spatial
 * reuse pass measures that change exactly (A-SVGF: it re-traces the previous
 * frame's sample under the current scene) and writes it to the gradient
 * texture; here that gradient is dilated and drives a per-pixel history reset,
 * so a changed shadow refreshes at once while static shadows keep their full,
 * clean history.
 */

PUSHCONSTANT(push, RESTIRDIPushConstants);

// (rgb = accumulated diffuse irradiance mean, a = history length). The mean is
// both this frame's denoised irradiance (the spatial pass filters it) and the
// temporal history the next frame reprojects.
Texture2D<float4> accumHistory : register(t0);
Texture2D<float3> rawIrradiance : register(t1);  // per-frame E (spatial pass)
Texture2D<float> gradientTexture : register(t2); // A-SVGF temporal gradient
Texture2D<float> moment2History : register(t3);  // accumulated luma 2nd moment

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

	const float depth = texture_depth[pixel];
	[branch]
	if (depth <= 0)
	{
		// Sky / background: nothing to denoise, reset history.
		accumOutput[pixel] = 0;
		moment2Output[pixel] = 0;
		return;
	}

	// This frame's single-sample diffuse irradiance. A legitimately black pixel
	// (a frame that selected an occluded light) accumulates 0 - that is exactly
	// the balancing sample the mean needs, so it is never rejected.
	const float3 freshE = rawIrradiance[pixel];
	const float freshLuma = RESTIRLuma(freshE);

	// Local (3x3) statistics of this frame's raw irradiance. The neighborhood
	// mean is a low-noise estimate of the true signal here, and the deviation
	// sets the color bounding box used to reject stale history below. Averaging
	// over the neighborhood is what makes both the clamp and the change
	// detection robust to the per-pixel stochastic light selection - a single
	// reservoir sample is far too noisy to compare against history directly
	// (this is why comparing raw per-frame irradiance across frames fails).
	float3 nSum = 0;
	float3 nSum2 = 0;
	float nWeight = 0;
	[unroll]
	for (int dy = -1; dy <= 1; ++dy)
	{
		[unroll]
		for (int dx = -1; dx <= 1; ++dx)
		{
			const int2 tap = clamp(
				int2(pixel) + int2(dx, dy),
				int2(0, 0), int2(push.resolution) - 1);
			const float3 c = rawIrradiance[(uint2)tap];
			nSum += c;
			nSum2 += c * c;
			nWeight += 1;
		}
	}
	const float3 neighborMean = nSum / nWeight;
	const float3 neighborDev =
		sqrt(max((float3)0, nSum2 / nWeight - neighborMean * neighborMean));

	// Exact shadow change from the spatial reuse pass (dilated). Reliable for a
	// change to the currently held light (an occluder crossing it), but blind
	// to a light that newly appears - it only re-traces the previously chosen
	// light. The color clamp + variance-normalized change below cover that case
	// (a shadow's trailing edge, or a fill-lit region a brighter light
	// reaches).
	const float resetGradient = saturate(
		DilateGradient(pixel) * RESTIR_DI_DENOISE_GRADIENT_SENSITIVITY);

	// Start from this frame's single sample (history length 1). The second
	// luminance moment tracks the temporal variance the spatial filter needs to
	// tell residual sampling noise (high variance -> blur) from a real, static
	// shadow edge (low variance -> keep sharp) even after history is long.
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
			const float prevMoment2 = moment2History[prevPixel];

			// Spatial mean of the reprojected history luma over the same 3x3
			// neighborhood as the current-frame mean above. The change test
			// must compare like with like - spatial mean vs spatial mean.
			// Comparing the history at a single point against a spatial mean
			// falsely fires at a STATIC lighting boundary (a shadow of one
			// light filled by another), where the neighborhood mean
			// legitimately differs from the center; the spurious reset would
			// then cut short the temporal averaging that turns the stochastic
			// single-light samples back into the correct multi-light total,
			// darkening the boundary (the multi-light energy bias). Two spatial
			// means are equal there, so the test stays quiet.
			float histSum = 0;
			[unroll]
			for (int hy = -1; hy <= 1; ++hy)
			{
				[unroll]
				for (int hx = -1; hx <= 1; ++hx)
				{
					const int2 htap = clamp(
						int2(prevPixel) + int2(hx, hy),
						int2(0, 0), int2(push.resolution) - 1);
					histSum += RESTIRLuma(accumHistory[(uint2)htap].rgb);
				}
			}
			const float histNeighborLuma = histSum / 9.0;

			// Variance-normalized temporal change: compare the two low-noise
			// neighborhood means, measured in units of the history's own
			// temporal deviation so ordinary sampling noise (which the history
			// already predicts) does not trigger a reset. A shadow's trailing
			// edge - a light newly lighting this surface - makes the current
			// neighborhood far brighter than the stale history neighborhood, a
			// large change even though the A-SVGF gradient (re-tracing only the
			// previously held light) sees nothing.
			const float histLuma = RESTIRLuma(prevAccum.rgb);
			const float meanLuma = RESTIRLuma(neighborMean);
			const float temporalChange =
				abs(histNeighborLuma - meanLuma) /
				max(1e-4, histNeighborLuma + meanLuma);
			const float temporalDev =
				sqrt(max(0.0, prevMoment2 - histLuma * histLuma));
			const float relDev = temporalDev / max(1e-4, histLuma);
			const float changeSignal =
				smoothstep(0.1, 0.4, temporalChange - relDev);

			// Combined reset: the exact gradient, or the appearance / vanishing of
			// light detected from the signal itself.
			const float reset = max(resetGradient, changeSignal);

			// Color bounding-box clamp (TAA-style): reject history that lies
			// outside this frame's local color distribution. This is the
			// primary anti-trail mechanism - at a moving shadow edge the stale
			// history (dark on the trailing edge, bright on the leading edge)
			// is snapped into the current box at once, independent of how long
			// the history is. The box tightens as the reset signal rises, so a
			// confirmed change clamps hard while a clean static surface keeps a
			// loose box that preserves its denoised history.
			//
			// The box is centered on this pixel's OWN sample (freshE), not the
			// neighborhood mean: at a high-contrast shadow edge the mean is
			// pulled toward the darker side, so centering on it would clamp a
			// lit edge pixel down toward that mean and leave a dark rim
			// (visible only with
			// >= 2 lights, where the shadowed side is fill-lit and the edge has
			// real spread). The neighborhood deviation still sets the box
			// width.
			const float boxSize = lerp(2.0, 0.5, reset);
			const float3 clampedHistory = clamp(
				prevAccum.rgb,
				freshE - neighborDev * boxSize,
				freshE + neighborDev * boxSize);

			// Adaptive history length: shorten it where the signal changed so
			// the clamped history is replaced quickly, keep it long where
			// stable so the estimate stays clean. The exact gradient still
			// forces length 1 (a single confident sample) when it fires.
			const float maxHistory =
				lerp(RESTIR_DI_DENOISE_MAX_HISTORY, 4.0, changeSignal);
			historyLength =
				lerp(min(prevAccum.a + 1, maxHistory), 1.0, resetGradient);

			const float alpha = 1.0 / historyLength;
			mean = lerp(clampedHistory, freshE, alpha);
			moment2 = lerp(prevMoment2, freshLuma * freshLuma, alpha);
		}
	}

	accumOutput[pixel] = float4(mean, historyLength);
	moment2Output[pixel] = moment2;
}
