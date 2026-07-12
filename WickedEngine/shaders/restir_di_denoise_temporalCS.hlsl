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
 * correct colour; a visibility channel cannot. Forward shading multiplies the
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

	// Exact shadow change here (dilated), measured by the spatial reuse pass.
	const float reset = saturate(
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

			// The exact temporal gradient resets the history where the shadow
			// actually changed; everywhere else it accumulates up to the cap.
			historyLength = lerp(
				min(prevAccum.a + 1, RESTIR_DI_DENOISE_MAX_HISTORY),
				1.0, reset);

			const float alpha = 1.0 / historyLength;
			mean = lerp(prevAccum.rgb, freshE, alpha);
			moment2 = lerp(
				moment2History[prevPixel], freshLuma * freshLuma, alpha);
		}
	}

	accumOutput[pixel] = float4(mean, historyLength);
	moment2Output[pixel] = moment2;
}
