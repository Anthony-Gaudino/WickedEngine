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
 * long history smears the shadow into a trail. A single binary visibility
 * sample cannot tell "the shadow moved" (mean shifted) from "a different light
 * was sampled this frame" (RIS noise). To separate them the pass keeps a
 * second, short-window running mean (the fast accumulator): a *persistent* gap
 * between the fast and slow means means the true mean shifted, whereas sampling
 * noise leaves them in agreement. That gap drives an antilag term that
 * collapses the slow history so the trail clears immediately, while static
 * penumbrae keep their long, well-denoised history.
 *
 * The accumulated mean is written back in place into the reservoir's visibility
 * field so forward shading needs no change; (second moment, history length,
 * fast-accumulator mean) are kept in a ping-pong texture for antilag and for
 * the variance the spatial filter uses.
 */

PUSHCONSTANT(push, RESTIRDIPushConstants);

// (moment2, historyLength, meanFast, slowMean). The slow mean is the temporal
// history the next frame reprojects; keeping it here (not in the reservoir)
// means the spatial pass, which overwrites reservoir.visibility with the
// spatially filtered result, cannot feed that blur back into temporal history.
Texture2D<float4> momentsHistory : register(t0);

RWByteAddressBuffer reservoirCurrent : register(u0); // reservoir_final[cur], in place
RWTexture2D<float4> momentsOutput : register(u1);

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

	// Start from this frame's single sample (history length 1).
	float mean = freshVis;
	float moment2 = freshVis * freshVis;
	float meanFast = freshVis;
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
			const float prevMeanFast = prevMoments.z;

			// Fast accumulator: short window, tracks change quickly.
			const float fastLength =
				min(prevMoments.y + 1, RESTIR_DI_DENOISE_FAST_HISTORY);
			meanFast = lerp(prevMeanFast, freshVis, 1.0 / fastLength);

			// Antilag: a persistent gap between the fast and slow means is real
			// motion (the shadow moved); sampling noise leaves them close. The
			// gap is already denoised over the two windows, so it barely reacts
			// to single-sample flips. Amplify it so a moderate gap drives a
			// strong response, then collapse the slow history all the way to
			// the fresh sample at full antilag: a moving shadow snaps clear
			// (accepting some transient noise while in motion) and re-denoises
			// the instant it stops, while a small penumbra gap keeps most of
			// the long, well-denoised history.
			const float antilag = saturate(abs(meanFast - prevMeanSlow) * 2.0);
			historyLength =
				lerp(min(prevMoments.y + 1, RESTIR_DI_DENOISE_MAX_HISTORY),
					1.0, antilag);

			const float alpha = 1.0 / historyLength;
			mean = lerp(prevMeanSlow, freshVis, alpha);
			moment2 = lerp(prevMoments.x, freshVis * freshVis, alpha);
		}
	}

	// Write the temporal mean into the reservoir's visibility (the spatial pass
	// reads it there, then overwrites it with the filtered result shading
	// uses). The moments carry the second moment (variance for the spatial
	// filter), the history length, the fast-accumulator mean (antilag), and the
	// slow mean (the temporal history the next frame reprojects).
	reservoir.visibility = mean;
	RESTIRDIReservoirStore(reservoirCurrent, flatIndex, reservoir);
	momentsOutput[pixel] = float4(moment2, historyLength, meanFast, mean);
}
