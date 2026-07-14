#include "globals.hlsli"
#include "restir_giHF.hlsli"

/**
 * ReSTIR GI - temporal resampling pass.
 *
 * Merges the previous frame's reservoir (the previous spatial output,
 * reprojected through the velocity buffer and rejected on disocclusion) into
 * this frame's initial reservoir. The result feeds the spatial pass, which
 * gathers neighbors and resolves the irradiance for the denoiser.
 *
 * Reuse across frames tracks the same surface point, so the reconnection
 * Jacobian is ~1 here (the source and target visible points coincide); the full
 * Jacobian is applied by the spatial pass, whose neighbors are genuinely
 * different points. The confidence is bounded (RESTIR_MAX_HISTORY_LENGTH), with
 * weightSum scaled together with M so the unbiased weight W = weightSum /
 * (M * targetPdf) is preserved - capping M alone would inflate W every frame
 * and blow up through the history.
 */

PUSHCONSTANT(push, RESTIRGIPushConstants);

ByteAddressBuffer reservoirInitial : register(t0);
ByteAddressBuffer reservoirHistory : register(t1);

RWByteAddressBuffer reservoirOutput : register(u0);

[numthreads(8, 8, 1)]
void main(uint3 DTid : SV_DispatchThreadID)
{
	if (DTid.x >= push.resolution.x || DTid.y >= push.resolution.y)
		return;

	const uint2 pixel = DTid.xy;
	const uint flatIndex = pixel.y * push.resolution.x + pixel.x;

	RESTIRGIReservoir reservoir =
		RESTIRGIReservoirLoad(reservoirInitial, flatIndex);

	const float depth = texture_depth[pixel];
	[branch]
	if (depth > 0)
	{
		const float2 uv = (pixel + 0.5) * push.resolutionRcp;
		const float3 P = reconstruct_position(uv, depth);
		const float3 N = decode_normal(texture_normal_roughness[pixel]);

		const float2 velocity = texture_velocity[pixel];
		const float2 prevUV = uv + velocity;

		[branch]
		if (all(prevUV >= 0) && all(prevUV <= 1))
		{
			const uint2 prevPixel =
				min(uint2(prevUV * push.resolution), push.resolution - 1);

			// Disocclusion test: reject history that reprojects onto a
			// different surface (depth edge crossing) or a differently oriented
			// one (rotation past an edge). No previous normal buffer is kept,
			// so the current normal at the reprojected pixel is used as a
			// geometric proxy - the same approximation the DI temporal pass
			// makes.
			const float prevDepth = texture_depth_history.SampleLevel(
				sampler_point_clamp, prevUV, 0);
			const float linearCur = compute_lineardepth(depth);
			const float linearPrev = compute_lineardepth(prevDepth);
			const float3 prevN = decode_normal(texture_normal_roughness[prevPixel]);
			const bool consistent =
				abs(linearCur - linearPrev) <= 0.05 * linearCur &&
				dot(prevN, N) >= 0.9;

			[branch]
			if (consistent)
			{
				const uint prevFlat =
					prevPixel.y * push.resolution.x + prevPixel.x;

				RESTIRGIReservoir history =
					RESTIRGIReservoirLoad(reservoirHistory, prevFlat);

				if (history.M > (float)RESTIR_MAX_HISTORY_LENGTH)
				{
					history.weightSum *=
						(float)RESTIR_MAX_HISTORY_LENGTH / history.M;
					history.M = (float)RESTIR_MAX_HISTORY_LENGTH;
				}

				[branch]
				if (history.M > 0)
				{
					RNG rng;
					rng.init(pixel, GetFrame().frame_count + 1u);

					const float targetAtSelf = RESTIRGITargetFunction(
						history.sampleRadiance, history.samplePosition, P, N);

					RESTIRGIReservoirMerge(
						reservoir, history, targetAtSelf,
						RESTIR_GI_MAX_W, rng);
				}
			}
		}
	}

	// The spatial pass gathers neighbors into this reservoir and resolves the
	// denoiser inputs (irradiance + gradient).
	RESTIRGIReservoirStore(reservoirOutput, flatIndex, reservoir);
}
