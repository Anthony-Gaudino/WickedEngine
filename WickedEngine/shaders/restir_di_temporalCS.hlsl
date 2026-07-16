#include "globals.hlsli"
#include "lightingHF.hlsli"
#include "restir_lightsamplingHF.hlsli"
#include "restir_diHF.hlsli"
#include "restir_di_reuseHF.hlsli"

// ReSTIR DI - temporal resampling pass.
//
// Reprojects the previous frame's final reservoir through the velocity buffer
// and merges it into this frame's initial reservoir, bounded by a confidence
// cap (RESTIR_MAX_HISTORY_LENGTH) and rejected on disocclusion (depth
// mismatch). The result feeds the spatial pass.

PUSHCONSTANT(push, RESTIRDIPushConstants);

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

	RESTIRDIReservoir reservoir = RESTIRDIReservoirLoad(reservoirInitial, flatIndex);

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
			const uint2 prevPixel = min(uint2(prevUV * push.resolution), push.resolution - 1);

			// Disocclusion test: reject history when the reprojected sample
			// lies on a different surface. Depth catches translation across
			// depth edges; the normal test catches reprojection onto a
			// differently oriented surface (e.g. rotating past an edge), which
			// reduces the smeared ghosting the depth-only test leaves behind.
			// No previous normal buffer is kept, so this compares the current
			// normal at the reprojected location - a good geometric proxy for
			// disocclusion.
			const float prevDepth = texture_depth_history.SampleLevel(sampler_point_clamp, prevUV, 0);
			const float linearCur = compute_lineardepth(depth);
			const float linearPrev = compute_lineardepth(prevDepth);
			const float3 prevN = decode_normal(texture_normal_roughness[prevPixel]);
			const bool consistent =
				abs(linearCur - linearPrev) <= 0.05 * linearCur &&
				dot(prevN, N) >= 0.9;

			[branch]
			if (consistent)
			{
				const uint prevFlat = prevPixel.y * push.resolution.x + prevPixel.x;

				RESTIRDIReservoir history = RESTIRDIReservoirLoad(reservoirHistory, prevFlat);

				// Bound the temporal confidence so stale samples cannot
				// dominate. weightSum must be scaled together with M so the
				// unbiased weight W = weightSum / (M * targetPdf) is preserved;
				// capping M alone would inflate W every frame and blow the
				// image up to white (positive feedback through the history).
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

					// Canonical (current) reservoir plus the reprojected
					// history, each tagged with its own surface so the merge is
					// the unbiased balance heuristic (same combiner as the
					// spatial pass). Post-gate the two surfaces match, so this
					// reduces to a confidence-weighted average; the balance
					// heuristic still down-weights a history whose sample no
					// longer illuminates the surface instead of letting its
					// stale M darken the pixel.
					RESTIRDIMISSource sources[RESTIR_DI_MIS_MAX_SOURCES];
					sources[0].samplePosition = reservoir.samplePosition;
					sources[0].sampleRadiance = reservoir.sampleRadiance;
					sources[0].P = P;
					sources[0].N = N;
					sources[0].M = reservoir.M;
					sources[0].weightSum = reservoir.weightSum;
					sources[0].visibility = reservoir.visibility;
					sources[0].lightIndex = reservoir.lightIndex;
					sources[0].uv = reservoir.uv;

					// The history source is tagged with the CURRENT surface (P,
					// N), not the reprojected one: the disocclusion gate above
					// already guarantees the history lies on the same surface,
					// while reconstruct_position with the current camera
					// returns a wrong world point under camera motion
					// - near a point light the position-sensitive (1/d^2)
					//   target would then corrupt the history weight and dim
					//   the pixel as the viewer moves. Same-surface tagging
					//   makes the temporal merge a correct confidence-weighted
					//   combine.
					sources[1].samplePosition = history.samplePosition;
					sources[1].sampleRadiance = history.sampleRadiance;
					sources[1].P = P;
					sources[1].N = N;
					sources[1].M = history.M;
					sources[1].weightSum = history.weightSum;
					sources[1].visibility = history.visibility;
					sources[1].lightIndex = history.lightIndex;
					sources[1].uv = history.uv;

					const float maxW = (float)lights().item_count();
					reservoir = RESTIRDIMergeBalanceHeuristic(
						sources, 2, P, N, maxW, rng);
				}
			}
		}
	}

	RESTIRDIReservoirStore(reservoirOutput, flatIndex, reservoir);
}
