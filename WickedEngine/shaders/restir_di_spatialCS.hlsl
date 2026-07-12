#define TEXTURE_SLOT_NONUNIFORM // final visibility ray alpha-tests materials
#include "globals.hlsli"
#include "raytracingHF.hlsli"
#include "lightingHF.hlsli"
#include "restir_lightsamplingHF.hlsli"
#include "restir_diHF.hlsli"
#include "restir_di_visibilityHF.hlsli"

/**
 * Re-resolves a reservoir sample's world position at the light's current
 * transform.
 *
 * The reservoir stores a light reference (index + uv); resolving it against the
 * live light entity each frame makes a moving light's sample track the light
 * instead of lagging behind the stale resolved position carried through reuse.
 *
 * @param[in] lightIndex - Chosen light index.
 * @param[in] uv - Stored sample parameterization on the light.
 * @param[in] P - World-space shading point.
 * @param[in] N - World-space shading normal.
 * @param[out] samplePosition - Resolved world-space point on the light.
 * @param[out] sampleRadiance - Resolved unshadowed incident radiance.
 */
void RESTIRDIResolveSample(
	uint lightIndex,
	float2 uv,
	float3 P,
	float3 N,
	out float3 samplePosition,
	out float3 sampleRadiance)
{
	const ShaderEntity light = load_entity(lightIndex);
	const RESTIRLightSample s = RESTIRResolveAnalyticLight(light, P, N, uv);
	// Directional samples have no finite position; anchor them far along the
	// sample direction (matches the initial pass).
	const float dist = (s.distance >= FLT_MAX * 0.5) ? 100000.0 : s.distance;
	samplePosition = P + s.direction * dist;
	sampleRadiance = s.radiance;
}

// ReSTIR DI - spatial resampling pass.
//
// Merges a handful of nearby reservoirs into each pixel's temporal reservoir,
// rejecting neighbors that lie on a different surface (depth or normal
// mismatch), then re-traces a fresh shadow ray for the final sample so the
// shadow reflects the current frame (the visibility cached through reuse is
// stale and would ghost). The result is the final reservoir: consumed by
// forward shading this frame and becomes the temporal history next frame.
//
// This pass also produces the A-SVGF temporal gradient for the denoiser: it
// stores this frame's raw shadow-ray result and, by re-tracing the PREVIOUS
// frame's sample under the current scene, measures exactly how much the shadow
// changed (same ray, so the difference is noiseless occluder motion).

PUSHCONSTANT(push, RESTIRDIPushConstants);

ByteAddressBuffer reservoirInput : register(t0);
ByteAddressBuffer reservoirHistory : register(t1);    // reservoir_final[prev]
Texture2D<float> rawVisibilityHistory : register(t2); // raw_visibility[prev]

RWByteAddressBuffer reservoirOutput : register(u0);
RWTexture2D<float> rawVisibilityOutput : register(u1); // raw_visibility[cur]
RWTexture2D<float> gradientOutput : register(u2);      // A-SVGF gradient

[numthreads(8, 8, 1)]
void main(uint3 DTid : SV_DispatchThreadID)
{
	if (DTid.x >= push.resolution.x || DTid.y >= push.resolution.y)
		return;

	const uint2 pixel = DTid.xy;
	const uint flatIndex = pixel.y * push.resolution.x + pixel.x;

	RESTIRDIReservoir reservoir = RESTIRDIReservoirLoad(reservoirInput, flatIndex);

	// A-SVGF temporal gradient and raw visibility for this pixel (written at
	// the end; default 0 = no change / no valid sample).
	float gradient = 0;
	float rawVisibility = 0;

	const float depth = texture_depth[pixel];
	[branch]
	if (depth > 0)
	{
		const float2 uv = (pixel + 0.5) * push.resolutionRcp;
		const float3 P = reconstruct_position(uv, depth);
		const float3 N = decode_normal(texture_normal_roughness[pixel]);
		const float linearCenter = compute_lineardepth(depth);

		RNG rng;
		rng.init(pixel, GetFrame().frame_count + 2u);

		for (uint i = 0; i < push.spatialSampleCount; ++i)
		{
			// Random neighbor within the spatial radius (uniform in a square).
			const float2 offset = (rng.next_float2() * 2 - 1) * push.spatialRadius;
			const int2 neighborPixel = int2(pixel) + int2(round(offset));

			[branch]
			if (any(neighborPixel < 0) || any(neighborPixel >= int2(push.resolution)))
				continue;

			const float neighborDepth = texture_depth[neighborPixel];
			if (neighborDepth <= 0)
				continue;

			// Reject neighbors on a different surface.
			const float linearNeighbor = compute_lineardepth(neighborDepth);
			if (abs(linearNeighbor - linearCenter) > 0.05 * linearCenter)
				continue;

			const float3 neighborN = decode_normal(texture_normal_roughness[neighborPixel]);
			if (dot(neighborN, N) < 0.9)
				continue;

			const uint neighborFlat = neighborPixel.y * push.resolution.x + neighborPixel.x;
			const RESTIRDIReservoir neighbor = RESTIRDIReservoirLoad(reservoirInput, neighborFlat);

			[branch]
			if (neighbor.M > 0)
			{
				float3 dir;
				float dist;
				const float targetAtSelf = RESTIRDITargetFunction(
					neighbor.samplePosition, neighbor.sampleRadiance, P, N, dir, dist);

				const float maxW = (float)lights().item_count();
				RESTIRDIReservoirMerge(reservoir, neighbor, targetAtSelf, maxW, rng);
			}
		}

		// Final visibility: re-trace a fresh shadow ray for the reused sample
		// so the shadow reflects the current frame instead of the stale
		// visibility carried through spatiotemporal reuse (which ghosts).
		const float W = RESTIRDIReservoirGetInvPdf(reservoir);
		[branch]
		if (W > 0 && reservoir.targetPdf > 0)
		{
			// Re-resolve the chosen sample at the light's CURRENT transform so
			// a moving light's shadow tracks it instead of lagging behind the
			// stale resolved point carried through spatiotemporal reuse.
			RESTIRDIResolveSample(
				reservoir.lightIndex, reservoir.uv, P, N,
				reservoir.samplePosition, reservoir.sampleRadiance);

			reservoir.visibility =
				RESTIRDITraceVisibility(P, N, reservoir.samplePosition, rng);
			rawVisibility = reservoir.visibility;

			// A-SVGF gradient: re-evaluate the previous frame's sample under
			// the current scene and compare to its stored raw visibility. The
			// ray is identical across the two frames, so any difference is an
			// occluder that crossed it (an exact, noiseless "the shadow changed
			// here").
			const float2 prevUV = uv + texture_velocity[pixel];
			[branch]
			if (all(prevUV >= 0) && all(prevUV <= 1))
			{
				const uint2 prevPixel =
					min(uint2(prevUV * push.resolution), push.resolution - 1);
				const float prevDepth = texture_depth_history.SampleLevel(
					sampler_point_clamp, prevUV, 0);
				const float linearPrev = compute_lineardepth(prevDepth);
				const float3 prevN =
					decode_normal(texture_normal_roughness[prevPixel]);

				// Only compare on the same surface (else it is a disocclusion the
				// temporal denoise resets anyway).
				[branch]
				if (abs(linearCenter - linearPrev) <= 0.05 * linearCenter &&
					dot(prevN, N) >= 0.9)
				{
					const uint prevFlat =
						prevPixel.y * push.resolution.x + prevPixel.x;
					const RESTIRDIReservoir prevRes =
						RESTIRDIReservoirLoad(reservoirHistory, prevFlat);
					[branch]
					if (prevRes.targetPdf > 0)
					{
						if (prevRes.lightIndex != reservoir.lightIndex)
						{
							// The chosen light switched between frames, so the
							// visibility channel now means a different light.
							// Force a reset to avoid averaging across the two
							// lights (which ghosts on multi-light shadows).
							gradient = 1;
						}
						else
						{
							// Same light: re-resolve its sample at the light's
							// current transform (so a moved light is detected
							// too) and re-trace under the current scene,
							// comparing to the stored raw visibility.
							float3 prevPos;
							float3 prevRad;
							RESTIRDIResolveSample(
								prevRes.lightIndex, prevRes.uv, P, N,
								prevPos, prevRad);
							const float vReeval =
								RESTIRDITraceVisibility(P, N, prevPos, rng);
							gradient =
								abs(vReeval - rawVisibilityHistory[prevPixel]);
						}
					}
				}
			}
		}
	}

	RESTIRDIReservoirStore(reservoirOutput, flatIndex, reservoir);
	rawVisibilityOutput[pixel] = rawVisibility;
	gradientOutput[pixel] = gradient;
}
