#define TEXTURE_SLOT_NONUNIFORM // final visibility ray alpha-tests materials
#include "globals.hlsli"
#include "raytracingHF.hlsli"
#include "lightingHF.hlsli"
#include "restir_lightsamplingHF.hlsli"
#include "restir_diHF.hlsli"
#include "restir_di_reuseHF.hlsli"
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
ByteAddressBuffer reservoirHistory : register(t1);      // reservoir_final[prev]
Texture2D<float3> rawIrradianceHistory : register(t2);  // raw_irradiance[prev]

RWByteAddressBuffer reservoirOutput : register(u0);
RWTexture2D<float> gradientOutput : register(u1);       // A-SVGF gradient
RWTexture2D<float3> irradianceOutput : register(u2);    // per-frame diffuse E

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
	const uint flatIndex = pixel.y * push.resolution.x + pixel.x;

	RESTIRDIReservoir reservoir = RESTIRDIReservoirLoad(reservoirInput, flatIndex);

	// A-SVGF temporal gradient and per-frame diffuse irradiance for this pixel
	// (written at the end; default 0 = no change / no valid sample).
	float gradient = 0;
	float3 irradiance = 0;

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

		// Source 0 is this pixel's own (canonical) reservoir; spatial neighbors
		// are gathered after it. The balance-heuristic merge tags each source
		// with the surface it was built on so a neighbor that cannot see the
		// chosen light no longer dilutes the estimate.
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
		uint sourceCount = 1;

		for (uint i = 0; i < push.spatialSampleCount &&
			sourceCount < RESTIR_DI_MIS_MAX_SOURCES; ++i)
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
			if (neighbor.M <= 0)
				continue;

			// Reconstruct the neighbor's own surface: the MIS denominator
			// evaluates each sample's target there, which is what removes the
			// biased-combiner darkening.
			const float2 neighborUV = (neighborPixel + 0.5) * push.resolutionRcp;
			const float3 neighborP = reconstruct_position(neighborUV, neighborDepth);

			sources[sourceCount].samplePosition = neighbor.samplePosition;
			sources[sourceCount].sampleRadiance = neighbor.sampleRadiance;
			sources[sourceCount].P = neighborP;
			sources[sourceCount].N = neighborN;
			sources[sourceCount].M = neighbor.M;
			sources[sourceCount].weightSum = neighbor.weightSum;
			sources[sourceCount].visibility = neighbor.visibility;
			sources[sourceCount].lightIndex = neighbor.lightIndex;
			sources[sourceCount].uv = neighbor.uv;
			++sourceCount;
		}

		const float maxW = (float)lights().item_count();
		reservoir = RESTIRDIMergeBalanceHeuristic(
			sources, sourceCount, P, N, maxW, rng);

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

			// Per-frame diffuse irradiance for the full-signal denoiser: the
			// chosen sample's contribution without albedo (radiance * W *
			// visibility * NdotL). Denoising this - rather than the bare
			// visibility - lets the stochastic per-frame light selection
			// (bright when the visible light wins, black when the occluded one
			// does) average to the correct colour, which the visibility channel
			// alone cannot represent. W is firefly-clamped to the light count
			// exactly as forward shading does.
			const float Wclamp = min(W, (float)lights().item_count());
			const float3 L = normalize(reservoir.samplePosition - P);
			const float NdotL = saturate(dot(N, L));
			float3 irr =
				reservoir.sampleRadiance * Wclamp * reservoir.visibility * NdotL;
			irradiance =
				(any(isnan(irr)) || any(isinf(irr))) ? (float3)0 : irr;

			// A-SVGF gradient: re-evaluate the previous frame's sample's full
			// diffuse contribution under the current scene and compare to its
			// stored value. Because the denoised signal is now irradiance (not
			// bare visibility), the change detector must compare irradiance too
			// - that is what makes a moving light update at once: its shadow
			//   may not change (visibility 1 -> 1) yet its radiance / NdotL /
			//   distance do, which a visibility-only gradient could never see
			//   (the lit region would fade over the whole history instead).
			//   Skip the whole gradient (and its shadow ray) when nothing that
			//   can change a shadow moved this frame: the gradient is then
			//   identically zero, so leaving it at its default 0 gives the same
			//   result while saving one ray/pixel. Camera-only motion keeps
			//   dynamicLighting at 0 and is handled by the denoiser's geometric
			//   disocclusion instead.
			const float2 prevUV = uv + texture_velocity[pixel];
			[branch]
			if (push.dynamicLighting != 0 && all(prevUV >= 0) && all(prevUV <= 1))
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
						// Re-resolve the previous frame's chosen sample at the
						// light's current transform (so a moved light is
						// detected too), re-trace it under the current scene,
						// and rebuild its full diffuse contribution. The sample
						// is identical across the two frames, so any difference
						// is exact (noiseless): an occluder that crossed the
						// ray, or the light itself moving. We do not force a
						// reset when the chosen light switches - standard
						// (unbiased) reuse reselects a light every frame, and
						// the denoiser accumulates the contribution (meaningful
						// across switches), so a switch is not a change to
						// react to.
						float3 prevPos;
						float3 prevRad;
						RESTIRDIResolveSample(
							prevRes.lightIndex, prevRes.uv, P, N,
							prevPos, prevRad);
						const float vReeval =
							RESTIRDITraceVisibility(P, N, prevPos, rng);
						const float NdotLprev =
							saturate(dot(N, normalize(prevPos - P)));
						const float Wprev = min(
							RESTIRDIReservoirGetInvPdf(prevRes),
							(float)lights().item_count());
						const float irrNow =
							RESTIRLuma(prevRad * Wprev * vReeval * NdotLprev);
						const float irrThen =
							RESTIRLuma(rawIrradianceHistory[prevPixel]);
						// Relative change in [0, 1] so the reset is scale-free
						// across dim and bright lighting alike.
						gradient = abs(irrNow - irrThen) /
							(max(irrNow, irrThen) + 1e-3);
					}
				}
			}
		}
	}

	RESTIRDIReservoirStore(reservoirOutput, flatIndex, reservoir);
	gradientOutput[pixel] = gradient;
	irradianceOutput[pixel] = irradiance;
}
