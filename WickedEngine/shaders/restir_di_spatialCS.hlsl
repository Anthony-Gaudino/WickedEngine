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
	// lightIndex is a STABLE scene index; map it to this frame's entity slot. A
	// light culled this frame reaches no visible surface, so yield a zero-
	// radiance sample (the callers weight by radiance, so it contributes
	// nothing).
	const uint slot = RESTIRDICurrentLightSlot(lightIndex);
	[branch]
	if (slot == RESTIR_INVALID_LIGHT_INDEX)
	{
		samplePosition = P;
		sampleRadiance = 0;
		return;
	}

	const ShaderEntity light = load_entity(slot);
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

#ifdef RDI_HASHGRID
#include "restir_di_hashgridHF.hlsli"

// Read-only views of the hash grid built this frame (see the build passes).
StructuredBuffer<uint> hg_checksums : register(t3);
StructuredBuffer<uint> hg_cellConfidences : register(t4);
StructuredBuffer<uint> hg_cellCounters : register(t5); // [cell] res, [+N] inv
StructuredBuffer<uint> hg_cellOffsets : register(t6);  // [cell] res, [+N] inv
StructuredBuffer<RESTIRDIGridReservoir> hg_resSorted : register(t7);
StructuredBuffer<uint2> hg_invSorted : register(t8);   // (flatPixel, asuint(M))

/** Division guarded against a non-positive denominator. */
float RESTIRHGSafeDiv(float a, float b)
{
	return abs(b) > 0 ? a / b : 0;
}

/** Confidence quantize/dequantize (must match the build pass's atomic sums). */
float RESTIRHGRequant(float c)
{
	return (float)((uint)(c * 1000.0)) / 1000.0;
}
float RESTIRHGDequant(uint q)
{
	return (float)q / 1000.0;
}

/**
 * Defensive pairwise-MIS weight of a candidate's own contribution.
 *
 * @param[in] ci - Candidate confidence (cell-scaled).
 * @param[in] cc - Canonical confidence.
 * @param[in] cSum - Cell confidence (scaled to the candidate budget).
 * @param[in] pOwn - Candidate sample's target at its own surface.
 * @param[in] pCenter - Candidate sample's target at the center surface.
 *
 * @return The MIS weight.
 */
float RESTIRHGNoncanonical(
	float ci, float cc, float cSum, float pOwn, float pCenter)
{
	return (cSum / (cSum + cc)) *
		RESTIRHGSafeDiv(ci * pOwn, cSum * pOwn + cc * pCenter);
}

/**
 * One candidate's contribution to the canonical's defensive pairwise-MIS
 * weight.
 *
 * @param[in] ci - Candidate confidence (cell-scaled).
 * @param[in] cc - Canonical confidence.
 * @param[in] cSum - Cell confidence (scaled to the candidate budget).
 * @param[in] pOwn - Canonical sample's target at the candidate's surface.
 * @param[in] pCenter - Canonical sample's target at the center surface.
 *
 * @return The additive contribution to the canonical MIS weight.
 */
float RESTIRHGCanonical(
	float ci, float cc, float cSum, float pOwn, float pCenter)
{
	return (ci / (cSum + cc)) *
		RESTIRHGSafeDiv(cc * pCenter, cc * pCenter + cSum * pOwn);
}

/**
 * Importance-samples one reservoir from a cell (RIS over 8 stratified sub-taps,
 * weighted by the stored resampling-weight sum = the cell importance).
 *
 * @param[in] rangeBase - First sorted index of the cell's reservoirs.
 * @param[in] rangeN - Number of reservoirs in the cell.
 * @param[in,out] rng - Random generator.
 * @param[out] selectionW - Unbiasing weight of the returned candidate (0 if
 *   none had positive importance).
 *
 * @return The chosen reservoir's offset within the cell.
 */
uint RESTIRHGSampleCell(
	uint rangeBase, uint rangeN, inout RNG rng, out float selectionW)
{
	float wsum = 0;
	uint selected = 0;
	float selectedp = 0;
	[loop]
	for (uint k = 0; k < 8; ++k)
	{
		const float u = (k + rng.next_float()) / 8.0;
		const uint j = min((uint)(u * rangeN), rangeN - 1);
		const float p = hg_resSorted[rangeBase + j].weightSum;
		const float w = p * rangeN;
		wsum += w;
		if (wsum * rng.next_float() < w)
		{
			selected = j;
			selectedp = p;
		}
	}
	selectionW = selectedp > 0 ? wsum / (8.0 * selectedp) : 0;
	return selected;
}

/**
 * Stochastic Pairwise MIS spatial reuse from the hash grid.
 *
 * Finds this pixel's surface cell, importance-samples
 * RESTIR_HASHGRID_CANDIDATES reservoirs from it, and combines them with the
 * canonical via defensive pairwise MIS - scaling every confidence by the cell
 * budget and using a single inverse-map sample for the canonical term, exactly
 * as the reference. Only the chosen light reference is written; the caller
 * re-resolves the sample and re-traces its shadow.
 *
 * References:
 * Hedstrom et al. 2026, "Stochastic Pairwise MIS for Unbiased Large-Kernel
 * Reuse in Real Time".
 *
 * @param[in] pixel - Screen pixel.
 * @param[in] P - Center world-space shading point.
 * @param[in] N - Center world-space shading normal.
 * @param[in] canonical - This pixel's own reservoir.
 * @param[in,out] rng - Random generator.
 *
 * @return The merged reservoir.
 */
RESTIRDIReservoir RESTIRDIMergeHashGrid(
	uint2 pixel, float3 P, float3 N, RESTIRDIReservoir canonical, inout RNG rng)
{
	RESTIRDIReservoir result = canonical;

	const uint cellCount = push.resolution.x * push.resolution.y;
	const float linearDepth = compute_lineardepth(texture_depth[pixel]);
	const RESTIRHashKey key = RESTIRHashGridKey(pixel, N, linearDepth);
	const uint cell = RESTIRHashSetFind(hg_checksums, cellCount, key);
	[branch]
	if (cell == RESTIR_HASHGRID_INVALID)
		return result;

	const uint invBase = hg_cellOffsets[cellCount + cell];
	const uint invN = hg_cellCounters[cellCount + cell];
	const uint resBase = hg_cellOffsets[cell];
	const uint resN = hg_cellCounters[cell];
	[branch]
	if (invN == 0)
		return result;

	const float invNf = (float)invN;
	const float candScale = (float)RESTIR_HASHGRID_CANDIDATES / invNf;
	const float cSum =
		RESTIRHGDequant(hg_cellConfidences[cell]) * candScale;

	// Canonical (this pixel's own reservoir).
	bool cValid = canonical.M > 0 &&
		canonical.lightIndex != RESTIR_INVALID_LIGHT_INDEX &&
		(canonical.lightIndex & RESTIR_LIGHT_FLAG_EMISSIVE_TRIANGLE) == 0;
	ShaderEntity cLight = (ShaderEntity)0;
	float cW = 0;
	float cPCenter = 0;
	float cc = 0;
	[branch]
	if (cValid)
	{
		const uint slot = RESTIRDICurrentLightSlot(canonical.lightIndex);
		[branch]
		if (slot == RESTIR_INVALID_LIGHT_INDEX)
		{
			cValid = false;
		}
		else
		{
			cLight = load_entity(slot);
			cPCenter = RESTIRDIReuseTarget(cLight, canonical.uv, P, N);
			[branch]
			if (cPCenter <= 0)
			{
				cValid = false;
			}
			else
			{
				const float maxW =
					GetFrame().totalLightPower / RESTIRLightPower(cLight);
				cW = min(
					canonical.weightSum / (canonical.M * cPCenter), maxW);
				cc = RESTIRHGRequant(canonical.M);
			}
		}
	}
	const float centerConfidence = cc;

	float wSum = 0;

	// Canonical MIS weight, estimated from one inverse-map (all-pixels) sample.
	[branch]
	if (cValid)
	{
		const uint j = min((uint)(rng.next_float() * invNf), invN - 1);
		const uint2 pc = hg_invSorted[invBase + j];
		const uint2 pj =
			uint2(pc.x % push.resolution.x, pc.x / push.resolution.x);
		const float ci = RESTIRHGRequant(asfloat(pc.y)) * candScale;

		const float pjDepth = texture_depth[pj];
		float shiftedTargetPdf = 0;
		[branch]
		if (pjDepth > 0)
		{
			const float3 Pj =
				reconstruct_position((pj + 0.5) * push.resolutionRcp, pjDepth);
			const float3 Nj = decode_normal(texture_normal_roughness[pj]);
			shiftedTargetPdf =
				RESTIRDIReuseTarget(cLight, canonical.uv, Pj, Nj);
		}

		const float selectionW = invNf;
		float mc = centerConfidence / (centerConfidence + cSum);
		mc += selectionW * RESTIRHGCanonical(
			ci, centerConfidence, cSum, shiftedTargetPdf, cPCenter);

		wSum = mc * cPCenter * cW;
		result.targetPdf = cPCenter;
		result.lightIndex = canonical.lightIndex;
		result.uv = canonical.uv;
	}

	// Non-canonical: importance-sampled cell reservoirs.
	[branch]
	if (resN > 0)
	{
		[loop]
		for (uint i = 0; i < RESTIR_HASHGRID_CANDIDATES; ++i)
		{
			float selectionW;
			const uint j = RESTIRHGSampleCell(resBase, resN, rng, selectionW);
			[branch]
			if (!(selectionW > 0))
				continue;

			const RESTIRDIGridReservoir cand = hg_resSorted[resBase + j];
			const uint slot = RESTIRDICurrentLightSlot(cand.lightIndex);
			[branch]
			if (slot == RESTIR_INVALID_LIGHT_INDEX)
				continue;

			const ShaderEntity candLight = load_entity(slot);
			const float2 candUV = unpack_half2(cand.uvPacked);
			const uint2 candPixel =
				uint2(cand.pixelPacked & 0xFFFF, cand.pixelPacked >> 16);
			const float candDepth = texture_depth[candPixel];
			[branch]
			if (candDepth <= 0)
				continue;

			const float3 Pc = reconstruct_position(
				(candPixel + 0.5) * push.resolutionRcp, candDepth);
			const float3 Nc = decode_normal(texture_normal_roughness[candPixel]);
			const float pOwn = RESTIRDIReuseTarget(candLight, candUV, Pc, Nc);
			[branch]
			if (pOwn <= 0)
				continue;

			const float pCenter = RESTIRDIReuseTarget(candLight, candUV, P, N);
			const float maxW =
				GetFrame().totalLightPower / RESTIRLightPower(candLight);
			const float candW = min(cand.weightSum / (cand.M * pOwn), maxW);
			const float ci = RESTIRHGRequant(cand.M) * candScale;

			const float mi =
				(selectionW / (float)RESTIR_HASHGRID_CANDIDATES) *
				RESTIRHGNoncanonical(ci, centerConfidence, cSum, pOwn, pCenter);
			const float wi = mi * pCenter * candW;
			wSum += wi;
			if (wSum > 0 && rng.next_float() * wSum < wi)
			{
				result.lightIndex = cand.lightIndex;
				result.uv = candUV;
				result.targetPdf = pCenter;
			}
		}
	}

	result.M = (cSum + centerConfidence) /
		(1.0 + (float)RESTIR_HASHGRID_CANDIDATES);
	result.weightSum = (result.targetPdf > 0) ? wSum * result.M : 0;
	result.visibility = 0; // re-traced by the caller
	return result;
}
#endif // RDI_HASHGRID

[numthreads(8, 8, 1)]
void main(uint3 DTid : SV_DispatchThreadID)
{
	if (DTid.x >= push.resolution.x || DTid.y >= push.resolution.y)
		return;

	// Stable light-index translation state (see restir_diHF.hlsli).
	RESTIRDILightMapBuffer = push.lightIndexMapBuffer;
	RESTIRDILightMapOffset = push.lightIndexMapOffset;
	RESTIRDISceneLightCount = push.sceneLightCount;

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

#ifndef RDI_HASHGRID
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

			// Visibility-aware reuse: only borrow a neighbor whose chosen light
			// is actually visible from THIS surface. A shadowed pixel's lit
			// neighbors hold a light that is occluded here; because the merge
			// weights by the UNSHADOWED target, borrowing it lets the merge
			// pick it, the final re-trace then zeroes it, and the resulting
			// black frames darken the shadow boundary - the initial pass does
			// shadowed RIS (no occluded samples, no compensating inflation), so
			// the reuse must reject occluded borrows to match, otherwise the
			// boundary reads dark. This costs one visibility ray per accepted
			// neighbor.
			[branch]
			if (neighbor.lightIndex != RESTIR_INVALID_LIGHT_INDEX &&
				(neighbor.lightIndex & RESTIR_LIGHT_FLAG_EMISSIVE_TRIANGLE) == 0)
			{
				float3 nPos;
				float3 nRad;
				RESTIRDIResolveSample(
					neighbor.lightIndex, neighbor.uv, P, N, nPos, nRad);
				[branch]
				if (RESTIRDITraceVisibility(P, N, nPos, rng) <= 0)
					continue;
			}

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

		[branch]
		if (push.biasCorrection == 1)
			reservoir = RESTIRDIMergePairwiseMIS(
				sources, sourceCount, push.spatialSampleCount, P, N, rng);
		else
			reservoir =
				RESTIRDIMergeBalanceHeuristic(sources, sourceCount, P, N, rng);
#else
		// Hash-grid Stochastic Pairwise MIS: reuse from many surface-similar
		// reservoirs gathered into this pixel's cell (see the build passes).
		reservoir = RESTIRDIMergeHashGrid(pixel, P, N, reservoir, rng);
#endif // RDI_HASHGRID

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
			// alone cannot represent. W is firefly-clamped per light to
			// totalLightPower / this light's power: under power-proportional
			// sampling a dim light's unbiased weight legitimately reaches that
			// bound, so the old light-count clamp (a uniform-sampling bound)
			// darkened dim lights. This caps each light's contribution at
			// ~totalLightPower while still bounding fireflies.
			const uint lightSlot =
				RESTIRDICurrentLightSlot(reservoir.lightIndex);
			const float lightPower =
				(lightSlot != RESTIR_INVALID_LIGHT_INDEX)
					? RESTIRLightPower(load_entity(lightSlot))
					: 1.0;
			const float maxW = GetFrame().totalLightPower / lightPower;
			const float Wclamp = min(W, maxW);
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
						const uint prevSlot =
							RESTIRDICurrentLightSlot(prevRes.lightIndex);
						const float prevPower =
							(prevSlot != RESTIR_INVALID_LIGHT_INDEX)
								? RESTIRLightPower(load_entity(prevSlot))
								: 1.0;
						const float Wprev = min(
							RESTIRDIReservoirGetInvPdf(prevRes),
							GetFrame().totalLightPower / prevPower);
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
