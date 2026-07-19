#ifndef WI_RESTIR_DI_HF
#define WI_RESTIR_DI_HF
#include "ShaderInterop_ReSTIR.h"

/**
 * ReSTIR DI screen-space reservoir helpers.
 *
 * Operates on the resolved-sample DI reservoir (RESTIRDIReservoir): weighted
 * reservoir streaming, packing, and the merge-source descriptor. Deliberately
 * free of any light-sampling dependency so the light-free consumers (the
 * denoise passes and forward shading) can include it without pulling in
 * lightingHF.hlsli. The spatiotemporal balance-heuristic merge, which
 * re-resolves the analytic light per surface, lives in restir_di_reuseHF.hlsli
 * and is included only by the temporal/spatial passes.
 *
 * Include order: "globals.hlsli" first, then this header.
 */

/*
################################################################################
Stable light-index translation
################################################################################
*/

/**
 * Bindless SRV of the packed stable light-index translation buffer, -1 for
 * identity.
 *
 * A reservoir stores a **stable scene light index**; the per-frame entity array
 * is reindexed by frustum culling, so this maps between the stable index and
 * the current entity slot (see RESTIRDIPushConstants::lightIndexMapBuffer for
 * the layout). Set once per thread at the top of each pass's `main()`; the
 * resolve helpers read it. Only touches `bindless_buffers`, never
 * `load_entity`, so the light-free consumers (denoise passes, forward shading)
 * can still include this header.
 */
static int RESTIRDILightMapBuffer = -1;

/** Byte offset of the translation buffer's data within its GPU buffer. */
static uint RESTIRDILightMapOffset = 0;

/** Number of scene lights: bound + split point of the two map sub-arrays. */
static uint RESTIRDISceneLightCount = 0;

/**
 * Maps a per-frame entity slot to its stable scene light index (store side).
 *
 * @param[in] slot - Entity-array slot of a currently-visible analytic light.
 *
 * @return The light's stable scene index, or `slot` unchanged when translation
 *   is disabled (RESTIRDILightMapBuffer < 0).
 */
inline uint RESTIRDIStableLightIndex(uint slot)
{
	[branch]
	if (RESTIRDILightMapBuffer < 0)
		return slot;
	return bindless_buffers[descriptor_index(RESTIRDILightMapBuffer)].Load(
		RESTIRDILightMapOffset + (RESTIRDISceneLightCount + slot) * 4);
}

/**
 * Maps a stable scene light index back to its current entity slot (reuse side).
 *
 * @param[in] sceneIndex - Stable scene light index stored in a reservoir.
 *
 * @return The current entity slot, RESTIR_INVALID_LIGHT_INDEX if the light is
 *   culled (or out of range) this frame, or `sceneIndex` unchanged when
 *   translation is disabled.
 */
inline uint RESTIRDICurrentLightSlot(uint sceneIndex)
{
	[branch]
	if (RESTIRDILightMapBuffer < 0)
		return sceneIndex;
	[branch]
	if (sceneIndex >= RESTIRDISceneLightCount)
		return RESTIR_INVALID_LIGHT_INDEX;
	return bindless_buffers[descriptor_index(RESTIRDILightMapBuffer)].Load(
		RESTIRDILightMapOffset + sceneIndex * 4);
}

/**
 * Streams one resolved candidate into a DI reservoir (weighted reservoir
 * sampling).
 *
 * @param[in,out] r - Reservoir to update.
 * @param[in] samplePosition - World-space point on the candidate light.
 * @param[in] sampleRadiance - Unshadowed incident radiance of the candidate.
 * @param[in] lightIndex - Index of the candidate light.
 * @param[in] uv - Sample parameterization on the light (for re-resolution).
 * @param[in] targetPdf - Target function value p_hat of the candidate.
 * @param[in] risWeight - Resampling weight p_hat / p_source of the candidate.
 * @param[in,out] rng - Random generator.
 *
 * @return true if the candidate became the reservoir's selected sample.
 */
inline bool RESTIRDIReservoirUpdate(
	inout RESTIRDIReservoir r,
	float3 samplePosition,
	float3 sampleRadiance,
	uint lightIndex,
	float2 uv,
	float targetPdf,
	float risWeight,
	inout RNG rng)
{
	r.weightSum += risWeight;
	r.M += 1;

	if (r.weightSum > 0 && rng.next_float() * r.weightSum < risWeight)
	{
		r.samplePosition = samplePosition;
		r.sampleRadiance = sampleRadiance;
		r.lightIndex = lightIndex;
		r.uv = uv;
		r.targetPdf = targetPdf;
		return true;
	}
	return false;
}

/**
 * Maximum number of reservoirs the balance-heuristic merge combines at once:
 * the pixel's own (canonical) reservoir plus its spatial neighbors.
 */
#define RESTIR_DI_MIS_MAX_SOURCES (RESTIR_DI_SPATIAL_SAMPLES + 1)

/**
 * One input reservoir of a spatiotemporal merge, tagged with the surface it was
 * created on.
 *
 * The surface (P, N) is what makes the merge unbiased: the MIS weight of each
 * source's sample is evaluated against every source's **own** surface (see
 * RESTIRDIMergePairwiseMIS), so a neighbor that cannot see the chosen
 * light no longer deflates the estimate.
 */
struct RESTIRDIMISSource
{
	/** World-space point on the chosen light. */
	float3 samplePosition;

	/** Unshadowed incident radiance of the sample. */
	float3 sampleRadiance;

	/** World-space shading point of the surface this reservoir was built on. */
	float3 P;

	/** World-space shading normal of that surface. */
	float3 N;

	/** Confidence (effective candidate count) of the reservoir. */
	float M;

	/** Running resampling-weight sum of the reservoir. */
	float weightSum;

	/** Cached visibility of the held sample, in [0, 1]. */
	float visibility;

	/** Index of the chosen light. */
	uint lightIndex;

	/** Sample parameterization on the light's area, in [0, 1)^2. */
	float2 uv;
};

/**
 * Loads a DI reservoir from a raw reservoir buffer.
 *
 * The buffers are raw (ByteAddressBuffer) so the same resource can be read both
 * by the DI passes and bindlessly by forward shading. Each entry is 48 bytes.
 *
 * @param[in] buf - Reservoir byte-address buffer.
 * @param[in] index - Flat pixel index (y * width + x).
 *
 * @return The unpacked reservoir at that index.
 */
inline RESTIRDIReservoir RESTIRDIReservoirLoad(ByteAddressBuffer buf, uint index)
{
	RESTIRDIReservoirPacked p;
	p.data0 = buf.Load4(index * 48);
	p.data1 = buf.Load4(index * 48 + 16);
	p.data2 = buf.Load4(index * 48 + 32);
	return p.load();
}

/**
 * Loads a DI reservoir from a writable reservoir buffer (in-place read).
 *
 * @param[in] buf - Reservoir byte-address buffer (UAV).
 * @param[in] index - Flat pixel index (y * width + x).
 *
 * @return The unpacked reservoir at that index.
 */
inline RESTIRDIReservoir RESTIRDIReservoirLoad(RWByteAddressBuffer buf, uint index)
{
	RESTIRDIReservoirPacked p;
	p.data0 = buf.Load4(index * 48);
	p.data1 = buf.Load4(index * 48 + 16);
	p.data2 = buf.Load4(index * 48 + 32);
	return p.load();
}

/**
 * Stores a DI reservoir into a raw reservoir buffer.
 *
 * @param[in,out] buf - Reservoir byte-address buffer (UAV).
 * @param[in] index - Flat pixel index (y * width + x).
 * @param[in] r - The reservoir to store.
 */
inline void RESTIRDIReservoirStore(
	RWByteAddressBuffer buf, uint index, RESTIRDIReservoir r)
{
	RESTIRDIReservoirPacked p;
	p.store(r);
	buf.Store4(index * 48, p.data0);
	buf.Store4(index * 48 + 16, p.data1);
	buf.Store4(index * 48 + 32, p.data2);
}

#endif // WI_RESTIR_DI_HF
