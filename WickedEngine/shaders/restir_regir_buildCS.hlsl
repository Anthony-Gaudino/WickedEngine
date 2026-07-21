#include "globals.hlsli"
#include "lightingHF.hlsli"
#include "restir_lightsamplingHF.hlsli"

/**
 * ReGIR - grid accumulation pass.
 *
 * The ReGIR grid is anchored to the scene (not the camera), so every cell maps
 * to a fixed world region and its light reservoir **persists across frames**.
 * Each frame this pass streams RESTIR_REGIR_BUILD_SAMPLES uniformly-drawn
 * lights into every cell slot's reservoir, weighting them by the light's
 * contribution to that slot's cell (power / distance^2), and caps the reservoir
 * confidence. Over a few dozen frames each slot converges to a stable,
 * low-variance light importance-sampled for its cell.
 *
 * The stored reciprocal source pdf (weightSum / (M * targetWeight)) is the
 * reservoir's unbiased contribution weight, so the per-pixel RIS in the DI
 * initial pass draws position-aware candidates from the cell containing the
 * pixel and stays unbiased - the RIS-of-RIS construction, position-aware and
 * temporally stable.
 *
 * References: Wyman & Panteleev 2021, "Rearchitecting Spatiotemporal Resampling
 *   for Production" (ReGIR).
 */

PUSHCONSTANT(push, RESTIRReGIRBuildPushConstants);

// Persistent storage: one RESTIRReGIRCell (16 bytes, uint4) per slot. Read and
// written in place - each thread only touches its own slot, so there is no
// race.
RWByteAddressBuffer regir : register(u0);

[numthreads(64, 1, 1)]
void main(uint3 DTid : SV_DispatchThreadID)
{
	const uint slot = DTid.x;
	if (slot >= RESTIR_REGIR_SLOT_TOTAL)
		return;

	const int cellIndex = (int)(slot / RESTIR_REGIR_LIGHTS_PER_CELL);
	int3 worldCell;
	float3 cellCenter;
	float cellRadius;
	RESTIRReGIRCellToWorld(
		cellIndex, push.gridOriginCell, push.cellSize,
		worldCell, cellCenter, cellRadius);

	// Load this slot's persistent reservoir (32-byte cell: reservoir in the
	// first uint4, toroidal world-cell key in the second).
	const uint4 raw0 = regir.Load4(slot * 32);
	const uint4 raw1 = regir.Load4(slot * 32 + 16);
	uint lightIndex = raw0.x;
	float weightSum = asfloat(raw0.y);
	float M = asfloat(raw0.z);
	float targetWeight = asfloat(raw0.w);
	const int3 storedCell = int3(asint(raw1.x), asint(raw1.y), asint(raw1.z));

	// Toroidal invalidation: if the camera-centered window has scrolled so this
	// buffer slot now maps to a different world cell than it last cached, the
	// reservoir belongs to the opposite edge's world region - reset it so the
	// slot starts empty (and the initial pass falls back to the tiles for it)
	// instead of injecting a stale, wrong-location light set.
	[branch]
	if (any(storedCell != worldCell))
	{
		lightIndex = RESTIR_INVALID_LIGHT_INDEX;
		weightSum = 0;
		M = 0;
		targetWeight = 0;
	}

	// A fresh (zero-initialized) slot holds no light yet.
	[branch]
	if (M <= 0)
		lightIndex = RESTIR_INVALID_LIGHT_INDEX;

	// Cap confidence so stale samples decay and a light that moved is followed
	// within a fraction of a second (scale weightSum with M so W is preserved).
	[branch]
	if (M > RESTIR_REGIR_MAX_M)
	{
		weightSum *= RESTIR_REGIR_MAX_M / M;
		M = RESTIR_REGIR_MAX_M;
	}

	const ShaderEntityIterator iterator = lights();
	const uint lightCount = iterator.item_count();

	[branch]
	if (lightCount > 0)
	{
		// RNG::init packs id as (id.x << 16) | id.y, so a slot index that needs
		// more than 16 bits (RESTIR_REGIR_SLOT_TOTAL is ~2^21 here) would lose
		// its high bits and collide. Those high bits are exactly the cell's Z
		// column (z contributes slot >> 16), so a naive uint2(slot, 0) seed
		// makes every cell in a Z column share the same candidate stream. Split
		// the full slot across both id components so all bits survive.
		RNG rng;
		rng.init(uint2(slot & 0xFFFFu, slot >> 16u), push.frameIndex);

		// Uniform source pmf for this frame's streamed candidates.
		const float invUniformPdf = (float)lightCount;

		for (uint i = 0; i < RESTIR_REGIR_BUILD_SAMPLES; ++i)
		{
			const uint candidate =
				iterator.first_item() + rng.next_uint(lightCount);
			const float weight = RESTIRReGIRVolumeWeight(
				load_entity(candidate), cellCenter, cellRadius);

			// risWeight = target / source = weight / (1 / lightCount).
			const float risWeight = weight * invUniformPdf;
			weightSum += risWeight;
			M += 1;

			if (weightSum > 0 && rng.next_float() * weightSum < risWeight)
			{
				lightIndex = candidate;
				targetWeight = weight;
			}
		}
	}

	regir.Store4(
		slot * 32,
		uint4(lightIndex, asuint(weightSum), asuint(M), asuint(targetWeight)));
	// Stamp the slot's current toroidal identity so a later scroll can detect
	// that it has gone stale.
	regir.Store4(
		slot * 32 + 16,
		uint4(asuint(worldCell.x), asuint(worldCell.y), asuint(worldCell.z), 0));
}
