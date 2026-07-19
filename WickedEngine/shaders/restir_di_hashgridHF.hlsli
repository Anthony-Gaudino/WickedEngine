#ifndef WI_RESTIR_DI_HASHGRID_HF
#define WI_RESTIR_DI_HASHGRID_HF
#include "ShaderInterop_ReSTIR.h"

/**
 * ReSTIR DI hash-grid reuse structures.
 *
 * A screen-space hash grid that bins each pixel's reservoir into a cell keyed
 * by its tile and quantized surface (normal + depth), so the spatial pass can
 * reuse from many surface-similar reservoirs (large-kernel Stochastic Pairwise
 * MIS) instead of a few random screen neighbors. Ported from the reference
 * implementation's HashMap.slang (Hedstrom et al. 2026).
 *
 * The grid is built each frame in three passes:
 *   1. create cells - each pixel inserts its reservoir into its cell and adds
 *      its confidence to the cell total (atomics),
 *   2. compute offsets - prefix-sum the per-cell counts into buffer offsets,
 *   3. sort - scatter reservoirs into contiguous per-cell ranges.
 *
 * HLSL has no resource-holding structs, so the map is a set of free functions
 * that take the backing buffers as parameters. The buffers are:
 *   - `checksums`      : one uint per cell (open-addressing hash set),
 *   - `cellCounters`   : one uint per cell (elements in the cell),
 *   - `cellOffsets`    : one uint per cell (start of the cell's sorted range),
 *   - `counters`       : [0] = total inserted, [1] = running offset cursor,
 *   - `data`/`dataCell`: unsorted reservoirs + their (cell, indexInCell),
 *   - `sortedData`     : reservoirs compacted into per-cell ranges.
 *
 * References:
 * Hedstrom et al. 2026, "Stochastic Pairwise MIS for Unbiased Large-Kernel
 * Reuse in Real Time". xxHash (Cyan4973); PCG (pcg-random.org).
 */

/*
################################################################################
Hashing
################################################################################
*/

/**
 * 32-bit xxHash of a single word (used for the collision-checking checksum).
 *
 * @param[in] p - Input word.
 *
 * @return The hashed word.
 */
inline uint RESTIRHashXX32(uint p)
{
	const uint PRIME32_2 = 2246822519u;
	const uint PRIME32_3 = 3266489917u;
	const uint PRIME32_4 = 668265263u;
	const uint PRIME32_5 = 374761393u;
	uint h32 = p + PRIME32_5;
	h32 = PRIME32_4 * ((h32 << 17) | (h32 >> (32 - 17)));
	h32 = PRIME32_2 * (h32 ^ (h32 >> 15));
	h32 = PRIME32_3 * (h32 ^ (h32 >> 13));
	return h32 ^ (h32 >> 16);
}

/**
 * PCG hash of a single word (used for the cell-address key).
 *
 * @param[in] v - Input word.
 *
 * @return The hashed word.
 */
inline uint RESTIRHashPCG(uint v)
{
	uint state = v * 747796405u + 2891336453u;
	uint word = ((state >> ((state >> 28u) + 4u)) ^ state) * 277803737u;
	return (word >> 22u) ^ word;
}

/**
 * A hash key: an address `key` (which cell) and a `checksum` (collision guard).
 *
 * Attributes are folded in one at a time with RESTIRHashKeyAdd so two different
 * surfaces that land on the same address are still distinguished by checksum.
 */
struct RESTIRHashKey
{
	uint key;
	uint checksum;
};

/**
 * Starts a hash key from a single attribute value.
 *
 * @param[in] value - The first attribute.
 *
 * @return A hash key seeded with the value.
 */
inline RESTIRHashKey RESTIRHashKeyInit(uint value)
{
	RESTIRHashKey k;
	k.key = value;
	k.checksum = value;
	return k;
}

/**
 * Folds another attribute into a hash key.
 *
 * @param[in] lhs - The accumulated key.
 * @param[in] value - The attribute to fold in.
 *
 * @return The combined key.
 */
inline RESTIRHashKey RESTIRHashKeyAdd(RESTIRHashKey lhs, uint value)
{
	RESTIRHashKey r;
	r.key = RESTIRHashPCG(lhs.key + RESTIRHashPCG(value));
	r.checksum = RESTIRHashXX32(lhs.checksum + RESTIRHashXX32(value));
	return r;
}

/*
################################################################################
Cell hash set (open addressing)
################################################################################
*/

/**
 * Inserts (or finds) the cell for a key, returning its slot index.
 *
 * Open-addressing linear probe with a checksum stored per slot: the first
 * inserter of a slot wins it (atomic compare-exchange from INVALID), later
 * lookups with the same checksum land on the same slot.
 *
 * @param[in,out] checksums - Per-cell checksum buffer (size `cellCount`).
 * @param[in] cellCount - Number of cells.
 * @param[in] key - The cell key.
 * @param[out] first - true if this call created the cell.
 *
 * @return The cell slot index, or RESTIR_HASHGRID_INVALID if the probe overran.
 */
inline uint RESTIRHashSetInsert(
	RWStructuredBuffer<uint> checksums,
	uint cellCount,
	RESTIRHashKey key,
	out bool first)
{
	first = false;
	const uint checksum = min(key.checksum, RESTIR_HASHGRID_INVALID - 1);
	uint cell = key.key % cellCount;
	[loop]
	for (uint i = 0; i < RESTIR_HASHGRID_PROBE_LENGTH; ++i)
	{
		uint current;
		InterlockedCompareExchange(
			checksums[cell], RESTIR_HASHGRID_INVALID, checksum, current);
		if (current == checksum || current == RESTIR_HASHGRID_INVALID)
		{
			first = (current == RESTIR_HASHGRID_INVALID);
			return cell;
		}
		cell = (cell + 1) % cellCount;
	}
	return RESTIR_HASHGRID_INVALID;
}

/**
 * Finds the cell slot for a key (read-only), or RESTIR_HASHGRID_INVALID.
 *
 * @param[in] checksums - Per-cell checksum buffer.
 * @param[in] cellCount - Number of cells.
 * @param[in] key - The cell key.
 *
 * @return The cell slot index, or RESTIR_HASHGRID_INVALID if not present.
 */
inline uint RESTIRHashSetFind(
	StructuredBuffer<uint> checksums, uint cellCount, RESTIRHashKey key)
{
	const uint checksum = min(key.checksum, RESTIR_HASHGRID_INVALID - 1);
	uint cell = key.key % cellCount;
	[loop]
	for (uint i = 0; i < RESTIR_HASHGRID_PROBE_LENGTH; ++i)
	{
		const uint current = checksums[cell];
		if (current == RESTIR_HASHGRID_INVALID)
			return RESTIR_HASHGRID_INVALID;
		if (current == checksum)
			return cell;
		cell = (cell + 1) % cellCount;
	}
	return RESTIR_HASHGRID_INVALID;
}

/*
################################################################################
Cell key from a pixel's surface
################################################################################
*/

/**
 * Builds the hash key of a pixel: its screen tile plus a quantized surface
 * (octahedral normal + linear depth bucket), so only reservoirs on a similar
 * surface share a cell.
 *
 * @param[in] pixel - Screen pixel.
 * @param[in] normal - World-space normal at the pixel.
 * @param[in] linearDepth - Linear (view-space) depth at the pixel.
 *
 * @return The pixel's cell key.
 */
inline RESTIRHashKey RESTIRHashGridKey(
	uint2 pixel, float3 normal, float linearDepth)
{
	const int3 qn = int3(floor(normal * 4.0));
	// Log2 depth buckets keep near geometry finely separated and far geometry
	// coarsely, so a cell never straddles a large depth discontinuity.
	const int qd = int(floor(log2(max(linearDepth, 1e-3)) * 8.0));

	RESTIRHashKey k = RESTIRHashKeyInit(pixel.x / RESTIR_HASHGRID_CELL_SIZE);
	k = RESTIRHashKeyAdd(k, pixel.y / RESTIR_HASHGRID_CELL_SIZE);
	k = RESTIRHashKeyAdd(k, (uint)(qn.x + 8));
	k = RESTIRHashKeyAdd(k, (uint)(qn.y + 8));
	k = RESTIRHashKeyAdd(k, (uint)(qn.z + 8));
	k = RESTIRHashKeyAdd(k, (uint)(qd + 4096));
	return k;
}

#endif // WI_RESTIR_DI_HASHGRID_HF
