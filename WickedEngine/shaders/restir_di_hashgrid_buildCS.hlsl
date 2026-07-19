#include "globals.hlsli"
#include "restir_diHF.hlsli"
#include "restir_di_hashgridHF.hlsli"

/**
 * ReSTIR DI hash-grid build - create reuse cells.
 *
 * First of the three hash-grid build passes. For each pixel it finds (or
 * creates) its surface cell, adds its reservoir's confidence to the cell total,
 * inserts the pixel into the cell's **inverse map** (all valid-surface pixels,
 * for the canonical MIS term), and - when the reservoir holds a usable analytic
 * light sample - inserts a compact reservoir into the cell's **reservoir map**
 * (the reuse candidates). The offsets and sort passes then compact both maps
 * into contiguous per-cell ranges.
 */

PUSHCONSTANT(push, RESTIRHashGridPushConstants);

// Source reservoirs (this frame's temporal result), read as raw.
ByteAddressBuffer reservoirInput : register(t0);

RWStructuredBuffer<uint> checksums : register(u0);       // hash set
RWStructuredBuffer<uint> cellConfidences : register(u1); // atomic M sum per cell
RWStructuredBuffer<uint> counters : register(u2);        // [0]resTot [2]invTot
RWStructuredBuffer<uint> cellCounters : register(u3);    // [cell] res, [+N] inv
RWStructuredBuffer<RESTIRDIGridReservoir> resData : register(u4);
RWStructuredBuffer<uint2> resDataCell : register(u5);    // (cell, indexInCell)
RWStructuredBuffer<uint2> invData : register(u6);        // (flatPixel, asuint(M))
RWStructuredBuffer<uint2> invDataCell : register(u7);    // (cell, indexInCell)

/**
 * Quantizes a confidence to an integer so it can be summed with atomics.
 *
 * @param[in] c - Confidence (M).
 *
 * @return The quantized confidence.
 */
inline uint RESTIRQuantizeConfidence(float c)
{
	return (uint)(c * 1000.0);
}

[numthreads(8, 8, 1)]
void main(uint3 DTid : SV_DispatchThreadID)
{
	if (DTid.x >= push.resolution.x || DTid.y >= push.resolution.y)
		return;

	const uint2 pixel = DTid.xy;

	const float depth = texture_depth[pixel];
	[branch]
	if (depth <= 0) // sky / background: not a surface, no cell
		return;

	const uint flatIndex = pixel.y * push.resolution.x + pixel.x;
	const RESTIRDIReservoir reservoir =
		RESTIRDIReservoirLoad(reservoirInput, flatIndex);

	const float3 N = decode_normal(texture_normal_roughness[pixel]);
	const float linearDepth = compute_lineardepth(depth);

	const RESTIRHashKey key = RESTIRHashGridKey(pixel, N, linearDepth);
	bool first;
	const uint cell =
		RESTIRHashSetInsert(checksums, push.cellCount, key, first);
	[branch]
	if (cell == RESTIR_HASHGRID_INVALID)
		return; // probe overflow: skip this pixel this frame

	const float confidence = reservoir.M;

	// Cell confidence total (all valid-surface pixels).
	InterlockedAdd(cellConfidences[cell], RESTIRQuantizeConfidence(confidence));

	// Inverse map: this pixel + its confidence (for the canonical MIS term).
	{
		uint slot;
		InterlockedAdd(counters[2], 1u, slot);
		if (slot < push.cellCount)
		{
			uint indexInCell;
			InterlockedAdd(cellCounters[push.cellCount + cell], 1u, indexInCell);
			invData[slot] = uint2(flatIndex, asuint(confidence));
			invDataCell[slot] = uint2(cell, indexInCell);
		}
	}

	// Reservoir map: only pixels whose reservoir holds a usable analytic light
	// sample become reuse candidates.
	const bool usable = reservoir.M > 0 && reservoir.targetPdf > 0 &&
		reservoir.weightSum > 0 &&
		reservoir.lightIndex != RESTIR_INVALID_LIGHT_INDEX &&
		(reservoir.lightIndex & RESTIR_LIGHT_FLAG_EMISSIVE_TRIANGLE) == 0;
	[branch]
	if (usable)
	{
		uint slot;
		InterlockedAdd(counters[0], 1u, slot);
		if (slot < push.cellCount)
		{
			uint indexInCell;
			InterlockedAdd(cellCounters[cell], 1u, indexInCell);

			RESTIRDIGridReservoir g;
			g.lightIndex = reservoir.lightIndex;
			g.uvPacked = pack_half2(reservoir.uv.x, reservoir.uv.y);
			g.pixelPacked = (pixel.x & 0xFFFF) | (pixel.y << 16);
			g.M = reservoir.M;
			g.weightSum = reservoir.weightSum;

			resData[slot] = g;
			resDataCell[slot] = uint2(cell, indexInCell);
		}
	}
}
