#include "globals.hlsli"
#include "restir_di_hashgridHF.hlsli"

/**
 * ReSTIR DI hash-grid build - compute per-cell offsets.
 *
 * Second of the three build passes. Assigns each cell a contiguous range in the
 * sorted buffers by adding its element count to a running cursor with an
 * atomic; cells therefore get disjoint ranges regardless of processing order.
 * Runs for both the reservoir map and the inverse map. The `counters` cursors
 * ([1] and [3]) must be zero before this pass (cleared each frame).
 */

PUSHCONSTANT(push, RESTIRHashGridPushConstants);

RWStructuredBuffer<uint> cellCounters : register(u0); // [cell] res, [+N] inv
RWStructuredBuffer<uint> cellOffsets : register(u1);  // [cell] res, [+N] inv
RWStructuredBuffer<uint> counters : register(u2);     // [1]resCursor [3]invCursor

[numthreads(64, 1, 1)]
void main(uint3 DTid : SV_DispatchThreadID)
{
	const uint cell = DTid.x;
	if (cell >= push.cellCount)
		return;

	uint resOffset;
	InterlockedAdd(counters[1], cellCounters[cell], resOffset);
	cellOffsets[cell] = resOffset;

	uint invOffset;
	InterlockedAdd(counters[3], cellCounters[push.cellCount + cell], invOffset);
	cellOffsets[push.cellCount + cell] = invOffset;
}
