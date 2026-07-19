#include "globals.hlsli"
#include "restir_di_hashgridHF.hlsli"

/**
 * ReSTIR DI hash-grid build - sort into per-cell ranges.
 *
 * Third of the three build passes. Scatters each inserted element into its
 * cell's contiguous sorted range (`cellOffsets[cell] + indexInCell`), for both
 * the reservoir map and the inverse map, so the spatial pass can read a cell's
 * candidates as a simple index range.
 */

PUSHCONSTANT(push, RESTIRHashGridPushConstants);

RWStructuredBuffer<RESTIRDIGridReservoir> resData : register(u0);
RWStructuredBuffer<uint2> resDataCell : register(u1);
RWStructuredBuffer<RESTIRDIGridReservoir> resSorted : register(u2);
RWStructuredBuffer<uint2> invData : register(u3);
RWStructuredBuffer<uint2> invDataCell : register(u4);
RWStructuredBuffer<uint2> invSorted : register(u5);
RWStructuredBuffer<uint> counters : register(u6);    // [0]resTot [2]invTot
RWStructuredBuffer<uint> cellOffsets : register(u7); // [cell] res, [+N] inv

[numthreads(64, 1, 1)]
void main(uint3 DTid : SV_DispatchThreadID)
{
	const uint i = DTid.x;

	[branch]
	if (i < counters[0]) // reservoir map
	{
		const uint2 cd = resDataCell[i];
		resSorted[cellOffsets[cd.x] + cd.y] = resData[i];
	}

	[branch]
	if (i < counters[2]) // inverse map
	{
		const uint2 cd = invDataCell[i];
		invSorted[cellOffsets[push.cellCount + cd.x] + cd.y] = invData[i];
	}
}
