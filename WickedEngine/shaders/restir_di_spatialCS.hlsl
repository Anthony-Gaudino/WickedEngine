#include "globals.hlsli"
#include "restir_diHF.hlsli"

// ReSTIR DI - spatial resampling pass.
//
// Merges a handful of nearby reservoirs into each pixel's temporal reservoir,
// rejecting neighbors that lie on a different surface (depth or normal
// mismatch). The result is the final reservoir: it is consumed by forward
// shading this frame and becomes the temporal history for the next frame.

PUSHCONSTANT(push, RESTIRDIPushConstants);

ByteAddressBuffer reservoirInput : register(t0);

RWByteAddressBuffer reservoirOutput : register(u0);

[numthreads(8, 8, 1)]
void main(uint3 DTid : SV_DispatchThreadID)
{
	if (DTid.x >= push.resolution.x || DTid.y >= push.resolution.y)
		return;

	const uint2 pixel = DTid.xy;
	const uint flatIndex = pixel.y * push.resolution.x + pixel.x;

	RESTIRDIReservoir reservoir = RESTIRDIReservoirLoad(reservoirInput, flatIndex);

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
	}

	RESTIRDIReservoirStore(reservoirOutput, flatIndex, reservoir);
}
