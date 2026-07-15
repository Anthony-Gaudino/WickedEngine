#define TEXTURE_SLOT_NONUNIFORM // shadow ray alpha-tests divergent materials
#include "globals.hlsli"
#include "raytracingHF.hlsli"
#include "lightingHF.hlsli"
#include "restir_lightsamplingHF.hlsli"
#include "restir_diHF.hlsli"
#include "restir_di_visibilityHF.hlsli"

// ReSTIR DI - initial candidate generation pass.
//
// Per screen pixel: reconstruct the primary surface from the depth + normal
// buffers, resample RESTIRDIPushConstants.candidateCount analytic lights into
// an initial reservoir, then trace a single shadow ray for the chosen sample
// and cache its visibility. The reservoir is written to the temporal history
// buffer, which the temporal pass consumes next.

PUSHCONSTANT(push, RESTIRDIPushConstants);

// Reservoir buffers are raw (ByteAddressBuffer) so the same resource can also
// be read bindlessly by the forward shading pass. Each reservoir is 32 bytes
// (RESTIRDIReservoirPacked = 2x uint4).
RWByteAddressBuffer reservoirOutput : register(u0);

/**
 * Builds an initial DI reservoir by resampling the analytic light list.
 *
 * Draws candidateCount lights (from a pre-sampled tile when available, else
 * uniformly), resolves each to a world-space sample via the light-sampling
 * core, and resamples by the diffuse target function. The chosen sample's
 * visibility is left unset for the caller to fill after a shadow ray.
 *
 * @param[in] P - World-space shading point.
 * @param[in] N - World-space shading normal.
 * @param[in] candidateCount - Number of RIS candidates (>= 1).
 * @param[in] pixel - Screen pixel (selects the pre-sampled tile).
 * @param[in,out] rng - Random generator.
 *
 * @return The initial reservoir (visibility unset).
 */
inline RESTIRDIReservoir RESTIRDISampleInitial(
	float3 P,
	float3 N,
	uint candidateCount,
	uint2 pixel,
	inout RNG rng)
{
	RESTIRDIReservoir r = RESTIRDIReservoirInit();

	RESTIRLightSample winning;
	const RESTIRReservoir lr = RESTIRSampleLights(
		P, N, candidateCount, push.lightTileBuffer, pixel, push.frameIndex,
		false, rng, winning);

	if (lr.lightIndex == RESTIR_INVALID_LIGHT_INDEX)
		return r;

	// Resolve the winner to a world-space point on the light. Directional
	// samples have no finite position, so anchor them far along the direction.
	const float sampleDist =
		(winning.distance >= FLT_MAX * 0.5) ? 100000.0 : winning.distance;

	r.samplePosition = P + winning.direction * sampleDist;
	r.sampleRadiance = winning.radiance;
	r.targetPdf = lr.targetPdf;
	r.weightSum = lr.weightSum;
	r.M = lr.M;
	r.lightIndex = lr.lightIndex;
	r.uv = lr.uv;
	return r;
}

[numthreads(8, 8, 1)]
void main(uint3 DTid : SV_DispatchThreadID)
{
	if (DTid.x >= push.resolution.x || DTid.y >= push.resolution.y)
		return;

	const uint2 pixel = DTid.xy;
	const uint flatIndex = pixel.y * push.resolution.x + pixel.x;

	RESTIRDIReservoir reservoir = RESTIRDIReservoirInit();

	const float depth = texture_depth[pixel];
	[branch]
	if (depth > 0) // skip the sky
	{
		const float2 uv = (pixel + 0.5) * push.resolutionRcp;
		const float3 P = reconstruct_position(uv, depth);
		const float3 N = decode_normal(texture_normal_roughness[pixel]);

		RNG rng;
		rng.init(pixel, GetFrame().frame_count);

		reservoir =
			RESTIRDISampleInitial(P, N, push.candidateCount, pixel, rng);

		// Cache visibility for the chosen sample. This is used directly when
		// spatiotemporal reuse is off; with reuse on, the spatial pass
		// re-traces a fresh shadow ray for the final reused sample so the
		// shadow reflects the current frame instead of a stale cached value.
		const float W = RESTIRDIReservoirGetInvPdf(reservoir);
		[branch]
		if (W > 0 && reservoir.targetPdf > 0)
		{
			reservoir.visibility =
				RESTIRDITraceVisibility(P, N, reservoir.samplePosition, rng);

			// RTXDI initial-visibility test: an occluded initial sample carries
			// no weight into reuse, so temporal/spatial resampling converges to
			// the visible light instead of keeping the occluded one and
			// swinging black<->over-bright each frame (see visibilityReject).
			// weightSum=0 makes the merge start empty so a visible
			// history/neighbor wins.
			[branch]
			if (push.visibilityReject != 0 && reservoir.visibility <= 0)
			{
				reservoir.weightSum = 0;
			}
		}
	}

	RESTIRDIReservoirStore(reservoirOutput, flatIndex, reservoir);
}
