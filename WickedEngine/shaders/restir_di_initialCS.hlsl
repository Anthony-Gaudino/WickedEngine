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
 * Draws candidateCount lights uniformly, resolves each to a world-space sample
 * (via the light-sampling core), and resamples by the diffuse target function.
 * The chosen sample's visibility is left unset for the caller to fill after a
 * shadow ray.
 *
 * @param[in] P - World-space shading point.
 * @param[in] N - World-space shading normal.
 * @param[in] candidateCount - Number of RIS candidates (>= 1).
 * @param[in,out] rng - Random generator.
 *
 * @return The initial reservoir (visibility unset).
 */
inline RESTIRDIReservoir RESTIRDISampleInitial(
	float3 P,
	float3 N,
	uint candidateCount,
	inout RNG rng)
{
	RESTIRDIReservoir r = RESTIRDIReservoirInit();

	const ShaderEntityIterator iterator = lights();
	const uint lightCount = iterator.item_count();
	if (lightCount == 0)
		return r;

	const float invSourcePdf = (float)lightCount;

	for (uint i = 0; i < candidateCount; ++i)
	{
		const uint lightIndex = iterator.first_item() + rng.next_uint(lightCount);
		const ShaderEntity light = load_entity(lightIndex);

		const RESTIRLightSample s = RESTIRResolveAnalyticLight(light, P, N, rng);

		// Resolve to a world-space point on the light. Directional samples have
		// no finite position, so anchor them far along the sample direction.
		const float sampleDist = (s.distance >= FLT_MAX * 0.5) ? 100000.0 : s.distance;
		const float3 samplePosition = P + s.direction * sampleDist;

		const float NdotL = saturate(dot(s.direction, N));
		const float luma = dot(s.radiance, float3(0.2126, 0.7152, 0.0722));
		const float targetPdf = luma * NdotL;
		const float risWeight = targetPdf * invSourcePdf;

		RESTIRDIReservoirUpdate(
			r, samplePosition, s.radiance, lightIndex, targetPdf, risWeight,
			rng);
	}

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

		reservoir = RESTIRDISampleInitial(P, N, push.candidateCount, rng);

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
		}
	}

	RESTIRDIReservoirPacked packed;
	packed.store(reservoir);
	reservoirOutput.Store4(flatIndex * 32, packed.data0);
	reservoirOutput.Store4(flatIndex * 32 + 16, packed.data1);
}
