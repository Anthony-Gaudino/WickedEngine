#define TEXTURE_SLOT_NONUNIFORM // bounce hit + shadow ray load divergent materials
#include "globals.hlsli"
#include "raytracingHF.hlsli"
#include "lightingHF.hlsli"
#include "restir_lightsamplingHF.hlsli"
#include "restir_di_visibilityHF.hlsli"
#include "restir_giHF.hlsli"

/**
 * ReSTIR GI - initial sample generation pass.
 *
 * Per screen pixel: reconstruct the primary surface (P from depth, N from the
 * normal buffer), trace one cosine-weighted hemisphere ray, and shade the hit
 * as the radiance leaving it toward P (its direct lighting via the shared RIS
 * core + one shadow ray, plus emission; a miss reads the sky). That sample -
 * position, normal and outgoing radiance - is streamed into an initial GI
 * reservoir. The reservoir is written for the (future) spatiotemporal reuse
 * passes, and this frame's resolved albedo-free indirect irradiance is written
 * for the denoiser.
 *
 * This is the initial stage: no spatiotemporal reuse yet, so one bounce is
 * sampled per pixel and the result is noisy before the denoiser. Temporal and
 * spatial reuse (with the reconnection Jacobian) are added next.
 */

PUSHCONSTANT(push, RESTIRGIPushConstants);

// RIS candidate count for shading the bounce hit's direct lighting.
static const uint RESTIR_GI_LIGHT_CANDIDATES = 8;

// Initial per-pixel reservoir; the temporal pass merges history into it and
// resolves the denoiser inputs.
RWByteAddressBuffer reservoirOutput : register(u0);

/**
 * Shades a bounce-ray hit surface: its direct lighting plus emission, i.e. the
 * radiance leaving the hit toward the ray origin.
 *
 * @param[in,out] surface - The loaded, updated hit surface.
 * @param[in,out] rng - Random generator.
 *
 * @return Outgoing radiance from the hit (RGB, >= 0).
 */
float3 RESTIRGIShadeHit(inout Surface surface, inout RNG rng)
{
	float3 direct = 0;

	// Direct lighting at the hit (ReSTIR RIS over the analytic light list) with
	// one shadow ray for the winner - mirrors the DDGI/Surfel bounce shading.
	RESTIRLightSample winning;
	const RESTIRReservoir reservoir = RESTIRSampleLightsUniform(
		surface.P, surface.N, RESTIR_GI_LIGHT_CANDIDATES, false, rng, winning);

	const float W = RESTIRReservoirGetInvPdf(reservoir);
	const float NdotL = saturate(dot(winning.direction, surface.N));

	[branch]
	if (W > 0 && NdotL > 0 && winning.distance > 0)
	{
		const float dist =
			(winning.distance >= FLT_MAX * 0.5) ? 100000.0 : winning.distance;
		const float3 samplePosition = surface.P + winning.direction * dist;
		const float shadow = RESTIRDITraceVisibility(
			surface.P, surface.N, samplePosition, rng);
		direct += shadow * winning.radiance * NdotL * W / PI;
	}

	// Lambertian outgoing radiance = albedo * irradiance + emission.
	return max(0, direct) * surface.albedo + surface.emissiveColor;
}

[numthreads(8, 8, 1)]
void main(uint3 DTid : SV_DispatchThreadID)
{
	if (DTid.x >= push.resolution.x || DTid.y >= push.resolution.y)
		return;

	const uint2 pixel = DTid.xy;
	const uint flatIndex = pixel.y * push.resolution.x + pixel.x;

	RESTIRGIReservoir reservoir = RESTIRGIReservoirInit();

	const float depth = texture_depth[pixel];
	[branch]
	if (depth > 0) // skip the sky
	{
		const float2 uv = (pixel + 0.5) * push.resolutionRcp;
		const float3 P = reconstruct_position(uv, depth);
		const float3 N = decode_normal(texture_normal_roughness[pixel]);

		RNG rng;
		rng.init(pixel, GetFrame().frame_count);

		// Normal-offset the bounce-ray origin (scaled with camera distance) to
		// avoid self-intersection at grazing angles, matching the shadow ray.
		const float distToCamera = length(P - GetCamera().position);
		const float bias = clamp(0.002 * distToCamera, 0.02, 0.5);
		const float3 rayOrigin = P + N * bias;

		for (uint i = 0; i < push.candidateCount; ++i)
		{
			// Cosine-weighted hemisphere direction around N.
			const float2 u = float2(rng.next_float(), rng.next_float());
			const float3 dir = RESTIRHemisphereCos(N, u);
			const float cosTheta = saturate(dot(N, dir));
			if (cosTheta <= 0)
				continue;

			RayDesc ray;
			ray.Origin = rayOrigin;
			ray.TMin = 0;
			ray.TMax = FLT_MAX;
			ray.Direction = dir;

			float3 sampleRadiance = 0;
			float3 samplePosition = rayOrigin + dir * 100000.0;
			float3 sampleNormal = -dir;

#ifdef RTAPI
			wiRayQuery q;
			q.TraceRayInline(
				scene_acceleration_structure,
				RAY_FLAG_SKIP_PROCEDURAL_PRIMITIVES | RAY_FLAG_FORCE_OPAQUE,
				push.instanceInclusionMask,
				ray);
			while (q.Proceed());
			const bool hit = q.CommittedStatus() == COMMITTED_TRIANGLE_HIT;
#else
			const RayHit rayHit =
				TraceRay_Closest(ray, push.instanceInclusionMask, rng);
			const bool hit = rayHit.distance < FLT_MAX - 1;
#endif // RTAPI

			[branch]
			if (!hit)
			{
				// Sky miss: the environment is the outgoing radiance.
				sampleRadiance = IsStaticSky()
					? GetStaticSkyColor(dir)
					: GetDynamicSkyColor(dir, true, false, true);
			}
			else
			{
				Surface surface;
				surface.init();

#ifdef RTAPI
				PrimitiveID prim;
				prim.init();
				prim.primitiveIndex = q.CommittedPrimitiveIndex();
				prim.instanceIndex = q.CommittedInstanceID();
				prim.subsetIndex = q.CommittedGeometryIndex();
				surface.SetBackface(!q.CommittedTriangleFrontFace());
				if (!surface.load(prim, q.CommittedTriangleBarycentrics()))
					continue;
				const float3 hitPos =
					q.WorldRayOrigin() + dir * q.CommittedRayT();
#else
				surface.SetBackface(rayHit.is_backface);
				if (!surface.load(rayHit.primitiveID, rayHit.bary))
					continue;
				const float3 hitPos = ray.Origin + dir * rayHit.distance;
#endif // RTAPI

				surface.P = hitPos;
				surface.V = -dir;
				surface.update();

				samplePosition = hitPos;
				sampleNormal = surface.N;
				sampleRadiance = RESTIRGIShadeHit(surface, rng);
			}

			// Clamp to the half-float storage range (R11G11B10) so a bright sky
			// or emitter cannot overflow to Inf and poison the reservoir.
			sampleRadiance = saturateMediump(sampleRadiance);

			// Cosine-weighted source pdf (solid angle): cosTheta / PI.
			const float sourcePdf = cosTheta / PI;
			const float targetPdf = RESTIRGITargetFunction(
				sampleRadiance, samplePosition, P, N);
			const float risWeight = (sourcePdf > 0) ? targetPdf / sourcePdf : 0;

			RESTIRGIReservoirUpdate(
				reservoir, samplePosition, sampleNormal, sampleRadiance,
				targetPdf, risWeight, rng);
		}
	}

	// The temporal pass merges history into this initial reservoir and resolves
	// the denoiser inputs (irradiance + gradient).
	RESTIRGIReservoirStore(reservoirOutput, flatIndex, reservoir);
}
