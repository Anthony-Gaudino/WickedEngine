#ifndef WI_RESTIR_DI_VISIBILITY_HF
#define WI_RESTIR_DI_VISIBILITY_HF

/**
 * Traces a shadow ray from a shading point toward a resolved light sample.
 *
 * Shared by the ReSTIR DI initial pass and the spatial pass's final visibility
 * re-trace. The ray:
 *   - offsets its origin along the normal (scaled with camera distance) to
       avoid
 *     grazing self-shadow without leaking on small near geometry,
 *   - uses the shadow inclusion mask (bit 0) so water and cast-shadow-disabled
 *     objects do not occlude (mirrors
       wi::renderer::raytracing_inclusion_mask_shadow),
 *   - alpha-tests non-opaque hits so cutout foliage casts a shaped shadow.
 *
 * Requires "globals.hlsli", "raytracingHF.hlsli" and "lightingHF.hlsli"
   included
 * first, and TEXTURE_SLOT_NONUNIFORM defined (for the divergent material load).
 *
 * @param[in] P - World-space shading point.
 * @param[in] N - World-space shading normal.
 * @param[in] samplePosition - World-space point on the chosen light.
 * @param[in,out] rng - Random generator (stochastic alpha test / software BVH).
 *
 * @return 1 if the light is visible, 0 if occluded.
 */
inline float RESTIRDITraceVisibility(
	float3 P, float3 N, float3 samplePosition, inout RNG rng)
{
	const float distToCamera = length(P - GetCamera().position);
	const float shadowBias = clamp(0.002 * distToCamera, 0.02, 0.5);
	const float3 rayOrigin = P + N * shadowBias;
	const float3 d = samplePosition - rayOrigin;
	const float dist = length(d);
	const float3 dir = d / max(dist, 1e-6);

	RayDesc ray;
	ray.Origin = rayOrigin;
	ray.TMin = 0.01;
	ray.TMax = max(0.02, dist - 0.02);
	ray.Direction = dir;

	const uint shadow_ray_mask = 1u;

	float shadow = 1;
#ifdef RTAPI
	wiRayQuery q;
	q.TraceRayInline(
		scene_acceleration_structure,	// RaytracingAccelerationStructure AccelerationStructure
		// No RAY_FLAG_FORCE_OPAQUE: non-opaque hits are alpha-tested below.
		RAY_FLAG_SKIP_PROCEDURAL_PRIMITIVES |
		RAY_FLAG_ACCEPT_FIRST_HIT_AND_END_SEARCH,	// uint RayFlags
		shadow_ray_mask,				// uint InstanceInclusionMask
		ray								// RayDesc Ray
	);
	while (q.Proceed())
	{
		// Non-opaque candidate: occlude only where the (stochastic) alpha test
		// passes; additive materials do not cast shadows.
		PrimitiveID prim;
		prim.init();
		prim.primitiveIndex = q.CandidatePrimitiveIndex();
		prim.instanceIndex = q.CandidateInstanceID();
		prim.subsetIndex = q.CandidateGeometryIndex();

		Surface hitSurface;
		hitSurface.init();
		hitSurface.V = ray.Direction;
		[branch]
		if (hitSurface.load(prim, q.CandidateTriangleBarycentrics()))
		{
			if (!hitSurface.material.IsAdditive() &&
				hitSurface.opacity - rng.next_float() >= 0)
			{
				q.CommitNonOpaqueTriangleHit();
			}
		}
	}
	shadow = q.CommittedStatus() == COMMITTED_TRIANGLE_HIT ? 0 : shadow;
#else
	shadow = TraceRay_Any(ray, shadow_ray_mask, rng) ? 0 : shadow;
#endif // RTAPI

	return shadow;
}

#endif // WI_RESTIR_DI_VISIBILITY_HF
