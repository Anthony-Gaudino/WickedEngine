#ifndef WI_RTSCENE_HF
#define WI_RTSCENE_HF

/**
 * Shared hardware ray-tracing scene-radiance gather.
 *
 * Provides `TraceSceneRadiance`, the reusable "trace a ray into the scene,
 * shade whatever it hits, or sample the sky on a miss" core. This is the
 * expensive, feature-agnostic part of any ray-traced lighting effect: the only
 * thing that differs between a reflection, a refraction, or an
 *ambient-occlusion trace is how the caller builds the ray and
 * what it does with the returned radiance.
 *
 * The gather here mirrors the closest-hit shading of the ray-traced reflection
 * pass (direct lights + environment reflection + DDGI + emissive, composited
 * with `ApplyLighting`), so consumers get a fully lit result for free.
 *
 * @note Hardware ray tracing only: the whole header is compiled out unless
 *       `RTAPI` is defined (the `_rtapi` shader permutation). Callers must
 *       provide a non-ray-traced fallback for devices without ray tracing.
 */

#include "globals.hlsli"
#include "raytracingHF.hlsli"
#include "lightingHF.hlsli"
#include "ShaderInterop_DDGI.h"

#ifdef RTAPI

/**
 * Result of a `TraceSceneRadiance` call.
 */
struct SceneRadiance
{
	/** Radiance gathered along the ray. Sky color when the ray missed. */
	float3 color;

	/** Hit distance along the ray, or `FLT_MAX` when the ray missed. */
	float distance;

	/** True when scene geometry was hit, false on a sky miss. */
	bool hit;
};

/**
 * Traces a ray into the scene and returns the radiance it gathers.
 *
 * Walks non-opaque candidate hits (alpha-tested with blue noise, accumulating
 * additive emissive as it goes), then either shades the closest opaque hit with
 * the full lighting model or, on a miss, samples the static/dynamic sky in the
 * ray direction.
 *
 * Example usage:
 * @code
 * RayDesc ray;
 * ray.Origin = surfacePoint;
 * ray.Direction = refractedDir;
 * ray.TMin = 0.01;
 * ray.TMax = FLT_MAX;
 * RayCone cone = RayCone::from_spread_angle(
 *     pixel_cone_spread_angle_from_image_height(resolution.y));
 * SceneRadiance rad = TraceSceneRadiance(
 *     ray, 0xFF,
 *     RAY_FLAG_SKIP_PROCEDURAL_PRIMITIVES | RAY_FLAG_CULL_BACK_FACING_TRIANGLES,
 *     cone, pixel, resolutionRcp);
 * float3 aboveWater = rad.color; // hit surface, or sky on miss
 * @endcode
 *
 * @param[in] ray - Ray to trace (origin, direction, TMin, TMax).
 * @param[in] instanceInclusionMask - TLAS instance inclusion mask.
 * @param[in] rayFlags - `RAY_FLAG_*` bits passed to the inline query.
 * @param[in] raycone - Ray cone used for mip selection at the hit.
 * @param[in] pixel - Dispatch pixel, used for blue-noise alpha test and the
 *                    hit surface's screen-space fields.
 * @param[in] resolutionRcp - Reciprocal render resolution for `screenUV`.
 *
 * @return The gathered radiance, hit distance, and hit/miss flag.
 *
 * @note The radiance is returned raw: no height fog, no aerial perspective
 *       and no water fog are applied over the traced distance. A consumer
 *       whose ray crosses any of those applies them itself, over `distance`
 *       from the ray origin - the camera's own path to that origin is a
 *       different column of air and cannot stand in for it.
 *
 * @note A miss returns the sky, which already carries the atmosphere and is
 *       not to be fogged again. It is only meaningful for a ray that can
 *       reach the sky at all; one travelling down into water tests `hit`.
 */
SceneRadiance TraceSceneRadiance(
	RayDesc ray,
	uint instanceInclusionMask,
	uint rayFlags,
	RayCone raycone,
	uint2 pixel,
	float2 resolutionRcp
)
{
	SceneRadiance result;
	result.color = 0;
	result.distance = FLT_MAX;
	result.hit = false;

	// Additive (unlit-additive) surfaces do not stop the ray; accumulate their
	// emissive and remember the nearest one so it is only added if it sits in
	// front of the committed opaque hit.
	float4 additive_dist = float4(0, 0, 0, FLT_MAX);

	wiRayQuery q;
	q.TraceRayInline(
		scene_acceleration_structure,	// RaytracingAccelerationStructure
		rayFlags,						// uint RayFlags
		instanceInclusionMask,			// uint InstanceInclusionMask
		ray								// RayDesc Ray
	);
	while (q.Proceed())
	{
		PrimitiveID prim;
		prim.init();
		prim.primitiveIndex = q.CandidatePrimitiveIndex();
		prim.instanceIndex = q.CandidateInstanceID();
		prim.subsetIndex = q.CandidateGeometryIndex();

		Surface surface;
		surface.init();
		surface.V = -ray.Direction;
#ifdef SURFACE_LOAD_MIPCONE
		surface.raycone = raycone;
#endif // SURFACE_LOAD_MIPCONE
		surface.hit_depth = q.CandidateTriangleRayT();
		if (!surface.load(prim, q.CandidateTriangleBarycentrics()))
			break;

		float alphatest =
			clamp(blue_noise(pixel, q.CandidateTriangleRayT()).r, 0, 0.99);

		if (surface.material.IsAdditive())
		{
			additive_dist.xyz += surface.emissiveColor;
			additive_dist.w = min(additive_dist.w, q.CandidateTriangleRayT());
		}
		else if (surface.opacity - alphatest >= 0)
		{
			q.CommitNonOpaqueTriangleHit();
		}
	}

	if (additive_dist.w <= q.CommittedRayT())
	{
		result.color += max(0, additive_dist.xyz);
	}

	if (q.CommittedStatus() != COMMITTED_TRIANGLE_HIT)
	{
		// Miss: sample the sky in the ray direction.
		[branch]
		if (IsStaticSky())
		{
			result.color += GetStaticSkyColor(q.WorldRayDirection());
		}
		else
		{
			result.color += GetDynamicSkyColor(q.WorldRayDirection());
		}
		result.distance = FLT_MAX;
		result.hit = false;
		return result;
	}

	// Closest hit: reconstruct and shade the surface.
	PrimitiveID prim;
	prim.init();
	prim.primitiveIndex = q.CommittedPrimitiveIndex();
	prim.instanceIndex = q.CommittedInstanceID();
	prim.subsetIndex = q.CommittedGeometryIndex();

	Surface surface;
	surface.init();
	surface.SetBackface(!q.CommittedTriangleFrontFace());
	surface.V = -ray.Direction;
#ifdef SURFACE_LOAD_MIPCONE
	surface.raycone = raycone;
#endif // SURFACE_LOAD_MIPCONE
	surface.hit_depth = q.CommittedRayT();
	if (!surface.load(prim, q.CommittedTriangleBarycentrics()))
	{
		// Degenerate hit: treat as a miss so the caller can fall back.
		return result;
	}

	surface.pixel = pixel;
	surface.screenUV = (surface.pixel + 0.5) * resolutionRcp;

	if (surface.material.IsUnlit())
	{
		result.color += surface.albedo + surface.emissiveColor;
	}
	else
	{
		surface.P = q.WorldRayOrigin() + q.WorldRayDirection() * q.CommittedRayT();
		surface.V = -q.WorldRayDirection();
		surface.update();

		if (!surface.IsGIApplied())
		{
			float3 ambient = GetAmbient(surface.N);
			surface.gi = lerp(ambient, ambient * surface.sss.rgb, saturate(surface.sss.a));
		}

		Lighting lighting;
		lighting.create(0, 0, surface.gi, 0);

		[loop]
		for (uint iterator = 0; iterator < lights().item_count(); iterator++)
		{
			ShaderEntity light = load_entity(lights().first_item() + iterator);
			if ((light.layerMask & surface.material.layerMask) == 0)
				continue;

			if (light.IsStaticLight())
			{
				continue; // static lights are baked into lightmaps
			}

			switch (light.GetType())
			{
			case ENTITY_TYPE_DIRECTIONALLIGHT:
			{
				light_directional(light, surface, lighting);
			}
			break;
			case ENTITY_TYPE_POINTLIGHT:
			{
				light_point(light, surface, lighting);
			}
			break;
			case ENTITY_TYPE_RECTLIGHT:
			{
				light_rect(light, surface, lighting);
			}
			break;
			case ENTITY_TYPE_SPOTLIGHT:
			{
				light_spot(light, surface, lighting);
			}
			break;
			}
		}

		lighting.indirect.specular += max(0, EnvironmentReflection_Global(surface));
		lighting.indirect.specular += surface.emissiveColor;

		[branch]
		if (GetScene().ddgi.probe_buffer >= 0)
		{
			lighting.indirect.diffuse = ddgi_sample_irradiance(
				surface.P, surface.N,
				surface.dominant_lightdir, surface.dominant_lightcolor
			);
		}

		float4 color = 0;
		ApplyLighting(surface, lighting, color);
		result.color += color.rgb;
	}

	result.distance = q.CommittedRayT();
	result.hit = true;
	return result;
}

#endif // RTAPI

#endif // WI_RTSCENE_HF
