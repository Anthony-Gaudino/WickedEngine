#define TEXTURE_SLOT_NONUNIFORM
#include "globals.hlsli"
#include "raytracingHF.hlsli"
#include "lightingHF.hlsli"
#include "restir_lightsamplingHF.hlsli"
#include "ShaderInterop_DDGI.h"

// This shader runs one probe per thread group and each thread will trace rays and write the trace result to a ray data buffer
//	ray data buffer will be later integrated by ddgi_updateCS shader which updates the DDGI irradiance and depth textures

PUSHCONSTANT(push, DDGIPushConstants);

// ReSTIR RIS candidate count per probe ray. One candidate reproduces the old
// uniform single-sample estimator exactly; more candidates importance-sample
// the lights (still one shadow ray), cutting noise when many lights are
// present.
static const uint RESTIR_DDGI_CANDIDATES = 8;

StructuredBuffer<uint> rayallocationBuffer : register(t0);
Buffer<uint> raycountBuffer : register(t1);

RWStructuredBuffer<DDGIRayDataPacked> rayBuffer : register(u0);

groupshared float shared_inconsistency[DDGI_COLOR_RESOLUTION * DDGI_COLOR_RESOLUTION];
groupshared uint shared_rayCount;

static const uint THREADCOUNT = 32;

// spherical fibonacci: https://github.com/diharaw/hybrid-rendering/blob/master/src/shaders/gi/gi_ray_trace.rgen
#define madfrac(A, B) ((A) * (B)-floor((A) * (B)))
static const float PHI = sqrt(5) * 0.5 + 0.5;
float3 spherical_fibonacci(float i, float n)
{
	float phi = 2.0 * PI * madfrac(i, PHI - 1);
	float cos_theta = 1.0 - (2.0 * i + 1.0) * (1.0 / n);
	float sin_theta = sqrt(clamp(1.0 - cos_theta * cos_theta, 0.0f, 1.0f));
	return float3(cos(phi) * sin_theta, sin(phi) * sin_theta, cos_theta);
}

[numthreads(THREADCOUNT, 1, 1)]
void main(uint3 DTid : SV_DispatchThreadID, uint3 Gid : SV_GroupID, uint groupIndex : SV_GroupIndex)
{
	const uint allocCount = rayallocationBuffer[3];
	if(DTid.x >= allocCount)
		return;
	const uint rayAlloc = rayallocationBuffer[4 + DTid.x];
	const uint probeIndex = rayAlloc & 0xFFFFF;
	const uint rayIndex = rayAlloc >> 20u;
	const uint rayCount = raycountBuffer[probeIndex] * DDGI_RAY_BUCKET_COUNT;
	const uint3 probeCoord = ddgi_probe_coord(probeIndex);
	const float3 probePos = ddgi_probe_position(probeCoord);

	RNG rng;
	rng.init(DTid.xx, GetFrame().frame_count);
	
	float3 radiance = 0;

	const float3x3 random_orientation = (float3x3)g_xTransform;
	const float3 raydir = normalize(mul(random_orientation, spherical_fibonacci(rayIndex, rayCount)));

#if 1
	// Light sampling - direct static (ReSTIR RIS over the analytic light list):
	// resample RESTIR_DDGI_CANDIDATES lights by the diffuse target function and
	// trace one shadow ray for the winner. staticOnly keeps the probe cache
	// limited to static lights (dynamic lights are applied in real time
	// elsewhere). With one candidate this matches the previous single-uniform-
	// sample estimator exactly. The probe ray direction acts as the shading
	// normal (probes are omnidirectional).
	{
		RESTIRLightSample winning;
		const RESTIRReservoir reservoir = RESTIRSampleLightsUniform(
			probePos, raydir, RESTIR_DDGI_CANDIDATES, true, rng, winning);

		const float W = RESTIRReservoirGetInvPdf(reservoir);
		const float NdotL = saturate(dot(winning.direction, raydir));

		[branch]
		if (W > 0 && NdotL > 0 && winning.distance > 0)
		{
			float3 shadow = 1;

			RayDesc newRay;
			newRay.Origin = probePos;
			newRay.TMin = 0;
			newRay.TMax = winning.distance;
			newRay.Direction = winning.direction;

	#ifdef RTAPI
			wiRayQuery q;
			q.TraceRayInline(
				scene_acceleration_structure,	// RaytracingAccelerationStructure AccelerationStructure
				//RAY_FLAG_CULL_FRONT_FACING_TRIANGLES |
				RAY_FLAG_SKIP_PROCEDURAL_PRIMITIVES |
				RAY_FLAG_ACCEPT_FIRST_HIT_AND_END_SEARCH,	// uint RayFlags
				0xFF,							// uint InstanceInclusionMask
				newRay							// RayDesc Ray
			);
			while (q.Proceed());
			shadow = q.CommittedStatus() == COMMITTED_TRIANGLE_HIT ? 0 : shadow;
	#else
			shadow = TraceRay_Any(newRay, push.instanceInclusionMask, rng) ? 0 : shadow;
	#endif // RTAPI
			if (any(shadow))
			{
				radiance += max(0, shadow * winning.radiance * NdotL * W / PI);
			}
		}
	}
#endif
	
	{
		RayDesc ray;
		ray.Origin = probePos;
		ray.TMin = 0; // don't need TMin because we are not tracing from a surface
		ray.TMax = FLT_MAX;
		ray.Direction = normalize(raydir);

#ifdef RTAPI
		wiRayQuery q;
		q.TraceRayInline(
			scene_acceleration_structure,	// RaytracingAccelerationStructure AccelerationStructure
			//RAY_FLAG_CULL_BACK_FACING_TRIANGLES |
			RAY_FLAG_SKIP_PROCEDURAL_PRIMITIVES |
			RAY_FLAG_FORCE_OPAQUE,			// uint RayFlags
			push.instanceInclusionMask,		// uint InstanceInclusionMask
			ray								// RayDesc Ray
		);
		while (q.Proceed());
		if (q.CommittedStatus() != COMMITTED_TRIANGLE_HIT)
#else
		RayHit hit = TraceRay_Closest(ray, push.instanceInclusionMask, rng);

		if (hit.distance >= FLT_MAX - 1)
#endif // RTAPI

		{
			float3 envColor;
			[branch]
			if (IsStaticSky())
			{
				// We have envmap information in a texture:
				envColor = GetStaticSkyColor(ray.Direction);
			}
			else
			{
				envColor = GetDynamicSkyColor(ray.Direction, true, false, true);
			}
			radiance += envColor;

			// Clamp to the half-float range before storing: rayData.radiance is
			// half3, and an unclamped bright sky (e.g. the sun disk) overflows
			// to +Inf, which the firefly suppression later turns into NaN (Inf
			// - Inf), poisoning the probe and blowing the whole frame to white.
			// The surface-hit path below already does this.
			radiance = saturateMediump(radiance);

			DDGIRayData rayData;
			rayData.direction = ray.Direction;
			rayData.depth = -1;
			rayData.radiance = radiance;
			rayBuffer[probeIndex * DDGI_MAX_RAYCOUNT + rayIndex].store(rayData);
		}
		else
		{

			Surface surface;
			surface.init();

			float hit_depth = 0;
			float3 hit_result = 0;

#ifdef RTAPI

			// ray origin updated for next bounce:
			ray.Origin = q.WorldRayOrigin() + q.WorldRayDirection() * q.CommittedRayT();
			hit_depth = q.CommittedRayT();

			PrimitiveID prim;
			prim.init();
			prim.primitiveIndex = q.CommittedPrimitiveIndex();
			prim.instanceIndex = q.CommittedInstanceID();
			prim.subsetIndex = q.CommittedGeometryIndex();

			surface.SetBackface(!q.CommittedTriangleFrontFace());

			if (!surface.load(prim, q.CommittedTriangleBarycentrics()))
				return;

#else

			// ray origin updated for next bounce:
			ray.Origin = ray.Origin + ray.Direction * hit.distance;
			hit_depth = hit.distance;

			surface.SetBackface(hit.is_backface);

			if (!surface.load(hit.primitiveID, hit.bary))
				return;

#endif // RTAPI

			if (surface.IsBackface())
			{
				hit_depth *= 0.9; // push inwards to help avoid shadow leaks from inwards to outside
			}

			surface.P = ray.Origin;
			surface.V = -ray.Direction;
			surface.update();

#if 1
			// Light sampling (ReSTIR RIS over the analytic light list):
			// resample RESTIR_DDGI_CANDIDATES lights by the diffuse target
			// function and trace one shadow ray for the winner. All lights
			// contribute here (this is the indirect bounce hit, not the
			// static-only cache block).
			{
				RESTIRLightSample winning;
				const RESTIRReservoir reservoir = RESTIRSampleLightsUniform(
					surface.P, surface.N, RESTIR_DDGI_CANDIDATES, false, rng, winning);

				const float W = RESTIRReservoirGetInvPdf(reservoir);
				const float NdotL = saturate(dot(winning.direction, surface.N));

				[branch]
				if (W > 0 && NdotL > 0 && winning.distance > 0)
				{
					float3 shadow = 1;

					RayDesc newRay;
					newRay.Origin = surface.P;
					newRay.TMin = 0.001;
					newRay.TMax = winning.distance;
					newRay.Direction = normalize(winning.direction + max3(surface.sss));

#ifdef RTAPI
					q.TraceRayInline(
						scene_acceleration_structure,	// RaytracingAccelerationStructure AccelerationStructure
						//RAY_FLAG_CULL_FRONT_FACING_TRIANGLES |
						RAY_FLAG_SKIP_PROCEDURAL_PRIMITIVES |
						RAY_FLAG_ACCEPT_FIRST_HIT_AND_END_SEARCH,	// uint RayFlags
						0xFF,							// uint InstanceInclusionMask
						newRay							// RayDesc Ray
					);
					while (q.Proceed());
					shadow = q.CommittedStatus() == COMMITTED_TRIANGLE_HIT ? 0 : shadow;
#else
					shadow = TraceRay_Any(newRay, push.instanceInclusionMask, rng) ? 0 : shadow;
#endif // RTAPI
					if (any(shadow))
					{
						hit_result += max(0, shadow * winning.radiance * NdotL * W / PI);
					}
				}
			}
#endif

			// Infinite bounces based on previous frame probe sampling:
			if (push.frameIndex > 0)
			{
				float energy_conservation = 0.95;
				energy_conservation /= PI; // one more divide by PI is inside the ddgi_sample_irradiance, with that we will have 2 PI divides, which is needed for hemishpere sampling
				float3 ddgi = ddgi_sample_irradiance(surface.P, surface.facenormal, surface.dominant_lightdir, surface.dominant_lightcolor);
				ddgi *= energy_conservation;
				hit_result += ddgi;
			}
			hit_result *= surface.albedo;
			hit_result += surface.emissiveColor;
			radiance += hit_result;

			radiance = saturateMediump(radiance);

			DDGIRayData rayData;
			rayData.direction = ray.Direction;
			rayData.depth = hit_depth;
			rayData.radiance = radiance;
			rayBuffer[probeIndex * DDGI_MAX_RAYCOUNT + rayIndex].store(rayData);
		}

	}
}
