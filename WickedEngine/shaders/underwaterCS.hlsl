#include "globals.hlsli"
#include "ShaderInterop_Postprocess.h"
#include "oceanSurfaceHF.hlsli"
#include "lightingHF.hlsli"
#include "fogHF.hlsli"
#include "waterFogHF.hlsli"
#include "volumetricFroxelHF.hlsli"
#ifdef RTAPI
	// Hardware ray-traced Snell's window: trace the refracted ray into the real
	// scene so above-water objects, and their occlusion of the window, are
	// exact. Compiled only in the underwaterCS_rtapi permutation; the base
	// permutation shows the analytic sky in the refracted direction instead.
	#include "rtsceneHF.hlsli"
#endif // RTAPI

#define INTERSECTION_DISTORT
//#define ANIMATED_DISTORT

PUSHCONSTANT(postprocess, PostProcess);

Texture2D<float4> input : register(t0);

RWTexture2D<float4> output : register(u0);


// https://www.shadertoy.com/view/MlSXR3
float2 brownConradyDistortion(float2 uv)
{
	// positive values of K1 give barrel distortion, negative give pincushion
	//float barrelDistortion1 = 0.15; // K1 in text books
	//float barrelDistortion2 = 0.0; // K2 in text books
	float barrelDistortion1 = 0.25; // K1 in text books
	float barrelDistortion2 = -0.34; // K2 in text books
	float r2 = uv.x*uv.x + uv.y*uv.y;
	uv *= 1.0 + barrelDistortion1 * r2 + barrelDistortion2 * r2 * r2;
	
	// tangential distortion (due to off center lens elements)
	// is not modeled in this function, but if it was, the terms would go here
	return uv;
}

[numthreads(POSTPROCESS_BLOCKSIZE, POSTPROCESS_BLOCKSIZE, 1)]
void main(uint3 DTid : SV_DispatchThreadID)
{
	float2 uv = (DTid.xy + 0.5f) * postprocess.resolution_rcp;

	// Unproject near plane and determine for every pixel if it's below water surface:
	float4 clipspace = float4(uv_to_clipspace(uv), ocean_waterline_handoff_depth(), 1); // push further away from near plane
	float4 unproj = mul(GetCamera().inverse_view_projection, clipspace);
	unproj.xyz /= unproj.w;
	float3 world_pos = unproj.xyz;
	const ShaderOcean ocean = GetWeather().ocean;
	float3 ocean_pos = float3(world_pos.x, ocean.water_height, world_pos.z);
	[branch]
	if (ocean.texture_displacementmap >= 0)
	{
		const float2 ocean_uv = ocean_pos.xz * ocean.patch_size_rcp;
		Texture2D texture_displacementmap = bindless_textures[descriptor_index(ocean.texture_displacementmap)];
		const float3 displacement = texture_displacementmap.SampleLevel(sampler_linear_wrap, ocean_uv, 0).xzy;
		ocean_pos += displacement;
	}

#ifdef INTERSECTION_DISTORT
	// Distort at intersection:
	float intersection_direction = world_pos.y - ocean_pos.y;
	float intersection_distance = abs(intersection_direction);
	float intersection_blend = saturate(exp(-intersection_distance * 50));
	clipspace.xy *= lerp(1, 0.5, intersection_blend);
	uv = clipspace_to_uv(clipspace.xy);

	// Recompute ray after distortion:
	unproj = mul(GetCamera().inverse_view_projection, clipspace);
	unproj.xyz /= unproj.w;
	world_pos = unproj.xyz;
	ocean_pos = float3(world_pos.x, ocean.water_height, world_pos.z);
	[branch]
	if (ocean.texture_displacementmap >= 0)
	{
		const float2 ocean_uv = ocean_pos.xz * ocean.patch_size_rcp;
		Texture2D texture_displacementmap = bindless_textures[descriptor_index(ocean.texture_displacementmap)];
		const float3 displacement = texture_displacementmap.SampleLevel(sampler_linear_wrap, ocean_uv, 0).xzy;
		ocean_pos += displacement;
	}
#endif // INTERSECTION_DISTORT

	float4 original_color = input.SampleLevel(sampler_linear_clamp, uv, 0);
	float4 color = original_color;

	// Apply effects on full screen:
	{

		// Some lens distortion uv modulation, as if seen through a dive mask.
		// Scaled rather than switched, so the setting can dial it back as well
		// as turn it off, and so 1 reproduces the engine's long-standing look
		// exactly.
		[branch]
		if (underwater_lens_distortion > 0)
		{
			const float2 distorted =
				clipspace_to_uv(brownConradyDistortion(uv_to_clipspace(uv) * 0.9));
			uv = lerp(uv, distorted, saturate(underwater_lens_distortion));
		}

#ifdef ANIMATED_DISTORT
		// It's not realistic to apply much refraction underwater to camera, but looks cool:
		uv += sin(uv * 10 + GetFrame().time * 5) * 0.005;
		uv += sin(-uv.y * 5 + GetFrame().time * 2) * 0.005;
#endif // ANIMATED_DISTORT

		// Underwater magnification: refraction at the eye/mask interface makes
		// submerged objects appear larger (and nearer). Model it as a radial
		// zoom toward the optical axis (screen center). A magnification of 1 is
		// a no-op. The final waterline blend below keeps the above-water part
		// of a partially submerged view at its true scale.
		uv = 0.5 + (uv - 0.5) / underwater_magnification;

		const float depth = texture_depth.SampleLevel(sampler_linear_clamp, uv, 0);
		float3 surface_position = reconstruct_position(uv, depth);

		float4 clipspace2 = float4(uv_to_clipspace(uv), 1, 1); // exact near plane
		float4 unproj2 = mul(GetCamera().inverse_view_projection, clipspace2);
		unproj2.xyz /= unproj2.w;
		float3 campos = unproj2.xyz;

		float3 V = campos - surface_position.xyz;
		float surface_dist = length(V);
		V /= surface_dist;

		// View ray from the eye into the scene:
		const float3 rayDir = -V;

		// The ocean is not rendered into the lineardepth unfortunately, so we
		// also trace it: Otherwise the ocean surface could be same as infinite
		// depth and incorrectly fogged
		//
		// Every quantity the window is built from is meaningless without a
		// crossing, so where the ray finds none the window is switched off
		// rather than handed the nearest wrong answer.
		float ocean_dist = 0;
		bool ocean_surface_ahead = false;
		float3 ocean_surface_pos = campos;
		float3 ocean_surface_normal = float3(0, 1, 0);

		[branch]
		if (GetCamera().IsWaterSegmentModel())
		{
			// Walked along the ray, which is the only thing that finds water
			// standing AHEAD of the eye rather than over it - the dry half of a
			// camera straddling a crest, looking through the wave it is sitting
			// in.
			const WaterSegment water = TraceWaterSegment(
				campos,
				mad(GetCamera().z_far, rayDir, campos),
				blue_noise(DTid.xy).x);

			ocean_surface_ahead = water.crosses;
			ocean_dist = water.entry;
			ocean_surface_pos = water.entryPoint;
			ocean_surface_normal = water.entryNormal;
		}
		else
		{
			// Two steps of a vertical height field solve: the first lands the
			// hit at the right horizontal position, the second corrects the
			// height there. The refinement is kept only while it still lands in
			// front of the eye, because a grazing ray re-intersected against a
			// height field walks away rather than converging.
			const float3 ocean_up = float3(0, 1, 0);
			float solved_dist = intersectPlaneClampInfiniteDist(
				campos, rayDir, ocean_up, ocean_drawn_surface_height(campos));
			[branch]
			if (solved_dist > 0)
			{
				const float refined_dist = intersectPlaneClampInfiniteDist(
					campos, rayDir, ocean_up,
					ocean_drawn_surface_height(campos + rayDir * solved_dist));
				if (refined_dist > 0)
				{
					solved_dist = refined_dist;
				}
			}

			ocean_surface_ahead = solved_dist > 0;
			ocean_dist = max(solved_dist, 0);
			ocean_surface_pos = campos + rayDir * ocean_dist;
			ocean_surface_normal =
				WaterSegmentSurfaceNormal(ocean_surface_pos);
		}

		// The medium, for Snell's window below. This pass no longer fogs
		// anything: every fragment applied the water's fog over its own view
		// path as it was drawn, which is the only way a transparent can be
		// fogged by its own distance rather than by whatever is behind it.
		const WaterVolumetrics medium = MakeWaterVolumetrics(1);

		color = input.SampleLevel(sampler_linear_mirror, uv, 0);

		// Snell's window: looking up from under water, refraction gathers the
		// whole above-water hemisphere into an overhead circular window of
		// fixed angular size - half-angle equal to the critical angle for water
		// (~48.6 deg, so ~97.2 deg across), independent of depth. What changes
		// with depth is how much of it survives the water column: the light
		// that forms the window is attenuated along its path to the surface, so
		// the rim (a longer path than the center) fades first and the bright
		// disc appears to shrink, and as the camera descends the crisp
		// above-water view washes out into a glowing patch of sunlight that
		// finally fades to nothing once the water swallows the light (the
		// light-reach depth). Only drawn for upward rays that reach the surface
		// through open water, so it never covers geometry and never appears
		// when looking down. References:
		// https://en.wikipedia.org/wiki/Snell%27s_window
		if (underwater_snell > 0)
		{
			// Surface normal where the ray pierces the water, taken from the
			// same wave gradient map the ocean surface is shaded with. Snell's
			// law is relative to the local surface normal, so testing the
			// critical angle against this wavy normal (instead of a flat up
			// vector) makes the window rim genuinely ondulate with the waves -
			// physically based, no procedural ripple:
			const float3 surfaceNormal = ocean_surface_normal;

			const float cosTheta = dot(rayDir, surfaceNormal);
			const float sinTheta = sqrt(saturate(1.0 - cosTheta * cosTheta));

			// Snell's law: sin(theta) * n_water reaches 1 at the critical
			// angle, so the window edge sits right at that angle at any depth:
			const float escape = sinTheta * 1.333;
			const float cone = 1.0 - smoothstep(0.92, 1.03, escape);

			// Only upward rays looking at open water (the surface nearer than
			// any geometry) form the window: never looking down, never over
			// objects. Gate on true world up so it stays robust to the wavy
			// normal above:
			const float upward = smoothstep(0.0, 0.1, rayDir.y);
			// The window only applies where the ocean surface is the nearest
			// thing along the ray (so geometry in front is never covered). Use
			// an ABSOLUTE gap, not a ratio: an object a fixed height above the
			// water must stay fully inside the window at every depth. A ratio
			// shrinks with depth and peels the object open as a dark silhouette
			// from the bottom up:
			const float openWater = smoothstep(0.0, 0.5, surface_dist - ocean_dist);

			// The window's light has to cross the whole water column from the
			// surface down to the eye, so the medium attenuates it exactly like
			// anything else - Beer-Lambert over ocean_dist with the real
			// per-channel extinction. The rim's path is longer than the
			// center's, so it fades first and the disc appears to shrink as the
			// camera descends; and because the extinction is spectral the
			// window also SHIFTS colour as it goes, losing red long before
			// blue, instead of just dimming. Turbidity closes it much sooner
			// than clear water, with no separate depth to author. The fade
			// param scales the medium so an artist can keep the window alive
			// deeper (<1) or shut it earlier (>1); 1 is physical:
			const float3 pathTransmittance =
				exp(-ocean_dist * medium.sigmaT * underwater_snell_fade);
			// Scalar stand-in for the activity gate below, taken from the
			// channel that survives furthest:
			const float pathFade = max(pathTransmittance.r,
				max(pathTransmittance.g, pathTransmittance.b));

			// Inside the window, refract the ray through the wave surface into
			// the air (water -> air) and read the above-water world in that
			// direction. Refraction squeezes the whole 180 deg hemisphere into
			// the ~97 deg cone, with the sun and objects landing at their true
			// refracted places; refract() returns 0 past the critical angle,
			// matching the cone edge:
			const float3 refractedDir = refract(rayDir, -surfaceNormal, 1.333);

			// Only rays that actually form the window need above-water content:
			// past the critical angle refract() returns 0 (no real refraction),
			// and once pathFade has decayed the window has washed out to water
			// color regardless of what is above. This also skips the ray-traced
			// window entirely at the light-reach depth, saving that cost.
			const bool windowActive =
				ocean_surface_ahead
				&& (cone * upward * openWater) > 0.0
				&& any(refractedDir)
				&& pathFade > 0.003;

			// The above-water view seen through the window. With ray tracing
			// the refracted ray is traced into the real scene, so above-water
			// objects appear and correctly occlude the window. The first hit is
			// that true surface and a miss returns the sky. Without ray tracing
			// (no capability, or the toggle off) it falls back to the analytic
			// sky. This is the shared gather that also drives RT reflections;
			// see rtsceneHF.hlsli.
			float3 aboveWater;
#ifdef RTAPI
			[branch]
			if (underwater_snell_rt != 0 && windowActive && GetScene().TLAS >= 0)
			{
				RayDesc ray;
				ray.Origin = ocean_surface_pos;
				ray.Direction = refractedDir;
				ray.TMin = 0.01;
				ray.TMax = FLT_MAX;

				RayCone raycone = RayCone::from_spread_angle(
					pixel_cone_spread_angle_from_image_height(
						postprocess.resolution.y));
				raycone = raycone.propagate(0, ocean_dist);

				const SceneRadiance rad = TraceSceneRadiance(
					ray,
					~0u,
					RAY_FLAG_SKIP_PROCEDURAL_PRIMITIVES |
					RAY_FLAG_CULL_BACK_FACING_TRIANGLES,
					raycone,
					DTid.xy,
					postprocess.resolution_rcp
				);
				aboveWater = rad.color;

				// The gather gives the radiance leaving the hit and nothing
				// of the air in front of it, so the window has to fog its own
				// content or a distant shore reads through it as crisply as
				// one at arm's length. Measured from the crossing along the
				// refracted ray, which is the air the light actually crossed;
				// the camera's own path to here is underwater and would be
				// the wrong column. A miss returns sky, which arrives with
				// the atmosphere already in it.
				[branch]
				if (rad.hit)
				{
					const half4 airFog = GetFog(
						rad.distance, ocean_surface_pos, refractedDir);
					aboveWater = lerp(
						aboveWater, (float3)airFog.rgb, (float)airFog.a);
				}
			}
			else
#endif // RTAPI
			{
				aboveWater = GetDynamicSkyColor(refractedDir);
			}

			// The window shows the refracted above-water view, extinguished by
			// the water it crossed. Per channel, so it loses its warm end first
			// and goes blue before it goes dark.
			//
			// It fades to nothing rather than to the fog's inscatter, which
			// carries the sun and the god-ray modulation and would draw a
			// god-ray texture inside the window as it hazed out. What replaces
			// it is the surrounding underwater view, through the blend below.
			// The window shows what lies beyond the surface, so it carries the
			// light scattered on THIS side of it - the water between the eye
			// and the crossing, which the window does not replace and which is
			// exactly where a submerged lamp near the surface puts its beams.
			//
			// Taken at the crossing itself, which this pass has already had to
			// find. That is the whole column the window covers and no more,
			// where a screen-space buffer could only offer the light gathered
			// out to whatever opaque surface stood behind the window.
			half4 windowContent = half4(aboveWater * pathTransmittance, 1);
			ApplyVolumetricLight(uv, ocean_surface_pos, windowContent);

			const float windowShape = ocean_surface_ahead
				? saturate(cone * upward * openWater)
				: 0;
			const float windowBlend = saturate(windowShape * underwater_snell);
			color.rgb = lerp(color.rgb, (float3)windowContent.rgb, windowBlend);
		}

		//color = float4(1, 0, 0, 1);
	}

	// Constrain the effect smoothly just below the water line:
	color = lerp(color, original_color, smoothstep(0.0, 0.025, saturate(world_pos.y - ocean_pos.y - 0.01)));
	//color = smoothstep(0.0, 0.05, saturate(intersection_direction));
	
	output[DTid.xy] = color;
}
