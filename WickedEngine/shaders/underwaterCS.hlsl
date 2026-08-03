#include "globals.hlsli"
#include "ShaderInterop_Postprocess.h"
#include "oceanSurfaceHF.hlsli"
#include "lightingHF.hlsli"
#include "fogHF.hlsli"
#include "waterFogHF.hlsli"
#ifdef RTAPI
	// Hardware ray-traced Snell's window: trace the refracted ray into the real
	// scene so above-water objects (and their occlusion of the window) are
	// exact, with no single-point probe to sit inside geometry. Compiled only
	// in the underwaterCS_rtapi permutation; the base permutation falls back to
	// the scene probe / analytic sky.
	#include "rtsceneHF.hlsli"
#endif // RTAPI

#define INTERSECTION_DISTORT
#define LENS_DISTORT
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
	float4 clipspace = float4(uv_to_clipspace(uv), OCEAN_NEARPLANE_CUTOFF, 1); // push further away from near plane
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
	//color.rgb = lerp(color.rgb, ocean.water_color.rgb, intersection_blend);

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

#ifdef LENS_DISTORT
		// Some lens distortion uv modulation:
		uv = clipspace_to_uv(brownConradyDistortion(uv_to_clipspace(uv) * 0.9));
#endif // LENS_DISTORT

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

		// The ocean is not rendered into the lineardepth unfortunately, so we also trace it:
		//	Otherwise the ocean surface could be same as infinite depth and incorrectly fogged
		float3 ocean_surface_pos = intersectPlaneClampInfinite(campos, -V, float3(0, 1, 0), ocean.water_height);
		float2 ocean_surface_uv = ocean_surface_pos.xz * ocean.patch_size_rcp;
		Texture2D texture_displacementmap = bindless_textures[descriptor_index(ocean.texture_displacementmap)];
		const float3 displacement = texture_displacementmap.SampleLevel(sampler_linear_wrap, ocean_surface_uv, 0).xzy;
		ocean_surface_pos += displacement;
		const float ocean_dist = length(ocean_surface_pos - campos);

		// How deep the EYE sits below the surface. Not a path: it is the water
		// column the sun had to come down through before it could scatter
		// towards the eye at all.
		const float camera_depth = max(0, ocean_pos.y - world_pos.y);

		// View ray from the eye into the scene:
		const float3 rayDir = -V;

		// Water the light crossed on its way to the eye: the stretch of the
		// segment from this pixel's near-plane point to whatever it sees that
		// lies below the surface.
		//
		// One expression for every case the reconstruction here used to sort
		// out by hand. Both ends under water gives the whole distance; the eye
		// under and the far end above - a shore, a boat, the sky - gives the
		// stretch up to the crossing; the eye above and the far end under gives
		// the stretch down from it. The far plane needs no special case: there
		// is no geometry there, the reconstruction simply lands a very long way
		// off, and a downward ray correctly takes all of it.
		//
		// That replaces an exit-distance trace, a smoothstep over the shaded
		// point's own depth, and a lerp between the two - machinery that existed
		// only because a screen-space pass had to guess which of those three
		// cases it was looking at from a depth buffer the ocean is not in.
		//
		// Clipped against the STILL water plane, not the displaced surface: the
		// crossing point's xz is not known until the crossing is found, so a
		// displaced clip would need a march. It is also the plane
		// SubmergedLightPath clips the light's own path against, so the fog and
		// the lighting it fogs agree about where the water stops.
		const float submerged_fraction = saturate(
			(ocean.water_height - min(campos.y, surface_position.y))
			/ max(abs(surface_position.y - campos.y), 0.00001));
		const float fog_distance = surface_dist * submerged_fraction;

		// The water's own fog over the view ray's path through it. Shared with
		// the draw shaders, so a pixel fogged here and a pixel fogged where it
		// was drawn cannot drift apart.
		//
		// The medium is taken unfaded: this whole block is already masked to
		// the submerged part of the screen by the waterline blend at the end.
		const WaterVolumetrics medium = MakeWaterVolumetrics(1);
		const WaterFog fog = MakeWaterFog(
			medium,
			fog_distance,
			V,
			camera_depth,
			uv,
			clipspace2.xy,
			underwater_godrays_procedural != 0
		);

		const float3 transmittance = fog.transmittance;

		color = input.SampleLevel(sampler_linear_mirror, uv, 0);

		// Single-scattering composite: what survives the water plus what the
		// water scatters back in.
		color.rgb = color.rgb * transmittance + fog.inscatter;

		// Fraction of the volumetric light already sitting in the input that
		// survives this pass.
		//
		// Volumetric lights are composited into the scene back in
		// RenderTransparents, long before this runs, so everything below treats
		// them as if they were surface radiance emitted at the far end of the
		// water column - attenuating them all over again, when the volumetric
		// march has already applied its own transmittance along the same ray.
		// Rather than reorder the frame (the composite is a blended draw inside
		// the transparent render pass, and the post chain has no render pass and
		// targets without RTVs), track what each stage takes away and give the
		// difference back at the end. That is exact, needs only additions, and
		// leaves the above-water path untouched.
		float3 volSurvival = transmittance;

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
			// Surface normal where the ray pierces the water, rebuilt from the
			// same wave gradient map the ocean surface is shaded with. Snell's
			// law is relative to the local surface normal, so testing the
			// critical angle against this wavy normal (instead of a flat up
			// vector) makes the window rim genuinely ondulate with the waves -
			// physically based, no procedural ripple:
			float3 surfaceNormal = float3(0, 1, 0);
			if (ocean.texture_gradientmap >= 0)
			{
				Texture2D texture_gradientmap =
					bindless_textures[descriptor_index(ocean.texture_gradientmap)];
				const float2 gradient = texture_gradientmap.SampleLevel(
					sampler_linear_wrap, ocean_surface_uv, 0
				).rg;
				surfaceNormal = normalize(
					float3(gradient.x, ocean.texel_length * 2.0, gradient.y)
				);
			}

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
			// matching the cone edge. When a realtime scene probe is available
			// it carries the real above-water geometry; otherwise use analytic
			// sky:
			const float3 refractedDir = refract(rayDir, -surfaceNormal, 1.333);

			// Only rays that actually form the window need above-water content:
			// past the critical angle refract() returns 0 (no real refraction),
			// and once pathFade has decayed the window has washed out to water
			// color regardless of what is above. This also skips the ray-traced
			// window entirely at the light-reach depth, saving that cost.
			const bool windowActive =
				(cone * upward * openWater) > 0.0
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
			}
			else
#endif // RTAPI
			{
				aboveWater = GetDynamicSkyColor(refractedDir);
			}

			// The window shows the refracted above-water view, hazing toward
			// the water with DEPTH (not strength) so it dissolves into the blue
			// as it attenuates. The blend is per channel, so the window loses
			// its warm end first and goes blue before it goes dark. Fade to the
			// flat base water color, NOT the fog's inscatter: that carries the
			// sun and the god-ray modulation, which would otherwise show as a
			// god-ray texture inside the window as it hazes out. Strength is
			// the overall opacity, blending back to the normal underwater view:
			const float3 windowContent =
				lerp(ocean.water_color.rgb, aboveWater, pathTransmittance);
			const float windowShape = saturate(cone * upward * openWater);
			const float windowBlend = saturate(windowShape * underwater_snell);
			color.rgb = lerp(color.rgb, windowContent, windowBlend);

			// The window replaces the scene, so it replaces the volumetric
			// light in it too. Looking up at a submerged lamp near the surface
			// is exactly where the window is, so without this the beams would
			// be erased precisely where they matter most.
			volSurvival *= 1 - windowBlend;
		}

		// Give back the share of the volumetric light the stages above took
		// away, restoring it in full. Addition only: the composited value was
		// multiplied by exactly volSurvival, so adding the rest cannot drive
		// anything negative the way subtracting the beams out and re-adding
		// them would - and the half resolution buffer disagreeing slightly with
		// the temporally filtered copy in the input can only soften the result,
		// never ring around it.
		//
		// Near the surface volSurvival tends to 1 and the fully filtered
		// version is what survives; only at depth does this converge on the raw
		// buffer, which is where the water has hidden the detail anyway.
		[branch]
		if (underwater_volumetrics_texture >= 0)
		{
			const float3 volumetrics =
				bindless_textures[descriptor_index(underwater_volumetrics_texture)]
					.SampleLevel(sampler_linear_clamp, uv, 0).rgb;

			color.rgb += volumetrics * (1 - volSurvival);
		}

		//color = float4(1, 0, 0, 1);
	}

	// Constrain the effect smoothly just below the water line:
	color = lerp(color, original_color, smoothstep(0.0, 0.025, saturate(world_pos.y - ocean_pos.y - 0.01)));
	//color = smoothstep(0.0, 0.05, saturate(intersection_direction));
	
	output[DTid.xy] = color;
}
