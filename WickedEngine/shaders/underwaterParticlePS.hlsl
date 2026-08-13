/**
 * Underwater particle pixel shader.
 *
 * Shades one mote of suspended detritus: a small, dull, rough scatterer lit by
 * whatever reaches it through the water, and fogged by the water between it
 * and the eye.
 *
 * Both of those come free from the shared lighting path. Every light already
 * arrives dimmed and shifted by the water it crossed, and the ambient term
 * already arrives through the column of water standing over the point, so a
 * mote beside a lamp is lit in the lamp's colour while one twenty metres off
 * is not - which is the whole reason the field is worth drawing.
 */

#define TRANSPARENT // uses transparent light lists

// Note the absence of DISABLE_TRANSPARENT_SHADOWMAP, which the emitters set.
// Their reason is that a camera facing sprite self-shadows badly against the
// layer it casts into; these particles are never drawn in a shadow pass, so
// they cast nothing and cannot self-shadow. Keeping the layer is the point:
// the ocean surface writes its caustics there rather than any occlusion, and
// motes lighting up as the filaments sweep over them is the clearest thing
// the field has to show.

#include "globals.hlsli"
// Not objectHF.hlsli, which the emitters use: it declares the object push
// constant block, and a shader may only have one. This field carries its own.
#include "shadingHF.hlsli"
#include "waterFogHF.hlsli"
#include "underwaterParticleHF.hlsli"
#include "volumetricFroxelHF.hlsli"

/**
 * Depth range over which a particle fades out against solid geometry, in
 * metres.
 *
 * A particle is a flat sprite, so where it meets the sea bed it would
 * otherwise cut a hard line across it. Measured in world space rather than
 * against the particle's own size, because a mote is millimetres across and
 * fading over millimetres would hide nothing.
 */
static const float SOFT_FADE_DISTANCE = 0.25;

/**
 * Reflectance of a particle.
 *
 * Marine snow is a loose aggregate of dead cells and mucus, which is a dull
 * off-white rather than anything bright.
 */
static const half3 PARTICLE_ALBEDO = half3(0.75, 0.73, 0.68);

/**
 * Specular reflectance of a particle at normal incidence.
 *
 * Water's own value. An aggregate this sodden is barely a separate substance
 * from the water around it, and it presents no smooth face to reflect from in
 * any case - the roughness below leaves almost nothing of this.
 */
static const half PARTICLE_F0 = 0.02;

half4 main(UnderwaterParticleVertexToPixel input) : SV_TARGET
{
	// The particles are suspended IN the water, so nothing above the surface
	// is one. This is a stronger statement than the water side test the other
	// transparents make, and it replaces it: those are drawn on both sides of
	// the water and only have to be told which half is theirs, whereas a
	// particle in the air is not a particle at all.
	//
	// The DRAWN surface, so that the waterline the particles stop at is the
	// one the ocean put on screen rather than waves the mesh flattened away.
	//
	// Per fragment rather than per particle: a sprite is a disc with extent,
	// so one straddling the waterline is partly submerged and cannot be
	// classified whole.
	clip(ocean_drawn_surface_height(input.P) - input.P.y);

	// A round mote out of a square sprite, soft edged so it does not alias
	// into a diamond as it shrinks.
	const float radialSq = dot(input.corner, input.corner);
	clip(1 - radialSq);
	half alpha = (half)(1 - radialSq) * input.opacity;

	const float2 screenUV =
		input.position.xy * GetCamera().internal_resolution_rcp;

	[branch]
	if (GetCamera().texture_depth_index >= 0)
	{
		const float sceneDepth = compute_lineardepth(
			texture_depth.SampleLevel(sampler_linear_clamp, screenUV, 0));
		const float particleDepth = compute_lineardepth(input.position.z);

		alpha *= (half)saturate(
			(sceneDepth - particleDepth) / SOFT_FADE_DISTANCE);
	}

	// Shade the sprite as the sphere it stands for rather than as the disc it
	// is, so a mote has a lit side and a dark one and the field gains some
	// shape under a nearby lamp. In view space the visible hemisphere is the
	// one facing the camera, which is the negative z half.
	const float3 viewNormal =
		float3(input.corner, -sqrt(saturate(1 - radialSq)));

	float3 toEye = GetCamera().position - input.P;
	const float distanceToEye = length(toEye);
	toEye /= max(distanceToEye, 0.00001);

	Surface surface;
	surface.init();
	surface.P = input.P;
	surface.N = (half3)normalize(
		mul((float3x3)GetCamera().inverse_view, viewNormal));
	surface.V = (half3)toEye;
	surface.albedo = PARTICLE_ALBEDO;
	surface.f0 = PARTICLE_F0;
	surface.roughness = 1;
	surface.opacity = alpha;
	surface.pixel = (min16uint2)input.position.xy;
	surface.screenUV = screenUV;
	surface.update();

	Lighting lighting;
	lighting.create(0, 0, GetAmbient(surface.N), 0);

	TiledLighting(surface, lighting, GetFlatTileIndex(surface.pixel));

	half4 color = half4(0, 0, 0, alpha);
	ApplyLighting(surface, lighting, color);

	// Premultiplied from here, matching the blend state. The water's veiling
	// light then has to be scaled by the coverage too, or a mote covering a
	// tenth of a pixel would paint haze across the nine tenths it never
	// touched - which is why this takes the premultiplied fog rather than the
	// plain one the opaque geometry uses.
	color.rgb *= color.a;
	ApplyWaterFogPremultiplied(screenUV, input.P, color);
	ApplyVolumetricLightPremultiplied(screenUV, input.P, color);

	return color;
}
