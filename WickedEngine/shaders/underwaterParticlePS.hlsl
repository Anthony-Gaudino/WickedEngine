/**
 * Underwater particle pixel shader.
 *
 * Shades one mote of suspended detritus. The particles are drawn on both sides
 * of the water surface - the transparent pass issues them once before the ocean
 * and once after, exactly as it does the emitters - so the first thing this
 * does is drop the half that does not belong to the side being drawn.
 */

#include "globals.hlsli"
#include "underwaterHF.hlsli"
#include "underwaterParticleHF.hlsli"

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

	[branch]
	if (GetCamera().texture_depth_index >= 0)
	{
		const float2 screenUV =
			input.position.xy * GetCamera().internal_resolution_rcp;
		const float sceneDepth = compute_lineardepth(
			texture_depth.SampleLevel(sampler_linear_clamp, screenUV, 0));
		const float particleDepth = compute_lineardepth(input.position.z);

		alpha *= (half)saturate(
			(sceneDepth - particleDepth) / SOFT_FADE_DISTANCE);
	}

	return half4(PARTICLE_ALBEDO * alpha, alpha);
}
