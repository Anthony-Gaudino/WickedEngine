#ifndef WI_WATER_SUN_SHAFTS_HF
#define WI_WATER_SUN_SHAFTS_HF
/**
 * The sun's shafts in the water, seen from an eye that is not in it.
 *
 * The froxel volume draws these wherever the eye is under the surface, and
 * cannot draw them at all where it is not. The volume is built along
 * **straight** view rays, and the leg below the surface of a ray that left the
 * water for an eye in the air is refracted - confined by Snell's window to at
 * most 1.512 times the depth, where the straight segment grows without limit at
 * a grazing angle. Cells below the crossing are therefore not on the light's
 * path, by a factor that goes to infinity at the horizon, and no amount of
 * depth resolution moves them onto it.
 *
 * So the shafts are marched here, by the fragment that wants them, over the
 * very leg of water it already measured for its own fog. What comes back is an
 * occlusion fraction rather than a radiance: how much of the sun was blocked,
 * and nothing about what the water then did to it - `MakeWaterFog` owns that
 * and would otherwise account for it twice. That is also what makes this safe
 * where a stored field is not: there is nothing to interpolate between, so no
 * discontinuity at the surface for a reconstruction to straddle.
 *
 * References:
 * https://en.wikipedia.org/wiki/Snell%27s_law
 *
 * Include after `globals.hlsli`.
 */
#include "globals.hlsli"
#include "shadowHF.hlsli"
#include "waterVolumetricsHF.hlsli"

/**
 * How many points along the submerged leg are asked whether the sun reached
 * them.
 *
 * They are spread over whichever is shorter of the leg and the depth the water
 * can be seen into, so they stay close together in exactly the water that
 * extinguishes fastest. What this count buys is the number of distinct edges a
 * single shaft can have across its own depth.
 */
static const uint WATER_SUN_SHAFT_STEPS = 6;

/**
 * How deep into the water the march looks, in optical depths.
 *
 * What the eye receives from a stretch of water falls off as
 * \f$e^{-\sigma_t d}\f$ from the end nearest it, so past a few of these there
 * is nothing left to weigh: at 4 the far end contributes under two percent.
 *
 * Bounding the march rather than the step count is what keeps the answer a
 * property of the water alone. Every fragment deeper than this bound marches
 * the same stretch of it and returns the same veil, so nothing beyond the depth
 * the water can be seen through leaves a trace of itself in the light.
 */
static const float WATER_SUN_SHAFT_OPTICAL_DEPTH = 4.0;

/**
 * Whether the sun reached a point below the surface, and in what colour.
 *
 * One comparison tap, never `shadow_2D`. A soft lookup would blur the very edge
 * the shaft is made of, and at sixteen taps a step it is not affordable in a
 * shader that runs over every submerged fragment on screen. The atlas
 * addressing is written out rather than borrowed for the same reason
 * `sample_shadow` cannot be called directly: its filtering is chosen by macros
 * the including shader sets for its own surface shading, and this wants its own
 * answer either way.
 *
 * **The opaque atlas only.** The transparent layer beside it carries the
 * ocean's caustics, which describe how the surface focuses light onto what it
 * lands on - a property of a lit surface, and already applied to submerged ones
 * by the shading path. Putting them in the water's own veil instead spreads
 * that pattern through the volume at a strength that never falls with depth, so
 * it survives water far too murky to see any distance into.
 *
 * @param[in] light - The directional light to test. Must be casting a shadow;
 *                    there is no atlas entry to read otherwise.
 * @param[in] position - World position to test, which for a submerged sample is
 *                       where the light crossed into the water rather than the
 *                       sample itself. See `WaterVolumetrics::SurfaceEntry`.
 *
 * @return Per-channel surviving fraction, in [0, 1]. Fully lit where the point
 *         falls outside every cascade.
 */
inline half3 WaterSunShaftShadow(in ShaderEntity light, in float3 position)
{
	for (uint cascade = 0; cascade < light.GetShadowCascadeCount(); ++cascade)
	{
		// Ortho matrix, so no divide by w.
		const float3 shadowPosition = mul(
			load_entitymatrix(light.GetMatrixIndex() + cascade),
			float4(position, 1)).xyz;
		float3 shadowUV = clipspace_to_uv(shadowPosition);

		if (!is_saturated(shadowUV))
		{
			continue;
		}

		// Held half a texel inside the cascade's own footprint. The cascades are
		// packed side by side in one atlas, so a tap on the very edge filters
		// into the neighbour and reads a depth from a different projection
		// entirely.
		const float2 cascadeResolution =
			light.shadowAtlasMulAdd.xy * GetFrame().shadow_atlas_resolution;
		shadowUV.xy = clamp(
			shadowUV.xy * cascadeResolution, 0.5, cascadeResolution - 0.5)
			/ cascadeResolution;

		float2 atlasUV = shadowUV.xy;
		atlasUV.x += cascade;
		atlasUV = mad(
			atlasUV, light.shadowAtlasMulAdd.xy, light.shadowAtlasMulAdd.zw);

		Texture2D<half4> shadowAtlas = bindless_textures_half4[
			descriptor_index(GetFrame().texture_shadowatlas_index)];

		half3 shadow = shadowAtlas.SampleCmpLevelZero(
			sampler_cmp_depth, atlasUV, shadowPosition.z).rrr;

		if (GetFrame().options & OPTION_BIT_VOLUMETRICCLOUDS_CAST_SHADOW)
		{
			shadow *= shadow_2D_volumetricclouds(position);
		}

		return shadow;
	}

	return 1;
}

/**
 * How much of the sun's glow in a stretch of water actually survived.
 *
 * Walks down the submerged leg of a view ray from the surface and asks, at each
 * point, whether anything stood between it and the sun - then averages those
 * answers weighted by how much each point contributes to what the eye receives:
 * \f[
 * M = \frac{\sum_i T_{view}(d_i)\,\mathrm{shadow}(d_i)}
 *          {\sum_i T_{view}(d_i)}
 * \f]
 * for \f$d_i\f$ measured from the surface, which is the end nearest the eye.
 *
 * **Occlusion only, and deliberately.** What the water itself does to sunlight
 * on the way down belongs to the medium, and `MakeWaterFog` accounts for it
 * already - over the eye's own depth, and through the share of the segment that
 * was extinguished. Attenuating here as well would apply the same absorption
 * twice, and in water murky enough to matter that is the difference between a
 * sun and no sun at all. So this answers one question: how much of the sun was
 * blocked.
 *
 * **The weighting is not decoration.** Water takes red out of a beam roughly an
 * order of magnitude faster than blue, so an even average would let a shadow
 * deep in the column darken a red channel that never reached the eye from down
 * there at all. Weighting by `ViewTransmittance` is what keeps each channel's
 * answer to the water it can actually see through, and it is why this is a
 * `float3` rather than a scalar. It also settles what a shaft looks like with
 * depth: the sun's entry points for successive samples spread out across the
 * surface, so a deeper fragment averages a wider patch of it and its shadows
 * arrive softer - the same thing many scattering events do to real water, out
 * of the geometry rather than out of a fudge factor.
 *
 * Sampling at step centres is not an approximation of the closed form
 * `StepIntegralWeight` gives: over equal steps the two differ by
 * \f$e^{-\sigma_t \Delta / 2}\f$, the same factor on every step, and this is a
 * ratio - so the constant divides out and the midpoint rule is exact here.
 *
 * Example usage:
 * @code
 * const float3 modulation = WaterSunShaftModulation(
 *     medium, surfaceCrossing, toEye, path, pixel);
 * @endcode
 *
 * @param[in] medium - The water the leg lies in.
 * @param[in] surfaceCrossing - World position where this view ray enters the
 *                              water, which is the NEAR end of the leg. Taken
 *                              from the ray rather than from the fragment so
 *                              that everything drawn along one ray marches the
 *                              same water.
 * @param[in] toEye - Normalized direction from the fragment towards the eye,
 *                    taken above the surface.
 * @param[in] path - Length of the leg lying in water (in metres), from
 *                   `SubmergedViewPath`. An upper bound on what is walked: past
 *                   a few optical depths there is nothing left to weigh, and
 *                   the march stops there instead.
 * @param[in] pixel - Pixel coordinate, for spreading the sample depths.
 *
 * @return Per-channel unoccluded fraction of the sun over this water, in
 *         [0, 1]. Returns 1 - nothing in the way - wherever there is no leg or
 *         no sun to march it against, so a caller needs no test of its own.
 *
 * @note For an eye **in the air**. The leg is taken along the refracted
 *       direction, which is the path light took to leave the water; an eye
 *       under the surface sees along a straight one, and the froxel volume
 *       carries that case with a shadow at every cell rather than one leg.
 */
inline float3 WaterSunShaftModulation(
	in WaterVolumetrics medium,
	in float3 surfaceCrossing,
	in float3 toEye,
	in float path,
	in uint2 pixel
)
{
	// Nothing to shape. Anything shorter than this is a fragment sitting
	// essentially on the surface, where the whole term is vanishing anyway.
	[branch]
	if (path < 0.01)
	{
		return 1;
	}

	// The same sun `MakeWaterFog` builds its glow from, which the renderer takes
	// from the dominant directional light - so the first one carrying
	// volumetrics is the entity behind that colour and direction, and the one
	// whose shadow map has anything to say about it.
	ShaderEntityIterator iterator = directional_lights();
	for (uint index = iterator.first_item();
		index <= iterator.last_item() && !iterator.empty();
		++index)
	{
		ShaderEntity light = load_entity(index);

		if (!light.IsVolumetricsEnabled())
		{
			continue;
		}

		// No atlas entry to read, so nothing can cut this light into shafts.
		// The glow it makes is still there, whole.
		if (!light.IsCastingShadow())
		{
			return 1;
		}

		// The direction light travelled through the water on its way up to the
		// eye, which is the leg this fog was measured over. Reciprocity: bending
		// the direction back towards the eye into the water gives the direction
		// the light came up along.
		const float3 legDirection = RefractIntoWater(toEye);

		// Where the sun points once the surface has bent it, which is fixed for
		// the whole scene - the interface is a plane and the source is at
		// infinity - so each step costs one lookup and no geometry.
		const float3 toSun = RefractIntoWater((float3)GetSunDirection());

		// How far down the leg is worth walking. Set by the channel that is
		// extinguished SLOWEST, since it is the last one still carrying
		// anything: cutting the march to a faster one would end it while that
		// channel still had something to say.
		const float slowestExtinction = min3(
			medium.sigmaT.r, medium.sigmaT.g, medium.sigmaT.b);

		const float marched = min(
			path,
			WATER_SUN_SHAFT_OPTICAL_DEPTH / max(slowestExtinction, 0.00001));

		const float stepSize = marched / (float)WATER_SUN_SHAFT_STEPS;

		// Where inside its step each sample sits. Occlusion is a step function
		// - the shadow map says lit or unlit and nothing between - so six
		// samples of it can only ever return seven values, and a shadow's edge
		// crossing the water arrives as that many terraces. Moving the sample
		// depths per pixel trades the terraces for a fine noise that the
		// temporal pass averages back out.
		//
		// One offset for the whole march, so the samples stay evenly spread
		// rather than clumping. Blue noise rather than an ordered pattern: this
		// carries its own per-frame phase, and an ordered one would hold still
		// and print itself on the screen as a grid instead of averaging away.
		const float jitter = blue_noise(pixel).x;

		float3 weightSum = 0;
		float3 litSum = 0;

		for (uint step = 0; step < WATER_SUN_SHAFT_STEPS; ++step)
		{
			// Measured DOWN from the surface, so a step stands at the same depth
			// whatever is drawn behind it.
			const float depthAlongLeg = ((float)step + jitter) * stepSize;

			const float3 samplePosition =
				surfaceCrossing - legDirection * depthAlongLeg;

			// The water this step is seen through, which is the stretch of it
			// between here and the surface. The near end counts for most,
			// steeply so once the water is at all murky.
			const float3 weight = medium.ViewTransmittance(depthAlongLeg);

			// Read where the sun entered the water, never at the sample. What
			// occludes a submerged point does so above the surface, on the part
			// of the path that is still straight, so a shadow map that knows
			// nothing of refraction is exactly right there.
			const half3 shadow = WaterSunShaftShadow(
				light,
				medium.SurfaceEntry(samplePosition, toSun, FLT_MAX));

			weightSum += weight;
			litSum += weight * (float3)shadow;
		}

		return litSum / max(weightSum, 0.00001);
	}

	return 1;
}

#endif // WI_WATER_SUN_SHAFTS_HF
