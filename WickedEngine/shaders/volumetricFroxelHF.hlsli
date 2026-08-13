#ifndef WI_VOLUMETRIC_FROXEL_HF
#define WI_VOLUMETRIC_FROXEL_HF
/**
 * Shared geometry of the volumetric froxel volume.
 *
 * Everything that builds the volume and everything that reads it derives its
 * coordinates from here. Two places deriving the same mapping is how a lookup
 * drifts off the thing it is looking up, and the symptom - a smooth,
 * depth-dependent offset - reads as a tuning problem rather than as a bug.
 *
 * Include after `globals.hlsli`.
 */
#include "ShaderInterop_VolumetricFroxels.h"

/**
 * Radial distance from the eye to a slice boundary, in metres.
 *
 * Squared rather than linear, matching the aerial perspective volume this sits
 * beside. It spends resolution near the eye, which is where a medium is thick
 * enough to matter and where the light that crossed water still survives,
 * without collapsing the far field the way an exponential distribution does:
 * over 500 m in 128 slices the first cell is 3 cm deep and the last 8 m,
 * where an exponential of the same span would end in cells tens of metres
 * across.
 *
 * **Radial distance, not view depth**, because the volume is built along
 * normalized view rays. Anything sampling it has to measure the same way or the
 * error grows towards the edges of the screen, where the two diverge most.
 *
 * @param[in] slice - Slice coordinate, 0 at the eye and
 *                    `VOLUMETRIC_FROXEL_SLICES` at the far end.
 * @param[in] range - How far the volume reaches, in metres.
 *
 * @return Radial distance from the eye to that boundary (in metres).
 */
inline float VolumetricFroxelSliceToDepth(in float slice, in float range)
{
	const float normalized = slice / (float)VOLUMETRIC_FROXEL_SLICES;

	return range * normalized * normalized;
}

/**
 * Texture W coordinate for a distance from the eye.
 *
 * The inverse of `VolumetricFroxelSliceToDepth`, in the 0-1 a sampler wants.
 * Saturated here rather than left to the sampler's clamp, so the far behaviour
 * is stated next to the mapping it belongs to: everything past the volume's
 * range reads the last slice, which is what lets that slice carry the whole
 * remaining column.
 *
 * @param[in] depth - Radial distance from the eye (in metres).
 * @param[in] range - How far the volume reaches, in metres.
 *
 * @return W coordinate in 0-1.
 */
inline float VolumetricFroxelDepthToW(in float depth, in float range)
{
	return sqrt(saturate(depth / max(range, 0.0001)));
}

/**
 * World position of a point inside the volume.
 *
 * Built by unprojecting the cell's screen position onto the far plane and
 * walking the resulting ray, so it inherits the camera's own projection rather
 * than reconstructing one - which is what keeps the volume aligned with the
 * fragments that will sample it.
 *
 * @param[in] uv - Screen space position of the cell centre (0-1).
 * @param[in] depth - Radial distance from the eye (in metres).
 *
 * @return World position of that point.
 */
inline float3 VolumetricFroxelPosition(in float2 uv, in float depth)
{
	float4 unprojected = mul(
		GetCamera().inverse_view_projection,
		float4(uv_to_clipspace(uv), 0, 1));
	unprojected.xyz /= unprojected.w;

	const float3 rayOrigin = GetCamera().position;
	const float3 rayDirection = normalize(unprojected.xyz - rayOrigin);

	return rayOrigin + rayDirection * depth;
}

/**
 * Where inside a cell to take that cell's single sample, in 0-1 per axis.
 *
 * A cell stands for the average of the medium across its whole volume, but only
 * one sample is affordable, so where that sample sits decides what the cell
 * reports. Take them all at cell centres and every cell in a slice tests the
 * shadow map at the same distance at once, and the slice boundary becomes a
 * shell the eye can see: a shaft that should be a straight line is cut into
 * repeating bands, one per slice. Spreading the sample points inside their
 * cells trades that shell for a fine pattern which the volume's own bilinear
 * upsampling then averages away.
 *
 * An ordered Bayer pattern rather than noise, and one that does not move
 * between frames, so what is left after the averaging is stationary. Random
 * offsets would need several frames of accumulation to be worth more than this,
 * and would crawl until they got it.
 *
 * The slice index is folded into the lookup on both axes with different
 * strides, so neighbouring slices draw different offsets instead of stacking
 * the same 8x8 tile through the depth of the volume.
 *
 * @param[in] cell - Integer cell coordinate within the volume.
 *
 * @return Offset within the cell, in 0-1 on each axis.
 */
inline float3 VolumetricFroxelCellJitter(in uint3 cell)
{
	const uint2 xyOffset = uint2(cell.z * 5, cell.z * 3);
	const uint2 zOffset = uint2(cell.z * 3 + 2, cell.z * 7 + 5);

	return float3(
		(float)dither((min16uint2)(cell.xy + xyOffset)),
		(float)dither((min16uint2)(cell.yx + xyOffset)),
		(float)dither((min16uint2)(cell.xy + zOffset)));
}

/**
 * Extinction of the air at a point, in 1/m.
 *
 * The height fog's own density, read as the **local** quantity it is built from
 * rather than as the integral along a ray that `GetFogAmount` returns. A cell
 * needs the density at a place; the integral over the path is what the
 * integration pass is for, and asking `GetFogAmount` here would count the same
 * air once per slice.
 *
 * Shared by both build passes deliberately. They must agree about where the air
 * is thick, or the light one pass injects is attenuated by a medium the other
 * does not think is there.
 *
 * References:
 * https://www.iquilezles.org/www/articles/fog/fog.htm
 *
 * @param[in] position - World position to evaluate at.
 * @param[in] distanceToEye - Radial distance from the eye (in metres), for the
 *                            authored near fade.
 *
 * @return Extinction coefficient there (1/m).
 */
inline float VolumetricFroxelAirExtinction(
	in float3 position, in float distanceToEye)
{
	const ShaderFog fog = GetWeather().fog;

	// Artistic, not physical: fog is held off until `start` and ramped in over
	// the same distance again. Evaluated at the point's own distance, which is
	// a question a cell can answer exactly - unlike a march, which has to pick
	// one representative point for a whole segment.
	const float startFalloff = saturate((distanceToEye - fog.start) / fog.start);

	float density = fog.density * startFalloff;

	[branch]
	if (GetFrame().options & OPTION_BIT_HEIGHT_FOG)
	{
		const float falloffScale =
			rcp(max(0.01, fog.height_end - fog.height_start));

		// Solve e^(-h * x) = 0.001 for x, matching `GetFogAmount` so the two
		// descriptions of the same fog cannot disagree about where it thins.
		const float falloff = 6.907755 * falloffScale;

		density *= exp(-max(position.y - fog.height_start, 0) * falloff);
	}

	return density;
}

/**
 * Adds the light scattered between the eye and a fragment.
 *
 * This is what the volume exists for. The inscatter reaching a fragment depends
 * only on where that fragment is, so every kind of fragment - opaque, glass,
 * particle, sprite, splat, ocean surface - takes exactly its own share, and
 * there is no compositing order left to choose because nothing is composited.
 *
 * **The weight is a coverage fraction, not a radiance.** `ApplyFog` and
 * `ApplyWaterFog` take a `background` radiance and use `min(background, color)`
 * to find the share to hold out of a *multiplicative* operation. That trick does
 * not carry over to an additive term, which needs to know what *fraction* of the
 * pixel came from behind. Wherever a shader composites a refraction it has that
 * fraction to hand - `objectHF.hlsli` builds its background as
 * `surface.refraction.rgb * (1 - surface.F) * surface.refraction.a`, so the
 * weight is `(1 - surface.F) * surface.refraction.a`.
 *
 * **Alpha blended draws pass 0.** There the hardware multiplies by the source
 * alpha and the destination carries its own share already, so weighting here as
 * well would count the coverage twice.
 *
 * **Additive draws must not call this at all.** The destination already holds
 * the inscatter for that pixel, and adding it once per additive draw multiplies
 * the haze by however many of them land there - the same reasoning
 * `ApplyWaterFogAdditive` is built on.
 *
 * A camera with no volume - a reflection, an environment probe, a shadow
 * cascade - carries -1 and returns here, so no call site needs a condition of
 * its own.
 *
 * Example usage:
 * @code
 * ApplyWaterFog(dist, color);
 * ApplyVolumetricLight(ScreenCoord, dist, 0, color);
 * @endcode
 *
 * Takes the world position rather than a distance, matching `ApplyWaterFog`
 * beside it. The volume is indexed by **radial** distance from the eye, which
 * differs from view depth by more and more towards the edges of the screen; a
 * caller handed that job would eventually pass the wrong one, and a smooth
 * error that grows towards the corners reads as a tuning problem rather than as
 * a bug.
 *
 * @param[in] screenUV - Screen space position of the fragment (0-1).
 * @param[in] fragmentPosition - World position of the fragment.
 * @param[in] backgroundWeight - Fraction of this fragment's colour that came
 *                               from behind it and already carries its own
 *                               inscatter (0-1). Per channel, because a
 *                               refraction's Fresnel term is, and zero unless
 *                               the shader composited a refraction itself.
 * @param[in,out] color - Fragment colour, lit in place. Alpha is untouched.
 */
inline void ApplyVolumetricLight(
	in float2 screenUV,
	in float3 fragmentPosition,
	in half3 backgroundWeight,
	inout half4 color
)
{
	[branch]
	if (GetCamera().texture_volumetricfroxels_index < 0)
	{
		return;
	}

	const float distanceToEye = length(fragmentPosition - GetCamera().position);

	const float w = VolumetricFroxelDepthToW(
		distanceToEye, GetFrame().volumetricfroxel_range);

	const half3 inscatter = (half3)texture_volumetricfroxels.SampleLevel(
		sampler_linear_clamp, float3(screenUV, w), 0).rgb;

	color.rgb += inscatter * (1 - backgroundWeight);
}

/**
 * Adds the light scattered between the eye and a fragment that has no
 * already-lit background of its own.
 *
 * @param[in] screenUV - Screen space position of the fragment (0-1).
 * @param[in] fragmentPosition - World position of the fragment.
 * @param[in,out] color - Fragment colour, lit in place.
 */
inline void ApplyVolumetricLight(
	in float2 screenUV, in float3 fragmentPosition, inout half4 color
)
{
	ApplyVolumetricLight(screenUV, fragmentPosition, 0, color);
}

/**
 * Adds the light scattered between the eye and a fragment drawn with
 * premultiplied alpha.
 *
 * The blend adds the source unweighted, so the coverage has to be applied here
 * instead - otherwise a mote covering a tenth of a pixel paints a shaft across
 * the nine tenths it never touched. Same reasoning as
 * `ApplyWaterFogPremultiplied`.
 *
 * @param[in] screenUV - Screen space position of the fragment (0-1).
 * @param[in] fragmentPosition - World position of the fragment.
 * @param[in,out] color - Premultiplied fragment colour, lit in place.
 */
inline void ApplyVolumetricLightPremultiplied(
	in float2 screenUV, in float3 fragmentPosition, inout half4 color
)
{
	ApplyVolumetricLight(screenUV, fragmentPosition, 1 - color.a, color);
}

#endif // WI_VOLUMETRIC_FROXEL_HF
