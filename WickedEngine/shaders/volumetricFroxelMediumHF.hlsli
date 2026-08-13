#ifndef WI_VOLUMETRIC_FROXEL_MEDIUM_HF
#define WI_VOLUMETRIC_FROXEL_MEDIUM_HF
/**
 * What fills one cell of the volumetric froxel volume.
 *
 * Air or water, and what either of them does to a light passing through it.
 *
 * Kept out of `volumetricFroxelHF.hlsli`, which every shader that draws
 * anything includes and which has to stay small. This pulls in the ocean and
 * the entire water medium behind it, and only the two passes that build the
 * volume need any of that.
 *
 * Include after `globals.hlsli` and `fogHF.hlsli`.
 */
#include "volumetricFroxelHF.hlsli"
#include "waterVolumetricsHF.hlsli"

/**
 * Extinction of the air at a point, in 1/m.
 *
 * The height fog's own density, read as the **local** quantity it is built from
 * rather than as the integral along a ray that `GetFogAmount` returns. A cell
 * needs the density at a place; the integral over the path is what the
 * integration pass is for, and asking `GetFogAmount` here would count the same
 * air once per slice.
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
 * The medium filling one cell, and what the lights do in it.
 *
 * Holds both rather than choosing between them, so that one expression covers
 * either and neither has to know which the cell turned out to be. Only one of
 * them is ever present at a time: air and water are alternatives at a point,
 * not a mixture.
 *
 * Example usage:
 * @code
 * const VolumetricFroxelMedium medium = VolumetricFroxelMediumAt(uv, depth);
 * if (medium.Scatters())
 * {
 *     scattered = lightColor * shadow
 *         * medium.Scatter(P, L, medium.Bend(L), toEye, FLT_MAX, L.y);
 * }
 * @endcode
 */
struct VolumetricFroxelMedium
{
	/**
	 * The ocean's water, inactive where the eye is not under it.
	 *
	 * Carried whole rather than reduced to coefficients, because what the water
	 * does to a light is not only extinction: it bends the sun on the way in,
	 * attenuates each channel over the path the light took to get here, and
	 * scatters through its own strongly forward phase function.
	 */
	WaterVolumetrics water;

	/**
	 * Height fog's extinction here, zero where the water has the cell (1/m).
	 *
	 * Non-absorbing, so this is its scattering coefficient too - which is what
	 * the height fog already assumes by having no absorption term of its own to
	 * separate out.
	 */
	float airSigmaT;

	/** Extinction of everything in the cell, per channel (1/m). */
	float3 sigmaT;

	/**
	 * Whether anything here turns light towards the eye.
	 *
	 * @return true when there is a medium worth evaluating the lights for.
	 */
	bool Scatters()
	{
		return airSigmaT > 0 || water.IsActive();
	}

	/**
	 * The direction a light's rays actually arrive from.
	 *
	 * Sunlight bends at the surface before it reaches anything below it, so a
	 * submerged cell is lit from steeper than the direction towards the sun.
	 * Shared with surface shading, so a shaft and the sea bed it lands on agree
	 * about where the sun is.
	 *
	 * @param[in] toLight - Normalized direction from the cell to the light.
	 *
	 * @return The direction to sample and shade along.
	 *
	 * @note For a **directional** light only. Snell's law at a plane has one
	 *       solution for a source at infinity and needs a per-cell one for a
	 *       source that is not; a local light is left unbent, and its path is
	 *       clipped at the surface either way.
	 */
	float3 Bend(float3 toLight)
	{
		return water.IsActive() ? RefractIntoWater(toLight) : toLight;
	}

	/**
	 * Where to test a shadow map for a light shining into this cell.
	 *
	 * Whatever occludes a submerged cell does so above the surface, on the part
	 * of the path that is still straight - so a shadow map that knows nothing
	 * of refraction is exactly right when it is read at the point the light
	 * enters the water. That is also where the ocean writes its caustics into
	 * the transparent shadow layer, so they land where they form.
	 *
	 * @param[in] samplePos - World position of the cell.
	 * @param[in] toLight - Normalized direction to the light, already bent.
	 * @param[in] distanceToLight - Distance to the light, or `FLT_MAX` for a
	 *                              directional light.
	 *
	 * @return World position to transform into the light's shadow space.
	 */
	float3 ShadowSamplePosition(
		float3 samplePos,
		float3 toLight,
		float distanceToLight
	)
	{
		return water.IsActive()
			? water.SurfaceEntry(samplePos, toLight, distanceToLight)
			: samplePos;
	}

	/**
	 * Light turned towards the eye here, per metre of view ray.
	 *
	 * **Each medium is given the direction the light reaches it from**, which
	 * is why both are asked for. They differ only under water, and only for the
	 * sun, but the difference is not small: the air's phase function is sharply
	 * forward peaked, so handing it a direction bent by the surface changes how
	 * bright the fog is by a large factor. Since the bending switches on the
	 * moment a cell holds any water at all, sharing one direction between them
	 * puts a step in the **air** fog wherever a cell so much as grazes the
	 * surface - which is a step in a term that should not have moved.
	 *
	 * @param[in] samplePos - World position of the cell.
	 * @param[in] toLight - Normalized direction to the light through air.
	 * @param[in] toLightBent - The same direction as it travels through the
	 *                          water, from `Bend`. Equal to `toLight` for a
	 *                          local light, which is not bent.
	 * @param[in] toEye - Normalized direction from the cell towards the eye.
	 * @param[in] distanceToLight - Distance to the light, or `FLT_MAX` for a
	 *                              directional light.
	 * @param[in] cosThetaAbove - Cosine of the angle to the light taken ABOVE
	 *                            the surface, before any bending. Ignored in
	 *                            air, and for a light that is itself submerged.
	 *
	 * @return Per-channel radiance scattered towards the eye per metre, before
	 *         the light's own colour, falloff and shadow term.
	 */
	float3 Scatter(
		float3 samplePos,
		float3 toLight,
		float3 toLightBent,
		float3 toEye,
		float distanceToLight,
		float cosThetaAbove
	)
	{
		// Nothing attenuates the light on its way to a cell of air: the height
		// fog carries no record of where it has been, and the shadow map is
		// what stands in for the occlusion that matters.
		float3 scattered =
			airSigmaT * ComputeScattering(saturate(dot(toLight, -toEye)));

		[branch]
		if (water.IsActive())
		{
			scattered += water.InScatter(
				samplePos, toLightBent, toEye, distanceToLight, cosThetaAbove);
		}

		return scattered;
	}
};

/**
 * Evaluates the medium filling a cell.
 *
 * **The water is decided once for a whole column, at the eye, and never
 * changes along it.** The volume cannot carry the surface itself, and this is a
 * limit of its depth resolution rather than of the arithmetic.
 *
 * Water scatters one to two orders of magnitude more per metre than height fog,
 * so where a column crosses the surface the light gathered along it has a slope
 * discontinuity of that size. Cells are metres deep - about four at 125 m and
 * eight at 490 - and the volume stores one running total per cell, which a
 * fragment reads by interpolating between neighbours. Interpolating across a
 * hundred-fold change in slope leaves an error of order the whole water term,
 * and, because the crossing sweeps outward through cell after cell as the view
 * angle changes, that error alternates in sign once per cell. It draws
 * concentric rings on the one surface that always sits exactly at the
 * discontinuity: the drawn ocean.
 *
 * No placement of the boundary within a cell escapes it. Splitting a cell at
 * the crossing makes the stored totals exact and leaves the interpolation
 * between them just as wrong, because the kink is still inside a cell the
 * fragment reads across.
 *
 * @note **What this gives up** is a column that straddles - a camera above the
 *       water sees no shafts below the surface, and one crossing the waterline
 *       changes medium over the span `ocean_underwater_factor` sweeps rather
 *       than exactly at the surface. What it keeps is the case the water needs:
 *       an eye under the surface has water the whole way and gets its shafts,
 *       its caustics and its per-channel attenuation with no discontinuity
 *       anywhere in the volume.
 *
 * @note Shares `GetWaterVolumetricsAtEye` with the surface shading and the
 *       underwater post pass, so all of them switch on the same pixels at the
 *       same instant. That is worth more here than being right about where the
 *       surface is: an effect that disagrees with the others about which side
 *       of the water a pixel is on is visible immediately, and being a metre
 *       out about a surface no fragment is drawn on is not visible at all.
 *
 * @param[in] uv - Screen space position of the cell (0-1).
 * @param[in] depth - Radial distance from the eye to the cell (in metres).
 *
 * @return The medium filling that cell.
 */
inline VolumetricFroxelMedium VolumetricFroxelMediumAt(
	in float2 uv, in float depth)
{
	VolumetricFroxelMedium medium;

	medium.water = GetWaterVolumetricsAtEye(uv);

	// The air's own density does vary along a column, unlike the water's, so it
	// is asked at the cell rather than at the eye. It fills the column whenever
	// the water does not, so between them they account for it once.
	const float3 position =
		GetCamera().position + VolumetricFroxelRayDirection(uv) * depth;

	medium.airSigmaT = VolumetricFroxelAirExtinction(position, depth)
		* (1 - (float)medium.water.submersion);

	medium.sigmaT = medium.water.sigmaT + medium.airSigmaT;

	return medium;
}

#endif // WI_VOLUMETRIC_FROXEL_MEDIUM_HF
