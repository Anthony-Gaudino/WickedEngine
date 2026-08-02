#ifndef WI_WATERVOLUMETRICS_HF
#define WI_WATERVOLUMETRICS_HF

/**
 * The ocean's water as a participating medium, for volumetric light marches.
 *
 * Provides `WaterVolumetrics`, which answers the three questions a march has to
 * ask of a medium at each sample: how much of the light arriving there is
 * turned towards the eye, how much of it was lost on the way in, and how much
 * of what was turned survives the trip back out. The engine's volumetric light
 * passes ask those of the height fog; below the waterline they ask them here.
 *
 * These are deliberately stateless leaf functions rather than a march of their
 * own. HLSL has no function pointers, and the four volumetric light shaders
 * differ in exactly the places a shared loop would have to abstract over -
 * cascades versus cube maps versus a spot atlas, and four different falloffs -
 * so each keeps its own loop and only the medium is shared.
 *
 * Deliberately free of any shadow sampling, so this stays includable from
 * `lightingHF.hlsli` - which the callers themselves include - without a cycle.
 *
 * @note Callers that sample the shadow cascades also read the *transparent*
 *       shadow layer, which is where the ocean surface writes its caustics
 *       (see the `SHADOWMAPRENDERING` branch of `oceanSurfacePS.hlsl`). Their
 *       shafts therefore pick up caustic modulation without asking for it.
 *       Only the first couple of cascades render the ocean, so that modulation
 *       fades with distance.
 */

#include "globals.hlsli"
#include "skyAtmosphere.hlsli"

/**
 * Depth over which the water medium fades in as the camera submerges (metres).
 *
 * The volumetric light passes clamp their march at the water plane, so a ray
 * looking down flips from a segment that ends at the surface to one running all
 * the way to the sea bed the moment the eye crosses it. Ramping the medium in
 * over the first metre keeps that from landing as a step change in brightness.
 */
static const float WATER_VOLUMETRICS_FADE_DEPTH = 1.0;

/**
 * The ocean's water medium, evaluated for one camera position.
 *
 * Single scattering: what a sample sends to the eye is the light that reaches
 * it, turned through the phase function and attenuated on both legs:
 * \[
 * dL = T_{view}(s)\, \sigma_s\, p(\theta)\, T_{light}(s)\, L_{light}\, ds
 * \]
 * `InScatter` returns the middle three factors and `ViewTransmittance` the
 * first; the shadow term, the distance falloff and the light's colour stay with
 * the caller, since those are what differ between light types.
 *
 * Example usage:
 * @code
 * const WaterVolumetrics water = GetWaterVolumetrics(GetCamera().position);
 * [branch]
 * if (water.IsActive())
 * {
 *     attenuation *= water.InScatter(P, L, V, distanceToLight)
 *         * water.ViewTransmittance(cameraDistance - marchedDistance)
 *         * stepSize;
 * }
 * @endcode
 *
 * References:
 * https://www.pbr-book.org/3ed-2018/Volume_Scattering
 */
struct WaterVolumetrics
{
	/** Scattering coefficient, per RGB channel (1/m). */
	float3 sigmaS;

	/** Extinction, absorption plus scattering, per RGB channel (1/m). */
	float3 sigmaT;

	/**
	 * Henyey-Greenstein asymmetry, already reduced by similarity theory.
	 *
	 * Marine particulates scatter far forward, but that is the phase of a
	 * single event; a single-scattering march standing in for a multiply
	 * scattering medium should use the reduced value.
	 */
	half phaseG;

	/** World height of the still water plane. */
	float waterHeight;

	/** 0 at or above the surface, reaching 1 once fully submerged. */
	half submersion;

	/**
	 * Whether the eye is under water, so this medium applies at all.
	 *
	 * @return true when the camera is below the surface.
	 */
	bool IsActive()
	{
		return submersion > 0;
	}

	/**
	 * Beer-Lambert transmittance from the eye out to a sample.
	 *
	 * Analytic rather than a running product, because the volumetric light
	 * passes march from the far end of the segment back towards the camera - a
	 * product accumulated in that order would attenuate from the wrong end, and
	 * would have to be advanced even on samples a cone or shadow test rejects.
	 *
	 * @param[in] distanceFromEye - Distance from the eye to the sample (in
	 *                              metres).
	 *
	 * @return Per-channel surviving fraction, in [0, 1].
	 */
	float3 ViewTransmittance(float distanceFromEye)
	{
		return exp(-distanceFromEye * sigmaT);
	}

	/**
	 * Length of the light's path to a sample that lies under water.
	 *
	 * One formula for every light type. A submerged light is reached before the
	 * path can leave the water, so the whole distance counts; one above the
	 * surface clips at the crossing and only the part below it attenuates. A
	 * directional light passes `FLT_MAX` and gets the slant path up to the
	 * surface, which is why no separate sun case is needed.
	 *
	 * @param[in] samplePos - Sample position, in world space.
	 * @param[in] toLight - Normalized direction from the sample to the light.
	 * @param[in] distanceToLight - Distance to the light, or `FLT_MAX` for a
	 *                              directional light.
	 *
	 * @return Length of the submerged part of the path (in metres).
	 */
	float SubmergedLightPath(
		float3 samplePos,
		float3 toLight,
		float distanceToLight
	)
	{
		// A path heading down or along never reaches the surface. Clamp to the
		// far plane rather than infinity so the transmittance below cannot
		// evaluate 0 * INF and produce a NaN.
		const float toSurface = (toLight.y > 0.0001)
			? (waterHeight - samplePos.y) / toLight.y
			: GetCamera().z_far;

		return min(distanceToLight, clamp(toSurface, 0, GetCamera().z_far));
	}

	/**
	 * Point at which the light's path to a sample crosses into the water.
	 *
	 * Occlusion of a submerged sample by anything above the surface happens on
	 * the straight part of the path, so a shadow map that knows nothing about
	 * refraction is still exactly right when it is sampled here.
	 *
	 * @param[in] samplePos - Sample position, in world space.
	 * @param[in] toLight - Normalized direction from the sample to the light,
	 *                      refracted if the light is above the surface.
	 * @param[in] distanceToLight - Distance to the light, or `FLT_MAX` for a
	 *                              directional light.
	 *
	 * @return World position where the path enters the water.
	 */
	float3 SurfaceEntry(
		float3 samplePos,
		float3 toLight,
		float distanceToLight
	)
	{
		return samplePos
			+ toLight * SubmergedLightPath(samplePos, toLight, distanceToLight);
	}

	/**
	 * Light turned towards the eye at a sample, per metre of ray.
	 *
	 * @param[in] samplePos - Sample position, in world space.
	 * @param[in] toLight - Normalized direction from the sample to the light.
	 * @param[in] toEye - Normalized direction from the sample to the eye.
	 * @param[in] distanceToLight - Distance to the light, or `FLT_MAX` for a
	 *                              directional light.
	 *
	 * @return Per-channel scattered radiance per metre, before the light's own
	 *         colour, falloff and shadow term.
	 */
	float3 InScatter(
		float3 samplePos,
		float3 toLight,
		float3 toEye,
		float distanceToLight
	)
	{
		const float3 lightTransmittance = exp(
			-SubmergedLightPath(samplePos, toLight, distanceToLight) * sigmaT
		);

		// This HgPhase peaks at cosTheta = -1, so dot(toLight, toEye) puts the
		// peak where the eye looks along the light's own direction of travel -
		// forward scattering, which is where a forward-peaked medium belongs.
		const half phase = HgPhase(phaseG, (half)dot(toLight, toEye));

		return sigmaS * phase * lightTransmittance;
	}
};

/**
 * Builds the water medium as seen from a camera position.
 *
 * @param[in] eyePos - Camera position, in world space.
 *
 * @return The medium. Inactive when there is no ocean or the camera is above
 *         the surface, in which case the caller keeps its own fog.
 */
WaterVolumetrics GetWaterVolumetrics(float3 eyePos)
{
	const ShaderOcean ocean = GetWeather().ocean;

	WaterVolumetrics water;
	water.waterHeight = ocean.water_height;
	water.submersion = ocean.IsValid()
		? (half)saturate(
			(ocean.water_height - eyePos.y) / WATER_VOLUMETRICS_FADE_DEPTH)
		: (half)0;

	const float3 sigmaS = ocean.scattering.rgb;
	const float3 sigmaT = max(ocean.absorption.rgb + sigmaS, 0.00001);

	// Fade both coefficients together, so the medium degenerates to "not there"
	// rather than to "clear but still absorbing" - the latter would leave the
	// beams black on the way in rather than simply absent.
	water.sigmaS = sigmaS * water.submersion;
	water.sigmaT = sigmaT * water.submersion;

	// Reduce the asymmetry from the UNFADED ratio, so crossing the surface
	// changes how much the medium scatters but not the shape of its lobe.
	const float3 albedo = sigmaS / sigmaT;
	water.phaseG = (half)ocean.scattering.a * (1 - (half)saturate(albedo.g));

	return water;
}

#endif // WI_WATERVOLUMETRICS_HF
