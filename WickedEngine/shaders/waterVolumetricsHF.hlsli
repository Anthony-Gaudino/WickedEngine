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
 * `MarchWaterSunShafts` below is the older directional-only march, still used
 * by the underwater post pass until that adopts the volumetric light path.
 *
 * @note Sampling the shadow cascades also reads the *transparent* shadow layer,
 *       which is where the ocean surface writes its caustics (see the
 *       `SHADOWMAPRENDERING` branch of `oceanSurfacePS.hlsl`). The shafts
 *       therefore pick up caustic modulation without asking for it. Only the
 *       first couple of cascades render the ocean, so that modulation fades
 *       with distance.
 */

#include "globals.hlsli"
#include "lightingHF.hlsli"
#include "skyAtmosphere.hlsli"

/**
 * Optical depth past which the march stops.
 *
 * \f$e^{-5} \approx 0.007\f$, so beyond this the water has swallowed all but a
 * fraction of a percent and further steps buy nothing.
 */
static const float WATER_SHAFT_OPTICAL_DEPTH_CUTOFF = 5.0;

/**
 * Largest distance the dithered ray start may be offset by, in metres.
 *
 * The offset is a fraction of a step, so it grows with the step size. Left
 * unbounded, a long segment scatters neighbouring pixels far enough apart that
 * the dither's Bayer matrix becomes visible as a pattern on screen rather than
 * as noise.
 */
static const float WATER_SHAFT_MAX_DITHER_OFFSET = 10.0;

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

/**
 * Accumulates sunlight scattered towards the eye along a submerged ray.
 *
 * Single scattering: each step gathers sunlight that reaches the sample, turns
 * it towards the eye through the medium's phase function, and attenuates it
 * back along the ray already walked:
 * \[
 * L = \int_0^d T_{view}(s)\, \sigma_s\, p(\theta)\, V(s)\, T_{sun}(s)\, ds
 * \]
 * with \f$T\f$ Beer-Lambert transmittance and \f$V\f$ the shadow term.
 *
 * Example usage:
 * @code
 * const float3 shafts = MarchWaterSunShafts(
 *     campos, rayDir, waterPathLength,
 *     refractedSunDir, GetSunColor(),
 *     sigmaS, sigmaT, ocean.scattering.a,
 *     16, DTid.xy);
 * color.rgb += shafts;
 * @endcode
 *
 * @param[in] rayOrigin - Start of the submerged segment, in world space.
 * @param[in] rayDir - Normalized direction from the eye into the scene.
 * @param[in] rayLength - Length of the submerged segment (in metres).
 * @param[in] sunDir - Sun direction *after* refraction at the surface, so
 *                     pointing down into the water.
 * @param[in] sunColor - Sun radiance before the water attenuates it.
 * @param[in] sigmaS - Scattering coefficient of the water (1/m per channel).
 * @param[in] sigmaT - Extinction coefficient of the water (1/m per channel).
 * @param[in] phaseG - Henyey-Greenstein asymmetry of the medium.
 * @param[in] sampleCount - Steps along the segment. Must be at least 1.
 * @param[in] pixel - Dispatch pixel, used to dither the start and to sample
 *                    the shadow cascades.
 *
 * @return Radiance scattered towards the eye. Zero when there is no usable
 *         shadow casting sun, so the caller can add it unconditionally.
 */
float3 MarchWaterSunShafts(
	float3 rayOrigin,
	float3 rayDir,
	float rayLength,
	float3 sunDir,
	float3 sunColor,
	float3 sigmaS,
	float3 sigmaT,
	float phaseG,
	uint sampleCount,
	uint2 pixel
)
{
	// Reuses the engine's own "is this light usable for volumetrics" test,
	// which also hands back the last cascade worth walking.
	ShaderEntity light = (ShaderEntity)0;
	uint furthestCascade = 0;
	if (!furthest_cascade_volumetrics(light, furthestCascade))
		return 0;

	const float waterHeight = GetWeather().ocean.water_height;

	// Marching the full ray is pointless and actively harmful: where the depth
	// buffer holds the far plane the segment is kilometres long, and dividing
	// that into a handful of steps puts the samples so far apart that the
	// dithered start alone lands neighbouring pixels in completely different
	// places - which reads as the dither's Bayer pattern printed on screen.
	// Water swallows everything well before then, so stop once the view
	// transmittance has essentially died.
	const float opticalRange = WATER_SHAFT_OPTICAL_DEPTH_CUTOFF
		/ max(0.00001, min(sigmaT.r, min(sigmaT.g, sigmaT.b)));

	// A ray heading up leaves the water entirely, and everything past that
	// point is air which the submersion test below discards anyway. Stopping
	// there matters most exactly where the segment is shortest - a camera just
	// under the surface looking up - because otherwise the step stays long
	// while the water is thin, and the dither alone decides whether a pixel
	// samples that thin layer at all. Neighbouring pixels then disagree, which
	// reads as a halftone pattern.
	//
	// Continuous across the horizontal: as rayDir.y falls to zero the exit
	// distance grows without bound and the optical range takes over instead.
	float exitDistance = rayLength;
	[branch]
	if (rayDir.y > 0.0001)
	{
		exitDistance = max(0, (waterHeight - rayOrigin.y) / rayDir.y);
	}

	const float marchLength = min(min(rayLength, opticalRange), exitDistance);

	const float stepSize = marchLength / (float)max(1u, sampleCount);
	if (stepSize <= 0)
		return 0;

	// The sun is directional, so the scattering angle does not change along the
	// ray and the phase function is a constant.
	//
	// Reduce the asymmetry by how much of the light has been scattered, exactly
	// as the veiling term in underwaterCS does. Marine particulates peak some
	// 300x above isotropic, and summing that raw over a march produces a spike
	// at the sun bright enough to blow the image out; similarity theory says a
	// single scattering march standing in for a multiply scattering medium
	// should use the reduced value instead.
	//
	// The trailing uniform phase keeps this on the same scale as the rest of
	// the underwater pass, which applies one too (as the engine's height fog
	// does). Without it these shafts would sit some 4*PI above the veiling
	// light they are meant to modulate.
	const float3 albedo = sigmaS / sigmaT;
	const half multiScatter = (half)saturate(albedo.g);
	const half reducedG = (half)phaseG * (1 - multiScatter);
	const half cosTheta = (half)dot(rayDir, sunDir);
	const half phase = HgPhase(reducedG, cosTheta) * UniformPhase();

	// Sunlight crosses a slant path from the surface down to a sample, not just
	// the vertical drop. Refraction bounds the ray to the critical angle, so
	// that path can never be more than about 1.5x the depth.
	const float sunSlant = max(abs(sunDir.y), 0.66);

	// Marching in fixed steps bands badly at low sample counts; offsetting the
	// start per pixel turns the banding into noise, as the directional
	// volumetric light pass does.
	float3 P = rayOrigin + rayDir * min(
		stepSize * dither((min16uint2)pixel), WATER_SHAFT_MAX_DITHER_OFFSET
	);

	const float3 stepTransmittance = exp(-stepSize * sigmaT);
	float3 viewTransmittance = 1;
	float3 scattered = 0;

	[loop]
	for (uint i = 0; i < sampleCount; ++i)
	{
		// Walk out to the first cascade that contains the sample: the nearest
		// one covering it is also the highest resolution one.
		float3 shadow = 1;
		for (uint cascade = 0; cascade <= furthestCascade; ++cascade)
		{
			const float3 shadow_pos = mul(
				load_entitymatrix(light.GetMatrixIndex() + cascade),
				float4(P, 1)
			).xyz; // ortho matrix, no divide by .w
			const float3 shadow_uv = clipspace_to_uv(shadow_pos);

			[branch]
			if (is_saturated(shadow_uv))
			{
				shadow = shadow_2D(
					light, shadow_pos.z, shadow_uv.xy, cascade,
					(min16uint2)pixel
				);
				break;
			}
		}

		// Whatever light reaches this sample had to come down through the water
		// above it first, which is what makes the shafts fade out with depth
		// and lose their warm end on the way.
		const float sunDepth = waterHeight - P.y;
		const float3 sunTransmittance =
			exp(-(max(0, sunDepth) / sunSlant) * sigmaT);

		// There is no medium above the surface, so a ray that leaves the water
		// part way along must stop contributing rather than keep scattering in
		// air. Tested per sample instead of by shortening the segment, because
		// the distance at which the ray exits jumps discontinuously at the
		// horizontal and would draw a seam across the view.
		const float submerged = sunDepth >= 0 ? 1.0 : 0.0;

		scattered += submerged
			* viewTransmittance * shadow * sunTransmittance * sigmaS * phase * stepSize;

		viewTransmittance *= stepTransmittance;
		P += rayDir * stepSize;
	}

	return max(0, scattered * sunColor);
}

#endif // WI_WATERVOLUMETRICS_HF
