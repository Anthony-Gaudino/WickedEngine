#ifndef WI_WATERVOLUMETRICS_HF
#define WI_WATERVOLUMETRICS_HF

/**
 * The ocean's water as a participating medium, for lighting under the surface.
 *
 * Provides `WaterVolumetrics`, which answers the three questions a march has to
 * ask of a medium at each sample: how much of the light arriving there is
 * turned towards the eye, how much of it was lost on the way in, and how much
 * of what was turned survives the trip back out. The engine's volumetric light
 * passes ask those of the height fog; below the waterline they ask them here.
 *
 * Surface shading asks only the second of those - through
 * `WaterLightTransmittance` for a light it can name, and
 * `WaterAmbientTransmittance` for the daylight it cannot - so a submerged
 * surface is lit in the colour that actually reached it.
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
#include "underwaterHF.hlsli"

/**
 * Depth over which the water medium fades in below the surface (metres).
 *
 * Applies to SHADED POINTS only. A wall running down into the water would
 * otherwise pick up the full attenuation of a distant submerged light on the
 * first pixel below the waterline, drawing a hard line across it, and geometry
 * straddling the surface is common enough to be worth a metre of softening.
 *
 * Deliberately NOT applied to the eye - see GetWaterVolumetricsAtEye.
 */
static const float WATER_VOLUMETRICS_FADE_DEPTH = 1.0;

/**
 * Refractive index of water for visible light, relative to air.
 *
 * Matches the value the engine's other water shaders use - fresh water at green
 * wavelengths. Sea water is nearer 1.34 and the index climbs towards the blue
 * end, but neither difference survives being rendered.
 */
static const float WATER_REFRACTIVE_INDEX = 1.333;

/**
 * Height above the still plane within which the eye may still be submerged.
 *
 * A crest reaches above the still water plane, so an eye a little above it can
 * be inside the water even though the plane says otherwise. Only
 * `ocean_underwater_factor` can settle that, and it costs an unprojection and a
 * displacement map sample - so this is the threshold for bothering to ask.
 *
 * Deliberately loose. It has to clear both the wave displacement and the span
 * of the test plane that check sits on, a metre ahead of the near plane; being
 * generous costs one scalar compare, being tight would cut the fog off early.
 */
static const float WATER_EYE_SUBMERSION_MARGIN = 10.0;

/**
 * Cosine of the critical angle, measured from the vertical.
 *
 * \[
 * \sin\theta_c = \frac{1}{n}, \qquad
 * \cos\theta_c = \sqrt{1 - \frac{1}{n^2}} \approx 0.6612
 * \]
 * for \( n = 1.333 \), so about 48.6 degrees. Nothing that came through the
 * surface travels further from the vertical than this, in either direction: it
 * is the edge of Snell's window seen from below, and the steepest a ray leaving
 * the water can have been seen from above.
 */
static const float WATER_CRITICAL_ANGLE_COSINE = 0.6612;

/**
 * Optical path of downwelling daylight, as a multiple of depth.
 *
 * Refraction packs the whole sky into Snell's window, so light from above
 * reaches a submerged point from within a cone of about 48.6 degrees rather
 * than from the full hemisphere. Averaging the slant \( \sec\theta \) over that
 * cone under the cosine weight a diffuse surface applies anyway:
 * \[
 * \frac{\int_0^{\theta_c}\sec\theta\,\cos\theta\,\sin\theta\,d\theta}
 *      {\int_0^{\theta_c}\cos\theta\,\sin\theta\,d\theta}
 * = \frac{2\,(1 - \cos\theta_c)}{\sin^2\theta_c} \approx 1.2
 * \]
 * with \( \sin\theta_c = 1/1.333 \). So a point ten metres down sits under
 * roughly twelve metres of water, not ten.
 */
static const float WATER_DOWNWELLING_SLANT = 1.204;

/**
 * Optical depths a march covers before treating the water as opaque.
 *
 * The in-scatter itself no longer needs this - `StepIntegralWeight` makes the
 * march exact however coarsely the segment is divided. What still needs it is
 * everything sampled per step rather than integrated: the light's path to the
 * sample and the shadow lookup are held constant across a step, so a segment
 * left to run to the far plane would resolve the shafts at one sample per
 * hundred metres and lose their shape even while carrying the right energy.
 *
 * Three optical depths of the thinnest channel keeps the step in metres and
 * still captures 95% of that channel's in-scatter, essentially all of the
 * others'. Running further buys the last few percent of blue at the cost of a
 * coarser shadow.
 */
static const float WATER_MARCH_OPTICAL_DEPTHS = 3.0;

/**
 * Fraction of light reflected by the surface at straight-down incidence.
 *
 * \[
 * R_0 = \left(\frac{n - 1}{n + 1}\right)^2
 * \]
 * for \( n = 1.333 \), so a sun directly overhead loses about 2% at the
 * surface. It climbs to everything at grazing incidence.
 */
static const float WATER_FRESNEL_NORMAL_REFLECTANCE = 0.020373;

/**
 * Fraction of a uniform sky's light that gets through the surface.
 *
 * The cosine-weighted mean of the Fresnel transmittance over the hemisphere,
 * which is what the downwelling irradiance below a flat surface works out to:
 * \[
 * \bar{T} = 2\int_0^{\pi/2} T(\theta)\,\cos\theta\,\sin\theta\,d\theta
 * \]
 * evaluated against the exact unpolarized Fresnel equations. So the sky loses
 * about 7% at the surface however it is oriented - a constant, unlike the sun,
 * which arrives from one angle that changes through the day.
 */
static const float WATER_FRESNEL_DIFFUSE_TRANSMITTANCE = 0.9336;

/**
 * Fraction of the light arriving from above that crosses into the water.
 *
 * Some of what reaches the surface reflects off it instead of refracting down,
 * and how much depends steeply on the angle: about 2% straight down, still only
 * 6% at 60 degrees, but 58% at 85 and everything at the horizon. This is what
 * puts out a setting sun for anything below the surface, smoothly, rather than
 * leaving it lit until the moment it dips below the horizon.
 *
 * Parameterised by the SUBMERGED angle, which is what the shading code has to
 * hand, and converted back through Snell's law - exact for a directional light,
 * whose direction really was refracted, and an approximation for a local one,
 * which is treated as though it had been.
 *
 * Uses Schlick's approximation, which for this index ratio matches the exact
 * unpolarized Fresnel equations to better than a part in a thousand on the
 * hemisphere average.
 *
 * Example usage:
 * @code
 * light_color *= FresnelTransmittanceIntoWater(Lwater.y);
 * @endcode
 *
 * References:
 * https://en.wikipedia.org/wiki/Fresnel_equations
 *
 * @param[in] cosThetaSubmerged - Cosine of the angle between straight up and
 *                                the direction towards the light, taken below
 *                                the surface.
 *
 * @return Transmitted fraction, in [0, 1]. 0 for a light at or below the
 *         horizon, which delivers nothing through the surface at all.
 */
float FresnelTransmittanceIntoWater(float cosThetaSubmerged)
{
	// Snell's law back to the air side. The reflectance is governed by the
	// angle the light struck the surface at, not by the steeper one it travels
	// at afterwards, and the two differ by a lot - the whole hemisphere above
	// maps into a cone of 48.6 degrees below.
	const float cosBelow = saturate(cosThetaSubmerged);
	const float sinBelow = sqrt(saturate(1 - cosBelow * cosBelow));
	const float sinAbove = saturate(sinBelow * WATER_REFRACTIVE_INDEX);
	const float cosAbove = sqrt(saturate(1 - sinAbove * sinAbove));

	const float grazing = pow5(1 - cosAbove);

	return (1 - WATER_FRESNEL_NORMAL_REFLECTANCE) * (1 - grazing);
}

/**
 * Bends a direction pointing at an above-water light down into the water.
 *
 * Snell's law at the flat surface. Seen from below, a source above the water is
 * pulled towards the vertical - the same compression that fits the whole sky
 * inside Snell's window - so its light both arrives from a steeper angle and
 * crosses less water than the straight line towards it suggests.
 *
 * Exact for a directional light, where the interface is a plane and the source
 * is at infinity, so one refraction serves the whole scene. A local light would
 * need a per-sample solution and does not get one here; near-overhead entry is
 * close to straight anyway, and the path is clipped at the surface either way.
 *
 * Example usage:
 * @code
 * const float3 Lwater = RefractIntoWater(GetSunDirection());
 * @endcode
 *
 * References:
 * https://en.wikipedia.org/wiki/Snell%27s_law
 *
 * @param[in] toLight - Normalized direction from the sample to the light.
 *
 * @return The refracted direction, still pointing towards the light. Returned
 *         unchanged for a light at or below the horizon, which has no downward
 *         crossing to refract.
 */
float3 RefractIntoWater(float3 toLight)
{
	[branch]
	if (toLight.y <= 0)
	{
		return toLight;
	}

	// refract() wants the direction of travel and a normal facing into it, so
	// the ray is flipped to head downwards and flipped back on the way out.
	// This enters the denser medium, so there is no internal reflection case.
	return -refract(-toLight, float3(0, 1, 0), 1.0 / WATER_REFRACTIVE_INDEX);
}

/**
 * Length of the in-water leg of a view ray that left the water for the eye.
 *
 * Light reaching an eye in the air bent at the surface on its way out, so the
 * water it crossed is not the straight line the eye looks along. Reciprocity
 * gives that leg without tracing anything: the direction back towards the eye,
 * bent into the water, is the direction the light travelled up through it, and
 * the leg is the depth over that direction's vertical cosine:
 * \[
 * s = \frac{d}{\cos\theta_{water}}
 * \]
 * Taken against the flat plane rather than the wave normal, so the slant stays
 * steady instead of chattering with the surface detail.
 *
 * **The result is bounded**, which is the whole reason to prefer it to the
 * geometric segment. Refraction confines the leg to Snell's window, so it lies
 * between 1.0 and about 1.512 times the depth however far away the eye stands,
 * while the straight segment grows without limit at a grazing angle and would
 * extinguish water that is plainly see-through. Flooring the cosine at the
 * critical angle is that same analytic bound, so the clamp only ever catches a
 * degenerate direction.
 *
 * Example usage:
 * @code
 * const float path = SubmergedViewPath(waterDepth, surface.V);
 * const float3 transmittance = exp(-path * medium.sigmaT);
 * @endcode
 *
 * @param[in] depth - Depth of the far end below the surface (in metres).
 * @param[in] toEye - Normalized direction from the far end towards the eye.
 *
 * @return Length of the segment lying in water (in metres).
 *
 * @note For an eye **in the air** only. Looking out from beneath the surface
 *       nothing refracts on the way to the eye, and the straight segment is
 *       already the right answer.
 */
float SubmergedViewPath(float depth, float3 toEye)
{
	return depth / max(RefractIntoWater(toEye).y, WATER_CRITICAL_ANGLE_COSINE);
}

/**
 * The ocean's water medium, evaluated at one world position.
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
 * const WaterVolumetrics water = GetWaterVolumetricsAtEye(ScreenCoord);
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
	 * Henyey-Greenstein asymmetry as authored, before any reduction.
	 *
	 * Marine particulates scatter extremely far forward - g around 0.9, see
	 * WaterMedium::PhaseAsymmetry - but that is the phase of a SINGLE
	 * scattering event, and feeding it straight into a term that stands for
	 * many puts a spike some 300x the isotropic value at the sun. Use one of
	 * the reduced forms below instead; this is here to derive them from.
	 */
	half phaseG0;

	/**
	 * Henyey-Greenstein asymmetry reduced by similarity theory.
	 *
	 * Reduced by the single-scattering albedo, which is the limit the
	 * path-dependent reduction tends to as the path grows. The right choice
	 * for a MARCH: it integrates a segment sample by sample and so has no one
	 * path length to reduce by.
	 */
	half phaseG;

	/** World height of the still water plane. */
	float waterHeight;

	/** 0 at or above the surface, reaching 1 once fully submerged. */
	half submersion;

	/**
	 * Whether this medium applies at all.
	 *
	 * @return true when the position it was built for is below the surface.
	 */
	bool IsActive()
	{
		return submersion > 0;
	}

	/**
	 * Distance past which nothing more can reach the eye through this water.
	 *
	 * Taken on the thinnest channel - whichever penetrates furthest, blue in
	 * most water - so the range is the one that still has something left to
	 * gather rather than the one that gave up first.
	 *
	 * Exists because the in-scatter march INTEGRATES over a fixed number of
	 * samples: its accuracy depends on the step, and the step depends on how
	 * far the segment runs. See the caller for what that costs.
	 *
	 * @return Distance in metres. Meaningless unless `IsActive()`, since the
	 *         coefficients fall to nothing above the surface.
	 */
	float VisibleRange()
	{
		return WATER_MARCH_OPTICAL_DEPTHS / min3(sigmaT.r, sigmaT.g, sigmaT.b);
	}

	/**
	 * Asymmetry reduced by how much of the light was actually scattered.
	 *
	 * For a single stretch of water of known length, where the veiling term is
	 * to hand. It already combines "how much was extinguished" with "how much
	 * of that was redirected rather than absorbed", which is exactly what the
	 * similarity reduction asks for - so clear water keeps a sharp sun glow
	 * while turbid water, having scattered far more, spreads it into an even
	 * haze.
	 *
	 * Tends to `phaseG` as the path grows and the veiling term saturates at
	 * the albedo, so the two agree in the limit and differ only where there
	 * genuinely is a path length to be had.
	 *
	 * @param[in] inscatterAmount - Veiling term, albedo * (1 - transmittance).
	 *
	 * @return The reduced asymmetry, for `HgPhase`.
	 */
	half ReducedPhaseG(float3 inscatterAmount)
	{
		return phaseG0 * (1 - (half)saturate(inscatterAmount.g));
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
	 * Integral of the view transmittance across one march step.
	 *
	 * A step contributes the in-scatter along a whole segment rather than at a
	 * point, and the view transmittance falls exponentially across it:
	 * \[
	 * \int_{d-\Delta}^{d} e^{-\sigma_t s}\,ds
	 *   = e^{-\sigma_t (d - \Delta)}\,
	 *     \frac{1 - e^{-\sigma_t \Delta}}{\sigma_t}
	 * \]
	 * so a step weighs this much, not its own length. Sampling the
	 * transmittance at one end and multiplying by the length is a left Riemann
	 * sum over an exponential: it misses whatever falls inside the first step,
	 * and that is where nearly all the light reaching the eye comes from. Red
	 * loses about half its in-scatter at a two metre step, being extinguished
	 * several times faster than blue.
	 *
	 * Summed over the march this telescopes back to the closed form, so the
	 * result is exact at any step size and in every channel - the accuracy no
	 * longer depends on how finely the segment happens to be divided.
	 *
	 * Depends only on the step, so callers hoist it out of their loop. Pair it
	 * with `ViewTransmittance` evaluated at the step's NEAR end, which is what
	 * the formula above is anchored to.
	 *
	 * @param[in] stepSize - Length of one march step (in metres).
	 *
	 * @return Per-channel weight for the step's in-scatter (in metres). Tends
	 *         to `stepSize` itself as the step shrinks.
	 */
	float3 StepIntegralWeight(float stepSize)
	{
		return (1 - exp(-stepSize * sigmaT)) / sigmaT;
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
	 * What survives of a light's journey to a sample: surface, then water.
	 *
	 * Water takes red out of a beam roughly an order of magnitude faster than
	 * blue, so this is what leaves a submerged red lamp red within arm's reach
	 * and blue a few metres out. A light that started in the air loses a second
	 * share of itself at the surface, to reflection.
	 *
	 * @param[in] samplePos - Sample position, in world space.
	 * @param[in] toLight - Normalized direction from the sample to the light.
	 * @param[in] distanceToLight - Distance to the light, or `FLT_MAX` for a
	 *                              directional light.
	 *
	 * @return Per-channel surviving fraction, in [0, 1].
	 */
	float3 LightTransmittance(
		float3 samplePos,
		float3 toLight,
		float distanceToLight
	)
	{
		const float path =
			SubmergedLightPath(samplePos, toLight, distanceToLight);

		// A path that had to be clipped is one that reached the surface before
		// it reached the light, so the light is in the air and part of it
		// reflected away rather than coming down. A submerged light crosses
		// nothing and keeps all of it.
		const float crossing = (path < distanceToLight)
			? FresnelTransmittanceIntoWater(toLight.y)
			: 1;

		return crossing * exp(-path * sigmaT);
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
		const float3 lightTransmittance =
			LightTransmittance(samplePos, toLight, distanceToLight);

		// This HgPhase peaks at cosTheta = -1, so dot(toLight, toEye) puts the
		// peak where the eye looks along the light's own direction of travel -
		// forward scattering, which is where a forward-peaked medium belongs.
		const half phase = HgPhase(phaseG, (half)dot(toLight, toEye));

		return sigmaS * phase * lightTransmittance;
	}
};

/**
 * Builds the water medium for a given degree of submersion.
 *
 * Shared by the two entry points below, which differ only in how they decide
 * that number.
 *
 * @param[in] submersion - 0 for no medium at all, 1 for the full authored
 *                         water, in between to fade it in.
 *
 * @return The medium. Inactive when there is no ocean, whatever was asked for.
 */
WaterVolumetrics MakeWaterVolumetrics(half submersion)
{
	const ShaderOcean ocean = GetWeather().ocean;

	WaterVolumetrics water;
	water.waterHeight = ocean.water_height;
	water.submersion = ocean.IsValid() ? submersion : (half)0;

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
	water.phaseG0 = (half)ocean.scattering.a;
	water.phaseG = water.phaseG0 * (1 - (half)saturate(albedo.g));

	return water;
}

/**
 * Builds the water medium as it applies at a shaded point.
 *
 * Fades in over the first metre of depth, because a point's own submersion is
 * a SPATIAL quantity: a wall running down into the water would otherwise take
 * the full attenuation of a distant submerged light on the first pixel below
 * the waterline, drawing a hard line across it.
 *
 * @param[in] position - World position of the shaded point.
 *
 * @return The medium. Inactive when there is no ocean or the point is at or
 *         above the surface, in which case nothing attenuates.
 */
WaterVolumetrics GetWaterVolumetrics(float3 position)
{
	const float waterHeight = GetWeather().ocean.water_height;

	return MakeWaterVolumetrics((half)saturate(
		(waterHeight - position.y) / WATER_VOLUMETRICS_FADE_DEPTH));
}

/**
 * Builds the water medium as it applies to one pixel's view ray.
 *
 * Per PIXEL, not per camera, and using the very test the underwater post pass
 * and the submerged depth of field and chromatic aberration already share. A
 * camera crossing the surface is submerged across part of the frame long before
 * its own centre goes under - the waterline sweeps down the screen - and asking
 * where the camera is instead answered for the whole frame at once. That left
 * the post pass painting water over the bottom of the view while this medium
 * was still inactive, so the fog appeared but the shafts did not, and they came
 * in together only once the camera itself was a good way under.
 *
 * Judged against the wave displaced surface, so a crest counts as water, and at
 * the same plane a metre ahead of the near plane, so all of them switch on the
 * same pixels at the same instant.
 *
 * There is no fade with depth here, unlike a shaded point. Fading the medium
 * thins the water rather than removing it, and both halves of that are wrong:
 * scaling sigmaS costs the in-scatter over any path bounded by geometry - a sea
 * bed forty metres out gathers under a third of what it should at ten
 * centimetres down - while scaling sigmaT makes the water CLEARER the nearer
 * the surface the eye is, seen as looking far too far while going under.
 *
 * @param[in] screenUV - Screen space UV coordinates (0-1) of the pixel.
 *
 * @return The medium. Inactive where there is no ocean or the pixel is above
 *         the waterline, in which case the caller keeps its own fog.
 */
WaterVolumetrics GetWaterVolumetricsAtEye(float2 screenUV)
{
	return MakeWaterVolumetrics((half)ocean_underwater_factor(screenUV));
}

/**
 * Transmittance of the water lying between a lit point and one of its lights.
 *
 * The surface-shading counterpart to the marches above, and the missing leg of
 * underwater light transport: the water a light crossed on its way in dims and
 * colours it exactly as the water between the lit point and the eye does. Every
 * light needs it, not only the sun. Without it a lamp on the sea bed lights its
 * surroundings in the colour it was authored and only the view ray tints the
 * result, so a red lamp reads blue everywhere - including right beside it.
 *
 * Attenuation follows PATH LENGTH, which is why no depth-based tint can stand
 * in for it: a lamp a metre away across the sea bed and the sun twenty metres
 * overhead are at the same depth and nothing alike.
 *
 * Free in a scene with no ocean, where the medium is inactive and the branch is
 * uniform.
 *
 * Example usage:
 * @code
 * light_color *= WaterLightTransmittance(surface.P, L, dist);
 * @endcode
 *
 * @param[in] samplePos - The lit point, in world space.
 * @param[in] toLight - Normalized direction from the lit point to the light.
 * @param[in] distanceToLight - Distance to the light, or `FLT_MAX` for a
 *                              directional light.
 *
 * @return Per-channel surviving fraction, in [0, 1]. 1 at or above the surface,
 *         so callers need no test of their own.
 */
half3 WaterLightTransmittance(
	float3 samplePos,
	float3 toLight,
	float distanceToLight
)
{
	const WaterVolumetrics water = GetWaterVolumetrics(samplePos);

	[branch]
	if (!water.IsActive())
	{
		return 1;
	}

	return (half3)water.LightTransmittance(samplePos, toLight, distanceToLight);
}

/**
 * Transmittance of the water above a point, for downwelling daylight.
 *
 * The ambient term and the sky probes stand in for light arriving from the
 * whole upper hemisphere, and none of it can reach a submerged point without
 * coming down through the water first. That is the one part of the water's
 * effect which genuinely is a function of DEPTH rather than of a path to some
 * particular light, because the source is the sky itself.
 *
 * Without it a sea bed twenty metres down is lit by an ambient term that never
 * got wet, which is enough on its own to hide the water's colour: at that depth
 * daylight is already down to a ten-thousandth of its red and a fifth of its
 * blue, so any unattenuated white washes the blue straight back out. Diffuse
 * ambient also enters the shading without the 1/pi that direct light carries,
 * which makes it the larger of the two long before the water alone would.
 *
 * Includes the share the surface reflects away before any of it gets in, which
 * for a whole sky is a constant rather than the steep function of angle a
 * single light sees.
 *
 * Applies to reflected daylight as well as to the diffuse term, so it is only
 * an approximation for a reflection of something nearby and submerged, whose
 * light did not come down from the sky at all. The sky probe dominates.
 *
 * Example usage:
 * @code
 * lighting.indirect.diffuse *= WaterAmbientTransmittance(surface.P);
 * @endcode
 *
 * @param[in] samplePos - The lit point, in world space.
 *
 * @return Per-channel surviving fraction, in [0, 1]. 1 at or above the surface,
 *         so callers need no test of their own.
 */
half3 WaterAmbientTransmittance(float3 samplePos)
{
	const WaterVolumetrics water = GetWaterVolumetrics(samplePos);

	[branch]
	if (!water.IsActive())
	{
		return 1;
	}

	const float depth = max(0, water.waterHeight - samplePos.y);

	return (half3)(
		WATER_FRESNEL_DIFFUSE_TRANSMITTANCE
		* exp(-depth * WATER_DOWNWELLING_SLANT * water.sigmaT)
	);
}

#endif // WI_WATERVOLUMETRICS_HF
