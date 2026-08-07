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
 * Height above the wave surface within which the eye may still be submerged.
 *
 * The eye's own height settles which side of the surface its CENTRE is on, but
 * not which side each PIXEL is on: `ocean_underwater_factor` judges that on a
 * plane a metre ahead of the near plane, which spans a good deal of world at a
 * wide field of view, and that sweep across the screen is the whole point of
 * it. So this is the threshold for bothering to ask - the span of that test
 * plane, not the waves, which the surface height has already accounted for.
 *
 * Deliberately loose: being generous costs one scalar compare, being tight
 * would cut the fog off while part of the frame is still under.
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
 * Raman scattering coefficient of pure water, by the channel it EMITS into.
 *
 * Water scatters a small share of what passes through it **inelastically**: the
 * photon gives up a fixed quantum of energy to an O-H stretch vibration of the
 * molecule and carries on at a longer wavelength. The shift is fixed in
 * wavenumber rather than in wavelength - about \f$3400\,\mathrm{cm}^{-1}\f$ -
 * which across these three reference bands works out to very nearly one band
 * each time: 460 nm re-emerges at 545 nm and 540 nm at 661 nm. So red is lit by
 * green, green by blue, and blue by an ultraviolet that never reaches the water
 * in any quantity worth counting - which is why the blue component here is
 * zero, and why nothing converts twice.
 *
 * Taken from the measured \f$b_R(488) = 2.7\times10^{-4}\,\mathrm{m}^{-1}\f$
 * and carried to each excitation band by the \f$\lambda^{-5.3}\f$ law.
 *
 * **A property of the water molecule**, so unlike every other coefficient the
 * medium carries it does not move with what is dissolved or suspended in the
 * water - particles do not Raman scatter appreciably. It is a constant here
 * rather than something the scene publishes.
 *
 * References:
 * Bartlett et al. 1998, *Raman scattering by pure water and seawater*;
 * Sathyendranath & Platt 1998, *Ocean-color model incorporating transspectral
 * processes*.
 */
static const float3 WATER_RAMAN_SCATTERING = float3(1.6e-4, 3.7e-4, 0.0);

/**
 * Share of the red channel's response that light at 685 nm produces.
 *
 * Chlorophyll fluoresces in a narrow line at 685 nm, and this render has no
 * channel there - it samples red at 620 nm. Handing a 685 nm line straight to
 * the red channel would claim it looks as red as 620 nm light of the same
 * radiance, and it does not: 685 nm is far out on the tail of both the eye's
 * response and the display's red primary, where a watt buys around a
 * thirteenth of the response it buys at 620 nm.
 * \[
 * \frac{R(685)}{R(620)}, \qquad
 * R(\lambda) \propto 3.2406\,\bar{x} - 1.5372\,\bar{y} - 0.4986\,\bar{z}
 * \]
 * evaluated on the CIE 1931 observer through the sRGB matrix.
 *
 * A statement about the primaries rather than about the water, which is why it
 * lives here and not in `WaterMedium::Fluorescence`.
 *
 * References:
 * https://en.wikipedia.org/wiki/CIE_1931_color_space
 */
static const float WATER_FLUORESCENCE_BAND_RESPONSE = 0.079;

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
	 * Backscattering coefficient in 1/m per channel.
	 *
	 * What comes back rather than carrying on forwards. See `EmergentAlbedo`
	 * for why `sigmaS` cannot answer that.
	 */
	float3 sigmaB;

	/**
	 * Share of the light one band up that comes back out at this one.
	 *
	 * The Raman counterpart to `EmergentAlbedo`, and dimensionless like it:
	 * where that returns a share of the light at its own wavelength, this
	 * returns a share of the light at the shorter one that excited it. Carried
	 * as a member rather than derived on demand because it has to be built
	 * from the UNFADED coefficients - see `MakeWaterVolumetrics`.
	 */
	float3 ramanAlbedo;

	/**
	 * Share of the light in each band that comes back out as red.
	 *
	 * Indexed by the band that EXCITED it, not the one it emerges in - every
	 * band chlorophyll absorbs feeds one line at 685 nm, so this gathers into
	 * red where `ramanAlbedo` shifts channel by channel. Carries the display
	 * response at 685 nm as well as the water's own conversion, so it is ready
	 * to be dotted with the downwelling and put straight into red.
	 */
	float3 fluorescenceAlbedo;

	/**
	 * Single-scattering albedo on the green channel, as authored.
	 *
	 * How much of what the medium extinguishes it merely redirects, which is
	 * to say how many times a photon is likely to scatter before it is
	 * absorbed. Governs how thoroughly multiple scattering has smoothed the
	 * upwelling light out - see `BidirectionalFactor`. Unfaded, for the reason
	 * `phaseG` is.
	 */
	half singleScatterAlbedo;

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

	/**
	 * World height of the water surface this medium was built for.
	 *
	 * The WAVE surface where the caller had a position to measure it over, and
	 * the still plane where it did not - see the two builders below for which
	 * is which and why. Everything that needs to know how far a point is under
	 * the water reads this rather than the plane, so a crest counts as the
	 * water it is.
	 */
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
	 * Colour the medium settles on once light has scattered many times in it.
	 *
	 * **Scattering is not what decides the colour of water - backscattering
	 * is.** Only light turned back towards the viewer is ever seen, and how
	 * much that is varies by two orders of magnitude between constituents.
	 * Silt and phytoplankton cells are large compared to the wavelength and
	 * scatter overwhelmingly forwards, returning under two percent and half a
	 * percent of what they scatter; pure water's Rayleigh-like scattering is
	 * nearly symmetric and returns half. So the enormous, nearly grey
	 * \f$\sigma_s\f$ of a turbid water counts for far less than its size
	 * suggests, and the water's own sharply selective absorption survives.
	 *
	 * The transport (similarity) approximation supplies the equivalent
	 * isotropic scattering: an isotropic scatterer returns half of what it
	 * scatters, so \f$2 b_b\f$ is the coefficient that sends back the same
	 * amount. Feeding its albedo to the two-stream reflectance of a
	 * semi-infinite medium,
	 * \f[
	 * R_\infty = \frac{1 - \sqrt{1 - \omega}}{1 + \sqrt{1 - \omega}},
	 * \qquad \omega = \frac{2 b_b}{\sigma_a + 2 b_b}
	 * \f]
	 * lands where measurements do: about 5% for the murkiest water and under
	 * 2% for the clearest.
	 *
	 * @note Carries no directionality, so callers must not weight it with a
	 *       phase function. Light that has scattered this many times has no
	 *       memory of the direction it arrived from; the phase function
	 *       belongs to the single-scattering term alongside it.
	 *
	 * References:
	 * Gordon et al. 1988, *A semianalytic radiance model of ocean color*;
	 * https://en.wikipedia.org/wiki/Kubelka%E2%80%93Munk_theory
	 *
	 * @return Emergent albedo per channel, 0 for a purely absorbing medium.
	 */
	float3 EmergentAlbedo()
	{
		const float3 reducedScattering = 2 * sigmaB;
		const float3 absorption = max(sigmaT - sigmaS, 0);

		const float3 albedo = reducedScattering
			/ max(absorption + reducedScattering, 1e-6);
		const float3 root = sqrt(saturate(1 - albedo));

		return (1 - root) / (1 + root);
	}

	/**
	 * How the emergent colour reshapes with the sun and the view direction.
	 *
	 * `EmergentAlbedo` is one number for a whole body of water, which quietly
	 * assumes what comes back out is spread evenly over the sky. It is not.
	 * A low sun drives its light in at a slant, giving it further to travel
	 * through the scattering layer and so more chance to be turned back, and
	 * light leaving the water is somewhat peaked about the vertical, so a
	 * slanted view collects less of it than a view straight down. Ocean optics
	 * carries these as the pair \f$f/Q\f$ - the first factor for the sun, the
	 * second for the view.
	 *
	 * Taken from the single-scattering reflectance of a semi-infinite medium,
	 * \f[
	 * R(\mu_s, \mu_v) \propto \frac{1}{\mu_s + \mu_v}
	 * \f]
	 * so what a viewer collects goes as \f$\mu_v / (\mu_s + \mu_v)\f$,
	 * normalized here to the overhead sun seen from straight above.
	 *
	 * **Refraction is what keeps this small.** Both directions are measured
	 * inside the water, and nothing crossing the surface bends further from
	 * the vertical than the critical angle, so neither cosine can fall below
	 * about 0.66 however low the sun sits or however grazing the view. The
	 * whole span is therefore only about \f$\pm 20\%\f$, where the published
	 * tables range much wider over in-water directions whose light could never
	 * have reached an eye in the air to begin with.
	 *
	 * Faded towards no reshaping at all as the single-scattering albedo rises,
	 * by the similarity argument `phaseG` already uses: the angular law above
	 * describes ONE scattering event, and water that scatters many times over
	 * before absorbing anything has smoothed most of the direction back out.
	 *
	 * @note Applies to the ELASTIC return only. The inelastic terms beside it
	 *       are isotropic emission - a molecule or a cell radiates the same in
	 *       every direction - so they carry no angular preference to reshape.
	 *
	 * @note Derived for light escaping upward out of deep water, which is the
	 *       geometry an eye in the air sees. A submerged eye looking level has
	 *       no such column in front of it; clamping at the critical angle
	 *       bounds what that case can do rather than describing it.
	 *
	 * References:
	 * Morel & Gentili 1996, *Diffuse reflectance of oceanic waters III*;
	 * Chandrasekhar 1960, *Radiative Transfer*.
	 *
	 * @param[in] cosSun - Cosine between straight up and the direction towards
	 *                     the sun, measured below the surface.
	 * @param[in] cosView - Cosine between straight up and the direction towards
	 *                      the eye, measured below the surface.
	 *
	 * @return Factor on the emergent albedo, in roughly [0.8, 1.2] and 1 at the
	 *         reference geometry. Tends to 1 as the water turns turbid.
	 */
	float BidirectionalFactor(float cosSun, float cosView)
	{
		// The floor is the analytic bound rather than a guard, so it joins up
		// smoothly: a ray arriving along the surface refracts to exactly the
		// critical angle, which is the floor itself. Replace it with a smaller
		// number and a camera looking level under water gets a hard band across
		// the screen where the two sides stop meeting.
		const float muSun = max(cosSun, WATER_CRITICAL_ANGLE_COSINE);
		const float muView = max(cosView, WATER_CRITICAL_ANGLE_COSINE);

		const float singleScattering = 2 * muView / (muSun + muView);

		return lerp(
			1, singleScattering, 1 - (float)saturate(singleScatterAlbedo));
	}

	/**
	 * Light the water gives back at a longer wavelength than it took in.
	 *
	 * Every other term here is elastic - a photon changes direction and keeps
	 * its colour. A small share instead loses a fixed quantum of energy to the
	 * water molecule and comes back **redder than it arrived**, which makes the
	 * water a source of light at wavelengths nothing in it reflected.
	 *
	 * That is what stops the clearest water being as purely blue as its
	 * absorption says it should be. A few metres down blue is very nearly all
	 * that is left, and a hundredth of a percent of a large number beats the
	 * elastic return of a small one: green comes out around a quarter brighter
	 * than backscatter alone would leave it, and past ten metres or so what red
	 * there is is mostly this rather than red that made it down. Turbid water
	 * drowns it - the elastic return is two orders of magnitude larger there and
	 * the blue driving this is gone - so it costs nothing to leave switched on.
	 *
	 * \f[
	 * R_{Raman}(\lambda_{em}) = \frac{b_R(\lambda_{ex})}
	 *   {2\,\left(K(\lambda_{ex}) + K(\lambda_{em})\right)}
	 * \f]
	 * from integrating an isotropic volume source down a semi-infinite column:
	 * half of what converts at each depth heads back up, having been attenuated
	 * on the way down at the wavelength that excited it and on the way back at
	 * the one it now carries.
	 *
	 * @note Carries no directionality, for the same reason `EmergentAlbedo`
	 *       does not: the emission is isotropic and remembers nothing of where
	 *       the exciting photon came from. Do not weight it with a phase
	 *       function.
	 *
	 * References:
	 * Sathyendranath & Platt 1998, *Ocean-color model incorporating
	 * transspectral processes*.
	 *
	 * @param[in] downwelling - Daylight reaching this depth, per channel, in
	 *                          whatever units the result is wanted in.
	 *
	 * @return Per-channel radiance the water re-emits, to be weighted by how
	 *         much of the segment it fills exactly as `EmergentAlbedo` is.
	 */
	float3 RamanEmission(float3 downwelling)
	{
		// Each channel is lit by the one above it in energy: red by green,
		// green by blue, and blue by an ultraviolet there is no channel for -
		// which the zero coefficient in the albedo takes care of.
		return ramanAlbedo * float3(downwelling.g, downwelling.b, 0);
	}

	/**
	 * Light the phytoplankton in the water give back after absorbing it.
	 *
	 * A cell that catches more light than its photosynthesis can consume sheds
	 * the excess a few nanoseconds later as a narrow line at **685 nm**. Like
	 * Raman this is inelastic and so puts light back at a wavelength that may
	 * have been extinguished long before it could reach anything; unlike Raman
	 * it is a property of what is living in the water rather than of the water
	 * itself, so it rises with the algae and vanishes in a sterile sea.
	 *
	 * Every band the chlorophyll absorbs feeds the same line - the cell relaxes
	 * to one state before emitting, however it was excited - so this **gathers
	 * into red** rather than shifting band by band as `RamanEmission` does.
	 * Blue does most of the exciting, chlorophyll's absorption peaking there.
	 *
	 * @note Radiometrically this is the larger of the two inelastic terms in
	 *       any water with much life in it, several times Raman's by the
	 *       turbid end. Almost none of it survives being shown: the line sits
	 *       far enough into the red that the display's primary answers it with
	 *       about a thirteenth of the response 620 nm would get, which
	 *       `WATER_FLUORESCENCE_BAND_RESPONSE` accounts for and which leaves
	 *       the visible result the smaller of the two.
	 *
	 * @note Carries no directionality, exactly as `EmergentAlbedo` and
	 *       `RamanEmission` do not - the emission is isotropic and has
	 *       forgotten where the exciting photon came from.
	 *
	 * References:
	 * Babin et al. 1996, *Remote sensing of sea surface sun-induced chlorophyll
	 * fluorescence*.
	 *
	 * @param[in] downwelling - Daylight reaching this depth, per channel, in
	 *                          whatever units the result is wanted in.
	 *
	 * @return Per-channel radiance the water re-emits, red only, to be weighted
	 *         by how much of the segment it fills exactly as `EmergentAlbedo`
	 *         is.
	 */
	float3 FluorescenceEmission(float3 downwelling)
	{
		return float3(dot(fluorescenceAlbedo, downwelling), 0, 0);
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
 * Shared by the entry points below, which differ only in how they decide that
 * number and which surface they measure it against.
 *
 * @param[in] submersion - 0 for no medium at all, 1 for the full authored
 *                         water, in between to fade it in.
 * @param[in] waterHeight - World height of the surface this medium sits under.
 *
 * @return The medium. Inactive when there is no ocean, whatever was asked for.
 */
WaterVolumetrics MakeWaterVolumetrics(half submersion, float waterHeight)
{
	const ShaderOcean ocean = GetWeather().ocean;

	WaterVolumetrics water;
	water.waterHeight = waterHeight;
	water.submersion = ocean.IsValid() ? submersion : (half)0;

	const float3 sigmaS = ocean.scattering.rgb;
	const float3 sigmaB = ocean.backscattering.rgb;
	const float3 sigmaT = max(ocean.absorption.rgb + sigmaS, 0.00001);

	// Fade both coefficients together, so the medium degenerates to "not there"
	// rather than to "clear but still absorbing" - the latter would leave the
	// beams black on the way in rather than simply absent.
	water.sigmaS = sigmaS * water.submersion;
	water.sigmaT = sigmaT * water.submersion;
	water.sigmaB = sigmaB * water.submersion;

	// Reduce the asymmetry from the UNFADED ratio, so crossing the surface
	// changes how much the medium scatters but not the shape of its lobe.
	const float3 albedo = sigmaS / sigmaT;
	water.singleScatterAlbedo = (half)saturate(albedo.g);
	water.phaseG0 = (half)ocean.scattering.a;
	water.phaseG = water.phaseG0 * (1 - water.singleScatterAlbedo);

	// Unfaded for a stronger reason than the asymmetry: this RISES as the
	// medium thins, since a longer mean free path converts more of what crosses
	// it, while the caller's weighting falls away. Faded, the two would cancel
	// and the water would go on emitting as it left the frame.
	//
	// The same reduced extinction EmergentAlbedo works from, so both describe
	// one medium. Each channel is paired with the extinction of the band above
	// it, which is where its light was taken from; blue has no band above it
	// and a zero coefficient, so what it pairs with never shows.
	const float3 reducedExtinction = max(sigmaT - sigmaS, 0) + 2 * sigmaB;
	water.ramanAlbedo = 0.5 * WATER_RAMAN_SCATTERING
		/ max(reducedExtinction + reducedExtinction.gbb, 0.00001);

	// The same integral, unfaded for the same reason, but every band escapes
	// against red's extinction rather than its own: whatever excited the cell,
	// what leaves it is one line at 685 nm.
	water.fluorescenceAlbedo =
		0.5 * WATER_FLUORESCENCE_BAND_RESPONSE * ocean.fluorescence
		/ max(reducedExtinction + reducedExtinction.r, 0.00001);

	return water;
}

/**
 * Builds the water medium for a given submersion, under the still plane.
 *
 * For callers that have already decided the medium is fully present and have
 * no position to measure a surface over - the ocean surface compositing its own
 * reflection, or a wave transmitting light through itself.
 *
 * @param[in] submersion - 0 for no medium at all, 1 for the full authored
 *                         water, in between to fade it in.
 *
 * @return The medium. Inactive when there is no ocean, whatever was asked for.
 */
WaterVolumetrics MakeWaterVolumetrics(half submersion)
{
	return MakeWaterVolumetrics(submersion, GetWeather().ocean.water_height);
}

/**
 * Builds the water medium as it applies at a shaded point.
 *
 * **Measured against the wave surface over the point**, so a shaded point
 * inside a crest is in the water. Against the still plane it was not: every
 * surface around a camera sitting in a crest reads as above the water, the
 * medium comes back inactive, and the whole scene is lit as though in air -
 * however much water the fog says is in front of it.
 *
 * Fades in over the first metre of depth, because a point's own submersion is
 * a SPATIAL quantity: a wall running down into the water would otherwise take
 * the full attenuation of a distant submerged light on the first pixel below
 * the waterline, drawing a hard line across it.
 *
 * @note **Always asks where the surface is** rather than rejecting cheaply
 *       against the still plane first. A cheap rejection needs a bound on how
 *       tall a crest can be, and there is none to have: the amplitude is
 *       authored, the displacement comes off an FFT on the GPU, and nothing
 *       reads its extent back. A guessed bound does not clamp the waves, it
 *       just stops asking about them - so a sea that exceeds it goes back to
 *       being lit as air above the guess, with a hard line across everything
 *       standing in it and no clue as to why. One fetch is the price of not
 *       having that failure mode.
 *
 * @note Uses the world's surface, NOT the drawn one, and so disagrees with the
 *       fog beyond the distance the waves flatten out over. Deliberate: this
 *       shades a point, and a shaded point may not move with the camera. The
 *       disagreement is a metre or two of depth on something far enough away
 *       that the fog dominates how it looks; the alternative was a waterline
 *       that slides up a tree trunk as you walk towards it, and cached
 *       world-space lighting that drifts with the view.
 *
 * @param[in] position - World position of the shaded point.
 *
 * @return The medium. Inactive when there is no ocean or the point is at or
 *         above the surface, in which case nothing attenuates.
 */
WaterVolumetrics GetWaterVolumetrics(float3 position, float waterHeight)
{
	return MakeWaterVolumetrics(
		(half)saturate(
			(waterHeight - position.y) / WATER_VOLUMETRICS_FADE_DEPTH),
		waterHeight);
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
 * **The surface height it carries stays the still plane, deliberately**, unlike
 * the shaded-point medium above. This one is handed to a MARCH, whose samples
 * are spread along the view ray and can be a long way off; the height is what
 * `SubmergedLightPath` clips a light against at each of them. The crest the
 * camera happens to be under says nothing about the surface over a sample a
 * hundred metres away, so using it would replace an unbiased answer - the plane
 * is the mean surface - with a biased one. Measuring per sample is the honest
 * alternative and costs a displacement fetch per step per light.
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
 * light_color *= WaterLightTransmittance(
 *     surface.P, surface.waterSurfaceHeight, L, dist);
 * @endcode
 *
 * @param[in] samplePos - The lit point, in world space.
 * @param[in] waterHeight - World height of the water surface over that point.
 * @param[in] toLight - Normalized direction from the lit point to the light.
 * @param[in] distanceToLight - Distance to the light, or `FLT_MAX` for a
 *                              directional light.
 *
 * @return Per-channel surviving fraction, in [0, 1]. 1 at or above the surface,
 *         so callers need no test of their own.
 */
half3 WaterLightTransmittance(
	float3 samplePos,
	float waterHeight,
	float3 toLight,
	float distanceToLight
)
{
	const WaterVolumetrics water = GetWaterVolumetrics(samplePos, waterHeight);

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
 * lighting.indirect.diffuse *=
 *     WaterAmbientTransmittance(surface.P, surface.waterSurfaceHeight);
 * @endcode
 *
 * @param[in] samplePos - The lit point, in world space.
 * @param[in] waterHeight - World height of the water surface over that point.
 *
 * @return Per-channel surviving fraction, in [0, 1]. 1 at or above the surface,
 *         so callers need no test of their own.
 */
half3 WaterAmbientTransmittance(float3 samplePos, float waterHeight)
{
	const WaterVolumetrics water = GetWaterVolumetrics(samplePos, waterHeight);

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
