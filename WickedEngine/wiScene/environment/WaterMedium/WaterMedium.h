#pragma once

#include <cstdint>
#if __has_include(<DirectXMath.h>)
    // In this case, DirectXMath is coming from Windows SDK.
    // It is better to use this on Windows as some Windows libraries could
    // depend on the same DirectXMath headers
	#include <DirectXMath.h>
#else
    // In this case, DirectXMath is coming from supplied source code
    // On platforms that don't have Windows SDK, the source code for DirectXMath
    // is provided as part of the engine utilities
	#include "Utility/DirectXMath/DirectXMath.h"
#endif

namespace wi::scene::environment
{
	/**
	 * Jerlov optical water types.
	 *
	 * The Jerlov classification sorts natural waters by how strongly they
	 * attenuate downwelling daylight, from the clearest open ocean (**oceanic
	 * type I**) through progressively more turbid coastal waters (**types 1C
	 * to 9C**). Each entry here selects a preset mix of the `algae`, `silt` and
	 * `stain` constituents in WaterMedium, chosen so the resulting visibility
	 * distance follows the classical ladder (roughly 50 m for oceanic type I
	 * down to well under a metre for coastal type 9C).
	 *
	 * References:
	 * https://en.wikipedia.org/wiki/Jerlov_water_type
	 * Ocean Optics Web Book, "Jerlov Water Types":
	 * https://www.oceanopticsbook.info
	 */
	enum class JerlovWaterType : uint32_t
	{
		/** Hand-authored coefficients; no preset is applied. */
		CUSTOM,

		/** Oceanic type I: clearest open ocean, deep blue. */
		OCEANIC_I,

		/** Oceanic type IA: very clear open ocean. */
		OCEANIC_IA,

		/** Oceanic type IB: clear open ocean. */
		OCEANIC_IB,

		/** Oceanic type II: moderately clear ocean, blue-green. */
		OCEANIC_II,

		/** Oceanic type III: turbid ocean / clear shelf water, green. */
		OCEANIC_III,

		/** Coastal type 1C: clearest coastal water. */
		COASTAL_1C,

		/** Coastal type 3C: typical coastal water, green. */
		COASTAL_3C,

		/** Coastal type 5C: turbid coastal water, green-brown. */
		COASTAL_5C,

		/** Coastal type 7C: very turbid harbour water. */
		COASTAL_7C,

		/** Coastal type 9C: murkiest coastal water, brown. */
		COASTAL_9C,

		/** Number of water types. Not a selectable type. */
		COUNT
	};

	class WaterMedium;
} // namespace wi::scene::environment

/**
 * Inherent optical properties of a body of water.
 *
 * Describes the water as a participating medium through the two coefficients
 * that radiative transfer actually needs, each per RGB channel and in units of
 * inverse metres:
 *
 * - **Absorption** \f$\sigma_a\f$ removes light from the beam. It is strongly
 *   wavelength selective in water (red dies within a couple of metres, blue
 *   travels furthest), and it is what makes deep clear water look dark blue.
 * - **Scattering** \f$\sigma_s\f$ redirects light rather than destroying it.
 *   Pure water barely scatters at all, so **scattering is essentially the
 *   definition of turbidity**: suspended particulate matter is what makes
 *   water look *milky*, collapsing contrast long before it collapses
 *   brightness.
 *
 * Together they give the extinction (beam attenuation) coefficient
 * \f$\sigma_t = \sigma_a + \sigma_s\f$, and the single-scattering albedo
 * \f$\omega_0 = \sigma_s / \sigma_t\f$ that decides whether the medium reads
 * as a dark blue filter or a bright haze.
 *
 * Rather than exposing six raw coefficients, the medium is authored as pure
 * water plus the three constituents that actually dirty natural water, and the
 * spectra are summed from them:
 *
 * \f[
 * \sigma_a(\lambda) = a_{water}(\lambda)
 *                   + \text{algae} \cdot a_{ph}(\lambda)
 *                   + \text{silt} \cdot a_{min}(\lambda)
 *                   + \text{stain} \cdot a_{CDOM}(\lambda)
 * \f]
 * \f[
 * \sigma_s(\lambda) = b_{water}(\lambda)
 *                   + \text{algae} \cdot b_{ph}(\lambda)
 *                   + \text{silt} \cdot b_{min}(\lambda)
 * \f]
 *
 * The three are kept apart because they are **not interchangeable amounts of
 * the same dirt** - they differ in single-scattering albedo, which is what the
 * eye actually reads:
 *
 * - **Silt** scatters hard and absorbs little, so it raises \f$\omega_0\f$ and
 *   gives a bright, milky, near-neutral haze that destroys contrast.
 * - **Algae** scatters *and* absorbs selectively (blue and red, sparing green),
 *   so \f$\omega_0\f$ stays lower and the water goes green and darker.
 * - **Stain** only absorbs, adding no haze at all - it just drives the hue
 *   through yellow to brown and takes light away.
 *
 * References:
 * Ocean Optics Web Book, "Inherent Optical Properties":
 * https://www.oceanopticsbook.info
 *
 * Example usage:
 * @code
 * WaterMedium medium(JerlovWaterType::COASTAL_3C);
 * const XMFLOAT3 extinction = medium.Extinction(); // sigma_t, 1/m per channel
 * const float visibility = medium.VisibilityDistance(); // ~2 m
 * @endcode
 *
 * @note All coefficients are given at the RGB reference wavelengths
 *       620 / 540 / 460 nm; this is a three-band approximation of a continuous
 *       spectral quantity, not a spectral renderer.
 */
class wi::scene::environment::WaterMedium final
{
	/*
	############################################################################
	Private
	############################################################################
	*/
	private:

	// Properties
	//==========================================================================

	/**
	 * Stores the authored description of the medium.
	 *
	 * The optical coefficients are derived from these on demand rather than
	 * cached, so the struct stays trivially copyable and serializable.
	 */
	struct WaterMediumProps
	{
		/**
		 * Living phytoplankton, as chlorophyll-a concentration in mg/m^3.
		 *
		 * Absorbs in the blue (~440 nm) and again in the red (~675 nm) while
		 * sparing green, and scatters moderately. That double-ended absorption
		 * is what makes productive water read as green and comparatively dark
		 * rather than merely hazy. Roughly 0.05 in open ocean, 1-5 in coastal
		 * water, tens in a bloom. Must not be negative.
		 *
		 * The defaults across all three constituents are Jerlov oceanic type
		 * III.
		 */
		float algae = 2.5F;

		/**
		 * Suspended mineral matter, in g/m^3.
		 *
		 * Scatters strongly and almost achromatically while absorbing very
		 * little, so it raises the single-scattering albedo: this is the
		 * constituent that makes water bright and milky and destroys contrast
		 * without making it dark. Its weak absorption does climb towards blue
		 * (iron oxides), which tints heavy loads brown. Roughly 0 offshore,
		 * 1-10 in coastal or river water, far more in a stirred estuary. Must
		 * not be negative.
		 */
		float silt = 0.25F;

		/**
		 * Dissolved organic staining, as CDOM absorption at 440 nm, in 1/m.
		 *
		 * CDOM (coloured dissolved organic matter, historically *gelbstoff* or
		 * "yellow substance") is dissolved rather than suspended, so it absorbs
		 * without scattering at all - it darkens and tints but never hazes. It
		 * absorbs blue far more strongly than red, which is what drives water
		 * through green to brown as it accumulates. Must not be negative.
		 */
		float stain = 0.18F;

		/**
		 * Preset the coefficients were last taken from.
		 *
		 * Purely descriptive: it records which JerlovWaterType was applied so
		 * a UI can show it. It never feeds into the coefficients itself.
		 */
		JerlovWaterType waterType = JerlovWaterType::OCEANIC_III;
	} props;

	/*
	############################################################################
	Public
	############################################################################
	*/
	public:

	// Constructors
	//==========================================================================

	/**
	 * Constructs a medium describing moderately turbid ocean water.
	 *
	 * @note Defaults to Jerlov oceanic type III (roughly 6 m of visibility).
	 */
	WaterMedium() noexcept = default;

	/**
	 * Constructs a medium from a Jerlov water type preset.
	 *
	 * @param[in] waterType - Preset to apply. JerlovWaterType::CUSTOM leaves
	 *                        the default (oceanic type III) constituents in
	 *                        place.
	 */
	explicit WaterMedium(JerlovWaterType waterType) noexcept;

	// Getters
	//==========================================================================

	/**
	 * Returns the phytoplankton concentration.
	 *
	 * @return Chlorophyll-a (in mg/m^3). Always non-negative.
	 */
	[[nodiscard]] float GetAlgae() const noexcept;

	/**
	 * Returns the suspended mineral load.
	 *
	 * @return Mineral matter (in g/m^3). Always non-negative.
	 */
	[[nodiscard]] float GetSilt() const noexcept;

	/**
	 * Returns the dissolved organic staining.
	 *
	 * @return CDOM absorption at 440 nm (in 1/m). Always non-negative.
	 */
	[[nodiscard]] float GetStain() const noexcept;

	/**
	 * Returns the preset the coefficients were last taken from.
	 *
	 * @return The recorded water type. JerlovWaterType::CUSTOM when the
	 *         coefficients were authored by hand.
	 */
	[[nodiscard]] JerlovWaterType GetWaterType() const noexcept;

	// Setters
	//==========================================================================

	/**
	 * Sets the phytoplankton concentration.
	 *
	 * @param[in] algae - Chlorophyll-a (in mg/m^3). Negative values are
	 *                    clamped to 0.
	 *
	 * @note This deliberately does **not** reset the recorded water type to
	 *       JerlovWaterType::CUSTOM, so that deserializing a preset followed
	 *       by its stored coefficients round-trips exactly. A UI that lets the
	 *       user drag this value should set the type to `CUSTOM` itself.
	 */
	void SetAlgae(float algae) noexcept;

	/**
	 * Sets the suspended mineral load.
	 *
	 * @param[in] silt - Mineral matter (in g/m^3). Negative values are clamped
	 *                   to 0.
	 *
	 * @note Does not reset the recorded water type; see SetAlgae().
	 */
	void SetSilt(float silt) noexcept;

	/**
	 * Sets the dissolved organic staining.
	 *
	 * @param[in] stain - CDOM absorption at 440 nm (in 1/m). Negative values
	 *                    are clamped to 0.
	 *
	 * @note Does not reset the recorded water type; see SetAlgae().
	 */
	void SetStain(float stain) noexcept;

	/**
	 * Applies a Jerlov water type preset.
	 *
	 * Overwrites both coefficients with the preset values and records the
	 * type.
	 *
	 * @param[in] waterType - Preset to apply. JerlovWaterType::CUSTOM only
	 *                        records the type and leaves the coefficients
	 *                        untouched.
	 */
	void SetWaterType(JerlovWaterType waterType) noexcept;

	/**
	 * Sets the coefficients to approximate a legacy scalar water density.
	 *
	 * Before this medium existed the water was described by a single density
	 * in 1/m, applied equally to every channel. This scales all three
	 * constituents together until they reproduce that attenuation, holding
	 * them in the proportion the Jerlov presets use so the result stays a
	 * plausible natural water rather than an arbitrary set of numbers.
	 *
	 * It is matched on the green channel, which is the least attenuated one
	 * across the range of densities this is meant to convert, and therefore
	 * the channel that decided how far the old density let you see.
	 *
	 * @param[in] density - Legacy water density (in 1/m). Values at or below
	 *                      the absorption floor of pure water give the
	 *                      clearest medium this model can represent.
	 *
	 * @note Intended for deserializing content authored before the medium
	 *       existed; records the water type as JerlovWaterType::CUSTOM.
	 */
	void SetFromLegacyDensity(float density) noexcept;

	// Methods
	//==========================================================================

	/**
	 * Computes the absorption coefficient of the medium.
	 *
	 * Pure water absorption plus every constituent that absorbs - all three of
	 * them do:
	 * \f[
	 * \sigma_a(\lambda) = a_{water}(\lambda)
	 *                   + \text{algae} \cdot a_{ph}(\lambda)
	 *                   + \text{silt} \cdot a_{min}(\lambda)
	 *                   + \text{stain} \cdot a_{CDOM}(\lambda)
	 * \f]
	 *
	 * References:
	 * Pope & Fry 1997, *Absorption spectrum (380-700 nm) of pure water*;
	 * Smith & Baker 1981, *Optical properties of the clearest natural waters*;
	 * Bricaud et al. 1995, *Variability in the chlorophyll-specific absorption
	 * coefficients of natural phytoplankton*.
	 *
	 * @return \f$\sigma_a\f$ per RGB channel (in 1/m).
	 */
	[[nodiscard]] DirectX::XMFLOAT3 Absorption() const noexcept;

	/**
	 * Computes the scattering coefficient of the medium.
	 *
	 * Molecular (Rayleigh) scattering of pure water plus the two *suspended*
	 * constituents. Staining is dissolved, so it contributes nothing here -
	 * which is precisely why it darkens water without hazing it:
	 * \f[
	 * \sigma_s(\lambda) = b_{water}(\lambda)
	 *                   + \text{algae} \cdot b_{ph}(\lambda)
	 *                   + \text{silt} \cdot b_{min}(\lambda)
	 * \f]
	 *
	 * References:
	 * Morel 1974, *Optical properties of pure water and pure sea water*;
	 * Gordon & Morel 1983, *Remote assessment of ocean color*.
	 *
	 * @return \f$\sigma_s\f$ per RGB channel (in 1/m).
	 */
	[[nodiscard]] DirectX::XMFLOAT3 Scattering() const noexcept;

	/**
	 * Computes the overall turbidity of the medium.
	 *
	 * The total particulate scattering at 550 nm. Useful as a readout: it says
	 * how hazy the water is overall, while saying nothing about whether that
	 * haze is bright mineral or dark algal.
	 *
	 * @return Particulate scattering at 550 nm (in 1/m). Excludes the
	 *         molecular scattering of the water itself.
	 */
	[[nodiscard]] float Turbidity() const noexcept;

	/**
	 * Computes the extinction (beam attenuation) coefficient of the medium.
	 *
	 * \f[
	 * \sigma_t = \sigma_a + \sigma_s
	 * \f]
	 *
	 * @return \f$\sigma_t\f$ per RGB channel (in 1/m).
	 */
	[[nodiscard]] DirectX::XMFLOAT3 Extinction() const noexcept;

	/**
	 * Computes the Henyey-Greenstein asymmetry parameter of the medium.
	 *
	 * Pure water scatters almost symmetrically (\f$g \approx 0\f$) while
	 * marine particulates are extremely forward peaked (Petzold's average
	 * phase function gives \f$g \approx 0.924\f$). The effective asymmetry is
	 * the scattering-weighted blend of the two, so it rises from 0 towards
	 * 0.924 as turbidity takes over:
	 * \f[
	 * g = \frac{b_{particulate} \cdot g_{particulate}}
	 *          {b_{water} + b_{particulate}}
	 * \f]
	 *
	 * References:
	 * Petzold 1972, *Volume scattering functions for selected ocean waters*.
	 *
	 * @return Asymmetry parameter \f$g\f$ in [0, 0.924]. Evaluated on the
	 *         green channel, since the phase function is treated as
	 *         achromatic.
	 */
	[[nodiscard]] float PhaseAsymmetry() const noexcept;

	/**
	 * Computes the horizontal visibility distance through the medium.
	 *
	 * Uses Duntley's contrast transmittance law: the apparent contrast of a
	 * target falls as \f$e^{-\sigma_t r}\f$, and a target is taken to
	 * disappear once it drops to the conventional 2% threshold, giving
	 * \f[
	 * r = \frac{\ln 50}{\sigma_t} \approx \frac{3.912}{\sigma_t}
	 * \f]
	 *
	 * Evaluated on the *least* attenuated channel, because that is the
	 * wavelength an observer still sees the target by. Which channel wins
	 * migrates from blue in clear ocean water to green and finally red as the
	 * water turns murky, matching real observation.
	 *
	 * References:
	 * Duntley 1963, *Light in the sea*.
	 *
	 * @return Visibility distance (in metres). Bounded above for a perfectly
	 *         transparent medium.
	 *
	 * @note This is a horizontal sighting range, not a Secchi depth; the two
	 *       differ because Secchi depth also involves vertical daylight
	 *       attenuation.
	 */
	[[nodiscard]] float VisibilityDistance() const noexcept;
};
