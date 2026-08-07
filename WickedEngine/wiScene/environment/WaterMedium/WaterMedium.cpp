/**
 * Implementation of the water inherent optical properties model.
 *
 * All spectra below are sampled at the RGB reference wavelengths
 * 620 / 540 / 460 nm.
 */

#include "wiScene/environment/WaterMedium/WaterMedium.h"

#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef>

using namespace DirectX;
using namespace wi::scene::environment;

namespace
{
	/**
	 * Absorption coefficient of pure sea water, in 1/m.
	 *
	 * Red is absorbed roughly fifteen times faster than blue, which is why
	 * depth alone desaturates a scene towards blue even in perfectly clear
	 * water.
	 *
	 * Sources:
	 * Pope & Fry 1997, *Absorption spectrum (380-700 nm) of pure water*;
	 * Smith & Baker 1981, *Optical properties of the clearest natural
	 * waters*.
	 */
	constexpr XMFLOAT3 PURE_WATER_ABSORPTION = XMFLOAT3(0.45F, 0.08F, 0.03F);

	/**
	 * Molecular (Rayleigh) scattering coefficient of pure sea water, in 1/m.
	 *
	 * Follows a \f$\lambda^{-4.3}\f$ law, so blue scatters most. The magnitude
	 * is tiny compared to any realistic particulate load, which is precisely
	 * why turbidity dominates the scattering term in practice.
	 *
	 * Sources:
	 * Morel 1974, *Optical properties of pure water and pure sea water*.
	 */
	constexpr XMFLOAT3 PURE_WATER_SCATTERING =
		XMFLOAT3(0.0007F, 0.0015F, 0.0030F);

	/**
	 * Relative absorption spectrum of CDOM, normalized to 1 at 440 nm.
	 *
	 * Coloured dissolved organic matter absorbs following
	 * \f$a_g(\lambda) = a_g(440) e^{-S(\lambda - 440)}\f$ with a slope
	 * \f$S \approx 0.014\,\text{nm}^{-1}\f$. Blue is hit hardest, so
	 * accumulating CDOM shifts water from blue through green to brown.
	 *
	 * Sources:
	 * Bricaud, Morel & Prieur 1981, *Absorption by dissolved organic matter of
	 * the sea (yellow substance) in the UV and visible domains*.
	 */
	constexpr XMFLOAT3 CDOM_ABSORPTION_SPECTRUM =
		XMFLOAT3(0.08F, 0.25F, 0.76F);

	/**
	 * Absorption spectrum of phytoplankton, per mg/m^3 of chlorophyll-a.
	 *
	 * Chlorophyll absorbs in two bands - a strong one in the blue near 440 nm
	 * and a weaker one in the red near 675 nm - with a pronounced minimum in
	 * the green between them. That double-ended shape is what leaves green as
	 * the surviving channel and makes productive water read green rather than
	 * merely hazy.
	 *
	 * Sources:
	 * Bricaud et al. 1995, *Variability in the chlorophyll-specific absorption
	 * coefficients of natural phytoplankton*.
	 */
	constexpr XMFLOAT3 ALGAE_ABSORPTION_SPECTRUM =
		XMFLOAT3(0.012F, 0.008F, 0.045F);

	/**
	 * Scattering spectrum of phytoplankton, per mg/m^3 of chlorophyll-a.
	 *
	 * Cells are large compared to the wavelength, so this is nearly grey with
	 * only a weak tilt towards blue.
	 *
	 * Sources:
	 * Gordon & Morel 1983, *Remote assessment of ocean color for
	 * interpretation of satellite visible imagery*.
	 */
	constexpr XMFLOAT3 ALGAE_SCATTERING_SPECTRUM =
		XMFLOAT3(0.140F, 0.150F, 0.160F);

	/**
	 * Absorption spectrum of suspended minerals, per g/m^3.
	 *
	 * Weak overall, and decaying exponentially from blue towards red much like
	 * CDOM does - iron oxide content is what gives heavy mineral loads their
	 * brown cast. It is far smaller than the mineral scattering below, which is
	 * exactly why silt brightens water instead of darkening it.
	 *
	 * Sources:
	 * Babin et al. 2003, *Light scattering properties of marine particles in
	 * coastal and open ocean waters*.
	 */
	constexpr XMFLOAT3 MINERAL_ABSORPTION_SPECTRUM =
		XMFLOAT3(0.004F, 0.010F, 0.024F);

	/**
	 * Scattering spectrum of suspended minerals, per g/m^3.
	 *
	 * Mineral grains are large compared to the wavelength, so their scattering
	 * is nearly achromatic; only a weak \f$\lambda^{-0.5}\f$ tilt remains. This
	 * near-grey, high-albedo response is what gives silty water its milky,
	 * contrast-destroying haze instead of a colour cast.
	 *
	 * Sources:
	 * Babin et al. 2003, *Light scattering properties of marine particles in
	 * coastal and open ocean waters*.
	 */
	constexpr XMFLOAT3 MINERAL_SCATTERING_SPECTRUM =
		XMFLOAT3(0.470F, 0.500F, 0.540F);

	/**
	 * Share of pure water's scattering that is turned backwards.
	 *
	 * Molecular scattering is Rayleigh-like, whose phase function is symmetric
	 * about the scattering plane, so very nearly half of it heads back the way
	 * it came. Two orders of magnitude more than the particles below, which is
	 * why the water's own sharply selective colour survives at all once they
	 * are present.
	 *
	 * Sources:
	 * Morel 1974, *Optical properties of pure water and pure sea water*.
	 */
	constexpr float WATER_BACKSCATTER_FRACTION = 0.500F;

	/**
	 * Share of phytoplankton scattering that is turned backwards.
	 *
	 * Cells are large compared to the wavelength, so their scattering is
	 * overwhelmingly forward and only a few parts per thousand return. This is
	 * why productive water can be very hazy to see through and still show
	 * comparatively little colour of its own.
	 *
	 * Sources:
	 * Petzold 1972, *Volume scattering functions for selected ocean waters*.
	 */
	constexpr float ALGAE_BACKSCATTER_FRACTION = 0.005F;

	/**
	 * Share of mineral scattering that is turned backwards.
	 *
	 * Grains are larger still, but rougher and more angular than cells, so
	 * they return somewhat more - a few percent rather than a few tenths of
	 * one. Silt therefore does brighten water, just far less than its
	 * enormous \f$\sigma_s\f$ taken at face value would suggest.
	 *
	 * Sources:
	 * Babin et al. 2003, *Light scattering properties of marine particles*.
	 */
	constexpr float MINERAL_BACKSCATTER_FRACTION = 0.018F;

	/**
	 * Share of the photons phytoplankton absorb that are re-emitted at 685 nm.
	 *
	 * A quantum yield, so photons out per photon in regardless of how much
	 * energy each carries. Around one percent in open water, climbing several
	 * fold under nutrient stress or low light, when the cell has caught more
	 * than its photosynthesis can consume and sheds the excess.
	 *
	 * Sources:
	 * Babin et al. 1996, *Remote sensing of sea surface sun-induced chlorophyll
	 * fluorescence*.
	 */
	constexpr float CHLOROPHYLL_FLUORESCENCE_QUANTUM_YIELD = 0.01F;

	/**
	 * Share of the ENERGY absorbed in each band re-emitted as fluorescence.
	 *
	 * The quantum yield above counts photons, and the emitted ones are longer
	 * than the ones that were caught, so each carries proportionally less
	 * energy. Scaling by \f$\lambda_{ex} / 685\f$ converts one to the other.
	 *
	 * The line sits at 685 nm however it was excited - the cell relaxes to the
	 * same state before emitting - so the shortest wavelengths lose the most in
	 * the conversion.
	 */
	constexpr XMFLOAT3 CHLOROPHYLL_FLUORESCENCE_ENERGY_YIELD = XMFLOAT3(
		CHLOROPHYLL_FLUORESCENCE_QUANTUM_YIELD * (620.0F / 685.0F),
		CHLOROPHYLL_FLUORESCENCE_QUANTUM_YIELD * (540.0F / 685.0F),
		CHLOROPHYLL_FLUORESCENCE_QUANTUM_YIELD * (460.0F / 685.0F)
	);

	/**
	 * Henyey-Greenstein asymmetry parameter of marine particulates.
	 *
	 * Petzold's measured average ocean phase function is extremely forward
	 * peaked; this is the equivalent-asymmetry fit commonly used in its place.
	 *
	 * Sources:
	 * Petzold 1972, *Volume scattering functions for selected ocean waters*.
	 */
	constexpr float PARTICULATE_PHASE_ASYMMETRY = 0.924F;

	/**
	 * Optical depth at which apparent contrast falls to the 2% threshold.
	 *
	 * \f$\ln(1 / 0.02) = \ln 50 \approx 3.912\f$.
	 */
	constexpr float CONTRAST_THRESHOLD_OPTICAL_DEPTH = 3.912F;

	/**
	 * Smallest extinction used when deriving a visibility distance, in 1/m.
	 *
	 * Guards the division for a hypothetical perfectly transparent medium,
	 * which would otherwise report an infinite sighting range.
	 */
	constexpr float MIN_EXTINCTION = 1.0e-4F;

	/**
	 * Authored constituent mix for one Jerlov water type.
	 */
	struct JerlovPreset
	{
		/** Chlorophyll-a, in mg/m^3. */
		float algae;

		/** Suspended minerals, in g/m^3. */
		float silt;

		/** CDOM absorption at 440 nm, in 1/m. */
		float stain;
	};

	/**
	 * Converts a water type to its index in the preset table.
	 *
	 * @param[in] waterType - Water type to convert.
	 *
	 * @return Zero-based index into JERLOV_PRESETS.
	 */
	[[nodiscard]] constexpr size_t PresetIndex(
		const JerlovWaterType waterType
	) noexcept
	{
		return static_cast<size_t>(waterType);
	}

	/**
	 * Number of entries in the Jerlov preset table.
	 */
	constexpr size_t JERLOV_PRESET_COUNT =
		PresetIndex(JerlovWaterType::COUNT);

	/**
	 * Constituent presets for every JerlovWaterType, indexed by its value.
	 *
	 * The mixes were fitted so that WaterMedium::VisibilityDistance()
	 * reproduces the sighting ranges the Jerlov ladder is known for, quoted per
	 * row below. They also follow the real Case 1 / Case 2 split: the oceanic
	 * types are phytoplankton dominated with essentially no mineral load, while
	 * the coastal types take on progressively more silt until it carries most
	 * of the scattering. The `CUSTOM` row is never read; it only keeps the
	 * indices aligned with the enum.
	 */
	constexpr std::array<JerlovPreset, JERLOV_PRESET_COUNT> JERLOV_PRESETS = {{
		//  algae  silt   stain
		{ 0.20F, 0.00F, 0.01F }, // CUSTOM (unused, keeps indices aligned)
		{ 0.20F, 0.00F, 0.01F }, // OCEANIC_I
		{ 0.33F, 0.00F, 0.02F }, // OCEANIC_IA
		{ 0.53F, 0.00F, 0.03F }, // OCEANIC_IB
		{ 1.50F, 0.05F, 0.08F }, // OCEANIC_II
		{ 2.50F, 0.25F, 0.18F }, // OCEANIC_III
		{ 3.00F, 0.90F, 0.30F }, // COASTAL_1C
		{ 4.00F, 2.00F, 0.55F }, // COASTAL_3C
		{ 5.00F, 3.70F, 0.90F }, // COASTAL_5C
		{ 6.00F, 6.20F, 1.40F }, // COASTAL_7C
		{ 8.00F, 9.60F, 2.00F }, // COASTAL_9C
	}};

	/**
	 * Reference constituent mix used to resolve a single authored number.
	 *
	 * Oceanic type III, the mid point of the ladder. Scaling this whole mix
	 * gives a one-parameter family of plausible natural waters, which is what
	 * lets a legacy water density be converted without having to invent a
	 * split between the constituents.
	 */
	constexpr JerlovPreset REFERENCE_MIX =
		JERLOV_PRESETS[PresetIndex(JerlovWaterType::OCEANIC_III)];

	// The in-class member initializers of WaterMediumProps have to spell out
	// the oceanic type III values so the type stays trivially default
	// constructible; keep them from silently drifting away from the table.
	static_assert(
		REFERENCE_MIX.algae == 2.50F
		&& REFERENCE_MIX.silt == 0.25F
		&& REFERENCE_MIX.stain == 0.18F,
		"WaterMedium defaults must match the oceanic type III preset"
	);
} // namespace

/*
################################################################################
Public
################################################################################
*/

// Constructors
//==============================================================================

WaterMedium::WaterMedium(const JerlovWaterType waterType) noexcept
{
	SetWaterType(waterType);
}

// Getters
//==============================================================================

float WaterMedium::GetAlgae() const noexcept
{
	return props.algae;
}

float WaterMedium::GetSilt() const noexcept
{
	return props.silt;
}

float WaterMedium::GetStain() const noexcept
{
	return props.stain;
}

JerlovWaterType WaterMedium::GetWaterType() const noexcept
{
	return props.waterType;
}

// Setters
//==============================================================================

void WaterMedium::SetAlgae(const float algae) noexcept
{
	props.algae = std::max(0.0F, algae);
}

void WaterMedium::SetSilt(const float silt) noexcept
{
	props.silt = std::max(0.0F, silt);
}

void WaterMedium::SetStain(const float stain) noexcept
{
	props.stain = std::max(0.0F, stain);
}

void WaterMedium::SetWaterType(const JerlovWaterType waterType) noexcept
{
	assert(waterType < JerlovWaterType::COUNT && "Unknown Jerlov water type");

	props.waterType = waterType;

	if (waterType == JerlovWaterType::CUSTOM
		|| waterType >= JerlovWaterType::COUNT)
	{
		return;
	}

	const JerlovPreset& preset = JERLOV_PRESETS[PresetIndex(waterType)];
	props.algae = preset.algae;
	props.silt = preset.silt;
	props.stain = preset.stain;
}

void WaterMedium::SetFromLegacyDensity(const float density) noexcept
{
	// Extinction on the green channel with nothing dissolved or suspended at
	// all: the floor the legacy density cannot go below.
	constexpr float pureWaterFloor =
		PURE_WATER_ABSORPTION.y + PURE_WATER_SCATTERING.y;

	// Green extinction contributed by one unit of the reference mix, counting
	// every constituent in it.
	constexpr float extinctionPerReferenceUnit =
		(REFERENCE_MIX.algae
			* (ALGAE_ABSORPTION_SPECTRUM.y + ALGAE_SCATTERING_SPECTRUM.y))
		+ (REFERENCE_MIX.silt
			* (MINERAL_ABSORPTION_SPECTRUM.y + MINERAL_SCATTERING_SPECTRUM.y))
		+ (REFERENCE_MIX.stain * CDOM_ABSORPTION_SPECTRUM.y);

	const float scale = std::max(0.0F, density - pureWaterFloor)
		/ extinctionPerReferenceUnit;

	props.algae = REFERENCE_MIX.algae * scale;
	props.silt = REFERENCE_MIX.silt * scale;
	props.stain = REFERENCE_MIX.stain * scale;
	props.waterType = JerlovWaterType::CUSTOM;
}

// Methods
//==============================================================================

XMFLOAT3 WaterMedium::Absorption() const noexcept
{
	const float chlorophyll = GetAlgae();
	const float minerals = GetSilt();
	const float cdom = GetStain();

	return XMFLOAT3(
		PURE_WATER_ABSORPTION.x
			+ (chlorophyll * ALGAE_ABSORPTION_SPECTRUM.x)
			+ (minerals * MINERAL_ABSORPTION_SPECTRUM.x)
			+ (cdom * CDOM_ABSORPTION_SPECTRUM.x),
		PURE_WATER_ABSORPTION.y
			+ (chlorophyll * ALGAE_ABSORPTION_SPECTRUM.y)
			+ (minerals * MINERAL_ABSORPTION_SPECTRUM.y)
			+ (cdom * CDOM_ABSORPTION_SPECTRUM.y),
		PURE_WATER_ABSORPTION.z
			+ (chlorophyll * ALGAE_ABSORPTION_SPECTRUM.z)
			+ (minerals * MINERAL_ABSORPTION_SPECTRUM.z)
			+ (cdom * CDOM_ABSORPTION_SPECTRUM.z)
	);
}

XMFLOAT3 WaterMedium::Scattering() const noexcept
{
	const float chlorophyll = GetAlgae();
	const float minerals = GetSilt();

	return XMFLOAT3(
		PURE_WATER_SCATTERING.x
			+ (chlorophyll * ALGAE_SCATTERING_SPECTRUM.x)
			+ (minerals * MINERAL_SCATTERING_SPECTRUM.x),
		PURE_WATER_SCATTERING.y
			+ (chlorophyll * ALGAE_SCATTERING_SPECTRUM.y)
			+ (minerals * MINERAL_SCATTERING_SPECTRUM.y),
		PURE_WATER_SCATTERING.z
			+ (chlorophyll * ALGAE_SCATTERING_SPECTRUM.z)
			+ (minerals * MINERAL_SCATTERING_SPECTRUM.z)
	);
}

XMFLOAT3 WaterMedium::Backscattering() const noexcept
{
	const float chlorophyll = GetAlgae();
	const float minerals = GetSilt();

	return XMFLOAT3(
		(PURE_WATER_SCATTERING.x * WATER_BACKSCATTER_FRACTION)
			+ (chlorophyll * ALGAE_SCATTERING_SPECTRUM.x
				* ALGAE_BACKSCATTER_FRACTION)
			+ (minerals * MINERAL_SCATTERING_SPECTRUM.x
				* MINERAL_BACKSCATTER_FRACTION),
		(PURE_WATER_SCATTERING.y * WATER_BACKSCATTER_FRACTION)
			+ (chlorophyll * ALGAE_SCATTERING_SPECTRUM.y
				* ALGAE_BACKSCATTER_FRACTION)
			+ (minerals * MINERAL_SCATTERING_SPECTRUM.y
				* MINERAL_BACKSCATTER_FRACTION),
		(PURE_WATER_SCATTERING.z * WATER_BACKSCATTER_FRACTION)
			+ (chlorophyll * ALGAE_SCATTERING_SPECTRUM.z
				* ALGAE_BACKSCATTER_FRACTION)
			+ (minerals * MINERAL_SCATTERING_SPECTRUM.z
				* MINERAL_BACKSCATTER_FRACTION)
	);
}

XMFLOAT3 WaterMedium::Fluorescence() const noexcept
{
	const float chlorophyll = GetAlgae();

	// Only what the phytoplankton themselves caught can be re-emitted, so this
	// is their share of the absorption alone - the water, the minerals and the
	// staining all absorb without fluorescing.
	return XMFLOAT3(
		chlorophyll * ALGAE_ABSORPTION_SPECTRUM.x
			* CHLOROPHYLL_FLUORESCENCE_ENERGY_YIELD.x,
		chlorophyll * ALGAE_ABSORPTION_SPECTRUM.y
			* CHLOROPHYLL_FLUORESCENCE_ENERGY_YIELD.y,
		chlorophyll * ALGAE_ABSORPTION_SPECTRUM.z
			* CHLOROPHYLL_FLUORESCENCE_ENERGY_YIELD.z
	);
}

XMFLOAT3 WaterMedium::Extinction() const noexcept
{
	const XMFLOAT3 absorption = Absorption();
	const XMFLOAT3 scattering = Scattering();

	return XMFLOAT3(
		absorption.x + scattering.x,
		absorption.y + scattering.y,
		absorption.z + scattering.z
	);
}

float WaterMedium::Turbidity() const noexcept
{
	// Green channel, the 550 nm reference the constituent loads are quoted
	// against; the molecular scattering of the water itself is not turbidity.
	return (GetAlgae() * ALGAE_SCATTERING_SPECTRUM.y)
		+ (GetSilt() * MINERAL_SCATTERING_SPECTRUM.y);
}

float WaterMedium::PhaseAsymmetry() const noexcept
{
	// Green channel: the phase function is treated as achromatic. Algae and
	// minerals are not told apart here either - both are large particles with
	// strongly forward peaked scattering, and Petzold's average covers both
	// closely enough that separating them would be inventing precision.
	const float particles = Turbidity();
	const float total = PURE_WATER_SCATTERING.y + particles;

	// Pure water contributes symmetric (g = 0) scattering, so only the
	// particulate lobe appears in the numerator.
	return (particles * PARTICULATE_PHASE_ASYMMETRY) / total;
}

float WaterMedium::VisibilityDistance() const noexcept
{
	const XMFLOAT3 extinction = Extinction();

	// The least attenuated channel is the one the target is still seen by.
	const float bestChannel =
		std::min({ extinction.x, extinction.y, extinction.z });

	return CONTRAST_THRESHOLD_OPTICAL_DEPTH
		/ std::max(MIN_EXTINCTION, bestChannel);
}
