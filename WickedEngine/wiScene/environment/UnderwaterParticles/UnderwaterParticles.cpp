/**
 * Description of the suspended detritus drifting in the water.
 *
 * The constants below are the whole of the field's physical description. They
 * are file-local on purpose: nothing outside should be able to set them
 * directly, because the interesting ones are about to stop being constants and
 * start being derived from the water medium.
 */

#include "wiScene/environment/UnderwaterParticles/UnderwaterParticles.h"

#include <algorithm>
#include <cmath>
#include <cstdint>

using namespace DirectX;
using namespace wi::scene::environment;

namespace
{
	/**
	 * Shortest distance the field is allowed to reach, in metres.
	 *
	 * The wrapping face sits one reach from the eye, and the murkiest water
	 * here can only be seen through for about 60 cm. Following that literally
	 * would recycle particles inside arm's reach, where the eye tracks
	 * individual motes and would catch them vanishing. Everything past the
	 * true sighting range is hidden by the haze regardless, so the extra
	 * volume is cheap insurance.
	 */
	constexpr float MIN_FIELD_REACH = 4.0F;

	/**
	 * Furthest the field is allowed to reach, in metres.
	 *
	 * Economy rather than physics. Past roughly this distance a particle is
	 * enlarged to the renderer's minimum sprite size and dimmed by the square
	 * of that enlargement, which leaves it contributing a hundredth of its
	 * light or less. Volume grows as the cube of the reach, so following the
	 * clearest ocean water out to its full fifty metres would multiply the
	 * particle count by two orders of magnitude to add nothing anyone can see.
	 */
	constexpr float MAX_FIELD_REACH = 25.0F;

	/**
	 * Visible flakes per cubic metre, per unit of turbidity.
	 *
	 * A **calibrated proportionality, not an inversion**: see the note on
	 * UnderwaterParticles::ParticleCount for why a scattering coefficient
	 * cannot be turned into a number density of visible flakes. Chosen so that
	 * moderately turbid water - Jerlov oceanic III, turbidity 0.5 - puts a
	 * little over a hundred motes inside the few metres where they read
	 * individually, which is populated without being soup.
	 */
	constexpr float PARTICLES_PER_TURBIDITY = 1.1F;

	/**
	 * Largest number of particles that may be drawn.
	 *
	 * The presets stay well under this; it exists so that a hand-authored
	 * medium cannot turn one draw call into an unbounded one.
	 */
	constexpr uint32_t MAX_PARTICLE_COUNT = 32768;

	/**
	 * Radius of a single particle, in metres.
	 *
	 * Marine snow aggregates run from well under a millimetre to about a
	 * centimetre across, the larger flakes being the ones the eye actually
	 * picks out. This sits near the top of that range for exactly that reason.
	 */
	constexpr float PARTICLE_RADIUS = 0.006F;

	/**
	 * Velocity the current carries the whole field at, in metres per second.
	 *
	 * A few centimetres per second, which is an unremarkable coastal current
	 * and is fast enough to read as motion without the field streaking.
	 *
	 * The downward component is much smaller than the horizontal one and is
	 * what makes this *snow*: aggregates settle, but only at metres per hour,
	 * so any real current overwhelms the sinking many times over.
	 */
	constexpr XMFLOAT3 DRIFT_VELOCITY = XMFLOAT3(0.05F, -0.008F, 0.03F);
} // namespace

/*
################################################################################
Public
################################################################################
*/

// Constructors
//==============================================================================

UnderwaterParticles::UnderwaterParticles(const WaterMedium& medium) noexcept
	: props{ medium }
{
}

// Methods
//==============================================================================

float UnderwaterParticles::FieldSize() const noexcept
{
	const float reach = std::clamp(
		props.medium.VisibilityDistance(), MIN_FIELD_REACH, MAX_FIELD_REACH);

	// The viewer stands at the centre, so the cube's edge is twice its reach.
	return 2.0F * reach;
}

uint32_t UnderwaterParticles::ParticleCount() const noexcept
{
	const float size = FieldSize();
	const float density = PARTICLES_PER_TURBIDITY * props.medium.Turbidity();
	const float count = density * size * size * size;

	return static_cast<uint32_t>(std::clamp(
		count, 0.0F, static_cast<float>(MAX_PARTICLE_COUNT)));
}

float UnderwaterParticles::ParticleRadius() const noexcept
{
	return PARTICLE_RADIUS;
}

XMFLOAT3 UnderwaterParticles::DriftOffset(const float time) const noexcept
{
	const float period = FieldSize();

	return XMFLOAT3(
		std::fmod(DRIFT_VELOCITY.x * time, period),
		std::fmod(DRIFT_VELOCITY.y * time, period),
		std::fmod(DRIFT_VELOCITY.z * time, period)
	);
}
