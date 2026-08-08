/**
 * Description of the suspended detritus drifting in the water.
 *
 * The constants below are the whole of the field's physical description. They
 * are file-local on purpose: nothing outside should be able to set them
 * directly, because the interesting ones are about to stop being constants and
 * start being derived from the water medium.
 */

#include "wiScene/environment/UnderwaterParticles/UnderwaterParticles.h"

#include <cmath>
#include <cstdint>

using namespace DirectX;
using namespace wi::scene::environment;

namespace
{
	/**
	 * Edge length of the cube the particles occupy, in metres.
	 *
	 * Large enough that the far face sits beyond the distance at which the
	 * water's own haze has swallowed everything in all but the clearest
	 * presets, so the boundary falls where there is nothing left to see.
	 */
	constexpr float FIELD_SIZE = 40.0F;

	/**
	 * Number of particles filling the field.
	 *
	 * Spread through 64000 cubic metres this is roughly one particle every two
	 * cubic metres, which puts a few hundred within the five metres where they
	 * are individually resolvable and leaves the rest to merge into texture.
	 * Six vertices each, so the draw is comparable in size to the ocean
	 * surface's own clipmap.
	 */
	constexpr uint32_t PARTICLE_COUNT = 32768;

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

// Methods
//==============================================================================

float UnderwaterParticles::FieldSize() const noexcept
{
	return FIELD_SIZE;
}

uint32_t UnderwaterParticles::ParticleCount() const noexcept
{
	return PARTICLE_COUNT;
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
