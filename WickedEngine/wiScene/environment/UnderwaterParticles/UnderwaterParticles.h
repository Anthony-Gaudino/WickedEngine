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

namespace wi::scene::environment { class UnderwaterParticles; }

/**
 * Suspended detritus drifting in the water around the viewer.
 *
 * Natural water is never empty. Below the surface it carries a continuous
 * rain of organic aggregates - **marine snow** - together with mineral grains
 * and bubbles, and it is those catching the light that read as *underwater*
 * far more strongly than the colour of the water does.
 *
 * This class holds only the **description of the field**: how far it reaches,
 * how many particles it holds, how large they are and how fast the current
 * carries them. It knows nothing about graphics; the renderer asks it these
 * questions once a frame and puts the answers in a push constant.
 *
 * The particles deliberately add **no attenuation of their own**. What the
 * suspended load does to light passing through the water is already accounted
 * for, in full and per wavelength, by WaterMedium - the scattering it reports
 * *is* the scattering off these very particles. Letting them dim the scene
 * again would count the same physics twice. They contribute structure, not
 * extinction: something for a light shaft to fall on, and a sense of scale and
 * motion that a smooth exponential haze cannot give.
 *
 * Example usage:
 * @code
 * const UnderwaterParticles particles;
 * const uint32_t count = particles.ParticleCount();
 * const float size = particles.FieldSize();
 * const XMFLOAT3 drift = particles.DriftOffset(sceneTime);
 * @endcode
 */
class wi::scene::environment::UnderwaterParticles final
{
	/*
	############################################################################
	Public
	############################################################################
	*/
	public:

	// Constructors
	//==========================================================================

	/**
	 * Constructs a particle field with the default marine snow description.
	 */
	UnderwaterParticles() noexcept = default;

	// Methods
	//==========================================================================

	/**
	 * Returns the edge length of the cube the particles occupy.
	 *
	 * The field is centred on the viewer and wraps at this period, so the same
	 * number both places the far boundary and decides when a particle recycles.
	 *
	 * @return Edge length of the field, in metres.
	 */
	[[nodiscard]] float FieldSize() const noexcept;

	/**
	 * Returns how many particles the field holds.
	 *
	 * Together with FieldSize() this fixes the number density. The particles
	 * are spread **uniformly** through the cube rather than concentrated near
	 * the viewer, which is what the physical distribution actually is; the
	 * price is that most of them land far away and cover well under a pixel,
	 * and it is the renderer's job to stop those from aliasing.
	 *
	 * @return Number of particles to draw.
	 */
	[[nodiscard]] uint32_t ParticleCount() const noexcept;

	/**
	 * Returns the radius of a single particle.
	 *
	 * @return Particle radius, in metres.
	 */
	[[nodiscard]] float ParticleRadius() const noexcept;

	/**
	 * Returns how far the current has carried the whole field by a given time.
	 *
	 * A single uniform translation, since a current large enough to see is
	 * coherent over far more than the few tens of metres this field spans. The
	 * per-particle motion that keeps the field from reading as a rigid body is
	 * a small wobble, applied in the shader where it costs nothing to vary per
	 * particle.
	 *
	 * Reduced modulo FieldSize(), which the toroidal wrapping makes exact
	 * rather than approximate: a translation by a whole period maps the field
	 * onto itself.
	 *
	 * @param[in] time - Seconds since the scene started.
	 *
	 * @return Displacement of the field, in metres, each component within one
	 *         field period of zero.
	 */
	[[nodiscard]] DirectX::XMFLOAT3 DriftOffset(float time) const noexcept;
};
