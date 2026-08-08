#ifndef WI_SHADERINTEROP_UNDERWATERPARTICLES_H
#define WI_SHADERINTEROP_UNDERWATERPARTICLES_H
#include "ShaderInterop.h"

/**
 * Parameters of the suspended detritus filling the water around the viewer.
 *
 * The field has **no storage and no simulation**: a particle's position is a
 * closed-form function of its index and the time, so the whole draw is
 * described by the handful of scalars below. Everything an artist or the
 * medium controls has already been resolved into them on the CPU.
 *
 * Positions come out of a hash of the index, so the particles occupy a cube of
 * side `fieldSize` centred on `fieldCenter`, wrapped toroidally: a particle
 * leaving one face reappears at the opposite one. That is what keeps the field
 * around the viewer without anchoring it to the viewer - wrapping moves one
 * particle at a time, where snapping the whole cube to a lattice would jerk
 * every particle at once each time the camera crossed a quantum.
 *
 * Example usage:
 * @code
 * UnderwaterParticlePushConstants push;
 * push.fieldCenter = camera.Eye;
 * push.fieldSize = particles.FieldSize();
 * device->PushConstants(&push, sizeof(push), cmd);
 * device->DrawInstanced(6, push.particleCount, 0, 0, cmd);
 * @endcode
 */
struct UnderwaterParticlePushConstants
{
	/** World position the field is centred on, i.e. the viewer. */
	float3 fieldCenter;

	/**
	 * Edge length of the cube the particles occupy, in metres.
	 *
	 * Doubles as the wrapping period, so it is the one number that decides
	 * both how far the field reaches and where a particle recycles.
	 */
	float fieldSize;

	/**
	 * Displacement the current has carried the whole field by, in metres.
	 *
	 * Reduced modulo `fieldSize` on the CPU. It could be accumulated in the
	 * shader from the time instead, but that grows without bound and loses the
	 * low bits of a metre-scale offset within an hour of play.
	 */
	float3 driftOffset;

	/** Radius of one particle, in metres. */
	float particleRadius;

	/** Number of particles to draw. */
	uint particleCount;
};

#endif // WI_SHADERINTEROP_UNDERWATERPARTICLES_H
