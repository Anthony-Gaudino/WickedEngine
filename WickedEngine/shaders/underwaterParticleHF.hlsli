#ifndef WI_UNDERWATERPARTICLE_HF
#define WI_UNDERWATERPARTICLE_HF
/**
 * Shared declarations for the underwater particle draw.
 *
 * Include after `globals.hlsli` (requires `GetCamera()`).
 */

#include "ShaderInterop_UnderwaterParticles.h"

PUSHCONSTANT(particles, UnderwaterParticlePushConstants);

/**
 * Interpolants passed from the particle vertex shader to the pixel shader.
 */
struct UnderwaterParticleVertexToPixel
{
	/** Clip space position. */
	float4 position : SV_POSITION;

	/**
	 * World space position of this corner of the sprite.
	 *
	 * Carried rather than reconstructed from the depth, because the sprite is
	 * a transparent and so has not written any depth to reconstruct from.
	 */
	float3 P : WORLDPOSITION;

	/** Position within the sprite, -1 to 1 on both axes from its centre. */
	float2 corner : PARTICLECORNER;

	/**
	 * Fraction of its light this particle still carries, 0 to 1.
	 *
	 * Below 1 only where the sprite had to be enlarged to stay resolvable, and
	 * exactly compensates that enlargement, so a particle emits the same total
	 * light however many pixels it was spread over.
	 */
	nointerpolation half opacity : PARTICLEOPACITY;

	/** Distance to the camera's clip plane, for planar reflection passes. */
	float clip : SV_ClipDistance0;
};

#endif // WI_UNDERWATERPARTICLE_HF
