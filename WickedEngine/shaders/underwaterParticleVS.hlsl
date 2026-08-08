/**
 * Underwater particle vertex shader.
 *
 * Builds a camera-facing sprite per particle straight out of `SV_InstanceID`
 * and `SV_VertexID`. There is no vertex buffer and no simulation: a particle's
 * position is a closed-form function of its index and the time, so the field
 * costs one draw call and nothing else.
 *
 * The field is a cube centred on the viewer and **wrapped toroidally** rather
 * than snapped to a lattice. Wrapping recycles one particle at a time, at the
 * far face where the water's haze has already hidden it; snapping the cube
 * would instead displace every particle at once each time the camera crossed a
 * quantum, which reads as the whole field twitching.
 */

#include "globals.hlsli"
#include "underwaterParticleHF.hlsli"

/**
 * Sprite corners, as two triangles.
 *
 * Wound for no particular facing: the field is rasterized with culling off,
 * since a sprite turned to face the camera has no meaningful back.
 */
static const float2 QUAD_CORNERS[6] =
{
	float2(-1, -1), float2( 1, -1), float2(-1,  1),
	float2(-1,  1), float2( 1, -1), float2( 1,  1),
};

/**
 * How far a particle strays from its drifting position, in metres.
 *
 * Small enough to be a jostle rather than a path. Without it every particle
 * shares one velocity exactly, and a field translating as a rigid body is read
 * by the eye as the camera moving rather than the water.
 */
static const float WOBBLE_AMPLITUDE = 0.06;

/**
 * Angular rate of the wobble, in radians per second.
 *
 * Slow: suspended matter is neutrally buoyant and moves with the water around
 * it, so anything brisk enough to notice as vibration is wrong.
 */
static const float WOBBLE_RATE = 0.8;

UnderwaterParticleVertexToPixel main(
	uint vertexID : SV_VertexID,
	uint instanceID : SV_InstanceID
)
{
	// Three decorrelated values per particle, placing it in the unit cube and
	// doubling as the phase of its wobble. Reusing them for the phase is safe
	// precisely because they are already the random part of the position:
	// neighbouring particles are no more alike in phase than they are in place.
	const float3 unit = float3(
		hash1(instanceID * 3 + 0),
		hash1(instanceID * 3 + 1),
		hash1(instanceID * 3 + 2)
	);

	float3 particleCenter = (unit - 0.5) * particles.fieldSize;

	particleCenter += WOBBLE_AMPLITUDE
		* sin((GetFrame().time * WOBBLE_RATE) + (unit * (2 * PI)));

	particleCenter += particles.driftOffset;

	// Into the cube around the viewer. `round` puts the particle in whichever
	// period is nearest, so the result is always within half a field of the
	// centre whatever the drift has accumulated.
	const float3 relative = particleCenter - particles.fieldCenter;
	particleCenter = particles.fieldCenter + relative
		- (particles.fieldSize * round(relative / particles.fieldSize));

	const float2 corner = QUAD_CORNERS[vertexID];

	// Camera facing: the sprite's plane is the view plane, so its corners are
	// the camera's right and up axes taken out of the inverse view matrix.
	const float3 worldPosition = particleCenter + mul(
		(float3x3)GetCamera().inverse_view,
		float3(corner * particles.particleRadius, 0)
	);

	UnderwaterParticleVertexToPixel output;
	output.position =
		mul(GetCamera().view_projection, float4(worldPosition, 1));
	output.P = worldPosition;
	output.corner = corner;
	output.clip = dot(float4(worldPosition, 1), GetCamera().clip_plane);

	return output;
}
