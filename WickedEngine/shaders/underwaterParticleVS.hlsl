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

/**
 * Smallest radius a particle may be drawn at, in pixels.
 *
 * A little over a pixel, so that a sprite always straddles enough of the grid
 * to be hit consistently as it moves. Larger is steadier still, but the
 * particle is spread over more of the image and reads as a soft blob rather
 * than a fleck.
 */
static const float MIN_PIXEL_RADIUS = 1.5;

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

	// Depth the perspective divide will use. A camera facing offset lies in the
	// view plane and so does not change it, which is why the whole sprite can
	// be sized from the centre alone.
	//
	// Floored at the near plane rather than at zero: behind the camera this
	// quantity goes negative, and a negative depth would ask for an enormous
	// sprite that, although centred behind the eye, would still reach across
	// the screen in front of it.
	const float viewDepth = max(
		dot(
			GetCamera().view_projection._41_42_43_44,
			float4(particleCenter, 1)
		),
		GetCamera().z_near
	);

	// Radius on screen, in pixels. The projection's vertical scale takes a
	// world length to normalized device coordinates, which span two units
	// across the viewport.
	const float pixelRadius = particles.particleRadius
		* GetCamera().projection._22 * 0.5
		* (float)GetCamera().internal_resolution.y / viewDepth;

	// Anything smaller than a pixel is missed by the rasterizer as often as it
	// is caught, and a mote that appears and disappears with sub-pixel camera
	// movement reads as a field of fireflies. So a distant particle is drawn
	// larger than it is and dimmed to match: spreading fixed light over
	// `scale^2` times the area costs it `scale^2` of its brightness, which is
	// the same total contribution the correctly sized sprite would have made.
	// The dimming is what keeps this from thickening the water - the field
	// must not fake attenuation the medium already accounts for.
	const float scale = max(1, MIN_PIXEL_RADIUS / max(pixelRadius, 0.00001));

	const float2 corner = QUAD_CORNERS[vertexID];

	// Camera facing: the sprite's plane is the view plane, so its corners are
	// the camera's right and up axes taken out of the inverse view matrix.
	const float3 worldPosition = particleCenter + mul(
		(float3x3)GetCamera().inverse_view,
		float3(corner * (particles.particleRadius * scale), 0)
	);

	UnderwaterParticleVertexToPixel output;
	output.opacity = (half)(1 / (scale * scale));
	output.position =
		mul(GetCamera().view_projection, float4(worldPosition, 1));
	output.P = worldPosition;
	output.corner = corner;
	output.clip = dot(float4(worldPosition, 1), GetCamera().clip_plane);

	return output;
}
