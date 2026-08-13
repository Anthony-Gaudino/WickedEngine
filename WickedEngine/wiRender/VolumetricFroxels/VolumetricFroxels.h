#pragma once
#include "wiGraphicsDevice.h"
#include "wiScene_Decl.h"
#include "shaders/ShaderInterop_VolumetricFroxels.h"

namespace wi::render { class VolumetricFroxels; } // namespace wi::render

/**
 * Light scattered by the air and water in front of the camera, held in a
 * volume.
 *
 * A view-frustum aligned 3D grid, rebuilt every frame, in which each cell holds
 * the light scattered towards the eye by the medium **between the eye and that
 * cell**. Because it is indexed by distance, a fragment reads exactly the share
 * of the column that stands in front of *it*.
 *
 * **That is the whole point of the structure.** The alternative - integrating
 * the column once per pixel out to the opaque surface and blending the result
 * over the frame - produces a single number that cannot be split. Blended in
 * before the transparents it loses the light in front of them; blended in
 * after, it paints the light from *behind* them over the top, and a pane of
 * glass is hazed by air it does not stand in. Neither is recoverable from that
 * number, because a transparent cannot occlude light gathered against a depth
 * buffer it was never written into.
 *
 * References:
 * https://advances.realtimerendering.com/s2014/index.html
 * (Wronski, *Volumetric Fog: Unified Compute Shader Based Solution to
 * Atmospheric Scattering*, SIGGRAPH 2014)
 *
 * Example usage:
 * @code
 * volumetricFroxels.Create(500.0F);
 * volumetricFroxels.AdvanceFrame();
 * volumetricFroxels.Build(camera, cameraPrevious, cameraReflection, cmd);
 * // Publish GetVolume() and GetTail() to shaders, and sample per fragment.
 * @endcode
 *
 * @note The volume belongs to one camera. A reflection, probe or impostor view
 *       must not sample a volume built for the main eye, which is why the
 *       descriptor index lives on the camera rather than the frame.
 */
class wi::render::VolumetricFroxels final
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
	 * Light scattered per metre at each cell, before it is integrated.
	 *
	 * Two of them, because this is the half that is filtered over time: it is a
	 * **local** property of a point in the world, so the same point seen from a
	 * moved camera holds the same value and last frame's estimate can be
	 * reused. The integrated volume below cannot be filtered that way - it is a
	 * path integral along one particular ray, and a ray from somewhere else is
	 * not a correction to it.
	 */
	wi::graphics::Texture injectVolume[2];

	/**
	 * Light gathered from the eye out to each cell.
	 *
	 * What the fragment shaders sample. Cumulative rather than per-slice so a
	 * lookup is one texture fetch at the fragment's own distance. Everything
	 * beyond the volume's range reads its last slice, and takes the rest of the
	 * column from the tail below.
	 */
	wi::graphics::Texture integratedVolume;

	/**
	 * Light gathered from the end of the volume out to the far plane.
	 *
	 * One value per screen column rather than per cell, because the stretch it
	 * describes has no slices: the volume's last texel spans a few metres and
	 * this spans thousands. Its alpha carries the extinction over that stretch,
	 * which is the shape a fragment out there fades it in along - written into
	 * the last texel instead, the whole far field would arrive across those few
	 * metres as a hard arc at exactly the range.
	 */
	wi::graphics::Texture tailTexture;

	/**
	 * Frames built since the last reset, or -1 before the first.
	 *
	 * Drives the ping-pong and, more importantly, says whether there is any
	 * history to reuse. A camera cut leaves the previous frame describing a
	 * different world, so this resets rather than blending across it.
	 */
	mutable int frame = -1;

	/**
	 * How far the volume reaches from the eye, in metres.
	 *
	 * Authored rather than derived: it trades resolution against reach, and the
	 * right answer for a room and for an ocean horizon are orders of magnitude
	 * apart. The view is not bounded by it - see
	 * `VOLUMETRIC_FROXEL_TAIL_STEPS`.
	 */
	float range = 0.0F;

	/*
	############################################################################
	Public
	############################################################################
	*/
	public:

	// Constructors
	//==========================================================================

	/**
	 * Constructs an empty volume that allocates nothing until `Create`.
	 */
	VolumetricFroxels() noexcept = default;

	// Getters
	//==========================================================================

	/**
	 * Returns whether the volume has been allocated and can be built.
	 *
	 * @return true once `Create` has succeeded, false before it or after
	 *         `Reset`.
	 */
	[[nodiscard]] bool IsValid() const noexcept;

	/**
	 * Returns the volume the fragment shaders sample.
	 *
	 * @return The integrated volume, or an invalid texture when there is none.
	 */
	[[nodiscard]] const wi::graphics::Texture* GetVolume() const noexcept;

	/**
	 * Returns the column beyond the volume.
	 *
	 * @return The tail, or an invalid texture when there is none.
	 */
	[[nodiscard]] const wi::graphics::Texture* GetTail() const noexcept;

	/**
	 * Returns how far the volume reaches from the eye.
	 *
	 * @return Range in metres. Shaders need it to turn a distance into a
	 *         texture coordinate, so it has to be published alongside the
	 *         volume itself.
	 */
	[[nodiscard]] float GetRange() const noexcept;

	// Methods
	//==========================================================================

	/**
	 * Loads the shaders this uses, and re-loads them on request.
	 *
	 * Call once at startup, before any volume is built.
	 */
	static void Initialize();

	/**
	 * Allocates the volume.
	 *
	 * Safe to call when already allocated at the same range, in which case it
	 * does nothing - which is what lets a caller treat it as "make sure this
	 * exists" once per frame rather than tracking the resize itself.
	 *
	 * @param[in] range - How far the volume should reach, in metres.
	 */
	void Create(float range);

	/**
	 * Releases the volume.
	 */
	void Reset() noexcept;

	/**
	 * Discards the temporal history without releasing anything.
	 *
	 * For a camera cut, where the previous frame describes a different world
	 * and blending towards it would drag the old scene across several frames.
	 */
	void ResetFrame() const noexcept;

	/**
	 * Moves to the next frame's half of the ping-pong.
	 *
	 * Call exactly once per frame, before `Build`.
	 */
	void AdvanceFrame() const noexcept;

	/**
	 * Builds the volume for a camera.
	 *
	 * Two dispatches: one filling every cell with what the medium there
	 * scatters towards the eye, and one accumulating those along each view ray.
	 *
	 * Rebinds the camera constant buffer with the temporal anti-aliasing jitter
	 * removed, and restores it before returning. That jitter moves the whole
	 * frustum by a fraction of a pixel each frame, and a volume rebuilt on a
	 * frustum that will not sit still cannot be told from one whose contents
	 * are moving - the sample point inside each cell has to be the only thing
	 * that changes between frames, or the accumulation is averaging two things
	 * at once.
	 *
	 * @param[in] camera - The camera the volume is built for. Its frustum
	 *                     defines the grid, so a volume built for one camera
	 *                     means nothing to another.
	 * @param[in] cameraPrevious - Where that camera stood last frame, for
	 *                             reprojecting into the history.
	 * @param[in] cameraReflection - Reflection camera, carried through the
	 *                               rebind untouched.
	 * @param[in] cmd - Command list to record into.
	 *
	 * @note Records nothing when the volume is invalid, so a caller need not
	 *       test first.
	 */
	void Build(
		const wi::scene::CameraComponent& camera,
		const wi::scene::CameraComponent& cameraPrevious,
		const wi::scene::CameraComponent& cameraReflection,
		wi::graphics::CommandList cmd
	) const;

	/*
	############################################################################
	Private
	############################################################################
	*/
	private:

	// Getters
	//==========================================================================

	/**
	 * Returns which half of the ping-pong this frame writes.
	 *
	 * @return Index into `injectVolume`.
	 */
	[[nodiscard]] int GetInjectOutputIndex() const noexcept;

	/**
	 * Returns which half of the ping-pong holds last frame's estimate.
	 *
	 * @return Index into `injectVolume`. Meaningless until a frame has been
	 *         built - ask `frame` before trusting it.
	 */
	[[nodiscard]] int GetInjectHistoryIndex() const noexcept;
};
