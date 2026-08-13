#ifndef WI_SHADERINTEROP_VOLUMETRICFROXELS_H
#define WI_SHADERINTEROP_VOLUMETRICFROXELS_H
/**
 * Shared description of the volumetric froxel volume.
 *
 * Compiled as both C++ and HLSL, so the passes that BUILD the volume and the
 * shaders that SAMPLE it read one description of its shape. Two places deriving
 * the same slice mapping is how a lookup drifts off the thing it is looking up,
 * and the error would be a smooth depth-dependent offset - the kind that looks
 * like a tuning problem rather than a bug.
 */
#include "ShaderInterop.h"

/**
 * Cells across the volume, and the thread group shape that fills it.
 *
 * Sized against what is at risk rather than against a memory budget: caustic
 * modulation inside an underwater shaft is the finest detail the volume has to
 * carry, and it is resolved at the volume's own resolution rather than the
 * screen's. Independent of output resolution on purpose - the volume is sampled
 * bilinearly, so a froxel spanning six pixels at 1080p and twelve at 4K is not
 * perceptible in a field this smooth, and a fixed size keeps the memory
 * predictable.
 */
static const uint VOLUMETRIC_FROXEL_WIDTH = 320;
static const uint VOLUMETRIC_FROXEL_HEIGHT = 180;
static const uint VOLUMETRIC_FROXEL_SLICES = 128;

static const uint VOLUMETRIC_FROXEL_BLOCKSIZE_X = 8;
static const uint VOLUMETRIC_FROXEL_BLOCKSIZE_Y = 8;
static const uint VOLUMETRIC_FROXEL_BLOCKSIZE_Z = 4;

/**
 * How far the volume reaches from the eye by default, in metres.
 *
 * Trades resolution against reach, and the right answer for a room and for an
 * ocean horizon are orders of magnitude apart - so this is a starting point, not
 * a law. The view is **not** bounded by it; see `VOLUMETRIC_FROXEL_TAIL_STEPS`.
 */
static const float VOLUMETRIC_FROXEL_DEFAULT_RANGE = 500.0F;

/**
 * Steps spent marching from the end of the volume out to the far plane.
 *
 * The volume ends at an authored range, but the view does not: without this a
 * mountain at three kilometres would collect only the light in front of the
 * range, and the sun shafts crossing a landscape would simply stop.
 *
 * Few, because the result is published as a single value per screen column
 * rather than as slices. There is no depth detail to resolve out there - only a
 * total and the rate it fades in at.
 */
static const uint VOLUMETRIC_FROXEL_TAIL_STEPS = 8;

/**
 * Largest radiance a cell may hold.
 *
 * `R11G11B10_FLOAT` carries a five bit exponent and stops just above 65024;
 * anything larger is stored as infinity. That is not a rounding error that
 * fades out - the integration pass carries a running sum along each column, so
 * one infinite cell makes every cell behind it infinite too, and a tone mapper
 * dividing by it turns the lot into black.
 *
 * Reachable in ordinary content rather than in abuse: a point light's falloff
 * is bounded only by `1 / max(0.0001, distanceSquared)`, so any cell landing
 * within a few centimetres of a bright lamp arrives here in the tens of
 * thousands.
 */
static const float VOLUMETRIC_FROXEL_MAX_RADIANCE = 60000.0F;

/**
 * Share of a cell taken from this frame's sample rather than from the history.
 *
 * One sample per cell per frame is far too few to describe what a bright local
 * light does across a cell, where the falloff goes as \f$1/d^2\f$; the sample
 * point is moved around inside the cell every frame so that the average over
 * several frames is what a great many samples in one frame would have given.
 * This is the length of that average.
 *
 * Small, because the quantity being averaged is a **local** property of a point
 * in the world: the same point seen from a moved camera holds the same value,
 * so an old estimate of it is not stale, only noisy. What does go stale is an
 * occluder moving through the light, and that is the cost this trades against.
 */
static const float VOLUMETRIC_FROXEL_HISTORY_BLEND = 0.05F;

/**
 * How far a reused estimate may stray from what its neighbours say now.
 *
 * A long average is only sound while the thing being averaged holds still. Move
 * an occluder through a shaft and the history describes light that is no longer
 * there, which arrives as a smear trailing the occluder for as many frames as
 * the average is long.
 *
 * The guard is that the cells around a given one are a fair sample of what it
 * should be: the medium is smooth, so the history is kept only where it is
 * consistent with the spread its neighbours show this frame, and pulled to the
 * edge of that spread where it is not. Measured in standard deviations rather
 * than clamped to the neighbours' outright range, which throws away most of the
 * averaging exactly where the noise is worst.
 */
static const float VOLUMETRIC_FROXEL_HISTORY_CLIP_SIGMA = 2.0F;

/**
 * Constants both build passes need.
 *
 * Carried as push constants rather than folded into the generic postprocess
 * struct, so every field is named for what it is at each use.
 */
struct VolumetricFroxelPush
{
	/** How far the volume reaches from the eye, in metres. */
	float range;

	/**
	 * Frames built since the last reset.
	 *
	 * Advances the jitter sequence, and says whether there is a history worth
	 * blending towards: zero means there is not.
	 */
	uint frame;

	float padding0;
	float padding1;

	/**
	 * Where the eye stood last frame, in world space.
	 *
	 * The volume is indexed by distance from the eye, so finding a point in
	 * last frame's volume needs last frame's eye - and a screen position alone
	 * cannot supply it. Passed in rather than recovered from the previous view
	 * matrix in the shader, where it would be three lines of transpose that
	 * quietly depend on which way round the engine stores its matrices.
	 */
	float3 previousEye;

	float padding2;
};

/**
 * The slice mapping itself lives in `volumetricFroxelHF.hlsli`.
 *
 * It is written in terms of HLSL intrinsics, and this header also compiles as
 * C++. Nothing on the C++ side needs to convert a distance into a slice - the
 * shape published here is all it has to agree on.
 */

#endif // WI_SHADERINTEROP_VOLUMETRICFROXELS_H
