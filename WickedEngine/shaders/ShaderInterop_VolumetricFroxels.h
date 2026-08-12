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
 * Steps spent carrying the last slice out to the far plane.
 *
 * The volume ends at an authored range, but the view does not. Everything
 * beyond that range samples the last slice, so that slice has to hold the
 * integral over the *whole* remaining column or a mountain at three kilometres
 * would collect only the light in front of the range - and the sun shafts that
 * cross a landscape would simply stop.
 */
static const uint VOLUMETRIC_FROXEL_TAIL_STEPS = 8;

/**
 * Constants both build passes need.
 *
 * Carried as push constants rather than folded into the generic postprocess
 * struct, so the two fields are named for what they are at every use.
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
};

/**
 * The slice mapping itself lives in `volumetricFroxelHF.hlsli`.
 *
 * It is written in terms of HLSL intrinsics, and this header also compiles as
 * C++. Nothing on the C++ side needs to convert a distance into a slice - the
 * shape published here is all it has to agree on.
 */

#endif // WI_SHADERINTEROP_VOLUMETRICFROXELS_H
