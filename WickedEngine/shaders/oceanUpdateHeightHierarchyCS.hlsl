/**
 * Builds one level of the ocean's height hierarchy.
 *
 * A min/max pyramid over the displacement map's height channel, so that a ray
 * can be tested against a whole region of sea at once instead of being sampled
 * along its length. Each texel holds the highest and lowest the drawn surface
 * reaches anywhere inside the patch region it covers, which is what lets a
 * trace reject a region outright - a ray passing entirely above the maximum, or
 * entirely below the minimum, cannot cross the surface there however far it
 * runs.
 *
 * That is the property a fixed-step march cannot have. Twelve samples spread
 * over a near-horizontal segment of several hundred metres land tens of metres
 * apart and step straight over a wave; the bound here is exact regardless of
 * how long the segment is.
 *
 * **The reduction is conservative, never averaged.** A mip chain built by
 * filtering answers "what is the height around here on average", which is the
 * wrong question and is not a bound at all. `GenerateMipChain` cannot be used
 * for this.
 *
 * References:
 * https://developer.nvidia.com/gpugems/gpugems3/part-iii-rendering/chapter-18-relaxed-cone-stepping-relief-mapping
 */

#include "globals.hlsli"
#include "ShaderInterop_Postprocess.h"
#include "ShaderInterop_Ocean.h"

PUSHCONSTANT(postprocess, PostProcess);

Texture2D<float4> displacementmap : register(t0);
Texture2D<float2> hierarchy_previous : register(t1);

RWTexture2D<float2> output : register(u0);

/**
 * Wraps a texel coordinate into a power-of-two map.
 *
 * The patch tiles, so the texel past the last edge is the first one again.
 * Sampling off the end without this would clamp and report a wall of water
 * that is not there.
 *
 * @param[in] coord - Texel coordinate, possibly outside the map.
 * @param[in] dim - Map size in texels. Must be a power of two.
 *
 * @return The coordinate wrapped into `[0, dim)`.
 */
uint2 WrapTexel(in int2 coord, in uint2 dim)
{
	return (uint2)coord & (dim - 1);
}

[numthreads(OCEAN_COMPUTE_TILESIZE, OCEAN_COMPUTE_TILESIZE, 1)]
void main(uint3 DTid : SV_DispatchThreadID)
{
	const uint2 outputDim = (uint2)postprocess.params0.xy;

	if (any(DTid.xy >= outputDim))
	{
		return;
	}

	float heightMax;
	float heightMin;

	[branch]
	if (postprocess.params0.z != 0)
	{
		// **The base level, reduced straight from the displacement map.**
		//
		// Four texels wide, from `2i - 1` to `2i + 2`, for a cell that only
		// spans `2i` and `2i + 1`. Both extra texels are needed and neither is
		// caution: the height field is read BILINEARLY, and a bilinear tap sits
		// half a texel back from the coordinate it is asked about. So the
		// surface just inside this cell's leading edge is interpolated from
		// `2i - 1`, and the surface just inside the next cell's is interpolated
		// from `2i + 2`. Bounding only the two texels the cell nominally covers
		// leaves a strip at each edge where the surface can stand outside the
		// bound - which is not a loose bound but a wrong one, and a trace built
		// on it steps through the surface exactly there.
		uint2 sourceDim;
		displacementmap.GetDimensions(sourceDim.x, sourceDim.y);

		const int2 base = int2(DTid.xy) * 2;

		heightMax = -FLT_MAX;
		heightMin = FLT_MAX;

		[unroll]
		for (int y = -1; y <= 2; ++y)
		{
			[unroll]
			for (int x = -1; x <= 2; ++x)
			{
				// The displacement map holds xzy, so height is `z`.
				const float height = displacementmap[
					WrapTexel(base + int2(x, y), sourceDim)].z;

				heightMax = max(heightMax, height);
				heightMin = min(heightMin, height);
			}
		}
	}
	else
	{
		// Each parent covers exactly the four children below it, and every
		// child already bounds its own span including the seams, so two by two
		// is enough here.
		uint2 sourceDim;
		hierarchy_previous.GetDimensions(sourceDim.x, sourceDim.y);

		const int2 base = int2(DTid.xy) * 2;

		heightMax = -FLT_MAX;
		heightMin = FLT_MAX;

		[unroll]
		for (int y = 0; y < 2; ++y)
		{
			[unroll]
			for (int x = 0; x < 2; ++x)
			{
				const float2 child = hierarchy_previous[
					WrapTexel(base + int2(x, y), sourceDim)];

				heightMax = max(heightMax, child.x);
				heightMin = min(heightMin, child.y);
			}
		}
	}

	output[DTid.xy] = float2(heightMax, heightMin);
}
