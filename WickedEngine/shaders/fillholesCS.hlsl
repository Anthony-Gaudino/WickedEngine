#include "globals.hlsli"
#include "ShaderInterop_Postprocess.h"

/**
 * Fills the empty regions of a mipped texture from the content that surrounds
 * them.
 *
 * Runs as two sweeps over one mip chain. The **pull** sweep (this file's
 * default permutation) builds each coarser level out of whichever of its four
 * children hold anything, so content spreads upwards until some level covers
 * the hole entirely. The **push** sweep (`PUSH`) then walks back down, writing
 * into empty texels alone and leaving everything already filled untouched.
 *
 * A texel counts as empty when all three channels are exactly zero, which is
 * the state a masked `downsample4xCS` leaves wherever its whole footprint was
 * rejected. Scene content that is exactly black is therefore treated as a hole
 * and filled from its neighbours; the substitution is between one black region
 * and another, so it cannot be seen.
 *
 * The push weights each of its four taps by whether that tap holds anything,
 * rather than assuming the level above came out full. An unweighted bilinear
 * tap would blend the zero of a still-empty neighbour into the result and
 * darken it, which reads as a shadow around every hole instead of a fill.
 *
 * References:
 * https://en.wikipedia.org/wiki/Inpainting
 *
 * @note Both permutations expect `postprocess.resolution` to describe the level
 *       being written and `postprocess.params0.xy` the level being read.
 */

PUSHCONSTANT(postprocess, PostProcess);

/**
 * The level being read: the finer one when pulling, the coarser when pushing.
 */
Texture2D<float4> input : register(t0);

/** The level being written. */
RWTexture2D<float4> output : register(u0);

/**
 * True when a texel holds no content at all.
 *
 * @param[in] color - The texel's colour.
 *
 * @return true when the texel is empty and may be filled.
 */
bool IsEmpty(float3 color)
{
	return !any(color > 0);
}

[numthreads(POSTPROCESS_BLOCKSIZE, POSTPROCESS_BLOCKSIZE, 1)]
void main(uint3 DTid : SV_DispatchThreadID)
{
	if (any(DTid.xy >= postprocess.resolution))
		return;

	const int2 read_dim = int2(postprocess.params0.xy);

#ifdef PUSH
	// Already answered by the level below, so the pull's coarser guess has
	// nothing to add here.
	if (!IsEmpty(output[DTid.xy].rgb))
		return;

	const float2 uv = (DTid.xy + 0.5f) * postprocess.resolution_rcp;
	const float2 coord = (uv * read_dim) - 0.5f;
	const int2 corner = int2(floor(coord));
	const float2 subtexel = coord - corner;

	float3 sum = 0;
	float weight = 0;

	[unroll]
	for (int y = 0; y < 2; ++y)
	{
		[unroll]
		for (int x = 0; x < 2; ++x)
		{
			const float3 color =
				input[clamp(corner + int2(x, y), 0, read_dim - 1)].rgb;

			[branch]
			if (!IsEmpty(color))
			{
				const float tap =
					(x == 0 ? 1 - subtexel.x : subtexel.x) *
					(y == 0 ? 1 - subtexel.y : subtexel.y);
				sum += color * tap;
				weight += tap;
			}
		}
	}

	// Still nothing within reach, so leave the hole for a coarser level that
	// has already been filled to answer on the next sweep down.
	if (weight > 0)
	{
		output[DTid.xy] = float4(sum / weight, 1);
	}
#else
	const int2 base = int2(DTid.xy) * 2;

	float3 sum = 0;
	float weight = 0;

	[unroll]
	for (int y = 0; y < 2; ++y)
	{
		[unroll]
		for (int x = 0; x < 2; ++x)
		{
			const float3 color =
				input[min(base + int2(x, y), read_dim - 1)].rgb;

			[branch]
			if (!IsEmpty(color))
			{
				sum += color;
				weight += 1;
			}
		}
	}

	output[DTid.xy] = float4(weight > 0 ? sum / weight : 0, 1);
#endif // PUSH
}
