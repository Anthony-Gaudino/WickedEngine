#include "globals.hlsli"
#include "ShaderInterop_Postprocess.h"

PUSHCONSTANT(postprocess, PostProcess);

Texture2D<float4> input : register(t0);

#ifdef REJECT_MASKED
/**
 * Coverage of whatever the caller wants the copy to describe.
 *
 * Source texels the mask leaves black are dropped from the average, so the
 * result is built out of the marked ones alone. The ocean uses it to keep the
 * scene copy free of everything standing in front of the water: refracting
 * those would draw them a second time, displaced by a surface they are nearer
 * than.
 */
Texture2D<float> reject_mask : register(t1);

// Throwaway. 0 off, 1 paints magenta wherever the whole footprint was
// rejected, 2 shows the mask itself.
#define REJECT_MASKED_DEBUG 0
#endif // REJECT_MASKED

RWTexture2D<float4> output : register(u0);

[numthreads(POSTPROCESS_BLOCKSIZE, POSTPROCESS_BLOCKSIZE, 1)]
void main(uint3 DTid : SV_DispatchThreadID)
{
	const float2 uv = (DTid.xy + 0.5f) * postprocess.resolution_rcp;

	const bool hdrToSRGB = postprocess.params0.x > 0;

	float4 color = 0;

	uint2 dim;
	input.GetDimensions(dim.x, dim.y);
	float2 dim_rcp = rcp(dim);

#ifdef REJECT_MASKED
	// One tap per source texel rather than four bilinear ones, because a
	// bilinear tap has already blended a rejected texel into its neighbours by
	// the time it is read, and no weight applied afterwards can separate them
	// again.
	const int2 base = int2(DTid.xy) * 4;
	float weight = 0;

	[unroll]
	for (int y = 0; y < 4; ++y)
	{
		[unroll]
		for (int x = 0; x < 4; ++x)
		{
			const int2 coord = min(base + int2(x, y), int2(dim) - 1);

			[branch]
			if (reject_mask[coord] > 0.5f)
			{
				color += tonemap(input[coord]);
				weight += 1;
			}
		}
	}

	if (weight > 0)
	{
		color /= weight;
	}
	else
	{
		// Nothing marked anywhere in the footprint. Left empty rather than
		// averaged, because the plain average here is made of nothing but the
		// texels the mask rejected - it would hand back exactly the content
		// this pass exists to remove. `Postprocess_FillHoles` answers for it
		// out of the surrounding content instead.
		color = 0;
	}

	float3 debug_color = 0;
	bool debug_override = false;
#if REJECT_MASKED_DEBUG == 1
	debug_color = float3(1, 0, 1);
	debug_override = weight <= 0;
#elif REJECT_MASKED_DEBUG == 2
	debug_color = weight / 16.0f;
	debug_override = true;
#endif // REJECT_MASKED_DEBUG
#else
	color += tonemap(input.SampleLevel(sampler_linear_clamp, uv + float2(-1, -1) * dim_rcp, 0));
	color += tonemap(input.SampleLevel(sampler_linear_clamp, uv + float2(1, -1) * dim_rcp, 0));
	color += tonemap(input.SampleLevel(sampler_linear_clamp, uv + float2(-1, 1) * dim_rcp, 0));
	color += tonemap(input.SampleLevel(sampler_linear_clamp, uv + float2(1, 1) * dim_rcp, 0));

	color /= 4.0f;
#endif // REJECT_MASKED

	color.rgb = inverse_tonemap(color.rgb);

#if defined(REJECT_MASKED) && REJECT_MASKED_DEBUG != 0
	if (debug_override)
	{
		color.rgb = debug_color;
	}
#endif // REJECT_MASKED_DEBUG

	if (hdrToSRGB)
	{
		// reverse of the IMAGE_FLAG_OUTPUT_COLOR_SPACE_LINEAR in imagePS
		//	this is used when scene HDR result is used for GUI background
		color.rgb /= 9.0f;
		color.rgb = ApplySRGBCurve_Fast(color.rgb);
	}

	output[DTid.xy] = color;
}
