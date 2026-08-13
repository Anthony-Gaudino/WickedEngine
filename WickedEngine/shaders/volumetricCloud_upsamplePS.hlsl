#include "globals.hlsli"
#include "ShaderInterop_Postprocess.h"
#include "waterFogHF.hlsli"
#include "volumetricFroxelHF.hlsli"

PUSHCONSTANT(postprocess, PostProcess);

Texture2D<float4> cloud_current : register(t0);
Texture2D<float2> cloud_depth_current : register(t1);

static const int UPSAMPLE_SAMPLE_RADIUS = 1;

#define GAUSSIAN_SIGMA_SPATIAL 0.5
#define GAUSSIAN_SIGMA_RANGE 1000.0

half Gaussian(half x, half sigma)
{
	return exp(-x * x / (2.0 * sigma * sigma));
}

float4 main(float4 pos : SV_Position, float2 uv : TEXCOORD) : SV_Target
{
	const uint2 pixel = pos.xy;
	const float depth = texture_depth[pixel];
	
	const uint2 reprojectionCoord = pixel / 2;
	
	const float2 reprojectionResolution = postprocess.params0.xy;
	const float2 reprojectionTexelSize = postprocess.params0.zw;

	float3 depthWorldPosition = reconstruct_position(uv, depth);
	float tToDepthBuffer = length(depthWorldPosition - GetCamera().position);
	
	// If sky, set distance to infinite
	tToDepthBuffer = depth == 0.0 ? FLT_MAX : tToDepthBuffer;

	// Adapted from upsample_bilateral_float4CS:
	
	const float2 uv00 = uv - reprojectionTexelSize * 0.5;
	const float2 uv10 = float2(uv00.x + reprojectionTexelSize.x, uv00.y);
	const float2 uv01 = float2(uv00.x, uv00.y + reprojectionTexelSize.y);
	const float2 uv11 = float2(uv00.x + reprojectionTexelSize.x, uv00.y + reprojectionTexelSize.y);

	const float4 lineardepth_lowres = float4(
		cloud_depth_current.SampleLevel(sampler_point_clamp, uv00, 0).g,
		cloud_depth_current.SampleLevel(sampler_point_clamp, uv10, 0).g,
		cloud_depth_current.SampleLevel(sampler_point_clamp, uv01, 0).g,
		cloud_depth_current.SampleLevel(sampler_point_clamp, uv11, 0).g
	);

	const float4 depthDiff = abs(tToDepthBuffer - lineardepth_lowres);	
	float depthDiffMax = max(max(depthDiff.x, depthDiff.y), max(depthDiff.z, depthDiff.w));
	
	half4 result = 0;
	
	[branch]
	if (depthDiffMax < tToDepthBuffer * 0.2)
	{
		// small error, take bilinear sample:
		result = cloud_current.SampleLevel(sampler_linear_clamp, uv, 0);
	}
	else
	{
		// large error, calculate weight and color depending on depth difference with gaussian configuration
		half4 color = 0;
		float weightSum = 0; // Note: weights need full precision on Nvidia Vulkan!
				
		[unroll]
		for (int y = -UPSAMPLE_SAMPLE_RADIUS; y <= UPSAMPLE_SAMPLE_RADIUS; y++)
		{
			[unroll]
			for (int x = -UPSAMPLE_SAMPLE_RADIUS; x <= UPSAMPLE_SAMPLE_RADIUS; x++)
			{
				int2 offset = int2(x, y);
			
				int2 neighborReprojectionCoord = reprojectionCoord + offset;
				float2 neighborReprojectionUV = (neighborReprojectionCoord + 0.5) / reprojectionResolution;
			
				half4 cloudResult = cloud_current.SampleLevel(sampler_linear_clamp, neighborReprojectionUV, 0);
				half cloudDepth = cloud_depth_current[neighborReprojectionCoord].g;
				
				float spatialWeight = Gaussian(length(float2(offset)), GAUSSIAN_SIGMA_SPATIAL);
				float rangeWeight = Gaussian(abs(tToDepthBuffer - cloudDepth), GAUSSIAN_SIGMA_RANGE);
				float weight = spatialWeight * rangeWeight;
				
				color += cloudResult * weight;
				weightSum += weight;
			}
		}

		if (weightSum > 0)
		{
			result = color / weightSum;
		}
	}

	// The water between the clouds and the eye. They are blended over the scene
	// after the sky, so without this they overwrite the sky's fog and leave a
	// bright unfogged band at the horizon when seen from under water.
	//
	// **A cloud is not under the water**, which is what the zero height says.
	// With the eye above the water there is nothing in front of it and this is a
	// no-op; with the eye below, the segment is fogged over its whole length,
	// which is right in the only case it is visible at all.
	//
	// It must be said explicitly, because the position handed over is the DEPTH
	// BUFFER's, not the cloud's. Letting that be measured against the surface
	// puts the cloud wherever the sea bed is, and flying over the ocean absorbed
	// every cloud in front of the water away to nothing.
	ApplyWaterFogPremultiplied(GetWaterFog(uv, depthWorldPosition, 0), result);

	// Sampled on the far plane, not at the depth buffer's position. A cloud is
	// kilometres away and the volume ends long before it, so the far plane and
	// the cloud read the same texel - the whole column the volume holds - while
	// the depth buffer's position would hand the cloud the light gathered as
	// far as whatever opaque thing happens to stand in front of it.
	ApplyVolumetricLightPremultiplied(uv, reconstruct_position(uv, 0), result);

	return result;
}
