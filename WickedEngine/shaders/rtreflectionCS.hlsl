#define RTAPI
#define DISABLE_SOFT_SHADOWMAP
#define DISABLE_TRANSPARENT_SHADOWMAP
#define SURFACE_LOAD_MIPCONE
#define TEXTURE_SLOT_NONUNIFORM

#include "globals.hlsli"
#include "ShaderInterop_Postprocess.h"
#include "raytracingHF.hlsli"
#include "stochasticSSRHF.hlsli"
#include "lightingHF.hlsli"
#include "ShaderInterop_DDGI.h"
#include "rtsceneHF.hlsli"

PUSHCONSTANT(postprocess, PostProcess);

RWTexture2D<float4> output_rayIndirectSpecular : register(u0);
RWTexture2D<float4> output_rayDirectionPDF : register(u1);
RWTexture2D<float> output_rayLengths : register(u2);

struct RayPayload
{
	float4 data;
};

[numthreads(8, 4, 1)]
void main(uint2 DTid : SV_DispatchThreadID)
{
	const float2 uv = ((float2)DTid.xy + 0.5) * postprocess.resolution_rcp;

	const float downsampleFactor = rtreflection_downscalefactor;

	// This is necessary for accurate upscaling. This is so we don't reuse the same half-res pixels
	uint2 screenJitter = floor(blue_noise(uint2(0, 0)).xy * downsampleFactor);
	uint2 jitterPixel = screenJitter + DTid.xy * downsampleFactor;
	float2 jitterUV = (screenJitter + DTid.xy + 0.5f) * postprocess.resolution_rcp;

	const float depth = texture_depth.SampleLevel(sampler_linear_clamp, jitterUV, 0);
	const float lineardepth = texture_lineardepth.SampleLevel(sampler_linear_clamp, jitterUV, 0);

	const half3 normal_roughness = texture_normal_roughness[jitterPixel].rgb;
	const half roughness = normal_roughness.b;

	if (!NeedReflection(roughness, depth, rtreflection_roughness_cutoff))
	{
		output_rayIndirectSpecular[DTid.xy] = 0;
		output_rayDirectionPDF[DTid.xy] = 0;
		output_rayLengths[DTid.xy] = FLT_MAX;
		return;
	}

	const float3 N = decode_normal(normal_roughness.rg);
	const float3 P = reconstruct_position(jitterUV, depth);
	const float3 V = normalize(GetCamera().frustum_corners.screen_to_nearplane(uv) - P); // ortho support

	const float4 GGX = ReflectionDir_GGX(V, N, roughness, blue_noise(DTid.xy).xy);
	const float3 R = GGX.xyz;
	const float PDF = GGX.w;

	RayDesc ray;
	ray.TMin = 0.01;
	ray.TMax = rtreflection_range;
	ray.Origin = P;
	ray.Direction = normalize(R);

	RayPayload payload;
	payload.data = 0;

	const float minraycone = 0.05;
	RayCone raycone = RayCone::from_spread_angle(pixel_cone_spread_angle_from_image_height(postprocess.resolution.y));
	raycone = raycone.propagate(sqr(max(minraycone, roughness)), lineardepth * GetCamera().z_far);

	// Trace + shade the reflection ray through the shared scene-radiance
	// gather:
	SceneRadiance rad = TraceSceneRadiance(
		ray,
		asuint(postprocess.params1.x),
		RAY_FLAG_SKIP_PROCEDURAL_PRIMITIVES |
		RAY_FLAG_CULL_BACK_FACING_TRIANGLES,
		raycone,
		DTid.xy,
		postprocess.resolution_rcp
	);
	payload.data.xyz = rad.color;
	payload.data.w = rad.distance;

	output_rayIndirectSpecular[DTid.xy] = float4(payload.data.xyz, 1);
	output_rayDirectionPDF[DTid.xy] = float4(R, PDF);
	output_rayLengths[DTid.xy] = payload.data.w;
}
