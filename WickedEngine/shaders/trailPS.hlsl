#include "globals.hlsli"
#include "waterFogHF.hlsli"

float4 main(float4 pos : SV_Position, float4 screen : SCREEN, float4 uv : TEXCOORD, float4 color : COLOR) : SV_TARGET
{
	color *= bindless_textures[descriptor_index(g_xTrailTextureIndex1)].Sample(sampler_linear_mirror, uv.xy);
	color *= bindless_textures[descriptor_index(g_xTrailTextureIndex2)].Sample(sampler_linear_mirror, uv.zw);

	[branch]
	if (g_xTrailDepthTextureIndex >= 0)
	{
		Texture2D depthTexture = bindless_textures[descriptor_index(g_xTrailDepthTextureIndex)];
		float2 screenUV = clipspace_to_uv(screen.xy / screen.w);
		float depthScene = compute_lineardepth(depthTexture.SampleLevel(sampler_linear_clamp, screenUV, 0).r, 0.1, g_xTrailCameraFar);
		float depthFragment = pos.w;
		color.a *= saturate(g_xTrailDepthSoften * (depthScene - depthFragment)); // soft depth fade
	}

	// The water between this trail and the eye. Without it a submerged trail
	// stays at full brightness however far away it is, and a trail crossing the
	// waterline shows no change at all as it goes under.
	//
	// Which form the fog takes depends on how this draw is combined with the
	// screen, which only the CPU knows - one pixel shader serves every blend
	// state. See TRAIL_WATERFOG_* in ShaderInterop_Renderer.h.
	[branch]
	if (g_xTrailWaterFogMode != TRAIL_WATERFOG_NONE)
	{
		const float2 fogUV = pos.xy * GetCamera().internal_resolution_rcp;
		const float3 P = reconstruct_position(fogUV, pos.z);

		half4 fogged = (half4)color;

		if (g_xTrailWaterFogMode == TRAIL_WATERFOG_ADDITIVE)
		{
			ApplyWaterFogAdditive(fogUV, P, fogged);
		}
		else if (g_xTrailWaterFogMode == TRAIL_WATERFOG_PREMULTIPLIED)
		{
			ApplyWaterFogPremultiplied(fogUV, P, fogged);
		}
		else
		{
			ApplyWaterFog(fogUV, P, fogged);
		}

		color = (float4)fogged;
	}

	return color;
}
