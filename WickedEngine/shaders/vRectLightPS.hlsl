#include "globals.hlsli"
#include "waterFogHF.hlsli"

float4 main(float4 pos : SV_Position, float2 uv : TEXCOORD) : SV_TARGET
{
	float4 color = float4(xLightColor.rgb * xLightEnerdis.x, 1);

	int maskTexDescriptor = (int)xLightEnerdis.y;
	if (maskTexDescriptor >= 0)
	{
		color *= bindless_textures[descriptor_index(maskTexDescriptor)].Sample(sampler_linear_clamp, uv);
	}

	// The water between this visualizer and the eye. Without it the gizmo for a
	// submerged light stays at full brightness however far away it is.
	const float2 screenUV = pos.xy * GetCamera().internal_resolution_rcp;
	half4 fogged = (half4)color;
	ApplyWaterFog(screenUV, reconstruct_position(screenUV, pos.z), fogged);

	return (float4)fogged;
}
