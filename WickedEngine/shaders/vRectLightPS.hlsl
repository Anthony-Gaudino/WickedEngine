#include "globals.hlsli"
#include "waterFogHF.hlsli"
#include "volumetricFroxelHF.hlsli"

float4 main(float4 pos : SV_Position, float2 uv : TEXCOORD) : SV_TARGET
{
	float4 color = float4(xLightColor.rgb * xLightEnerdis.x, 1);

	int maskTexDescriptor = (int)xLightEnerdis.y;
	if (maskTexDescriptor >= 0)
	{
		color *= bindless_textures[descriptor_index(maskTexDescriptor)].Sample(sampler_linear_clamp, uv);
	}

	const float2 screenUV = pos.xy * GetCamera().internal_resolution_rcp;
	const float3 P = reconstruct_position(screenUV, pos.z);

	// Half of this visualizer may belong on the far side of the water surface,
	// where the transparent pass issued it as a separate draw.
	ClipToWaterSide(
		P,
		GetCamera().IsWaterSideBeyond(),
		GetCamera().IsWaterSideNear()
	);

	// The water between this visualizer and the eye. Without it the gizmo for a
	// submerged light stays at full brightness however far away it is.
	half4 fogged = (half4)color;
	ApplyWaterFog(screenUV, P, fogged);
	ApplyVolumetricLight(screenUV, P, fogged);

	return (float4)fogged;
}
