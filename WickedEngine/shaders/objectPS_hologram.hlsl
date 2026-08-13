#define OBJECTSHADER_LAYOUT_COMMON
#include "objectHF.hlsli"

[earlydepthstencil]
float4 main(PixelInput input) : SV_TARGET
{
	ShaderMaterial material = GetMaterial();

	float4 uvsets = input.GetUVSets();
	
	write_mipmap_feedback(push.materialIndex, ddx_coarse(uvsets), ddy_coarse(uvsets));
	
	float4 color;
	[branch]
	if (material.textures[BASECOLORMAP].IsValid() && (GetFrame().options & OPTION_BIT_DISABLE_ALBEDO_MAPS) == 0)
	{
		color = material.textures[BASECOLORMAP].Sample(sampler_objectshader, uvsets);
		color.rgb = max(color.r, max(color.g, color.b));
	}
	else
	{
		color = 1;
	}
	color *= input.color;

	float3 emissiveColor = material.GetEmissive();
	[branch]
	if (any(emissiveColor) && material.textures[EMISSIVEMAP].IsValid())
	{
		float4 emissiveMap = material.textures[EMISSIVEMAP].Sample(sampler_objectshader, uvsets);
		emissiveColor *= emissiveMap.rgb * emissiveMap.a;
	}
	color.rgb += emissiveColor;

	float time = GetFrame().time;
	float2 uv = input.pos.xy * GetCamera().internal_resolution_rcp;
	uv.y *= GetCamera().internal_resolution.y * GetCamera().internal_resolution_rcp.x;

	// wave:
	color.a *= sin(input.GetPos3D().y * 30 + time * 10) * 0.5 + 0.5;

	// rim:
	color *= lerp(0.3, 6, pow(1 - saturate(dot(normalize(input.nor), normalize(input.GetViewVector()))), 2));

	// keep some base color
	color.a += 0.2;
	color.a = saturate(color.a);

	// grain:
	float noise = 1;
	noise *= texture_random64x64.SampleLevel(sampler_linear_wrap, uv * 4 + time, 0).r * 2 - 1;
	noise *= texture_random64x64.SampleLevel(sampler_linear_wrap, uv * float2(-2, 3) + time, 0).g * 2 - 1;
	noise = noise * 0.5 + 0.5;

	color *= noise;

	// The water between this fragment and the eye. Its own main(), so it never
	// reached the fog that objectHF applies for every other object shader - a
	// hologram is stylised, but it is still something placed in the scene
	// rather than a diagnostic overlay, so the water reaches it like anything
	// else.
	// NOT the uv above: that one has had its y scaled by the aspect ratio for
	// the grain, so it is no longer a screen coordinate.
	const float2 screenUV = input.pos.xy * GetCamera().internal_resolution_rcp;

	half4 fogged = (half4)color;
	[branch]
	if (material.IsAdditive())
	{
		ApplyWaterFogAdditive(screenUV, input.GetPos3D(), fogged);
	}
	else
	{
		ApplyWaterFog(screenUV, input.GetPos3D(), fogged);
		ApplyVolumetricLight(screenUV, input.GetPos3D(), fogged);
	}

	return saturateMediump((float4)fogged);
}
