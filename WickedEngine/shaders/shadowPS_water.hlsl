#define OBJECTSHADER_LAYOUT_SHADOW_TEX
#define OBJECTSHADER_USE_COLOR
#include "objectHF.hlsli"

[earlydepthstencil]
float4 main(PixelInput input) : SV_TARGET
{
	ShaderMaterial material = GetMaterial();
	
	float4 uvsets = input.GetUVSets();

	half4 color;
	[branch]
	if (material.textures[BASECOLORMAP].IsValid())
	{
		color = material.textures[BASECOLORMAP].Sample(sampler_objectshader, uvsets);
	}
	else
	{
		color = 1;
	}
	color *= input.color;

	half opacity = color.a;

	color.rgb = 1; // disable water shadow because it has already fog

	const ShaderOcean ocean = GetWeather().ocean;
	if (ocean.texture_displacementmap >= 0)
	{
		float2 caustics_uv = uvsets.xy * ocean.caustics_scale;
		color.rgb += texture_caustics.SampleLevel(sampler_linear_mirror, caustics_uv, 0).rgb * ocean.caustics_intensity;
	}

	color.a = input.pos.z; // secondary depth

	return color;
}
