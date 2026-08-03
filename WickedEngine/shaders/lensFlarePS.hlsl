#include "globals.hlsli"

Texture2D<float4> texture_flare : register(t1);

struct VertexOut
{
	float4 pos : SV_POSITION;
	float2 uv : TEXCOORD0;
	nointerpolation float opacity : TEXCOORD1;
	nointerpolation float3 waterTransmittance : TEXCOORD2;
};

float4 main(VertexOut input) : SV_TARGET
{
	float4 color = texture_flare.SampleLevel(sampler_linear_clamp, input.uv, 0);
	color.a *= input.opacity;

	// Extinction only: the flare is drawn additively, so the water's veiling
	// light is already in whatever it is drawn over.
	color.rgb *= input.waterTransmittance;

	return color;
}
