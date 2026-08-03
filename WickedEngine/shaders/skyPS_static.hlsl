#include "skyHF.hlsli"
#include "waterFogHF.hlsli"

struct PushData
{
	uint clouds_enabled;
};
PUSHCONSTANT(push, PushData);

float4 main(float4 pos : SV_Position) : SV_Target
{
	const float3 V = normalize(-GetCamera().screen_to_view(pos));

	float4 color = float4(GetStaticSkyColor(V, push.clouds_enabled), 1);

	// The water this ray crossed before it left the surface. The sky covers
	// every pixel no geometry wrote, so fogging it here is what gives the
	// "infinite water" backdrop that the post pass used to paint from the far
	// plane. Composited by hand rather than through ApplyWaterFog because the
	// colour here is full precision and that applier works in half.
	const float2 waterFogUV = pos.xy * GetCamera().internal_resolution_rcp;
	const WaterFog waterFog =
		GetWaterFog(waterFogUV, reconstruct_position(waterFogUV, pos.z));
	color.rgb = color.rgb * waterFog.transmittance + waterFog.inscatter;

	color = saturateMediump(color);
	return color;
}

