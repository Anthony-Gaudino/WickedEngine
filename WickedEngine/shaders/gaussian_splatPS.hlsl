#include "globals.hlsli"
#include "ShaderInterop_GaussianSplat.h"
#include "waterFogHF.hlsli"

float4 main(float4 pos : SV_Position, half4 color : COLOR, half2 localPos : LOCALPOS) : SV_Target
{
	const float A = dot(localPos, localPos);
	const half opacity = exp(-0.5 * A) * color.a;
	clip(opacity - 1.0 / 255.0);

	half4 result = color * opacity;

	// The water between this splat and the eye. Premultiplied, matching the
	// blend. The whole billboard takes its centre's depth - reconstruct_position
	// builds its own clip vector and never reads SV_Position.w, which the vertex
	// shader hardcodes to 1, so pos.z is usable as it stands.
	const float2 screenUV = pos.xy * GetCamera().internal_resolution_rcp;
	ApplyWaterFogPremultiplied(
		screenUV, reconstruct_position(screenUV, pos.z), result);

	return result;
}
