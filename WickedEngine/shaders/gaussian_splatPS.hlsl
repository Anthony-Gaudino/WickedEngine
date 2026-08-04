#include "globals.hlsli"
#include "ShaderInterop_GaussianSplat.h"
#include "waterFogHF.hlsli"

float4 main(float4 pos : SV_Position, half4 color : COLOR, half2 localPos : LOCALPOS) : SV_Target
{
	const float A = dot(localPos, localPos);
	const half opacity = exp(-0.5 * A) * color.a;
	clip(opacity - 1.0 / 255.0);

	// The whole billboard takes its centre's depth - reconstruct_position
	// builds its own clip vector and never reads SV_Position.w, which the
	// vertex shader hardcodes to 1, so pos.z is usable as it stands. The
	// position still varies across the billboard, tracing a plane at that one
	// depth, which is what lets a splat straddling the waterline be cut rather
	// than snapped.
	const float2 screenUV = pos.xy * GetCamera().internal_resolution_rcp;
	const float3 P = reconstruct_position(screenUV, pos.z);

	// Half of this model may belong on the far side of the water surface, where
	// the transparent pass issued it as a separate draw.
	ClipToWaterSide(
		P,
		GetCamera().IsWaterSideSubmerged(),
		GetCamera().IsWaterSideAbove()
	);

	half4 result = color * opacity;

	// The water between this splat and the eye. Premultiplied, matching the
	// blend.
	ApplyWaterFogPremultiplied(screenUV, P, result);

	return result;
}
