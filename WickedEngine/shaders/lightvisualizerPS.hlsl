/**
 * Point and spot light visualizers.
 *
 * Identical to `vertexcolorPS` except that it fogs itself with the water. It
 * exists only so that it can: `vertexcolorPS` is shared with the wireframe
 * overlay and the debug renderer, which are diagnostic and must stay legible
 * whatever they are drawn through.
 *
 * Drawn additively, so the water only takes light away here - the veiling light
 * is already in whatever the visualizer is drawn over.
 */

#include "globals.hlsli"
#include "waterFogHF.hlsli"

float4 main(float4 pos : SV_Position, half4 col : COLOR) : SV_Target
{
	const float2 screenUV = pos.xy * GetCamera().internal_resolution_rcp;
	const float3 P = reconstruct_position(screenUV, pos.z);

	// Half of this visualizer may belong on the far side of the water surface,
	// where the transparent pass issued it as a separate draw.
	ClipToWaterSide(
		P,
		GetCamera().IsWaterSideBeyond(),
		GetCamera().IsWaterSideNear()
	);

	half4 color = col;
	ApplyWaterFogAdditive(screenUV, P, color);

	return color;
}
