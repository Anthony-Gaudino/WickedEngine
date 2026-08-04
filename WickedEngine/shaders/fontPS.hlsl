#include "globals.hlsli"
#include "ShaderInterop_Font.h"
#include "waterFogHF.hlsli"

struct VertextoPixel
{
	float4 pos : SV_Position;
	float2 uv : TEXCOORD0;
	float2 bary : TEXCOORD1;
};

float4 main(VertextoPixel input) : SV_TARGET
{
	Texture2D<half4> tex = bindless_textures_half4[descriptor_index(font.texture_index)];
	half value = tex.SampleLevel(sampler_linear_clamp, input.uv, 0).r;
	half4 color = unpack_half4(font.color);
	
	const half3 softness_bolden_hdrscaling = unpack_half3(font.softness_bolden_hdrscaling);
	const half softness = softness_bolden_hdrscaling.x;
	const half bolden = softness_bolden_hdrscaling.y;
	const half hdr_scaling = softness_bolden_hdrscaling.z;
	const min16uint flags = font.softness_bolden_hdrscaling.y >> 16u;

	// Per pixel clipping, for the cases where this shader is drawing text that
	// lives in the world. Both share one world position.
	[branch]
	if (flags & (
		FONT_FLAG_CLIP_PLANE |
		FONT_FLAG_WATERSIDE_SUBMERGED |
		FONT_FLAG_WATERSIDE_ABOVE))
	{
		const float2 screenUV = input.pos.xy * GetCamera().internal_resolution_rcp;
		const float3 P = reconstruct_position(screenUV, input.pos.z);

		// A planar reflection camera carries a clip plane, and its view must
		// not show what is behind the mirror. fontVS has no world position to
		// output a clip distance from, so it is tested here - which also means
		// text crossing the water has exactly its dry part reflected.
		[branch]
		if (flags & FONT_FLAG_CLIP_PLANE)
		{
			clip(dot(float4(P, 1), GetCamera().clip_plane));
		}

		// ...and half of it may belong on the far side of the water surface,
		// where the transparent pass issued it as a separate draw.
		ClipToWaterSide(
			P,
			(flags & FONT_FLAG_WATERSIDE_SUBMERGED) != 0,
			(flags & FONT_FLAG_WATERSIDE_ABOVE) != 0
		);
	}

	[branch]
	if (flags & FONT_FLAG_SDF_RENDERING)
	{
		float2 bary_fw = fwidth(input.bary);
		half w = max(bary_fw.x, bary_fw.y); // screen coverage dependency
		w = max(w, 1.0 / 255.0); // min softness to avoid pixelated hard edge in magnification
		w += softness;
		w = saturate(w);
		half mid = lerp((half)SDF::onedge_value_unorm, 0, bolden);
		color.a *= smoothstep(saturate(mid - w), saturate(mid + w), value);
	}
	else
	{
		color.a *= value;
	}

	color.rgb *= color.a; // NOTE: premultiplied blending! This is important for blending in HDR linear space

	// The water between this glyph and the eye. Premultiplied, matching the
	// line above: the veiling light has to be scaled by the coverage too, or a
	// glyph's antialiased edge would paint haze where it barely covers.
	[branch]
	if (flags & FONT_FLAG_UNDERWATER_FOG)
	{
		const float2 screenUV = input.pos.xy * GetCamera().internal_resolution_rcp;
		ApplyWaterFogPremultiplied(
			screenUV, reconstruct_position(screenUV, input.pos.z), color);
	}

	[branch]
	if (flags & FONT_FLAG_OUTPUT_COLOR_SPACE_LINEAR)
	{
		color.rgb = RemoveSRGBCurve_Fast(color.rgb);
		color.rgb *= hdr_scaling;
	}
	
	[branch]
	if (flags & FONT_FLAG_OUTPUT_COLOR_SPACE_HDR10_ST2084)
	{
		// https://github.com/microsoft/DirectX-Graphics-Samples/blob/master/Samples/Desktop/D3D12HDR/src/presentPS.hlsl
		const half referenceWhiteNits = 80.0;
		const half st2084max = 10000.0;
		const half hdrScalar = referenceWhiteNits / st2084max;
		// The input is in Rec.709, but the display is Rec.2020
		color.rgb = REC709toREC2020(color.rgb);
		// Apply the ST.2084 curve to the result.
		color.rgb = ApplyREC2084Curve(color.rgb * hdrScalar);
	}

	return color;
}
