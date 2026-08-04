#ifndef WI_SHADERINTEROP_FONT_H
#define WI_SHADERINTEROP_FONT_H
#include "ShaderInterop.h"

static const uint FONT_FLAG_SDF_RENDERING = 1u << 0u;
static const uint FONT_FLAG_OUTPUT_COLOR_SPACE_HDR10_ST2084 = 1u << 1u;
static const uint FONT_FLAG_OUTPUT_COLOR_SPACE_LINEAR = 1u << 2u;
// Fog this draw with the ocean's water. Set only for text drawn into the
// scene: this same shader draws all UI text and the debug overlay, whose
// transform is a 2D canvas projection with no world position behind it.
static const uint FONT_FLAG_UNDERWATER_FOG = 1u << 3u;
// Keep only the half of this draw on one side of the ocean surface. The
// transparent pass issues scene text twice, once on each side of the water, so
// that text crossing the waterline is refracted below it and dry above it.
// Exactly one of the two is set, or neither.
static const uint FONT_FLAG_WATERSIDE_SUBMERGED = 1u << 4u;
static const uint FONT_FLAG_WATERSIDE_ABOVE = 1u << 5u;
// Clip this draw against the camera's clip plane, for a planar reflection. Set
// only for text drawn into the scene, which is the only text with a world
// position behind its transform.
static const uint FONT_FLAG_CLIP_PLANE = 1u << 6u;

struct FontVertex
{
	float2 pos;
	float2 uv;
};

namespace SDF
{
	static const uint padding = 5;
	static const uint onedge_value = 127;
	static const float onedge_value_unorm = float(onedge_value) / 255.0f;
	static const float pixel_dist_scale = float(onedge_value) / float(padding);
}
struct FontConstants
{
	int buffer_index;
	uint buffer_offset;
	int texture_index;
	int padding0;

	uint2 color; // packed half4
	uint2 softness_bolden_hdrscaling; // packed half3 | uint16 flags

	float4x4 transform;
};
CONSTANTBUFFER(font, FontConstants, CBSLOT_FONT);

#endif // WI_SHADERINTEROP_FONT_H
