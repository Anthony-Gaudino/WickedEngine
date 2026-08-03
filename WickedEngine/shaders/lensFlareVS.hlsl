#include "globals.hlsli"
#include "waterVolumetricsHF.hlsli"

PUSHCONSTANT(push, LensFlarePush);

Texture2D<float4> texture_occlusion : register(t0);

static const float2 BILLBOARD[] = {
	float2(-1, -1),
	float2(1, -1),
	float2(-1, 1),
	float2(1, 1),
};

struct VertexOut
{
	float4 pos : SV_POSITION;
	float2 uv : TEXCOORD0;
	nointerpolation float opacity : TEXCOORD1;
	nointerpolation float3 waterTransmittance : TEXCOORD2;
};

VertexOut main(uint vertexID : SV_VertexID)
{
	VertexOut Out;
	Out.uv = BILLBOARD[vertexID] * 0.5f + 0.5f;

	// determine the flare opacity based on depth buffer occlusion and occlusion mask:
	float referenceDepth = saturate(1 - push.xLensFlarePos.z);
	const float2 step = 1.0 / (GetInternalResolution() * push.xLensFlarePos.z);
	const float2 range = 10.5 * step;
	float samples = 0;
	float visibility = 0;
	for (float y = -range.y; y <= range.y; y += step.y)
	{
		for (float x = -range.x; x <= range.x; x += step.x)
		{
			samples++;
			visibility += texture_depth.SampleLevel(sampler_point_clamp, push.xLensFlarePos.xy + float2(x, y), 0).r <= referenceDepth ? 1 : 0;
		}
	}
	visibility /= samples;

	if (any(push.xLensFlareDirectionalLight) && GetScene().texture_cloudmap >= 0)
	{
		// Cloud mask directional occlusion:
		half3 V = -unpack_half3(push.xLensFlareDirectionalLight);
		float4 cloudmap = bindless_textures[descriptor_index(GetScene().texture_cloudmap)].SampleLevel(sampler_linear_clamp, encode_hemioct(V.xzy) * 0.5 + 0.5, 0);
		visibility *= 1 - cloudmap.a;
	}

	Out.opacity = visibility;

	float2 pos = (push.xLensFlarePos.xy - 0.5) * float2(2, -2);
	float2 moddedpos = pos * push.xLensFlareOffset;
	Out.opacity *= saturate(1 - length(pos - moddedpos));

	Out.pos = float4(moddedpos + BILLBOARD[vertexID] * push.xLensFlareSize * GetCamera().canvas_size_rcp, 0, 1);

	// What the water leaves of the light before it ever reaches the lens.
	//
	// A flare is made inside the lens, so being under water does not suppress
	// it - a submerged camera pointed at a lamp flares just as one in air does.
	// What changes is the light arriving to make it: dimmed, and shifted blue,
	// by however much water it crossed. A red lamp thirty metres off flares
	// faintly and blue, not brightly and red.
	//
	// Spectral, so it cannot ride on the scalar opacity above. Per flare rather
	// than per pixel - one light, one path, and the quad is a screen space
	// artefact rather than a place in the world.
	Out.waterTransmittance = 1;

	[branch]
	if (GetCamera().IsUnderwaterFog())
	{
		// The same submersion test the rest of the frame uses, taken where the
		// light appears on screen. Saturated because a directional light
		// projects to wherever the sun is, which may be off screen entirely.
		const WaterVolumetrics medium =
			GetWaterVolumetricsAtEye(saturate(push.xLensFlarePos.xy));

		[branch]
		if (medium.IsActive())
		{
			const float3 toLight =
				push.xLensFlareWorldPos - GetCamera().position;
			const float distanceToLight = length(toLight);

			// LightTransmittance clips the path at the surface and takes the
			// Fresnel loss when it had to, so a lamp above the water and one
			// below it are both handled - and the sun's flare dies out as it
			// drops, most of a low sun reflecting off the surface rather than
			// entering the water at all.
			Out.waterTransmittance = medium.LightTransmittance(
				GetCamera().position,
				toLight / max(distanceToLight, 0.00001),
				distanceToLight
			);
		}
	}

	return Out;
}
