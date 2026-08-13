#define TRANSPARENT // uses transparent light lists
#define LIGHTING_SCATTER
#define DISABLE_TRANSPARENT_SHADOWMAP // particle casts transparent shadow layer, but particle self-shadow doesn't look good, there is mismatch due to camera-facingness
#include "globals.hlsli"
#include "emittedparticleHF.hlsli"
#include "ShaderInterop_EmittedParticle.h"
#include "objectHF.hlsli"

#ifdef EMITTEDPARTICLE_DISTORTION
static const uint SLOT = NORMALMAP;
#else
static const uint SLOT = BASECOLORMAP;
#endif // EMITTEDPARTICLE_DISTORTION

[earlydepthstencil]
float4 main(VertextoPixel input) : SV_TARGET
{
	float3 pos3D = input.GetPos3D();

	// Half of this emitter may belong on the far side of the water surface,
	// where the transparent pass issued it as a separate draw. Discard that
	// half up front. Per fragment rather than per emitter: an emitter is a
	// volume, so classifying the whole of it by its centre would put spray
	// straddling the waterline entirely under the water or entirely clear of
	// it.
	ClipToWaterSide(
		pos3D,
		GetCamera().IsWaterSideSubmerged(),
		GetCamera().IsWaterSideAbove()
	);

	// Blocker shadow map check:
	[branch]
	if ((xEmitterOptions & EMITTER_OPTION_BIT_USE_RAIN_BLOCKER) && rain_blocker_check(pos3D))
	{
		return 0;
	}

	ShaderCamera camera = GetCamera();
	float2 ScreenCoord = input.pos.xy * camera.internal_resolution_rcp; // use pixel center!
	min16uint2 pixel = input.pos.xy; // no longer pixel center!
	
	write_mipmap_feedback(EmitterGetGeometry().materialIndex, ddx_coarse(input.tex.xyxy), ddy_coarse(input.tex.xyxy));
	
	ShaderMaterial material = EmitterGetMaterial();

	half4 color = 1;

	[branch]
	if (material.textures[SLOT].IsValid())
	{
		color = material.textures[SLOT].Sample(sampler_linear_clamp, input.tex.xyxy);

		[branch]
		if (xEmitterOptions & EMITTER_OPTION_BIT_FRAME_BLENDING_ENABLED)
		{
			half4 color2 = material.textures[SLOT].Sample(sampler_linear_clamp, input.tex.zwzw);
			color = lerp(color, color2, input.frameBlend);
		}
	}
	
	half4 inputColor;
	inputColor.r = ((input.color >> 0)  & 0xFF) / 255.0f;
	inputColor.g = ((input.color >> 8)  & 0xFF) / 255.0f;
	inputColor.b = ((input.color >> 16) & 0xFF) / 255.0f;
	inputColor.a = ((input.color >> 24) & 0xFF) / 255.0f;

	half opacity = color.a * inputColor.a;

	half3 normal = 0;
#ifndef EMITTEDPARTICLE_DISTORTION // the "distortion" shader is just using normal map as the color map and uses it for signed blending, this normal logic won't be used for that
	[branch]
	if (material.textures[NORMALMAP].IsValid())
	{
		normal = material.textures[NORMALMAP].Sample(sampler_linear_clamp, input.tex.xyxy).rgb;

		[branch]
		if (xEmitterOptions & EMITTER_OPTION_BIT_FRAME_BLENDING_ENABLED)
		{
			half3 normal2 = material.textures[NORMALMAP].Sample(sampler_linear_clamp, input.tex.zwzw).rgb;
			normal = lerp(normal, normal2, input.frameBlend);
		}

		normal -= 0.5;
		normal *= material.GetNormalMapStrength();
		normal *= opacity;
	}
#endif // EMITTEDPARTICLE_DISTORTION
	
	[branch]
	if (camera.texture_depth_index >= 0)
	{
		float4 depthScene = compute_lineardepth(texture_depth.GatherRed(sampler_linear_clamp, ScreenCoord));
		float depthFragment = input.pos.w;
		opacity *= saturate(1.0 / input.size * (max(max(depthScene.x, depthScene.y), max(depthScene.z, depthScene.w)) - depthFragment));
	}

	opacity = saturate(opacity);

	color.rgb *= inputColor.rgb * (1 + material.GetEmissive());
	color.a = opacity;

#ifdef EMITTEDPARTICLE_DISTORTION
	// just make normal maps blendable:
	color.rgb = color.rgb - 0.5;
#endif // EMITTEDPARTICLE_DISTORTION

#ifdef EMITTEDPARTICLE_LIGHTING

	[branch]
	if (color.a > 0)
	{
		float3 N;
		N.x = -cos(PI * input.unrotated_uv.x);
		N.y = cos(PI * input.unrotated_uv.y);
		N.z = -sin(PI * length(input.unrotated_uv));
		N.xz += normal.rg;
		N = mul((float3x3)camera.inverse_view, N);
		N = normalize(N);
		
		float3 V = input.GetViewVector();
		float dist = length(V);
		V /= dist;

		Lighting lighting;
		lighting.create(0, 0, GetAmbient(N), 0);

		Surface surface;
		surface.init();
		surface.create(material, color, surfacemap_simple);
		surface.P = pos3D;
		surface.N = N;
		surface.V = V;
		surface.pixel = pixel;
		surface.sss = material.GetSSS();
		surface.sss_inv = material.GetSSSInverse();
		surface.extinction = 0;
		surface.update();

		TiledLighting(surface, lighting, GetFlatTileIndex(pixel, camera), camera);

		color.rgb *= lighting.direct.diffuse + lighting.indirect.diffuse;
		color.rgb += lighting.indirect.specular;

		//color.rgb = float3(unrotated_uv, 0);
		//color.rgb = float3(input.tex, 0);

		color = max(0, color);
	}

#endif // EMITTEDPARTICLE_LIGHTING

#ifndef EMITTEDPARTICLE_DISTORTION
	// Height fog.
	//
	// Outside the lighting block on purpose. It used to sit inside it, and the
	// DEFAULT emitter type is SOFT, which does not define EMITTEDPARTICLE_
	// LIGHTING - so an ordinary emitter was never fogged at all, above water or
	// below. Skipped for the distortion permutation, which writes a normal map
	// into the distortion buffer rather than radiance into the scene.
	const float3 toEye = input.GetViewVector();
	const float distToEye = length(toEye);
	ApplyFog(distToEye, toEye / max(distToEye, 0.00001), color);

	// The water between this particle and the eye.
	//
	// Outside the lighting block on purpose: the default emitter type is SOFT,
	// which does not define it, so a particle would otherwise never be fogged
	// at all - and it is a particle's own distance that has to fog it, which is
	// the whole reason this moved out of the post pass.
	//
	// An additive particle takes the extinction only. The water's veiling light
	// is already in the destination, put there by whatever drew it, so adding
	// it again would thicken the haze once per additive particle stacked on the
	// pixel. An alpha blended one does take it, or a distant particle would
	// fade towards black instead of towards the colour of the water.
	[branch]
	if (material.IsAdditive())
	{
		ApplyWaterFogAdditive(ScreenCoord, pos3D, color);
	}
	else
	{
		ApplyWaterFog(ScreenCoord, pos3D, color);

		// Straight alpha here, so the blend applies the coverage itself and
		// this must not apply it again.
		ApplyVolumetricLight(ScreenCoord, pos3D, color);
	}
#endif // EMITTEDPARTICLE_DISTORTION

	color.rgb = mul(saturationMatrix(material.GetSaturation()), color.rgb);

	return color;
}
