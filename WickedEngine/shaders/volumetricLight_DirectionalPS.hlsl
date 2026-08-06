#define DISABLE_SOFT_SHADOWMAP
#include "volumetricLightHF.hlsli"
#include "volumetricCloudsHF.hlsli"
#include "fogHF.hlsli"
#include "oceanSurfaceHF.hlsli"
#include "waterVolumetricsHF.hlsli"

float4 main(VertexToPixel input) : SV_Target
{
	ShaderEntity light = load_entity(directional_lights().first_item() + (uint)g_xColor.x);
	
	float2 ScreenCoord = input.pos2D.xy / input.pos2D.w * float2(0.5, -0.5) + 0.5;
	float4 depths = texture_depth.GatherRed(sampler_point_clamp, ScreenCoord);
	float depth = max(depths.x, max(depths.y, max(depths.z, depths.w)));
	float3 P = reconstruct_position(ScreenCoord, depth);
	float3 V = GetCamera().frustum_corners.screen_to_nearplane(ScreenCoord) - P; // ortho support
	float cameraDistance = length(V);
	V /= cameraDistance;

	// Fix for ocean: because ocean is not in linear depth, we trace it instead
	const ShaderOcean ocean = GetWeather().ocean;
	if (ocean.IsValid())
	{
		// Traced from this pixel's own submersion test point against the
		// wave displaced surface, matching GetWaterVolumetricsAtEye below.
		// Using the camera and the flat plane here would clamp a downward ray
		// to nothing for the very pixels that test as submerged, so the medium
		// would be active over a segment too short to scatter anything.
		const float3 waterEye = ocean_underwater_test_position(ScreenCoord);
		float3 rayDirection = -V;
		float dist = intersectPlaneClampInfiniteDist(waterEye, rayDirection, float3(0, 1, 0), ocean_surface_height(waterEye));
		if (dist > 0 && dist < cameraDistance)
		{
			P = waterEye + rayDirection * dist;
			cameraDistance = dist;
		}
	}

	float marchedDistance = 0;
	half3 accumulation = 0;

	const half3 L = light.GetDirection();
	const half scattering = ComputeScattering(saturate(dot(L, -V)));

	// The segment is entirely on one side of the surface, because the clamp
	// above already ends an upward ray there and a ray heading down or along
	// never reaches it from below. So one test decides the medium for the whole
	// loop, with no straddling - but it has to be the PIXEL's test, not the
	// camera's, or this disagrees with every other underwater effect while the
	// camera is crossing.
	const WaterVolumetrics water = GetWaterVolumetricsAtEye(ScreenCoord);

	// Sunlight bends at the surface before it ever reaches a submerged sample.
	// Shared with surface shading, so a shaft and the sea bed it lands on agree
	// about where the sun is.
	const float3 Lwater = water.IsActive()
		? RefractIntoWater((float3)L)
		: (float3)L;

	// Cut the segment down to the distance this water is transparent over.
	//
	// The water branch INTEGRATES - every sample carries a stepSize - so its
	// accuracy rests on the near end of the segment being sampled finely. All
	// the in-scatter that reaches the eye comes from the last few metres, yet
	// the closest sample sits a whole step away. Let the segment run to the far
	// plane, which is what happens over open water with nothing behind it to
	// put in the depth buffer, and that step becomes hundreds of metres: every
	// sample lands where the water is already opaque. A thousand metre segment
	// returns three percent of the true integral in blue and none of the red.
	//
	// Beyond a few optical depths there is nothing left out there to gather, so
	// discarding that stretch costs nothing and buys back a step measured in
	// metres. Safe to move P here because this pass marches from the camera, so
	// shortening the far end cannot cross the near one.
	[branch]
	if (water.IsActive())
	{
		const float visibleRange = water.VisibleRange();
		if (cameraDistance > visibleRange)
		{
			P = GetCamera().position - V * visibleRange;
			cameraDistance = visibleRange;
		}
	}

	float3 rayEnd = GetCamera().position;

	const uint sampleCount = 16;
	const half sampleCount_rcp = rcp((half)sampleCount);
	const float stepSize = length(P - rayEnd) / sampleCount;

	// dither ray start to help with undersampling:
	P = P + V * min(stepSize * dither(input.pos.xy), 10); // limit dithering step to 10 to avoid very large dither on sky

	// Weight of one march step's in-scatter, hoisted out of the loop because it
	// depends only on the step. Guarded, because above the waterline the
	// medium's coefficients are zero and this divides by them.
	float3 stepIntegralWeight = 0;
	[branch]
	if (water.IsActive())
	{
		stepIntegralWeight = water.StepIntegralWeight(stepSize);
	}

	// Perform ray marching to integrate light volume along view ray:
	[loop]
	for (uint i = 0; i < sampleCount; ++i)
	{
		bool valid = false;

		half3 shadow = 1;
		// A submerged sample is occluded above the surface, on the part of the
		// path that is still straight, so the cascades - which know nothing of
		// refraction - are exactly right when sampled at the point the light
		// enters the water. That is also where the ocean writes its caustics
		// into the transparent shadow layer, so they land where they form.
		const float3 shadowSample = water.IsActive()
			? water.SurfaceEntry(P, Lwater, FLT_MAX)
			: P;
		for (uint cascade = 0; cascade < light.GetShadowCascadeCount(); ++cascade)
		{
			float3 shadow_pos = mul(load_entitymatrix(light.GetMatrixIndex() + cascade), float4(shadowSample, 1)).xyz; // ortho matrix, no divide by .w
			float3 shadow_uv = clipspace_to_uv(shadow_pos.xyz);

			[branch]
			if (is_saturated(shadow_uv))
			{
				shadow *= shadow_2D(light, shadow_pos.z, shadow_uv.xy, cascade, input.pos.xy);
				break;
			}
		}

		if (GetFrame().options & OPTION_BIT_VOLUMETRICCLOUDS_CAST_SHADOW)
		{
			shadow *= shadow_2D_volumetricclouds(P);
		}

		[branch]
		if (water.IsActive())
		{
			// Physical single scattering, so this integrates over the segment
			// rather than averaging across it: the beam now scales with how
			// much water the ray actually crosses, and with the medium's own
			// scattering coefficient instead of an authored fog density.
			shadow *= (half3)(
				water.InScatter(P, Lwater, V, FLT_MAX)
				* water.ViewTransmittance(
					cameraDistance - marchedDistance - stepSize)
				* stepIntegralWeight
			);
		}
		else
		{
			// Evaluate sample height for height fog calculation, given 0 for V:
			shadow *= (half)g_xColor.y + GetFogAmount(cameraDistance - marchedDistance, P, 0);
			shadow *= scattering;
		}

		accumulation += shadow;

		marchedDistance += stepSize;
		P = P + V * stepSize;
	}
	// The water march integrates, so it must not be divided by the sample
	// count. That turns the per-light boost from a density floor for the fog to
	// scatter from into a gain over the physical result, 0 meaning physical.
	accumulation *= water.IsActive()
		? (half)(1 + g_xColor.y)
		: sampleCount_rcp;

	half3 atmosphere_transmittance = 1;
	if (GetFrame().options & OPTION_BIT_REALISTIC_SKY)
	{
		atmosphere_transmittance = GetAtmosphericLightTransmittance(GetWeather().atmosphere, P, L, texture_transmittancelut);
	}

	return saturateMediump(max(0, half4(accumulation * light.GetColor().rgb * atmosphere_transmittance, 1)));
}
