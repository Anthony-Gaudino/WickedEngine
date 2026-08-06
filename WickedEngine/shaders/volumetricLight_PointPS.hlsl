#define DISABLE_SOFT_SHADOWMAP
#include "volumetricLightHF.hlsli"
#include "fogHF.hlsli"
#include "oceanSurfaceHF.hlsli"
#include "waterVolumetricsHF.hlsli"

float4 main(VertexToPixel input) : SV_TARGET
{
	ShaderEntity light = load_entity(pointlights().first_item() + (uint)g_xColor.x);

	float2 ScreenCoord = input.pos2D.xy / input.pos2D.w * float2(0.5, -0.5) + 0.5;
	float4 depths = texture_depth.GatherRed(sampler_point_clamp, ScreenCoord);
	float depth = max(input.pos.z, max(depths.x, max(depths.y, max(depths.z, depths.w))));
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
		float dist = intersectPlaneClampInfiniteDist(waterEye, rayDirection, float3(0, 1, 0), ocean_drawn_surface_height(waterEye));
		if (dist > 0 && dist < cameraDistance)
		{
			P = waterEye + rayDirection * dist;
			cameraDistance = dist;
		}
	}

	float marchedDistance = 0;
	half3 accumulation = 0;

	float3 rayEnd = GetCamera().position;
	if (length(rayEnd - light.position) > light.GetRange())
	{
		// if we are outside the light volume, then rayEnd will be the traced sphere frontface:
		float t = trace_sphere(rayEnd, -V, light.position, light.GetRange());
		rayEnd = rayEnd - t * V;
	}

	const uint sampleCount = 16;
	const half sampleCount_rcp = rcp((half)sampleCount);
	const float stepSize = distance(rayEnd, P) / sampleCount;

	// dither ray start to help with undersampling:
	P = P + V * stepSize * dither(input.pos.xy);
	
	const uint maskTex = light.GetTextureIndex();

	// The whole marched segment is on one side of the surface, because the
	// clamp above ends an upward ray there and a ray heading down or along
	// never reaches it from below. One test, no straddling - the PIXEL's test
	// rather than the camera's, so this agrees with the other underwater
	// effects while the camera is crossing.
	const WaterVolumetrics water = GetWaterVolumetricsAtEye(ScreenCoord);

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
	for(uint i = 0; i < sampleCount; ++i)
	{
		float3 L = light.position - P;
		const float3 Lunnormalized = L;
		const half dist2 = dot(L, L);
		const half dist = sqrt(dist2);
		L /= dist;

		const half range = light.GetRange();
		const half range2 = range * range;
		half3 attenuation = attenuation_pointlight(dist2, range, light.GetRange2Rcp());

		[branch]
		if (light.IsCastingShadow())
		{
			attenuation *= shadow_cube(light, Lunnormalized, input.pos.xy);
		}
		
		[branch]
		if (maskTex > 0)
		{
			half4 mask = bindless_cubemaps_half4[descriptor_index(maskTex)].Sample(sampler_linear_clamp, -Lunnormalized);
			attenuation *= mask.rgb * mask.a;
		}

		[branch]
		if (water.IsActive())
		{
			// Physical single scattering: the beam is now attenuated on both
			// legs, so a submerged lamp loses its warm end within a metre or
			// two and fades with distance instead of carrying forever.
			attenuation *= (half3)(
				water.InScatter(P, L, V, dist)
				* water.ViewTransmittance(
					cameraDistance - marchedDistance - stepSize)
				* stepIntegralWeight
			);
		}
		else
		{
			// Evaluate sample height for height fog calculation, given 0 for V:
			attenuation *= (half)g_xColor.y + GetFogAmount(cameraDistance - marchedDistance, P, 0);
			attenuation *= ComputeScattering(saturate(dot(L, -V)));
		}

		accumulation += attenuation;

		marchedDistance += stepSize;
		P = P + V * stepSize;
	}

	// The water march integrates rather than averages, so the per-light boost
	// becomes a gain over the physical result, 0 meaning physical.
	accumulation *= water.IsActive()
		? (half)(1 + g_xColor.y)
		: sampleCount_rcp;

	return saturateMediump(max(0, half4(accumulation * light.GetColor().rgb, 1)));
}
