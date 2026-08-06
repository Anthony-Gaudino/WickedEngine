#define DISABLE_SOFT_SHADOWMAP
#include "volumetricLightHF.hlsli"
#include "fogHF.hlsli"
#include "oceanSurfaceHF.hlsli"
#include "waterVolumetricsHF.hlsli"

float4 main(VertexToPixel input) : SV_TARGET
{
	//return float4(1,0,0,1);
	ShaderEntity light = load_entity(rectlights().first_item() + (uint)g_xColor.x);

	float2 ScreenCoord = input.pos2D.xy / input.pos2D.w * float2(0.5, -0.5) + 0.5;
	float4 depths = texture_depth.GatherRed(sampler_point_clamp, ScreenCoord);
	float depth = max(input.pos.z, max(depths.x, max(depths.y, max(depths.z, depths.w))));
	float3 P = reconstruct_position(ScreenCoord, depth);
	float3 nearP = GetCamera().frustum_corners.screen_to_nearplane(ScreenCoord);
	float3 V = nearP - P; // ortho support
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

	float3 rayEnd = nearP;
	
	const uint sampleCount = 16;
	const half sampleCount_rcp = rcp((half)sampleCount);
	const float stepSize = distance(rayEnd, P) / sampleCount;

	// dither ray start to help with undersampling:
	P = P + V * stepSize * dither(input.pos.xy);
	
	const uint maskTex = light.GetTextureIndex();
	
	const half4 quaternion = light.GetQuaternion();
	const half3 right = rotate_vector(half3(1, 0, 0), quaternion);
	const half3 up = rotate_vector(half3(0, 1, 0), quaternion);
	const half3 forward = cross(up, right);
	const half light_length = max(0.01, light.GetLength());
	const half light_height = max(0.01, light.GetHeight());
	const float3 p0 = light.position - right * light_length * 0.5 + up * light_height * 0.5;
	const float3 p1 = light.position + right * light_length * 0.5 + up * light_height * 0.5;
	const float3 p2 = light.position + right * light_length * 0.5 - up * light_height * 0.5;
	const float3 p3 = light.position - right * light_length * 0.5 - up * light_height * 0.5;

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
	//
	// The ray advances in the increment expression rather than at the foot of
	// the body, so that the "behind light" early-out below still moves it. As a
	// `continue` over an advance written at the end, that test froze P: every
	// remaining sample re-tested the same point, found it behind the light too,
	// and the rest of the segment - including whatever lay in front of the
	// light - was never marched at all.
	[loop]
	for (uint i = 0; i < sampleCount; ++i, marchedDistance += stepSize, P += V * stepSize)
	{
		float3 L = light.position - P;
		const float3 Lunnormalized = L;
		const half dist2 = dot(L, L);
		const half dist = sqrt(dist2);
		L /= dist;

		const half range = light.GetRange();
		const half range2 = range * range;
		half3 attenuation = attenuation_pointlight(dist2, range, light.GetRange2Rcp());
	
		if (dot(P - light.position, forward) <= 0)
			continue; // behind light
		
		// Solid angle based on the Frostbite presentation: Moving Frostbite to Physically Based Rendering by Sebastien Lagarde, Charles de Rousiers, Siggraph 2014
		float3 v0 = normalize(p0 - P);
		float3 v1 = normalize(p1 - P);
		float3 v2 = normalize(p2 - P);
		float3 v3 = normalize(p3 - P);
		float3 n0 = normalize(cross(v0, v1));
		float3 n1 = normalize(cross(v1, v2));
		float3 n2 = normalize(cross(v2, v3));
		float3 n3 = normalize(cross(v3, v0));
		float g0 = acos(dot(-n0, n1));
		float g1 = acos(dot(-n1, n2));
		float g2 = acos(dot(-n2, n3));
		float g3 = acos(dot(-n3, n0));
		const float solid_angle = saturate(g0 + g1 + g2 + g3 - 2 * PI);
	
		attenuation *= solid_angle * 0.25 * (
			saturate(dot(v0, L)) +
			saturate(dot(v1, L)) +
			saturate(dot(v2, L)) +
			saturate(dot(v3, L))
		);
		
		[branch]
		if (light.IsCastingShadow())
		{
			float4 shadow_pos = mul(load_entitymatrix(light.GetMatrixIndex() + 0), float4(P, 1));
			shadow_pos.xyz /= shadow_pos.w;
			float2 shadow_uv = clipspace_to_uv(shadow_pos.xy);
			[branch]
			if (is_saturated(shadow_uv))
			{
				attenuation *= shadow_2D(light, shadow_pos.z, shadow_uv.xy, 0, input.pos.xy);
			}
		}
			
		[branch]
		if (maskTex > 0)
		{
			float4 shadow_pos = mul(load_entitymatrix(light.GetMatrixIndex() + 0), float4(P, 1));
			shadow_pos.xyz /= shadow_pos.w;
			float2 shadow_uv = shadow_pos.xy * float2(0.5f, -0.5f) + float2(0.5f, 0.5f);
			half4 mask = bindless_textures_half4[descriptor_index(maskTex)].Sample(sampler_linear_clamp, shadow_uv);
			attenuation *= mask.rgb * mask.a;
		}

		[branch]
		if (water.IsActive())
		{
			// Physical single scattering, attenuated on both legs.
			attenuation *= (half3)(
				water.InScatter(P, L, V, dist)
				* water.ViewTransmittance(
					cameraDistance - marchedDistance - stepSize)
				* stepIntegralWeight
			);
		}
		else
		{
			// Evaluate sample height for exponential fog calculation, given 0 for V:
			attenuation *= g_xColor.y + GetFogAmount(cameraDistance - marchedDistance, P, 0);
			attenuation *= ComputeScattering(saturate(dot(L, -V)));
		}

		accumulation += attenuation;
	}

	// The water march integrates rather than averages, so the per-light boost
	// becomes a gain over the physical result, 0 meaning physical.
	accumulation *= water.IsActive()
		? (half)(1 + g_xColor.y)
		: sampleCount_rcp;

	return saturateMediump(max(0, half4(accumulation * light.GetColor().rgb, 1)));
}
