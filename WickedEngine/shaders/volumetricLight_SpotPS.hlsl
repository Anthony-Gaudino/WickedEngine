#define DISABLE_SOFT_SHADOWMAP
#include "volumetricLightHF.hlsli"
#include "fogHF.hlsli"
#include "oceanSurfaceHF.hlsli"
#include "waterVolumetricsHF.hlsli"

// https://github.com/dong-zhan/ray-tracer/blob/master/ray%20cone%20intersect.hlsl
//p: ray o
//v: ray dir
//pa: cone apex
//va: cone dir
//sina2: sin(half apex angle) squared
//cosa2 : cos(half apex angle) squared
bool intersectInfiniteCone(float3 p, float3 v, float3 pa, float3 va, float sina2, float cosa2, out float tnear, out float tfar)
{
	tnear = 0;
	tfar = 0;

	float3 dp = p - pa;
	float vva = dot(v, va);
	float dpva = dot(dp, va);
	
	float3 v_vvava = v - vva * va;
	float3 dpvava = dp - dot(dp, va) * va;

	float A = cosa2 * dot(v_vvava, v_vvava) - sina2 * vva * vva;
	float B = 2 * cosa2 * dot(v_vvava, dpvava) - 2 * sina2 * vva * dpva;
	float C = cosa2 * dot(dpvava, dpvava) - sina2 * dpva * dpva;

	float disc = B*B - 4 * A*C;
	if (disc < 0.0001)
		return false;

	float sqrtDisc = sqrt(disc);
	tnear = (-B - sqrtDisc) / (2*A);
	tfar = (-B + sqrtDisc) / (2*A);
	return true;
}

float4 main(VertexToPixel input) : SV_TARGET
{
	ShaderEntity light = load_entity(spotlights().first_item() + (uint)g_xColor.x);

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

	if(g_xColor.w > 0)
	{
		// If camera is outside light volume, do a cone trace to determine closest point on cone:
		float tnear = 0;
		float tfar = 0;
		float2 sina2_cosa2 = unpack_half2(asuint(g_xColor.z));
		if (intersectInfiniteCone(GetCamera().position, -V, light.position, light.GetDirection(), sina2_cosa2.x, sina2_cosa2.y, tnear, tfar))
		{
			rayEnd -= V * max(0, tnear);
			//return float4(1,0,0,1);
		}
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
	for (uint i = 0; i < sampleCount; ++i)
	{
		float3 L = light.position - P;
		const half dist2 = dot(L, L);
		const half dist = sqrt(dist2);
		L /= dist;

		const half spot_factor = dot(L, light.GetDirection());
		const half spot_cutoff = light.GetConeAngleCos();

		[branch]
		if (spot_factor > spot_cutoff)
		{
			const half range = light.GetRange();
			const half range2 = range * range;
			half3 attenuation = attenuation_spotlight(dist2, range, light.GetRange2Rcp(), spot_factor, light.GetAngleScale(), light.GetAngleOffset());

			[branch]
			if (light.IsCastingShadow())
			{
				float4 shadow_pos = mul(load_entitymatrix(light.GetMatrixIndex() + 0), float4(P, 1));
				shadow_pos.xyz /= shadow_pos.w;
				float2 shadow_uv = shadow_pos.xy * float2(0.5f, -0.5f) + float2(0.5f, 0.5f);
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
				// Physical single scattering, attenuated on both legs. Safe
				// inside the cone gate because the transmittance is analytic -
				// there is no running state to advance on rejected samples.
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
