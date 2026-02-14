#ifndef WI_SKY_HF
#define WI_SKY_HF
#include "globals.hlsli"
#include "skyAtmosphere.hlsli"
#include "fogHF.hlsli"

// Custom Atmosphere based on: https://www.shadertoy.com/view/Ml2cWG
// Cloud noise based on: https://www.shadertoy.com/view/4tdSWr

float3 AccurateAtmosphericScattering(float2 pixelPosition, float3 rayOrigin, float3 rayDirection, float3 sunDirection, float3 sunColor,
	bool enableSun, bool darkMode, bool stationary, bool highQuality, bool perPixelNoise, bool receiveShadow)
{
    AtmosphereParameters atmosphere = GetWeather().atmosphere;

	// We get compiler warnings: "floating point division by zero" when stationary is true, but it gets handled by GetCameraPlanetPos anyway
    float3 skyRelativePosition = stationary ? float3(0.00001, 0.00001, 0.00001) : rayOrigin;
	
    float3 worldPosition = GetCameraPlanetPos(atmosphere, skyRelativePosition);
    float3 worldDirection = rayDirection;

    float viewHeight = length(worldPosition);

    float3 luminance = 0;

	// Switch to high quality when above atmosphere layer, if high quality is not already enabled:
	if (viewHeight < atmosphere.topRadius && !highQuality)
    {
        float2 uv;
        float3 upVector = normalize(worldPosition);
        float viewZenithCosAngle = dot(worldDirection, upVector);

        float3 sideVector = normalize(cross(upVector, worldDirection)); // Assumes non parallel vectors
        float3 forwardVector = normalize(cross(sideVector, upVector)); // Aligns toward the sun light but perpendicular to up vector
        float2 lightOnPlane = float2(dot(sunDirection, forwardVector), dot(sunDirection, sideVector));
        lightOnPlane = normalize(lightOnPlane);
        float lightViewCosAngle = lightOnPlane.x;

        bool intersectGround = RaySphereIntersectNearest(worldPosition, worldDirection, float3(0, 0, 0), atmosphere.bottomRadius) >= 0.0f;

        SkyViewLutParamsToUv(atmosphere, intersectGround, viewZenithCosAngle, lightViewCosAngle, viewHeight, uv);

		luminance = texture_skyviewlut.SampleLevel(sampler_linear_clamp, uv, 0).rgb;						// Darken scattered sky luminance near the sun if the moon covers it (improves visual eclipse)
						if (GetSunEclipseStrength() > 0.0f)
						{
							float3 moonDir_local = GetMoonDirection();
							float moonLenSq_local = dot(moonDir_local, moonDir_local);
							if (moonLenSq_local > 1e-6f)
							{
								float3 moonDirNorm_local = moonDir_local * rsqrt(moonLenSq_local);
								float3 referenceUp_local = abs(sunDirection.y) > 0.95f ? float3(1,0,0) : float3(0,1,0);
								float3 sunRight_local = normalize(cross(referenceUp_local, sunDirection));
								float3 sunUp_local = normalize(cross(sunDirection, sunRight_local));
								float2 sunPlane_local = float2(dot(worldDirection, sunRight_local), dot(worldDirection, sunUp_local));
								float2 moonPlane_local = float2(dot(moonDirNorm_local, sunRight_local), dot(moonDirNorm_local, sunUp_local));
					// `GetMoonSize()` is angular diameter (rad) in Weather; use half-angle for spherical radius
					float moonRadius_local = sin(0.5f * GetMoonSize());
					float2 delta_local = sunPlane_local - moonPlane_local;
								float distSq_local = dot(delta_local, delta_local);
								float moonRadiusSq_local = moonRadius_local * moonRadius_local;
								if (moonRadiusSq_local > 1e-8f)
								{
									float coverage_local = saturate(1.0f - distSq_local / moonRadiusSq_local);
									coverage_local = pow(coverage_local, 1.5f);
									float silhouette_local = saturate(GetSunEclipseStrength() * coverage_local);
									luminance = lerp(luminance, float3(0,0,0), silhouette_local * 0.9f);
								}
							}
						}	}
    else
    {
        // Move to top atmosphere as the starting point for ray marching.
		// This is critical to be after the above to not disrupt above atmosphere tests and voxel selection.
		if (MoveToTopAtmosphere(worldPosition, worldDirection, atmosphere.topRadius))
		{
			// Apply the start offset after moving to the top of atmosphere to avoid black pixels
			worldPosition += worldDirection * AP_START_OFFSET_KM;

			float3 sunIlluminance = sunColor;

			const float tDepth = 0.0;
			const float sampleCountIni = 0.0;
			const bool variableSampleCount = true;
			const bool opaque = false;
			const bool ground = false;
			const bool mieRayPhase = true;
			const bool multiScatteringApprox = true;
			const bool volumetricCloudShadow = receiveShadow;
			const bool opaqueShadow = receiveShadow;
			SingleScatteringResult ss = IntegrateScatteredLuminance(
				atmosphere, pixelPosition, worldPosition, worldDirection, sunDirection, sunIlluminance, tDepth, sampleCountIni, variableSampleCount,
				perPixelNoise, opaque, ground, mieRayPhase, multiScatteringApprox, volumetricCloudShadow, opaqueShadow, texture_transmittancelut, texture_multiscatteringlut);

			luminance = ss.L;						// Also darken scattered luminance near the sun when the moon covers it
						if (GetSunEclipseStrength() > 0.0f)
						{
							float3 moonDir_local = GetMoonDirection();
							float moonLenSq_local = dot(moonDir_local, moonDir_local);
							if (moonLenSq_local > 1e-6f)
							{
								float3 moonDirNorm_local = moonDir_local * rsqrt(moonLenSq_local);
								float3 referenceUp_local = abs(sunDirection.y) > 0.95f ? float3(1,0,0) : float3(0,1,0);
								float3 sunRight_local = normalize(cross(referenceUp_local, sunDirection));
								float3 sunUp_local = normalize(cross(sunDirection, sunRight_local));
								float2 sunPlane_local = float2(dot(worldDirection, sunRight_local), dot(worldDirection, sunUp_local));
								float2 moonPlane_local = float2(dot(moonDirNorm_local, sunRight_local), dot(moonDirNorm_local, sunUp_local));
					// `GetMoonSize()` is angular diameter (rad) in Weather; use half-angle for spherical radius
					float moonRadius_local = sin(0.5f * GetMoonSize());
								float2 delta_local = sunPlane_local - moonPlane_local;
								float distSq_local = dot(delta_local, delta_local);
								float moonRadiusSq_local = moonRadius_local * moonRadius_local;
								if (moonRadiusSq_local > 1e-8f)
								{
									float coverage_local = saturate(1.0f - distSq_local / moonRadiusSq_local);
									coverage_local = pow(coverage_local, 1.5f);
									float silhouette_local = saturate(GetSunEclipseStrength() * coverage_local);
									luminance = lerp(luminance, float3(0,0,0), silhouette_local * 0.9f);
								}
							}
						}		}
    }

    float3 totalColor = 0;

    if (enableSun)
    {
        float3 sunIlluminance = sunColor;
		totalColor = luminance + GetSunLuminance(worldPosition, worldDirection, sunDirection, sunIlluminance, atmosphere, texture_transmittancelut);
	}
    else
    {
        totalColor = luminance; // We cant really seperate mie from luminance due to precomputation, todo?
    }

    if (darkMode)
    {
        totalColor = max(pow(saturate(dot(sunDirection, rayDirection)), 64) * sunColor, 0) * luminance * 1.0;
    }

	if (GetFrame().options & OPTION_BIT_HEIGHT_FOG)
	{
		const float3 planet_center = atmosphere.planetCenter * SKY_UNIT_TO_M;
		const float bottom_radius = atmosphere.bottomRadius * SKY_UNIT_TO_M;
		const float top_radius = bottom_radius + GetWeather().fog.height_end;
		float dist = RaySphereIntersectNearest(rayOrigin, rayDirection, planet_center, top_radius);

		if(dist >= 0)
		{
			// Offset origin with fog start value.
			// We can't do this with normal distance due to infinite distance.
			const float3 offsetO = rayOrigin + rayDirection * GetWeather().fog.start;
			float4 fog = GetFog(dist, offsetO, rayDirection);
			if (fog.a > 0) // this check avoids switching to fully fogged above fog level camera at zero density
			{
				if (length(rayOrigin - planet_center) > top_radius) // check if we are above fog height sphere
				{
					// hack: flip fog when camera is above
					fog.a = 1 - fog.a; // this only supports non-premultiplied fog
				}
				totalColor = lerp(totalColor, fog.rgb, fog.a); // non-premultiplied fog
			}
		}
	}

	return totalColor;
}

// Returns sky color modulated by the sun and clouds
//	pixel	: screen pixel position
//	V		: view direction
float3 GetDynamicSkyColor(in float2 pixel, in float3 V, bool sun_enabled = true, bool dark_enabled = false, bool stationary = false, bool highQuality = false, bool perPixelNoise = false, bool receiveShadow = false, bool clouds_enabled = true)
{
    float3 sky = 0;

    if (GetFrame().options & OPTION_BIT_REALISTIC_SKY)
    {
        sky = AccurateAtmosphericScattering
        (
			pixel,
            GetCamera().position,       // Ray origin
            V,                          // Ray direction
			GetSunDirection(),          // Position of the sun
			GetSunColor(),              // Sun Color
            sun_enabled,                // Use sun and total
            dark_enabled,               // Enable dark mode for light shafts etc.
            stationary,					// Fixed position for ambient and environment capture.
			highQuality,				// Skip color lookup from lut
			perPixelNoise,				// Vary sampling position with TAA
			receiveShadow				// Atmosphere to use pre-rendered shadow data
        );
    }
    else
    {
		sky = lerp(GetHorizonColor(), GetZenithColor(), saturate(V.y * 0.5f + 0.5f));
    }

	if (!dark_enabled)
	{
		float moonSize = GetMoonSize();
		float moonHalfSize = 0.5f * moonSize; // shader expects half-angle for trig/radius
		float3 moonColor = GetMoonColor();
		if (moonSize > 0 && dot(moonColor, moonColor) > 0)
		{
			float3 moonDir = GetMoonDirection();
			float3 sunDir = (float3)GetSunDirection();
			float sunLenSq = dot(sunDir, sunDir);
			float3 sunToMoonDir = moonDir;
			float phaseVisibility = 1.0f;
			if (sunLenSq > 1e-6f)
			{
				float invSunLen = rsqrt(sunLenSq);
				sunToMoonDir = -sunDir * invSunLen;
				phaseVisibility = saturate(0.5f * (1.0f + dot(sunToMoonDir, moonDir)));
			}
			float eclipseStrength = GetMoonEclipseStrength();
			phaseVisibility *= saturate(1.0f - eclipseStrength);
			float cosAngle = dot(V, moonDir);
			float innerEdge = cos(moonHalfSize);
			float core = smoothstep(innerEdge, cos(moonHalfSize * 0.8f), cosAngle);
			float3 diskColor = moonColor;
			// Always mask the moon disk by geometry (core) so it occludes correctly even when unlit.
			// Lighting/texture only affect the lit hemisphere; unlit hemisphere will render as silhouette.
			float diskMask = 0.0f;
			if (phaseVisibility > 0.0f)
			{
				float3 referenceUp = abs(moonDir.y) > 0.95f ? float3(1, 0, 0) : float3(0, 1, 0);
				float3 moonRight = normalize(cross(referenceUp, moonDir));
				float3 moonUp = normalize(cross(moonDir, moonRight));
				float2 local = float2(dot(V, moonRight), dot(V, moonUp));
				// map angular offset -> [0,1] UV : use sin(half-angle) as spherical radius
				float invRadius = 0.5f / max(sin(moonHalfSize), 0.0001f);
				float2 moonUV = local * invRadius + 0.5f;
				if (all(moonUV >= 0.0f) && all(moonUV <= 1.0f))
				{
					float2 diskCoord = (moonUV - 0.5f) * 2.0f;
					float radialSq = dot(diskCoord, diskCoord);
					float localZ = sqrt(saturate(1.0f - radialSq));
					float3 localNormal = normalize(moonRight * diskCoord.x + moonUp * diskCoord.y + moonDir * localZ);
					float lit = saturate(dot(localNormal, sunToMoonDir));

					if (lit > 0.0f)
					{
						// If this moon pixel is also over the sun disk and eclipse is active,
						// render the moon fully opaque black to form a clean silhouette.
						// Use the physical apparent sun size here so the moon can occlude the actual bright solar disk.
						float sunApexAngleDegree_local = 0.545; // physical sun angular diameter (deg) used for eclipse silhouette
						float sunHalfApexAngleRadian_local = 0.5 * sunApexAngleDegree_local * PI / 180.0;
						float sunCosHalfApex_local = cos(sunHalfApexAngleRadian_local);
						float3 sunDir_local = (float3)GetSunDirection();
						float sunDotPixel_local = dot(V, sunDir_local);
							// compute sun disk mask using angular distance (cos-space offset was too wide for small sun)
							float sunAngle_local = acos(saturate(sunDotPixel_local));
							const float sunEdgeSoftRad_local = 0.0025; // radians (small soft edge)
							float sunDiskMask_local = 1.0 - smoothstep(sunHalfApexAngleRadian_local - sunEdgeSoftRad_local, sunHalfApexAngleRadian_local + sunEdgeSoftRad_local, sunAngle_local);
						float overlapMask = core * sunDiskMask_local * saturate(GetSunEclipseStrength());
									if (overlapMask > 1e-4)
						{
								diskMask = overlapMask; // only overlapping portion becomes silhouette
							diskColor = float3(0,0,0); // pure black silhouette
							// avoid lit-shading on this pixel
							lit = 0.0f;
						}
						else
						{
							diskMask = core * lit;
							if (HasMoonTexture())
							{
								float4 tex = bindless_textures[NonUniformResourceIndex(descriptor_index(GetWeather().moon_texture))].SampleLevel(sampler_linear_clamp, moonUV, GetMoonTextureMipBias());
								diskMask *= tex.a;
								diskColor *= tex.rgb;
							}
							float phaseBlend = saturate(phaseVisibility);
							float3 shadingNormal = normalize(lerp(moonDir, localNormal, phaseBlend));
							float diffuse = saturate(dot(shadingNormal, sunToMoonDir));
							float rim = pow(saturate(1.0f - dot(localNormal, V)), 4.0f);
							float wave = dot(diskCoord, float2(21.7f, -14.9f));
							float time = GetTime();
							float anim0 = sin(wave + time * 0.35f);
							float anim1 = sin((diskCoord.x * -18.3f + diskCoord.y * 11.1f) - time * 0.22f);
							float craterMask = 0.85f + 0.15f * (anim0 * 0.6f + anim1 * 0.4f);
							float detailMask = lerp(1.0f, craterMask, phaseBlend);
							float shading = lerp(0.35f, 1.0f, diffuse);
							diskColor *= shading * detailMask;
							diskColor += moonColor * rim * 0.08f * phaseBlend;
						}
					}
				}
			}

			float haloContribution = 0;
			float haloIntensity = GetMoonHaloIntensity();
			if (haloIntensity > 0)
			{
				float haloSize = max(GetMoonHaloSize(), 0.0001f);
			float haloRadius = moonHalfSize + haloSize; // haloSize is an additional radius (rad)
				float halo = smoothstep(cos(haloRadius), innerEdge, cosAngle);
				halo = pow(saturate(halo), max(GetMoonHaloSharpness(), 0.0001f));

				// Mask halo to the lit hemisphere so crescents don't get a full circular glow.
				// Use the view direction V and sunToMoonDir (direction toward lit side).
				float litHaloMask = saturate(dot(V, sunToMoonDir));
							haloContribution = halo * haloIntensity * litHaloMask; // show halo even for unlit/new moon (so UI halo controls are effective)
			}
			sky += moonColor * haloContribution;
		sky = lerp(sky, diskColor, diskMask * saturate(1.0f - eclipseStrength * 0.9f));
		}
	}

	sky *= GetWeather().sky_exposure;
	
	if (clouds_enabled && V.y > 0 && GetScene().texture_cloudmap >= 0 && !stationary)
	{
		float4 cloudmap = bindless_textures[descriptor_index(GetScene().texture_cloudmap)].SampleLevel(sampler_linear_clamp, encode_hemioct(V.xzy) * 0.5 + 0.5, 0);
		if (dark_enabled)
			cloudmap.rgb = 0;
		sky.rgb = sky.rgb * (1.0 - cloudmap.a) + cloudmap.rgb;
	}

	return sky;
}
float3 GetDynamicSkyColor(in float3 V, bool sun_enabled = true, bool dark_enabled = false, bool stationary = false, bool highQuality = false, bool perPixelNoise = false, bool receiveShadow = false, bool clouds_enabled = true)
{
	return GetDynamicSkyColor(float2(0.0f, 0.0f), V, sun_enabled, dark_enabled, stationary, highQuality, perPixelNoise, receiveShadow, clouds_enabled);
}

float3 GetStaticSkyColor(in float3 V, bool clouds_enabled = true)
{
	ShaderWeather weather = GetWeather();
	float2x2 rot = float2x2(
		weather.sky_rotation_cos, -weather.sky_rotation_sin,
		weather.sky_rotation_sin, weather.sky_rotation_cos
	);
	V.xz = mul(V.xz, rot);

	float3 sky = 0;
	if (GetFrame().options & OPTION_BIT_STATIC_SKY_SPHEREMAP)
	{
		float2 uv = (float2(-atan2(V.z, V.x) / PI, -V.y) + 1.0) * 0.5;
		sky = bindless_textures[descriptor_index(GetScene().globalenvmap)].SampleLevel(sampler_linear_clamp, uv, 0).rgb;
	}
	else
	{
		sky = bindless_cubemaps[descriptor_index(GetScene().globalenvmap)].SampleLevel(sampler_linear_clamp, V, 0).rgb;
	}
	
	sky *= GetWeather().sky_exposure;

	if (clouds_enabled && V.y > 0 && GetScene().texture_cloudmap >= 0)
	{
		float4 cloudmap = bindless_textures[descriptor_index(GetScene().texture_cloudmap)].SampleLevel(sampler_linear_clamp, encode_hemioct(V.xzy) * 0.5 + 0.5, 0);
		sky.rgb = sky.rgb * (1.0 - cloudmap.a) + cloudmap.rgb;
	}

	return sky;
}


#endif // WI_SKY_HF
