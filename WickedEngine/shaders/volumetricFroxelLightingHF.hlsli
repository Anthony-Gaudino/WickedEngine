#ifndef WI_VOLUMETRIC_FROXEL_LIGHTING_HF
#define WI_VOLUMETRIC_FROXEL_LIGHTING_HF
/**
 * What the lights scatter towards the eye at a point in the volume.
 *
 * Separate from `volumetricFroxelHF.hlsli` on purpose. That header is included
 * by every shader that draws anything and must stay small; this one pulls in
 * the whole lighting and shadow machinery and belongs only to the passes that
 * build the volume.
 *
 * Include after `globals.hlsli`, `lightingHF.hlsli` and `fogHF.hlsli`.
 */
#include "volumetricFroxelHF.hlsli"

/**
 * Reads one bucket of the transparent entity tile list.
 *
 * **The transparent list, never the opaque one.** The opaque list is clipped
 * against the depth buffer, so a lamp standing behind a wall is dropped from
 * it. That is right for shading the wall and wrong here: the air in *front* of
 * the wall is still lit by that lamp, and dropping it makes a shaft vanish the
 * moment the geometry behind it comes closer than the light.
 *
 * Addressed here rather than through `load_entitytile`, which picks its list
 * from whether `TRANSPARENT` happens to be defined. Defining that in a compute
 * shader to steer a macro would also silently change four unrelated decisions
 * in `lightingHF.hlsli`.
 *
 * @param[in] camera - Camera whose tile list is being read.
 * @param[in] tileIndex - Flattened tile index plus bucket.
 *
 * @return The bucket's 32 entity bits.
 */
inline uint VolumetricFroxelEntityTile(in ShaderCamera camera, in uint tileIndex)
{
	return bindless_structured_uint[
		descriptor_index(camera.buffer_entitytiles_index)][
			camera.entity_culling_tile_offset_transparent + tileIndex];
}

/**
 * Light scattered towards the eye at a point, summed over every light.
 *
 * **Per metre of path, and not yet multiplied by the medium's scattering
 * coefficient.** What comes back is the radiance a scattering event at this
 * point would send towards the eye; how much scattering happens there is the
 * medium's business and is applied by the caller. That split is what lets one
 * light evaluation serve air and water alike.
 *
 * The tile list is two-dimensional, so every light in a screen column is listed
 * for every cell of that column however far apart they are. Each light is
 * therefore range- and cone-tested against the cell's own position before
 * anything expensive happens: there is one cell per slice, and a shadow lookup
 * for a lamp a hundred metres away is pure waste repeated 128 times.
 *
 * References:
 * https://advances.realtimerendering.com/s2014/index.html
 * (Wronski, *Volumetric Fog*, SIGGRAPH 2014)
 *
 * @param[in] position - World position of the point.
 * @param[in] toEye - Normalized direction from the point towards the eye.
 * @param[in] screenUV - Screen position of the cell (0-1), for the tile lookup.
 * @param[in] pixel - Cell coordinate the shadow lookups dither against.
 *
 * @return Radiance scattered towards the eye, per metre of path.
 */
inline float3 VolumetricFroxelScatteredLight(
	in float3 position,
	in float3 toEye,
	in float2 screenUV,
	in min16uint2 pixel
)
{
	const ShaderCamera camera = GetCamera();

	float3 scattered = 0;

	// Directional lights
	//==========================================================================
	//
	// Never culled, so they are swept rather than looked up in a tile.
	{
		ShaderEntityIterator iterator = directional_lights();
		for (uint index = iterator.first_item();
			index <= iterator.last_item() && !iterator.empty();
			++index)
		{
			ShaderEntity light = load_entity(index);

			if (!light.IsVolumetricsEnabled())
			{
				continue;
			}

			const half3 L = light.GetDirection();

			half3 shadow = 1;
			for (uint cascade = 0;
				cascade < light.GetShadowCascadeCount();
				++cascade)
			{
				// Ortho matrix, so no divide by w.
				const float3 shadowPosition = mul(
					load_entitymatrix(light.GetMatrixIndex() + cascade),
					float4(position, 1)).xyz;
				const float3 shadowUV = clipspace_to_uv(shadowPosition);

				[branch]
				if (is_saturated(shadowUV))
				{
					shadow *= shadow_2D(
						light, shadowPosition.z, shadowUV.xy, cascade, pixel);
					break;
				}
			}

			if (GetFrame().options & OPTION_BIT_VOLUMETRICCLOUDS_CAST_SHADOW)
			{
				shadow *= shadow_2D_volumetricclouds(position);
			}

			half3 lightColor = light.GetColor().rgb;
			if (GetFrame().options & OPTION_BIT_REALISTIC_SKY)
			{
				lightColor *= GetAtmosphericLightTransmittance(
					GetWeather().atmosphere,
					position,
					L,
					texture_transmittancelut);
			}

			scattered += (float3)(lightColor * shadow)
				* ComputeScattering(saturate(dot(L, -toEye)))
				* (1 + GetFrame().GetVolumetricBoost(index));
		}
	}

	[branch]
	if (camera.buffer_entitytiles_index < 0)
	{
		return scattered;
	}

	const uint2 screenPixel = uint2(screenUV * camera.internal_resolution);
	const uint flatTileIndex = flatten2D(
		screenPixel / TILED_CULLING_BLOCKSIZE,
		camera.entity_culling_tilecount.xy) * SHADER_ENTITY_TILE_BUCKET_COUNT;

	// Spot lights
	//==========================================================================
	[branch]
	if (!spotlights().empty())
	{
		ShaderEntityIterator iterator = spotlights();
		for (uint bucket = iterator.first_bucket();
			bucket <= iterator.last_bucket();
			++bucket)
		{
			uint bucketBits = iterator.mask_entity(
				bucket, VolumetricFroxelEntityTile(camera, flatTileIndex + bucket));

			// Bucket scalarizer - Siggraph 2017 - Improved Culling
			// [Michal Drobot]
			bucketBits = WaveReadLaneFirst(WaveActiveBitOr(bucketBits));

			[loop]
			while (bucketBits != 0)
			{
				const uint bitIndex = firstbitlow(bucketBits);
				bucketBits ^= 1u << bitIndex;

				const uint index = bucket * 32 + bitIndex;
				ShaderEntity light = load_entity(index);

				if (!light.IsVolumetricsEnabled())
				{
					continue;
				}

				float3 L = light.position - position;
				const float distanceSq = dot(L, L);
				const float distance = sqrt(distanceSq);
				L /= distance;

				if (distance > light.GetRange())
				{
					continue;
				}

				const half spotFactor = dot((half3)L, light.GetDirection());
				if (spotFactor <= light.GetConeAngleCos())
				{
					continue;
				}

				float3 attenuation = attenuation_spotlight(
					(half)distanceSq,
					light.GetRange(),
					light.GetRange2Rcp(),
					spotFactor,
					light.GetAngleScale(),
					light.GetAngleOffset());

				[branch]
				if (light.IsCastingShadow())
				{
					float4 shadowPosition = mul(
						load_entitymatrix(light.GetMatrixIndex()),
						float4(position, 1));
					shadowPosition.xyz /= shadowPosition.w;
					const float2 shadowUV = clipspace_to_uv(shadowPosition.xy);

					[branch]
					if (is_saturated(shadowUV))
					{
						attenuation *= shadow_2D(
							light, shadowPosition.z, shadowUV, 0, pixel);
					}
				}

				const uint maskTexture = light.GetTextureIndex();
				[branch]
				if (maskTexture > 0)
				{
					float4 shadowPosition = mul(
						load_entitymatrix(light.GetMatrixIndex()),
						float4(position, 1));
					shadowPosition.xyz /= shadowPosition.w;
					const half4 mask =
						bindless_textures_half4[descriptor_index(maskTexture)]
							.SampleLevel(
								sampler_linear_clamp,
								clipspace_to_uv(shadowPosition.xy),
								0);
					attenuation *= mask.rgb * mask.a;
				}

				scattered += (float3)light.GetColor().rgb * attenuation
					* ComputeScattering(saturate(dot(L, -toEye)))
					* (1 + GetFrame().GetVolumetricBoost(index));
			}
		}
	}

	// Point lights
	//==========================================================================
	[branch]
	if (!pointlights().empty())
	{
		ShaderEntityIterator iterator = pointlights();
		for (uint bucket = iterator.first_bucket();
			bucket <= iterator.last_bucket();
			++bucket)
		{
			uint bucketBits = iterator.mask_entity(
				bucket, VolumetricFroxelEntityTile(camera, flatTileIndex + bucket));

			bucketBits = WaveReadLaneFirst(WaveActiveBitOr(bucketBits));

			[loop]
			while (bucketBits != 0)
			{
				const uint bitIndex = firstbitlow(bucketBits);
				bucketBits ^= 1u << bitIndex;

				const uint index = bucket * 32 + bitIndex;
				ShaderEntity light = load_entity(index);

				if (!light.IsVolumetricsEnabled())
				{
					continue;
				}

				const float3 toLight = light.position - position;
				const float distanceSq = dot(toLight, toLight);
				const float distance = sqrt(distanceSq);
				const float3 L = toLight / distance;

				if (distance > light.GetRange())
				{
					continue;
				}

				float3 attenuation = attenuation_pointlight(
					(half)distanceSq, light.GetRange(), light.GetRange2Rcp());

				[branch]
				if (light.IsCastingShadow())
				{
					attenuation *= shadow_cube(light, toLight, pixel);
				}

				const uint maskTexture = light.GetTextureIndex();
				[branch]
				if (maskTexture > 0)
				{
					const half4 mask =
						bindless_cubemaps_half4[descriptor_index(maskTexture)]
							.SampleLevel(sampler_linear_clamp, -toLight, 0);
					attenuation *= mask.rgb * mask.a;
				}

				scattered += (float3)light.GetColor().rgb * attenuation
					* ComputeScattering(saturate(dot(L, -toEye)))
					* (1 + GetFrame().GetVolumetricBoost(index));
			}
		}
	}

	// Rectangle lights
	//==========================================================================
	[branch]
	if (!rectlights().empty())
	{
		ShaderEntityIterator iterator = rectlights();
		for (uint bucket = iterator.first_bucket();
			bucket <= iterator.last_bucket();
			++bucket)
		{
			uint bucketBits = iterator.mask_entity(
				bucket, VolumetricFroxelEntityTile(camera, flatTileIndex + bucket));

			bucketBits = WaveReadLaneFirst(WaveActiveBitOr(bucketBits));

			[loop]
			while (bucketBits != 0)
			{
				const uint bitIndex = firstbitlow(bucketBits);
				bucketBits ^= 1u << bitIndex;

				const uint index = bucket * 32 + bitIndex;
				ShaderEntity light = load_entity(index);

				if (!light.IsVolumetricsEnabled())
				{
					continue;
				}

				const float3 toLight = light.position - position;
				const float distanceSq = dot(toLight, toLight);
				const float distance = sqrt(distanceSq);
				const float3 L = toLight / distance;

				if (distance > light.GetRange())
				{
					continue;
				}

				const half4 quaternion = light.GetQuaternion();
				const half3 right = rotate_vector(half3(1, 0, 0), quaternion);
				const half3 up = rotate_vector(half3(0, 1, 0), quaternion);
				const half3 forward = cross(up, right);

				// Only the lit face of the panel emits.
				if (dot(position - light.position, forward) <= 0)
				{
					continue;
				}

				const half halfLength = max(0.01, light.GetLength()) * 0.5;
				const half halfHeight = max(0.01, light.GetHeight()) * 0.5;

				const float3 p0 =
					light.position - right * halfLength + up * halfHeight;
				const float3 p1 =
					light.position + right * halfLength + up * halfHeight;
				const float3 p2 =
					light.position + right * halfLength - up * halfHeight;
				const float3 p3 =
					light.position - right * halfLength - up * halfHeight;

				// Solid angle of the panel as seen from this cell, after
				// Lagarde and de Rousiers, *Moving Frostbite to Physically
				// Based Rendering*, Siggraph 2014.
				const float3 v0 = normalize(p0 - position);
				const float3 v1 = normalize(p1 - position);
				const float3 v2 = normalize(p2 - position);
				const float3 v3 = normalize(p3 - position);
				const float3 n0 = normalize(cross(v0, v1));
				const float3 n1 = normalize(cross(v1, v2));
				const float3 n2 = normalize(cross(v2, v3));
				const float3 n3 = normalize(cross(v3, v0));
				const float solidAngle = saturate(
					acos(dot(-n0, n1)) +
					acos(dot(-n1, n2)) +
					acos(dot(-n2, n3)) +
					acos(dot(-n3, n0)) -
					(2 * PI));

				float3 attenuation = attenuation_pointlight(
					(half)distanceSq, light.GetRange(), light.GetRange2Rcp());

				attenuation *= solidAngle * 0.25 * (
					saturate(dot(v0, L)) +
					saturate(dot(v1, L)) +
					saturate(dot(v2, L)) +
					saturate(dot(v3, L)));

				[branch]
				if (light.IsCastingShadow())
				{
					float4 shadowPosition = mul(
						load_entitymatrix(light.GetMatrixIndex()),
						float4(position, 1));
					shadowPosition.xyz /= shadowPosition.w;
					const float2 shadowUV = clipspace_to_uv(shadowPosition.xy);

					[branch]
					if (is_saturated(shadowUV))
					{
						attenuation *= shadow_2D(
							light, shadowPosition.z, shadowUV, 0, pixel);
					}
				}

				scattered += (float3)light.GetColor().rgb * attenuation
					* ComputeScattering(saturate(dot(L, -toEye)))
					* (1 + GetFrame().GetVolumetricBoost(index));
			}
		}
	}

	return scattered;
}

#endif // WI_VOLUMETRIC_FROXEL_LIGHTING_HF
