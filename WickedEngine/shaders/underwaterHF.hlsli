#ifndef WI_UNDERWATER_HF
#define WI_UNDERWATER_HF
/**
 * Shared underwater helpers.
 *
 * Provides utilities to determine pixels that are submerged for effect masking.
 * Include after `globals.hlsli` (requires `GetWeather()`, `GetCamera()`,
 * `ShaderOcean`, bindless textures, and `reconstruct_position()`).
 */

/**
 * Computes the fraction of a pixel submerged below the ocean surface.
 *
 * Reconstructs the world position at the near plane and compares against the
 * displaced ocean surface height to determine if the pixel is underwater.
 * Used to mask screen-space effects (depth-of-field, chromatic aberration) to
 * only the submerged region, preserving the waterline boundary.
 *
 * @param[in] uv - Screen space UV coordinates (0-1).
 *
 * @return Submersion factor: 1 for fully submerged pixels, 0 for above water,
 *         smooth transition near the waterline.
 */
inline float ocean_underwater_factor(in float2 uv)
{
	const ShaderOcean ocean = GetWeather().ocean;
	// Reconstruct a world position just in front of the near plane for this
	// pixel; its height vs. the ocean surface decides submersion.
	const float cutoff = compute_inverse_lineardepth(max(GetCamera().z_near + 1, 1.0));
	float4 unproj = mul(GetCamera().inverse_view_projection, float4(uv_to_clipspace(uv), cutoff, 1));
	const float3 world_pos = unproj.xyz / unproj.w;
	float3 ocean_pos = float3(world_pos.x, ocean.water_height, world_pos.z);

	[branch]
	if (ocean.texture_displacementmap >= 0)
	{
		const float2 ocean_uv = ocean_pos.xz * ocean.patch_size_rcp;
		Texture2D texture_displacementmap = bindless_textures[descriptor_index(ocean.texture_displacementmap)];
		ocean_pos += texture_displacementmap.SampleLevel(sampler_linear_wrap, ocean_uv, 0).xzy;
	}

	return 1 - smoothstep(0.0, 0.025, saturate(world_pos.y - ocean_pos.y - 0.01));
}

#endif // WI_UNDERWATER_HF
