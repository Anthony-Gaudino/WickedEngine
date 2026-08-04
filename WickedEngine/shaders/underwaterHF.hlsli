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
 * Reconstructs the world position this pixel's submersion is judged at.
 *
 * Every underwater effect has to agree on where the waterline falls, or one
 * engages while another has not. Sharing this is what keeps them together.
 *
 * @param[in] uv - Screen space UV coordinates (0-1).
 *
 * @return World position on the test plane for this pixel.
 */
inline float3 ocean_underwater_test_position(in float2 uv)
{
	// A metre in front of the near plane rather than on it, so the plane this
	// reconstructs spans enough of the world for the waterline to cross the
	// screen instead of snapping over a few centimetres of camera travel.
	const float cutoff = compute_inverse_lineardepth(max(GetCamera().z_near + 1, 1.0));
	float4 unproj = mul(GetCamera().inverse_view_projection, float4(uv_to_clipspace(uv), cutoff, 1));
	return unproj.xyz / unproj.w;
}

/**
 * Height of the ocean surface above a world space xz position.
 *
 * The displaced surface, not the still water plane, so a wave crest counts as
 * water and a trough does not.
 *
 * @param[in] worldXZ - Horizontal world position to sample under.
 *
 * @return World height of the surface there.
 */
inline float ocean_surface_height(in float2 worldXZ)
{
	const ShaderOcean ocean = GetWeather().ocean;
	float height = ocean.water_height;

	[branch]
	if (ocean.texture_displacementmap >= 0)
	{
		const float2 ocean_uv = worldXZ * ocean.patch_size_rcp;
		Texture2D texture_displacementmap = bindless_textures[descriptor_index(ocean.texture_displacementmap)];
		height += texture_displacementmap.SampleLevel(sampler_linear_wrap, ocean_uv, 0).z;
	}

	return height;
}

/**
 * Discards a fragment lying on the far side of the ocean surface.
 *
 * The ocean is the only transparent that writes depth, so the transparent pass
 * draws the two sides of the water around it: the far side first, where the
 * surface's screen space refraction can pick it up, and the near side after,
 * where the surface's depth cannot reject it. A billboard straddling the
 * waterline belongs to **both**, so it is issued in both passes and each keeps
 * only its own half. That is why the test is per fragment rather than a
 * per-object classification - classifying by origin makes a sprite snap between
 * fully wet and fully dry as it crosses, and it can never be partly submerged.
 *
 * The two sides are strictly complementary, so a fragment exactly on the
 * surface is kept by one pass and dropped by the other, never by both (which
 * would blend it twice along the waterline) and never by neither (which would
 * leave a transparent seam).
 *
 * Example usage:
 * @code
 * ClipToWaterSide(P, flags & FLAG_SUBMERGED, flags & FLAG_ABOVE);
 * @endcode
 *
 * @param[in] position - World space position of the fragment.
 * @param[in] keepSubmerged - Keep only what is at or below the surface.
 * @param[in] keepAbove - Keep only what is strictly above the surface.
 *
 * @note Does nothing unless exactly one side is selected, so a draw that is not
 *       split pays only the flag test.
 */
inline void ClipToWaterSide(
	in float3 position,
	in bool keepSubmerged,
	in bool keepAbove
)
{
	[branch]
	if (keepSubmerged != keepAbove)
	{
		const float side = position.y - ocean_surface_height(position.xz);

		if (keepSubmerged ? (side > 0) : (side <= 0))
		{
			discard;
		}
	}
}

/**
 * Computes the fraction of a pixel submerged below the ocean surface.
 *
 * Used to mask screen-space effects - depth of field, chromatic aberration,
 * the underwater composite and the volumetric light marches - to the submerged
 * region, preserving the waterline boundary. They all share this one test so
 * that a camera crossing the surface brings them in together.
 *
 * @param[in] uv - Screen space UV coordinates (0-1).
 *
 * @return Submersion factor: 1 for fully submerged pixels, 0 for above water,
 *         smooth transition near the waterline.
 */
inline float ocean_underwater_factor(in float2 uv)
{
	const float3 world_pos = ocean_underwater_test_position(uv);
	const float surface = ocean_surface_height(world_pos.xz);

	return 1 - smoothstep(0.0, 0.025, saturate(world_pos.y - surface - 0.01));
}

#endif // WI_UNDERWATER_HF
