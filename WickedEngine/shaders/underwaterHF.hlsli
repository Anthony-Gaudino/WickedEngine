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
 * Fraction of the wave displacement flattened away at a distance from the eye.
 *
 * The surface mesh cannot carry wave detail finer than its cells, so the
 * displacement is faded out over the band where the cells grow past the wave
 * patch. That band is derived from the mesh on the CPU
 * (`Ocean::GetDisplacementFadeBand`) and published in `ShaderOcean`, so it
 * follows the patch size and the mesh resolution instead of being a pair of
 * authored metres that quietly stops being true when either changes.
 *
 * **Everything that fades with the drawn waves calls this**, the surface vertex
 * shader included. Two places computing the same curve is exactly how the
 * surface and the tests against it drift apart, and a test that classifies
 * against waves that are not on screen paints wave-shaped noise across the
 * whole distance.
 *
 * @param[in] distanceToEye - Distance from the eye to the point on the still
 *                            water plane, in metres. Taken as a parameter
 *                            rather than derived because the eye is not the
 *                            same camera in every pass.
 *
 * @return 0 where the waves are drawn at full height, 1 where the surface is
 *         the still plane, blended in between.
 */
inline float ocean_displacement_fade(in float distanceToEye)
{
	const float2 band = GetWeather().ocean.displacement_fade;

	// A degenerate band would divide by zero inside smoothstep and take every
	// caller with it; the ocean can be disabled, which zeroes this.
	return smoothstep(band.x, max(band.y, band.x + 1), distanceToEye);
}

/**
 * How far the waves reach from the still water plane, in metres.
 *
 * The tallest crest and the deepest trough are taken together, as one magnitude
 * either side of the plane: a wave train is very nearly symmetric about it, and
 * one bound covering both directions is one reduction instead of two.
 *
 * **Measured, not derived.** The value is reduced over the displacement map by
 * the pass that writes it (`oceanUpdateDisplacementMapCS`), so it describes the
 * sea that is on screen this frame rather than bounding it from the parameters.
 * A bound guessed from the spectrum has to be generous to stay safe, and a
 * generous bound gives away exactly the work it was introduced to save.
 *
 * Everything that asks whether a ray can meet the water at all tests against
 * this: a straight segment's lowest point is at one of its ends, so a segment
 * with both ends further than this from the plane cannot cross the surface, and
 * that is exact rather than conservative.
 *
 * Example usage:
 * @code
 * const float reach = ocean_max_displacement();
 * if (min(eyeHeight, targetHeight) > reach) { return; } // cannot cross
 * @endcode
 *
 * @return Largest displacement from the still plane, in metres. Zero where
 *         there is no displacement map, which is exact - the surface is then
 *         the still plane everywhere.
 */
inline float ocean_max_displacement()
{
	const int buffer = GetWeather().ocean.buffer_max_displacement;

	[branch]
	if (buffer < 0)
	{
		return 0;
	}

	return asfloat(bindless_buffers[descriptor_index(buffer)].Load(0));
}

/**
 * Distance in front of the near plane where the drawn ocean surface hands over
 * to the analytic waterline, in metres.
 *
 * A metre rather than the near plane itself, so the plane the waterline is
 * reconstructed on spans enough of the world to cross the screen instead of
 * snapping over a few centimetres of camera travel.
 */
static const float OCEAN_WATERLINE_HANDOFF = 1.0;

/**
 * Device depth of the plane where the surface hands over to the waterline.
 *
 * **Two things read this and they must agree.** The ocean surface shader
 * discards fragments nearer than this plane, and the underwater effects
 * reconstruct their waterline exactly on it, so between them the nearest metre
 * is drawn analytically and everything beyond it is drawn as geometry. Split
 * the two apart and either the mesh puts a hard geometric edge over the
 * waterline band, or a metre of water goes missing with nothing drawing it.
 *
 * @return Reverse-Z device depth of the handoff plane.
 */
inline float ocean_waterline_handoff_depth()
{
	return compute_inverse_lineardepth(
		max(GetCamera().z_near + OCEAN_WATERLINE_HANDOFF, 1.0));
}

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
	float4 unproj = mul(GetCamera().inverse_view_projection,
		float4(uv_to_clipspace(uv), ocean_waterline_handoff_depth(), 1));
	return unproj.xyz / unproj.w;
}

/**
 * Height of the ocean surface over a world space position, as the WORLD has it.
 *
 * The displaced surface, not the still water plane, so a wave crest counts as
 * water and a trough does not.
 *
 * **Depends on nothing but the position**, which is what makes it the right
 * question for anything that shades: whether a tree trunk stands in water is a
 * property of the world, and an answer that moved with the camera would slide
 * the waterline up and down the trunk as you walked towards it. Cached
 * world-space lighting - surfels, DDGI probes - would carry that drift into
 * everything they later light.
 *
 * Use `ocean_drawn_surface_height` instead wherever the answer has to agree
 * with the surface actually on screen.
 *
 * @param[in] position - World position to judge the surface height over.
 *
 * @return World height of the surface there.
 */
/**
 * The highest and lowest the drawn surface reaches over a region of the patch.
 *
 * Reads the ocean's min/max height pyramid. Level 0 covers two displacement
 * texels across, and each level after it covers twice as much ground, up to a
 * single texel bounding the whole patch. The bound is exact in the sense that
 * matters: the surface anywhere inside the region lies between the two, so a
 * ray passing entirely above the first or below the second cannot cross it
 * there, however long it runs.
 *
 * Point sampled, and it must stay that way - a filtered read would blend two
 * neighbouring bounds into something that bounds neither.
 *
 * @param[in] position - World position; only `xz` is used.
 * @param[in] level - Mip level of the pyramid to ask.
 *
 * @return `x` the highest and `y` the lowest world height, or the still plane
 *         in both where the pyramid is unavailable.
 */
inline float2 ocean_height_bounds(in float3 position, in uint level)
{
	const ShaderOcean ocean = GetWeather().ocean;

	[branch]
	if (ocean.texture_heightHierarchy < 0)
	{
		return ocean.water_height.xx;
	}

	Texture2D<float2> hierarchy =
		bindless_textures_float2[
			descriptor_index(ocean.texture_heightHierarchy)];

	const float2 patchUV = position.xz * ocean.patch_size_rcp;

	const float2 bounds = hierarchy.SampleLevel(
		sampler_point_wrap, patchUV, level);

	return ocean.water_height + bounds;
}

/**
 * The surface's bounds over a whole stretch of the patch, not just a point.
 *
 * Two point taps, at the stretch's own ends, unioned. The caller asks for a
 * level whose texel is wider than the stretch, so the two texels sampled
 * between them hold nearly all of it.
 *
 * **Nearly, not provably.** A stretch running diagonally can clip the corner of
 * a texel neither end sits in, and that corner is unbounded here. The slack
 * that covers it is the base level's own: it is built from a four by four
 * window for a cell spanning two texels, so each cell already reports a good
 * deal more than it strictly owns. Four taps would close the gap outright and
 * cost twice as much on the hottest path in the water - and this is asked at
 * every node of every descent.
 *
 * A bilinear tap would be cheaper still and simply wrong: it averages
 * neighbouring bounds into a value that bounds neither.
 *
 * @param[in] pointNear - World position at one end of the stretch; `xz` used.
 * @param[in] pointFar - World position at the other end.
 * @param[in] level - Mip level of the pyramid to ask.
 *
 * @return `x` the highest and `y` the lowest world height across the stretch.
 */
inline float2 ocean_height_bounds_region(
	in float3 pointNear, in float3 pointFar, in uint level
)
{
	const ShaderOcean ocean = GetWeather().ocean;

	[branch]
	if (ocean.texture_heightHierarchy < 0)
	{
		return ocean.water_height.xx;
	}

	Texture2D<float2> hierarchy =
		bindless_textures_float2[
			descriptor_index(ocean.texture_heightHierarchy)];

	uint width;
	uint height;
	uint levels;
	hierarchy.GetDimensions(0, width, height, levels);

	const uint clamped = min(level, levels - 1);
	const float2 uvNear = pointNear.xz * ocean.patch_size_rcp;
	const float2 uvFar = pointFar.xz * ocean.patch_size_rcp;

	const float2 near = hierarchy.SampleLevel(
		sampler_point_wrap, uvNear, clamped);
	const float2 far = hierarchy.SampleLevel(
		sampler_point_wrap, uvFar, clamped);

	return ocean.water_height + float2(
		max(near.x, far.x), min(near.y, far.y));
}

/**
 * The coarsest pyramid level whose texel is at least as wide as a stretch.
 *
 * `ocean_height_bounds_region` bounds half a texel either side of its centre,
 * so a stretch is covered once the texel is as wide as the stretch itself.
 * Asking a finer level than that would bound only part of it.
 *
 * @param[in] extent - Length of the stretch across the patch, in metres.
 *
 * @return Level to ask for. 0 where the stretch is shorter than a base texel.
 */
inline uint ocean_height_bounds_level(in float extent)
{
	const ShaderOcean ocean = GetWeather().ocean;

	[branch]
	if (ocean.texture_heightHierarchy < 0)
	{
		return 0;
	}

	Texture2D<float2> hierarchy =
		bindless_textures_float2[
			descriptor_index(ocean.texture_heightHierarchy)];

	uint width;
	uint height;
	uint levels;
	hierarchy.GetDimensions(0, width, height, levels);

	// Metres one base texel covers: the pyramid spans exactly one patch.
	const float baseTexel = rcp(max(ocean.patch_size_rcp, 0.000001))
		/ (float)max(width, 1u);

	const float wanted = log2(max(extent / max(baseTexel, 0.000001), 1.0));

	return min((uint)ceil(wanted), levels - 1);
}

inline float ocean_surface_height(in float3 position)
{
	const ShaderOcean ocean = GetWeather().ocean;
	float height = ocean.water_height;

	[branch]
	if (ocean.texture_displacementmap >= 0)
	{
		const float2 ocean_uv = position.xz * ocean.patch_size_rcp;
		Texture2D texture_displacementmap = bindless_textures[descriptor_index(ocean.texture_displacementmap)];

		height += texture_displacementmap.SampleLevel(
			sampler_linear_wrap, ocean_uv, 0).z;
	}

	return height;
}

/**
 * Height of the ocean surface over a position, as it is DRAWN.
 *
 * The vertex shader flattens the displacement towards the horizon
 * (`ocean_displacement_fade`), so past that band the surface on screen *is* the
 * still plane. A test that has to agree with what was rendered - which side of
 * the water a fragment was drawn on, or how much water a view ray crossed -
 * has to flatten with it, or it classifies against waves that were never there.
 *
 * **View dependent, and legitimately so**: every caller of this is answering a
 * question about a particular view. Anything shading a point in the world wants
 * `ocean_surface_height` instead, whose answer does not move with the camera.
 *
 * Measured to the point on the still plane above or below the position, which
 * is what the vertex shader measures to.
 *
 * @param[in] position - World position to judge the surface height over.
 *
 * @return World height of the drawn surface there.
 */
inline float ocean_drawn_surface_height(in float3 position)
{
	const ShaderOcean ocean = GetWeather().ocean;

	const float fade = ocean_displacement_fade(distance(
		GetCamera().position,
		float3(position.x, ocean.water_height, position.z)));

	// Returned before sampling where the displacement is faded away entirely,
	// rather than sampled and multiplied out, because most of what asks this is
	// far off - the sky covers every pixel no geometry wrote.
	[branch]
	if (fade >= 1)
	{
		return ocean.water_height;
	}

	return lerp(ocean_surface_height(position), ocean.water_height, fade);
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
	const float surface = ocean_drawn_surface_height(world_pos);

	return 1 - smoothstep(0.0, 0.025, saturate(world_pos.y - surface - 0.01));
}

#endif // WI_UNDERWATER_HF
