#ifndef WI_VOLUMETRIC_FROXEL_HF
#define WI_VOLUMETRIC_FROXEL_HF
/**
 * Shared geometry of the volumetric froxel volume.
 *
 * Everything that builds the volume and everything that reads it derives its
 * coordinates from here. Two places deriving the same mapping is how a lookup
 * drifts off the thing it is looking up, and the symptom - a smooth,
 * depth-dependent offset - reads as a tuning problem rather than as a bug.
 *
 * Include after `globals.hlsli`.
 */
#include "ShaderInterop_VolumetricFroxels.h"
#include "underwaterHF.hlsli"
#include "waterSegmentHF.hlsli"

/**
 * Radial distance from the eye to a slice boundary, in metres.
 *
 * Squared rather than linear, matching the aerial perspective volume this sits
 * beside. It spends resolution near the eye, which is where a medium is thick
 * enough to matter and where the light that crossed water still survives,
 * without collapsing the far field the way an exponential distribution does:
 * over 500 m in 128 slices the first cell is 3 cm deep and the last 8 m,
 * where an exponential of the same span would end in cells tens of metres
 * across.
 *
 * **Radial distance, not view depth**, because the volume is built along
 * normalized view rays. Anything sampling it has to measure the same way or the
 * error grows towards the edges of the screen, where the two diverge most.
 *
 * @param[in] slice - Slice coordinate, 0 at the eye and
 *                    `VOLUMETRIC_FROXEL_SLICES` at the far end.
 * @param[in] range - How far the volume reaches, in metres.
 *
 * @return Radial distance from the eye to that boundary (in metres).
 */
inline float VolumetricFroxelSliceToDepth(in float slice, in float range)
{
	const float normalized = slice / (float)VOLUMETRIC_FROXEL_SLICES;

	return range * normalized * normalized;
}

/**
 * Texture W coordinate for a distance from the eye.
 *
 * The inverse of `VolumetricFroxelSliceToDepth`, in the 0-1 a sampler wants.
 * Saturated here rather than left to the sampler's clamp, so the far behaviour
 * is stated next to the mapping it belongs to: everything past the volume's
 * range reads the last slice, and what lies beyond that is the tail's business
 * rather than this mapping's.
 *
 * @param[in] depth - Radial distance from the eye (in metres).
 * @param[in] range - How far the volume reaches, in metres.
 *
 * @return W coordinate in 0-1.
 */
inline float VolumetricFroxelDepthToW(in float depth, in float range)
{
	return sqrt(saturate(depth / max(range, 0.0001)));
}

/**
 * Direction of the view ray a column of the volume is built along.
 *
 * Interpolated between the camera's own frustum corners, which is how the rest
 * of the engine reconstructs a view ray, so the volume lines up with the
 * fragments that will sample it.
 *
 * **Not unprojected through `inverse_view_projection`.** Depth is reversed, so
 * the far plane sits at clip \f$z = 0\f$, and recovering a point there divides
 * by a \f$w\f$ that has gone to nearly nothing - the direction comes back
 * quantised. The corners are unprojected once on the CPU and interpolated,
 * which carries no such division.
 *
 * @param[in] uv - Screen space position of the cell (0-1).
 *
 * @return Normalized direction from the eye through that point.
 */
inline float3 VolumetricFroxelRayDirection(in float2 uv)
{
	return normalize(
		GetCamera().frustum_corners.screen_to_farplane(uv)
			- GetCamera().position);
}

/**
 * World position of a point inside the volume.
 *
 * @param[in] uv - Screen space position of the cell centre (0-1).
 * @param[in] depth - Radial distance from the eye (in metres).
 *
 * @return World position of that point.
 */
inline float3 VolumetricFroxelPosition(in float2 uv, in float depth)
{
	return GetCamera().position + VolumetricFroxelRayDirection(uv) * depth;
}

/**
 * Where inside a cell to take that cell's single sample, in 0-1 per axis.
 *
 * A cell stands for the average of the medium across its whole volume, but only
 * one sample is affordable, so where that sample sits decides what the cell
 * reports. Take them all at cell centres and every cell in a slice tests the
 * shadow map at the same distance at once, and the slice boundary becomes a
 * shell the eye can see: a shaft that should be a straight line is cut into
 * repeating bands, one per slice. Spreading the sample points inside their
 * cells trades that shell for a fine pattern which the volume's own bilinear
 * upsampling then averages away.
 *
 * An ordered Bayer pattern spatially, so that within any one frame the offsets
 * are spread evenly over each 8x8 block of cells rather than clumping the way
 * random ones would. The slice index is folded into the lookup on both axes with
 * different strides, so neighbouring slices draw different offsets instead of
 * stacking the same tile through the depth of the volume.
 *
 * **Advanced every frame by an irrational step**, which is what makes the
 * temporal accumulation worth having: a pattern that stood still would sample
 * the same 64 points forever and average to a fixed, biased answer. Stepping by
 * the golden ratio and wrapping spreads each cell's successive offsets evenly
 * over its own extent for any number of frames, with no period to line up
 * against - the Cranley-Patterson rotation.
 *
 * @param[in] cell - Integer cell coordinate within the volume.
 * @param[in] frame - Frames built since the last reset.
 *
 * @return Offset within the cell, in 0-1 on each axis.
 */
inline float3 VolumetricFroxelCellJitter(in uint3 cell, in uint frame)
{
	const uint2 xyOffset = uint2(cell.z * 5, cell.z * 3);
	const uint2 zOffset = uint2(cell.z * 3 + 2, cell.z * 7 + 5);

	const float3 spatial = float3(
		(float)dither((min16uint2)(cell.xy + xyOffset)),
		(float)dither((min16uint2)(cell.yx + xyOffset)),
		(float)dither((min16uint2)(cell.xy + zOffset)));

	// Three different irrationals, so the axes advance independently instead of
	// the sample point walking one straight line through the cell.
	const float3 rotation = frac(
		(float)frame * float3(0.6180340, 0.7548777, 0.8191725));

	return frac(spatial + rotation);
}

/**
 * Whether the volume carries the sun's in-scattering for this camera.
 *
 * The fog holds an analytic sun term of its own, and the volume supplies the
 * same light with the shadows on it, so exactly one of them may contribute.
 * Asked here rather than tested at each site, because the answer is a policy
 * about which pass owns the sun and two copies of it would drift apart.
 *
 * **Both halves are needed.** The camera settles whether a volume was built at
 * all - a reflection, an environment probe or a shadow cascade has none, and
 * keeps its analytic term. The frame flag settles whether the SUN is in that
 * volume, which is a different question: the volume is built for any light with
 * volumetrics, so a scene lit by a single lamp would otherwise lose its fog
 * brightening entirely.
 *
 * @return true when the fog should leave the sun to the volume.
 */
inline bool VolumetricFroxelCarriesTheSun()
{
	return GetCamera().texture_volumetricfroxels_index >= 0
		&& (GetFrame().options & OPTION_BIT_VOLUMETRIC_SUN) != 0;
}

/**
 * How far along a view ray the volume still describes what is there.
 *
 * The volume holds **one medium per column, chosen at the eye**
 * (`VolumetricFroxelMediumAt`), so its running total stops being a description
 * of this ray at the surface the ray crosses: past that point it has gone on
 * gathering air through water, or water through air. Reading it at the
 * fragment's own distance hands a fragment the whole of that, and because the
 * volume is added rather than composited, the water's own fog cannot take any
 * of it back out.
 *
 * Visible wherever a dry eye looks at something under the sea: the sea bed
 * arrives fogged as though the air reached it, and where nothing was drawn at
 * all the sky - which the water correctly extinguishes to nothing - comes back
 * carrying an entire atmosphere's haze.
 *
 * **Which side the ray starts on comes from the camera's own height; the fade
 * between the two answers comes from the pixel.** They are different questions
 * and only one of them the pixel can answer. Where a ray starts is a fact about
 * the eye, and the pixel's probe reads a point a hand's breadth ahead of the
 * near plane - so a ray pointing up out of a trough reports air while the eye
 * is under it. The fade is the opposite: it is what sweeps the waterline down
 * the screen instead of snapping it, and only the pixel knows where that line
 * falls.
 *
 * **Where the water begins is traced, not inferred from the two ends.** Whether
 * a fragment is itself submerged is a different question from whether the ray
 * reaching it crossed water, and only the second one bounds the air: a hillside
 * standing in the open, seen through a wave, is dry at both ends of a ray that
 * spends metres inside the sea.
 *
 * @param[in] screenUV - Screen space position of the fragment (0-1).
 * @param[in] fragmentPosition - World position of the fragment.
 * @param[in] distanceToEye - Radial distance from the eye to the fragment, in
 *                            metres.
 *
 * @return Distance the volume may be read out to, in metres. The whole of
 *         `distanceToEye` wherever the ray stays in air.
 *
 * @note Costs a water segment trace for a dry eye whose ray dips into the wave
 *       slab, which is the only case that can have water on it. A segment clear
 *       of the slab at both ends is rejected on two compares, and a submerged
 *       eye on one height sample.
 *
 * @note The same segment `ApplyWaterFog` traces for this fragment immediately
 *       beforehand. Sharing one trace between them would halve the cost, at
 *       the price of a distance threaded through every call site.
 */
inline float VolumetricFroxelColumnLength(
	in float2 screenUV, in float3 fragmentPosition, in float distanceToEye
)
{
	[branch]
	if (!GetCamera().IsWaterFog())
	{
		return distanceToEye;
	}

	const float3 eye = GetCamera().position;
	const float eyeHeight = eye.y - ocean_drawn_surface_height(eye);

	// **An eye under the surface reads the whole of its column, whichever way
	// the ray points.** The volume fills a column with one medium and describes
	// it from the eye onwards: a column it filled with water is water all the
	// way down, and one it filled with air is the air a ray leaving the surface
	// looks up through. Neither stops being described part way along.
	//
	// **Asked of the camera's own height, never of the pixel's probe.** That
	// probe samples a point a hand's breadth ahead of the near plane, so a ray
	// pointing upwards out of a trough reports air while the eye is plainly
	// under - and cutting on it discards the water that ray really crosses,
	// which is most of what it was looking through.
	[branch]
	if (eyeHeight < 0)
	{
		return distanceToEye;
	}

	// Above every crest the sea can raise at both ends, so the straight line
	// between them is above it too and no water can stand on this ray. Costs
	// two compares and spares the trace below to every fragment in a scene the
	// waves are not in front of, which is most of them.
	const float slabTop =
		GetWeather().ocean.water_height + ocean_max_displacement();

	[branch]
	if (eye.y >= slabTop && fragmentPosition.y >= slabTop)
	{
		return distanceToEye;
	}

	// **Where the air ends, found by following the ray.** The volume's column
	// is filled with air for a dry eye, so it describes this ray only as far as
	// the first water on it.
	//
	// The two ends cannot stand in for the trace, cheap as they are: a fragment
	// standing in air is still seen THROUGH water wherever a crest lies between
	// it and the eye, and both ends of that ray are dry. Since the volume is
	// ADDED after the water fog rather than composited into it, haze handed out
	// past the water is haze the water can no longer take back.
	const WaterSegment water = TraceWaterSegment(eye, fragmentPosition, 0);

	[branch]
	if (!water.crosses)
	{
		return distanceToEye;
	}

	// Faded on the pixel's own submersion, which is what sweeps the waterline
	// down the screen as the camera comes down to meet it.
	return lerp(water.entry, distanceToEye, ocean_underwater_factor(screenUV));
}

/**
 * Adds the light scattered between the eye and a fragment.
 *
 * This is what the volume exists for. The inscatter reaching a fragment depends
 * only on where that fragment is, so every kind of fragment - opaque, glass,
 * particle, sprite, splat, ocean surface - takes exactly its own share, and
 * there is no compositing order left to choose because nothing is composited.
 *
 * **The weight is a coverage fraction, not a radiance.** `ApplyFog` and
 * `ApplyWaterFog` take a `background` radiance and use `min(background, color)`
 * to find the share to hold out of a *multiplicative* operation. That trick does
 * not carry over to an additive term, which needs to know what *fraction* of the
 * pixel came from behind. Wherever a shader composites a refraction it has that
 * fraction to hand - `objectHF.hlsli` builds its background as
 * `surface.refraction.rgb * (1 - surface.F) * surface.refraction.a`, so the
 * weight is `(1 - surface.F) * surface.refraction.a`.
 *
 * **Alpha blended draws pass 0.** There the hardware multiplies by the source
 * alpha and the destination carries its own share already, so weighting here as
 * well would count the coverage twice.
 *
 * **Additive draws must not call this at all.** The destination already holds
 * the inscatter for that pixel, and adding it once per additive draw multiplies
 * the haze by however many of them land there - the same reasoning
 * `ApplyWaterFogAdditive` is built on.
 *
 * A camera with no volume - a reflection, an environment probe, a shadow
 * cascade - carries -1 and returns here, so no call site needs a condition of
 * its own.
 *
 * Example usage:
 * @code
 * ApplyWaterFog(dist, color);
 * ApplyVolumetricLight(ScreenCoord, dist, 0, color);
 * @endcode
 *
 * Takes the world position rather than a distance, matching `ApplyWaterFog`
 * beside it. The volume is indexed by **radial** distance from the eye, which
 * differs from view depth by more and more towards the edges of the screen; a
 * caller handed that job would eventually pass the wrong one, and a smooth
 * error that grows towards the corners reads as a tuning problem rather than as
 * a bug.
 *
 * @param[in] screenUV - Screen space position of the fragment (0-1).
 * @param[in] fragmentPosition - World position of the fragment.
 * @param[in] backgroundWeight - Fraction of this fragment's colour that came
 *                               from behind it and already carries its own
 *                               inscatter (0-1). Per channel, because a
 *                               refraction's Fresnel term is, and zero unless
 *                               the shader composited a refraction itself.
 * @param[in,out] color - Fragment colour, lit in place. Alpha is untouched.
 */
inline void ApplyVolumetricLight(
	in float2 screenUV,
	in float3 fragmentPosition,
	in half3 backgroundWeight,
	inout half4 color
)
{
	[branch]
	if (GetCamera().texture_volumetricfroxels_index < 0)
	{
		return;
	}

	const float distanceToEye = VolumetricFroxelColumnLength(
		screenUV,
		fragmentPosition,
		length(fragmentPosition - GetCamera().position));

	float slice = VolumetricFroxelDepthToW(
		distanceToEye, GetFrame().volumetricfroxel_range)
		* VOLUMETRIC_FROXEL_SLICES;

	// Nothing is stored nearer than the first texel's own depth, and the
	// sampler's clamp would hand a fragment closer than that the whole of the
	// first slice - light gathered over a stretch it does not stand behind. The
	// first texel is held and faded in instead, which is exact at zero and at
	// the texel, and the only thing it can be in between.
	half weight = 1;
	if (slice < 0.5)
	{
		weight = (half)saturate(slice * 2);
		slice = 0.5;
	}

	const half3 inscatter = (half3)texture_volumetricfroxels.SampleLevel(
		sampler_linear_clamp,
		float3(screenUV, slice / VOLUMETRIC_FROXEL_SLICES),
		0).rgb;

	half3 total = inscatter * weight;

	// The column beyond the volume, faded in over the distance it occupies
	// rather than delivered at the boundary. It is zero at the last texel's own
	// depth, so the two meet without a step, and complete at the far plane.
	[branch]
	if (GetFrame().texture_volumetricfroxeltail_index >= 0)
	{
		const float lastTexelDepth = VolumetricFroxelSliceToDepth(
			(float)VOLUMETRIC_FROXEL_SLICES - 0.5,
			GetFrame().volumetricfroxel_range);

		const float tailSpan = GetCamera().z_far - lastTexelDepth;
		const float intoTail = distanceToEye - lastTexelDepth;

		[branch]
		if (intoTail > 0 && tailSpan > 0)
		{
			const half4 tail = texture_volumetricfroxeltail.SampleLevel(
				sampler_linear_clamp, screenUV, 0);

			// Shaped by the medium's own extinction over that stretch, so the
			// far field arrives the way a medium delivers it - most of it
			// early, tailing off - instead of climbing in a straight line.
			// Falls back to that straight line where there is no medium to
			// shape it.
			const float extinction = (float)tail.a;
			const float ramp = extinction > 0.000001
				? (1 - exp(-extinction * intoTail))
					/ max(1 - exp(-extinction * tailSpan), 0.000001)
				: (intoTail / tailSpan);

			total += tail.rgb * (half)saturate(ramp);
		}
	}

	color.rgb += total * (1 - backgroundWeight);
}

/**
 * Adds the light scattered between the eye and a fragment that has no
 * already-lit background of its own.
 *
 * @param[in] screenUV - Screen space position of the fragment (0-1).
 * @param[in] fragmentPosition - World position of the fragment.
 * @param[in,out] color - Fragment colour, lit in place.
 */
inline void ApplyVolumetricLight(
	in float2 screenUV, in float3 fragmentPosition, inout half4 color
)
{
	ApplyVolumetricLight(screenUV, fragmentPosition, 0, color);
}

/**
 * Adds the light scattered between the eye and a fragment drawn with
 * premultiplied alpha.
 *
 * The blend adds the source unweighted, so the coverage has to be applied here
 * instead - otherwise a mote covering a tenth of a pixel paints a shaft across
 * the nine tenths it never touched. Same reasoning as
 * `ApplyWaterFogPremultiplied`.
 *
 * @param[in] screenUV - Screen space position of the fragment (0-1).
 * @param[in] fragmentPosition - World position of the fragment.
 * @param[in,out] color - Premultiplied fragment colour, lit in place.
 */
inline void ApplyVolumetricLightPremultiplied(
	in float2 screenUV, in float3 fragmentPosition, inout half4 color
)
{
	ApplyVolumetricLight(screenUV, fragmentPosition, 1 - color.a, color);
}

#endif // WI_VOLUMETRIC_FROXEL_HF
