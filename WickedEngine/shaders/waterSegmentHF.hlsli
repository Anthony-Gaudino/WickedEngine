#ifndef WI_WATER_SEGMENT_HF
#define WI_WATER_SEGMENT_HF
/**
 * Where a straight segment meets the ocean surface, and how much water it
 * crossed.
 *
 * One question, asked along the ray, in place of the several that used to be
 * asked about a point or about the vertical. A point test - "is this fragment
 * above or below the surface" - and a vertical solve - "how far up is the
 * surface from here" - are both right on a flat sea and both wrong on a crest,
 * because a crest is neither at an end of the segment nor overhead. Everything
 * that has to know about water in front of something asks this instead.
 *
 * **Measured against the surface as DRAWN.** `ocean_drawn_surface_height`
 * flattens the waves towards the horizon exactly as the surface vertex shader
 * does, so the answer agrees with what was rasterized. Measuring against
 * `ocean_surface_height` instead would classify against waves that are not on
 * screen and paint their shape across the distance.
 *
 * Include after `globals.hlsli` and `underwaterHF.hlsli`.
 */

#include "globals.hlsli"
#include "underwaterHF.hlsli"

/**
 * How many points along the segment are asked which side of the surface they
 * are on.
 *
 * What this count buys is the thinnest crest that can still be found: a wave
 * narrower than one step can fall between two samples and be missed. It is not
 * the precision of the crossing, which comes from the refinement below.
 */
static const uint WATER_SEGMENT_STEPS = 12;

/**
 * How many times a bracketed crossing is halved.
 *
 * Each one doubles the precision of the crossing's position, so three take a
 * step of a metre down to about twelve centimetres. Cheap next to the search:
 * the bracket is already known, so these only ever refine it.
 */
static const uint WATER_SEGMENT_REFINEMENTS = 3;

/**
 * Sample phase for a caller whose answer is a decision rather than a quantity.
 *
 * The step centre, fixed. A moving phase is what keeps a continuous result from
 * terracing, but a fragment that is classified rather than measured takes one
 * branch or the other, and the two look nothing alike - a moving phase would
 * flip it between them from frame to frame, which arrives as flicker along a
 * crest's edge rather than as a dither the temporal pass can resolve.
 */
static const float WATER_SEGMENT_STABLE_PHASE = 0.5;

/**
 * Where a view segment meets the ocean, and how much water it crossed.
 */
struct WaterSegment
{
	/**
	 * Whether the segment meets the surface at all.
	 *
	 * False both for a segment entirely in the air and for one entirely under
	 * water - it says a *surface* was crossed, not that water was present. Use
	 * `submerged` for the latter.
	 */
	bool crosses;

	/** Distance from the origin to the first crossing (in metres). */
	float entry;

	/** Length of the segment lying in water (in metres). */
	float submerged;

	/** World position of the first crossing. */
	float3 entryPoint;

	/**
	 * Wave surface normal at the first crossing, pointing up out of the water.
	 */
	float3 entryNormal;
};

/**
 * Height of a point above the drawn ocean surface, signed.
 *
 * @param[in] position - World position to measure.
 *
 * @return Metres above the surface; negative under water, zero on it.
 */
inline float WaterSegmentHeightAbove(in float3 position)
{
	return position.y - ocean_drawn_surface_height(position);
}

/**
 * Surface normal of the waves over a world position.
 *
 * Rebuilt from the same gradient map the ocean surface is shaded with, so a
 * test against this normal ondulates with the waves that were drawn rather
 * than with a flat plane.
 *
 * @param[in] position - World position to take the normal over.
 *
 * @return Normalized surface normal, pointing up out of the water.
 */
inline float3 WaterSegmentSurfaceNormal(in float3 position)
{
	const ShaderOcean ocean = GetWeather().ocean;

	[branch]
	if (ocean.texture_gradientmap < 0)
	{
		return float3(0, 1, 0);
	}

	Texture2D texture_gradientmap =
		bindless_textures[descriptor_index(ocean.texture_gradientmap)];

	const float2 gradient = texture_gradientmap.SampleLevel(
		sampler_linear_wrap, position.xz * ocean.patch_size_rcp, 0).rg;

	return normalize(
		float3(gradient.x, ocean.texel_length * 2.0, gradient.y));
}

/**
 * Traces a straight segment against the drawn ocean surface.
 *
 * The surface is a height field, so the trace is a search for sign changes in
 * the height above it. Three things bound the work, and each of them is exact
 * rather than a tuning choice:
 *
 * - **The slab.** No crossing exists further from the still plane than the
 *   tallest wave (`ocean_max_displacement`), so the segment is clipped to that
 *   band before anything is sampled. A steep ray leaves it almost at once.
 * - **The fade.** Past `ocean_displacement_fade` the drawn surface *is* the
 *   still plane, so beyond that reach the crossing is a plane intersection and
 *   needs no samples. This is what bounds a near-horizontal ray, which can lie
 *   inside the slab for kilometres.
 * - **The ends.** A straight segment's extremes are its endpoints, so both ends
 *   clear of the slab on the same side settles the answer with no search at
 *   all.
 *
 * On a flat sea every sample of the height field returns the still plane, so
 * the search reduces to the plane intersection and the result is the same one
 * the plane test gives - exactly, not nearly.
 *
 * The first sample sits wherever the caller's `jitter` puts it, because a step
 * function sampled N times returns N+1 distinct answers and a fixed phase
 * prints them as terraces along a crest. A caller measuring a quantity passes
 * `blue_noise`, which carries its own per-frame phase for the temporal pass to
 * average out; a caller taking a decision passes `WATER_SEGMENT_STABLE_PHASE`
 * instead, for the reason recorded there.
 *
 * **At most two crossings are resolved**, the first and the last. A ray
 * crossing three or more times - two separate crests seen edge on - counts the
 * dry gap between them as water. The error is bounded by the gap and always
 * over-counts, which fogs slightly too much rather than showing water that is
 * not there.
 *
 * Example usage:
 * @code
 * const WaterSegment water = TraceWaterSegment(eye, P, blue_noise(pixel).x);
 * if (water.crosses) { ... } // the surface stands between the eye and P
 * const float3 transmittance = exp(-water.submerged * medium.sigmaT);
 * @endcode
 *
 * @param[in] origin - World position the segment starts at, normally the eye.
 *                     The wave flattening is measured from the camera, so an
 *                     origin far from it is answered against the surface as the
 *                     CAMERA sees it drawn.
 * @param[in] target - World position the segment ends at.
 * @param[in] jitter - Where inside its step the first sample sits, in [0, 1).
 *
 * @return What the segment met. All fields are zero where there is no ocean.
 */
inline WaterSegment TraceWaterSegment(
	in float3 origin,
	in float3 target,
	in float jitter
)
{
	WaterSegment result;
	result.crosses = false;
	result.entry = 0;
	result.submerged = 0;
	result.entryPoint = origin;
	result.entryNormal = float3(0, 1, 0);

	const ShaderOcean ocean = GetWeather().ocean;

	[branch]
	if (!ocean.IsValid())
	{
		return result;
	}

	const float3 delta = target - origin;
	const float segment = length(delta);

	[branch]
	if (segment < 0.0001)
	{
		return result;
	}

	const float3 direction = delta / segment;

	// Signed heights against the STILL plane, which is what the slab below is
	// measured in. The drawn surface is asked for separately, only where it can
	// differ from the plane.
	const float originPlane = origin.y - ocean.water_height;
	const float targetPlane = target.y - ocean.water_height;

	// Where the drawn surface can stand, either side of the plane.
	const float reach = ocean_max_displacement();

	// Both ends clear of the slab on the same side: the segment is straight, so
	// its extremes are its ends and nothing between them can reach the surface.
	[branch]
	if (min(originPlane, targetPlane) >= reach)
	{
		return result;
	}

	[branch]
	if (max(originPlane, targetPlane) <= -reach)
	{
		result.submerged = segment;
		return result;
	}

	// The stretch of the segment lying within the slab, which is the only place
	// a crossing can be. Solved rather than searched: the segment's height is
	// linear in the distance along it.
	float slabNear = 0;
	float slabFar = segment;

	[branch]
	if (abs(targetPlane - originPlane) > 0.0001)
	{
		const float toTop = (reach - originPlane) / (targetPlane - originPlane);
		const float toBottom =
			(-reach - originPlane) / (targetPlane - originPlane);

		slabNear = clamp(min(toTop, toBottom), 0, 1) * segment;
		slabFar = clamp(max(toTop, toBottom), 0, 1) * segment;
	}

	// A near-horizontal segment can lie inside the slab for kilometres, and
	// beyond the fade band there is nothing out there to find: the drawn
	// surface is the still plane, and a plane crossing is already covered by
	// the linear solve above. Walking further only spreads the samples thinner
	// over the near water, which is where the waves that matter are.
	const float marchFar = min(
		slabFar,
		slabNear + max(GetWeather().ocean.displacement_fade.y, 1.0));

	const float marched = marchFar - slabNear;

	[branch]
	if (marched <= 0)
	{
		return result;
	}

	// Which side each end is on, against the surface as drawn. The origin's
	// side says whether the spans between crossings hold water or air, and the
	// two together say whether the number of crossings is odd or even.
	const bool originSubmerged = WaterSegmentHeightAbove(origin) <= 0;
	const bool targetSubmerged = WaterSegmentHeightAbove(target) <= 0;

	// Spread so that the last step lands exactly on the far end: the stretch
	// has to be covered whole, or a crossing in the final fraction of it is
	// missed however many samples are spent on the rest.
	const float lastStep = (float)(WATER_SEGMENT_STEPS - 1);
	const float stepSize = marched / lastStep;

	// Brackets around the first and the last sign change found. Held as a
	// distance pair so the refinement below has something to halve.
	float firstLo = -1;
	float firstHi = -1;
	float lastLo = -1;
	float lastHi = -1;

	float previousDistance = slabNear;
	bool previousSubmerged = WaterSegmentHeightAbove(
		mad(previousDistance, direction, origin)) <= 0;

	for (uint step = 0; step < WATER_SEGMENT_STEPS; ++step)
	{
		const float distanceAlong =
			slabNear + min((float)step + jitter, lastStep) * stepSize;

		const bool sampleSubmerged = WaterSegmentHeightAbove(
			mad(distanceAlong, direction, origin)) <= 0;

		if (sampleSubmerged != previousSubmerged)
		{
			if (firstLo < 0)
			{
				firstLo = previousDistance;
				firstHi = distanceAlong;
			}

			lastLo = previousDistance;
			lastHi = distanceAlong;
		}

		previousDistance = distanceAlong;
		previousSubmerged = sampleSubmerged;
	}

	// Nothing was crossed inside the marched stretch. Beyond it the surface is
	// the still plane, so the ends decide: they disagree only if the plane was
	// crossed out there, and where they agree the whole segment is on one side.
	[branch]
	if (firstLo < 0)
	{
		[branch]
		if (originSubmerged == targetSubmerged)
		{
			result.submerged = originSubmerged ? segment : 0;
			return result;
		}

		// Solved against the plane, which is what the surface is out where this
		// crossing has to be - the marched stretch covers everywhere it is not.
		//
		// The difference keeps its sign. Flooring it would turn a segment
		// running upwards out of the water into one running down, and put the
		// crossing at the eye instead of at the horizon.
		const float heightDrop = originPlane - targetPlane;
		const float crossing = abs(heightDrop) > 0.0001
			? segment * saturate(originPlane / heightDrop)
			: 0;

		result.crosses = true;
		result.entry = crossing;
		result.submerged = originSubmerged ? crossing : segment - crossing;
		result.entryPoint = mad(crossing, direction, origin);
		result.entryNormal = WaterSegmentSurfaceNormal(result.entryPoint);

		return result;
	}

	// Halve each bracket down onto its crossing. The sample at the low end is
	// known to be on the origin's side of that crossing, so only the midpoint
	// has to be asked.
	float firstLow = firstLo;
	float firstHigh = firstHi;
	const bool lowSideSubmerged = WaterSegmentHeightAbove(
		mad(firstLow, direction, origin)) <= 0;

	for (uint refine = 0; refine < WATER_SEGMENT_REFINEMENTS; ++refine)
	{
		const float middle = 0.5 * (firstLow + firstHigh);
		const bool middleSubmerged = WaterSegmentHeightAbove(
			mad(middle, direction, origin)) <= 0;

		if (middleSubmerged == lowSideSubmerged)
		{
			firstLow = middle;
		}
		else
		{
			firstHigh = middle;
		}
	}

	const float firstCrossing = 0.5 * (firstLow + firstHigh);

	result.crosses = true;
	result.entry = firstCrossing;
	result.entryPoint = mad(firstCrossing, direction, origin);
	result.entryNormal = WaterSegmentSurfaceNormal(result.entryPoint);

	// One crossing, so the segment changes side exactly once and the first is
	// the only one there is.
	[branch]
	if (originSubmerged != targetSubmerged)
	{
		result.submerged = originSubmerged
			? firstCrossing
			: segment - firstCrossing;

		return result;
	}

	// Both ends on the same side, and something was crossed: the segment goes
	// out and comes back. Refine the last bracket as well and take the span
	// between the two for whichever medium is NOT the one the ends are in.
	float lastLow = lastLo;
	float lastHigh = lastHi;
	const bool lastLowSubmerged = WaterSegmentHeightAbove(
		mad(lastLow, direction, origin)) <= 0;

	for (uint refineLast = 0; refineLast < WATER_SEGMENT_REFINEMENTS;
		++refineLast)
	{
		const float middle = 0.5 * (lastLow + lastHigh);
		const bool middleSubmerged = WaterSegmentHeightAbove(
			mad(middle, direction, origin)) <= 0;

		if (middleSubmerged == lastLowSubmerged)
		{
			lastLow = middle;
		}
		else
		{
			lastHigh = middle;
		}
	}

	const float span = max(0, 0.5 * (lastLow + lastHigh) - firstCrossing);

	result.submerged = originSubmerged ? segment - span : span;

	return result;
}

/**
 * Discards a fragment that belongs to the other side of the water.
 *
 * The ocean is the only transparent that writes depth, so the transparent pass
 * draws around it: what the water stands in front of goes down first, where the
 * surface's screen space refraction picks it up, and everything else after,
 * where the surface's depth cannot reject it. A billboard can belong to
 * **both** - part of it behind a wave, part of it clear - so it is issued in
 * both passes and each keeps only its own share. That is why this is per
 * fragment and not a per-object classification: classifying by origin makes a
 * sprite snap whole from one side to the other as it crosses, and it can never
 * be partly behind anything.
 *
 * **The question is about the segment, not the fragment.** A point above the
 * surface over its own position can still have a crest between it and the eye,
 * and what stands between them is what decides which pass can draw it: only the
 * first pass runs early enough for the water to refract and fog what it covers.
 * Asking whether the fragment is wet answers a different question, and answers
 * it correctly only on a flat sea.
 *
 * The two sides are strictly complementary - one boolean out of one function,
 * given the same two points in either pass - so no fragment is kept by both
 * (which would blend it twice) or dropped by both (which would leave a hole).
 * The phase is fixed for the same reason: see `WATER_SEGMENT_STABLE_PHASE`.
 *
 * Example usage:
 * @code
 * ClipToWaterSide(P, camera.IsWaterSideBeyond(), camera.IsWaterSideNear());
 * @endcode
 *
 * @param[in] position - World space position of the fragment.
 * @param[in] keepBeyond - Keep only what the water stands in front of.
 * @param[in] keepNear - Keep only what nothing but air separates from the eye.
 *
 * @note Does nothing unless exactly one side is selected, so a draw that is not
 *       split pays only the flag test.
 */
inline void ClipToWaterSide(
	in float3 position,
	in bool keepBeyond,
	in bool keepNear
)
{
	[branch]
	if (keepBeyond != keepNear)
	{
		bool beyond;

		[branch]
		if (GetCamera().IsWaterSegmentModel())
		{
			beyond = TraceWaterSegment(
				GetCamera().position,
				position,
				WATER_SEGMENT_STABLE_PHASE
			).crosses;
		}
		else
		{
			// The DRAWN surface over the fragment's own position. Judged
			// against waves the grid flattened away, a fragment would be
			// dropped by both sides or kept by both.
			beyond = WaterSegmentHeightAbove(position) <= 0;
		}

		if (keepBeyond != beyond)
		{
			discard;
		}
	}
}

#endif // WI_WATER_SEGMENT_HF
