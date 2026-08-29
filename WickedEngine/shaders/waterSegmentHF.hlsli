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
 * Finds crossings by walking the height pyramid instead of sampling a fixed
 * number of times along the segment.
 *
 * The march samples a set number of times whatever the stretch, so a wave
 * narrower than the spacing is missed - and a near-horizontal ray lies inside
 * the slab for hundreds of metres, which is where it fails. The walk skips a
 * cell only when the ray passes wholly above or wholly below everything the
 * surface reaches inside it, which holds at any length.
 *
 * Set to 0 to fall back to the march, for comparing the two.
 */
#define WATER_SEGMENT_HIERARCHICAL 1

/**
 * Cells the walk may enter before giving up.
 *
 * A ceiling on the work, not a resolution. Open water is crossed in a handful
 * of patch-sized strides because every skip widens the next one; the budget
 * only binds where a ray runs along many crests at once, which is the
 * near-horizontal case. Reaching it leaves whatever was found so far.
 */
static const uint WATER_SEGMENT_WALK_BUDGET = 128;

/**
 * Finest pyramid level the walk descends to.
 *
 * **Zero, and it cannot be raised to save work.** A leaf is answered by testing
 * which side of the surface each of its ends is on, so a sheet of water thinner
 * than the leaf is entered and left inside one cell, puts both ends on the same
 * side, and is not seen at all. That is the very failure this walk exists to
 * remove - a feature narrower than the sampling - reintroduced at the size of a
 * cell instead of the size of a step, and it draws as holes quantised to the
 * cell grid.
 *
 * Crests really do throw sheets that thin: a camera placed between two of them
 * is underwater, so they are water volume and not an artifact to be tuned away.
 *
 * The cost of walking to level zero is bounded by the bound itself, which is
 * the point of having one - a cell is only entered when the surface could
 * genuinely be crossed inside it.
 */
static const uint WATER_SEGMENT_WALK_MIN_LEVEL = 0;

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

	/**
	 * Cells the hierarchical walk entered, for diagnosis. 0 on the march.
	 *
	 * A search that reaches its budget returns what it had, which is a hole in
	 * the water - so knowing how close a view runs to the ceiling says whether
	 * the ceiling is the problem before it is raised.
	 */
	uint cellsUsed;
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
/**
 * Brackets the first and last surface crossing on a stretch of the segment.
 *
 * Walks the ocean's height pyramid cell by cell rather than sampling the
 * stretch a set number of times. A cell is skipped whole when the ray passes
 * above everything the surface reaches inside it or below everything - which is
 * true however long the stretch is, and is the guarantee a fixed sample count
 * cannot give. Twelve samples over the four hundred metres a near-horizontal
 * ray spends inside the wave slab land thirty-three metres apart, and step over
 * crests entirely.
 *
 * **Walked, not subdivided.** Halving a stretch and testing each half re-tests
 * ground already covered and costs by the LENGTH of the stretch; a walk costs
 * by the number of cells actually entered, which is what the sea contains
 * rather than how far away it is. It also knows exactly which cell it stands
 * in, so one point tap bounds that cell exactly where a subdivision must sample
 * several to cover a stretch that could straddle any of them.
 *
 * Climbing a level after every skip is what makes open water nearly free: an
 * empty stretch is crossed in a handful of patch-sized steps. Descending on a
 * possible crossing is what keeps it exact.
 *
 * References:
 * https://www.nvidia.com/docs/IO/8228/GDC2003_SummedAreaVariance.pdf
 *
 * @param[in] origin - Start of the whole segment, in world space.
 * @param[in] direction - Normalized direction along it.
 * @param[in] searchNear - Distance along the segment where the search starts.
 * @param[in] searchFar - Distance where it ends.
 * @param[out] firstLo - Near side of the bracket around the first crossing.
 * @param[out] firstHi - Far side of it. Negative where nothing was crossed.
 * @param[out] lastLo - Near side of the bracket around the last crossing.
 * @param[out] lastHi - Far side of it.
 */
inline void WaterSegmentBracketCrossings(
	in float3 origin,
	in float3 direction,
	in float searchNear,
	in float searchFar,
	out float firstLo,
	out float firstHi,
	out float lastLo,
	out float lastHi,
	out uint cellsUsed
)
{
	firstLo = -1;
	firstHi = -1;
	lastLo = -1;
	lastHi = -1;

	const ShaderOcean ocean = GetWeather().ocean;

	uint baseWidth;
	uint baseHeight;
	uint levelCount;
	bindless_textures_float2[
		descriptor_index(ocean.texture_heightHierarchy)].GetDimensions(
			0, baseWidth, baseHeight, levelCount);

	// Patch space, where one unit is one patch and the pyramid is one texel
	// wide at its coarsest. The walk is done unwrapped so the geometry stays
	// continuous; only the lookups wrap.
	const float2 originUV = origin.xz * ocean.patch_size_rcp;
	const float2 directionUV = direction.xz * ocean.patch_size_rcp;

	// **Started at a level that fits the segment, not at the coarsest.**
	//
	// The coarsest cell covers the whole patch, so for any segment shorter than
	// that it always straddles the ray and the walk spends its first steps
	// descending without advancing - one per level, every call, before the
	// search begins. Beginning where the cell is about the size of the segment
	// skips that descent entirely.
	//
	// Never finer than the segment needs either: a coarser start is free, since
	// an empty stretch is skipped whole and the level climbs back up anyway.
	const float3 searchStart = mad(searchNear, direction, origin);
	const float3 searchEnd = mad(searchFar, direction, origin);

	uint level = min(
		ocean_height_bounds_level(length(searchEnd.xz - searchStart.xz)),
		levelCount - 1);

	float along = searchNear;

	cellsUsed = 0;

	uint budget = WATER_SEGMENT_WALK_BUDGET;

	[loop]
	while (along < searchFar && budget > 0)
	{
		--budget;
		++cellsUsed;

		const uint dimension = max(baseWidth >> level, 1u);
		const float2 hereUV = originUV + directionUV * along;

		// How far to the cell's far side, on whichever axis leaves first. A
		// component going nowhere never leaves, hence the guard.
		const float2 cellNext = (floor(hereUV * dimension)
			+ select(directionUV > 0, 1.0, 0.0)) / dimension;

		const float2 toEdge = select(
			abs(directionUV) > 1e-9,
			(cellNext - hereUV) / directionUV,
			FLT_MAX);

		// **Nudged past the boundary by a fraction of the cell, never by a
		// fixed distance.** A ray running almost parallel to a cell edge
		// leaves it a hair beyond where it entered, and a fixed nudge that is
		// small next to a coarse cell is not enough to carry such a ray across
		// at all - it stalls on the seam, taking near-zero steps until the
		// budget runs out, which draws as a line down the screen at whichever
		// direction happens to be axis aligned in patch space.
		//
		// A minimum step sized by the cell removes the stall outright: no
		// matter how grazing the crossing, the walk always leaves.
		const float cellWorld = rcp(max(ocean.patch_size_rcp, 1e-6))
			/ (float)dimension;

		const float cellSpan = max(min(toEdge.x, toEdge.y), 0.0);
		const float exitAlong = min(
			along + max(cellSpan, cellWorld * 0.01) + cellWorld * 1e-3,
			searchFar);

		const float3 pointNear = mad(along, direction, origin);
		const float3 pointFar = mad(exitAlong, direction, origin);

		const float2 bounds = ocean_height_bounds_uv(hereUV, level);

		const bool clearOfCell =
			min(pointNear.y, pointFar.y) > bounds.x
			|| max(pointNear.y, pointFar.y) < bounds.y;

		[branch]
		if (clearOfCell)
		{
			// Nothing here. Step past it and widen the stride, so an empty sea
			// is crossed in patch-sized strides instead of texel-sized ones.
			along = exitAlong;
			level = min(level + 1, levelCount - 1);
			continue;
		}

		[branch]
		if (level > WATER_SEGMENT_WALK_MIN_LEVEL)
		{
			// A crossing is possible somewhere in here. Ask a finer cell the
			// same question without moving.
			--level;
			continue;
		}

		// Finest cell, and the ray could be crossing inside it. Only the ends
		// can say whether it does; the bound said it was possible, not that it
		// happened.
		const bool nearSubmerged = WaterSegmentHeightAbove(pointNear) <= 0;
		const bool farSubmerged = WaterSegmentHeightAbove(pointFar) <= 0;

		if (nearSubmerged != farSubmerged)
		{
			if (firstLo < 0)
			{
				firstLo = along;
				firstHi = exitAlong;
			}

			lastLo = along;
			lastHi = exitAlong;
		}
		else
		{
			// The bound allowed a crossing and the ends disagree with it, so
			// there is nothing here after all - which is the same news as a
			// skip and earns the same widening. Without this the walk grinds
			// along at the finest level for as far as the bound stays
			// hopeful, which along a crest is most of the way.
			level = min(level + 1, levelCount - 1);
		}

		along = exitAlong;
	}
}

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
	result.cellsUsed = 0;

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

#if WATER_SEGMENT_HIERARCHICAL
	WaterSegmentBracketCrossings(
		origin, direction, slabNear, marchFar,
		firstLo, firstHi, lastLo, lastHi, result.cellsUsed);
#else
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
#endif // WATER_SEGMENT_HIERARCHICAL

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
		// Walked along the segment from the eye. Asking instead which side of
		// the surface the fragment's own position lies on cannot see a crest
		// standing between the two, so a fragment behind one is judged dry and
		// sorted into the near pass, where nothing can refract it.
		const bool beyond = TraceWaterSegment(
			GetCamera().position,
			position,
			WATER_SEGMENT_STABLE_PHASE
		).crosses;

		if (keepBeyond != beyond)
		{
			discard;
		}
	}
}

#endif // WI_WATER_SEGMENT_HF
