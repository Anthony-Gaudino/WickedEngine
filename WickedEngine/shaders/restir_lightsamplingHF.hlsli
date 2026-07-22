#ifndef WI_RESTIR_LIGHTSAMPLING_HF
#define WI_RESTIR_LIGHTSAMPLING_HF
#include "ShaderInterop_ReSTIR.h"

/**
 * ReSTIR light-sampling library.
 *
 * Shared, pass-agnostic building blocks used by every ReSTIR consumer:
 *   - weighted reservoir sampling (update / merge / unbiased weight),
 *   - resolving an analytic ShaderEntity light into a point sample,
 *   - the diffuse target function used to importance-sample lights,
 *   - a Resampled Importance Sampling (RIS) driver over the analytic light
 *     list.
 *
 * Include order (matches surfel_raytraceCS / ddgi_raytraceCS): include
 * "globals.hlsli" and "lightingHF.hlsli" first, then this header. Those provide
 * load_entity(), lights(), the RNG type, sample_hemisphere_cos(),
 * attenuation_pointlight()/attenuation_spotlight(), rotate_vector() and the
 * atmosphere helpers used below.
 *
 * The RIS driver here samples uniformly from the analytic light list and keeps
 * the winning sample in registers (no reservoir storage). With one candidate it
 * reduces exactly to the estimator SurfelGI/DDGI use today (light_count *
 * radiance * NdotL), so it is a safe drop-in that only improves as the
 * candidate count grows. Emissive-triangle lights and pre-sampled light tiles
 * plug into the same reservoir via RESTIRReservoirUpdate() and are added
 * alongside the light pre-sampling pass.
 */

/*
################################################################################
Reservoir operations
################################################################################
*/

/**
 * Streams one candidate into a reservoir using weighted reservoir sampling.
 *
 * @param[in,out] r - Reservoir to update.
 * @param[in] lightIndex - Candidate light reference (see RESTIRLightRef).
 * @param[in] uv - Sample parameterization to store for later reconstruction.
 * @param[in] targetPdf - Target function value p_hat of the candidate.
 * @param[in] risWeight - Resampling weight p_hat / p_source of the candidate.
 * @param[in,out] rng - Random generator.
 *
 * @return true if the candidate became the reservoir's selected sample.
 */
inline bool RESTIRReservoirUpdate(
	inout RESTIRReservoir r,
	uint lightIndex,
	float2 uv,
	float targetPdf,
	float risWeight,
	inout RNG rng)
{
	r.weightSum += risWeight;
	r.M += 1;

	// Select this candidate with probability risWeight / weightSum. Guarded so a
	// zero-weight run never selects and never divides by zero.
	if (r.weightSum > 0 && rng.next_float() * r.weightSum < risWeight)
	{
		r.lightIndex = lightIndex;
		r.uv = uv;
		r.targetPdf = targetPdf;
		return true;
	}
	return false;
}

/**
 * Merges reservoir 'other' into 'self' (spatial/temporal reuse).
 *
 * The other reservoir's sample is re-weighted by its target function evaluated
 * at self's surface, so the merge stays unbiased across differing shading
 * points.
 *
 * @param[in,out] self - Reservoir receiving the merge (M and weightSum grow).
 * @param[in] other - Reservoir to merge in.
 * @param[in] targetPdfAtSelf - other's sample target function evaluated at
 *   self's surface point/normal.
 * @param[in,out] rng - Random generator.
 */
inline void RESTIRReservoirMerge(
	inout RESTIRReservoir self,
	RESTIRReservoir other,
	float targetPdfAtSelf,
	inout RNG rng)
{
	const float risWeight =
		targetPdfAtSelf * RESTIRReservoirGetInvPdf(other) * other.M;

	self.M += other.M;
	self.weightSum += risWeight;

	if (self.weightSum > 0 && rng.next_float() * self.weightSum < risWeight)
	{
		self.lightIndex = other.lightIndex;
		self.uv = other.uv;
		self.targetPdf = targetPdfAtSelf;
	}
}

/*
################################################################################
Analytic light sampling
################################################################################
*/

/**
 * Diffuse target function: luminance of the unshadowed diffuse contribution.
 *
 * Albedo is deliberately factored out (it is applied later, in the consumer's
 * BRDF), which keeps the reservoir reusable across surfaces and matches the
 * forward+ integration. An unsamplable sample (radiance 0 or back-facing)
 * yields 0 and can never win the reservoir.
 *
 * @param[in] s - The resolved light sample.
 * @param[in] N - World-space shading normal.
 *
 * @return The scalar target function value p_hat (>= 0).
 */
inline float RESTIRTargetFunction(RESTIRLightSample s, float3 N)
{
	const float NdotL = saturate(dot(s.direction, N));
	const float luma = dot(s.radiance, float3(0.2126, 0.7152, 0.0722));
	return luma * NdotL;
}

/**
 * Emitted-power proxy of an analytic light, used to importance-sample the
 * pre-sampled light tiles.
 *
 * A cheap, view-independent scalar proportional to the light's contribution:
 * the luminance of its (energy-scaled) color. It is only a *sampling* weight -
 * the per-pixel estimator stays unbiased for any positive proxy because the
 * tile carries the exact reciprocal source pdf - so an approximate proxy just
 * trades variance, never correctness. A small floor keeps every light reachable
 * (a zero would make a light unsamplable and bias the estimate dark).
 *
 * @param[in] light - The analytic light entity.
 *
 * @return A positive power proxy (>= a small epsilon).
 */
inline float RESTIRLightPower(ShaderEntity light)
{
	const float luma = dot(light.GetColor().rgb, float3(0.2126, 0.7152, 0.0722));
	return max(luma, 1e-3);
}

/**
 * Cosine-weighted hemisphere sample from an explicit 2D sample point.
 *
 * The deterministic counterpart of sample_hemisphere_cos: reproducible from a
 * stored uv so a light sample can be re-resolved at a moving light's current
 * transform.
 *
 * @param[in] normal - Hemisphere axis.
 * @param[in] uv - Sample point in [0, 1)^2.
 *
 * @return A direction in the hemisphere around normal.
 */
inline float3 RESTIRHemisphereCos(float3 normal, float2 uv)
{
	return mul(hemispherepoint_cos(uv.x, uv.y), get_tangentspace(normal));
}

/**
 * Resolves an analytic ShaderEntity light into a point sample at point P.
 *
 * Deterministic form: the area extent (sphere radius, rectangle) is sampled
 * from the explicit uv, so the same sample can be reproduced at the light's
 * current transform (used to re-resolve a reused sample so a moving light's
 * shadow does not lag). The tube-length dimension is not parameterized by uv
 * and is sampled at the center. The returned radiance includes distance/spot
 * attenuation but NOT visibility. solidAnglePdf is left at 1: the
 * light-selection probability is carried by the RIS source pdf in the driver,
 * not here.
 *
 * @param[in] light - The analytic light entity.
 * @param[in] P - World-space shading point.
 * @param[in] N - World-space shading normal (unused here; kept for symmetry).
 * @param[in] uv - Sample parameterization on the light's area, in [0, 1)^2.
 *
 * @return The resolved sample; radiance is 0 when the light does not reach P.
 */
inline RESTIRLightSample RESTIRResolveAnalyticLight(
	ShaderEntity light,
	float3 P,
	float3 N,
	float2 uv)
{
	RESTIRLightSample s;
	s.direction = float3(0, 0, 1);
	s.distance = 0;
	s.radiance = 0;
	s.solidAnglePdf = 1;

	switch (light.GetType())
	{
	case ENTITY_TYPE_DIRECTIONALLIGHT:
	{
		float3 L = light.GetDirection().xyz;
		L += RESTIRHemisphereCos(L, uv) * light.GetRadius();
		L = normalize(L);
		s.direction = L;
		s.distance = FLT_MAX;

		float3 color = light.GetColor().rgb;
		[branch]
		if (GetFrame().options & OPTION_BIT_REALISTIC_SKY)
		{
			color *= GetAtmosphericLightTransmittance(GetWeather().atmosphere, P, L, texture_transmittancelut);
		}
		s.radiance = color;
	}
	break;
	case ENTITY_TYPE_POINTLIGHT:
	{
		float3 pos = light.position;
		pos += RESTIRHemisphereCos(normalize(pos - P), uv) * light.GetRadius();
		float3 L = pos - P;
		const float dist2 = dot(L, L);
		const float range = light.GetRange();

		[branch]
		if (dist2 > 0 && dist2 < range * range)
		{
			const float dist = sqrt(dist2);
			s.direction = L / dist;
			s.distance = dist;
			s.radiance = light.GetColor().rgb * attenuation_pointlight(dist2, range, light.GetRange2Rcp());
		}
	}
	break;
	case ENTITY_TYPE_SPOTLIGHT:
	{
		const float3 Loriginal = normalize(light.position - P);
		float3 pos = light.position;
		pos += RESTIRHemisphereCos(normalize(light.position - P), uv) * light.GetRadius();
		float3 L = pos - P;
		const float dist2 = dot(L, L);
		const float range = light.GetRange();

		[branch]
		if (dist2 > 0 && dist2 < range * range)
		{
			const float dist = sqrt(dist2);
			const float3 Ln = L / dist;
			const float spot_factor = dot(Loriginal, light.GetDirection());
			const float spot_cutoff = light.GetConeAngleCos();

			[branch]
			if (spot_factor > spot_cutoff)
			{
				s.direction = Ln;
				s.distance = dist;
				s.radiance = light.GetColor().rgb * attenuation_spotlight(dist2, range, light.GetRange2Rcp(), spot_factor, light.GetAngleScale(), light.GetAngleOffset());
			}
		}
	}
	break;
	case ENTITY_TYPE_RECTLIGHT:
	{
		const half4 quaternion = light.GetQuaternion();
		const half3 right = rotate_vector(half3(1, 0, 0), quaternion);
		const half3 up = rotate_vector(half3(0, 1, 0), quaternion);
		const half3 forward = cross(up, right);

		[branch]
		if (dot(P - light.position, forward) > 0) // in front of the light
		{
			float3 pos = light.position;
			pos += right * (uv.x - 0.5) * light.GetLength();
			pos += up * (uv.y - 0.5) * light.GetHeight();
			float3 L = pos - P;
			const float dist2 = dot(L, L);
			const float range = light.GetRange();

			[branch]
			if (dist2 > 0 && dist2 < range * range)
			{
				const float dist = sqrt(dist2);
				s.direction = L / dist;
				s.distance = dist;
				s.radiance = light.GetColor().rgb * attenuation_pointlight(dist2, range, light.GetRange2Rcp());
			}
		}
	}
	break;
	}

	return s;
}

/**
 * Draws a fresh uv and resolves the light (stochastic convenience overload).
 *
 * Kept for the light-only consumers (SurfelGI / DDGI) that do not persist a
 * sample and just need a random one each call.
 *
 * @param[in] light - The analytic light entity.
 * @param[in] P - World-space shading point.
 * @param[in] N - World-space shading normal.
 * @param[in,out] rng - Random generator.
 *
 * @return The resolved sample; radiance is 0 when the light does not reach P.
 */
inline RESTIRLightSample RESTIRResolveAnalyticLight(
	ShaderEntity light,
	float3 P,
	float3 N,
	inout RNG rng)
{
	const float2 uv = float2(rng.next_float(), rng.next_float());
	return RESTIRResolveAnalyticLight(light, P, N, uv);
}

/*
################################################################################
RIS driver (analytic light list)
################################################################################
*/

/**
 * Resamples the analytic light list into a reservoir at surface point (P, N).
 *
 * Draws candidateCount lights uniformly, resamples them by the diffuse target
 * function, and returns the reservoir together with the resolved winning
 * sample. The winning sample is kept in registers (nothing is stored), so the
 * caller must use it immediately:
 *
 * @code
 * RESTIRLightSample winning;
 * RESTIRReservoir r = RESTIRSampleLightsUniform(P, N, 8, false, rng, winning);
 * // trace one visibility ray for winning.direction up to winning.distance ...
 * float NdotL = saturate(dot(winning.direction, N));
 * float3 diffuse =
 *     winning.radiance * NdotL * visibility * RESTIRReservoirGetInvPdf(r);
 * @endcode
 *
 * @param[in] P - World-space shading point.
 * @param[in] N - World-space shading normal.
 * @param[in] candidateCount - Number of RIS candidates (>= 1).
 * @param[in] staticOnly - When true, only lights flagged static contribute
 *   (matches SurfelGI's direct-cache block, which bakes static lights only).
 *   Non-static candidates are still counted (uniform selection is preserved, so
 *   the estimator stays unbiased) but carry zero weight.
 * @param[in,out] rng - Random generator.
 * @param[out] winning - Resolved sample of the reservoir's selected light.
 *
 * @return The reservoir; RESTIRReservoirGetInvPdf(r) is the unbiased weight W.
 */
inline RESTIRReservoir RESTIRSampleLightsUniform(
	float3 P,
	float3 N,
	uint candidateCount,
	bool staticOnly,
	inout RNG rng,
	out RESTIRLightSample winning)
{
	RESTIRReservoir r = RESTIRReservoirInit();

	winning.direction = float3(0, 0, 1);
	winning.distance = 0;
	winning.radiance = 0;
	winning.solidAnglePdf = 1;

	const ShaderEntityIterator iterator = lights();
	const uint lightCount = iterator.item_count();
	if (lightCount == 0)
		return r;

	// Uniform selection over the analytic light list.
	const float invSourcePdf = (float)lightCount;

	for (uint i = 0; i < candidateCount; ++i)
	{
		const uint lightIndex = iterator.first_item() + rng.next_uint(lightCount);
		const ShaderEntity light = load_entity(lightIndex);

		RESTIRLightSample s;
		s.direction = float3(0, 0, 1);
		s.distance = 0;
		s.radiance = 0;
		s.solidAnglePdf = 1;
		float targetPdf = 0;

		[branch]
		if (!staticOnly || light.IsStaticLight())
		{
			s = RESTIRResolveAnalyticLight(light, P, N, rng);
			targetPdf = RESTIRTargetFunction(s, N);
		}

		const float risWeight = targetPdf * invSourcePdf;

		if (RESTIRReservoirUpdate(r, lightIndex, float2(0, 0), targetPdf, risWeight, rng))
		{
			winning = s;
		}
	}

	return r;
}

/*
################################################################################
RIS driver (pre-sampled light tiles)
################################################################################
*/

/**
 * Loads one pre-sampled light reference from the tile buffer (raw storage).
 *
 * @param[in] tileBuffer - Bindless SRV index of the RESTIRLightRef tile buffer.
 * @param[in] index - Flat entry index in [0, RESTIR_LIGHT_TILE_TOTAL).
 *
 * @return The stored light reference (index + reciprocal source pdf).
 */
inline RESTIRLightRef RESTIRLoadLightRef(int tileBuffer, uint index)
{
	// Each RESTIRLightRef is 8 bytes: uint lightIndex, float invSourcePdf.
	const uint2 raw =
		bindless_buffers[descriptor_index(tileBuffer)].Load2(index * 8);

	RESTIRLightRef ref;
	ref.lightIndex = raw.x;
	ref.invSourcePdf = asfloat(raw.y);
	return ref;
}

/**
 * Picks which pre-sampled tile a pixel draws its candidates from.
 *
 * An 8x8 pixel block shares one tile so a wave's per-pixel candidate reads stay
 * within a single cached tile (the point of pre-sampling), and the choice
 * rotates per frame so temporal reuse still sees the whole light set over time.
 *
 * @param[in] pixel - Screen pixel (primary-surface pixel for GI bounce
 *                    shading).
 * @param[in] frameIndex - Current frame counter.
 *
 * @return A tile index in [0, RESTIR_LIGHT_TILE_COUNT).
 */
inline uint RESTIRSelectTile(uint2 pixel, uint frameIndex)
{
	const uint block = (pixel.x >> 3) + (pixel.y >> 3) * 71u;
	return (block + frameIndex * 37u) % RESTIR_LIGHT_TILE_COUNT;
}

/**
 * Resamples analytic lights into a reservoir at (P, N), drawing candidates from
 * a pre-sampled light tile when one is available (falls back to uniform).
 *
 * With a valid tile the candidates are already importance-sampled by emitted
 * power, so fewer of them land on lights that cannot brighten this surface; the
 * tile also keeps the candidate reads cache-coherent across a wave. Each tile
 * entry carries the reciprocal of the (power-proportional) probability with
 * which it was pre-selected, so the estimator stays unbiased. With tileBuffer <
 * 0 this is exactly RESTIRSampleLightsUniform.
 *
 * The winning sample is kept in registers (nothing is stored); use it
 * immediately (see RESTIRSampleLightsUniform for the pattern).
 *
 * @param[in] P - World-space shading point.
 * @param[in] N - World-space shading normal.
 * @param[in] candidateCount - Number of RIS candidates (>= 1).
 * @param[in] tileBuffer - Bindless SRV index of the tile buffer (-1 = uniform).
 * @param[in] pixel - Screen pixel selecting the tile (ignored when uniform).
 * @param[in] frameIndex - Current frame counter (ignored when uniform).
 * @param[in] staticOnly - When true, only static lights carry weight (see
 *   RESTIRSampleLightsUniform).
 * @param[in,out] rng - Random generator.
 * @param[out] winning - Resolved sample of the reservoir's selected light.
 *
 * @return The reservoir; RESTIRReservoirGetInvPdf(r) is the unbiased weight W.
 */
inline RESTIRReservoir RESTIRSampleLights(
	float3 P,
	float3 N,
	uint candidateCount,
	int tileBuffer,
	uint2 pixel,
	uint frameIndex,
	bool staticOnly,
	inout RNG rng,
	out RESTIRLightSample winning)
{
	[branch]
	if (tileBuffer < 0)
	{
		return RESTIRSampleLightsUniform(
			P, N, candidateCount, staticOnly, rng, winning);
	}

	RESTIRReservoir r = RESTIRReservoirInit();

	winning.direction = float3(0, 0, 1);
	winning.distance = 0;
	winning.radiance = 0;
	winning.solidAnglePdf = 1;

	const uint tileBase =
		RESTIRSelectTile(pixel, frameIndex) * RESTIR_LIGHT_TILE_SIZE;

	for (uint i = 0; i < candidateCount; ++i)
	{
		const uint slot = rng.next_uint(RESTIR_LIGHT_TILE_SIZE);
		const RESTIRLightRef ref = RESTIRLoadLightRef(tileBuffer, tileBase + slot);

		RESTIRLightSample s;
		s.direction = float3(0, 0, 1);
		s.distance = 0;
		s.radiance = 0;
		s.solidAnglePdf = 1;
		float targetPdf = 0;
		float2 uv = float2(0, 0);

		// Only analytic lights are pre-sampled for now; an emissive-triangle
		// reference (flag set) or an empty slot contributes nothing but is
		// still counted, keeping the reservoir statistics consistent.
		[branch]
		if (ref.lightIndex != RESTIR_INVALID_LIGHT_INDEX &&
			(ref.lightIndex & RESTIR_LIGHT_FLAG_EMISSIVE_TRIANGLE) == 0)
		{
			const ShaderEntity light = load_entity(ref.lightIndex);
			if (!staticOnly || light.IsStaticLight())
			{
				uv = float2(rng.next_float(), rng.next_float());
				s = RESTIRResolveAnalyticLight(light, P, N, uv);
				targetPdf = RESTIRTargetFunction(s, N);
			}
		}

		const float risWeight = targetPdf * ref.invSourcePdf;

		if (RESTIRReservoirUpdate(
				r, ref.lightIndex, uv, targetPdf, risWeight, rng))
		{
			winning = s;
		}
	}

	return r;
}

/*
################################################################################
ReGIR (Reservoir-based Grid Importance Resampling)
################################################################################
*/

/**
 * Importance weight of a light for a whole grid cell (a sphere of the given
 * radius about its center), used to build the ReGIR cells.
 *
 * Unlike the global power proxy (light color luminance, position-blind), this
 * is the light's approximate irradiance AT the cell: a local light is weighted
 * by power / distance^2 so a dim light near the cell dominates that cell's set,
 * while a bright light far away contributes little - exactly the mismatch that
 * makes global power sampling noisy for many disparate lights. The final
 * per-pixel target function still resolves the exact contribution; this only
 * steers which lights are cached in the cell.
 *
 * **Directional lights are excluded** (weight 0). They reach every cell
 * equally, so a position-aware cell buys them nothing - but a bright
 * directional (a sun) would then win a constant, high weight in *every* cell
 * and crowd the local lights out of the far cells where their own weight is
 * small. The compensating reservoir weight the local light then needs exceeds
 * the power-based firefly clamp and gets chopped, so the local light reads dim
 * away from its source. Directionals are sampled directly from the global tiles
 * instead (see the DI initial pass), where their weight is power-bounded and
 * never clamped.
 *
 * @param[in] light - Analytic light entity.
 * @param[in] cellCenter - World-space cell center.
 * @param[in] cellRadius - Cell bounding-sphere radius (world units).
 *
 * @return A non-negative importance weight (0 if the light cannot reach the
 *   cell, or is directional).
 */
inline float RESTIRReGIRVolumeWeight(
	ShaderEntity light, float3 cellCenter, float cellRadius)
{
	[branch]
	if (light.GetType() == ENTITY_TYPE_DIRECTIONALLIGHT)
		return 0; // sampled via the global tiles, not position-aware cells

	const float luma =
		max(dot(light.GetColor().rgb, float3(0.2126, 0.7152, 0.0722)), 1e-3);

	// Local light: irradiance ~ power / distance^2, range-limited.
	const float3 toCell = cellCenter - light.position;
	const float dist = sqrt(dot(toCell, toCell));

	[branch]
	if (dist - cellRadius > light.GetRange())
		return 0;

	// Never treat the light as closer than the cell radius, so a cell
	// containing the light does not get an unbounded (firefly) weight.
	const float avgDist = max(dist, cellRadius);
	return luma / (avgDist * avgDist);
}

/**
 * Positive modulo for each component of an int3 (result in [0, n)).
 *
 * HLSL `%` follows the sign of the dividend, so a plain `a % n` yields negative
 * results for negative world-cell coordinates; this wraps them into [0, n) as
 * the toroidal buffer addressing requires.
 *
 * @param[in] a - Dividend (may be negative).
 * @param[in] n - Positive modulus.
 *
 * @return `a` reduced modulo `n`, per component, in [0, n).
 */
inline int3 RESTIRReGIRPosMod(int3 a, int n)
{
	return ((a % n) + n) % n;
}

/**
 * Maps a world position to its (nearest) camera-centered grid cell.
 *
 * A world position maps to world cell floor(worldPos / cellSize); only cells
 * inside the window [gridOriginCell, gridOriginCell + RESTIR_REGIR_GRID_RES)
   are cached. The returned index is the **toroidal** buffer slot (world cell
 * modulo the grid resolution), so a world cell keeps the same slot as the
 * window scrolls and its reservoir accumulates across frames.
 *
 * The caller (RESTIRReGIRSampleCell) jitters `worldPos` by ~a cell before this
 * lookup, which is what dithers the hard cell boundaries into
 * denoiser-removable noise instead of a low-frequency blob (the same role the
 * old stochastic trilinear pick served, but with a tunable, unbounded spread).
 *
 * @param[in] worldPos - World-space shading point (already jittered).
 * @param[in] gridOriginCell - World-cell coordinate of buffer cell (0,0,0).
 * @param[in] cellSize - Cell edge length (world units).
 *
 * @return Toroidal cell index in [0, RESTIR_REGIR_CELL_COUNT), or -1 if outside
 *   the camera-centered window.
 */
inline int RESTIRReGIRGridWorldToCell(
	float3 worldPos, int3 gridOriginCell, float cellSize)
{
	const int3 worldCell = int3(floor(worldPos / cellSize));
	const int res = (int)RESTIR_REGIR_GRID_RES;

	[branch]
	if (any(worldCell < gridOriginCell) ||
		any(worldCell >= gridOriginCell + res))
		return -1;

	const int3 cell = RESTIRReGIRPosMod(worldCell, res);
	return cell.x + cell.y * res + cell.z * (res * res);
}

/**
 * Returns the world cell and world-space geometry of a toroidal ReGIR slot.
 *
 * Given a buffer cell index and the current window origin, recovers the world
 * cell that slot currently caches (the unique cell congruent to the buffer
 * coordinate modulo the resolution inside
 * [gridOriginCell, gridOriginCell + RESTIR_REGIR_GRID_RES)), and from it the
 * cell center and bounding-sphere radius. The build pass also uses the returned
 * world cell as the slot's toroidal identity for staleness detection.
 *
 * @param[in] cellIndex - Buffer cell index in [0, RESTIR_REGIR_CELL_COUNT).
 * @param[in] gridOriginCell - World-cell coordinate of buffer cell (0,0,0).
 * @param[in] cellSize - Cell edge length (world units).
 * @param[out] worldCell - World-cell coordinate this slot represents.
 * @param[out] center - World-space cell center.
 * @param[out] radius - Cell bounding-sphere radius (half the cell diagonal).
 */
inline void RESTIRReGIRCellToWorld(
	int cellIndex, int3 gridOriginCell, float cellSize,
	out int3 worldCell, out float3 center, out float radius)
{
	const int res = (int)RESTIR_REGIR_GRID_RES;
	const int3 bufferCoord = int3(
		cellIndex % res,
		(cellIndex / res) % res,
		cellIndex / (res * res));
	// Invert the toroidal map: the world cell is the representative of
	// bufferCoord (mod res) that lies inside the current window.
	worldCell =
		gridOriginCell + RESTIRReGIRPosMod(bufferCoord - gridOriginCell, res);
	center = (float3(worldCell) + 0.5) * cellSize;
	radius = cellSize * 0.8660254; // sqrt(3)/2 * edge
}

/*
--- ReGIR onion mode --------------------------------------------------------
 *
 * Camera-centered concentric spherical shells. Shell radii grow geometrically
 * (r_i = cellSize * g^i, g = 1 + dTheta, dTheta = PI / elevationBands), so a
 * shell's radial thickness is proportional to its radius. Every shell shares
 * one spherical angular grid: `elevationBands` equal-height latitude bands,
 * each split into a cos(latitude)-reduced number of longitude cells so cells
 * stay roughly square. With a constant angular grid and radius-proportional
 * radial thickness, every cell subtends about the same solid angle and
 * therefore projects to a roughly constant screen size at any depth.
 *
 * The onion buffer is the prefix sum of per-band longitude counts:
 * bandOffset[b] = sum of longitude counts of bands < b, with
 * bandOffset[elevationBands] = cellsPerShell. A (band, longitude) pair maps to
 * a flat angular index bandOffset[band] + longitude, and back.
 *
 * Independent construction (geometric shells + reduced lat-long sphere
 * binning).
 */

/**
 * Reads a band prefix-sum entry (cumulative longitude count of bands < band).
 *
 * @param[in] onionBuffer - Bindless SRV of the band prefix-sum buffer (raw).
 * @param[in] band - Band index in [0, elevationBands].
 *
 * @return Cumulative longitude-cell count before this band.
 */
inline uint RESTIRReGIROnionBandOffset(int onionBuffer, uint band)
{
	return bindless_buffers[descriptor_index(onionBuffer)].Load(band * 4);
}

/**
 * Spherical decomposition used by the onion
 * (x = r cos(az) cos(el), y = r sin(el), z = r sin(az) cos(el)).
 *
 * @param[in] p - Point relative to the onion center.
 * @param[out] r - Radius.
 * @param[out] azimuth - Azimuth angle (radians).
 * @param[out] elevation - Elevation angle (radians, in [-pi/2, pi/2]).
 */
inline void RESTIRReGIROnionCartesianToSpherical(
	float3 p, out float r, out float azimuth, out float elevation)
{
	r = length(p);
	elevation = (r > 0) ? asin(clamp(p.y / r, -1.0, 1.0)) : 0.0;
	azimuth = atan2(p.z, p.x);
}

/**
 * Inverse of RESTIRReGIROnionCartesianToSpherical.
 *
 * @param[in] r - Radius.
 * @param[in] azimuth - Azimuth angle (radians).
 * @param[in] elevation - Elevation angle (radians).
 *
 * @return The Cartesian point (onion units).
 */
inline float3 RESTIRReGIROnionSphericalToCartesian(
	float r, float azimuth, float elevation)
{
	return float3(
		r * cos(azimuth) * cos(elevation),
		r * sin(elevation),
		r * sin(azimuth) * cos(elevation));
}

/**
 * Maps a world position to its onion cell index.
 *
 * Takes the position relative to the onion center, decomposes it into radius /
 * latitude / longitude, picks the geometric shell from the radius and the
 * latitude band + cos-reduced longitude cell from the shared angular grid, then
 * combines them into a flat cell index.
 *
 * @param[in] worldPos - World-space point.
 * @param[in] center - World-space onion center.
 * @param[in] cellSize - Inner shell radius / world scale.
 * @param[in] onionBuffer - Bindless SRV of the band prefix-sum buffer.
 * @param[in] elevationBands - Latitude band count.
 * @param[in] shells - Concentric shell count.
 * @param[in] cellsPerShell - Angular cells per shell.
 *
 * @return Onion cell index, or -1 if beyond the outermost shell.
 */
inline int RESTIRReGIROnionWorldToCell(
	float3 worldPos, float3 center, float cellSize,
	int onionBuffer, uint elevationBands, uint shells, uint cellsPerShell)
{
	float r, azimuth, elevation;
	RESTIRReGIROnionCartesianToSpherical(
		worldPos - center, r, azimuth, elevation);
	if (azimuth < 0)
		azimuth += 2 * PI; // normalize to [0, 2*PI)

	const float dTheta = PI / (float)elevationBands;
	const float g = 1.0 + dTheta; // shell growth ratio

	// Inner sphere is a single cell; beyond the outermost shell, fall back.
	[branch]
	if (r <= cellSize)
		return 0;
	const int shell = (int)floor(log(r / cellSize) / log(g));
	[branch]
	if (shell >= (int)shells)
		return -1;

	// Latitude band, then the cos-reduced longitude cell within it.
	const uint band = min(
		(uint)floor((elevation + 0.5 * PI) / dTheta), elevationBands - 1);
	const uint bandStart = RESTIRReGIROnionBandOffset(onionBuffer, band);
	const uint bandEnd = RESTIRReGIROnionBandOffset(onionBuffer, band + 1);
	const uint bandCount = max(bandEnd - bandStart, 1u);
	const uint longitude =
		min((uint)floor(azimuth / (2 * PI) * bandCount), bandCount - 1);

	return (int)(1 + (uint)shell * cellsPerShell + bandStart + longitude);
}

/**
 * Returns the world-space center and bounding radius of an onion cell.
 *
 * The inverse of RESTIRReGIROnionWorldToCell, used by the build to weight
 * lights by their contribution to the cell.
 *
 * @param[in] cellIndex - Onion cell index.
 * @param[in] center - World-space onion center.
 * @param[in] cellSize - Inner shell radius / world scale.
 * @param[in] onionBuffer - Bindless SRV of the band prefix-sum buffer.
 * @param[in] elevationBands - Latitude band count.
 * @param[in] cellsPerShell - Angular cells per shell.
 * @param[out] cellCenter - World-space cell center.
 * @param[out] cellRadius - Cell bounding-sphere radius (world units).
 *
 * @return true if the cell index is valid.
 */
inline bool RESTIRReGIROnionCellToWorld(
	int cellIndex, float3 center, float cellSize,
	int onionBuffer, uint elevationBands, uint cellsPerShell,
	out float3 cellCenter, out float cellRadius)
{
	cellCenter = center;
	cellRadius = cellSize;

	// Cell 0 is the inner sphere; a negative index is invalid.
	[branch]
	if (cellIndex <= 0)
		return cellIndex == 0;

	const float dTheta = PI / (float)elevationBands;
	const float g = 1.0 + dTheta;

	const uint c = (uint)(cellIndex - 1);
	const uint shell = c / cellsPerShell;
	const uint angular = c % cellsPerShell;

	// Find the latitude band whose prefix range contains this angular index.
	uint band = elevationBands - 1;
	[loop]
	for (uint bi = 0; bi < elevationBands; ++bi)
	{
		[branch]
		if (RESTIRReGIROnionBandOffset(onionBuffer, bi + 1) > angular)
		{
			band = bi;
			break;
		}
	}
	const uint bandStart = RESTIRReGIROnionBandOffset(onionBuffer, band);
	const uint bandEnd = RESTIRReGIROnionBandOffset(onionBuffer, band + 1);
	const uint bandCount = max(bandEnd - bandStart, 1u);
	const uint longitude = angular - bandStart;

	const float elevation = -0.5 * PI + ((float)band + 0.5) * dTheta;
	const float azimuth =
		((float)longitude + 0.5) * (2 * PI / (float)bandCount);

	const float rInner = cellSize * pow(g, (float)shell);
	const float rOuter = rInner * g;
	const float rMid = sqrt(rInner * rOuter);

	cellCenter =
		center + RESTIRReGIROnionSphericalToCartesian(rMid, azimuth, elevation);

	// Bounding radius ~ half the larger of the radial thickness and the local
	// tangential cell extent (radius * band angle).
	cellRadius = 0.5 * max(rOuter - rInner, rMid * dTheta);
	return true;
}

/**
 * Jitter magnitude (world units) for onion cell sampling at a point.
 *
 * Onion cells grow with distance, so the jitter that hides cell boundaries must
 * scale with radius; the local tangential cell extent is about radius * dTheta
 * (floored at the inner shell so the near field still jitters).
 *
 * @param[in] worldPos - World-space point.
 * @param[in] center - World-space onion center.
 * @param[in] cellSize - Inner shell radius / world scale.
 * @param[in] elevationBands - Latitude band count.
 *
 * @return A world-space jitter radius.
 */
inline float RESTIRReGIROnionJitterScale(
	float3 worldPos, float3 center, float cellSize, uint elevationBands)
{
	const float dTheta = PI / (float)elevationBands;
	const float r = max(length(worldPos - center), cellSize);
	return r * dTheta;
}

/**
 * Selects a ReGIR cell for a shading point in the active mode.
 *
 * Jitters the lookup position by a random offset scaled to the local cell size
 * (uniform for the grid, distance-growing for the onion) and then takes the
 * nearest cell. Jittering the position, rather than reading the exact cell,
 * dithers the hard cell/shell boundaries into per-pixel noise the denoiser can
 * remove - otherwise each cell's slightly different estimate prints through as
 * a low-frequency, cell-sized brightness blob (most visible when the camera is
 * close and a cell covers many pixels). `samplingJitter` scales the spread: 0
 * gives hard cells, 1 spreads over ~one cell, larger values reach further
 * neighbors to break up bigger blobs (at the cost of position accuracy).
 *
 * @param[in] mode - RESTIR_REGIR_MODE_GRID or RESTIR_REGIR_MODE_ONION.
 * @param[in] P - World-space shading point.
 * @param[in] gridOriginCell - Grid window origin (grid mode).
 * @param[in] onionCenter - Onion center (onion mode).
 * @param[in] cellSize - Cell size / onion inner-cell size.
 * @param[in] onionBuffer - Bindless SRV of the onion band table (onion mode).
 * @param[in] elevationBands - Onion latitude bands.
 * @param[in] shells - Onion concentric shell count.
 * @param[in] cellsPerShell - Onion angular cells per shell.
 * @param[in] samplingJitter - Lookup-jitter width in local cell sizes.
 * @param[in,out] rng - Random generator.
 *
 * @return Cell index, or -1 if the point is outside the structure.
 */
inline int RESTIRReGIRSampleCell(
	uint mode, float3 P, int3 gridOriginCell, float3 onionCenter, float cellSize,
	int onionBuffer, uint elevationBands, uint shells, uint cellsPerShell,
	float samplingJitter, inout RNG rng)
{
	// Local cell extent at P: uniform for the grid, growing with distance for
	// the onion (its cells are bigger far from the eye).
	const float localCell = (mode == RESTIR_REGIR_MODE_ONION)
		? RESTIRReGIROnionJitterScale(P, onionCenter, cellSize, elevationBands)
		: cellSize;

	// Random offset in [-0.5, 0.5]^3 cells, scaled by samplingJitter.
	const float3 jitter =
		(float3(rng.next_float(), rng.next_float(), rng.next_float()) - 0.5)
		* (samplingJitter * localCell);
	const float3 Pj = P + jitter;

	[branch]
	if (mode == RESTIR_REGIR_MODE_ONION)
		return RESTIRReGIROnionWorldToCell(Pj, onionCenter, cellSize,
			onionBuffer, elevationBands, shells, cellsPerShell);
	return RESTIRReGIRGridWorldToCell(Pj, gridOriginCell, cellSize);
}

/**
 * Loads a ReGIR cell reservoir slot as a light reference for the per-pixel RIS.
 *
 * Unpacks the accumulated RESTIRReGIRCell and forms its reciprocal source pdf
 * weightSum / (M * targetWeight) - the reservoir's unbiased contribution 
 * weight- so the per-pixel RIS can weight the candidate without re-deriving its
 * selection probability. An empty / not-yet-converged slot yields an invalid
 * reference.
 *
 * @param[in] regirBuffer - Bindless SRV index of the ReGIR cell buffer (raw).
 * @param[in] slot - Flat slot index in [0, RESTIR_REGIR_SLOT_TOTAL).
 *
 * @return The cached light reference (index + reciprocal source pdf).
 */
inline RESTIRLightRef RESTIRReGIRLoadSlot(int regirBuffer, uint slot)
{
	// Cells are 32 bytes (RESTIRReGIRCell); the reservoir fields the sampler
	// needs are the first uint4 (the second holds the toroidal world-cell key).
	const uint4 raw =
		bindless_buffers[descriptor_index(regirBuffer)].Load4(slot * 32);

	RESTIRLightRef ref;
	ref.lightIndex = raw.x;

	const float weightSum = asfloat(raw.y);
	const float m = asfloat(raw.z);
	const float targetWeight = asfloat(raw.w);
	ref.invSourcePdf =
		(m > 0 && targetWeight > 0) ? weightSum / (m * targetWeight) : 0;

	return ref;
}

#endif // WI_RESTIR_LIGHTSAMPLING_HF
