#define TEXTURE_SLOT_NONUNIFORM // shadow ray alpha-tests divergent materials
#include "globals.hlsli"
#include "raytracingHF.hlsli"
#include "lightingHF.hlsli"
#include "restir_lightsamplingHF.hlsli"
#include "restir_diHF.hlsli"
#include "restir_di_visibilityHF.hlsli"

// ReSTIR DI - initial candidate generation pass.
//
// Per screen pixel: reconstruct the primary surface from the depth + normal
// buffers, resample RESTIRDIPushConstants.candidateCount analytic lights into
// an initial reservoir, then trace a single shadow ray for the chosen sample
// and cache its visibility. The reservoir is written to the temporal history
// buffer, which the temporal pass consumes next.

PUSHCONSTANT(push, RESTIRDIPushConstants);

// Reservoir buffers are raw (ByteAddressBuffer) so the same resource can also
// be read bindlessly by the forward shading pass. Each reservoir is 32 bytes
// (RESTIRDIReservoirPacked = 2x uint4).
RWByteAddressBuffer reservoirOutput : register(u0);

/**
 * Streams one resolved analytic-light candidate into an initial DI reservoir
 * (shadowed RIS with a balance-heuristic resampling weight).
 *
 * Resolves the light at (P, N), evaluates the diffuse target, traces one shadow
 * ray, and - if the sample is visible - streams it with resampling weight
 * \[ w = \frac{\hat p \, \mathrm{vis}}{\mathrm{denom}} \], where `denom` is the
 * candidate's **combined sampling density** across all strategies. Encoding the
 * density here (rather than a single strategy's reciprocal pdf) is what lets
 * the local tile pool and the exhaustive directional sweep share one reservoir
 * without over-counting a light both can draw. Occluded or zero-target
 * candidates contribute nothing.
 *
 * @param[in] light - The candidate light entity (loaded from `slot`).
 * @param[in] slot - The light's current entity slot (stored as a stable index).
 * @param[in] uv - Area-sample parameterization on the light.
 * @param[in] denom - Combined sampling density of the candidate (must be > 0).
 * @param[in] P - World-space shading point.
 * @param[in] N - World-space shading normal.
 * @param[in,out] weightSum - Running balance-heuristic weight sum.
 * @param[in,out] r - Reservoir being built.
 * @param[in,out] rng - Random generator.
 */
inline void RESTIRDIStreamInitialCandidate(
	ShaderEntity light,
	uint slot,
	float2 uv,
	float denom,
	float3 P,
	float3 N,
	inout float weightSum,
	inout RESTIRDIReservoir r,
	inout RNG rng)
{
	[branch]
	if (denom <= 0)
		return;

	const RESTIRLightSample s = RESTIRResolveAnalyticLight(light, P, N, uv);
	const float unshadowedTarget = RESTIRTargetFunction(s, N);
	[branch]
	if (unshadowedTarget <= 0)
		return;

	// One shadow ray per candidate: occluded -> zero target -> excluded.
	const float sampleDist =
		(s.distance >= FLT_MAX * 0.5) ? 100000.0 : s.distance;
	const float3 samplePos = P + s.direction * sampleDist;
	const float vis = RESTIRDITraceVisibility(P, N, samplePos, rng);
	[branch]
	if (vis <= 0)
		return;

	const float risWeight = unshadowedTarget * vis / denom;
	weightSum += risWeight;

	if (weightSum > 0 && rng.next_float() * weightSum < risWeight)
	{
		r.samplePosition = samplePos;
		r.sampleRadiance = s.radiance;
		// Visible winner: shadowed target == unshadowed target.
		r.targetPdf = unshadowedTarget;
		// Store the light's STABLE scene index (not the volatile per-frame
		// entity slot) so reuse resolves the same physical light across frames.
		r.lightIndex = RESTIRDIStableLightIndex(slot);
		r.uv = uv;
		r.visibility = vis;
	}
}

/**
 * Combined sampling density of a light across the initial pass's strategies.
 *
 * The balance-heuristic denominator: the expected number of times the light is
 * drawn across the power-weighted tile pool (candidateCount draws), the uniform
 * pool (uniformCount draws), and - for a directional light - the exhaustive
 * sweep (one draw). Every candidate of a given light uses this same density, so
 * the mixture stays unbiased regardless of which strategy actually drew it.
 *
 * @param[in] light - The candidate light.
 * @param[in] candidateCount - Number of power-weighted tile draws.
 * @param[in] uniformCount - Number of uniform draws.
 * @param[in] lightCount - Total analytic light count.
 * @param[in] totalPower - Sum of RESTIRLightPower over all lights.
 * @param[in] hasTiles - Whether the power-weighted tiles are available.
 *
 * @return The combined sampling density (> 0 when any strategy can draw it).
 */
inline float RESTIRDIInitialDenom(
	ShaderEntity light,
	uint candidateCount,
	uint uniformCount,
	uint lightCount,
	float totalPower,
	bool hasTiles)
{
	// Tile selection pmf: power-proportional when tiles exist, else uniform.
	const float pTile = hasTiles
		? RESTIRLightPower(light) / totalPower
		: 1.0 / (float)lightCount;
	float d = candidateCount * pTile + uniformCount * (1.0 / (float)lightCount);
	[branch]
	if (light.GetType() == ENTITY_TYPE_DIRECTIONALLIGHT)
		d += 1.0;
	return d;
}

/**
 * Builds an initial DI reservoir by resampling the analytic light list.
 *
 * Draws candidateCount lights (from a pre-sampled tile when available, else
 * uniformly), resolves each to a world-space sample via the light-sampling
 * core, and resamples by the diffuse target function. The chosen sample's
 * visibility is left unset for the caller to fill after a shadow ray.
 *
 * @param[in] P - World-space shading point.
 * @param[in] N - World-space shading normal.
 * @param[in] candidateCount - Number of RIS candidates (>= 1).
 * @param[in] pixel - Screen pixel (selects the pre-sampled tile).
 * @param[in,out] rng - Random generator.
 *
 * @return The initial reservoir (visibility unset).
 */
inline RESTIRDIReservoir RESTIRDISampleInitial(
	float3 P,
	float3 N,
	uint candidateCount,
	uint2 pixel,
	inout RNG rng)
{
	RESTIRDIReservoir r = RESTIRDIReservoirInit();

	const ShaderEntityIterator iterator = lights();
	const uint lightCount = iterator.item_count();
	if (lightCount == 0)
		return r;

	// Shadowed RIS: draw candidates (from a pre-sampled power-weighted tile
	// when available) and resample by the **shadowed** target
	// target*visibility. A candidate that is occluded gets zero target, so it
	// can neither win nor add to weightSum - which is what stops an occluded
	// but high-target light (a bright point light behind a wall) from inflating
	// a *visible* light's weight and brightening its region. The winner is
	// therefore always visible, so its stored target equals its unshadowed
	// target and the balance-heuristic reuse stays consistent without
	// re-tracing every reused sample.
	const int tileBuffer = push.lightTileBuffer;
	const uint tileBase = (tileBuffer >= 0)
		? RESTIRSelectTile(pixel, push.frameIndex) * RESTIR_LIGHT_TILE_SIZE
		: 0u;

	// Directional (sun / infinite) lights are swept exhaustively below, so
	// their combined sampling density needs the analytic local-pool selection
	// pmf.
	const ShaderEntityIterator dirIterator = directional_lights();
	const uint dirCount = dirIterator.item_count();
	const float totalPower = max(GetFrame().totalLightPower, 1e-3);

	float weightSum = 0;
	float M = 0;

	// Local-light candidates from the power-weighted tile (or uniform
	// fallback). A directional light can also surface here, but because it is
	// additionally swept below, its resampling weight uses the density of BOTH
	// strategies (balance heuristic), so a rare tile draw of the dim sun no
	// longer spikes.
	for (uint i = 0; i < candidateCount; ++i)
	{
		// Draw one candidate light from the power-weighted tile (its resampling
		// weight is reweighted analytically by RESTIRDIInitialDenom, so only
		// the tile's light selection is used here, not its stored
		// invSourcePdf).
		uint lightIndex;
		[branch]
		if (tileBuffer >= 0)
		{
			const uint slot = rng.next_uint(RESTIR_LIGHT_TILE_SIZE);
			const RESTIRLightRef ref =
				RESTIRLoadLightRef(tileBuffer, tileBase + slot);
			lightIndex = ref.lightIndex;
		}
		else
		{
			lightIndex = iterator.first_item() + rng.next_uint(lightCount);
		}

		M += 1;

		// Only analytic lights are sampled here (skip empty / emissive-tri
		// refs).
		[branch]
		if (lightIndex == RESTIR_INVALID_LIGHT_INDEX ||
			(lightIndex & RESTIR_LIGHT_FLAG_EMISSIVE_TRIANGLE) != 0)
			continue;

		const ShaderEntity light = load_entity(lightIndex);

		// Combined sampling density across the tile, uniform and directional
		// strategies (balance heuristic). Analytic (power-proportional) so it
		// is identical whether this light is drawn from the tile or the uniform
		// pool below - the tile's stored invSourcePdf is only its sampling
		// mechanism.
		const float denom = RESTIRDIInitialDenom(
			light, candidateCount, RESTIR_DI_UNIFORM_CANDIDATES,
			lightCount, totalPower, tileBuffer >= 0);

		const float2 uv = float2(rng.next_float(), rng.next_float());
		RESTIRDIStreamInitialCandidate(
			light, lightIndex, uv, denom, P, N, weightSum, r, rng);
	}

	// Uniform-light candidates: every light gets a power-independent chance, so
	// a dim light is sampled even when a much brighter light floods the tiles
	// (which importance-sample by global power). MIS-combined with the tiles
	// via the shared density above, so it stays unbiased and the tiles still
	// handle the many-lights case.
	for (uint u = 0; u < RESTIR_DI_UNIFORM_CANDIDATES; ++u)
	{
		const uint lightIndex =
			iterator.first_item() + rng.next_uint(lightCount);

		M += 1;

		const ShaderEntity light = load_entity(lightIndex);
		const float denom = RESTIRDIInitialDenom(
			light, candidateCount, RESTIR_DI_UNIFORM_CANDIDATES,
			lightCount, totalPower, tileBuffer >= 0);

		const float2 uv = float2(rng.next_float(), rng.next_float());
		RESTIRDIStreamInitialCandidate(
			light, lightIndex, uv, denom, P, N, weightSum, r, rng);
	}

	// Directional (sun / infinite) lights: sampled exhaustively every frame,
	// separate from the power-weighted local pool. A directional light reaches
	// every pixel but carries a small power proxy (irradiance, not intensity),
	// so the local pool tiles it rarely and the estimator would only
	// occasionally select the sun - with a large compensating weight, a
	// per-frame brightness spike (the "whitening" shimmer) the denoiser cannot
	// hide under motion. Sweeping the (very few) directional lights every frame
	// makes the sun's contribution near-deterministic and low-variance. The two
	// strategies are combined by the balance heuristic (the `denom` term).
	//
	// References: Veach & Guibas 1995, "Optimally Combining Sampling Techniques
	// for Monte Carlo Rendering".
	for (uint d = 0; d < dirCount; ++d)
	{
		const uint lightIndex = dirIterator.first_item() + d;

		M += 1;

		const ShaderEntity light = load_entity(lightIndex);

		// Same combined density as every other strategy, so a light that more
		// than one strategy can draw is weighted consistently (unbiased).
		const float denom = RESTIRDIInitialDenom(
			light, candidateCount, RESTIR_DI_UNIFORM_CANDIDATES,
			lightCount, totalPower, tileBuffer >= 0);

		const float2 uv = float2(rng.next_float(), rng.next_float());
		RESTIRDIStreamInitialCandidate(
			light, lightIndex, uv, denom, P, N, weightSum, r, rng);
	}

	// Pre-multiply by M so the downstream unbiased weight W = weightSum / (M *
	// targetPdf) recovers the balance-heuristic W = weightSum / targetPdf
	// (balance-heuristic RIS does not divide by the sample count). Mirrors the
	// same trick in RESTIRDIMergeBalanceHeuristic.
	r.weightSum = M * weightSum;
	r.M = M;
	return r;
}

[numthreads(8, 8, 1)]
void main(uint3 DTid : SV_DispatchThreadID)
{
	if (DTid.x >= push.resolution.x || DTid.y >= push.resolution.y)
		return;

	// Stable light-index translation state (see restir_diHF.hlsli).
	RESTIRDILightMapBuffer = push.lightIndexMapBuffer;
	RESTIRDILightMapOffset = push.lightIndexMapOffset;
	RESTIRDISceneLightCount = push.sceneLightCount;

	const uint2 pixel = DTid.xy;
	const uint flatIndex = pixel.y * push.resolution.x + pixel.x;

	RESTIRDIReservoir reservoir = RESTIRDIReservoirInit();

	const float depth = texture_depth[pixel];
	[branch]
	if (depth > 0) // skip the sky
	{
		const float2 uv = (pixel + 0.5) * push.resolutionRcp;
		const float3 P = reconstruct_position(uv, depth);
		const float3 N = decode_normal(texture_normal_roughness[pixel]);

		RNG rng;
		rng.init(pixel, GetFrame().frame_count);

		// Shadowed RIS already traced a visibility ray per candidate, so the
		// selected sample is guaranteed visible (visibility = 1) and its weight
		// excludes occluded candidates - no separate initial-visibility test is
		// needed. The spatial pass still re-traces the final reused sample so
		// the shadow reflects the current frame.
		reservoir =
			RESTIRDISampleInitial(P, N, push.candidateCount, pixel, rng);
	}

	RESTIRDIReservoirStore(reservoirOutput, flatIndex, reservoir);
}
