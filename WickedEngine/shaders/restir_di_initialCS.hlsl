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

	// ReGIR: when the pixel lies inside the camera-centered grid, its
	// candidates come from the cell built for this location (lights
	// importance-sampled by their contribution HERE, not by a global power
	// proxy) - this is what stops a dim light near the surface from being
	// crowded out by a brighter light far away. Outside the grid it falls back
	// to the global tiles / uniform.
	const int regirBuffer = push.regirBuffer;

	float weightSum = 0;
	float M = 0;

	for (uint i = 0; i < candidateCount; ++i)
	{
		// Draw one candidate light + its reciprocal source pdf.
		uint lightIndex;
		float invSourcePdf;
		bool drawn = false;

		[branch]
		if (regirBuffer >= 0)
		{
			// Stochastic trilinear cell pick blends the eight surrounding
			// cells' light sets so the grid does not print through as
			// low-frequency blobs on nearby surfaces, while staying weighted
			// toward the pixel's own cell (position-aware).
			const int cell = RESTIRReGIRSampleCellTrilinear(
				P, push.regirGridMin, push.regirCellSize, rng);

			[branch]
			if (cell >= 0)
			{
				const uint slot = (uint)cell * RESTIR_REGIR_LIGHTS_PER_CELL +
					rng.next_uint(RESTIR_REGIR_LIGHTS_PER_CELL);
				const RESTIRLightRef ref =
					RESTIRReGIRLoadSlot(regirBuffer, slot);

				// Only consume the ReGIR draw when the cell slot actually holds
				// a light. A cell that has not converged yet (or whose slot
				// caches no in-range light) returns an invalid reference with a
				// zero source pdf; treating that as a drawn candidate would
				// suppress the tile / uniform fallback below and leave the
				// pixel with no candidate at all - i.e. an unlit region
				// wherever the grid is not yet populated. Falling through keeps
				// coverage correct.
				[branch]
				if (ref.lightIndex != RESTIR_INVALID_LIGHT_INDEX &&
					ref.invSourcePdf > 0)
				{
					lightIndex = ref.lightIndex;
					invSourcePdf = ref.invSourcePdf;
					drawn = true;
				}
			}
		}

		[branch]
		if (!drawn && tileBuffer >= 0)
		{
			const uint slot = rng.next_uint(RESTIR_LIGHT_TILE_SIZE);
			const RESTIRLightRef ref =
				RESTIRLoadLightRef(tileBuffer, tileBase + slot);
			lightIndex = ref.lightIndex;
			invSourcePdf = ref.invSourcePdf;
			drawn = true;
		}

		[branch]
		if (!drawn)
		{
			lightIndex = iterator.first_item() + rng.next_uint(lightCount);
			invSourcePdf = (float)lightCount;
		}

		M += 1;

		// Only analytic lights are sampled here (skip empty / emissive-tri
		// refs).
		[branch]
		if (lightIndex == RESTIR_INVALID_LIGHT_INDEX ||
			(lightIndex & RESTIR_LIGHT_FLAG_EMISSIVE_TRIANGLE) != 0)
			continue;

		const ShaderEntity light = load_entity(lightIndex);
		const float2 uv = float2(rng.next_float(), rng.next_float());
		const RESTIRLightSample s = RESTIRResolveAnalyticLight(light, P, N, uv);
		const float unshadowedTarget = RESTIRTargetFunction(s, N);
		[branch]
		if (unshadowedTarget <= 0)
			continue;

		// One shadow ray per candidate: occluded -> zero target -> excluded.
		const float sampleDist =
			(s.distance >= FLT_MAX * 0.5) ? 100000.0 : s.distance;
		const float3 samplePos = P + s.direction * sampleDist;
		const float vis = RESTIRDITraceVisibility(P, N, samplePos, rng);
		[branch]
		if (vis <= 0)
			continue;

		const float risWeight = unshadowedTarget * vis * invSourcePdf;
		weightSum += risWeight;

		if (weightSum > 0 && rng.next_float() * weightSum < risWeight)
		{
			r.samplePosition = samplePos;
			r.sampleRadiance = s.radiance;
			// Visible winner: shadowed target == unshadowed target.
			r.targetPdf = unshadowedTarget;
			r.lightIndex = lightIndex;
			r.uv = uv;
			r.visibility = vis;
		}
	}

	r.weightSum = weightSum;
	r.M = M;
	return r;
}

[numthreads(8, 8, 1)]
void main(uint3 DTid : SV_DispatchThreadID)
{
	if (DTid.x >= push.resolution.x || DTid.y >= push.resolution.y)
		return;

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
