#include "globals.hlsli"
#include "lightingHF.hlsli"
#include "restir_lightsamplingHF.hlsli"

/**
 * ReSTIR - light tile pre-sampling pass.
 *
 * Fills the pre-sampled light-tile buffer once per frame (RTXDI's light
 * pre-sampling). Each thread owns one tile slot and resamples
 * RESTIR_PRESAMPLE_CANDIDATES uniformly-drawn analytic lights into a single
 * winner, importance-sampled by emitted power. The winner is stored together
 * with the reciprocal of its (power-proportional) selection probability, so the
 * per-pixel RIS in the DI initial / GI bounce passes can draw its candidates
 * from a small cached tile - coherent memory access and fewer wasted candidates
 * - while staying unbiased.
 *
 * The stored reciprocal source pdf is the reservoir's unbiased contribution
 * weight \[ W = \frac{1}{M}\sum_i \frac{p_{power}(y_i)}{p_{uniform}} \big/
 * p_{power}(Y) \]. Feeding samples whose reciprocal source pdf is this unbiased
 * weight into a second RIS layer is the standard RIS-of-RIS construction and
 * leaves the final estimator unbiased for any positive power proxy.
 *
 * References: Wyman & Panteleev 2021, "Rearchitecting Spatiotemporal Resampling
 *   for Production" (RTXDI light pre-sampling).
 */

PUSHCONSTANT(push, RESTIRPresamplePushConstants);

// Raw storage: RESTIRLightRef is 8 bytes (uint lightIndex, float invSourcePdf).
RWByteAddressBuffer lightTiles : register(u0);

[numthreads(64, 1, 1)]
void main(uint3 DTid : SV_DispatchThreadID)
{
	const uint slot = DTid.x;
	if (slot >= RESTIR_LIGHT_TILE_TOTAL)
		return;

	// Default to an empty slot (no light reaches this entry).
	uint winnerIndex = RESTIR_INVALID_LIGHT_INDEX;
	float invSourcePdf = 0;

	const ShaderEntityIterator iterator = lights();
	const uint lightCount = iterator.item_count();

	[branch]
	if (lightCount > 0)
	{
		RNG rng;
		rng.init(uint2(slot, 0), push.frameIndex);

		// Uniform source pmf for the mini-RIS candidates.
		const float invUniformPdf = (float)lightCount;

		// Weighted reservoir sampling with the light power as the RIS target.
		float weightSum = 0;
		float winnerPower = 0;
		for (uint i = 0; i < RESTIR_PRESAMPLE_CANDIDATES; ++i)
		{
			const uint lightIndex =
				iterator.first_item() + rng.next_uint(lightCount);
			const ShaderEntity light = load_entity(lightIndex);
			const float power = RESTIRLightPower(light);

			// risWeight = target / source = power / (1 / lightCount).
			const float risWeight = power * invUniformPdf;
			weightSum += risWeight;

			if (weightSum > 0 && rng.next_float() * weightSum < risWeight)
			{
				winnerIndex = lightIndex;
				winnerPower = power;
			}
		}

		// Unbiased contribution weight W = weightSum / (M * targetPdf), stored
		// as the reciprocal source pdf of this pre-sampled candidate.
		const float denom = RESTIR_PRESAMPLE_CANDIDATES * winnerPower;
		invSourcePdf = denom > 0 ? weightSum / denom : 0;
	}

	lightTiles.Store2(slot * 8, uint2(winnerIndex, asuint(invSourcePdf)));
}
