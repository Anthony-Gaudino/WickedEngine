#ifndef WI_RESTIR_DI_REUSE_HF
#define WI_RESTIR_DI_REUSE_HF

/**
 * ReSTIR DI spatiotemporal reuse: balance-heuristic MIS merge.
 *
 * Split out from restir_diHF.hlsli because it **re-resolves the analytic light
 * per surface** (needs `load_entity` / `RESTIRResolveAnalyticLight`), which the
 * light-free consumers of restir_diHF.hlsli (the denoise passes, forward
 * shading) must not pull in. Only the temporal and spatial passes include this.
 *
 * Include order: "globals.hlsli", "lightingHF.hlsli",
 * "restir_lightsamplingHF.hlsli" and "restir_diHF.hlsli" first, then this
 * header.
 */

/**
 * Diffuse target of a reservoir's chosen light, re-resolved against the live
 * light at an arbitrary surface.
 *
 * The reservoir stores a light **reference** (index + area sample `uv`), so the
 * exact incident radiance - with the light's real distance attenuation, spot
 * cone and directional transmittance - is recomputed at whatever surface the
 * merge needs it, rather than reusing the radiance that was baked for the
 * sample's originating surface. That baked radiance is only correct at its own
 * distance; the balance-heuristic MIS weights compare a sample across surfaces
 * at different distances from the light, so an approximation there biases the
 * estimate (most visibly, a dark ring around a point light). Re-resolving is
 * exact and matches how the spatial pass re-resolves the final chosen sample.
 *
 * @param[in] light - The chosen analytic light entity.
 * @param[in] uv - Stored area-sample parameterization on the light.
 * @param[in] P - World-space shading point to evaluate the target at.
 * @param[in] N - World-space shading normal at that point.
 *
 * @return The diffuse target function value p_hat (>= 0): luminance of the
 *   attenuated incident radiance times N·L. Albedo is factored out.
 */
inline float RESTIRDIReuseTarget(
	ShaderEntity light, float2 uv, float3 P, float3 N)
{
	const RESTIRLightSample s = RESTIRResolveAnalyticLight(light, P, N, uv);
	const float NdotL = saturate(dot(s.direction, N));
	const float luma = dot(s.radiance, float3(0.2126, 0.7152, 0.0722));
	return luma * NdotL;
}

/**
 * Combines several reservoirs at a surface using the **generalized balance
 * heuristic** MIS (unbiased spatiotemporal reuse).
 *
 * Each source `s` streams its sample into the result with resampling weight
 * \[
 * w_s = m_s(X_s) \; \hat p_c(X_s) \; W_s, \qquad
 * m_s(X_s) = \frac{M_s \, \hat p_s(X_s)}{\sum_j M_j \, \hat p_j(X_s)},
 * \]
 * where \f$\hat p_c\f$ is the target at the center (canonical) surface,
 * \f$\hat p_j(X_s)\f$ is the sample's target re-resolved at source `j`'s
   **own**
 * surface, \f$W_s\f$ is source `s`'s (firefly-clamped) unbiased contribution
 * weight and \f$M_j\f$ are the confidences. The MIS denominator counts a
 * source's confidence only where that source could actually have produced the
 * sample (its own-surface target is non-zero); the final weight is formed from
 * this MIS-normalized sum, **not** by dividing by the total sample count.
 *
 * The \f$O(N^2)\f$ cost (each source re-resolved at every other source's
 * surface) is fine at the small spatial sample counts used here, and it is the
 * default reuse combiner: it reproduces correct soft shadows where the linear
 * pairwise variant (RESTIRDIMergePairwiseMIS) does not.
 *
 * The result is packed so RESTIRDIReservoirGetInvPdf still returns the correct
 * weight: `weightSum` is stored as (Σ w_s)·M so the consumer's `÷M` cancels,
 * leaving W = Σ w_s / targetPdf.
 *
 * References:
 * Talbot 2005, "Importance Resampling for Global Illumination". Lin et al.
 * 2022, "Generalized Resampled Importance Sampling". Wyman 2023 SIGGRAPH
 *course, "A Gentle Introduction to ReSTIR".
 *
 * @param[in] sources - Input reservoirs; index 0 must be the canonical (center)
 *   reservoir. Entries with M <= 0 or a non-analytic light ref are skipped.
 * @param[in] count - Number of valid entries (1 .. RESTIR_DI_MIS_MAX_SOURCES).
 * @param[in] centerP - Center (canonical) surface shading point.
 * @param[in] centerN - Center (canonical) surface shading normal.
 * @param[in,out] rng - Random generator.
 *
 * @return The merged reservoir (visibility of the selected sample carried
 *         over).
 */
inline RESTIRDIReservoir RESTIRDIMergeBalanceHeuristic(
	RESTIRDIMISSource sources[RESTIR_DI_MIS_MAX_SOURCES],
	uint count,
	float3 centerP,
	float3 centerN,
	inout RNG rng)
{
	RESTIRDIReservoir result = RESTIRDIReservoirInit();

	float wSum = 0;
	float mSum = 0;

	for (uint s = 0; s < count; ++s)
	{
		[branch]
		if (sources[s].M <= 0)
			continue;

		// DI resolves analytic lights only; skip empty or emissive-triangle
		// refs.
		[branch]
		if (sources[s].lightIndex == RESTIR_INVALID_LIGHT_INDEX ||
			(sources[s].lightIndex & RESTIR_LIGHT_FLAG_EMISSIVE_TRIANGLE) != 0)
			continue;

		// The reservoir stores the light's STABLE scene index; map it to this
		// frame's entity slot. A light culled this frame no longer reaches any
		// visible surface, so drop the source entirely (before it counts toward
		// the confidence sum).
		const uint slot = RESTIRDICurrentLightSlot(sources[s].lightIndex);
		[branch]
		if (slot == RESTIR_INVALID_LIGHT_INDEX)
			continue;

		mSum += sources[s].M;

		// Load the chosen light once, then re-resolve its target at each
		// surface (attenuation/cone/transmittance recomputed per surface -
		// exact, no baked-radiance approximation).
		const ShaderEntity light = load_entity(slot);

		// Target at the sample's own surface: the pdf that produced this
		// reservoir's weight, so W = weightSum / (M * pSelf) is its unbiased
		// contribution weight.
		const float pSelf =
			RESTIRDIReuseTarget(light, sources[s].uv, sources[s].P, sources[s].N);

		[branch]
		if (pSelf <= 0)
			continue;

		// Per-light firefly clamp. Under power-proportional light sampling a
		// light's unbiased weight legitimately reaches totalLightPower / its
		// own power (a dim light is sampled proportionally less, so each hit
		// carries proportionally more); clamping at the light count (a
		// uniform-sampling bound) would lop that off and darken dim lights.
		// This bound instead caps each light's contribution at
		// ~totalLightPower, so dim lights stay correct while true fireflies are
		// still bounded.
		const float maxW = GetFrame().totalLightPower / RESTIRLightPower(light);
		const float wS = min(sources[s].weightSum / (sources[s].M * pSelf), maxW);

		// Target of this source's sample at the center (canonical) surface.
		const float pCenter =
			RESTIRDIReuseTarget(light, sources[s].uv, centerP, centerN);

		// Balance-heuristic denominator: confidence times own-surface target,
		// summed over all sources. A source whose surface cannot see this
		// sample's light adds nothing, so it no longer dilutes the weight.
		float denom = 0;
		for (uint j = 0; j < count; ++j)
		{
			const float pj = RESTIRDIReuseTarget(
				light, sources[s].uv, sources[j].P, sources[j].N);
			denom += sources[j].M * pj;
		}

		[branch]
		if (denom <= 0)
			continue;

		const float misWeight = (sources[s].M * pSelf) / denom;
		const float risWeight = misWeight * pCenter * wS;
		wSum += risWeight;

		if (wSum > 0 && rng.next_float() * wSum < risWeight)
		{
			result.samplePosition = sources[s].samplePosition;
			result.sampleRadiance = sources[s].sampleRadiance;
			result.lightIndex = sources[s].lightIndex;
			result.uv = sources[s].uv;
			result.targetPdf = pCenter;
			result.visibility = sources[s].visibility;
		}
	}

	result.M = mSum;
	// Pre-multiply by M so the consumer's RESTIRDIReservoirGetInvPdf (which
	// divides by M) cancels it, leaving the MIS-normalized W = wSum /
	// targetPdf.
	result.weightSum = (mSum > 0) ? wSum * mSum : 0;
	return result;
}

/**
 * Balance-heuristic MIS weight for one pairwise comparison.
 *
 * Weighs technique 0 (target `w0`, confidence `m0`) against technique 1 (target
 * `w1`, confidence `m1`): \f$\tfrac{m_0 w_0}{m_0 w_0 + m_1 w_1}\f$.
 *
 * @param[in] w0 - Target of the sample in technique 0's domain.
 * @param[in] w1 - Target of the sample in technique 1's domain.
 * @param[in] m0 - Confidence of technique 0.
 * @param[in] m1 - Confidence of technique 1.
 *
 * @return The balance-heuristic weight in [0, 1] (0 if the denominator is 0).
 */
inline float RESTIRDIPairwiseMisWeight(float w0, float w1, float m0, float m1)
{
	const float denom = m0 * w0 + m1 * w1;
	return denom > 0 ? max(0.0, m0 * w0) / denom : 0;
}

/**
 * Combines several reservoirs at a surface using **pairwise MIS** (unbiased
 * spatiotemporal reuse), faithfully mirroring RTXDI's PairwiseStreaming.
 *
 * Each non-canonical neighbor is paired with the canonical (center) reservoir
 * alone, so the merge costs \f$O(N)\f$ target evaluations rather than the
 * generalized balance heuristic's \f$O(N^2)\f$. In each pairwise step the
 * neighbor's confidence is scaled by `neighborBudget` (the intended neighbor
 * count, **not** how many turned out valid - so a neighbor is not
 * under-weighted when others are rejected, which otherwise over-favors the
 * canonical and hardens shadow penumbrae). The canonical, touched once per
 * neighbor, is streamed a single time afterwards with the accumulated
 * compensation weight  \f$\sum(1 - w_1)\f$, and the estimate is normalized by
 * the number of valid neighbors.
 *
 * The result is packed so RESTIRDIReservoirGetInvPdf returns the correct
 * weight: `weightSum` is stored as (Σ w)·M / N so the consumer's `÷M` leaves
 * W = Σ w / (N · targetPdf).
 *
 * Each source's contribution weight \f$W\f$ carries the same per-light firefly
 * clamp (totalLightPower / the light's own power). A source whose light was
 * culled this frame, or whose sample cannot illuminate the surface, contributes
 * nothing; the canonical always keeps a share, so a disocclusion falls back to
 * the pixel's own sample instead of darkening.
 *
 * References:
 * Bitterli 2022 (Chapter 9.1). Hedstrom et al. 2026, "Stochastic Pairwise MIS
 * for Unbiased Large-Kernel Reuse in Real Time". Wyman 2023 SIGGRAPH course,
 * "A Gentle Introduction to ReSTIR".
 *
 * @param[in] sources - Input reservoirs; index 0 must be the canonical (center)
 *   reservoir. Entries with M <= 0 or a non-analytic light ref are skipped.
 * @param[in] count - Number of valid entries (1 .. RESTIR_DI_MIS_MAX_SOURCES).
 * @param[in] neighborBudget - Intended neighbor count (the pairwise weight
 *   scales each neighbor's confidence by this).
 * @param[in] centerP - Center (canonical) surface shading point.
 * @param[in] centerN - Center (canonical) surface shading normal.
 * @param[in,out] rng - Random generator.
 *
 * @return The merged reservoir (visibility of the selected sample carried
 *         over).
 */
inline RESTIRDIReservoir RESTIRDIMergePairwiseMIS(
	RESTIRDIMISSource sources[RESTIR_DI_MIS_MAX_SOURCES],
	uint count,
	uint neighborBudget,
	float3 centerP,
	float3 centerN,
	inout RNG rng)
{
	RESTIRDIReservoir result = RESTIRDIReservoirInit();
	[branch]
	if (count == 0)
		return result;

	// Resolve the canonical (source 0): its light, its target at the center
	// surface, and its firefly-clamped contribution weight.
	bool cValid = sources[0].M > 0 &&
		sources[0].lightIndex != RESTIR_INVALID_LIGHT_INDEX &&
		(sources[0].lightIndex & RESTIR_LIGHT_FLAG_EMISSIVE_TRIANGLE) == 0;
	ShaderEntity cLight = (ShaderEntity)0;
	float cW = 0;
	float cPCenter = 0;
	float cc = 0;
	[branch]
	if (cValid)
	{
		const uint slot = RESTIRDICurrentLightSlot(sources[0].lightIndex);
		[branch]
		if (slot == RESTIR_INVALID_LIGHT_INDEX)
		{
			cValid = false;
		}
		else
		{
			cLight = load_entity(slot);
			cPCenter =
				RESTIRDIReuseTarget(cLight, sources[0].uv, centerP, centerN);
			[branch]
			if (cPCenter <= 0)
			{
				cValid = false;
			}
			else
			{
				const float maxW =
					GetFrame().totalLightPower / RESTIRLightPower(cLight);
				cW = min(
					sources[0].weightSum / (sources[0].M * cPCenter), maxW);
				cc = sources[0].M;
			}
		}
	}

	// Resolve each neighbor once: cache its own-surface and center-surface
	// targets and its clamped contribution weight, and count the valid
	// neighbors (the normalization denominator). A neighbor with a
	// culled/invalid light or a sample that cannot illuminate its own surface
	// is dropped.
	bool nValid[RESTIR_DI_MIS_MAX_SOURCES];
	float nW[RESTIR_DI_MIS_MAX_SOURCES];
	float nPOwn[RESTIR_DI_MIS_MAX_SOURCES];
	float nPCenter[RESTIR_DI_MIS_MAX_SOURCES];
	uint validNeighbors = 0;
	for (uint i = 1; i < count; ++i)
	{
		nValid[i] = false;
		[branch]
		if (sources[i].M <= 0 ||
			sources[i].lightIndex == RESTIR_INVALID_LIGHT_INDEX ||
			(sources[i].lightIndex & RESTIR_LIGHT_FLAG_EMISSIVE_TRIANGLE) != 0)
			continue;

		const uint slot = RESTIRDICurrentLightSlot(sources[i].lightIndex);
		[branch]
		if (slot == RESTIR_INVALID_LIGHT_INDEX)
			continue;

		const ShaderEntity light = load_entity(slot);
		const float pOwn = RESTIRDIReuseTarget(
			light, sources[i].uv, sources[i].P, sources[i].N);
		[branch]
		if (pOwn <= 0)
			continue;

		const float maxW =
			GetFrame().totalLightPower / RESTIRLightPower(light);
		nW[i] = min(sources[i].weightSum / (sources[i].M * pOwn), maxW);
		nPOwn[i] = pOwn;
		nPCenter[i] =
			RESTIRDIReuseTarget(light, sources[i].uv, centerP, centerN);
		nValid[i] = true;
		validNeighbors += 1;
	}

	// The pairwise weight scales the neighbor confidence by the intended
	// budget, clamped so it is never below the number actually streamed.
	const float budget = max((float)neighborBudget, (float)validNeighbors);

	float wSum = 0;
	float mAccum = 0;

	// Pairwise MIS touches the canonical once per neighbor (overweighting it);
	// each step records how much (1 - w1) so the single canonical stream after
	// the loop compensates exactly.
	float canonicalWeight = 0;

	for (uint i = 1; i < count; ++i)
	{
		[branch]
		if (!nValid[i])
			continue;

		const float m0 = sources[i].M * budget;
		const float m1 = cc;

		// Neighbor MIS weight: its own-surface target vs the canonical's; then
		// stream its center-surface contribution.
		const float w0 =
			RESTIRDIPairwiseMisWeight(nPOwn[i], nPCenter[i], m0, m1);
		const float risWeight = nPCenter[i] * nW[i] * w0;

		// Canonical's own sample re-resolved at THIS neighbor's surface, for
		// the canonical's pairwise compensation weight.
		const float cAtI = cValid
			? RESTIRDIReuseTarget(
				cLight, sources[0].uv, sources[i].P, sources[i].N)
			: 0;
		const float w1 = RESTIRDIPairwiseMisWeight(cAtI, cPCenter, m0, m1);
		canonicalWeight += 1.0 - w1;

		// Carry the full confidence (matching the balance heuristic's mSum), so
		// the reservoir's temporal weight next frame is not disrupted at shadow
		// boundaries.
		wSum += risWeight;
		mAccum += sources[i].M;
		if (wSum > 0 && rng.next_float() * wSum < risWeight)
		{
			result.samplePosition = sources[i].samplePosition;
			result.sampleRadiance = sources[i].sampleRadiance;
			result.lightIndex = sources[i].lightIndex;
			result.uv = sources[i].uv;
			result.targetPdf = nPCenter[i];
			result.visibility = sources[i].visibility;
		}
	}

	// With no usable neighbor the canonical carries the full weight.
	if (validNeighbors == 0)
		canonicalWeight = 1;

	// Stream the canonical once, compensating for the pairwise overweighting.
	[branch]
	if (cValid)
	{
		const float risWeight = cPCenter * cW * canonicalWeight;
		wSum += risWeight;
		mAccum += cc;
		if (wSum > 0 && rng.next_float() * wSum < risWeight)
		{
			result.samplePosition = sources[0].samplePosition;
			result.sampleRadiance = sources[0].sampleRadiance;
			result.lightIndex = sources[0].lightIndex;
			result.uv = sources[0].uv;
			result.targetPdf = cPCenter;
			result.visibility = sources[0].visibility;
		}
	}

	// Finalize: pairwise MIS normalizes by the valid-neighbor count. Store
	// weightSum pre-multiplied by M so the consumer's
	// RESTIRDIReservoirGetInvPdf (which divides by M) leaves W = wSum /
	// (validNeighbors * targetPdf).
	const float denom = max(1.0, (float)validNeighbors);
	result.M = mAccum;
	result.weightSum =
		(mAccum > 0 && result.targetPdf > 0) ? wSum * mAccum / denom : 0;
	return result;
}

#endif // WI_RESTIR_DI_REUSE_HF
