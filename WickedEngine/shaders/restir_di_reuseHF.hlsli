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
 * This replaces the biased `÷M` combiner, which added every neighbor's full `M`
 * to the denominator even for samples the neighbor's surface could not
 * illuminate - underestimating W and darkening lighting/shadow discontinuities.
 * Because RIS still uses the unshadowed target, the estimator stays unbiased in
 * expectation for visibility (occluded samples stay eligible and contribute the
 * black frames that balance the bright ones); that per-frame visibility swing
 *is removed by the full-signal diffuse denoiser. See [[wickedengine-restir]].
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
 * @param[in] maxW - Firefly clamp for each source's unbiased weight (the light
 *   count: uniform RIS cannot legitimately exceed it).
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
	float maxW,
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

		mSum += sources[s].M;

		// Load the chosen light once, then re-resolve its target at each
		// surface (attenuation/cone/transmittance recomputed per surface -
		// exact, no baked-radiance approximation).
		const ShaderEntity light = load_entity(sources[s].lightIndex);

		// Target at the sample's own surface: the pdf that produced this
		// reservoir's weight, so W = weightSum / (M * pSelf) is its unbiased
		// contribution weight.
		const float pSelf =
			RESTIRDIReuseTarget(light, sources[s].uv, sources[s].P, sources[s].N);

		[branch]
		if (pSelf <= 0)
			continue;

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

#endif // WI_RESTIR_DI_REUSE_HF
