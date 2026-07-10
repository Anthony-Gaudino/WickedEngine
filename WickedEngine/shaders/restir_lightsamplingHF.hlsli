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
 * Resolves an analytic ShaderEntity light into a point sample at point P.
 *
 * Area extent (sphere radius, tube length, rectangle) is stochastically sampled
 * with rng, mirroring the sampling used by surfel_raytraceCS so results match.
 * The returned radiance includes distance/spot attenuation but NOT visibility.
 * solidAnglePdf is left at 1: the light-selection probability is carried by the
 * RIS source pdf in the driver, not here.
 *
 * @param[in] light - The analytic light entity.
 * @param[in] P - World-space shading point.
 * @param[in] N - World-space shading normal (unused here; kept for symmetry).
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
		L += sample_hemisphere_cos(L, rng) * light.GetRadius();
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
		pos += light.GetDirection() * (rng.next_float() - 0.5) * light.GetLength();
		pos += sample_hemisphere_cos(normalize(pos - P), rng) * light.GetRadius();
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
		pos += sample_hemisphere_cos(normalize(light.position - P), rng) * light.GetRadius();
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
			pos += right * (rng.next_float() - 0.5) * light.GetLength();
			pos += up * (rng.next_float() - 0.5) * light.GetHeight();
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

#endif // WI_RESTIR_LIGHTSAMPLING_HF
