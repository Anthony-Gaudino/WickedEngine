#ifndef WI_RESTIR_DI_HF
#define WI_RESTIR_DI_HF
#include "ShaderInterop_ReSTIR.h"

/**
 * ReSTIR DI screen-space reservoir helpers.
 *
 * Operates on the resolved-sample DI reservoir (RESTIRDIReservoir): weighted
 * reservoir streaming, spatial/temporal merge, and the diffuse target function
 * evaluated from a stored sample. Deliberately free of any light-sampling
 * dependency so the temporal/spatial passes (which never sample lights) do not
 * need lightingHF.hlsli. The initial pass, which does resample the light list,
 * includes restir_lightsamplingHF.hlsli itself.
 *
 * Include order: "globals.hlsli" first, then this header.
 */

/**
 * Diffuse target function for a stored DI sample at a surface.
 *
 * Computes the direction and distance from P to the sample point and returns
 * p_hat = luminance(sampleRadiance) * NdotL. Albedo is factored out (applied
 * later by the consumer's BRDF), matching the light-sampling core.
 *
 * @param[in] samplePosition - World-space point on the light.
 * @param[in] sampleRadiance - Unshadowed incident radiance of the sample.
 * @param[in] P - World-space shading point.
 * @param[in] N - World-space shading normal.
 * @param[out] dirToLight - Normalized direction from P toward the sample.
 * @param[out] distToLight - Distance from P to the sample.
 *
 * @return The target function value p_hat (>= 0); 0 if the sample is
 * degenerate.
 */
inline float RESTIRDITargetFunction(
	float3 samplePosition,
	float3 sampleRadiance,
	float3 P,
	float3 N,
	out float3 dirToLight,
	out float distToLight)
{
	dirToLight = float3(0, 0, 1);
	distToLight = 0;

	const float3 d = samplePosition - P;
	const float dist2 = dot(d, d);
	if (dist2 <= 0)
		return 0;

	distToLight = sqrt(dist2);
	dirToLight = d / distToLight;

	const float NdotL = saturate(dot(dirToLight, N));
	const float luma = dot(sampleRadiance, float3(0.2126, 0.7152, 0.0722));
	return luma * NdotL;
}

/**
 * Streams one resolved candidate into a DI reservoir (weighted reservoir
 * sampling).
 *
 * @param[in,out] r - Reservoir to update.
 * @param[in] samplePosition - World-space point on the candidate light.
 * @param[in] sampleRadiance - Unshadowed incident radiance of the candidate.
 * @param[in] lightIndex - Index of the candidate light.
 * @param[in] uv - Sample parameterization on the light (for re-resolution).
 * @param[in] targetPdf - Target function value p_hat of the candidate.
 * @param[in] risWeight - Resampling weight p_hat / p_source of the candidate.
 * @param[in,out] rng - Random generator.
 *
 * @return true if the candidate became the reservoir's selected sample.
 */
inline bool RESTIRDIReservoirUpdate(
	inout RESTIRDIReservoir r,
	float3 samplePosition,
	float3 sampleRadiance,
	uint lightIndex,
	float2 uv,
	float targetPdf,
	float risWeight,
	inout RNG rng)
{
	r.weightSum += risWeight;
	r.M += 1;

	if (r.weightSum > 0 && rng.next_float() * r.weightSum < risWeight)
	{
		r.samplePosition = samplePosition;
		r.sampleRadiance = sampleRadiance;
		r.lightIndex = lightIndex;
		r.uv = uv;
		r.targetPdf = targetPdf;
		return true;
	}
	return false;
}

/**
 * Merges DI reservoir 'other' into 'self' (spatial/temporal reuse).
 *
 * The merged sample carries over its cached visibility, and is re-weighted by
 * its target function evaluated at self's surface.
 *
 * Standard (unbiased) reuse: other's full confidence M is combined, so occluded
 * candidates stay eligible for selection. This is deliberate. RIS uses the
 * unshadowed target function, so it is only unbiased in expectation: a region
 * lit by one light through another light's shadow is correct only because the
 * occluded light is sometimes selected and contributes the black frames that
 * balance the over-bright frames (whose W folds in every candidate's unshadowed
 * target). Down-weighting the occluded reused samples by their visibility
 * (which an earlier revision did, to hide the per-frame flicker) removes those
 * black frames and leaves a net over-estimate - the shadow reads brighter than
 * the fully-shadowed surface, which is wrong. The per-frame variance this
 * reintroduces is instead removed by denoising the full diffuse contribution
 * (radiance x W x visibility x NdotL), where the bright and black frames
 * average to the correct colour. See [[wickedengine-restir]].
 *
 * @param[in,out] self - Reservoir receiving the merge.
 * @param[in] other - Reservoir to merge in.
 * @param[in] targetPdfAtSelf - other's sample target function at self's
 * surface.
 * @param[in] maxW - Firefly clamp for other's unbiased weight (typically the
 *   light count: uniform RIS cannot legitimately exceed it, so this bounds a
 *   near-zero-target sample from injecting spurious energy through reuse).
 * @param[in,out] rng - Random generator.
 */
inline void RESTIRDIReservoirMerge(
	inout RESTIRDIReservoir self,
	RESTIRDIReservoir other,
	float targetPdfAtSelf,
	float maxW,
	inout RNG rng)
{
	const float otherW = min(RESTIRDIReservoirGetInvPdf(other), maxW);
	const float risWeight = targetPdfAtSelf * otherW * other.M;

	self.M += other.M;
	self.weightSum += risWeight;

	if (self.weightSum > 0 && rng.next_float() * self.weightSum < risWeight)
	{
		self.samplePosition = other.samplePosition;
		self.sampleRadiance = other.sampleRadiance;
		self.lightIndex = other.lightIndex;
		self.uv = other.uv;
		self.targetPdf = targetPdfAtSelf;
		self.visibility = other.visibility;
	}
}

/**
 * Loads a DI reservoir from a raw reservoir buffer.
 *
 * The buffers are raw (ByteAddressBuffer) so the same resource can be read both
 * by the DI passes and bindlessly by forward shading. Each entry is 48 bytes.
 *
 * @param[in] buf - Reservoir byte-address buffer.
 * @param[in] index - Flat pixel index (y * width + x).
 *
 * @return The unpacked reservoir at that index.
 */
inline RESTIRDIReservoir RESTIRDIReservoirLoad(ByteAddressBuffer buf, uint index)
{
	RESTIRDIReservoirPacked p;
	p.data0 = buf.Load4(index * 48);
	p.data1 = buf.Load4(index * 48 + 16);
	p.data2 = buf.Load4(index * 48 + 32);
	return p.load();
}

/**
 * Loads a DI reservoir from a writable reservoir buffer (in-place read).
 *
 * @param[in] buf - Reservoir byte-address buffer (UAV).
 * @param[in] index - Flat pixel index (y * width + x).
 *
 * @return The unpacked reservoir at that index.
 */
inline RESTIRDIReservoir RESTIRDIReservoirLoad(RWByteAddressBuffer buf, uint index)
{
	RESTIRDIReservoirPacked p;
	p.data0 = buf.Load4(index * 48);
	p.data1 = buf.Load4(index * 48 + 16);
	p.data2 = buf.Load4(index * 48 + 32);
	return p.load();
}

/**
 * Stores a DI reservoir into a raw reservoir buffer.
 *
 * @param[in,out] buf - Reservoir byte-address buffer (UAV).
 * @param[in] index - Flat pixel index (y * width + x).
 * @param[in] r - The reservoir to store.
 */
inline void RESTIRDIReservoirStore(
	RWByteAddressBuffer buf, uint index, RESTIRDIReservoir r)
{
	RESTIRDIReservoirPacked p;
	p.store(r);
	buf.Store4(index * 48, p.data0);
	buf.Store4(index * 48 + 16, p.data1);
	buf.Store4(index * 48 + 32, p.data2);
}

#endif // WI_RESTIR_DI_HF
