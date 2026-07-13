#ifndef WI_RESTIR_GI_HF
#define WI_RESTIR_GI_HF
#include "ShaderInterop_ReSTIR.h"

/**
 * ReSTIR GI screen-space reservoir helpers.
 *
 * Operates on the GI reservoir (RESTIRGIReservoir), which holds an indirect
 * sample point (hit position + normal) and the radiance leaving it toward the
 * pixel. Provides weighted reservoir streaming, the diffuse target function,
 * the albedo-free irradiance resolve, and raw byte-address load/store (buffers
 * are ByteAddressBuffer so the same resource is read by both the GI passes and,
 * eventually, spatial reuse).
 *
 * Include order: "globals.hlsli" first (for PI, pack helpers), then this
 * header.
 */

/**
 * Diffuse target function for a GI sample evaluated at a surface (P, N).
 *
 * \[ \hat p = \mathrm{luminance}(L_o)\, \frac{\cos\theta_v}{\pi} \]
 *
 * where \( \cos\theta_v \) is the cosine at the receiver toward the sample
 * point. Albedo is factored out (applied later by the consumer), matching the
 * DI target and the forward+ integration.
 *
 * @param[in] sampleRadiance - Radiance leaving the sample point (L_o).
 * @param[in] samplePosition - World-space sample (hit) point.
 * @param[in] P - World-space receiver point.
 * @param[in] N - World-space receiver normal.
 *
 * @return The target function value p_hat (>= 0); 0 if degenerate/back-facing.
 */
inline float RESTIRGITargetFunction(
	float3 sampleRadiance,
	float3 samplePosition,
	float3 P,
	float3 N)
{
	const float3 d = samplePosition - P;
	const float dist2 = dot(d, d);
	if (dist2 <= 0)
		return 0;

	const float cosTheta = saturate(dot(N, d * rsqrt(dist2)));
	const float luma = dot(sampleRadiance, float3(0.2126, 0.7152, 0.0722));
	return luma * cosTheta / PI;
}

/**
 * Albedo-free indirect-diffuse irradiance contributed by a GI reservoir at
 * (P, N).
 *
 * \[ E = L_o\, \frac{\cos\theta_v}{\pi}\, W \]
 *
 * Forward shading multiplies this by the receiver albedo. For a single
 * cosine-weighted initial sample this reduces exactly to \( L_o \) (the
 * cosine and \(\pi\) cancel the source pdf), as expected.
 *
 * @param[in] r - The GI reservoir.
 * @param[in] P - World-space receiver point.
 * @param[in] N - World-space receiver normal.
 *
 * @return The albedo-free indirect irradiance (RGB, >= 0).
 */
inline float3 RESTIRGIResolve(RESTIRGIReservoir r, float3 P, float3 N)
{
	const float W = RESTIRGIReservoirGetInvPdf(r);
	if (!(W > 0))
		return float3(0, 0, 0);

	const float3 d = r.samplePosition - P;
	const float dist2 = dot(d, d);
	if (dist2 <= 0)
		return float3(0, 0, 0);

	const float cosTheta = saturate(dot(N, d * rsqrt(dist2)));
	return r.sampleRadiance * (cosTheta / PI) * W;
}

/**
 * Streams one indirect sample into a GI reservoir (weighted reservoir
 * sampling).
 *
 * @param[in,out] r - Reservoir to update.
 * @param[in] samplePosition - World-space hit point of the candidate.
 * @param[in] sampleNormal - World-space normal at the hit point.
 * @param[in] sampleRadiance - Radiance leaving the hit toward the pixel (L_o).
 * @param[in] targetPdf - Target function value p_hat of the candidate.
 * @param[in] risWeight - Resampling weight p_hat / p_source of the candidate.
 * @param[in,out] rng - Random generator.
 *
 * @return true if the candidate became the reservoir's selected sample.
 */
inline bool RESTIRGIReservoirUpdate(
	inout RESTIRGIReservoir r,
	float3 samplePosition,
	float3 sampleNormal,
	float3 sampleRadiance,
	float targetPdf,
	float risWeight,
	inout RNG rng)
{
	r.weightSum += risWeight;
	r.M += 1;

	if (r.weightSum > 0 && rng.next_float() * r.weightSum < risWeight)
	{
		r.samplePosition = samplePosition;
		r.sampleNormal = sampleNormal;
		r.sampleRadiance = sampleRadiance;
		r.targetPdf = targetPdf;
		return true;
	}
	return false;
}

/**
 * Reconnection-shift Jacobian for reusing a sample point across visible points.
 *
 * A sample point x_s was drawn from a source visible point (by hemisphere,
 * i.e. solid-angle, sampling). Reusing it at a different target visible point
 * changes the solid angle x_s subtends, so its measure must be corrected:
 *
 * \[
 * J = \frac{|\cos\theta_r|\; \lVert x_q - x_s \rVert^2}
 *          {|\cos\theta_q|\; \lVert x_r - x_s \rVert^2}
 * \]
 *
 * where \(\theta_q,\theta_r\) are the angles at the sample-point normal toward
 * the source and target visible points. Clamped to a bounded range so a grazing
 * source connection (tiny \(\cos\theta_q\)) cannot spike the reused weight into
 * a firefly.
 *
 * References: Ouyang et al. 2021, "ReSTIR GI", eq. (11).
 *
 * @param[in] samplePosition - World-space sample (hit) point x_s.
 * @param[in] sampleNormal - World-space normal at x_s.
 * @param[in] sourceP - Visible point the sample was drawn from.
 * @param[in] targetP - Visible point the sample is reused at.
 *
 * @return The (clamped) reconnection Jacobian (>= 0).
 */
inline float RESTIRGIReconnectionJacobian(
	float3 samplePosition,
	float3 sampleNormal,
	float3 sourceP,
	float3 targetP)
{
	const float3 dq = sourceP - samplePosition;
	const float3 dr = targetP - samplePosition;
	const float dq2 = dot(dq, dq);
	const float dr2 = dot(dr, dr);
	if (dq2 <= 0 || dr2 <= 0)
		return 0;

	const float cosQ = abs(dot(sampleNormal, dq * rsqrt(dq2)));
	const float cosR = abs(dot(sampleNormal, dr * rsqrt(dr2)));
	if (cosQ <= 1e-4)
		return 0;

	const float j = (cosR * dq2) / (cosQ * dr2);
	return clamp(j, 0.0, RESTIR_GI_JACOBIAN_CLAMP);
}

/**
 * Merges reservoir 'other' into 'self' with the reconnection Jacobian
 * (spatial/temporal reuse).
 *
 * The other reservoir's sample is re-weighted by its target function evaluated
 * at self's surface and by the shift-map Jacobian, keeping the merge unbiased
 * across differing visible points.
 *
 * @param[in,out] self - Reservoir receiving the merge (M and weightSum grow).
 * @param[in] other - Reservoir to merge in.
 * @param[in] targetPdfAtSelf - other's sample target function evaluated at
 *   self's surface (see RESTIRGITargetFunction).
 * @param[in] jacobian - Reconnection Jacobian (1 for temporal reuse, where the
 *   visible point is tracked across frames).
 * @param[in,out] rng - Random generator.
 */
inline void RESTIRGIReservoirMerge(
	inout RESTIRGIReservoir self,
	RESTIRGIReservoir other,
	float targetPdfAtSelf,
	float jacobian,
	inout RNG rng)
{
	const float risWeight = targetPdfAtSelf
		* RESTIRGIReservoirGetInvPdf(other) * other.M * jacobian;

	self.M += other.M;
	self.weightSum += risWeight;

	if (self.weightSum > 0 && rng.next_float() * self.weightSum < risWeight)
	{
		self.samplePosition = other.samplePosition;
		self.sampleNormal = other.sampleNormal;
		self.sampleRadiance = other.sampleRadiance;
		self.targetPdf = targetPdfAtSelf;
	}
}

/**
 * Loads a GI reservoir from a raw reservoir buffer. Each entry is 48 bytes.
 *
 * @param[in] buf - Reservoir byte-address buffer.
 * @param[in] index - Flat pixel index (y * width + x).
 *
 * @return The unpacked reservoir at that index.
 */
inline RESTIRGIReservoir RESTIRGIReservoirLoad(ByteAddressBuffer buf, uint index)
{
	RESTIRGIReservoirPacked p;
	p.data0 = buf.Load4(index * 48);
	p.data1 = buf.Load4(index * 48 + 16);
	p.data2 = buf.Load4(index * 48 + 32);
	return p.load();
}

/**
 * Loads a GI reservoir from a writable reservoir buffer (in-place read).
 *
 * @param[in] buf - Reservoir byte-address buffer (UAV).
 * @param[in] index - Flat pixel index (y * width + x).
 *
 * @return The unpacked reservoir at that index.
 */
inline RESTIRGIReservoir RESTIRGIReservoirLoad(
	RWByteAddressBuffer buf, uint index)
{
	RESTIRGIReservoirPacked p;
	p.data0 = buf.Load4(index * 48);
	p.data1 = buf.Load4(index * 48 + 16);
	p.data2 = buf.Load4(index * 48 + 32);
	return p.load();
}

/**
 * Stores a GI reservoir into a raw reservoir buffer.
 *
 * @param[in,out] buf - Reservoir byte-address buffer (UAV).
 * @param[in] index - Flat pixel index (y * width + x).
 * @param[in] r - The reservoir to store.
 */
inline void RESTIRGIReservoirStore(
	RWByteAddressBuffer buf, uint index, RESTIRGIReservoir r)
{
	RESTIRGIReservoirPacked p;
	p.store(r);
	buf.Store4(index * 48, p.data0);
	buf.Store4(index * 48 + 16, p.data1);
	buf.Store4(index * 48 + 32, p.data2);
}

#endif // WI_RESTIR_GI_HF
