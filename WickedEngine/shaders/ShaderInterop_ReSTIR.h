#ifndef WI_SHADERINTEROP_RESTIR_H
#define WI_SHADERINTEROP_RESTIR_H
#include "ShaderInterop.h"
#include "ShaderInterop_Renderer.h"

/**
 * ReSTIR (Reservoir-based Spatiotemporal Importance Resampling) shared interop.
 *
 * This header defines the data shared between the CPU (light-buffer build,
 * resource allocation, dispatch) and the GPU (light sampling, resampling). It
 * is the foundation of the ReSTIR feature and is intentionally kept free of any
 * pass-specific state so the same reservoir/RIS core can be reused by:
 *   - ReSTIR DI (screen-space direct lighting),
 *   - ReSTIR GI (screen-space indirect lighting),
 *   - the path tracer, DDGI and SurfelGI ray-hit shading (light sampling only).
 *
 * References: Bitterli et al. 2020, "Spatiotemporal reservoir resampling for
 *   real-time ray tracing with dynamic direct lighting" (ReSTIR). Wyman &
 *   Panteleev 2021, "Rearchitecting Spatiotemporal Resampling for Production"
 *   (RTXDI, light presampling tiles).
 */

/**
 * Number of candidate lights pre-sampled into a single light tile.
 *
 * Per-pixel RIS draws its candidates from one randomly chosen tile instead of
 * the full light set, which keeps memory access coherent when the scene has
 * many lights.
 */
static const uint RESTIR_LIGHT_TILE_SIZE = 1024;

/** Number of independently pre-sampled light tiles kept per frame. */
static const uint RESTIR_LIGHT_TILE_COUNT = 128;

/** Total number of pre-sampled light references (one structured buffer). */
static const uint RESTIR_LIGHT_TILE_TOTAL =
	RESTIR_LIGHT_TILE_SIZE * RESTIR_LIGHT_TILE_COUNT;

/**
 * Confidence cap for temporal reuse.
 *
 * A reservoir's history weight (M) is clamped to this so stale samples cannot
 * dominate and the estimator stays responsive to lighting changes.
 */
static const uint RESTIR_MAX_HISTORY_LENGTH = 20;

/**
 * Marks a light reference that addresses an emissive mesh triangle (in the
 * emissive-triangle light buffer) rather than an analytic ShaderEntity light.
 */
static const uint RESTIR_LIGHT_FLAG_EMISSIVE_TRIANGLE = 0x80000000u;

/** Index value denoting "no light" (an empty reservoir / failed sample). */
static const uint RESTIR_INVALID_LIGHT_INDEX = 0x7FFFFFFFu;

/*
################################################################################
Emissive triangle light
################################################################################
*/

/**
 * One emissive mesh triangle promoted to an area light source.
 *
 * The CPU rebuilds this buffer each frame from meshes whose material is
 * emissive; the GPU samples it uniformly by area within a triangle and by flux
 * across triangles (via the pre-sampling pass). Positions are world space, and
 * each is padded to a 16-byte boundary for GPU alignment.
 */
struct RESTIREmissiveTriangle
{
	float3 p0;
	float pad0;
	float3 p1;
	float pad1;
	float3 p2;
	float pad2;

	/**
	 * Emitted radiance toward the front face (material emissive color scaled by
	 * intensity), constant over the triangle. RGB, linear.
	 */
	float3 radiance;

	/**
	 * Luminance-weighted flux (radiance luminance * world-space area). Used as
	 * the importance weight when building the alias/CDF sampling table.
	 */
	float flux;
};

/*
################################################################################
Pre-sampled light reference
################################################################################
*/

/**
 * One entry of a pre-sampled light tile.
 *
 * A light chosen with some source pdf, stored together with the reciprocal of
 * that pdf so per-pixel RIS can weight it without re-deriving the selection
 * probability.
 */
struct RESTIRLightRef
{
	/**
	 * Analytic ShaderEntity index, or emissive-triangle index OR'd with
	 * RESTIR_LIGHT_FLAG_EMISSIVE_TRIANGLE.
	 */
	uint lightIndex;

	/**
	 * 1 / p_source(lightIndex): reciprocal of the probability with which this
	 * light was selected into the tile.
	 */
	float invSourcePdf;
};

/*
################################################################################
Light sample
################################################################################
*/

#ifndef __cplusplus
/**
 * A resolved point sample on a light, ready for target-function and BRDF
 * evaluation.
 *
 * Produced by the light-sampling library from a RESTIRLightRef plus a surface
 * point. GPU-only (uses shader vector types directly).
 */
struct RESTIRLightSample
{
	/**
	 * World-space direction from the shading point toward the light
	 * (normalized).
	 */
	float3 direction;

	/**
	 * Distance from the shading point to the sample (FLT_MAX for directional).
	 */
	float distance;

	/**
	 * Incident radiance arriving from the light along -direction, already
	 * including analytic attenuation / spot falloff (but NOT visibility).
	 */
	float3 radiance;

	/**
	 * Probability density of this sample in solid-angle measure at the shading
	 * point. Zero marks an invalid (unsamplable) sample.
	 */
	float solidAnglePdf;
};
#endif // __cplusplus

/*
################################################################################
Reservoir
################################################################################
*/

#ifndef __cplusplus
/**
 * A working reservoir (kept in registers during a pass).
 *
 * Holds the single surviving light sample selected by weighted reservoir
 * sampling over all seen candidates, plus the statistics needed to form the
 * unbiased contribution weight and to merge with other reservoirs.
 */
struct RESTIRReservoir
{
	/** Chosen light reference (see RESTIRLightRef.lightIndex encoding). */
	uint lightIndex;

	/**
	 * Sample parameterization on the chosen light: barycentric (x,y) for an
	 * emissive triangle, or the canonical area/disk sample for an analytic
	 * light.
	 */
	float2 uv;

	/** Target function value p_hat of the currently held sample. */
	float targetPdf;

	/** Running sum of resampling weights of all candidates streamed in. */
	float weightSum;

	/**
	 * Confidence: effective number of candidates represented (float to allow
	 * clamping and MIS-weighted merges).
	 */
	float M;
};

/**
 * Resets a reservoir to empty.
 *
 * @return A reservoir holding no sample (lightIndex invalid, all stats zero).
 */
inline RESTIRReservoir RESTIRReservoirInit()
{
	RESTIRReservoir r;
	r.lightIndex = RESTIR_INVALID_LIGHT_INDEX;
	r.uv = float2(0, 0);
	r.targetPdf = 0;
	r.weightSum = 0;
	r.M = 0;
	return r;
}

/**
 * Unbiased contribution weight W = weightSum / (M * targetPdf).
 *
 * Guarded against division by zero (an empty or degenerate reservoir yields 0).
 * Multiplying the evaluated (visible) sample contribution by W gives the ReSTIR
 * estimate.
 *
 * @param[in] r - The reservoir to weigh.
 *
 * @return The unbiased contribution weight W (>= 0).
 */
inline float RESTIRReservoirGetInvPdf(RESTIRReservoir r)
{
	const float denom = r.M * r.targetPdf;
	return denom > 0 ? r.weightSum / denom : 0;
}
#endif // __cplusplus

/**
 * GPU-storable packed reservoir (32 bytes).
 *
 * Used by the screen-space DI/GI passes for their temporal history and spatial
 * reuse buffers; the light-only consumers (path tracer, DDGI, SurfelGI) never
 * store reservoirs and can ignore this. The layout keeps everything needed to
 * reconstruct the estimator across frames plus an optional cached visibility.
 */
struct RESTIRReservoirPacked
{
	uint4 data0;
	uint4 data1;

#ifndef __cplusplus
	/**
	 * Packs a reservoir and its cached visibility into this storage slot.
	 *
	 * @param[in] r - Reservoir to store.
	 * @param[in] cachedVisibility - Visibility of the selected sample (RGB),
	 *   packed into R11G11B10 for cheap reuse without re-tracing.
	 */
	inline void store(RESTIRReservoir r, float3 cachedVisibility)
	{
		data0.x = r.lightIndex;
		data0.y = pack_half2(r.uv.x, r.uv.y);
		data0.z = asuint(r.targetPdf);
		data0.w = asuint(r.weightSum);
		data1.x = asuint(r.M);
		data1.y = Pack_R11G11B10_FLOAT(cachedVisibility);
		data1.z = 0;
		data1.w = 0;
	}

	/**
	 * Unpacks the stored reservoir (without the cached visibility).
	 *
	 * @return The reconstructed reservoir.
	 */
	inline RESTIRReservoir load()
	{
		RESTIRReservoir r;
		r.lightIndex = data0.x;
		r.uv = unpack_half2(data0.y);
		r.targetPdf = asfloat(data0.z);
		r.weightSum = asfloat(data0.w);
		r.M = asfloat(data1.x);
		return r;
	}

	/**
	 * Unpacks the cached visibility of the stored sample.
	 *
	 * @return The cached visibility (RGB).
	 */
	inline float3 loadVisibility()
	{
		return Unpack_R11G11B10_FLOAT(data1.y);
	}
#endif // __cplusplus
};

/*
################################################################################
Screen-space DI reservoir
################################################################################
*/

/** ReSTIR DI RIS candidate count for the initial pass. */
static const uint RESTIR_DI_INITIAL_CANDIDATES = 8;

/** Spatial-reuse neighbor count per pixel. */
static const uint RESTIR_DI_SPATIAL_SAMPLES = 4;

/** Spatial-reuse search radius, in pixels. */
static const float RESTIR_DI_SPATIAL_RADIUS = 16.0;

/**
 * Visibility denoiser history cap. The temporal accumulation blends the fresh
 * per-frame visibility toward a running mean with weight 1 / historyLength;
 * capping the length floors that weight so the shadow stays responsive to
 * change (a higher cap denoises more but reacts slower).
 *
 * Kept high because the A-SVGF temporal gradient adapts the effective history
 * per pixel: static shadows accumulate the full length (very clean), while a
 * pixel whose shadow actually changed this frame resets its history to ~1
 * (ghost-free), so a long cap no longer means long ghosting.
 */
static const float RESTIR_DI_DENOISE_MAX_HISTORY = 32.0;

/**
 * Dilation radius (in pixels) of the A-SVGF temporal gradient before it drives
 * the history reset.
 *
 * The exact gradient is measured per pixel by re-tracing the previous frame's
 * sample under the current scene, so it fires only on the thin band of pixels
 * whose occlusion actually flipped. Dilating (a max over a
 * (2r+1)x(2r+1) window) spreads that reset across the moving edge's
 * neighborhood so the whole affected region refreshes together, instead of
 * leaving a one-pixel seam.
 */
static const int RESTIR_DI_DENOISE_GRADIENT_DILATE = 2;

/**
 * Sensitivity of the history reset to the (dilated) temporal gradient.
 *
 * The gradient is in [0, 1] (fraction of the neighborhood whose shadow
 * changed); the reset applied to the history is saturate(gradient *
 * sensitivity). Larger resets more eagerly on partial change (less ghost,
 * slightly more edge flicker); 1 resets only in proportion to the measured
 * change.
 */
static const float RESTIR_DI_DENOISE_GRADIENT_SENSITIVITY = 4.0;

/** Number of edge-aware a-trous passes in the spatial visibility denoiser. */
static const uint RESTIR_DI_DENOISE_ATROUS_PASSES = 4;

/**
 * Depth edge-stop scale for the spatial a-trous filter.
 *
 * The depth weight is \[ \exp(-|z_c - z_t| / (s\, z_c + \epsilon)) \], so a
 * smaller scale keeps the blur tighter across depth discontinuities.
 */
static const float RESTIR_DI_DENOISE_DEPTH_SCALE = 0.5;

/** Normal edge-stop exponent: higher preserves creases/silhouettes more. */
static const float RESTIR_DI_DENOISE_NORMAL_POWER = 64.0;

/**
 * Visibility (luminance) edge-stop scale for the spatial a-trous filter.
 *
 * Scales the standard deviation \[ \sqrt{\mathrm{Var}} \] used as the tolerance
 * on visibility differences: the filter blurs freely where the temporal
 * variance is high (noise) but keeps crisp penumbra gradients where it is low.
 */
static const float RESTIR_DI_DENOISE_LUMINANCE_SCALE = 4.0;

#ifndef __cplusplus
/**
 * A screen-space ReSTIR DI reservoir.
 *
 * Unlike the light-ref reservoir used by the GI ray-hit consumers, this one
 * carries the *resolved* light sample (a world-space point on the light plus
 * its unshadowed incident radiance) so temporal/spatial reuse can re-evaluate
 * the target function at a neighbor surface without re-sampling the light. The
 * carried radiance is the value at the sample's originating surface; reuse
 * currently ignores the attenuation/Jacobian change across pixels, which is a
 * small-radius approximation to be tightened later (proper light-ref
 * re-evaluation with the reconnection Jacobian).
 */
struct RESTIRDIReservoir
{
	/** World-space point on the chosen light (directional: P + dir * large). */
	float3 samplePosition;

	/**
	 * Unshadowed incident radiance from the sample (attenuated, no visibility).
	 */
	float3 sampleRadiance;

	/**
	 * Target function value p_hat of the held sample at this reservoir's
	 * surface.
	 */
	float targetPdf;

	/** Running sum of resampling weights. */
	float weightSum;

	/** Confidence (effective candidate count). */
	float M;

	/** Cached visibility of the held sample, in [0, 1]. */
	float visibility;

	/**
	 * Index of the chosen light (ShaderEntity index, or an emissive-triangle
	 * index OR'd with RESTIR_LIGHT_FLAG_EMISSIVE_TRIANGLE). Lets the denoiser
	 * reset history when a pixel's chosen light switches between frames, so the
	 * single visibility channel is never temporally averaged across two
	 * different lights (which would ghost). RESTIR_INVALID_LIGHT_INDEX if
	 * empty.
	 */
	uint lightIndex;

	/**
	 * Sample parameterization on the light's area, in [0, 1)^2. Together with
	 * lightIndex it re-resolves samplePosition/sampleRadiance at the light's
	 * current transform each frame, so a moving light's shadow tracks instead
	 * of lagging behind a stale resolved position.
	 */
	float2 uv;
};

/**
 * Resets a DI reservoir to empty.
 *
 * @return A DI reservoir holding no sample (all fields zero).
 */
inline RESTIRDIReservoir RESTIRDIReservoirInit()
{
	RESTIRDIReservoir r;
	r.samplePosition = float3(0, 0, 0);
	r.sampleRadiance = float3(0, 0, 0);
	r.targetPdf = 0;
	r.weightSum = 0;
	r.M = 0;
	r.visibility = 0;
	r.lightIndex = RESTIR_INVALID_LIGHT_INDEX;
	r.uv = float2(0, 0);
	return r;
}

/**
 * Unbiased contribution weight W = weightSum / (M * targetPdf).
 *
 * Guarded against division by zero (an empty reservoir yields 0).
 *
 * @param[in] r - The DI reservoir to weigh.
 *
 * @return The unbiased contribution weight W (>= 0).
 */
inline float RESTIRDIReservoirGetInvPdf(RESTIRDIReservoir r)
{
	const float denom = r.M * r.targetPdf;
	return denom > 0 ? r.weightSum / denom : 0;
}
#endif // __cplusplus

/**
 * GPU-storable packed DI reservoir (48 bytes), one per screen pixel.
 *
 * Backs the temporal history and spatial reuse buffers of the ReSTIR DI passes.
 */
struct RESTIRDIReservoirPacked
{
	uint4 data0;
	uint4 data1;
	uint4 data2;

#ifndef __cplusplus
	/**
	 * Packs a DI reservoir into this storage slot.
	 *
	 * @param[in] r - The DI reservoir to store.
	 */
	inline void store(RESTIRDIReservoir r)
	{
		data0.xyz = asuint(r.samplePosition);
		data0.w = Pack_R11G11B10_FLOAT(r.sampleRadiance);
		data1.x = asuint(r.targetPdf);
		data1.y = asuint(r.weightSum);
		// M (confidence, small) and visibility ([0, 1]) share one 16:16 word.
		data1.z = pack_half2(r.M, r.visibility);
		data1.w = r.lightIndex;
		// Sample uv (re-resolves the light each frame for moving-light
		// tracking).
		data2.x = pack_half2(r.uv.x, r.uv.y);
		data2.yzw = uint3(0, 0, 0);
	}

	/**
	 * Unpacks the stored DI reservoir.
	 *
	 * @return The reconstructed DI reservoir.
	 */
	inline RESTIRDIReservoir load()
	{
		RESTIRDIReservoir r;
		r.samplePosition = asfloat(data0.xyz);
		r.sampleRadiance = Unpack_R11G11B10_FLOAT(data0.w);
		r.targetPdf = asfloat(data1.x);
		r.weightSum = asfloat(data1.y);
		const float2 mv = unpack_half2(data1.z);
		r.M = mv.x;
		r.visibility = mv.y;
		r.lightIndex = data1.w;
		r.uv = unpack_half2(data2.x);
		return r;
	}
#endif // __cplusplus
};

/*
################################################################################
Screen-space GI reservoir
################################################################################
*/

/** ReSTIR GI initial-sample count per pixel (hemisphere rays traced). */
static const uint RESTIR_GI_INITIAL_CANDIDATES = 1;

/** Spatial-reuse neighbor count per pixel. */
static const uint RESTIR_GI_SPATIAL_SAMPLES = 4;

/** Spatial-reuse search radius, in pixels. */
static const float RESTIR_GI_SPATIAL_RADIUS = 16.0;

/**
 * Firefly clamp for the per-frame indirect irradiance (max luminance).
 *
 * The unbiased weight W = pi / cos(theta) blows up at grazing hemisphere
 * angles, so a single sample that escapes to the bright sky or lands on a
 * faintly lit surface can spike to a huge value the denoiser cannot absorb (a
 * flickering white speckle). Clamping the sample luminance suppresses those
 * fireflies at the cost of a small, bounded energy loss on genuinely bright
 * indirect bounces.
 */
static const float RESTIR_GI_FIREFLY_CLAMP = 8.0;

/**
 * Upper bound on the unbiased contribution weight W of a GI sample.
 *
 * The cosine-weighted source pdf gives W = pi / cos(theta), which explodes for
 * grazing hemisphere samples: a single such sample (or a neighbor carrying one
 * through reuse) becomes a bright spike the temporal denoiser cannot absorb,
 * seen as flickering blobs and an over-bright estimate on freshly disoccluded
 * pixels. Clamping W bounds the spike (a small dark bias at grazing angles),
 * and unlike the post-resolve luminance clamp it also bounds the weight a
 * reused neighbor contributes during resampling, before it can dominate.
 */
static const float RESTIR_GI_MAX_W = 8.0;

/**
 * SSAO-similarity falloff for the spatial resolve filter (per unit AO
 * difference, in exp2 units).
 *
 * A neighbor's weight is scaled by \[ 2^{-k\,|ao_c - ao_n|} \]: reuse across an
 * occlusion boundary (a corner, where ambient occlusion differs) is strongly
 * down-weighted, so a bright bounce sample in a tight corner is not spread onto
 * a differently occluded neighbor - the main source of the residual corner
 * fireflies. Applied only when an AO texture is available.
 *
 * References: EmbarkStudios/kajiya - rtdgi/restir_resolve.hlsl.
 */
static const float RESTIR_GI_SSAO_WEIGHT = 20.0;

/**
 * Spatial outlier-clamp factor for the resolve filter.
 *
 * Indirect diffuse irradiance is spatially smooth, so a pixel whose own
 * contribution is more than this many times its neighborhood mean is almost
 * certainly a firefly (an isolated bright bounce sample that the resolve
 * average did not dilute because its neighbors were rejected). Its contribution
 * is scaled back toward the neighbor mean, which removes single-pixel fireflies
 * without touching spatially coherent regions.
 */
static const float RESTIR_GI_OUTLIER_CLAMP = 4.0;

/**
 * Neighborhood color-clamp width for the GI temporal denoiser (in std devs).
 *
 * The reprojected history is clamped to the current frame's local irradiance
 * color box \[ [\mu - k\sigma,\ \mu + k\sigma] \] (a 3x3 neighborhood of the
 * resolved irradiance) before it is blended. A stale value - a dark trail that
 * a camera move drags across an area, or a lingering bright one - cannot then
 * ghost, because it is pulled back into the range the current frame actually
 * observes. Larger keeps more history (smoother, but more ghosting); smaller
 * reacts faster (less ghosting, more noise).
 *
 * A shared denoiser constant: both ReSTIR DI and GI clamp their reprojected
 * history to this color box. It catches temporal bleed across a lighting-
 * contrast edge (e.g. a multi-light shadow overlap) that the geometric
 * disocclusion and the A-SVGF gradient both miss, because the surface, depth
 * and normal are continuous there - only the lighting changes.
 *
 * References: EmbarkStudios/kajiya - rtdgi/temporal_filter.hlsl (color bbox).
 */
static const float RESTIR_DENOISE_HISTORY_CLAMP_K = 2.0;

/**
 * Sensitivity of the history-confidence reset to color staleness (shared by
 * ReSTIR DI and GI).
 *
 * The neighborhood color clamp bounds a ghost's magnitude, but the slow
 * 1/historyLength blend would still drag a clamped (stale) value out over many
 * frames. This measures how far the reprojected history lay outside the current
 * color box (relative to the local brightness) and resets the pixel's history
 * length toward 1 in proportion - a soft, color-based disocclusion on top of
 * the geometric one, so a detected ghost fades in about a frame instead of the
 * whole history. Larger resets more eagerly (less ghosting, more noise).
 */
static const float RESTIR_DENOISE_STALE_SENSITIVITY = 4.0;

/**
 * Number of edge-aware a-trous passes for the GI spatial denoiser.
 *
 * More than the DI pass count: a freshly disoccluded region has no temporal
 * history, so the only way to fill it is spatially, from converged neighbors.
 * Each extra pass doubles the a-trous stride, widening the reach so a
 * disocclusion band survives fast camera motion. Indirect diffuse is
 * low-frequency, so the extra blur costs little detail, and the variance /
 * depth / normal edge-stops keep converged regions sharp.
 */
static const uint RESTIR_GI_DENOISE_ATROUS_PASSES = 6;

/**
 * Low-history variance boost for the a-trous spatial denoiser (shared by ReSTIR
 * DI and GI).
 *
 * On a disoccluded pixel (short history) the a-trous edge-stop must open so the
 * pixel can be filled from its neighbors. The luminance-difference tolerance is
 * scaled up by up to this factor as history length falls toward 1, and relaxes
 * to 1 as the history reaches the denoise cap - the SVGF low-history spatial
 * fallback, made explicit.
 *
 * For DI this is what averages a freshly disoccluded multi-light shadow overlap
 * across its black/over-bright per-frame swing frames, so the region resolves
 * to the correct mean at once instead of flashing.
 *
 * References: Schied et al. 2017, "Spatiotemporal Variance-Guided Filtering".
 */
static const float RESTIR_DENOISE_LOWHISTORY_BOOST = 8.0;

#ifndef __cplusplus
/**
 * A screen-space ReSTIR GI reservoir.
 *
 * Holds one indirect-lighting *sample point*: the world-space point a
 * hemisphere ray hit, its normal, and the radiance leaving it toward this
 * pixel's surface. Unlike the DI reservoir (which references a light), GI
 * resamples surface points in the scene, so reuse across pixels needs the hit
 * position and normal to re-evaluate the target function and the reconnection
 * Jacobian at a different visible point.
 *
 * References: Ouyang et al. 2021, "ReSTIR GI: Path Resampling for Real-Time
 * Path Tracing".
 */
struct RESTIRGIReservoir
{
	/** World-space indirect hit point (the sample point). */
	float3 samplePosition;

	/**
	 * World-space normal at the hit point (for reuse Jacobian + validation).
	 */
	float3 sampleNormal;

	/**
	 * Radiance leaving the hit point toward this pixel's surface (view-
	 * independent diffuse outgoing radiance; no receiver albedo).
	 */
	float3 sampleRadiance;

	/** Target function value p_hat of the held sample at the owning pixel. */
	float targetPdf;

	/** Running sum of resampling weights. */
	float weightSum;

	/** Confidence (effective candidate count). */
	float M;
};

/**
 * Resets a GI reservoir to empty.
 *
 * @return A GI reservoir holding no sample (all fields zero).
 */
inline RESTIRGIReservoir RESTIRGIReservoirInit()
{
	RESTIRGIReservoir r;
	r.samplePosition = float3(0, 0, 0);
	r.sampleNormal = float3(0, 0, 0);
	r.sampleRadiance = float3(0, 0, 0);
	r.targetPdf = 0;
	r.weightSum = 0;
	r.M = 0;
	return r;
}

/**
 * Unbiased contribution weight W = weightSum / (M * targetPdf).
 *
 * @param[in] r - The GI reservoir to weigh.
 *
 * @return The unbiased contribution weight W (>= 0).
 */
inline float RESTIRGIReservoirGetInvPdf(RESTIRGIReservoir r)
{
	const float denom = r.M * r.targetPdf;
	return denom > 0 ? r.weightSum / denom : 0;
}
#endif // __cplusplus

/**
 * GPU-storable packed GI reservoir (48 bytes), one per screen pixel.
 *
 * Backs the trace output, temporal history and spatial reuse buffers of the
 * ReSTIR GI passes.
 */
struct RESTIRGIReservoirPacked
{
	uint4 data0;
	uint4 data1;
	uint4 data2;

#ifndef __cplusplus
	/**
	 * Packs a GI reservoir into this storage slot.
	 *
	 * @param[in] r - The GI reservoir to store.
	 */
	inline void store(RESTIRGIReservoir r)
	{
		data0.xyz = asuint(r.samplePosition);
		data0.w = Pack_R11G11B10_FLOAT(r.sampleRadiance);
		data1.xyz = asuint(r.sampleNormal);
		data1.w = asuint(r.targetPdf);
		data2.x = asuint(r.weightSum);
		data2.y = asuint(r.M);
		data2.zw = uint2(0, 0);
	}

	/**
	 * Unpacks the stored GI reservoir.
	 *
	 * @return The reconstructed GI reservoir.
	 */
	inline RESTIRGIReservoir load()
	{
		RESTIRGIReservoir r;
		r.samplePosition = asfloat(data0.xyz);
		r.sampleRadiance = Unpack_R11G11B10_FLOAT(data0.w);
		r.sampleNormal = asfloat(data1.xyz);
		r.targetPdf = asfloat(data1.w);
		r.weightSum = asfloat(data2.x);
		r.M = asfloat(data2.y);
		return r;
	}
#endif // __cplusplus
};

/*
################################################################################
Push constants
################################################################################
*/

/**
 * Parameters shared by the ReSTIR DI compute passes
 * (initial / temporal / spatial).
 */
struct RESTIRDIPushConstants
{
	uint2 resolution;
	float2 resolutionRcp;
	uint frameIndex;
	uint candidateCount;
	uint spatialSampleCount;
	float spatialRadius;
};

/**
 * Parameters shared by the ReSTIR GI compute passes
 * (trace / temporal / spatial).
 */
struct RESTIRGIPushConstants
{
	uint2 resolution;
	float2 resolutionRcp;
	uint frameIndex;
	uint candidateCount;
	uint spatialSampleCount;
	float spatialRadius;

	/** Ray-tracing instance inclusion mask for the indirect bounce ray. */
	uint instanceInclusionMask;
	uint pad0;
	uint pad1;
	uint pad2;
};

/**
 * Parameters for one pass of the spatial visibility a-trous denoiser.
 *
 * The denoiser runs RESTIR_DI_DENOISE_ATROUS_PASSES passes with a doubling tap
 * stride. The first pass reads the temporally accumulated visibility straight
 * from the reservoir buffer, the last writes the filtered result back into it;
 * intermediate passes ping-pong two (visibility, variance) scratch textures.
 */
struct RESTIRDIDenoisePushConstants
{
	uint2 resolution;
	float2 resolutionRcp;

	/** Tap stride for this a-trous pass (1, 2, 4, ...). */
	uint stepSize;

	/** 1 on the first pass: read the input from the reservoir + moments. */
	uint inputFromReservoir;

	/** 1 on the last pass: write the filtered visibility into the reservoir. */
	uint outputToReservoir;

	/**
	 * Low-history variance boost (0 disables it). Set by ReSTIR GI to
	 * RESTIR_GI_DENOISE_LOWHISTORY_BOOST so disoccluded pixels blur wider; left
	 * 0 by ReSTIR DI.
	 */
	float lowHistoryBoost;
};

/**
 * Parameters for the light pre-sampling pass that fills the light tiles.
 */
struct RESTIRPresamplePushConstants
{
	uint frameIndex;

	/**
	 * Number of emissive triangles currently in the light buffer (0 = none).
	 */
	uint emissiveTriangleCount;

	/** Bindless index of the RESTIRLightRef tile buffer (UAV). */
	int lightTileBuffer;

	/**
	 * Bindless index of the RESTIREmissiveTriangle buffer (SRV; -1 if none).
	 */
	int emissiveTriangleBuffer;
};

#endif // WI_SHADERINTEROP_RESTIR_H
