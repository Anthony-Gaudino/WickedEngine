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
Push constants
################################################################################
*/

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
