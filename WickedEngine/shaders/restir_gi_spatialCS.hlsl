#include "globals.hlsli"
#include "restir_giHF.hlsli"

/**
 * ReSTIR GI - spatial resolve filter.
 *
 * Resolves each pixel's indirect irradiance as a geometry-weighted average of
 * its own and a few neighbors' temporal-reservoir contributions. Averaging (as
 * in kajiya's resolve) DILUTES a bright, localized bounce sample across the
 * neighborhood, whereas a reservoir resample ADOPTS it wholesale - which is
 * what concentrated such samples into flickering firefly blobs. This still
 * fills freshly disoccluded pixels from their neighbors and lowers variance -
 * the benefits of spatial reuse - without spreading fireflies.
 *
 * Neighbors are rejected on depth and normal mismatch (different surface) and
 * weighted by their geometric similarity to the center pixel; the per-frame
 * random offsets let the temporal accumulation average different footprints.
 *
 * References: EmbarkStudios/kajiya - rtdgi/restir_resolve.hlsl.
 */

PUSHCONSTANT(push, RESTIRGIPushConstants);

ByteAddressBuffer reservoirInput : register(t0); // temporal reservoirs

RWTexture2D<float3> irradianceOutput : register(u0); // per-frame indirect E

[numthreads(8, 8, 1)]
void main(uint3 DTid : SV_DispatchThreadID)
{
	if (DTid.x >= push.resolution.x || DTid.y >= push.resolution.y)
		return;

	const uint2 pixel = DTid.xy;
	const uint flatIndex = pixel.y * push.resolution.x + pixel.x;

	float3 irradiance = 0;

	const float depth = texture_depth[pixel];
	[branch]
	if (depth > 0)
	{
		const float2 uv = (pixel + 0.5) * push.resolutionRcp;
		const float3 P = reconstruct_position(uv, depth);
		const float3 N = decode_normal(texture_normal_roughness[pixel]);
		const float linearCenter = compute_lineardepth(depth);

		RNG rng;
		rng.init(pixel, GetFrame().frame_count + 2u);

		// SSAO similarity down-weights reuse across occlusion boundaries
		// (corners), where a bright bounce sample must not spread onto a
		// differently occluded neighbor - the main source of the residual
		// corner fireflies. Only available when an AO texture is bound.
		const int aoIndex = GetCamera().texture_ao_index;
		const float centerAO = (aoIndex >= 0)
			? bindless_textures_half4[descriptor_index(aoIndex)].SampleLevel(
				sampler_linear_clamp, uv, 0).r
			: 1.0;

		// Accumulate the neighbors separately from the center so the center can
		// be outlier-clamped against the neighbor mean below.
		float3 neighborWeighted = 0;
		float neighborWsum = 0;

		for (uint i = 0; i < push.spatialSampleCount; ++i)
		{
			const float2 offset =
				(rng.next_float2() * 2 - 1) * push.spatialRadius;
			const int2 neighborPixel = int2(pixel) + int2(round(offset));

			[branch]
			if (any(neighborPixel < 0) ||
				any(neighborPixel >= int2(push.resolution)))
				continue;

			const float neighborDepth = texture_depth[neighborPixel];
			if (neighborDepth <= 0)
				continue;

			// Reject neighbors on a different surface; weight the rest by their
			// geometric similarity so a mildly different neighbor contributes
			// less than a coplanar one.
			const float linearNeighbor = compute_lineardepth(neighborDepth);
			if (abs(linearNeighbor - linearCenter) > 0.05 * linearCenter)
				continue;

			const float3 neighborN =
				decode_normal(texture_normal_roughness[neighborPixel]);
			const float ndot = dot(neighborN, N);
			if (ndot < 0.9)
				continue;

			const float wNormal = pow(saturate(ndot),
				RESTIR_DI_DENOISE_NORMAL_POWER);
			const float wDepth = exp(-abs(linearNeighbor - linearCenter) /
				(RESTIR_DI_DENOISE_DEPTH_SCALE * linearCenter + 1e-3));
			float w = wNormal * wDepth;

			[branch]
			if (aoIndex >= 0)
			{
				const float2 neighborUV =
					(neighborPixel + 0.5) * push.resolutionRcp;
				const float neighborAO =
					bindless_textures_half4[descriptor_index(aoIndex)].SampleLevel(
						sampler_linear_clamp, neighborUV, 0).r;
				w *= exp2(-RESTIR_GI_SSAO_WEIGHT * abs(centerAO - neighborAO));
			}

			const uint neighborFlat =
				neighborPixel.y * push.resolution.x + neighborPixel.x;
			neighborWeighted += RESTIRGIResolve(
				RESTIRGIReservoirLoad(reservoirInput, neighborFlat), P, N) * w;
			neighborWsum += w;
		}

		float3 centerContribution = RESTIRGIResolve(
			RESTIRGIReservoirLoad(reservoirInput, flatIndex), P, N);

		// Spatial outlier clamp: an isolated pixel whose contribution is far
		// brighter than its neighbor mean is a firefly (indirect diffuse is
		// smooth); scale it back toward the neighbor mean. Coherent regions
		// (where the neighbors are equally bright) are untouched.
		[branch]
		if (neighborWsum > 0)
		{
			const float3 neighborMean = neighborWeighted / neighborWsum;
			const float centerLuma =
				dot(centerContribution, float3(0.2126, 0.7152, 0.0722));
			const float maxLuma = dot(neighborMean,
				float3(0.2126, 0.7152, 0.0722)) * RESTIR_GI_OUTLIER_CLAMP + 1e-3;
			if (centerLuma > maxLuma)
				centerContribution *= maxLuma / centerLuma;
		}

		float3 irr = (centerContribution + neighborWeighted) /
			(1.0 + neighborWsum);
		irr = (any(isnan(irr)) || any(isinf(irr))) ? (float3)0 : irr;

		// Firefly clamp: bound the resolved irradiance luminance (preserving
		// hue) as a final guard.
		const float lum = dot(irr, float3(0.2126, 0.7152, 0.0722));
		if (lum > RESTIR_GI_FIREFLY_CLAMP)
			irr *= RESTIR_GI_FIREFLY_CLAMP / lum;

		irradiance = irr;
	}

	irradianceOutput[pixel] = irradiance;
}
