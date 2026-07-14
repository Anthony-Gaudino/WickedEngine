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
RWTexture2D<float> gradientOutput : register(u1);     // antilag (0 for now)

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

		// Center reservoir contributes with full weight.
		float3 weighted = RESTIRGIResolve(
			RESTIRGIReservoirLoad(reservoirInput, flatIndex), P, N);
		float wsum = 1;

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
			const float w = wNormal * wDepth;

			const uint neighborFlat =
				neighborPixel.y * push.resolution.x + neighborPixel.x;
			weighted += RESTIRGIResolve(
				RESTIRGIReservoirLoad(reservoirInput, neighborFlat), P, N) * w;
			wsum += w;
		}

		float3 irr = weighted / max(wsum, 1e-4);
		irr = (any(isnan(irr)) || any(isinf(irr))) ? (float3)0 : irr;

		// Firefly clamp: bound the resolved irradiance luminance (preserving
		// hue) as a final guard.
		const float lum = dot(irr, float3(0.2126, 0.7152, 0.0722));
		if (lum > RESTIR_GI_FIREFLY_CLAMP)
			irr *= RESTIR_GI_FIREFLY_CLAMP / lum;

		irradiance = irr;
	}

	irradianceOutput[pixel] = irradiance;
	// No antilag gradient yet; a zero gradient makes the denoiser rely on pure
	// temporal accumulation with depth/normal disocclusion.
	gradientOutput[pixel] = 0;
}
