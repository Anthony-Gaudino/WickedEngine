#include "globals.hlsli"
#include "restir_giHF.hlsli"

/**
 * ReSTIR GI - temporal resampling pass.
 *
 * Merges the previous frame's reservoir (reprojected through the velocity
 * buffer, rejected on disocclusion) into this frame's initial reservoir, then
 * resolves the temporally accumulated indirect irradiance for the denoiser.
 *
 * Reuse across frames tracks the same surface point, so the reconnection
 * Jacobian is ~1 here (the source and target visible points coincide); the full
 * Jacobian is applied by the spatial pass, whose neighbors are genuinely
 * different points. The confidence is bounded (RESTIR_MAX_HISTORY_LENGTH), with
 * weightSum scaled together with M so the unbiased weight W = weightSum /
 * (M * targetPdf) is preserved - capping M alone would inflate W every frame
 * and blow up through the history.
 */

PUSHCONSTANT(push, RESTIRGIPushConstants);

ByteAddressBuffer reservoirInitial : register(t0);
ByteAddressBuffer reservoirHistory : register(t1);

RWByteAddressBuffer reservoirOutput : register(u0);
RWTexture2D<float3> irradianceOutput : register(u1); // per-frame indirect E
RWTexture2D<float> gradientOutput : register(u2);     // antilag (0 for now)

[numthreads(8, 8, 1)]
void main(uint3 DTid : SV_DispatchThreadID)
{
	if (DTid.x >= push.resolution.x || DTid.y >= push.resolution.y)
		return;

	const uint2 pixel = DTid.xy;
	const uint flatIndex = pixel.y * push.resolution.x + pixel.x;

	RESTIRGIReservoir reservoir =
		RESTIRGIReservoirLoad(reservoirInitial, flatIndex);
	float3 irradiance = 0;

	const float depth = texture_depth[pixel];
	[branch]
	if (depth > 0)
	{
		const float2 uv = (pixel + 0.5) * push.resolutionRcp;
		const float3 P = reconstruct_position(uv, depth);
		const float3 N = decode_normal(texture_normal_roughness[pixel]);

		const float2 velocity = texture_velocity[pixel];
		const float2 prevUV = uv + velocity;

		[branch]
		if (all(prevUV >= 0) && all(prevUV <= 1))
		{
			const uint2 prevPixel =
				min(uint2(prevUV * push.resolution), push.resolution - 1);

			// Disocclusion test: reject history that reprojects onto a
			// different surface (depth edge crossing) or a differently oriented
			// one (rotation past an edge). No previous normal buffer is kept,
			// so the current normal at the reprojected pixel is used as a
			// geometric proxy - the same approximation the DI temporal pass
			// makes.
			const float prevDepth = texture_depth_history.SampleLevel(
				sampler_point_clamp, prevUV, 0);
			const float linearCur = compute_lineardepth(depth);
			const float linearPrev = compute_lineardepth(prevDepth);
			const float3 prevN = decode_normal(texture_normal_roughness[prevPixel]);
			const bool consistent =
				abs(linearCur - linearPrev) <= 0.05 * linearCur &&
				dot(prevN, N) >= 0.9;

			[branch]
			if (consistent)
			{
				const uint prevFlat =
					prevPixel.y * push.resolution.x + prevPixel.x;

				RESTIRGIReservoir history =
					RESTIRGIReservoirLoad(reservoirHistory, prevFlat);

				if (history.M > (float)RESTIR_MAX_HISTORY_LENGTH)
				{
					history.weightSum *=
						(float)RESTIR_MAX_HISTORY_LENGTH / history.M;
					history.M = (float)RESTIR_MAX_HISTORY_LENGTH;
				}

				[branch]
				if (history.M > 0)
				{
					RNG rng;
					rng.init(pixel, GetFrame().frame_count + 1u);

					const float targetAtSelf = RESTIRGITargetFunction(
						history.sampleRadiance, history.samplePosition, P, N);

					// Temporal reuse tracks the same surface point, so J = 1.
					RESTIRGIReservoirMerge(
						reservoir, history, targetAtSelf, 1.0, rng);
				}
			}
		}

		float3 irr = RESTIRGIResolve(reservoir, P, N);
		irr = (any(isnan(irr)) || any(isinf(irr))) ? (float3)0 : irr;

		// Firefly clamp: bound the resolved irradiance luminance (preserving
		// hue) so a grazing-angle weight cannot spike a speckle the denoiser
		// cannot absorb.
		const float lum = dot(irr, float3(0.2126, 0.7152, 0.0722));
		if (lum > RESTIR_GI_FIREFLY_CLAMP)
			irr *= RESTIR_GI_FIREFLY_CLAMP / lum;

		irradiance = irr;
	}

	RESTIRGIReservoirStore(reservoirOutput, flatIndex, reservoir);
	irradianceOutput[pixel] = irradiance;
	// No antilag gradient yet; a zero gradient makes the denoiser rely on pure
	// temporal accumulation with depth/normal disocclusion.
	gradientOutput[pixel] = 0;
}
