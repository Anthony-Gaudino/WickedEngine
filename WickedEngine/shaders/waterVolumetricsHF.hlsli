#ifndef WI_WATERVOLUMETRICS_HF
#define WI_WATERVOLUMETRICS_HF

/**
 * Shared sun shaft march through water.
 *
 * Provides `MarchWaterSunShafts`, the reusable "walk a submerged ray segment
 * and accumulate the sunlight scattered towards the eye" core. It is the
 * physical counterpart to the procedural god rays in `underwaterCS.hlsl`: the
 * shafts here come from the sun's own shadow cascades, so real geometry casts
 * them and real geometry blocks them.
 *
 * The segment and the medium are both parameters rather than read from the
 * ocean, so the same march can serve the underwater post pass (whole view ray)
 * and, later, the directional volumetric light pass (only the part of its ray
 * below the waterline).
 *
 * @note Sampling the shadow cascades also reads the *transparent* shadow layer,
 *       which is where the ocean surface writes its caustics (see the
 *       `SHADOWMAPRENDERING` branch of `oceanSurfacePS.hlsl`). The shafts
 *       therefore pick up caustic modulation without asking for it. Only the
 *       first couple of cascades render the ocean, so that modulation fades
 *       with distance.
 */

#include "globals.hlsli"
#include "lightingHF.hlsli"
#include "skyAtmosphere.hlsli"

/**
 * Optical depth past which the march stops.
 *
 * \f$e^{-5} \approx 0.007\f$, so beyond this the water has swallowed all but a
 * fraction of a percent and further steps buy nothing.
 */
static const float WATER_SHAFT_OPTICAL_DEPTH_CUTOFF = 5.0;

/**
 * Largest distance the dithered ray start may be offset by, in metres.
 *
 * The offset is a fraction of a step, so it grows with the step size. Left
 * unbounded, a long segment scatters neighbouring pixels far enough apart that
 * the dither's Bayer matrix becomes visible as a pattern on screen rather than
 * as noise.
 */
static const float WATER_SHAFT_MAX_DITHER_OFFSET = 10.0;

/**
 * Accumulates sunlight scattered towards the eye along a submerged ray.
 *
 * Single scattering: each step gathers sunlight that reaches the sample, turns
 * it towards the eye through the medium's phase function, and attenuates it
 * back along the ray already walked:
 * \[
 * L = \int_0^d T_{view}(s)\, \sigma_s\, p(\theta)\, V(s)\, T_{sun}(s)\, ds
 * \]
 * with \f$T\f$ Beer-Lambert transmittance and \f$V\f$ the shadow term.
 *
 * Example usage:
 * @code
 * const float3 shafts = MarchWaterSunShafts(
 *     campos, rayDir, waterPathLength,
 *     refractedSunDir, GetSunColor(),
 *     sigmaS, sigmaT, ocean.scattering.a,
 *     16, DTid.xy);
 * color.rgb += shafts;
 * @endcode
 *
 * @param[in] rayOrigin - Start of the submerged segment, in world space.
 * @param[in] rayDir - Normalized direction from the eye into the scene.
 * @param[in] rayLength - Length of the submerged segment (in metres).
 * @param[in] sunDir - Sun direction *after* refraction at the surface, so
 *                     pointing down into the water.
 * @param[in] sunColor - Sun radiance before the water attenuates it.
 * @param[in] sigmaS - Scattering coefficient of the water (1/m per channel).
 * @param[in] sigmaT - Extinction coefficient of the water (1/m per channel).
 * @param[in] phaseG - Henyey-Greenstein asymmetry of the medium.
 * @param[in] sampleCount - Steps along the segment. Must be at least 1.
 * @param[in] pixel - Dispatch pixel, used to dither the start and to sample
 *                    the shadow cascades.
 *
 * @return Radiance scattered towards the eye. Zero when there is no usable
 *         shadow casting sun, so the caller can add it unconditionally.
 */
float3 MarchWaterSunShafts(
	float3 rayOrigin,
	float3 rayDir,
	float rayLength,
	float3 sunDir,
	float3 sunColor,
	float3 sigmaS,
	float3 sigmaT,
	float phaseG,
	uint sampleCount,
	uint2 pixel
)
{
	// Reuses the engine's own "is this light usable for volumetrics" test,
	// which also hands back the last cascade worth walking.
	ShaderEntity light = (ShaderEntity)0;
	uint furthestCascade = 0;
	if (!furthest_cascade_volumetrics(light, furthestCascade))
		return 0;

	const float waterHeight = GetWeather().ocean.water_height;

	// Marching the full ray is pointless and actively harmful: where the depth
	// buffer holds the far plane the segment is kilometres long, and dividing
	// that into a handful of steps puts the samples so far apart that the
	// dithered start alone lands neighbouring pixels in completely different
	// places - which reads as the dither's Bayer pattern printed on screen.
	// Water swallows everything well before then, so stop once the view
	// transmittance has essentially died.
	const float opticalRange = WATER_SHAFT_OPTICAL_DEPTH_CUTOFF
		/ max(0.00001, min(sigmaT.r, min(sigmaT.g, sigmaT.b)));

	// A ray heading up leaves the water entirely, and everything past that
	// point is air which the submersion test below discards anyway. Stopping
	// there matters most exactly where the segment is shortest - a camera just
	// under the surface looking up - because otherwise the step stays long
	// while the water is thin, and the dither alone decides whether a pixel
	// samples that thin layer at all. Neighbouring pixels then disagree, which
	// reads as a halftone pattern.
	//
	// Continuous across the horizontal: as rayDir.y falls to zero the exit
	// distance grows without bound and the optical range takes over instead.
	float exitDistance = rayLength;
	[branch]
	if (rayDir.y > 0.0001)
	{
		exitDistance = max(0, (waterHeight - rayOrigin.y) / rayDir.y);
	}

	const float marchLength = min(min(rayLength, opticalRange), exitDistance);

	const float stepSize = marchLength / (float)max(1u, sampleCount);
	if (stepSize <= 0)
		return 0;

	// The sun is directional, so the scattering angle does not change along the
	// ray and the phase function is a constant.
	//
	// Reduce the asymmetry by how much of the light has been scattered, exactly
	// as the veiling term in underwaterCS does. Marine particulates peak some
	// 300x above isotropic, and summing that raw over a march produces a spike
	// at the sun bright enough to blow the image out; similarity theory says a
	// single scattering march standing in for a multiply scattering medium
	// should use the reduced value instead.
	//
	// The trailing uniform phase keeps this on the same scale as the rest of
	// the underwater pass, which applies one too (as the engine's height fog
	// does). Without it these shafts would sit some 4*PI above the veiling
	// light they are meant to modulate.
	const float3 albedo = sigmaS / sigmaT;
	const half multiScatter = (half)saturate(albedo.g);
	const half reducedG = (half)phaseG * (1 - multiScatter);
	const half cosTheta = (half)dot(rayDir, sunDir);
	const half phase = HgPhase(reducedG, cosTheta) * UniformPhase();

	// Sunlight crosses a slant path from the surface down to a sample, not just
	// the vertical drop. Refraction bounds the ray to the critical angle, so
	// that path can never be more than about 1.5x the depth.
	const float sunSlant = max(abs(sunDir.y), 0.66);

	// Marching in fixed steps bands badly at low sample counts; offsetting the
	// start per pixel turns the banding into noise, as the directional
	// volumetric light pass does.
	float3 P = rayOrigin + rayDir * min(
		stepSize * dither((min16uint2)pixel), WATER_SHAFT_MAX_DITHER_OFFSET
	);

	const float3 stepTransmittance = exp(-stepSize * sigmaT);
	float3 viewTransmittance = 1;
	float3 scattered = 0;

	[loop]
	for (uint i = 0; i < sampleCount; ++i)
	{
		// Walk out to the first cascade that contains the sample: the nearest
		// one covering it is also the highest resolution one.
		float3 shadow = 1;
		for (uint cascade = 0; cascade <= furthestCascade; ++cascade)
		{
			const float3 shadow_pos = mul(
				load_entitymatrix(light.GetMatrixIndex() + cascade),
				float4(P, 1)
			).xyz; // ortho matrix, no divide by .w
			const float3 shadow_uv = clipspace_to_uv(shadow_pos);

			[branch]
			if (is_saturated(shadow_uv))
			{
				shadow = shadow_2D(
					light, shadow_pos.z, shadow_uv.xy, cascade,
					(min16uint2)pixel
				);
				break;
			}
		}

		// Whatever light reaches this sample had to come down through the water
		// above it first, which is what makes the shafts fade out with depth
		// and lose their warm end on the way.
		const float sunDepth = waterHeight - P.y;
		const float3 sunTransmittance =
			exp(-(max(0, sunDepth) / sunSlant) * sigmaT);

		// There is no medium above the surface, so a ray that leaves the water
		// part way along must stop contributing rather than keep scattering in
		// air. Tested per sample instead of by shortening the segment, because
		// the distance at which the ray exits jumps discontinuously at the
		// horizontal and would draw a seam across the view.
		const float submerged = sunDepth >= 0 ? 1.0 : 0.0;

		scattered += submerged
			* viewTransmittance * shadow * sunTransmittance * sigmaS * phase * stepSize;

		viewTransmittance *= stepTransmittance;
		P += rayDir * stepSize;
	}

	return max(0, scattered * sunColor);
}

#endif // WI_WATERVOLUMETRICS_HF
