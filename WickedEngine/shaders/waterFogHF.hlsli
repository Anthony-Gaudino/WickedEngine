#ifndef WI_WATERFOG_HF
#define WI_WATERFOG_HF

/**
 * The water's fog over one view segment, per channel.
 *
 * What the eye receives from beyond a stretch of water is what survived it plus
 * what the water scattered back in its place:
 * \[
 * L = T\,L_{beyond} + \bar{L}\,\omega_0\,(1 - T),
 * \qquad T = e^{-\sigma_t d}
 * \]
 * The first term darkens and the second hazes, and they are separated by the
 * single-scattering albedo \( \omega_0 = \sigma_s / \sigma_t \) - absorption
 * only destroys, scattering redirects. Clear deep water has a low albedo and
 * reads as a dark filter; turbid water has a high one and reads as a bright
 * milky haze, which is why murk kills contrast long before it kills brightness.
 *
 * Kept apart from the engine's height fog rather than folded into it. `GetFog`
 * (`fogHF.hlsli`) returns a colour and a SCALAR weight, and a scalar cannot say
 * that red is gone at two metres while blue is still going at twelve - which is
 * the whole character of water. The two are also physically exclusive media, so
 * they compose rather than merge: height fog first, then this.
 *
 * Deliberately free of anything screen-space beyond the god rays, so the same
 * fog can be evaluated from a post pass over a reconstructed depth or from a
 * draw shader over its own fragment.
 */

#include "waterVolumetricsHF.hlsli"

/**
 * Screen-space radial stripes, swept around a point.
 *
 * A stylised god ray pattern: it is anchored to where the sun projects on
 * screen rather than to the world, and nothing occludes it, so it is a look
 * rather than a simulation. Kept because it costs only arithmetic.
 *
 * References:
 * https://www.shadertoy.com/view/XdyfR1
 *
 * @param[in] lightScreen - Normalized device position of the light.
 * @param[in] ndc - Normalized device position of the pixel.
 * @param[in] uv - Screen space UV coordinates (0-1) of the pixel.
 * @param[in] time - Animation time, in seconds.
 * @param[in] rayLength - Exponent on the vertical falloff.
 * @param[in] rayFrequency - Angular frequency of the stripes.
 *
 * @return Stripe intensity, nominally in [0, 1].
 */
float GodRays(
	float2 lightScreen,
	float2 ndc,
	float2 uv,
	float time,
	float rayLength = 0.1,
	float rayFrequency = 48.0
)
{
	const float2 godRayOrigin = ndc - lightScreen;
	const float rayInputFunc = atan2(godRayOrigin.y, godRayOrigin.x);

	float light =
		(sin(rayInputFunc * rayFrequency + time * -2.25) * 0.5 + 0.5);
	light = 0.5 * (light + (sin(rayInputFunc * 13.0 + time) * 0.5 + 0.5));
	light *= pow(uv.y, rayLength);

	return pow(light, 1.75);
}

/**
 * What the water does to one segment of a view ray.
 *
 * Example usage:
 * @code
 * const WaterFog fog = MakeWaterFog(
 *     MakeWaterVolumetrics(1), path, V, eyeDepth, uv, ndc, true);
 * color.rgb = color.rgb * fog.transmittance + fog.inscatter;
 * @endcode
 */
struct WaterFog
{
	/**
	 * Per-channel survival of whatever lies beyond the segment, in [0, 1].
	 *
	 * Multiply the incoming radiance by this.
	 */
	float3 transmittance;

	/**
	 * Veiling radiance the water puts in its place.
	 *
	 * Already weighted by the albedo and by how much was extinguished, so it
	 * is added rather than blended.
	 *
	 * @note Full precision deliberately. The colour is authored in half, but
	 *       the albedo weighting that scales it is not, and rounding their
	 *       product back down loses the low end of the haze where the water is
	 *       nearly clear.
	 */
	float3 inscatter;

	/**
	 * Whether this fog does anything at all.
	 *
	 * @return true when the water is present and the segment has length.
	 */
	bool IsActive()
	{
		return any(transmittance < 1);
	}
};

/**
 * Builds the water's fog over a segment of known length.
 *
 * The medium is passed in rather than read here, so that a caller which
 * already has one - a march, or a fragment that has decided how submerged it
 * is - cannot end up with a second, subtly different set of coefficients.
 *
 * The inscattered colour is the water's authored base plus the sun, and the sun
 * is put through everything that happened to it on the way: the atmosphere it
 * crossed, the share the surface reflected away instead of admitting, the water
 * column above the eye, the phase function, and the god rays.
 *
 * Example usage:
 * @code
 * const WaterFog fog = MakeWaterFog(
 *     MakeWaterVolumetrics(1), path, V, eyeDepth, uv, ndc, true);
 * @endcode
 *
 * @param[in] medium - The water, supplying sigmaS and sigmaT.
 * @param[in] path - Length of the segment lying in water (in metres).
 * @param[in] toEye - Normalized direction from the far end towards the eye.
 * @param[in] eyeDepth - Depth of the EYE below the surface (in metres), which
 *                       is the column the sun crossed before it could scatter
 *                       towards the eye at all.
 * @param[in] screenUV - Screen space UV coordinates (0-1), for the god rays.
 * @param[in] ndc - Normalized device position of the pixel, for the god rays.
 * @param[in] godRays - Whether to modulate the inscatter with god rays.
 *
 * @return The fog for this segment.
 */
WaterFog MakeWaterFog(
	WaterVolumetrics medium,
	float path,
	float3 toEye,
	float eyeDepth,
	float2 screenUV,
	float2 ndc,
	bool godRays
)
{
	const ShaderOcean ocean = GetWeather().ocean;

	WaterFog fog;
	fog.transmittance = exp(-path * medium.sigmaT);

	// Only the scattered share of what was extinguished comes back.
	const float3 scatterAlbedo = medium.sigmaS / medium.sigmaT;
	const float3 inscatterAmount = scatterAlbedo * (1 - fog.transmittance);

	const half3 L = GetSunDirection();
	const half3 refractedLightDir =
		refract(-L, float3(0, 1, 0), 1.0 / WATER_REFRACTIVE_INDEX);

	half3 inscatteringColor = GetSunColor();

	// Apply atmosphere transmittance:
	if (GetFrame().options & OPTION_BIT_REALISTIC_SKY)
	{
		// 0 for position since fog is centered around world center
		inscatteringColor *= GetAtmosphericLightTransmittance(
			GetWeather().atmosphere, 0, L, texture_transmittancelut);
	}

	// Phase function of the medium, at a REDUCED asymmetry.
	//
	// Marine particulates scatter extremely far forward (g around 0.9, see
	// WaterMedium::PhaseAsymmetry), but that is the phase of a SINGLE
	// scattering event. What is being modulated here is the whole veiling
	// radiance integrated along the path, and by the time light has bounced
	// several times it has lost its original direction entirely. Feeding the
	// raw single-scattering g into an integrated term puts an enormous spike at
	// the sun: its peak is some 300x the isotropic value, so even a little of
	// it swamps everything else.
	//
	// Similarity theory says to use a reduced asymmetry instead, scaled down by
	// how much of the light has been scattered at all - which is exactly the
	// veiling term computed above, since it already combines "how much was
	// extinguished" with "how much of that was redirected rather than
	// absorbed". Clear water keeps a sharp sun glow; turbid water, having
	// scattered far more, spreads it into an even haze.
	//
	// Taken from the raw authored asymmetry rather than the medium's own
	// phaseG, which is reduced by the albedo instead - the right choice for a
	// march, which has no single path length to reduce by, and the wrong one
	// here where there is exactly one.
	const half multiScatter = (half)saturate(inscatterAmount.g);
	const half phaseG = (half)ocean.scattering.a * (1 - multiScatter);
	const half cosTheta = dot(toEye, -refractedLightDir);
	inscatteringColor *= HgPhase(phaseG, cosTheta);

	// Uniform base phase, matching what the engine's height fog does for a
	// constant medium (see GetFog in fogHF.hlsli):
	inscatteringColor *= UniformPhase();

	// Sunlight has already crossed the water column above the eye before it can
	// scatter towards it, so it arrives filtered by the same medium. Red goes
	// first, which is what tints the scatter green-blue with depth without any
	// authored colour:
	inscatteringColor *= exp(-eyeDepth * medium.sigmaT);

	// And some of it never got in: the surface reflects a share away instead of
	// refracting it down, nearly all of it once the sun is low. The same term
	// the lit surfaces use, so the haze and the sea bed agree about how much
	// sun reached the water - without it the water would still glow at sunset
	// while everything in it went dark.
	//
	// RefractIntoWater rather than refractedLightDir above: the two agree while
	// the sun is up, but refract() keeps bending a sun that has already set
	// into a downward ray, whose cosine would read as light still getting in.
	inscatteringColor *=
		(half)FresnelTransmittanceIntoWater(RefractIntoWater(L).y);

	[branch]
	if (godRays)
	{
		float4 lightScreen = mul(
			GetCamera().view_projection,
			float4(GetCamera().position + refractedLightDir * 10000, 1));
		lightScreen.xy /= lightScreen.w;

		float godray =
			GodRays(lightScreen.xy, ndc, screenUV, GetTime(), 0.1, 32.0)
			+ GodRays(lightScreen.xy, ndc, screenUV, -GetTime() * 0.5, 0.1, 20.0);

		// Blend out at the light's centre and again on the far side of it.
		godray *= 1 - pow(abs(dot(refractedLightDir, toEye)), 2);
		godray *= pow(1 - saturate(-dot(refractedLightDir, toEye)), 1);

		inscatteringColor *= 1 - godray;
	}

	const half3 fogColor = ocean.water_color.rgb + inscatteringColor;
	fog.inscatter = fogColor * inscatterAmount;

	return fog;
}

#endif // WI_WATERFOG_HF
