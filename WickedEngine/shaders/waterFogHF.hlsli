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
 * The inscattered colour is the water's authored base plus the daylight that
 * reached this depth, from the sun and from the sky. Each is put through
 * everything that happened to it on the way down: the atmosphere it crossed,
 * the share the surface reflected away instead of admitting, and the water
 * column above the eye. The sun additionally carries the phase function and the
 * god rays, both of which need a direction the sky has not got.
 *
 * On top of what the water hands back at the wavelength it received, a little
 * comes back at a longer one - see `WaterVolumetrics::RamanEmission` and
 * `FluorescenceEmission`, which between them are most of what clear water shows
 * outside the blue.
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

	// Only the scattered share of what was extinguished comes back. This is
	// the SINGLE scattering veiling term, and it stays that way: it feeds the
	// phase reduction below, which asks how much this one path scattered.
	const float3 scatterAlbedo = medium.sigmaS / medium.sigmaT;
	const float3 inscatterAmount = scatterAlbedo * (1 - fog.transmittance);

	const half3 L = GetSunDirection();
	const half3 refractedLightDir =
		refract(-L, float3(0, 1, 0), 1.0 / WATER_REFRACTIVE_INDEX);

	// What that veil is coloured, which is a different question. Light leaves
	// turbid water having bounced many times, and the colour it settles on is
	// far more saturated than one bounce keeps - see EmergentAlbedo.
	//
	// Reshaped by where the sun is and where this is being seen from, which the
	// emergent colour on its own has no way to express: it brightens under a
	// low sun and dims towards a slanted view. Both directions are taken below
	// the surface, since that is where the light being described is.
	const float3 inscatterColorAmount = medium.EmergentAlbedo()
		* medium.BidirectionalFactor(
			RefractIntoWater(L).y, RefractIntoWater(toEye).y)
		* (1 - fog.transmittance);

	half3 sunLight = GetSunColor();

	// Apply atmosphere transmittance:
	if (GetFrame().options & OPTION_BIT_REALISTIC_SKY)
	{
		// 0 for position since fog is centered around world center
		sunLight *= GetAtmosphericLightTransmittance(
			GetWeather().atmosphere, 0, L, texture_transmittancelut);
	}

	// Sunlight has already crossed the water column above the eye before it can
	// scatter towards it, so it arrives filtered by the same medium. Red goes
	// first, which is what tints the scatter green-blue with depth without any
	// authored colour.
	//
	// Over the sun's own slant path, not the vertical drop: the sun enters at
	// whatever angle Snell's law bends it to and crosses more water than the
	// depth alone, which is the same reciprocity the view path uses.
	sunLight *= exp(-SubmergedViewPath(eyeDepth, L) * medium.sigmaT);

	// And some of it never got in: the surface reflects a share away instead of
	// refracting it down, nearly all of it once the sun is low. The same term
	// the lit surfaces use, so the haze and the sea bed agree about how much
	// sun reached the water - without it the water would still glow at sunset
	// while everything in it went dark.
	//
	// RefractIntoWater rather than refractedLightDir above: the two agree while
	// the sun is up, but refract() keeps bending a sun that has already set
	// into a downward ray, whose cosine would read as light still getting in.
	sunLight *=
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

		sunLight *= 1 - godray;
	}

	// The sky scatters in the water as well as the sun does, and unlike the sun
	// it is still there when the sun is not. Without this the water is lit by a
	// single directional source: it goes black at night, black under overcast,
	// and black in the shadow of anything - while the surfaces inside it stay
	// visibly lit, because they get the same daylight through
	// WaterAmbientTransmittance.
	//
	// Isotropic, so no phase function. A uniform field integrates the phase
	// over the whole sphere, and any phase function integrates to one - the
	// directionality the sun gets from HgPhase has nothing to weight here.
	// What reaches the water is the COSINE weighted mean of sky radiance over
	// the hemisphere above it: a patch of sky overhead delivers its full
	// radiance downwards while one at the horizon delivers none, so the zenith
	// counts for far more than an even average of the two would give it.
	//
	// The global probe's most convolved mip is that integral already, so
	// sampling it straight up asks the question directly - and asks it of the
	// same source the submerged surfaces are lit from, which is what keeps the
	// haze and the things seen through it agreeing about the daylight.
	half3 ambientColor;

	[branch]
	if (GetScene().globalprobe >= 0)
	{
		TextureCube<half4> globalProbe =
			bindless_cubemaps_half4[descriptor_index(GetScene().globalprobe)];
		uint2 probeDim;
		uint probeMips;
		globalProbe.GetDimensions(0, probeDim.x, probeDim.y, probeMips);

		ambientColor = globalProbe.SampleLevel(
			sampler_linear_clamp, float3(0, 1, 0), probeMips).rgb;
	}
	else
	{
		// Without a probe, the best available. The sky luminance LUT is an
		// average over the whole sphere rather than the hemisphere overhead,
		// so it under-reads what comes down; the authored sky can at least be
		// weighted properly, taking it as linear in \f$\cos\theta\f$ between
		// its ends:
		// \f[
		// \frac{\int_0^1 L(u)\,u\,du}{\int_0^1 u\,du}
		//   = \tfrac{1}{3} L_{horizon} + \tfrac{2}{3} L_{zenith}
		// \f]
		// which returns the sky's own radiance when it is uniform, as it must.
		ambientColor = (GetFrame().options & OPTION_BIT_REALISTIC_SKY)
			? texture_skyluminancelut.SampleLevel(
				sampler_point_clamp, float2(0.5, 0.5), 0).rgb
			: (GetHorizonColor() + GetZenithColor() * 2) / 3;
	}

	// The same two losses the sun takes, in the forms a hemispherical source
	// takes them: the mean Fresnel transmittance over the sky rather than the
	// angle-dependent one, and the cosine-weighted mean slant through Snell's
	// window rather than one refracted direction. Matches
	// WaterAmbientTransmittance exactly, so the haze and the surfaces lit
	// through it agree about how much daylight got down here.
	ambientColor *= WATER_FRESNEL_DIFFUSE_TRANSMITTANCE;
	ambientColor *= (half3)exp(
		-eyeDepth * WATER_DOWNWELLING_SLANT * medium.sigmaT);

	// Light leaving after many scattering events has no memory of the
	// direction it arrived from, so the phase function has nothing to weight -
	// the sun enters this the same way the sky does, as downwelling
	// irradiance. This term carries the water's colour.
	const half3 downwelling = ambientColor + sunLight * (half)saturate(L.y);

	// A share of that daylight comes back out redder than it went in, rather
	// than being reflected at the colour it arrived: converted by the water
	// molecules themselves, and re-emitted by the phytoplankton living in
	// them. Both put light back at wavelengths that may have been extinguished
	// long before, which is why they are what the clearest water shows outside
	// the blue and what the red owes itself to at any depth.
	const half3 inelastic = (half3)(
		medium.RamanEmission(downwelling)
		+ medium.FluorescenceEmission(downwelling));

	// One scattering event, which does remember: this is the glow that
	// sharpens towards the sun and draws the shafts. The medium's own albedo
	// weights it, not the emergent colour, because a single bounce keeps
	// exactly what a single bounce keeps.
	//
	// Only the part of that bounce which is DIRECTIONAL belongs here. Its
	// isotropic share carries no direction and is already accounted for in the
	// term above, so what is left is the excess over an even spread. That
	// excess falls away as the reduced asymmetry approaches isotropic, which
	// is the behaviour the water should have: a turbid sea has no sun glow to
	// show, only an even haze.
	const half cosTheta = dot(toEye, -refractedLightDir);
	const half phaseExcess = max(
		(half)0,
		HgPhase(medium.ReducedPhaseG(inscatterAmount), cosTheta)
			- UniformPhase());

	const half3 directional = sunLight * phaseExcess * UniformPhase();

	fog.inscatter =
		((half3)ocean.water_color.rgb + downwelling)
			* (half3)inscatterColorAmount
		+ inelastic * (half3)(1 - fog.transmittance)
		+ directional * (half3)inscatterAmount;

	return fog;
}

/**
 * The water's fog between the eye and a fragment.
 *
 * The whole point of applying this where a fragment is drawn rather than over
 * the finished frame: a post pass reads one depth per pixel, so anything
 * blended in front of geometry - a particle, a sprite, a light's visualizer -
 * gets attenuated over the distance to whatever opaque surface is BEHIND it.
 * A fragment knows how far away it actually is.
 *
 * **Both sides of the surface**, which is what makes the water absorb
 * consistently. An eye in the air looking down gets the refracted leg, bounded
 * by Snell's window; an eye in the water gets the straight segment clipped at
 * the surface; and the two are faded across as the camera crosses. Because every
 * fragment carries its own, nothing downstream has to reconstruct how much
 * water is in front of it from a depth buffer - which is impossible to do for
 * a transparent, since transparents are not in one.
 *
 * **Clipped against the WAVE surface, not the still plane.** Each end knows its
 * own signed height above the surface over it, and the water between them is
 * taken to be bounded by the straight line joining those two heights. That
 * keeps the submerged share of the segment a linear clip, exactly as a flat
 * plane made it - only tilted - so it costs a sample at each end and nothing
 * else. It is what stops an eye sitting in a trough, which is in the air but
 * below the still plane, from fogging the whole world as though submerged, and
 * an object inside a crest, which is in the water but above the still plane,
 * from being seen through as if the crest were glass.
 *
 * @note A crest BETWEEN the two ends is still missed - a straight line joining
 *       the endpoints cannot know about it. Seeing through the shoulder of a
 *       wave to the sky beyond needs the surface sampled along the segment,
 *       which is a march, not a clip.
 *
 * **Rejected cheapest first**, because this is called from every draw shader in
 * the frame. The camera flag is uniform across the draw and costs nothing in a
 * scene with no ocean; past it each end is measured against the surface, and
 * only past that exact test does anything unproject a pixel.
 *
 * There is deliberately no cheap height rejection before those measurements.
 * One would need a bound on how tall a crest can be, and there is none to have
 * - the amplitude is authored and the displacement comes off an FFT on the GPU.
 * A guessed bound does not clamp the waves, it just stops asking about them, so
 * a sea that exceeds it silently loses its fog above the guess.
 *
 * @param[in] screenUV - Screen space UV coordinates (0-1) of the fragment.
 * @param[in] fragmentPosition - World position of the fragment.
 * @param[in] fragmentHeight - Signed height of the fragment above the water
 *                             surface over it (in metres): negative under the
 *                             water, positive in the air, zero on the interface
 *                             itself.
 *
 * @return The fog over that segment. Transparent and unlit when the water does
 *         not apply, so callers need no test of their own.
 */
WaterFog GetWaterFog(
	float2 screenUV, float3 fragmentPosition, float fragmentHeight
)
{
	WaterFog fog;
	fog.transmittance = 1;
	fog.inscatter = 0;

	[branch]
	if (!GetCamera().IsWaterFog())
	{
		return fog;
	}

	const float3 eye = GetCamera().position;

	// The same quantity the caller supplied for the far end, for the near one.
	const float eyeHeight = eye.y - ocean_drawn_surface_height(eye);

	// Nothing on this segment is under water, so there is nothing to fog. Exact
	// now that both ends have been measured against the surface rather than the
	// plane: the submerged share below is clipped between these same two
	// heights, so this is the case where it would come out zero.
	[branch]
	if (min(eyeHeight, fragmentHeight) >= 0)
	{
		return fog;
	}

	// How submerged the EYE is at this pixel: 0 in the air, 1 under water, and
	// graded in between so the waterline sweeps the screen instead of snapping
	// as the camera crosses. Its height settles where the camera's centre is but
	// not where each pixel is, and this is the expensive way to ask - it
	// unprojects the pixel - so it is only asked when the answer is in doubt.
	const half eyeSubmersion =
		eyeHeight < WATER_EYE_SUBMERSION_MARGIN
			? (half)ocean_underwater_factor(screenUV)
			: (half)0;

	// The water is there if EITHER end of the segment is in it. Taking the
	// larger is what admits the case this whole function used to miss: an eye
	// in the air looking at something under water. The medium is fully present
	// along that segment even though the eye is not in it.
	const WaterVolumetrics medium = MakeWaterVolumetrics(
		max(eyeSubmersion, fragmentHeight < 0 ? (half)1 : (half)0));

	[branch]
	if (!medium.IsActive())
	{
		return fog;
	}

	const float3 towardsEye = eye - fragmentPosition;
	const float segment = length(towardsEye);
	const float3 toEye = towardsEye / max(segment, 0.00001);

	// The straight segment, clipped at the surface. Right whenever the eye is in
	// the water: between two points in one medium nothing bends. Both ends on
	// the same side saturate to the whole segment or to none of it, so only a
	// segment that actually crosses divides.
	const float straightPath = segment * saturate(
		-min(eyeHeight, fragmentHeight)
		/ max(abs(fragmentHeight - eyeHeight), 0.00001));

	// The refracted leg, for an eye in the air. This is not a refinement of the
	// straight path, it REPLACES it - light leaving the water bends towards the
	// horizontal, so a grazing view crosses far less water than the line of
	// sight suggests. An eye two metres up looking at a sea bed one metre down
	// three hundred metres out has a hundred metres of straight segment below
	// the plane and about one and a half metres of actual water. Using the
	// segment would extinguish water you can plainly see through.
	const float refractedPath =
		SubmergedViewPath(max(0, -fragmentHeight), toEye);

	// Blend on the EYE's submersion, not the medium's: the medium is fully
	// present in both cases and cannot tell them apart. The two formulas differ
	// by a lot at the crossing - which is genuine, and why the surface is a
	// mirror from just above and a window from just below - so they are faded
	// across rather than switched, on the same graded test that sweeps the
	// waterline down the screen.
	const float path = lerp(refractedPath, straightPath, eyeSubmersion);

	return MakeWaterFog(
		medium,
		path,
		toEye,
		max(0, -eyeHeight),
		screenUV,
		uv_to_clipspace(screenUV),
		// Screen space stripes swept around where the sun projects. They read
		// as shafts seen from inside the water and as stripes painted over the
		// sea seen from a boat, so they need the eye to be under, not merely
		// the water to be present.
		GetCamera().IsUnderwaterGodRays() && eyeSubmersion > 0
	);
}

/**
 * The water's fog between the eye and a fragment somewhere in the scene.
 *
 * Measures the fragment against the wave surface over it, which is what every
 * ordinary draw wants. The camera flag is tested here as well as inside, so a
 * scene with no ocean never pays the displacement map sample that measuring
 * would take.
 *
 * @param[in] screenUV - Screen space UV coordinates (0-1) of the fragment.
 * @param[in] fragmentPosition - World position of the fragment.
 *
 * @return The fog over that segment.
 */
WaterFog GetWaterFog(float2 screenUV, float3 fragmentPosition)
{
	WaterFog fog;
	fog.transmittance = 1;
	fog.inscatter = 0;

	[branch]
	if (!GetCamera().IsWaterFog())
	{
		return fog;
	}

	return GetWaterFog(
		screenUV,
		fragmentPosition,
		fragmentPosition.y - ocean_drawn_surface_height(fragmentPosition)
	);
}

/**
 * Applies an already-built fog to a fragment, passing a background through.
 *
 * The background excludes radiance that was sampled from the scene behind this
 * fragment - a refraction, most often - because that was already fogged over
 * its own longer path when it was drawn. Fogging it a second time would darken
 * a pane of glass by the water on both sides of it.
 *
 * @param[in] fog - The fog over the segment.
 * @param[in] background - Radiance already fogged elsewhere, to pass through
 *                         untouched.
 * @param[in,out] color - Fragment colour, fogged in place.
 */
void ApplyWaterFog(WaterFog fog, half3 background, inout half4 color)
{
	// Never take out more than is there. A shader that overwrites its colour
	// after the refraction was composited - UNLIT, interior mapping, forced
	// unlit - would otherwise drive the fogged term negative.
	const half3 alreadyFogged = min(background, color.rgb);

	color.rgb = (half3)(
		(color.rgb - alreadyFogged) * fog.transmittance
		+ fog.inscatter
		+ alreadyFogged
	);
}

/**
 * Fogs a fragment by the water between it and the eye.
 *
 * Example usage:
 * @code
 * const half3 background =
 *     surface.refraction.rgb * (1 - surface.F) * surface.refraction.a;
 * ApplyWaterFog(ScreenCoord, surface.P, background, color);
 * @endcode
 *
 * @param[in] screenUV - Screen space UV coordinates (0-1) of the fragment.
 * @param[in] fragmentPosition - World position of the fragment.
 * @param[in] background - Radiance already fogged elsewhere, to pass through
 *                         untouched.
 * @param[in,out] color - Fragment colour, fogged in place.
 */
void ApplyWaterFog(
	float2 screenUV,
	float3 fragmentPosition,
	half3 background,
	inout half4 color
)
{
	ApplyWaterFog(GetWaterFog(screenUV, fragmentPosition), background, color);
}

/**
 * Fogs a fragment that lies ON the water surface.
 *
 * The interface is never submerged with respect to itself: its height above the
 * water surface is zero by definition, so what separates it from the eye is
 * water exactly when the eye is under, and air when the eye is not.
 *
 * Said rather than measured, because measuring would be slightly wrong. The
 * surface is sampled at the position's xz, but choppiness moves a point
 * horizontally as well as vertically, so a crest's fragment does not sit over
 * the part of the map that produced it. The answer would come back near zero
 * and drift either side of it, and the drift is what shows - as banding along
 * the wave shape, trough by trough.
 *
 * @param[in] screenUV - Screen space UV coordinates (0-1) of the fragment.
 * @param[in] fragmentPosition - World position of the fragment.
 * @param[in] background - Radiance already fogged elsewhere, to pass through
 *                         untouched.
 * @param[in,out] color - Fragment colour, fogged in place.
 */
void ApplyWaterFogAtSurface(
	float2 screenUV,
	float3 fragmentPosition,
	half3 background,
	inout half4 color
)
{
	ApplyWaterFog(
		GetWaterFog(screenUV, fragmentPosition, 0), background, color
	);
}

/**
 * Fogs a fragment that carries no already-fogged background.
 *
 * @param[in] screenUV - Screen space UV coordinates (0-1) of the fragment.
 * @param[in] fragmentPosition - World position of the fragment.
 * @param[in,out] color - Fragment colour, fogged in place.
 */
void ApplyWaterFog(
	float2 screenUV,
	float3 fragmentPosition,
	inout half4 color
)
{
	ApplyWaterFog(screenUV, fragmentPosition, 0, color);
}

/**
 * Applies an already-built fog to a premultiplied fragment.
 *
 * Fonts, gaussian splats and clouds output premultiplied radiance, so the
 * veiling light has to be scaled the same way. Adding the full inscatter to a
 * fragment that only covers a tenth of the pixel would paint haze into the nine
 * tenths it never touched.
 *
 * @param[in] fog - The fog over the segment.
 * @param[in,out] color - Premultiplied fragment colour, fogged in place.
 */
void ApplyWaterFogPremultiplied(WaterFog fog, inout half4 color)
{
	color.rgb = (half3)(
		color.rgb * fog.transmittance + fog.inscatter * color.a
	);
}

/**
 * Fogs a premultiplied fragment by the water between it and the eye.
 *
 * @param[in] screenUV - Screen space UV coordinates (0-1) of the fragment.
 * @param[in] fragmentPosition - World position of the fragment.
 * @param[in,out] color - Premultiplied fragment colour, fogged in place.
 */
void ApplyWaterFogPremultiplied(
	float2 screenUV,
	float3 fragmentPosition,
	inout half4 color
)
{
	ApplyWaterFogPremultiplied(GetWaterFog(screenUV, fragmentPosition), color);
}

/**
 * Fogs a fragment that is added to what is already on screen.
 *
 * Transmittance only, deliberately: the destination has already had the water's
 * veiling light put into it by whatever drew it. An additive draw contributes
 * nothing of its own to that haze, and adding it again would multiply the haze
 * by the number of additive draws stacked on the pixel - so a lamp's glow and
 * its visualizer would each thicken the water.
 *
 * @param[in] screenUV - Screen space UV coordinates (0-1) of the fragment.
 * @param[in] fragmentPosition - World position of the fragment.
 * @param[in,out] color - Fragment colour, attenuated in place.
 */
void ApplyWaterFogAdditive(
	float2 screenUV,
	float3 fragmentPosition,
	inout half4 color
)
{
	color.rgb = (half3)(
		color.rgb * GetWaterFog(screenUV, fragmentPosition).transmittance
	);
}

#endif // WI_WATERFOG_HF
