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

	// Phase function of the medium, reduced by how much of the light this
	// particular path actually scattered - see ReducedPhaseG. There is exactly
	// one path length here, unlike in a march, so the medium's own albedo
	// reduced phaseG would be throwing that away.
	const half cosTheta = dot(toEye, -refractedLightDir);
	inscatteringColor *= HgPhase(medium.ReducedPhaseG(inscatterAmount), cosTheta);

	// Uniform base phase, matching what the engine's height fog does for a
	// constant medium (see GetFog in fogHF.hlsli):
	inscatteringColor *= UniformPhase();

	// Sunlight has already crossed the water column above the eye before it can
	// scatter towards it, so it arrives filtered by the same medium. Red goes
	// first, which is what tints the scatter green-blue with depth without any
	// authored colour.
	//
	// Over the sun's own slant path, not the vertical drop: the sun enters at
	// whatever angle Snell's law bends it to and crosses more water than the
	// depth alone, which is the same reciprocity the view path uses.
	inscatteringColor *= exp(-SubmergedViewPath(eyeDepth, L) * medium.sigmaT);

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
	half3 ambientColor = (GetFrame().options & OPTION_BIT_REALISTIC_SKY)
		? texture_skyluminancelut.SampleLevel(
			sampler_point_clamp, float2(0.5, 0.5), 0).rgb
		: (GetHorizonColor() + GetZenithColor()) * 0.5;

	// The same two losses the sun takes, in the forms a hemispherical source
	// takes them: the mean Fresnel transmittance over the sky rather than the
	// angle-dependent one, and the cosine-weighted mean slant through Snell's
	// window rather than one refracted direction. Matches
	// WaterAmbientTransmittance exactly, so the haze and the surfaces lit
	// through it agree about how much daylight got down here.
	ambientColor *= WATER_FRESNEL_DIFFUSE_TRANSMITTANCE;
	ambientColor *= (half3)exp(
		-eyeDepth * WATER_DOWNWELLING_SLANT * medium.sigmaT);

	const half3 fogColor =
		ocean.water_color.rgb + ambientColor + inscatteringColor;
	fog.inscatter = fogColor * inscatterAmount;

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
 * the plane; and the two are faded across as the camera crosses. Because every
 * fragment carries its own, nothing downstream has to reconstruct how much
 * water is in front of it from a depth buffer - which is impossible to do for
 * a transparent, since transparents are not in one.
 *
 * Clipped against the still water plane, matching `SubmergedLightPath`, so the
 * fog agrees with the lighting it fogs.
 *
 * **Rejected cheapest first**, because this is called from every draw shader in
 * the frame. The camera option is uniform across the draw; the still plane
 * test after it is two scalar compares and is exact, since a segment with both
 * ends above the plane has no submerged share to fog whatever the waves are
 * doing; and only past those does anything unproject a pixel or sample the
 * displacement map to find out where the eye sits.
 *
 * @param[in] screenUV - Screen space UV coordinates (0-1) of the fragment.
 * @param[in] fragmentPosition - World position of the fragment.
 *
 * @return The fog over that segment. Transparent and unlit when the water does
 *         not apply, so callers need no test of their own.
 */
WaterFog GetWaterFog(
	float2 screenUV, float3 fragmentPosition, half fragmentSubmersion
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
	const float waterHeight = GetWeather().ocean.water_height;

	// Nothing on this segment is under water, so there is nothing to fog. Exact
	// rather than conservative: the submerged share below is clipped against
	// this same plane, so this is the case where it would come out zero.
	[branch]
	if (min(eye.y, fragmentPosition.y) >= waterHeight)
	{
		return fog;
	}

	// How submerged the EYE is at this pixel: 0 in the air, 1 under water, and
	// graded in between so the waterline sweeps the screen instead of snapping
	// as the camera crosses. This is the expensive question - it unprojects the
	// pixel and samples the displacement map - so it is only asked when the eye
	// is near enough to the surface for the answer to be in doubt.
	const half eyeSubmersion =
		eye.y < waterHeight + WATER_EYE_SUBMERSION_MARGIN
			? (half)ocean_underwater_factor(screenUV)
			: (half)0;

	// The water is there if EITHER end of the segment is in it. Taking the
	// larger is what admits the case this whole function used to miss: an eye
	// in the air looking at something under water. The medium is fully present
	// along that segment even though the eye is not in it.
	const WaterVolumetrics medium =
		MakeWaterVolumetrics(max(eyeSubmersion, fragmentSubmersion));

	[branch]
	if (!medium.IsActive())
	{
		return fog;
	}

	const float3 towardsEye = eye - fragmentPosition;
	const float segment = length(towardsEye);
	const float3 toEye = towardsEye / max(segment, 0.00001);

	// The straight segment, clipped at the still plane. Right whenever the eye
	// is in the water: between two points in one medium nothing bends.
	const float straightPath = segment * saturate(
		(waterHeight - min(eye.y, fragmentPosition.y))
		/ max(abs(fragmentPosition.y - eye.y), 0.00001));

	// The refracted leg, for an eye in the air. This is not a refinement of the
	// straight path, it REPLACES it - light leaving the water bends towards the
	// horizontal, so a grazing view crosses far less water than the line of
	// sight suggests. An eye two metres up looking at a sea bed one metre down
	// three hundred metres out has a hundred metres of straight segment below
	// the plane and about one and a half metres of actual water. Using the
	// segment would extinguish water you can plainly see through.
	const float refractedPath =
		SubmergedViewPath(max(0, waterHeight - fragmentPosition.y), toEye);

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
		max(0, waterHeight - eye.y),
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
 * Takes the fragment to be in the water whenever it is below the still plane,
 * which is what every ordinary draw wants.
 *
 * @param[in] screenUV - Screen space UV coordinates (0-1) of the fragment.
 * @param[in] fragmentPosition - World position of the fragment.
 *
 * @return The fog over that segment.
 */
WaterFog GetWaterFog(float2 screenUV, float3 fragmentPosition)
{
	return GetWaterFog(
		screenUV,
		fragmentPosition,
		fragmentPosition.y < GetWeather().ocean.water_height
			? (half)1
			: (half)0
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
 * The interface is never submerged with respect to itself: what separates it
 * from the eye is water exactly when the eye is under, and air when the eye is
 * not. Every other fragment can be classified by its height, but this one sits
 * on the **displaced** surface, which dips below the still plane in a trough -
 * so the ordinary test would have it fog itself over a sliver of water that is
 * not there, trough by trough, as banding along the wave shape.
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
