#ifndef WI_FOG_HF
#define WI_FOG_HF
#include "globals.hlsli"
#include "skyAtmosphere.hlsli"
#include "underwaterHF.hlsli"
#include "volumetricFroxelHF.hlsli"

// [-0.999; 0.999] Describes how the lighting is destributed across sky
#define FOG_INSCATTERING_PHASE_G 0.6

/**
 * The stretch of a view segment that is filled with air, in metres.
 *
 * Height fog is a description of the air, and an eye under the water has none
 * of it in front of anything it can see. Without this the whole submerged view
 * is hazed towards the fog colour over its full length, which is the one place
 * the fog is applied to a medium that is not there.
 *
 * One start and one length describe the result: the segment is a straight line
 * and it is judged by its two ends alone, so the air it crosses comes out as a
 * single contiguous stretch.
 *
 * **An eye in the air is left alone**, and the straight line's dip below the
 * still plane must NOT be taken for water. Snell confines the leg a submerged
 * thing's light travels in the water to the critical angle, so it left the
 * water within about 1.5 times its own depth - the rest of the way to a dry eye
 * is air, however far below the plane the line of sight passes. Cutting the fog
 * at the plane would take hundreds of metres of that air at a grazing angle,
 * and take it as a function of the wave height at the far end: the end of the
 * air stretch moves by `distance / eyeHeight` metres per metre of wave, which
 * paints the swell onto the fog as hard-edged rings. What is genuinely
 * over-fogged is the wet leg itself, and it is bounded, and wherever it is long
 * enough to matter the water fog has already extinguished the fragment - so the
 * air counted with it cannot be seen.
 *
 * **The two ends are measured against different surfaces, deliberately.** The
 * far end is measured against the drawn waves, because that end is often a
 * fragment ON the surface - the ocean's own - and the plane would read a crest
 * as standing in air and a trough as drowned, which is the swell painted into
 * the water's underside. The near end is measured against the plane, because it
 * is what rejects a dry camera, and rejecting it there is what keeps the
 * displacement sample off every fragment of every scene whose camera is not in
 * the water. A camera inside a crest but above the plane therefore counts as
 * dry, which is what the rest of the engine already makes of it.
 *
 * @param[in] distance - Length of the segment (in metres).
 * @param[in] O - World position the segment starts at.
 * @param[in] V - Normalized direction the segment runs in.
 *
 * @return Distance from `O` to where the air begins in `x`, and the length of
 *         the air stretch in `y`, both in metres. A zero length means the whole
 *         segment is under water.
 *
 * @note Switches on the eye's own height, for the whole screen at once, where
 *       the water's fog sweeps the waterline down the screen on a per-pixel
 *       test it cannot share (this has no pixel to test). The two therefore
 *       disagree for as long as the camera straddles the surface.
 */
inline float2 GetFogAirSegment(float distance, float3 O, float3 V)
{
	// Nothing more than "does this scene have an ocean". With no water there is
	// no surface to clip against and the whole segment is air.
	[branch]
	if (!GetCamera().IsWaterFog())
	{
		return float2(0, distance);
	}

	const float originAbove = O.y - GetWeather().ocean.water_height;

	[branch]
	if (originAbove >= 0)
	{
		return float2(0, distance);
	}

	const float3 end = mad(distance, V, O);
	const float endAbove = end.y - ocean_drawn_surface_height(end);

	[branch]
	if (endAbove <= 0)
	{
		return float2(0, 0);
	}

	// The two ends straddle the surface here, so the difference below cannot
	// vanish and the saturate only guards the arithmetic.
	const float crossing =
		distance * saturate(-originAbove / (endAbove - originAbove));

	return float2(crossing, distance - crossing);
}

// Exponential height fog based on: https://www.iquilezles.org/www/articles/fog/fog.htm
// Non constant density function
//	distance	: sample to point distance
//	O			: sample position
//	V			: sample to point vector
inline half GetFogAmount(float distance, float3 O, float3 V)
{
	ShaderFog fog = GetWeather().fog;
	
	// The near fade keeps the fragment's own distance rather than the air's:
	// it is authored against how far away a thing is, not against how much air
	// stands in front of it.
	float startDistanceFalloff = saturate((distance - fog.start) / fog.start);

	const float2 air = GetFogAirSegment(distance, O, V);

	[branch]
	if (air.y <= 0)
	{
		return 0;
	}

	const float3 airOrigin = mad(air.x, V, O);
	const float airDistance = air.y;

	if (GetFrame().options & OPTION_BIT_HEIGHT_FOG)
	{
		float fogFalloffScale = rcp(max(0.01, fog.height_end - fog.height_start));

		// solve for x, e^(-h * x) = 0.001
		// x = 6.907755 * h^-1
		float fogFalloff = 6.907755 * fogFalloffScale;
		
		float originHeight = airOrigin.y;
		float Z = V.y;
		float effectiveZ = max(abs(Z), 0.001);

		float endLineHeight = mad(airDistance, Z, originHeight); // Isolated vector equation for y
		float minLineHeight = min(originHeight, endLineHeight);
		float heightLineFalloff = max(minLineHeight - fog.height_start, 0);
		
		float baseHeightFogDistance = clamp((fog.height_start - minLineHeight) / effectiveZ, 0, airDistance);
		float exponentialFogDistance = airDistance - baseHeightFogDistance; // Exclude distance below base height
		float exponentialHeightLineIntegral = exp(-heightLineFalloff * fogFalloff) * (1.0 - exp(-exponentialFogDistance * effectiveZ * fogFalloff)) / (effectiveZ * fogFalloff);
		
		float opticalDepthHeightFog = fog.density * startDistanceFalloff * (baseHeightFogDistance + exponentialHeightLineIntegral);
		float transmittanceHeightFog = exp(-opticalDepthHeightFog);
		
		float fogAmount = transmittanceHeightFog;
		return 1.0 - fogAmount;
	}
	else
	{
		// Height fog algorithm (above) reduced with infinity start and end heights:
		
		float opticalDepthHeightFog = fog.density * startDistanceFalloff * airDistance;
		float transmittanceHeightFog = exp(-opticalDepthHeightFog);
		
		float fogAmount = transmittanceHeightFog;
		return 1.0 - fogAmount;
	}
}

inline half4 GetFog(float distance, float3 O, float3 V)
{
	half3 fogColor = 0;
	
	if ((GetFrame().options & OPTION_BIT_REALISTIC_SKY) && (GetFrame().options & OPTION_BIT_OVERRIDE_FOG_COLOR) == 0)
	{
		// Sample captured ambient color from realisitc sky:
		fogColor = texture_skyluminancelut.SampleLevel(sampler_point_clamp, float2(0.5, 0.5), 0).rgb;
	}
	else
	{
		fogColor = GetHorizonColor();
	}

	// Sample inscattering color:
	//
	// Left to the froxel volume where it carries the sun, because that version
	// is SHADOWED and this one is not - the two together brighten the fog twice
	// over and put the second helping inside every shaft, where there is no
	// sunlight to scatter at all.
	[branch]
	if (!VolumetricFroxelCarriesTheSun())
	{
		const half3 L = GetSunDirection();
		
		half3 inscatteringColor = GetSunColor();

		// Apply atmosphere transmittance:
		if (GetFrame().options & OPTION_BIT_REALISTIC_SKY)
		{
			// 0 for position since fog is centered around world center
			inscatteringColor *= GetAtmosphericLightTransmittance(GetWeather().atmosphere, 0, L, texture_transmittancelut);
		}
		
		// Apply phase function solely for directionality:
		const half cosTheta = dot(-V, L);
		inscatteringColor *= HgPhase(FOG_INSCATTERING_PHASE_G, cosTheta);

		// Apply uniform phase since this medium is constant:
		inscatteringColor *= UniformPhase();
		
		fogColor += inscatteringColor;
	}
	
	return half4(fogColor, GetFogAmount(distance, O, V));
}


#endif // WI_FOG_HF
