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
 * **Both ends are measured against the drawn waves.** The plane reads a crest
 * as standing in air and a trough as drowned, and it is wrong about the eye in
 * exactly the same way: an eye sitting in a trough is below the plane and
 * plainly dry, and taking it for submerged invents a wet leg that is not there,
 * cutting the air fog short over the whole view - which under a realistic sky
 * removes a bright haze and reads as the far water going dark.
 *
 * The plane still rejects a dry camera cheaply, offset by the tallest wave the
 * sea can raise, so the displacement sample is paid for only inside the band
 * where the answer is genuinely in doubt.
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
 *       still disagree while the camera straddles the surface, though now over
 *       the wave it is actually in rather than over the whole swell.
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

	// Above every crest the sea can raise, nothing needs measuring: the eye is
	// dry and the whole segment is air. Costs no sample, which is what keeps
	// this affordable in every draw shader that fogs anything.
	[branch]
	if (O.y >= GetWeather().ocean.water_height + ocean_max_displacement())
	{
		return float2(0, distance);
	}

	// Inside that band the plane cannot answer and the wave over the eye has to
	// be measured, against the same drawn surface the far end is measured
	// against.
	const float originAbove = O.y - ocean_drawn_surface_height(O);

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
// Analytic optical depth of the air over a segment, which is the quantity the
// closed form below is really built out of. Exposed on its own because a caller
// integrating a column in pieces needs each piece's own optical depth: sampling
// an extinction at one point and calling it constant over the piece is exactly
// wrong for height fog, whose density falls off exponentially, and a piece
// spanning any real height is then out by however far the density moved across
// it.
//
//	distance			: length of the segment
//	O					: segment start
//	V					: segment direction, normalized
//	nearFadeDistance	: distance the artistic near fade is judged at, which
//						  for a piece of a column is the piece's own middle
//						  rather than its length
inline float GetFogOpticalDepth(
	float distance, float3 O, float3 V, float nearFadeDistance
)
{
	ShaderFog fog = GetWeather().fog;
	
	// The near fade keeps the fragment's own distance rather than the air's:
	// it is authored against how far away a thing is, not against how much air
	// stands in front of it.
	float startDistanceFalloff =
		saturate((nearFadeDistance - fog.start) / fog.start);

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
		
		return fog.density * startDistanceFalloff * (baseHeightFogDistance + exponentialHeightLineIntegral);
	}
	else
	{
		// Height fog algorithm (above) reduced with infinity start and end heights:
		
		return fog.density * startDistanceFalloff * airDistance;
	}
}

//	distance	: sample to point distance
//	O			: sample position
//	V			: sample to point vector
inline half GetFogAmount(float distance, float3 O, float3 V)
{
	return (half)(1.0 - exp(-GetFogOpticalDepth(distance, O, V, distance)));
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
