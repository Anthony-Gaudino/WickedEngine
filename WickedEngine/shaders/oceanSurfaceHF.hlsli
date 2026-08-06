#ifndef WI_OCEAN_SURFACE_HF
#define WI_OCEAN_SURFACE_HF
#include "globals.hlsli"
#include "underwaterHF.hlsli"
#include "ShaderInterop_Ocean.h"

//static const float OCEAN_NEARPLANE_CUTOFF = 0.1;
#define OCEAN_NEARPLANE_CUTOFF compute_inverse_lineardepth(max(GetCamera().z_near + 1, 1.0))

/**
 * Water a wave presents to light coming from behind it, in metres.
 *
 * `x` is flat water, `y` a fully folding crest. The medium extinguishes over
 * roughly a metre at the default turbidity, so the first number has to be deep
 * enough that open water stays dark and only the crests light up - otherwise
 * the whole sea glows and the effect reads as a wash rather than as thin water.
 *
 * **This mapping is calibrated, not derived**, and it is the only part of the
 * transmission that is. Everything else - the colour, the falloff, the phase,
 * both Fresnel interfaces - comes from the water medium, which is why the
 * strength control defaults to 1: that is the unmodified physical result.
 *
 * A derived thickness is not available. The Jacobian fold is a dimensionless
 * measure of how far the surface has compressed, with no length in it, and on
 * a deep ocean there is no thin slab to measure anyway - the glow comes from
 * the lip of a steep or folding crest, which is exactly what the fold detects
 * and exactly what a wave height would miss.
 */
static const float2 OCEAN_SUBSURFACE_THICKNESS = float2(4.0, 0.2);

/**
 * Distance band over which the FFT gradient gives way to the perlin
 * substitute, in metres.
 *
 * Shading only, and far longer than the geometric band
 * `ocean_displacement_fade` covers - the two used to share one curve, which
 * cost real wave shading to buy nothing.
 *
 * Note the two bands now overlap: the perlin starts taking over here while the
 * displacement is still being drawn, so between them the geometry repeats the
 * FFT patch while the shading has already begun hiding that repetition. That
 * is the seam to look at first if tiling shows up at range.
 *
 * The perlin octaves are all LOWER frequency than the base patch, so their job
 * is to hide the patch tiling at range, not to supply detail. Handing over
 * early therefore discards genuine FFT slopes and replaces them with something
 * coarser. That trade only made sense while those slopes aliased once they went
 * sub-pixel; the slope variance folded into roughness handles that properly
 * now, so the gradient can survive an order of magnitude further and the perlin
 * takes over only where tiling would actually start to show.
 */
static const float2 OCEAN_GRADIENT_FADE = float2(160, 1200);

struct PSIn
{
	float4 pos : SV_POSITION;
	float2 uv : TEXCOORD0;
	min16uint cameraIndex : CAMERAINDEX;
	
	inline float3 GetPos3D()
	{
		return GetCameraIndexed(cameraIndex).screen_to_world(pos);
	}

	inline float3 GetViewVector()
	{
		return GetCameraIndexed(cameraIndex).screen_to_view(pos);
	}
};

float intersectPlaneClampInfiniteDist(in float3 rayOrigin, in float3 rayDirection, in float3 planeNormal, float planeHeight)
{
	return (planeHeight - dot(planeNormal, rayOrigin)) / dot(planeNormal, rayDirection);
}
float3 intersectPlaneClampInfinite(in float3 rayOrigin, in float3 rayDirection, in float3 planeNormal, float planeHeight)
{
	float dist = intersectPlaneClampInfiniteDist(rayOrigin, rayDirection, planeNormal, planeHeight);
	if (dist > 0.0)
		return rayOrigin + rayDirection * dist;

	// The ray points away from the plane, so push it out to the horizon along
	// its horizontal heading instead. A ray aimed straight up or straight down
	// has no horizontal heading to normalize - that is normalize(0), which is
	// NaN, and a NaN vertex takes every triangle sharing it with it. Fall back
	// to an arbitrary but finite heading; the point is z_far away either way,
	// so which direction it leaves in is not observable.
	const float2 horizontal = rayDirection.xz;
	const float horizontalLength = length(horizontal);
	const float2 heading = horizontalLength > 1e-6 ? horizontal / horizontalLength : float2(0, 1);

	return float3(rayOrigin.x, planeHeight, rayOrigin.z) + float3(heading.x, 0, heading.y) * GetCamera().z_far;
}

#endif // WI_OCEAN_SURFACE_HF
