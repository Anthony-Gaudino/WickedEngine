#ifndef WI_VOLUMETRIC_CLOUD_APPLY_HF
#define WI_VOLUMETRIC_CLOUD_APPLY_HF

/**
 * The clouds standing in front of a fragment, applied where it is drawn.
 *
 * The cloud march resolves to one screen space buffer, so the frame can only
 * composite it at a single moment - after the opaque pass, before the
 * transparents. Everything drawn later paints straight over it: an ocean or a
 * transmissive mesh stays plainly visible from inside a cloud that has already
 * swallowed the land around it.
 *
 * A fragment that asks for itself is not bound to that moment. This is the same
 * answer `ApplyVolumetricLight` is built on, and for the same reason.
 *
 * @note **All or nothing, by construction.** The buffer holds one integrated
 *       result and one representative distance for the whole march, so it can
 *       say whether the cloud is in front of a fragment but not how much of it
 *       is. A fragment standing INSIDE the layer takes all of the cloud or none
 *       of it. Answering that properly needs the march to publish a running
 *       total the way the froxel volume does.
 *
 * @note **A refracted share receives the clouds twice**: once in the scene copy
 *       it was sampled from, which is taken after they are composited, and once
 *       here. It shows as a depth edge behind the water drawing itself while
 *       the layer is being entered or left. Weighing it out takes the share the
 *       copy actually supplied - `refraction.a` alone is not that, because a
 *       shader standing its own content in where the copy could not answer
 *       still reports it as fully refracted.
 */

#include "globals.hlsli"

/**
 * Blends the clouds standing in front of a fragment over it.
 *
 * The buffer is premultiplied - radiance already scaled by its own coverage,
 * with that coverage in alpha - so it composites as `over` and never needs the
 * colour dividing back out.
 *
 * Sampled at half resolution with a single bilinear tap rather than the
 * bilateral upsample the composite pass uses. That upsample exists to keep a
 * cloud edge off the geometry behind it, which is a hard edge on an opaque
 * surface; here it is being laid over something transparent and soft, and nine
 * taps at every such fragment would cost more than the halo it avoids.
 *
 * Example usage:
 * @code
 * ApplyVolumetricClouds(ScreenCoord, surface.P, color);
 * @endcode
 *
 * @param[in] screenUV - Screen space UV coordinates (0-1) of the fragment.
 * @param[in] fragmentPosition - World position of the fragment.
 * @param[in,out] color - Radiance to blend the clouds over.
 *
 * @note Does nothing for a camera without clouds, which reads -1 and skips.
 */
inline void ApplyVolumetricClouds(
	in float2 screenUV, in float3 fragmentPosition, inout half4 color
)
{
	[branch]
	if (GetCamera().texture_volumetricclouds_index < 0)
	{
		return;
	}

	// **The distance the march settled on**, which is `tDepth` - the depth
	// weighted mean of where the cloud actually is along this ray. `.g` is the
	// distance to whatever the depth buffer held when the clouds were marched,
	// and is the bilateral upsample's business rather than this one's.
	const float cloudDistance = texture_volumetricclouds_depth.SampleLevel(
		sampler_linear_clamp, screenUV, 0).r;

	// Behind this fragment, so the fragment covers it and there is nothing to
	// lay over. `FLT_MAX` where the march found no cloud at all falls out here
	// without a test of its own.
	[branch]
	if (cloudDistance >= distance(fragmentPosition, GetCamera().position))
	{
		return;
	}

	const half4 cloud = texture_volumetricclouds.SampleLevel(
		sampler_linear_clamp, screenUV, 0);

	color.rgb = mad(color.rgb, 1 - cloud.a, cloud.rgb);
}

/**
 * The same, for a draw whose colour is already scaled by its own coverage.
 *
 * A premultiplied draw carries its coverage in alpha and its blend does not
 * scale what it writes, so the cloud laid over it has to be scaled by that
 * coverage here - otherwise a fragment covering a tenth of a pixel hands over a
 * whole pixel of cloud.
 *
 * @param[in] screenUV - Screen space UV coordinates (0-1) of the fragment.
 * @param[in] fragmentPosition - World position of the fragment.
 * @param[in,out] color - Premultiplied radiance to blend the clouds over.
 */
inline void ApplyVolumetricCloudsPremultiplied(
	in float2 screenUV, in float3 fragmentPosition, inout half4 color
)
{
	half4 opaqueColor = half4(color.rgb, 1);
	ApplyVolumetricClouds(screenUV, fragmentPosition, opaqueColor);

	color.rgb = lerp(color.rgb, opaqueColor.rgb, color.a);
}

#endif // WI_VOLUMETRIC_CLOUD_APPLY_HF
