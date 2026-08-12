/**
 * Volumetric froxels - screen space apply pass.
 *
 * Adds the volume's light to every pixel at the depth of whatever opaque
 * surface stands there.
 *
 * **Scaffolding, not the destination.** It reproduces the shape of the pass it
 * replaces - one lookup per pixel against the opaque depth buffer - so the two
 * can be compared side by side while the volume itself is being got right. It
 * inherits that pass's defect exactly: a transparent drawn in front of the
 * surface is painted with the light from behind it, because a depth buffer
 * holding only opaque geometry cannot say a transparent was ever there. Fixing
 * that is the whole point of the volume and is the applier's job, not this
 * pass's.
 */
#include "globals.hlsli"
#include "ShaderInterop_Postprocess.h"
#include "volumetricFroxelHF.hlsli"

PUSHCONSTANT(froxels, VolumetricFroxelPush);

Texture3D<float3> input : register(t0);

RWTexture2D<float4> output : register(u0);

[numthreads(POSTPROCESS_BLOCKSIZE, POSTPROCESS_BLOCKSIZE, 1)]
void main(uint3 DTid : SV_DispatchThreadID)
{
	const float2 uv = (DTid.xy + 0.5) * GetCamera().internal_resolution_rcp;

	if (uv.x > 1 || uv.y > 1)
	{
		return;
	}

	const float depth = texture_depth.SampleLevel(sampler_point_clamp, uv, 0);
	const float3 P = reconstruct_position(uv, depth);

	// Radial distance, matching how the volume was built. Measuring view depth
	// instead would agree down the centre of the screen and drift towards the
	// edges, which is the hardest kind of error to notice and the easiest to
	// mistake for a tuning problem.
	const float distanceToEye = length(P - GetCamera().position);

	const float w = VolumetricFroxelDepthToW(distanceToEye, froxels.range);

	output[DTid.xy] += float4(
		input.SampleLevel(sampler_linear_clamp, float3(uv, w), 0), 0);
}
