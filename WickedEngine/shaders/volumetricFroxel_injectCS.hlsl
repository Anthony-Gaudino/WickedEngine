/**
 * Volumetric froxels - injection pass.
 *
 * Fills every cell with the light the medium there scatters towards the eye,
 * **per metre of path**. Nothing is accumulated here: a cell's value depends
 * only on where it is and what shines on it, which is what makes it a local
 * quantity - and so the half of the work that a later frame can reuse.
 *
 * The integration pass turns these into what a fragment reads.
 */
#include "globals.hlsli"
#include "lightingHF.hlsli"
#include "fogHF.hlsli"
#include "volumetricFroxelLightingHF.hlsli"

// One comparison tap per cell rather than a filtered kernel. The volume is
// already a spatial average over a cell far larger than a shadow texel, so a
// soft lookup would blur something that is about to be blurred anyway.
#define DISABLE_SOFT_SHADOWMAP

PUSHCONSTANT(froxels, VolumetricFroxelPush);

RWTexture3D<float3> output : register(u0);

[numthreads(
	VOLUMETRIC_FROXEL_BLOCKSIZE_X,
	VOLUMETRIC_FROXEL_BLOCKSIZE_Y,
	VOLUMETRIC_FROXEL_BLOCKSIZE_Z)]
void main(uint3 DTid : SV_DispatchThreadID)
{
	// Offset within the cell rather than its centre, on all three axes. The
	// depth axis is what breaks the slice boundaries up; the screen axes stop
	// the shadow test landing on the same grid the volume is stored on, which
	// would otherwise print that grid along the edge of every shaft.
	const float3 jitter = VolumetricFroxelCellJitter(DTid);

	const float2 uv = (DTid.xy + jitter.xy) / float2(
		VOLUMETRIC_FROXEL_WIDTH, VOLUMETRIC_FROXEL_HEIGHT);

	const float depth =
		VolumetricFroxelSliceToDepth((float)DTid.z + jitter.z, froxels.range);
	const float3 P = VolumetricFroxelPosition(uv, depth);
	const float3 toEye = normalize(GetCamera().position - P);

	const float sigmaT = VolumetricFroxelAirExtinction(P, depth);

	// Non-absorbing air: everything the medium takes out of a beam it scatters
	// somewhere, which is what the height fog already assumes by having no
	// absorption term of its own to separate out.
	const float sigmaS = sigmaT;

	float3 scattered = 0;

	[branch]
	if (sigmaS > 0)
	{
		scattered = VolumetricFroxelScatteredLight(
			P, toEye, uv, (min16uint2)DTid.xy) * sigmaS;
	}

	output[DTid] = min(scattered, VOLUMETRIC_FROXEL_MAX_RADIANCE);
}
