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

/**
 * Last frame's injection volume.
 *
 * Reusable precisely because this pass stores a **local** quantity: what the
 * medium at a world position scatters towards the eye does not depend on how
 * the volume happens to be sliced this frame, so a cell can be looked up
 * wherever that same world position sat last frame. The integrated volume could
 * never be filtered this way - it is a path integral along one particular ray,
 * and a ray from somewhere else is not a correction to it.
 */
Texture3D<float3> history : register(t0);

RWTexture3D<float3> output : register(u0);

/**
 * This block's freshly injected values, for bounding the history against.
 *
 * A thread has only its own sample, and one sample says nothing about how much
 * a cell's value should be trusted. Its neighbours in the same thread group are
 * free to reach, and the medium is smooth enough that they are a fair sample of
 * what this cell ought to hold.
 */
groupshared float3 blockScattered[
	VOLUMETRIC_FROXEL_BLOCKSIZE_X *
	VOLUMETRIC_FROXEL_BLOCKSIZE_Y *
	VOLUMETRIC_FROXEL_BLOCKSIZE_Z];

/**
 * Flattens a thread's position within its group into `blockScattered`.
 *
 * @param[in] threadID - Thread position within the group.
 *
 * @return Index into `blockScattered`.
 */
inline uint VolumetricFroxelBlockIndex(in uint3 threadID)
{
	return
		(threadID.z * VOLUMETRIC_FROXEL_BLOCKSIZE_Y + threadID.y)
			* VOLUMETRIC_FROXEL_BLOCKSIZE_X
		+ threadID.x;
}

[numthreads(
	VOLUMETRIC_FROXEL_BLOCKSIZE_X,
	VOLUMETRIC_FROXEL_BLOCKSIZE_Y,
	VOLUMETRIC_FROXEL_BLOCKSIZE_Z)]
void main(uint3 DTid : SV_DispatchThreadID, uint3 GTid : SV_GroupThreadID)
{
	// Offset within the cell rather than its centre, on all three axes. The
	// depth axis is what breaks the slice boundaries up; the screen axes stop
	// the shadow test landing on the same grid the volume is stored on, which
	// would otherwise print that grid along the edge of every shaft.
	const float3 jitter = VolumetricFroxelCellJitter(DTid, froxels.frame);

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

	blockScattered[VolumetricFroxelBlockIndex(GTid)] = scattered;
	GroupMemoryBarrierWithGroupSync();

	// Spread of this cell's neighbourhood, for bounding the history. Neighbours
	// outside the group are simply not counted: the block is 8x8x4 cells, so
	// most threads have a full set, and the ones on a face lose a little
	// confidence rather than correctness.
	float3 neighbourSum = 0;
	float3 neighbourSumSquared = 0;
	uint neighbourCount = 0;

	[unroll]
	for (int offsetZ = -1; offsetZ <= 1; ++offsetZ)
	{
		[unroll]
		for (int offsetY = -1; offsetY <= 1; ++offsetY)
		{
			[unroll]
			for (int offsetX = -1; offsetX <= 1; ++offsetX)
			{
				const int3 neighbour =
					int3(GTid) + int3(offsetX, offsetY, offsetZ);

				if (any(neighbour < 0) ||
					neighbour.x >= (int)VOLUMETRIC_FROXEL_BLOCKSIZE_X ||
					neighbour.y >= (int)VOLUMETRIC_FROXEL_BLOCKSIZE_Y ||
					neighbour.z >= (int)VOLUMETRIC_FROXEL_BLOCKSIZE_Z)
				{
					continue;
				}

				const float3 value =
					blockScattered[VolumetricFroxelBlockIndex((uint3)neighbour)];

				neighbourSum += value;
				neighbourSumSquared += value * value;
				++neighbourCount;
			}
		}
	}

	const float3 neighbourMean = neighbourSum / neighbourCount;
	const float3 neighbourDeviation = sqrt(max(
		0,
		(neighbourSumSquared / neighbourCount)
			- (neighbourMean * neighbourMean)));

	const float3 historyMinimum = neighbourMean
		- (VOLUMETRIC_FROXEL_HISTORY_CLIP_SIGMA * neighbourDeviation);
	const float3 historyMaximum = neighbourMean
		+ (VOLUMETRIC_FROXEL_HISTORY_CLIP_SIGMA * neighbourDeviation);

	// Temporal accumulation
	//==========================================================================
	//
	// One sample per cell is nowhere near enough to describe a cell that a
	// bright local light falls across, where the falloff goes as 1/d^2 and two
	// neighbouring sample points can differ by an order of magnitude. Averaging
	// over frames is what buys the sample count back, and the jitter above is
	// what gives each frame something new to say.
	float3 result = scattered;

	[branch]
	if (froxels.frame > 0)
	{
		// Reprojected from the cell's CENTRE, never from the jittered sample
		// point. The history holds one value per cell, and this cell's
		// predecessor is wherever this cell's own extent stood last frame -
		// which does not move just because the sample inside it did. Carrying
		// the jitter in here instead would resample the history at a different
		// sub-cell offset every frame, blending neighbouring cells together by
		// weights that change as fast as the jitter does: the average blurs
		// sideways a little more each frame, and the whole volume crawls.
		const float2 centreUV = (DTid.xy + 0.5) / float2(
			VOLUMETRIC_FROXEL_WIDTH, VOLUMETRIC_FROXEL_HEIGHT);
		const float centreDepth = VolumetricFroxelSliceToDepth(
			(float)DTid.z + 0.5, froxels.range);
		const float3 centreP =
			VolumetricFroxelPosition(centreUV, centreDepth);

		const float4 previousClip = mul(
			GetCamera().previous_view_projection, float4(centreP, 1));

		[branch]
		if (previousClip.w > 0)
		{
			const float2 previousUV =
				clipspace_to_uv(previousClip.xy / previousClip.w);

			// Distance from where the eye stood then, not where it stands now.
			// The volume is indexed radially, so a screen position alone does
			// not locate a point in it.
			const float previousDistance =
				length(centreP - froxels.previousEye);

			// Everything past the range collapses onto the last slice, where
			// one texel stands for the whole remaining column. Blending towards
			// that would drag the far field's answer into cells that have a
			// real one.
			[branch]
			if (is_saturated(previousUV) && previousDistance < froxels.range)
			{
				const float previousW = VolumetricFroxelDepthToW(
					previousDistance, froxels.range);

				const float3 previous = history.SampleLevel(
					sampler_linear_clamp,
					float3(previousUV, previousW),
					0);

				result = lerp(
					clamp(previous, historyMinimum, historyMaximum),
					scattered,
					VOLUMETRIC_FROXEL_HISTORY_BLEND);
			}
		}
	}

	output[DTid] = min(result, VOLUMETRIC_FROXEL_MAX_RADIANCE);
}
