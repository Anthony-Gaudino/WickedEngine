/**
 * Volumetric froxels - integration pass.
 *
 * Walks each view ray from the eye outwards, turning the per-metre scattering
 * the injection pass left in every cell into the light gathered between the eye
 * and that cell. One thread owns a whole column, because the running
 * transmittance is carried from one slice to the next and there is nothing to
 * share sideways.
 *
 * **Near to far, and it has to stay that way.** The medium's own extinction is
 * carried forward as a running product once there is one, and that is the only
 * form that is right for a column whose medium changes along it - a ray leaving
 * water into air, for instance. Marching the other way forces the transmittance
 * to be re-derived from the eye at every step, which can only be done in closed
 * form by assuming one medium the whole way.
 */
#include "globals.hlsli"
#include "fogHF.hlsli"
#include "volumetricFroxelHF.hlsli"

PUSHCONSTANT(froxels, VolumetricFroxelPush);

Texture3D<float3> input : register(t0);

RWTexture3D<float3> output : register(u0);

[numthreads(VOLUMETRIC_FROXEL_BLOCKSIZE_X, VOLUMETRIC_FROXEL_BLOCKSIZE_Y, 1)]
void main(uint3 DTid : SV_DispatchThreadID)
{
	const float2 uv = (DTid.xy + 0.5) / float2(
		VOLUMETRIC_FROXEL_WIDTH, VOLUMETRIC_FROXEL_HEIGHT);

	float3 accumulated = 0;
	float transmittance = 1;
	float previousDepth = 0;

	for (uint slice = 0; slice < VOLUMETRIC_FROXEL_SLICES; ++slice)
	{
		const float nearDepth = previousDepth;
		const float farDepth =
			VolumetricFroxelSliceToDepth((float)slice + 1, froxels.range);
		const float thickness = farDepth - nearDepth;
		previousDepth = farDepth;

		const float3 P = VolumetricFroxelPosition(
			uv, (nearDepth + farDepth) * 0.5);
		const float sigmaT = VolumetricFroxelAirExtinction(
			P, (nearDepth + farDepth) * 0.5);

		const float sliceTransmittance = exp(-thickness * sigmaT);

		// The exact integral of the transmittance across the slice, not the
		// transmittance at one end times its thickness. Slices are deliberately
		// far from uniform - the near ones are centimetres and the far ones
		// metres - and the closed form is what keeps a thick slice from being
		// wrong by however much the beam decayed while crossing it.
		//
		// Reduces to the thickness as the medium thins, which is why the
		// vanishing-extinction case needs no branch of its own.
		const float weight = sigmaT > 0.00001
			? (1 - sliceTransmittance) / sigmaT
			: thickness;

		accumulated += transmittance * input[uint3(DTid.xy, slice)] * weight;
		transmittance *= sliceTransmittance;

		output[uint3(DTid.xy, slice)] = accumulated;
	}
}
