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
 * carried forward as a running product, which is the only form that is right
 * for a column the medium varies along - height fog thins with altitude, so a
 * ray climbing away from the ground crosses less of it the further it goes.
 * Marching the other way forces the transmittance to be re-derived from the eye
 * at every step, which can only be done in closed form by assuming the medium
 * is uniform the whole way.
 *
 * The world beyond the volume is not left unlit: a second, much shorter march
 * carries the column from the last slice out to the far plane, and is published
 * separately for the fragments out there to fade in along. See the tail at the
 * foot of this file.
 */

// Matching the injection pass, whose samples this one continues past the range.
// Ahead of the includes for the same reason it is there.
#define DISABLE_SOFT_SHADOWMAP

#include "globals.hlsli"
#include "lightingHF.hlsli"
#include "fogHF.hlsli"
#include "volumetricFroxelLightingHF.hlsli"

PUSHCONSTANT(froxels, VolumetricFroxelPush);

Texture3D<float3> input : register(t0);

RWTexture3D<float3> output : register(u0);

/**
 * The column beyond the volume, one value per screen column.
 *
 * `rgb` is the light gathered from the last slice out to the far plane; `a` is
 * the extinction the medium averages over that stretch, which is what gives the
 * far field a shape to fade in along rather than a step to jump at.
 */
RWTexture2D<float4> tail : register(u1);

[numthreads(VOLUMETRIC_FROXEL_BLOCKSIZE_X, VOLUMETRIC_FROXEL_BLOCKSIZE_Y, 1)]
void main(uint3 DTid : SV_DispatchThreadID)
{
	const float2 uv = (DTid.xy + 0.5) / float2(
		VOLUMETRIC_FROXEL_WIDTH, VOLUMETRIC_FROXEL_HEIGHT);

	float3 accumulated = 0;

	// Per channel, because water is not grey: it takes red out of a beam
	// roughly an order of magnitude faster than blue, which is what leaves a
	// submerged lamp warm within arm's reach and blue a few metres out.
	float3 transmittance = 1;

	// The state at the last texel's own depth, which is where the tail has to
	// start if the two are to meet without a step.
	float3 lastTexelTransmittance = 1;
	float lastTexelDepth = 0;

	for (uint slice = 0; slice < VOLUMETRIC_FROXEL_SLICES; ++slice)
	{
		// A texel does not stand for the far edge of its slice, it stands for
		// the MIDDLE: `VolumetricFroxelDepthToW` sends the depth at slice+0.5
		// to texel `slice`, and the injection samples there too. So the sum has
		// to be split, banked at the middle, and carried on - accumulating the
		// whole slice before storing would file each value half a slice further
		// out than the depth it is read back at. The error grows with slice
		// thickness, so the reconstruction overshoots at every texel and sags
		// between them: one ripple per slice, seen as arcs centred on the eye.
		const float nearDepth =
			VolumetricFroxelSliceToDepth((float)slice, froxels.range);
		const float centreDepth =
			VolumetricFroxelSliceToDepth((float)slice + 0.5, froxels.range);
		const float farDepth =
			VolumetricFroxelSliceToDepth((float)slice + 1, froxels.range);

		// Not equal halves. The slices are squared in depth, so the near half of
		// one is always the thinner.
		const float frontThickness = centreDepth - nearDepth;
		const float backThickness = farDepth - centreDepth;

		// Taken at the middle, which is both where the texel sits and where the
		// injection's jitter averages out to.
		//
		// **Recomputed rather than stored beside the injected light.** The
		// medium is a smooth analytic field with no noise in it to filter, so
		// it gains nothing from the temporal average the injection volume
		// carries - and it would lose by it, smearing the waterline across the
		// frames a camera takes to cross the surface.
		const float3 sigmaT =
			VolumetricFroxelMediumAt(uv, centreDepth).sigmaT;

		const float3 injected = input[uint3(DTid.xy, slice)];

		// The exact integral of the transmittance across the segment, not the
		// transmittance at one end times its length. Slices are deliberately far
		// from uniform - the near ones are centimetres and the far ones metres -
		// and the closed form is what keeps a thick one from being wrong by
		// however much the beam decayed while crossing it.
		//
		// Reduces to the length as the medium thins, which is why the
		// vanishing-extinction case needs no branch of its own: the floor is
		// small enough that the ratio is the length to within a rounding error
		// over any slice, and a medium that thin scatters nothing to weigh.
		const float3 safeSigmaT = max(sigmaT, 0.00001);

		const float3 frontTransmittance = exp(-frontThickness * safeSigmaT);
		const float3 frontWeight = (1 - frontTransmittance) / safeSigmaT;

		accumulated += transmittance * injected * frontWeight;
		transmittance *= frontTransmittance;

		// Clamped on the way out as well as on the way in. This is a running
		// sum over 128 slices, so it can reach the format's ceiling even where
		// no single cell came close to it.
		output[uint3(DTid.xy, slice)] =
			min(accumulated, VOLUMETRIC_FROXEL_MAX_RADIANCE);

		if (slice == VOLUMETRIC_FROXEL_SLICES - 1)
		{
			lastTexelTransmittance = transmittance;
			lastTexelDepth = centreDepth;
		}

		const float3 backTransmittance = exp(-backThickness * safeSigmaT);
		const float3 backWeight = (1 - backTransmittance) / safeSigmaT;

		accumulated += transmittance * injected * backWeight;
		transmittance *= backTransmittance;
	}

	// The tail
	//==========================================================================
	//
	// The volume stops at an authored range; the view does not. Marched here
	// rather than injected, because it has no cells: one segment per column
	// instead of one per slice, so it is a few hundred thousand samples against
	// the volume's seven million.
	//
	// **Kept out of the volume.** Writing the whole remaining column into the
	// last texel is the obvious thing and it is wrong: that texel spans about
	// six metres of depth, so the entire stretch from the range to the far plane
	// arrives over six metres and lands as a hard arc across the scene at
	// exactly the range. Handed to the fragments as a total and an extinction
	// instead, it can be faded in over the thousands of metres it actually
	// occupies.
	//
	// What is lost is any structure within the tail - a shaft crossing it is one
	// smooth ramp, not a shape. That is the honest end of the trade: the
	// alternative is a range large enough to reach a horizon, which spends the
	// whole volume's resolution out where nothing needs it.
	const float tailStart = lastTexelDepth;
	const float tailEnd = GetCamera().z_far;

	float3 tailInscatter = 0;
	float tailExtinction = 0;

	[branch]
	if (tailEnd > tailStart)
	{
		const float tailSpan = tailEnd - tailStart;
		const float tailStep = tailSpan / (float)VOLUMETRIC_FROXEL_TAIL_STEPS;

		float3 accumulatedTail = 0;
		float3 tailTransmittance = 1;

		for (uint step = 0; step < VOLUMETRIC_FROXEL_TAIL_STEPS; ++step)
		{
			const float stepDepth = tailStart + ((float)step + 0.5) * tailStep;

			const float3 P = VolumetricFroxelPosition(uv, stepDepth);
			const VolumetricFroxelMedium medium =
				VolumetricFroxelMediumAt(uv, stepDepth);

			float3 scattered = 0;

			[branch]
			if (medium.Scatters())
			{
				// One position for both, the step's own middle: the tail is
				// marched rather than jittered, so there is no separate point
				// to read the occlusion at.
				scattered = VolumetricFroxelScatteredLight(
					medium,
					P,
					P,
					normalize(GetCamera().position - P),
					uv,
					(min16uint2)DTid.xy);
			}

			const float3 safeSigmaT = max(medium.sigmaT, 0.00001);

			const float3 stepTransmittance = exp(-tailStep * safeSigmaT);
			const float3 weight = (1 - stepTransmittance) / safeSigmaT;

			accumulatedTail += tailTransmittance * scattered * weight;
			tailTransmittance *= stepTransmittance;
		}

		// Attenuated by everything the volume already crossed, since this
		// stretch is seen through all of it.
		tailInscatter = lastTexelTransmittance * accumulatedTail;

		// One extinction standing for the whole stretch, recovered from what it
		// did to the beam. It is the shape the fade follows, so that light
		// arriving from the far field comes in the way a medium delivers it -
		// most of it early, tailing off - rather than in a straight line.
		//
		// Scalar where the medium is not, taken on the channel that carries
		// furthest: what is still being delivered at the far end of the tail is
		// whatever has not been extinguished yet, so that channel is the one
		// with a shape left to describe. Nothing is lost in air, where all
		// three are the same number, and nothing is lost under water either -
		// several hundred metres of it stand in front of the tail, and no
		// channel arrives through that.
		const float3 perChannel =
			-log(max(tailTransmittance, 0.0001)) / tailSpan;

		tailExtinction = min3(perChannel.r, perChannel.g, perChannel.b);
	}

	tail[DTid.xy] = float4(
		min(tailInscatter, VOLUMETRIC_FROXEL_MAX_RADIANCE), tailExtinction);
}
