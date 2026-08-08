#ifndef WI_SHADERINTEROP_OCEAN_H
#define WI_SHADERINTEROP_OCEAN_H
#include "ShaderInterop.h"

static const uint OCEAN_COMPUTE_TILESIZE = 8;

CBUFFER(OceanCB, CBSLOT_OTHER_OCEAN)
{
	// Foam colour in rgb, how opaque the foam is in w.
	float4 xOceanFoamColor;

	/**
	 * World XZ the whole clipmap hierarchy is centred on.
	 *
	 * Snapped on the CPU to the cell size of the last level that must be
	 * exactly world-aligned, and shared by **every** level. Sharing it is what
	 * makes each level's hole line up exactly with the extent of the level
	 * inside it, so the rings are watertight without any trim geometry.
	 */
	float2 xOceanMeshCenter;

	/** Cell size of the innermost level, in metres. */
	float xOceanMeshCellSize;

	/** Vertices along one side of a level, i.e. cells per side + 1. */
	uint xOceanMeshVertsPerSide;

	/**
	 * Width of the morph band at the outer edge of a level, in cells.
	 *
	 * Across this band a level's vertices slide onto the lattice of the level
	 * outside it, so the 2:1 density step at the boundary does not leave
	 * T-junctions.
	 */
	uint xOceanMeshMorphCells;

	/**
	 * First clipmap level covered by the current draw.
	 *
	 * 0 for the solid centre patch, 1 for the instanced ring draw. Levels per
	 * instance follow from this and `xOceanMeshLevelCount`, so no separate
	 * span is needed.
	 */
	uint xOceanMeshLevelBase;

	/** Total number of clipmap levels, counting the centre patch. */
	uint xOceanMeshLevelCount;

	float xOceanTexelLength;

	float xOceanPatchSizeRecip;
	float xOceanMapHalfTexel;
	float xOceanWaterHeight;
	uint xOceanActualDim;

	uint xOceanInWidth;
	uint xOceanOutWidth;
	uint xOceanOutHeight;
	// Element offset to the second packed FFT field (Dy). The first field packs
	// height (real output) + Dx (imaginary output) via the two-for-one real FFT.
	uint xOceanSecondFieldOffset;

	float xOceanTimeScale;
	float xOceanChoppyScale;
	float xOceanGridLen;
	float xOceanWaveAmplitude;

	/**
	 * Scales the distance this mesh reports when asking how far the waves have
	 * flattened.
	 *
	 * The published fade band is derived for the main surface mesh. The coarse
	 * mesh used for shadow maps, cubemaps and the occlusion test has the same
	 * world extent from far fewer cells, so its cells reach any given size
	 * much sooner and it has to start flattening correspondingly sooner - a
	 * mesh with one cell per wave patch displaced at full amplitude is
	 * sampling noise.
	 *
	 * **1 for the main mesh**, so the surface and everything testing against
	 * it read the band identically and cannot drift.
	 */
	float xOceanMeshFadeScale;

	float3 xOcean_padding;
};

#endif // WI_SHADERINTEROP_OCEAN_H
