/**
 * Ocean surface tessellation - a camera-centred geometry clipmap.
 *
 * The mesh is a set of nested square levels centred on the viewer. Level 0 is
 * a solid patch of `cells x cells` quads; every level above it is a square
 * annulus whose cell size doubles and whose hole is exactly the extent of the
 * level inside it. Vertex positions are generated here from `SV_VertexID`
 * alone - there is no vertex buffer, only two index buffers (patch and ring).
 *
 * Here a vertex depends on world position and camera *position* only - never
 * orientation, never a trace - so the near field is always covered at the
 * finest cell size and nothing can collapse.
 *
 * **Stability.** `xOceanMeshCenter` is snapped on the CPU, so vertices land on
 * a fixed world lattice and jump by whole cells instead of sliding. That is
 * what keeps a vertex sampling the same part of the wave field frame to frame.
 *
 * **Watertightness.** Every level shares one centre, so level `k`'s hole
 * coincides exactly with level `k-1`'s extent - no trim pieces, no overlap,
 * and no double blending of a transparent surface. The 2:1 density step at
 * each boundary is closed by the morph below rather than by skirt geometry.
 *
 * References:
 * https://hhoppe.com/geomclipmap.pdf
 * https://crest.readthedocs.io/en/4.16/user/technical-information.html
 */

#include "globals.hlsli"
#include "oceanSurfaceHF.hlsli"

Texture2D<float4> texture_displacementmap : register(t0);
Texture2D<float4> texture_perlin : register(t2);

PSIn main(
	uint vertexID : SV_VertexID,
	uint instanceID : SV_InstanceID,
	out uint RTIndex : SV_RenderTargetArrayIndex
)
{
	PSIn Out;

	// The instance index carries both the clipmap level and the camera. The
	// centre patch draw covers one level, the ring draw covers every level
	// above it, so the span follows from the base and needs no constant of its
	// own - the two draws only have to agree on this convention.
	const uint levelSpan = xOceanMeshLevelBase == 0
		? 1
		: xOceanMeshLevelCount - xOceanMeshLevelBase;
	const uint level = xOceanMeshLevelBase + (instanceID % levelSpan);
	Out.cameraIndex = (min16uint)(instanceID / levelSpan);

	const ShaderCamera camera = GetCameraIndexed(Out.cameraIndex);
	RTIndex = camera.output_index;

	// Vertex position within the level, as a signed offset in cells from the
	// level centre. Cells per side is even, so the outer boundary sits at an
	// even offset and the morph below cannot move it off the boundary.
	const int halfCells = (int)((xOceanMeshVertsPerSide - 1) / 2);
	const uint2 gridCoord = unflatten2D(vertexID, xOceanMeshVertsPerSide);
	const int2 cellOffset = (int2)gridCoord - halfCells;

	const float cellSize = xOceanMeshCellSize * (float)(1u << level);

	// Chebyshev distance from the centre, in cells: constant along the square
	// contour a level's edges follow, which is what the bands below key off.
	const int ringCoord = max(abs(cellOffset.x), abs(cellOffset.y));

	// Morph the outer band onto the lattice of the level outside this one, by
	// sliding odd vertices onto their even neighbour. At the boundary itself
	// the level therefore has the coarser level's vertex spacing and the two
	// edges share every vertex, so the 2:1 density step leaves no T-junction
	// to crack open. Collapsed pairs become zero-area triangles, which is the
	// normal cost of this technique.
	const float morphBand = (float)max(1u, xOceanMeshMorphCells);
	const float morphAlpha = saturate(
		((float)ringCoord - ((float)halfCells - morphBand)) / morphBand);
	const float2 parity = (float2)(uint2(abs(cellOffset)) & 1u);

	float2 surfaceXZ = xOceanMeshCenter
		+ ((float2)cellOffset - parity * morphAlpha) * cellSize;

	// A ring mesh is finite, so the outermost boundary is pushed out to the
	// far plane. This is the one thing the projected grid genuinely did well -
	// it reached the horizon by construction - and the water out there is flat
	// anyway, since the displacement has long since faded out.
	[branch]
	if (level + 1 == xOceanMeshLevelCount && ringCoord == halfCells)
	{
		const float outerExtent = (float)halfCells * cellSize;
		surfaceXZ = xOceanMeshCenter + (surfaceXZ - xOceanMeshCenter)
			* max(1.0, camera.z_far / outerExtent);
	}

	float3 surfacePos =
		float3(surfaceXZ.x, xOceanWaterHeight, surfaceXZ.y);

	// Sampled at the undisplaced position, and handed to the pixel shader as
	// is: the displacement map is a mapping from this parameter to a displaced
	// point, so the gradient the pixel shader fetches has to be read at the
	// same place the displacement was.
	const float2 uv = surfacePos.xz * xOceanPatchSizeRecip;

	// Both sides of every level boundary compute this from the same world
	// position, so they displace identically and the seam stays closed.
	const float3 displacement = texture_displacementmap.SampleLevel(
		sampler_linear_wrap, uv, 0).xzy;

	// Waves are flattened once the cells grow past the wave patch, since the
	// mesh cannot carry detail finer than a cell. `ocean_drawn_surface_height`
	// calls the same function, which is what keeps every test against the
	// water agreeing with the surface that was actually drawn.
	const float displacementFade = ocean_displacement_fade(
		distance(camera.position, surfacePos) * xOceanMeshFadeScale);
	surfacePos += lerp(displacement, 0, displacementFade);

	Out.pos = mul(camera.view_projection, float4(surfacePos, 1));
	Out.uv = uv;

	return Out;
}
