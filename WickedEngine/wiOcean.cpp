#include "wiOcean.h"
#include "wiRenderer.h"
#include "wiResourceManager.h"
#include "shaders/ShaderInterop_Ocean.h"
#include "wiScene.h"
#include "wiBacklog.h"
#include "wiEventHandler.h"
#include "wiTimer.h"
#include "wiVector.h"

#include "perlin.h"

#include <algorithm>

using namespace wi::graphics;
using namespace wi::scene;

namespace wi
{
	namespace ocean_internal
	{
		Shader		updateSpectrumCS;
		Shader		updateDisplacementMapCS;
		Shader		updateGradientFoldingCS;
		Shader		oceanSurfVS;
		Shader		wireframePS;
		Shader		oceanSurfPS;
		Shader		oceanSurfPS_envmap;
		Shader		oceanSurfPS_shadowmap;

		RasterizerState		rasterizerState;
		RasterizerState		wireRS;
		DepthStencilState	depthStencilState, depthStencilState_occlusionTest, depthStencilState_shadowmap;
		BlendState			blendState, blendState_occlusionTest, blendState_shadowmap;

		PipelineState PSO, PSO_envmap, PSO_shadowmap, PSO_wire, PSO_occlusionTest;
		Texture perlinTex;

		/*
		########################################################################
		Clipmap mesh
		########################################################################
		*/

		/**
		 * Cells along one side of a clipmap level, for the main surface draw.
		 *
		 * Must be a multiple of 4: the outer boundary has to sit at an even
		 * cell offset so the vertex shader's morph cannot slide it off the
		 * boundary, and the centre hole boundary must land on a cell edge.
		 *
		 * `(cells + 1)^2` must stay under 65536 so the index buffers can be
		 * 16-bit; 128 gives 16641 vertices.
		 */
		constexpr uint32_t MESH_CELLS_PER_SIDE = 128;

		/**
		 * Cells along one side of a clipmap level, for the cubemap, shadow map
		 * and occlusion draws.
		 *
		 * These only need the ocean's silhouette, not its detail, so they use
		 * a much coarser mesh. The base cell size is scaled up to compensate,
		 * which keeps every level's world extent identical between the two
		 * resolutions - including level 0's, which the snap quantum below is
		 * sized against.
		 */
		constexpr uint32_t MESH_CELLS_PER_SIDE_COARSE = 16;

		/** Cell size of the innermost level at `MESH_CELLS_PER_SIDE`, in metres. */
		constexpr float MESH_BASE_CELL_SIZE = 0.25F;

		/** Number of clipmap levels, counting the solid centre patch. */
		constexpr uint32_t MESH_LEVEL_COUNT = 10;

		/** Width of the morph band at the outer edge of a level, in cells. */
		constexpr uint32_t MESH_MORPH_CELLS = MESH_CELLS_PER_SIDE / 8;

		/**
		 * Cells per wave patch below which the displacement starts flattening.
		 *
		 * A patch spanned by this many cells still carries its shape; below it
		 * the mesh is throwing wave detail away whether or not the
		 * displacement is applied, so applying it only buys noise.
		 */
		constexpr float MESH_FADE_START_CELLS_PER_PATCH = 32.0F;

		/**
		 * Cells per wave patch at which the displacement is fully flattened.
		 *
		 * The hard limit here is **2** - one patch is the longest wavelength
		 * the FFT produces, so below two cells per patch the mesh cannot
		 * represent even the fundamental and the displacement degenerates into
		 * aliasing. Stopping well short of it is deliberate rather than
		 * timid: there is only one FFT patch, so pushing displacement further
		 * out mostly buys a visible repetition of the same 50 m tile. The
		 * cascaded FFT is what would earn the rest of the range.
		 */
		constexpr float MESH_FADE_END_CELLS_PER_PATCH = 8.0F;

		/**
		 * Grid the shared mesh centre is snapped to, in metres.
		 *
		 * Snapping is what stops the waves swimming: vertices land on a fixed
		 * world lattice and jump by whole cells as the viewer moves, instead
		 * of sliding across the displacement map.
		 *
		 * One quantum is shared by **every** level, which is what makes each
		 * level's hole coincide exactly with the extent of the level inside it
		 * - the alternative, snapping each level to its own cell size, leaves
		 * the two disagreeing by up to a cell and needs trim geometry to
		 * close. Levels whose cell size exceeds this quantum can slide by at
		 * most one quantum, which is under half their cell and is out where
		 * the displacement has faded to nothing anyway.
		 *
		 * The bound that matters: half of this is how far the viewer can drift
		 * from the mesh centre, and it must stay well inside level 0's half
		 * extent (`MESH_CELLS_PER_SIDE / 2 * MESH_BASE_CELL_SIZE`, 16 m) or
		 * the viewer ends up near the edge of the finest ring.
		 */
		constexpr float MESH_SNAP_QUANTUM =
			MESH_BASE_CELL_SIZE * float(1u << 5); // 8 m

		GPUBuffer indexBuffer_patch;
		GPUBuffer indexBuffer_ring;
		GPUBuffer indexBuffer_patch_coarse;
		GPUBuffer indexBuffer_ring_coarse;

		/**
		 * Builds the index buffer for one clipmap level shape.
		 *
		 * The vertex positions themselves are generated in the vertex shader
		 * from `SV_VertexID`, so this only describes connectivity over a
		 * `(cellsPerSide + 1)^2` vertex grid and is independent of the ocean
		 * parameters - hence built once at startup and shared by every ocean.
		 *
		 * @param[in] cellsPerSide - Cells along one side of the level.
		 * @param[in] ring - When true, omit the centre quarter of the grid,
		 *                   which the level inside this one covers exactly.
		 * @param[out] indexBuffer - Buffer to create.
		 * @param[in] debugName - Name to tag the buffer with.
		 */
		void CreateClipmapIndexBuffer(
			const uint32_t cellsPerSide,
			const bool ring,
			GPUBuffer& indexBuffer,
			const char* debugName
		)
		{
			const uint32_t vertsPerSide = cellsPerSide + 1;
			const uint32_t holeMin = cellsPerSide / 4;
			const uint32_t holeMax = cellsPerSide - holeMin;

			wi::vector<uint16_t> indexData;
			indexData.reserve(size_t(cellsPerSide) * cellsPerSide * 6);

			for (uint32_t x = 0; x < cellsPerSide; ++x)
			{
				for (uint32_t y = 0; y < cellsPerSide; ++y)
				{
					const bool inHole =
						x >= holeMin && (x + 1) <= holeMax &&
						y >= holeMin && (y + 1) <= holeMax;
					if (ring && inHole)
					{
						continue;
					}

					const uint16_t lowerLeft =
						uint16_t(x + y * vertsPerSide);
					const uint16_t lowerRight =
						uint16_t((x + 1) + y * vertsPerSide);
					const uint16_t topLeft =
						uint16_t(x + (y + 1) * vertsPerSide);
					const uint16_t topRight =
						uint16_t((x + 1) + (y + 1) * vertsPerSide);

					indexData.push_back(topLeft);
					indexData.push_back(lowerLeft);
					indexData.push_back(lowerRight);

					indexData.push_back(topLeft);
					indexData.push_back(lowerRight);
					indexData.push_back(topRight);
				}
			}

			GPUBufferDesc desc;
			desc.bind_flags = BindFlag::INDEX_BUFFER;
			desc.size = indexData.size() * sizeof(uint16_t);

			GraphicsDevice* device = wi::graphics::GetDevice();
			device->CreateBuffer(&desc, indexData.data(), &indexBuffer);
			device->SetName(&indexBuffer, debugName);
		}

		/**
		 * Number of indices in a clipmap level's index buffer.
		 *
		 * @param[in] indexBuffer - Buffer built by `CreateClipmapIndexBuffer`.
		 *
		 * @return Index count to pass to the draw.
		 */
		uint32_t GetIndexCount(const GPUBuffer& indexBuffer)
		{
			return uint32_t(indexBuffer.GetDesc().size / sizeof(uint16_t));
		}

		/**
		 * Issues the two draws that make up one clipmap.
		 *
		 * The centre patch and the rings have different connectivity, so they
		 * cannot share an index buffer and this is two draws rather than one.
		 * Which level a vertex belongs to rides in `SV_InstanceID` alongside
		 * the camera index; `xOceanMeshLevelBase` tells the vertex shader how
		 * to unpack it, and the convention set here has to match the one it
		 * reads.
		 *
		 * @param[in,out] constants - Ocean constants. `xOceanMeshLevelBase`
		 *                            is overwritten and the buffer rebound per
		 *                            draw.
		 * @param[in] coarse - Use the low resolution mesh, which carries the
		 *                     ocean's silhouette but not its detail.
		 * @param[in] cameraCount - Cameras to replicate the mesh across; 6 for
		 *                          a cubemap, 1 otherwise.
		 * @param[in] cmd - Command list to record into.
		 */
		void DrawClipmap(
			OceanCB& constants,
			const bool coarse,
			const uint32_t cameraCount,
			const CommandList cmd
		)
		{
			GraphicsDevice* device = wi::graphics::GetDevice();

			const GPUBuffer& patch =
				coarse ? indexBuffer_patch_coarse : indexBuffer_patch;
			const GPUBuffer& ring =
				coarse ? indexBuffer_ring_coarse : indexBuffer_ring;

			constants.xOceanMeshLevelBase = 0;
			device->BindDynamicConstantBuffer(
				constants, CB_GETBINDSLOT(OceanCB), cmd);
			device->BindIndexBuffer(&patch, IndexBufferFormat::UINT16, 0, cmd);
			device->DrawIndexedInstanced(
				GetIndexCount(patch), cameraCount, 0, 0, 0, cmd);

			if constexpr (MESH_LEVEL_COUNT > 1)
			{
				constants.xOceanMeshLevelBase = 1;
				device->BindDynamicConstantBuffer(
					constants, CB_GETBINDSLOT(OceanCB), cmd);
				device->BindIndexBuffer(
					&ring, IndexBufferFormat::UINT16, 0, cmd);
				device->DrawIndexedInstanced(
					GetIndexCount(ring),
					(MESH_LEVEL_COUNT - 1) * cameraCount,
					0, 0, 0, cmd);
			}
		}

		void LoadShaders()
		{
			wi::renderer::LoadShader(ShaderStage::CS, updateSpectrumCS, "oceanSimulatorCS.cso");
			wi::renderer::LoadShader(ShaderStage::CS, updateDisplacementMapCS, "oceanUpdateDisplacementMapCS.cso");
			wi::renderer::LoadShader(ShaderStage::CS, updateGradientFoldingCS, "oceanUpdateGradientFoldingCS.cso");

			wi::renderer::LoadShader(ShaderStage::VS, oceanSurfVS, "oceanSurfaceVS.cso");

			wi::renderer::LoadShader(ShaderStage::PS, oceanSurfPS, "oceanSurfacePS.cso");
			wi::renderer::LoadShader(ShaderStage::PS, oceanSurfPS_envmap, "oceanSurfacePS_envmap.cso");
			wi::renderer::LoadShader(ShaderStage::PS, oceanSurfPS_shadowmap, "oceanSurfacePS_shadowmap.cso");
			wi::renderer::LoadShader(ShaderStage::PS, wireframePS, "oceanSurfaceSimplePS.cso");


			GraphicsDevice* device = wi::graphics::GetDevice();

			{
				PipelineStateDesc desc;
				desc.vs = &oceanSurfVS;
				desc.ps = &oceanSurfPS;
				desc.bs = &blendState;
				desc.rs = &rasterizerState;
				desc.dss = &depthStencilState;
				device->CreatePipelineState(&desc, &PSO);

				desc.ps = &oceanSurfPS_envmap;
				device->CreatePipelineState(&desc, &PSO_envmap);

				desc.ps = &oceanSurfPS_shadowmap;
				desc.bs = &blendState_shadowmap;
				desc.dss = &depthStencilState_shadowmap;
				device->CreatePipelineState(&desc, &PSO_shadowmap);

				desc.bs = &blendState;
				desc.ps = &wireframePS;
				desc.rs = &wireRS;
				desc.dss = &depthStencilState;
				device->CreatePipelineState(&desc, &PSO_wire);

				desc.ps = {};
				desc.rs = &rasterizerState;
				desc.bs = &blendState_occlusionTest;
				desc.dss = &depthStencilState_occlusionTest;
				device->CreatePipelineState(&desc, &PSO_occlusionTest);
			}
		}
	}
	using namespace ocean_internal;


#define HALF_SQRT_2	0.7071068f
#define GRAV_ACCEL	981.0f	// The acceleration of gravity, cm/s^2

	// Generating gaussian random number with mean 0 and standard deviation 1.
	float Gauss()
	{
		float u1 = rand() / (float)RAND_MAX;
		float u2 = rand() / (float)RAND_MAX;
		if (u1 < 1e-6f)
			u1 = 1e-6f;
		return std::sqrt(-2 * logf(u1)) * cosf(2 * XM_PI * u2);
	}

	// Phillips Spectrum
	// K: normalized wave vector, W: wind direction, v: wind velocity, a: amplitude constant
	float Phillips(XMFLOAT2 K, XMFLOAT2 W, float v, float a, float dir_depend)
	{
		// largest possible wave from constant wind of velocity v
		float l = v * v / GRAV_ACCEL;
		// damp out waves with very small length w << l
		float w = l / 1000;

		float Ksqr = K.x * K.x + K.y * K.y;
		float Kcos = K.x * W.x + K.y * W.y;
		float phillips = a * expf(-1 / (l * l * Ksqr)) / (Ksqr * Ksqr * Ksqr) * (Kcos * Kcos);

		// filter out waves moving opposite to wind
		if (Kcos < 0)
			phillips *= dir_depend;

		// damp out waves with very small length w << l
		return phillips * expf(-Ksqr * w * w);
	}



	void Ocean::Create(const OceanParameters& params)
	{
		this->params = params;
		for (int i = 0; i < arraysize(occlusionQueries); ++i)
		{
			occlusionQueries[i] = -1;
		}

		GraphicsDevice* device = wi::graphics::GetDevice();

		// Height map H(0)
		int height_map_size = (params.dmap_dim + 4) * (params.dmap_dim + 1);
		wi::vector<XMFLOAT2> h0_data(height_map_size);
		wi::vector<float> omega_data(height_map_size);
		initHeightMap(h0_data.data(), omega_data.data());

		int hmap_dim = params.dmap_dim;
		int input_full_size = (hmap_dim + 4) * (hmap_dim + 1);
		// This value should be (hmap_dim / 2 + 1) * hmap_dim, but we use full sized buffer here for simplicity.
		int input_half_size = hmap_dim * hmap_dim;
		int output_size = hmap_dim * hmap_dim;

		GPUBufferDesc buf_desc;
		buf_desc.usage = Usage::DEFAULT;
		buf_desc.bind_flags = BindFlag::UNORDERED_ACCESS | BindFlag::SHADER_RESOURCE;
		buf_desc.misc_flags = ResourceMiscFlag::BUFFER_STRUCTURED;

		// RW buffer allocations
		// H0
		buf_desc.stride = sizeof(float2);
		buf_desc.size = buf_desc.stride * input_full_size;
		device->CreateBuffer(&buf_desc, h0_data.data(), &buffer_Float2_H0);

		// Two packed FFT fields instead of three: the two-for-one real FFT packs
		// H(t) and Dx(t) into the real/imaginary parts of one complex field, and
		// Dy(t) into a second. See oceanSimulatorCS / oceanUpdateDisplacementMapCS.
		buf_desc.stride = sizeof(float2);
		buf_desc.size = buf_desc.stride * 2 * input_half_size;
		device->CreateBufferZeroed(&buf_desc, &buffer_Float2_Ht);

		// omega
		buf_desc.stride = sizeof(float);
		buf_desc.size = buf_desc.stride * input_full_size;
		device->CreateBuffer(&buf_desc, omega_data.data(), &buffer_Float_Omega);

		// Two packed output fields (see above): field 0 holds Dz (real) + Dx
		// (imaginary), field 1 holds Dy. The C2C FFT output is complex, so the
		// stride stays sizeof(float2).
		buf_desc.stride = sizeof(float2);
		buf_desc.size = buf_desc.stride * 2 * output_size;
		device->CreateBufferZeroed(&buf_desc, &buffer_Float_Dxyz);

		TextureDesc tex_desc;
		tex_desc.width = hmap_dim;
		tex_desc.height = hmap_dim;
		tex_desc.array_size = 1;
		tex_desc.sample_count = 1;
		tex_desc.usage = Usage::DEFAULT;
		tex_desc.bind_flags = BindFlag::SHADER_RESOURCE | BindFlag::UNORDERED_ACCESS;

		tex_desc.format = Format::R16G16B16A16_FLOAT;
		tex_desc.mip_levels = GetMipCount(tex_desc.width, tex_desc.height);
		tex_desc.layout = ResourceState::SHADER_RESOURCE_COMPUTE;
		device->CreateTexture(&tex_desc, nullptr, &gradientMap);
		device->SetName(&gradientMap, "gradientMap");

		for (uint32_t i = 0; i < gradientMap.GetDesc().mip_levels; ++i)
		{
			int subresource_index;
			subresource_index = device->CreateSubresource(&gradientMap, SubresourceType::SRV, 0, 1, i, 1);
			assert(subresource_index == i);
			subresource_index = device->CreateSubresource(&gradientMap, SubresourceType::UAV, 0, 1, i, 1);
			assert(subresource_index == i);
		}

		// Half-float displacement: the values are zero-mean wave offsets, so 16-bit
		// float's ~0.05% relative precision is ample for both rendering and the CPU
		// buoyancy readback, while halving the texture's bandwidth and memory.
		tex_desc.format = Format::R16G16B16A16_FLOAT;
		tex_desc.mip_levels = 1;
		tex_desc.layout = ResourceState::COPY_SRC | ResourceState::SHADER_RESOURCE_COMPUTE;
		device->CreateTextureZeroed(&tex_desc, &displacementMap); // zero init the heightmap to be valid before first simulation
		device->SetName(&displacementMap, "displacementMap");

		for (uint32_t i = 0; i < arraysize(displacementMap_readback); ++i)
		{
			tex_desc.usage = Usage::READBACK;
			tex_desc.bind_flags = BindFlag::NONE;
			tex_desc.layout = ResourceState::COPY_DST;
			device->CreateTexture(&tex_desc, nullptr, &displacementMap_readback[i]);
			device->SetName(&displacementMap_readback[i], "displacementMap_readback[i]");
		}
	}

	// Keep the displacement map readback alive for this many frames after the last
	// GetDisplacedPosition() query. Must comfortably exceed the GPU->CPU readback
	// latency (GetBufferCount() frames) so that continuous but irregular queries
	// don't repeatedly drop and re-warm the readback chain.
	static constexpr uint32_t DISPLACEMENT_READBACK_KEEPALIVE_FRAMES = 8;

	XMFLOAT3 Ocean::GetDisplacedPosition(const XMFLOAT3& worldPosition) const
	{
		// Signal that ocean height is being queried on the CPU, so the readback copy
		// in CopyDisplacementMapReadback() stays enabled.
		displacement_readback_request = DISPLACEMENT_READBACK_KEEPALIVE_FRAMES;

		XMFLOAT3 ocean_pos = XMFLOAT3(worldPosition.x, params.waterHeight, worldPosition.z);
		if (displacement_readback_valid[displacement_readback_index])
		{
			const Texture& tex = displacementMap_readback[displacement_readback_index];
			const uint32_t row_pitch = tex.mapped_subresources[0].row_pitch;
			const uint8_t* bytedata = (const uint8_t*)tex.mapped_data;
			if (bytedata != nullptr)
			{
				const float patch_size_rcp = 1.0f / params.patch_length;
				const XMFLOAT2 uv = XMFLOAT2(frac(ocean_pos.x * patch_size_rcp), frac(ocean_pos.z * patch_size_rcp));
				const XMFLOAT2 fpixel = XMFLOAT2(uv.x * (tex.desc.width - 1), uv.y * (tex.desc.height - 1));
				const XMFLOAT2 pixel_frac = XMFLOAT2(frac(fpixel.x), frac(fpixel.y));
				const XMUINT2 pixel_tl = XMUINT2((uint32_t)std::floor(fpixel.x), (uint32_t)std::floor(fpixel.y));
				const XMUINT2 pixel_tr = XMUINT2((uint32_t)std::ceil(fpixel.x), (uint32_t)std::floor(fpixel.y));
				const XMUINT2 pixel_bl = XMUINT2((uint32_t)std::floor(fpixel.x), (uint32_t)std::ceil(fpixel.y));
				const XMUINT2 pixel_br = XMUINT2((uint32_t)std::ceil(fpixel.x), (uint32_t)std::ceil(fpixel.y));
				const XMFLOAT4 displacement_tl = wi::math::unpack_half4(((const XMUINT2*)(bytedata + pixel_tl.y * row_pitch))[pixel_tl.x]);
				const XMFLOAT4 displacement_tr = wi::math::unpack_half4(((const XMUINT2*)(bytedata + pixel_tr.y * row_pitch))[pixel_tr.x]);
				const XMFLOAT4 displacement_bl = wi::math::unpack_half4(((const XMUINT2*)(bytedata + pixel_bl.y * row_pitch))[pixel_bl.x]);
				const XMFLOAT4 displacement_br = wi::math::unpack_half4(((const XMUINT2*)(bytedata + pixel_br.y * row_pitch))[pixel_br.x]);
				const XMFLOAT4 xxxx = XMFLOAT4(displacement_bl.x, displacement_br.x, displacement_tr.x, displacement_tl.x);
				const XMFLOAT4 yyyy = XMFLOAT4(displacement_bl.y, displacement_br.y, displacement_tr.y, displacement_tl.y);
				const XMFLOAT4 zzzz = XMFLOAT4(displacement_bl.z, displacement_br.z, displacement_tr.z, displacement_tl.z);
				const XMFLOAT4 wwww = XMFLOAT4(displacement_bl.w, displacement_br.w, displacement_tr.w, displacement_tl.w);
				const XMFLOAT4 displacement = XMFLOAT4(
					bilinear(xxxx, pixel_frac),
					bilinear(yyyy, pixel_frac),
					bilinear(zzzz, pixel_frac),
					bilinear(wwww, pixel_frac)
				);
				// xzy swizzle is on purpose, that's how the data is generated on GPU:
				ocean_pos.x += displacement.x;
				ocean_pos.y += displacement.z;
				ocean_pos.z += displacement.y;
			}
		}
		return ocean_pos;
	}

	// Initialize the vector field.
	// wlen_x: width of wave tile, in meters
	// wlen_y: length of wave tile, in meters
	void Ocean::initHeightMap(XMFLOAT2* out_h0, float* out_omega)
	{
		int i, j;
		XMFLOAT2 K;

		XMFLOAT2 wind_dir;
		XMStoreFloat2(&wind_dir, XMVector2Normalize(XMLoadFloat2(&params.wind_dir)));
		float a = params.wave_amplitude * 1e-7f;	// It is too small. We must scale it for editing.
		float v = params.wind_speed;
		float dir_depend = params.wind_dependency;

		int height_map_dim = params.dmap_dim;
		float patch_length = params.patch_length;


		for (i = 0; i <= height_map_dim; i++)
		{
			// K is wave-vector, range [-|DX/W, |DX/W], [-|DY/H, |DY/H]
			K.y = (-height_map_dim / 2.0f + i) * (2 * XM_PI / patch_length);

			for (j = 0; j <= height_map_dim; j++)
			{
				K.x = (-height_map_dim / 2.0f + j) * (2 * XM_PI / patch_length);

				float phil = (K.x == 0 && K.y == 0) ? 0 : std::sqrt(Phillips(K, wind_dir, v, a, dir_depend));

				out_h0[i * (height_map_dim + 4) + j].x = float(phil * Gauss() * HALF_SQRT_2);
				out_h0[i * (height_map_dim + 4) + j].y = float(phil * Gauss() * HALF_SQRT_2);

				// The angular frequency is following the dispersion relation:
				//            out_omega^2 = g*k
				// The equation of Gerstner wave:
				//            x = x0 - K/k * A * sin(dot(K, x0) - sqrt(g * k) * t), x is a 2D vector.
				//            z = A * cos(dot(K, x0) - sqrt(g * k) * t)
				// Gerstner wave shows that a point on a simple sinusoid wave is doing a uniform circular
				// motion with the center (x0, y0, z0), radius A, and the circular plane is parallel to
				// vector K.
				out_omega[i * (height_map_dim + 4) + j] = std::sqrt(GRAV_ACCEL * sqrtf(K.x * K.x + K.y * K.y));
			}
		}
	}

	/**
	 * Fills the ocean constants for a draw or simulation dispatch.
	 *
	 * @param[in] params - Ocean parameters to publish.
	 * @param[in] viewerPosition - World position the clipmap centres on. This
	 *                             is the **main** camera even for the shadow
	 *                             map, which renders the same stretch of sea
	 *                             from the light.
	 * @param[in] cellsPerSide - Cells along one side of a clipmap level. The
	 *                           base cell size is scaled so that every level's
	 *                           world extent is the same whichever resolution
	 *                           is used.
	 *
	 * @return Constants ready to bind. `xOceanMeshLevelBase` still has to be
	 *         set per draw: 0 for the centre patch, 1 for the ring instances.
	 */
	OceanCB GetOceanCB(
		const Ocean::OceanParameters& params,
		const XMFLOAT3& viewerPosition,
		const uint32_t cellsPerSide
	)
	{
		OceanCB cb = {};

		cb.xOceanMeshCenter = XMFLOAT2(
			std::round(viewerPosition.x / MESH_SNAP_QUANTUM)
				* MESH_SNAP_QUANTUM,
			std::round(viewerPosition.z / MESH_SNAP_QUANTUM)
				* MESH_SNAP_QUANTUM
		);
		cb.xOceanMeshCellSize = MESH_BASE_CELL_SIZE
			* float(MESH_CELLS_PER_SIDE) / float(cellsPerSide);
		cb.xOceanMeshVertsPerSide = cellsPerSide + 1;
		cb.xOceanMeshMorphCells = std::max(1u,
			MESH_MORPH_CELLS * cellsPerSide / MESH_CELLS_PER_SIDE);
		cb.xOceanMeshLevelBase = 0;
		cb.xOceanMeshLevelCount = MESH_LEVEL_COUNT;
		// Cells reach a given size this much sooner than on the main mesh, so
		// the waves have to flatten this much sooner too. Exactly 1 for the
		// main mesh, which is what keeps the surface and the tests against it
		// reading one band.
		cb.xOceanMeshFadeScale =
			float(MESH_CELLS_PER_SIDE) / float(cellsPerSide);

		uint32_t actual_dim = params.dmap_dim;
		uint32_t input_width = actual_dim + 4;
		uint32_t output_width = actual_dim;
		uint32_t output_height = actual_dim;
		// Offset to the second packed FFT field. Must equal one FFT slice stride
		// (actual_dim * actual_dim), which coincides with the hardcoded 512*512
		// stride of fft_512x512_c2c only at actual_dim == 512.
		uint32_t second_field_offset = actual_dim * actual_dim;
		cb.xOceanActualDim = actual_dim;
		cb.xOceanInWidth = input_width;
		cb.xOceanOutWidth = output_width;
		cb.xOceanOutHeight = output_height;
		cb.xOceanSecondFieldOffset = second_field_offset;

		cb.xOceanTimeScale = params.time_scale;
		cb.xOceanChoppyScale = params.choppy_scale;
		cb.xOceanGridLen = params.dmap_dim / params.patch_length;

		cb.xOceanWaterColor = params.waterColor;
		cb.xOceanExtinctionColor = XMFLOAT4(1 - params.extinctionColor.x, 1 - params.extinctionColor.y, 1 - params.extinctionColor.z, 1);
		cb.xOceanTexelLength = params.patch_length / params.dmap_dim;
		cb.xOceanPatchSizeRecip = 1.0f / params.patch_length;
		cb.xOceanMapHalfTexel = 0.5f / params.dmap_dim;
		cb.xOceanWaterHeight = params.waterHeight;
		cb.xOceanWaveAmplitude = params.wave_amplitude;

		return cb;
	}

	void Ocean::UpdateDisplacementMap(CommandList cmd) const
	{
		GraphicsDevice* device = wi::graphics::GetDevice();

		device->EventBegin("Ocean Simulation", cmd);

		// The simulation is camera independent, so the mesh fields this fills
		// in are simply unused here.
		const OceanCB cb =
			GetOceanCB(params, XMFLOAT3(0, 0, 0), MESH_CELLS_PER_SIDE);

		device->BindDynamicConstantBuffer(cb, CB_GETBINDSLOT(OceanCB), cmd);

		// ---------------------------- H(0) -> H(t), D(x, t), D(y, t) --------------------------------

		device->BindComputeShader(&updateSpectrumCS, cmd);

		// Buffers
		const GPUResource* cs0_srvs[2] = {
			&buffer_Float2_H0,
			&buffer_Float_Omega
		};
		device->BindResources(cs0_srvs, 0, arraysize(cs0_srvs), cmd);

		const GPUResource* cs0_uavs[1] = { &buffer_Float2_Ht };
		device->BindUAVs(cs0_uavs, 0, arraysize(cs0_uavs), cmd);

		uint32_t group_count_x = (params.dmap_dim + OCEAN_COMPUTE_TILESIZE - 1) / OCEAN_COMPUTE_TILESIZE;
		uint32_t group_count_y = (params.dmap_dim + OCEAN_COMPUTE_TILESIZE - 1) / OCEAN_COMPUTE_TILESIZE;
		device->Dispatch(group_count_x, group_count_y, 1, cmd);
		device->Barrier(GPUBarrier::Memory(), cmd);



		// ------------------------------------ Perform FFT -------------------------------------------
		wi::fftgenerator::fft_512x512_c2c(buffer_Float_Dxyz, buffer_Float_Dxyz, buffer_Float2_Ht, cmd);



		// The FFT binds its own constant buffer to the same slot on some backends,
		// so rebind the ocean constants before the displacement/gradient passes.
		device->BindDynamicConstantBuffer(cb, CB_GETBINDSLOT(OceanCB), cmd);

		GPUBarrier barriers[] = {
			GPUBarrier::Image(&displacementMap, displacementMap.desc.layout, ResourceState::UNORDERED_ACCESS),
			GPUBarrier::Image(&gradientMap, gradientMap.desc.layout, ResourceState::UNORDERED_ACCESS),
		};
		device->Barrier(barriers, arraysize(barriers), cmd);

		// Update displacement map:
		device->BindComputeShader(&updateDisplacementMapCS, cmd);
		device->BindUAV(&displacementMap, 0, cmd);
		device->BindResource(&buffer_Float_Dxyz, 0, cmd);
		device->Dispatch(
			(params.dmap_dim + OCEAN_COMPUTE_TILESIZE - 1) / OCEAN_COMPUTE_TILESIZE,
			(params.dmap_dim + OCEAN_COMPUTE_TILESIZE - 1) / OCEAN_COMPUTE_TILESIZE,
			1,
			cmd
		);
		device->Barrier(GPUBarrier::Image(&displacementMap, ResourceState::UNORDERED_ACCESS, displacementMap.desc.layout), cmd);

		// Update gradient map:
		device->BindComputeShader(&updateGradientFoldingCS, cmd);
		device->BindUAV(&gradientMap, 0, cmd);
		device->BindResource(&displacementMap, 0, cmd);
		device->Dispatch(
			(params.dmap_dim + OCEAN_COMPUTE_TILESIZE - 1) / OCEAN_COMPUTE_TILESIZE,
			(params.dmap_dim + OCEAN_COMPUTE_TILESIZE - 1) / OCEAN_COMPUTE_TILESIZE,
			1,
			cmd
		);
		device->Barrier(GPUBarrier::Image(&gradientMap, ResourceState::UNORDERED_ACCESS, gradientMap.desc.layout), cmd);

		wi::renderer::GenerateMipChain(gradientMap, wi::renderer::MIPGENFILTER_LINEAR, cmd);

		device->EventEnd(cmd);
	}

	void Ocean::CopyDisplacementMapReadback(wi::graphics::CommandList cmd) const
	{
		// Skip the full-texture GPU->CPU copy unless ocean height was queried on the
		// CPU recently (see GetDisplacedPosition). The readback is only needed for
		// CPU-side queries like buoyancy, so when nothing reads it this saves a
		// per-frame copy of the whole RGBA32F displacement map.
		if (displacement_readback_request == 0)
			return;
		displacement_readback_request--;

		GraphicsDevice* device = wi::graphics::GetDevice();
		device->EventBegin("Ocean Readback Copy", cmd);
		device->CopyResource(&displacementMap_readback[displacement_readback_index], &displacementMap, cmd);
		displacement_readback_valid[displacement_readback_index] = true;
		displacement_readback_index = (displacement_readback_index + 1) % device->GetBufferCount();
		device->EventEnd(cmd);
	}

	void Ocean::RenderForOcclusionTest(const CameraComponent& camera, CommandList cmd) const
	{
		GraphicsDevice* device = wi::graphics::GetDevice();

		device->EventBegin("Ocean Occlusion Test", cmd);

		device->BindPipelineState(&PSO_occlusionTest, cmd);

		device->BindResource(&displacementMap, 0, cmd);
		device->BindResource(&perlinTex, 2, cmd);

		OceanCB cb = GetOceanCB(params, camera.Eye, MESH_CELLS_PER_SIDE_COARSE);
		DrawClipmap(cb, true, 1, cmd);

		device->EventEnd(cmd);
	}

	void Ocean::RenderForCubemap(const XMFLOAT3& viewerPosition, CommandList cmd) const
	{
		GraphicsDevice* device = wi::graphics::GetDevice();

		device->EventBegin("Ocean Rendering into Cubemap", cmd);

		device->BindPipelineState(&PSO_envmap, cmd);

		device->BindResource(&displacementMap, 0, cmd);
		device->BindResource(&gradientMap, 1, cmd);
		device->BindResource(&perlinTex, 2, cmd);

		OceanCB cb =
			GetOceanCB(params, viewerPosition, MESH_CELLS_PER_SIDE_COARSE);
		DrawClipmap(cb, true, 6, cmd); // 6 instances, one per cube side

		device->EventEnd(cmd);
	}

	void Ocean::RenderForShadowmap(const XMFLOAT3& viewerPosition, CommandList cmd) const
	{
		GraphicsDevice* device = wi::graphics::GetDevice();

		device->EventBegin("Ocean Rendering into shadow map", cmd);

		device->BindPipelineState(&PSO_shadowmap, cmd);

		device->BindResource(&displacementMap, 0, cmd);
		device->BindResource(&gradientMap, 1, cmd);
		device->BindResource(&perlinTex, 2, cmd);

		OceanCB cb =
			GetOceanCB(params, viewerPosition, MESH_CELLS_PER_SIDE_COARSE);
		DrawClipmap(cb, true, 1, cmd);

		device->EventEnd(cmd);
	}

	void Ocean::Render(const CameraComponent& camera, CommandList cmd) const
	{
		GraphicsDevice* device = wi::graphics::GetDevice();

		device->EventBegin("Ocean Rendering", cmd);

		const bool wire = wi::renderer::IsWireRender();

		if (wire)
		{
			device->BindPipelineState(&PSO_wire, cmd);
		}
		else
		{
			device->BindPipelineState(&PSO, cmd);
		}

		device->BindResource(&displacementMap, 0, cmd);
		device->BindResource(&gradientMap, 1, cmd);
		device->BindResource(&perlinTex, 2, cmd);

		OceanCB cb = GetOceanCB(params, camera.Eye, MESH_CELLS_PER_SIDE);
		DrawClipmap(cb, false, 1, cmd);

		device->EventEnd(cmd);
	}


	void Ocean::Initialize()
	{
		wi::Timer timer;

		RasterizerState ras_desc;
		ras_desc.fill_mode = FillMode::SOLID;
		ras_desc.cull_mode = CullMode::NONE;
		ras_desc.front_counter_clockwise = false;
		ras_desc.depth_bias = 0;
		ras_desc.slope_scaled_depth_bias = 0.0f;
		ras_desc.depth_bias_clamp = 0.0f;
		ras_desc.depth_clip_enable = true;
		ras_desc.multisample_enable = false;
		ras_desc.antialiased_line_enable = false;

		rasterizerState = ras_desc;

		ras_desc.fill_mode = FillMode::WIREFRAME;
		wireRS = ras_desc;

		DepthStencilState depth_desc;
		depth_desc.depth_enable = true;
		depth_desc.depth_write_mask = DepthWriteMask::ALL;
		depth_desc.depth_func = ComparisonFunc::GREATER;
		depth_desc.stencil_enable = false;
		depthStencilState = depth_desc;

		depth_desc.depth_write_mask = DepthWriteMask::ZERO;
		depthStencilState_occlusionTest = depth_desc;
		depthStencilState_shadowmap = depth_desc;

		BlendState blend_desc;
		blend_desc.alpha_to_coverage_enable = false;
		blend_desc.independent_blend_enable = false;
		blend_desc.render_target[0].blend_enable = true;
		blend_desc.render_target[0].src_blend = Blend::SRC_ALPHA;
		blend_desc.render_target[0].dest_blend = Blend::INV_SRC_ALPHA;
		blend_desc.render_target[0].blend_op = BlendOp::ADD;
		blend_desc.render_target[0].src_blend_alpha = Blend::ONE;
		blend_desc.render_target[0].dest_blend_alpha = Blend::ZERO;
		blend_desc.render_target[0].blend_op_alpha = BlendOp::ADD;
		blend_desc.render_target[0].render_target_write_mask = ColorWrite::ENABLE_ALL;
		blendState = blend_desc;

		blend_desc.render_target[0].blend_enable = false;
		blendState_occlusionTest = blend_desc;

		blend_desc.render_target[0].src_blend = Blend::ZERO;
		blend_desc.render_target[0].dest_blend = Blend::SRC_COLOR;
		blend_desc.render_target[0].blend_op = BlendOp::ADD;
		blend_desc.render_target[0].src_blend_alpha = Blend::ONE;
		blend_desc.render_target[0].dest_blend_alpha = Blend::ONE;
		blend_desc.render_target[0].blend_op_alpha = BlendOp::MAX;
		blend_desc.render_target[0].blend_enable = true;
		blend_desc.render_target[0].render_target_write_mask = ColorWrite::ENABLE_ALL;
		blend_desc.alpha_to_coverage_enable = false;
		blend_desc.independent_blend_enable = false;
		blendState_shadowmap = blend_desc;


		static wi::eventhandler::Handle handle = wi::eventhandler::Subscribe(wi::eventhandler::EVENT_RELOAD_SHADERS, [](uint64_t userdata) { LoadShaders(); wi::fftgenerator::LoadShaders(); });

		LoadShaders();
		wi::fftgenerator::LoadShaders();

		// Connectivity only - the vertex shader generates the positions - so
		// these depend on nothing per-ocean and are built once here rather
		// than lazily behind a lock on every draw.
		CreateClipmapIndexBuffer(
			MESH_CELLS_PER_SIDE, false,
			indexBuffer_patch, "Ocean::indexBuffer_patch");
		CreateClipmapIndexBuffer(
			MESH_CELLS_PER_SIDE, true,
			indexBuffer_ring, "Ocean::indexBuffer_ring");
		CreateClipmapIndexBuffer(
			MESH_CELLS_PER_SIDE_COARSE, false,
			indexBuffer_patch_coarse, "Ocean::indexBuffer_patch_coarse");
		CreateClipmapIndexBuffer(
			MESH_CELLS_PER_SIDE_COARSE, true,
			indexBuffer_ring_coarse, "Ocean::indexBuffer_ring_coarse");

		TextureDesc desc;
		desc.bind_flags = BindFlag::SHADER_RESOURCE;
		desc.width = 64;
		desc.height = 64;
		desc.mip_levels = 7;
		desc.format = Format::R16G16B16A16_FLOAT;
		const uint32_t data_stride = GetFormatStride(desc.format);
		const uint32_t block_size = GetFormatBlockSize(desc.format);
		SubresourceData initdata[7] = {};
		const uint8_t* src = perlin;
		for (uint32_t mip = 0; mip < desc.mip_levels; ++mip)
		{
			const uint32_t num_blocks_x = std::max(1u, desc.width >> mip) / block_size;
			const uint32_t num_blocks_y = std::max(1u, desc.height >> mip) / block_size;
			initdata[mip].data_ptr = src;
			initdata[mip].row_pitch = num_blocks_x * data_stride;
			src += num_blocks_x * num_blocks_y * data_stride;
		}
		GraphicsDevice* device = GetDevice();
		device->CreateTexture(&desc, initdata, &perlinTex);
		device->SetName(&perlinTex, "ocean_internal::perlinTex");

		wilog("wi::Ocean Initialized (%d ms)", (int)std::round(timer.elapsed()));
	}

	const Texture* Ocean::getDisplacementMap() const
	{
		return &displacementMap;
	}

	const Texture* Ocean::getGradientMap() const
	{
		return &gradientMap;
	}

	XMFLOAT2 Ocean::GetDisplacementFadeBand() const noexcept
	{
		// A level's outer boundary sits at (cells / 2) * cellSize, so the
		// cells at distance d are 2d / cells across. Inverting that turns a
		// cell size into the distance at which the mesh reaches it.
		const float distancePerCellSize = float(MESH_CELLS_PER_SIDE) * 0.5F;

		return XMFLOAT2(
			params.patch_length / MESH_FADE_START_CELLS_PER_PATCH
				* distancePerCellSize,
			params.patch_length / MESH_FADE_END_CELLS_PER_PATCH
				* distancePerCellSize
		);
	}

	const wi::primitive::AABB Ocean::GetAABB(const XMFLOAT3& camera_pos) const
	{
		wi::primitive::AABB aabb;
		aabb.createFromHalfWidth(XMFLOAT3(camera_pos.x, params.waterHeight, camera_pos.z), XMFLOAT3(1000, 10, 1000));
		return aabb;
	}

}
