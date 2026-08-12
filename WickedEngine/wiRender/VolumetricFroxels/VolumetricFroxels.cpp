/**
 * @file
 * Builds the volumetric froxel volume - see `VolumetricFroxels.h`.
 */
#include "VolumetricFroxels.h"

#include "wiEventHandler.h"
#include "wiProfiler.h"
#include "wiRenderer.h"
#include "wiScene_Components.h"
#include "shaders/ShaderInterop_Postprocess.h"
#include "shaders/ShaderInterop_VolumetricFroxels.h"

#include <algorithm>
#include <cstdint>

using namespace wi::graphics;
using namespace wi::render;

namespace
{
	/** Fills every cell with what the medium there scatters towards the eye. */
	Shader injectCS;

	/** Gathers those along each view ray into what a fragment reads. */
	Shader integrateCS;

	/** Adds the result to a render target against the opaque depth. */
	Shader applyCS;

	/**
	 * Loads the compute shaders.
	 *
	 * Separate from `VolumetricFroxels::Initialize` because it also runs on a
	 * shader reload, where the event subscription must not be repeated.
	 */
	void LoadShaders()
	{
		wi::renderer::LoadShader(
			ShaderStage::CS, injectCS, "volumetricFroxel_injectCS.cso");
		wi::renderer::LoadShader(
			ShaderStage::CS, integrateCS, "volumetricFroxel_integrateCS.cso");
		wi::renderer::LoadShader(
			ShaderStage::CS, applyCS, "volumetricFroxel_applyCS.cso");
	}

	/**
	 * Thread groups needed to cover an extent.
	 *
	 * @param[in] extent - Cells to cover.
	 * @param[in] blockSize - Cells one group covers.
	 *
	 * @return Group count, rounded up so the last partial group still runs.
	 */
	[[nodiscard]] uint32_t GroupCount(
		const uint32_t extent, const uint32_t blockSize
	) noexcept
	{
		return (extent + blockSize - 1) / blockSize;
	}
} // namespace

/*
################################################################################
Public
################################################################################
*/

// Getters
//==============================================================================

bool VolumetricFroxels::IsValid() const noexcept
{
	return integratedVolume.IsValid();
}

const Texture* VolumetricFroxels::GetVolume() const noexcept
{
	return &integratedVolume;
}

float VolumetricFroxels::GetRange() const noexcept
{
	return range;
}

// Methods
//==============================================================================

void VolumetricFroxels::Initialize()
{
	static wi::eventhandler::Handle handle = wi::eventhandler::Subscribe(
		wi::eventhandler::EVENT_RELOAD_SHADERS,
		[](uint64_t userdata) { LoadShaders(); });

	LoadShaders();
}

void VolumetricFroxels::Create(const float range)
{
	// Nothing about the volume's shape depends on the range, so a change to it
	// only has to be recorded - the grid it applies to is the same one.
	if (IsValid())
	{
		this->range = range;
		return;
	}

	this->range = range;
	ResetFrame();

	GraphicsDevice* device = wi::graphics::GetDevice();

	TextureDesc desc;
	desc.type = TextureDesc::Type::TEXTURE_3D;
	desc.width = VOLUMETRIC_FROXEL_WIDTH;
	desc.height = VOLUMETRIC_FROXEL_HEIGHT;
	desc.depth = VOLUMETRIC_FROXEL_SLICES;
	// Radiance is never negative and never needs an alpha here, so eleven bits
	// of mantissa per channel is ample and halves what a volume this size costs
	// against a 16-bit format.
	desc.format = Format::R11G11B10_FLOAT;
	desc.bind_flags = BindFlag::SHADER_RESOURCE | BindFlag::UNORDERED_ACCESS;
	desc.layout = ResourceState::SHADER_RESOURCE_COMPUTE;

	for (uint32_t i = 0; i < arraysize(injectVolume); ++i)
	{
		device->CreateTextureZeroed(&desc, &injectVolume[i]);
		device->SetName(
			&injectVolume[i], "VolumetricFroxels::injectVolume");
	}

	device->CreateTextureZeroed(&desc, &integratedVolume);
	device->SetName(&integratedVolume, "VolumetricFroxels::integratedVolume");
}

void VolumetricFroxels::Reset() noexcept
{
	for (uint32_t i = 0; i < arraysize(injectVolume); ++i)
	{
		injectVolume[i] = {};
	}

	integratedVolume = {};
	range = 0.0F;
	ResetFrame();
}

void VolumetricFroxels::ResetFrame() const noexcept
{
	frame = -1;
}

void VolumetricFroxels::AdvanceFrame() const noexcept
{
	frame++;
}

void VolumetricFroxels::Build(
	const wi::scene::CameraComponent& camera, const CommandList cmd
) const
{
	if (!IsValid())
	{
		return;
	}

	GraphicsDevice* device = wi::graphics::GetDevice();

	device->EventBegin("Volumetric Froxels", cmd);
	auto profilerRange =
		wi::profiler::BeginRangeGPU("Volumetric Froxels", cmd);

	// Bound here rather than assumed of the caller. The injection reads the
	// light entities and the shadow atlas, and without these it finds neither -
	// which looks exactly like a scene with no volumetric lights in it rather
	// than like a missing binding.
	wi::renderer::BindCommonResources(cmd);

	VolumetricFroxelPush push;
	push.range = range;
	push.frame = (uint32_t)std::max(frame, 0);
	push.padding0 = 0;
	push.padding1 = 0;

	const Texture& injectTarget = injectVolume[GetInjectOutputIndex()];

	// Injection
	//==========================================================================
	{
		device->EventBegin("Inject", cmd);

		device->Barrier(
			GPUBarrier::Image(
				&injectTarget,
				injectTarget.desc.layout,
				ResourceState::UNORDERED_ACCESS),
			cmd);

		device->BindComputeShader(&injectCS, cmd);
		device->PushConstants(&push, sizeof(push), cmd);
		device->BindUAV(&injectTarget, 0, cmd);

		device->Dispatch(
			GroupCount(
				VOLUMETRIC_FROXEL_WIDTH, VOLUMETRIC_FROXEL_BLOCKSIZE_X),
			GroupCount(
				VOLUMETRIC_FROXEL_HEIGHT, VOLUMETRIC_FROXEL_BLOCKSIZE_Y),
			GroupCount(
				VOLUMETRIC_FROXEL_SLICES, VOLUMETRIC_FROXEL_BLOCKSIZE_Z),
			cmd);

		device->Barrier(
			GPUBarrier::Image(
				&injectTarget,
				ResourceState::UNORDERED_ACCESS,
				injectTarget.desc.layout),
			cmd);

		device->EventEnd(cmd);
	}

	// Integration
	//==========================================================================
	{
		device->EventBegin("Integrate", cmd);

		device->Barrier(
			GPUBarrier::Image(
				&integratedVolume,
				integratedVolume.desc.layout,
				ResourceState::UNORDERED_ACCESS),
			cmd);

		device->BindComputeShader(&integrateCS, cmd);
		device->PushConstants(&push, sizeof(push), cmd);
		device->BindResource(&injectTarget, 0, cmd);
		device->BindUAV(&integratedVolume, 0, cmd);

		// One thread per column, so the slice axis is not dispatched over.
		device->Dispatch(
			GroupCount(
				VOLUMETRIC_FROXEL_WIDTH, VOLUMETRIC_FROXEL_BLOCKSIZE_X),
			GroupCount(
				VOLUMETRIC_FROXEL_HEIGHT, VOLUMETRIC_FROXEL_BLOCKSIZE_Y),
			1,
			cmd);

		device->Barrier(
			GPUBarrier::Image(
				&integratedVolume,
				ResourceState::UNORDERED_ACCESS,
				integratedVolume.desc.layout),
			cmd);

		device->EventEnd(cmd);
	}

	wi::profiler::EndRange(profilerRange);
	device->EventEnd(cmd);
}

void VolumetricFroxels::ApplyScreenSpace(
	const Texture& target, const CommandList cmd
) const
{
	if (!IsValid())
	{
		return;
	}

	GraphicsDevice* device = wi::graphics::GetDevice();

	device->EventBegin("Volumetric Froxels Apply", cmd);

	wi::renderer::BindCommonResources(cmd);

	VolumetricFroxelPush push;
	push.range = range;
	push.frame = (uint32_t)std::max(frame, 0);
	push.padding0 = 0;
	push.padding1 = 0;

	device->Barrier(
		GPUBarrier::Image(
			&target, target.desc.layout, ResourceState::UNORDERED_ACCESS),
		cmd);

	device->BindComputeShader(&applyCS, cmd);
	device->PushConstants(&push, sizeof(push), cmd);
	device->BindResource(&integratedVolume, 0, cmd);
	device->BindUAV(&target, 0, cmd);

	device->Dispatch(
		GroupCount(target.desc.width, POSTPROCESS_BLOCKSIZE),
		GroupCount(target.desc.height, POSTPROCESS_BLOCKSIZE),
		1,
		cmd);

	device->Barrier(
		GPUBarrier::Image(
			&target, ResourceState::UNORDERED_ACCESS, target.desc.layout),
		cmd);

	device->EventEnd(cmd);
}

/*
################################################################################
Private
################################################################################
*/

// Getters
//==============================================================================

int VolumetricFroxels::GetInjectOutputIndex() const noexcept
{
	// Clamped rather than taken modulo directly: before the first frame the
	// counter is negative, and a negative modulo would index outside the pair.
	return std::max(frame, 0) % 2;
}

int VolumetricFroxels::GetInjectHistoryIndex() const noexcept
{
	return 1 - GetInjectOutputIndex();
}
