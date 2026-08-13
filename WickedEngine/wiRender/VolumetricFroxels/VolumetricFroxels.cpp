/**
 * @file
 * Builds the volumetric froxel volume - see `VolumetricFroxels.h`.
 */
#include "VolumetricFroxels.h"

#include "wiEventHandler.h"
#include "wiProfiler.h"
#include "wiRenderer.h"
#include "wiScene_Components.h"
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

const Texture* VolumetricFroxels::GetTail() const noexcept
{
	return &tailTexture;
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
	desc.bind_flags = BindFlag::SHADER_RESOURCE | BindFlag::UNORDERED_ACCESS;

	// Radiance is never negative and never needs an alpha here, so eleven bits
	// of mantissa per channel is ample and halves what a volume this size costs
	// against a 16-bit format.
	//
	// The injection pair is read back into itself, and a sixteen bit format was
	// measured there on the theory that a twentieth-of-the-way blend step falls
	// below what eleven bits can resolve. It changed nothing on screen, so the
	// precision floor is not what limits this.
	desc.format = Format::R11G11B10_FLOAT;

	// Nothing outside the integration pass ever reads these, and that is a
	// compute pass.
	desc.layout = ResourceState::SHADER_RESOURCE_COMPUTE;

	// Left with undefined contents deliberately. Both build passes write every
	// cell they own before anything reads one, so zeroing them buys nothing and
	// costs a copy of the whole volume through the staging allocator.
	for (uint32_t i = 0; i < arraysize(injectVolume); ++i)
	{
		device->CreateTexture(&desc, nullptr, &injectVolume[i]);
		device->SetName(
			&injectVolume[i], "VolumetricFroxels::injectVolume");
	}

	// Read by the fragment shaders, so it must be visible to the pixel stage
	// too - the compute-only state is exactly right for the pair above and
	// silently wrong here.
	desc.layout = ResourceState::SHADER_RESOURCE;

	device->CreateTexture(&desc, nullptr, &integratedVolume);
	device->SetName(&integratedVolume, "VolumetricFroxels::integratedVolume");

	// One column of the screen, not one cell of the volume, and sixteen bits
	// because the alpha carries an extinction rather than a colour.
	TextureDesc tailDesc;
	tailDesc.width = VOLUMETRIC_FROXEL_WIDTH;
	tailDesc.height = VOLUMETRIC_FROXEL_HEIGHT;
	tailDesc.format = Format::R16G16B16A16_FLOAT;
	tailDesc.bind_flags =
		BindFlag::SHADER_RESOURCE | BindFlag::UNORDERED_ACCESS;
	tailDesc.layout = ResourceState::SHADER_RESOURCE;

	device->CreateTexture(&tailDesc, nullptr, &tailTexture);
	device->SetName(&tailTexture, "VolumetricFroxels::tailTexture");
}

void VolumetricFroxels::Reset() noexcept
{
	for (uint32_t i = 0; i < arraysize(injectVolume); ++i)
	{
		injectVolume[i] = {};
	}

	integratedVolume = {};
	tailTexture = {};
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
	const wi::scene::CameraComponent& camera,
	const wi::scene::CameraComponent& cameraPrevious,
	const wi::scene::CameraComponent& cameraReflection,
	const CommandList cmd
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

	// The grid must stand still between frames while the sample point inside
	// each cell moves, so the sub-pixel jitter that temporal anti-aliasing puts
	// on the frustum is taken off both cameras for the build. Left in, it would
	// shift the whole volume by a fraction of a cell every frame, and the
	// accumulation would be averaging the grid's own movement along with the
	// light.
	wi::scene::CameraComponent cameraStill = camera;
	cameraStill.jitter = XMFLOAT2(0, 0);
	cameraStill.UpdateCamera();

	wi::scene::CameraComponent cameraPreviousStill = cameraPrevious;
	cameraPreviousStill.jitter = XMFLOAT2(0, 0);
	cameraPreviousStill.UpdateCamera();

	wi::renderer::BindCameraCB(
		cameraStill, cameraPreviousStill, cameraReflection, cmd);

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
	push.previousEye = cameraPrevious.Eye;
	push.padding2 = 0;

	const Texture& injectTarget = injectVolume[GetInjectOutputIndex()];
	const Texture& injectHistory = injectVolume[GetInjectHistoryIndex()];

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
		device->BindResource(&injectHistory, 0, cmd);
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

		const GPUBarrier toWrite[] = {
			GPUBarrier::Image(
				&integratedVolume,
				integratedVolume.desc.layout,
				ResourceState::UNORDERED_ACCESS),
			GPUBarrier::Image(
				&tailTexture,
				tailTexture.desc.layout,
				ResourceState::UNORDERED_ACCESS),
		};
		device->Barrier(toWrite, arraysize(toWrite), cmd);

		device->BindComputeShader(&integrateCS, cmd);
		device->PushConstants(&push, sizeof(push), cmd);
		device->BindResource(&injectTarget, 0, cmd);
		device->BindUAV(&integratedVolume, 0, cmd);
		device->BindUAV(&tailTexture, 1, cmd);

		// One thread per column, so the slice axis is not dispatched over.
		device->Dispatch(
			GroupCount(
				VOLUMETRIC_FROXEL_WIDTH, VOLUMETRIC_FROXEL_BLOCKSIZE_X),
			GroupCount(
				VOLUMETRIC_FROXEL_HEIGHT, VOLUMETRIC_FROXEL_BLOCKSIZE_Y),
			1,
			cmd);

		const GPUBarrier toRead[] = {
			GPUBarrier::Image(
				&integratedVolume,
				ResourceState::UNORDERED_ACCESS,
				integratedVolume.desc.layout),
			GPUBarrier::Image(
				&tailTexture,
				ResourceState::UNORDERED_ACCESS,
				tailTexture.desc.layout),
		};
		device->Barrier(toRead, arraysize(toRead), cmd);

		device->EventEnd(cmd);
	}

	// Restored, because the rest of this command list was recorded against the
	// jittered camera and everything after this point expects it back.
	wi::renderer::BindCameraCB(camera, cameraPrevious, cameraReflection, cmd);

	wi::profiler::EndRange(profilerRange);
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
