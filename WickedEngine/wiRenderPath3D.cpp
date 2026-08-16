#include "wiRenderPath3D.h"
#include "wiRenderer.h"
#include "wiImage.h"
#include "wiHelper.h"
#include "wiTextureHelper.h"
#include "wiProfiler.h"
#include "wiScene/environment/UnderwaterParticles/UnderwaterParticles.h"

using namespace wi::graphics;
using namespace wi::enums;
using namespace wi::scene;
using namespace wi::ecs;

namespace wi
{
	static constexpr float foreground_depth_range = 0.01f;

	void RenderPath3D::DeleteGPUResources()
	{
		rtMain = {};
		rtMain_render = {};
		rtPrimitiveID = {};
		rtPrimitiveID_render = {};
		rtVelocity = {};
		rtReflection = {};
		rtReflection_render = {};
		rtRaytracedDiffuse = {};
		rtSSR = {};
		rtSceneCopy = {};
		rtSceneCopy_tmp = {};
		rtWaterRipple = {};
		rtParticleDistortion_render = {};
		rtParticleDistortion = {};
		rtBloom = {};
		rtBloom_tmp = {};
		rtAO = {};
		rtShadow = {};
		rtSun[0] = {};
		rtSun[1] = {};
		rtSun[2] = {};
		rtSun_resolved = {};
		rtGUIBlurredBackground[0] = {};
		rtGUIBlurredBackground[1] = {};
		rtGUIBlurredBackground[2] = {};
		rtShadingRate = {};
		rtFSR[0] = {};
		rtFSR[1] = {};
		rtOutlineSource = {};

		rtPostprocess = {};

		depthBuffer_Main = {};
		depthBuffer_Ocean = {};
		depthBuffer_Copy = {};
		depthBuffer_Copy1 = {};
		depthBuffer_Reflection = {};
		depthBuffer_Reflection_render = {};
		reprojectedDepth = {};

		debugUAV = {};
		tiledLightResources = {};
		tiledLightResources_planarReflection = {};
		luminanceResources = {};
		ssaoResources = {};
		msaoResources = {};
		rtaoResources = {};
		rtdiffuseResources = {};
		rtreflectionResources = {};
		ssrResources = {};
		rtshadowResources = {};
		screenspaceshadowResources = {};
		depthoffieldResources = {};
		motionblurResources = {};
		aerialperspectiveResources = {};
		aerialperspectiveResources_reflection = {};
		volumetriccloudResources = {};
		bloomResources = {};
		surfelGIResources = {};
		temporalAAResources = {};
		visibilityResources = {};
		fsr2Resources = {};
		vxgiResources = {};
		meshblendResources = {};

		RenderPath2D::DeleteGPUResources();
	}

	void RenderPath3D::ResizeBuffers()
	{
		first_frame = true;
		DeleteGPUResources();

		GraphicsDevice* device = wi::graphics::GetDevice();

		XMUINT2 internalResolution = GetInternalResolution();
		camera->width = (float)internalResolution.x;
		camera->height = (float)internalResolution.y;

		// Render targets:

		{
			TextureDesc desc;
			desc.format = wi::renderer::format_rendertarget_main;
			desc.bind_flags = BindFlag::RENDER_TARGET | BindFlag::SHADER_RESOURCE | BindFlag::UNORDERED_ACCESS;
			desc.width = internalResolution.x;
			desc.height = internalResolution.y;
			desc.sample_count = 1;
			device->CreateTexture(&desc, nullptr, &rtMain);
			device->SetName(&rtMain, "renderpath3D.rtMain");

			if (getMSAASampleCount() > 1)
			{
				desc.sample_count = getMSAASampleCount();
				desc.bind_flags = BindFlag::RENDER_TARGET | BindFlag::SHADER_RESOURCE;

				device->CreateTexture(&desc, nullptr, &rtMain_render);
				device->SetName(&rtMain_render, "renderpath3D.rtMain_render");
				
				// Note: graphics API can downgrade sample count for last supported value, this will be reflected in the renderpath setting too
				msaaSampleCount = std::min(msaaSampleCount, rtMain_render.desc.sample_count);
			}
			else
			{
				rtMain_render = rtMain;
			}
		}
		{
			TextureDesc desc;
			desc.format = wi::renderer::format_idbuffer;
			desc.bind_flags = BindFlag::RENDER_TARGET | BindFlag::SHADER_RESOURCE;
			if (getMSAASampleCount() > 1)
			{
				desc.bind_flags |= BindFlag::UNORDERED_ACCESS;
			}
			desc.width = internalResolution.x;
			desc.height = internalResolution.y;
			desc.sample_count = 1;
			desc.layout = ResourceState::SHADER_RESOURCE_COMPUTE;
			desc.misc_flags = ResourceMiscFlag::ALIASING_TEXTURE_RT_DS;
			device->CreateTexture(&desc, nullptr, &rtPrimitiveID);
			device->SetName(&rtPrimitiveID, "renderpath3D.rtPrimitiveID");

			if (getMSAASampleCount() > 1)
			{
				desc.sample_count = getMSAASampleCount();
				desc.bind_flags = BindFlag::RENDER_TARGET | BindFlag::SHADER_RESOURCE;
				desc.misc_flags = ResourceMiscFlag::NONE;
				device->CreateTexture(&desc, nullptr, &rtPrimitiveID_render);
				device->SetName(&rtPrimitiveID_render, "renderpath3D.rtPrimitiveID_render");
			}
			else
			{
				rtPrimitiveID_render = rtPrimitiveID;
			}
		}
		{
			TextureDesc desc;
			desc.bind_flags = BindFlag::RENDER_TARGET | BindFlag::SHADER_RESOURCE;
			desc.format = Format::R16G16_FLOAT;
			desc.width = internalResolution.x;
			desc.height = internalResolution.y;
			desc.sample_count = 1;
			desc.misc_flags = ResourceMiscFlag::ALIASING_TEXTURE_RT_DS;
			device->CreateTexture(&desc, nullptr, &rtParticleDistortion);
			device->SetName(&rtParticleDistortion, "renderpath3D.rtParticleDistortion");
			if (getMSAASampleCount() > 1)
			{
				desc.sample_count = getMSAASampleCount();
				desc.misc_flags = ResourceMiscFlag::NONE;
				device->CreateTexture(&desc, nullptr, &rtParticleDistortion_render);
				device->SetName(&rtParticleDistortion_render, "renderpath3D.rtParticleDistortion_render");
			}
			else
			{
				rtParticleDistortion_render = rtParticleDistortion;
			}
		}
		{
			TextureDesc desc;
			desc.bind_flags = BindFlag::SHADER_RESOURCE | BindFlag::UNORDERED_ACCESS;
			desc.format = wi::renderer::format_rendertarget_main;
			desc.width = internalResolution.x / 4;
			desc.height = internalResolution.y / 4;
			desc.mip_levels = std::min(8u, (uint32_t)std::log2(std::max(desc.width, desc.height)));
			device->CreateTextureZeroed(&desc, &rtSceneCopy);
			device->SetName(&rtSceneCopy, "renderpath3D.rtSceneCopy");
			desc.bind_flags = BindFlag::SHADER_RESOURCE | BindFlag::UNORDERED_ACCESS | BindFlag::RENDER_TARGET; // render target for aliasing
			device->CreateTexture(&desc, nullptr, &rtSceneCopy_tmp, &rtPrimitiveID);
			device->SetName(&rtSceneCopy_tmp, "renderpath3D.rtSceneCopy_tmp");

			device->CreateMipgenSubresources(rtSceneCopy);
			device->CreateMipgenSubresources(rtSceneCopy_tmp);
		}
		{
			TextureDesc desc;
			desc.bind_flags = BindFlag::RENDER_TARGET | BindFlag::SHADER_RESOURCE | BindFlag::UNORDERED_ACCESS;
			desc.format = wi::renderer::format_rendertarget_main;
			desc.width = internalResolution.x;
			desc.height = internalResolution.y;
			assert(ComputeTextureMemorySizeInBytes(desc) <= ComputeTextureMemorySizeInBytes(rtPrimitiveID.desc)); // Aliased check
			device->CreateTexture(&desc, nullptr, &rtPostprocess, &rtPrimitiveID); // Aliased!
			device->SetName(&rtPostprocess, "renderpath3D.rtPostprocess");
		}
		{
			TextureDesc desc;
			desc.bind_flags = BindFlag::RENDER_TARGET | BindFlag::SHADER_RESOURCE | BindFlag::UNORDERED_ACCESS;
			desc.format = Format::R10G10B10A2_UNORM;
			desc.width = internalResolution.x / 4;
			desc.height = internalResolution.y / 4;
			desc.bind_flags = BindFlag::UNORDERED_ACCESS | BindFlag::SHADER_RESOURCE;
			device->CreateTexture(&desc, nullptr, &rtGUIBlurredBackground[0]);
			device->SetName(&rtGUIBlurredBackground[0], "renderpath3D.rtGUIBlurredBackground[0]");

			desc.width /= 4;
			desc.height /= 4;
			device->CreateTexture(&desc, nullptr, &rtGUIBlurredBackground[1]);
			device->SetName(&rtGUIBlurredBackground[1], "renderpath3D.rtGUIBlurredBackground[1]");
			device->CreateTexture(&desc, nullptr, &rtGUIBlurredBackground[2]);
			device->SetName(&rtGUIBlurredBackground[2], "renderpath3D.rtGUIBlurredBackground[2]");
		}
		if (device->CheckCapability(GraphicsDeviceCapability::VARIABLE_RATE_SHADING_TIER2) &&
			wi::renderer::GetVariableRateShadingClassification())
		{
			uint32_t tileSize = device->GetVariableRateShadingTileSize();

			TextureDesc desc;
			desc.layout = ResourceState::UNORDERED_ACCESS;
			desc.bind_flags = BindFlag::UNORDERED_ACCESS | BindFlag::SHADING_RATE;
			desc.format = Format::R8_UINT;
			desc.width = (internalResolution.x + tileSize - 1) / tileSize;
			desc.height = (internalResolution.y + tileSize - 1) / tileSize;

			device->CreateTexture(&desc, nullptr, &rtShadingRate);
			device->SetName(&rtShadingRate, "renderpath3D.rtShadingRate");
		}

		// Depth buffers:
		{
			TextureDesc desc;
			desc.width = internalResolution.x;
			desc.height = internalResolution.y;

			desc.sample_count = getMSAASampleCount();
			desc.layout = ResourceState::DEPTHSTENCIL;
			desc.format = wi::renderer::format_depthbuffer_main;
			desc.bind_flags = BindFlag::DEPTH_STENCIL;
			device->CreateTexture(&desc, nullptr, &depthBuffer_Main);
			device->SetName(&depthBuffer_Main, "renderpath3D.depthBuffer_Main");

			desc.layout = ResourceState::SHADER_RESOURCE_COMPUTE;
			desc.format = Format::R32_FLOAT;
			desc.bind_flags = BindFlag::SHADER_RESOURCE | BindFlag::UNORDERED_ACCESS;
			desc.sample_count = 1;
			desc.mip_levels = 5;
			device->CreateTexture(&desc, nullptr, &depthBuffer_Copy);
			device->SetName(&depthBuffer_Copy, "renderpath3D.depthBuffer_Copy");
			device->CreateTexture(&desc, nullptr, &depthBuffer_Copy1);
			device->SetName(&depthBuffer_Copy1, "renderpath3D.depthBuffer_Copy1");

			for (uint32_t i = 0; i < depthBuffer_Copy.desc.mip_levels; ++i)
			{
				int subresource = 0;
				subresource = device->CreateSubresource(&depthBuffer_Copy, SubresourceType::SRV, 0, 1, i, 1);
				assert(subresource == i);
				subresource = device->CreateSubresource(&depthBuffer_Copy, SubresourceType::UAV, 0, 1, i, 1);
				assert(subresource == i);
				subresource = device->CreateSubresource(&depthBuffer_Copy1, SubresourceType::SRV, 0, 1, i, 1);
				assert(subresource == i);
				subresource = device->CreateSubresource(&depthBuffer_Copy1, SubresourceType::UAV, 0, 1, i, 1);
				assert(subresource == i);
			}
		}

		// Other resources:
		{
			TextureDesc desc;
			desc.width = internalResolution.x;
			desc.height = internalResolution.y;
			desc.mip_levels = 1;
			desc.array_size = 1;
			desc.format = Format::R8G8B8A8_UNORM;
			desc.sample_count = 1;
			desc.usage = Usage::DEFAULT;
			desc.bind_flags = BindFlag::SHADER_RESOURCE | BindFlag::UNORDERED_ACCESS;
			desc.layout = ResourceState::SHADER_RESOURCE;

			device->CreateTexture(&desc, nullptr, &debugUAV);
			device->SetName(&debugUAV, "renderpath3D.debugUAV");
		}
		wi::renderer::CreateTiledLightResources(tiledLightResources, internalResolution);
		wi::renderer::CreateScreenSpaceShadowResources(screenspaceshadowResources, internalResolution);

		// These can trigger resource creations if needed:
		setAO(ao);
		setSSREnabled(ssrEnabled);
		setSSGIEnabled(ssgiEnabled);
		setRaytracedReflectionsEnabled(raytracedReflectionsEnabled);
		setRaytracedDiffuseEnabled(raytracedDiffuseEnabled);
		setFSREnabled(fsrEnabled);
		setFSR2Enabled(fsr2Enabled);
		setEyeAdaptionEnabled(eyeAdaptionEnabled);
		setReflectionsEnabled(reflectionsEnabled);
		setBloomEnabled(bloomEnabled);
		setVolumeLightsEnabled(volumeLightsEnabled);
		setLightShaftsEnabled(lightShaftsEnabled);
		setOutlineEnabled(outlineEnabled);

		RenderPath2D::ResizeBuffers();
	}

	void RenderPath3D::PreUpdate()
	{
		camera_previous = *camera;
		camera_reflection_previous = camera_reflection;
	}

	void RenderPath3D::Update(float dt)
	{
		GraphicsDevice* device = wi::graphics::GetDevice();

		RenderPath2D::Update(dt);

		wi::renderer::SetShadowsEnabled(getShadowsEnabled());

		const bool hw_raytrace = device->CheckCapability(GraphicsDeviceCapability::RAYTRACING);
		if (getSceneUpdateEnabled())
		{
			if (wi::renderer::GetSurfelGIEnabled() ||
				wi::renderer::GetDDGIEnabled() ||
				(hw_raytrace && wi::renderer::GetRaytracedShadowsEnabled()) ||
				(hw_raytrace && getAO() == AO_RTAO) ||
				(hw_raytrace && getRaytracedReflectionEnabled()) ||
				(hw_raytrace && getRaytracedDiffuseEnabled()) ||
				(hw_raytrace && getUnderwaterSnellEnabled() && getUnderwaterSnellRTEnabled())
				)
			{
				scene->SetAccelerationStructureUpdateRequested(true);
			}
			scene->camera = *camera;
			scene->Update(dt * wi::renderer::GetGameSpeed());
		}
	}

	void RenderPath3D::PreRender()
	{
		GraphicsDevice* device = wi::graphics::GetDevice();

		if (rtMain_render.desc.sample_count != msaaSampleCount)
		{
			ResizeBuffers();
		}

		// Frustum culling for main camera:
		visibility_main.layerMask = getLayerMask();
		visibility_main.scene = scene;
		visibility_main.camera = camera;
		visibility_main.flags = wi::renderer::Visibility::ALLOW_EVERYTHING;
		if (!getOcclusionCullingEnabled())
		{
			visibility_main.flags &= ~wi::renderer::Visibility::ALLOW_OCCLUSION_CULLING;
		}
		wi::renderer::UpdateVisibility(visibility_main);

		if (visibility_main.planar_reflection_visible)
		{
			// Frustum culling for planar reflections:
			camera_reflection = *camera;
			camera_reflection.jitter = XMFLOAT2(0, 0);
			camera_reflection.Reflect(visibility_main.reflectionPlane);
			visibility_reflection.layerMask = getLayerMask();
			visibility_reflection.scene = scene;
			visibility_reflection.camera = &camera_reflection;
			visibility_reflection.flags =
				wi::renderer::Visibility::ALLOW_OBJECTS |
				wi::renderer::Visibility::ALLOW_EMITTERS |
				wi::renderer::Visibility::ALLOW_HAIRS |
				wi::renderer::Visibility::ALLOW_LIGHTS
				;
			wi::renderer::UpdateVisibility(visibility_reflection);
		}

		XMUINT2 internalResolution = GetInternalResolution();

		wi::renderer::UpdatePerFrameData(
			*scene,
			visibility_main,
			frameCB,
			getSceneUpdateEnabled() ? scene->dt : 0
		);

		if (getFSR2Enabled())
		{
			camera->jitter = fsr2Resources.GetJitter();
		}
		else if (wi::renderer::GetTemporalAAEnabled() && wi::renderer::GetWireframeMode() == wi::renderer::WIREFRAME_DISABLED)
		{
			const XMFLOAT4& halton = wi::math::GetHaltonSequence(wi::graphics::GetDevice()->GetFrameCount() % 256);
			camera->jitter.x = (halton.x * 2 - 1) / (float)internalResolution.x;
			camera->jitter.y = (halton.y * 2 - 1) / (float)internalResolution.y;
			if (!temporalAAResources.IsValid())
			{
				wi::renderer::CreateTemporalAAResources(temporalAAResources, internalResolution);
			}
		}
		else
		{
			camera->jitter = XMFLOAT2(0, 0);
			temporalAAResources = {};
		}

		camera_reflection.jitter = XMFLOAT2(0, 0);

		camera->UpdateCamera();
		if (visibility_main.planar_reflection_visible)
		{
			camera_reflection.UpdateCamera();
		}

		if (getAO() != AO_RTAO)
		{
			rtaoResources.frame = 0;
		}
		if (!wi::renderer::GetRaytracedShadowsEnabled())
		{
			rtshadowResources.frame = 0;
		}
		if (!getSSREnabled() && !getRaytracedReflectionEnabled())
		{
			rtSSR = {};
		}
		if (!getSSGIEnabled())
		{
			rtSSGI = {};
		}
		if (!getRaytracedDiffuseEnabled())
		{
			rtRaytracedDiffuse = {};
		}
		if (getAO() == AO_DISABLED)
		{
			rtAO = {};
		}

		if (wi::renderer::GetRaytracedShadowsEnabled() && device->CheckCapability(GraphicsDeviceCapability::RAYTRACING))
		{
			if (!rtshadowResources.denoised.IsValid())
			{
				wi::renderer::CreateRTShadowResources(rtshadowResources, internalResolution);
			}
		}
		else
		{
			rtshadowResources = {};
		}

		if (scene->weather.IsRealisticSky() && scene->weather.IsRealisticSkyAerialPerspective())
		{
			if (!aerialperspectiveResources.texture_output.IsValid())
			{
				wi::renderer::CreateAerialPerspectiveResources(aerialperspectiveResources, internalResolution);
			}
			if (getReflectionsEnabled() && depthBuffer_Reflection.IsValid())
			{
				if (!aerialperspectiveResources_reflection.texture_output.IsValid())
				{
					wi::renderer::CreateAerialPerspectiveResources(aerialperspectiveResources_reflection, XMUINT2(depthBuffer_Reflection.desc.width, depthBuffer_Reflection.desc.height));
				}
			}
			else
			{
				aerialperspectiveResources_reflection = {};
			}
		}
		else
		{
			aerialperspectiveResources = {};
		}

		if (scene->weather.IsVolumetricClouds())
		{
			if (!volumetriccloudResources.texture_cloudRender.IsValid())
			{
				wi::renderer::CreateVolumetricCloudResources(volumetriccloudResources, internalResolution);
			}
			volumetriccloudResources.AdvanceFrame();
		}
		else
		{
			volumetriccloudResources = {};
		}

		// Held to the same gate as the volumetric lights it carries, so a scene
		// with none of them pays neither the memory nor the two dispatches.
		if (getVolumeLightsEnabled() &&
			visibility_main.IsRequestedVolumetricLights())
		{
			// Never past the far plane: the last slice already carries
			// everything from the range outwards, so reach bought beyond what
			// the camera can see costs resolution for nothing.
			volumetricFroxels.Create(
				std::min(getVolumetricFroxelRange(), camera->zFarP));
			volumetricFroxels.AdvanceFrame();
		}
		else
		{
			volumetricFroxels.Reset();
		}

		// Published to the shaders that sample it. The range travels with the
		// volume because a lookup needs both: the index says where it is and the
		// range says how to turn a distance into a coordinate in it.
		camera->texture_volumetricfroxels_index = volumetricFroxels.IsValid()
			? device->GetDescriptorIndex(
				volumetricFroxels.GetVolume(), SubresourceType::SRV)
			: -1;
		wi::renderer::SetVolumetricFroxelParameters(
			volumetricFroxels.GetRange(),
			volumetricFroxels.IsValid()
				? device->GetDescriptorIndex(
					volumetricFroxels.GetTail(), SubresourceType::SRV)
				: -1);

		if (!scene->waterRipples.empty() && rtParticleDistortion.IsValid())
		{
			if (!rtWaterRipple.IsValid())
			{
				TextureDesc desc;
				desc.bind_flags = BindFlag::RENDER_TARGET | BindFlag::SHADER_RESOURCE;
				desc.format = Format::R16G16_FLOAT;
				desc.width = internalResolution.x / 8;
				desc.height = internalResolution.y / 8;
				assert(ComputeTextureMemorySizeInBytes(desc) <= ComputeTextureMemorySizeInBytes(rtParticleDistortion.desc)); // aliasing check
				device->CreateTexture(&desc, nullptr, &rtWaterRipple, &rtParticleDistortion); // aliased!
				device->SetName(&rtWaterRipple, "renderpath3D.rtWaterRipple");
			}
		}
		else
		{
			rtWaterRipple = {};
		}

		// The ocean surface draws into a depth buffer of its own rather than
		// the main one. It has to write depth: its displaced grid folds over
		// itself wherever choppy waves overhang, so two water fragments land on
		// one pixel and only depth picks the nearer. But the transparent pass
		// is now ordered by water side, so nothing drawn after the ocean should
		// be depth-tested against it - that is what used to erase the editor
		// overlays and the gaussian splats across the waterline. Handing it a
		// copy of the opaque depth satisfies both: the water is still occluded
		// by geometry in front of it and by itself, and its writes are thrown
		// away with the copy.
		//
		// Allocated only while there is an ocean, because it is a full
		// resolution depth buffer and matches the main one sample for sample.
		if (scene->weather.IsOceanEnabled() && depthBuffer_Main.IsValid())
		{
			if (depthBuffer_Ocean.desc.width != depthBuffer_Main.desc.width ||
				depthBuffer_Ocean.desc.height != depthBuffer_Main.desc.height ||
				depthBuffer_Ocean.desc.sample_count != depthBuffer_Main.desc.sample_count)
			{
				TextureDesc desc = depthBuffer_Main.desc;
				device->CreateTexture(&desc, nullptr, &depthBuffer_Ocean);
				device->SetName(&depthBuffer_Ocean, "renderpath3D.depthBuffer_Ocean");
			}
		}
		else
		{
			depthBuffer_Ocean = {};
		}

		if (wi::renderer::GetSurfelGIEnabled())
		{
			if (!surfelGIResources.result.IsValid())
			{
				wi::renderer::CreateSurfelGIResources(surfelGIResources, internalResolution);
			}
		}

		if (wi::renderer::GetVXGIEnabled())
		{
			if (!vxgiResources.IsValid())
			{
				wi::renderer::CreateVXGIResources(vxgiResources, internalResolution);
			}
		}
		else
		{
			vxgiResources = {};
		}

		// Check whether reprojected depth is required:
		if (!first_frame && wi::renderer::IsMeshShaderAllowed() && wi::renderer::IsMeshletOcclusionCullingEnabled())
		{
			TextureDesc desc;
			desc.bind_flags = BindFlag::SHADER_RESOURCE | BindFlag::UNORDERED_ACCESS;
			desc.format = Format::R16_UNORM;
			desc.width = internalResolution.x;
			desc.height = internalResolution.y;
			desc.mip_levels = GetMipCount(desc.width, desc.height, 1, 4);
			desc.layout = ResourceState::SHADER_RESOURCE_COMPUTE;
			device->CreateTexture(&desc, nullptr, &reprojectedDepth);
			device->SetName(&reprojectedDepth, "renderpath3D.reprojectedDepth");

			for (uint32_t i = 0; i < reprojectedDepth.desc.mip_levels; ++i)
			{
				int subresource_index;
				subresource_index = device->CreateSubresource(&reprojectedDepth, SubresourceType::SRV, 0, 1, i, 1);
				assert(subresource_index == i);
				subresource_index = device->CreateSubresource(&reprojectedDepth, SubresourceType::UAV, 0, 1, i, 1);
				assert(subresource_index == i);
			}
		}
		else
		{
			reprojectedDepth = {};
		}

		// Check whether velocity buffer is required:
		if (
			getMotionBlurEnabled() ||
			wi::renderer::GetTemporalAAEnabled() ||
			getSSREnabled() ||
			getSSGIEnabled() ||
			getRaytracedReflectionEnabled() ||
			getRaytracedDiffuseEnabled() ||
			wi::renderer::GetRaytracedShadowsEnabled() ||
			getAO() == AO::AO_RTAO ||
			wi::renderer::GetVariableRateShadingClassification() ||
			getFSR2Enabled() ||
			reprojectedDepth.IsValid()
			)
		{
			if (!rtVelocity.IsValid())
			{
				TextureDesc desc;
				desc.format = Format::R16G16_FLOAT;
				desc.bind_flags = BindFlag::SHADER_RESOURCE | BindFlag::UNORDERED_ACCESS | BindFlag::RENDER_TARGET;
				desc.width = internalResolution.x;
				desc.height = internalResolution.y;
				desc.layout = ResourceState::SHADER_RESOURCE_COMPUTE;
				device->CreateTexture(&desc, nullptr, &rtVelocity);
				device->SetName(&rtVelocity, "renderpath3D.rtVelocity");
			}
		}
		else
		{
			rtVelocity = {};
		}

		// Check whether shadow mask is required:
		if (wi::renderer::GetScreenSpaceShadowsEnabled() || wi::renderer::GetRaytracedShadowsEnabled())
		{
			if (!rtShadow.IsValid())
			{
				TextureDesc desc;
				desc.bind_flags = BindFlag::SHADER_RESOURCE | BindFlag::UNORDERED_ACCESS;
				desc.format = Format::R8_UNORM;
				desc.array_size = 16;
				desc.width = internalResolution.x;
				desc.height = internalResolution.y;
				desc.layout = ResourceState::SHADER_RESOURCE_COMPUTE;
				device->CreateTexture(&desc, nullptr, &rtShadow);
				device->SetName(&rtShadow, "renderpath3D.rtShadow");
			}
		}
		else
		{
			rtShadow = {};
		}

		if (getFSR2Enabled())
		{
			// FSR2 also acts as a temporal AA, so we inform the shaders about it here
			//	This will allow improved stochastic alpha test transparency
			frameCB.options |= OPTION_BIT_TEMPORALAA_ENABLED;
			uint x = frameCB.frame_count % 4;
			uint y = frameCB.frame_count / 4;
			frameCB.temporalaa_samplerotation = (x & 0x000000FF) | ((y & 0x000000FF) << 8);
		}

		// Check whether visibility resources are required:
		if (
			visibility_shading_in_compute ||
			getSSREnabled() ||
			getSSGIEnabled() ||
			getRaytracedReflectionEnabled() ||
			getRaytracedDiffuseEnabled() ||
			wi::renderer::GetScreenSpaceShadowsEnabled() ||
			wi::renderer::GetRaytracedShadowsEnabled() ||
			wi::renderer::GetVXGIEnabled()
			)
		{
			if (!visibilityResources.IsValid())
			{
				wi::renderer::CreateVisibilityResources(visibilityResources, internalResolution);
			}
		}
		else
		{
			visibilityResources.DeleteOptionalResources();
			if (!visibilityResources.IsValidSimple())
			{
				wi::renderer::CreateVisibilityResourcesSimple(visibilityResources, internalResolution);
			}
		}

		// Check for depth of field allocation:
		//	The underwater vision blur reuses the same DoF pipeline, so allocate
		//	the resources when the ocean underwater blur is enabled too.
		const bool underwaterBlurNeedsDof =
			scene->weather.IsOceanEnabled() && getUnderwaterBlurEnabled();
		if ((getDepthOfFieldEnabled() &&
			getDepthOfFieldStrength() > 0 &&
			camera->aperture_size > 0) ||
			underwaterBlurNeedsDof
			)
		{
			if (!depthoffieldResources.IsValid())
			{
				XMUINT2 resolution = GetInternalResolution();
				if (getFSR2Enabled())
				{
					resolution = XMUINT2(GetPhysicalWidth(), GetPhysicalHeight());
				}
				wi::renderer::CreateDepthOfFieldResources(depthoffieldResources, resolution);
			}
		}
		else
		{
			depthoffieldResources = {};
		}

		// Check for motion blur allocation:
		if (getMotionBlurEnabled() && getMotionBlurStrength() > 0)
		{
			if (!motionblurResources.IsValid())
			{
				XMUINT2 resolution = GetInternalResolution();
				if (getFSR2Enabled())
				{
					resolution = XMUINT2(GetPhysicalWidth(), GetPhysicalHeight());
				}
				wi::renderer::CreateMotionBlurResources(motionblurResources, resolution);
			}
		}
		else
		{
			motionblurResources = {};
		}

		// Keep a copy of last frame's depth buffer for temporal disocclusion checks, so swap with current one every frame:
		std::swap(depthBuffer_Copy, depthBuffer_Copy1);

		visibilityResources.depthbuffer = &depthBuffer_Copy;
		if (getMSAASampleCount() > 1)
		{
			visibilityResources.primitiveID_resolved = &rtPrimitiveID;
		}
		else
		{
			visibilityResources.primitiveID_resolved = nullptr;
		}

		camera->canvas.init(*this);
		camera->width = (float)internalResolution.x;
		camera->height = (float)internalResolution.y;
		camera->scissor = GetScissorInternalResolution();
		camera->sample_count = depthBuffer_Main.desc.sample_count;
		camera->shadercamera_options = SHADERCAMERA_OPTION_NONE;

		// Let the draw shaders fog themselves with the water. Nothing more than
		// "does this scene have an ocean": GetWaterFog rejects a segment that
		// cannot touch the water with two scalar compares against the still
		// plane, which is cheaper than anything that could be decided here and
		// is exact, whereas a test on the eye alone has to be loose enough to
		// clear the wave displacement.
		//
		// Set on camera_reflection below only while the real eye is submerged -
		// see the reasoning there.
		if (scene->weather.IsOceanEnabled())
		{
			camera->shadercamera_options |= SHADERCAMERA_OPTION_WATER_FOG;
			if (getUnderwaterGodRaysProceduralEnabled())
			{
				camera->shadercamera_options |= SHADERCAMERA_OPTION_UNDERWATER_GODRAYS;
			}
			if (getWaterSunShaftsEnabled())
			{
				camera->shadercamera_options |= SHADERCAMERA_OPTION_WATER_SUN_SHAFTS;
			}
		}

		camera->texture_primitiveID_index = device->GetDescriptorIndex(&rtPrimitiveID, SubresourceType::SRV);
		camera->texture_depth_index = device->GetDescriptorIndex(&depthBuffer_Copy, SubresourceType::SRV);
		camera->texture_velocity_index = device->GetDescriptorIndex(&rtVelocity, SubresourceType::SRV);
		camera->texture_normal_roughness_index = device->GetDescriptorIndex(&visibilityResources.texture_normal_roughness, SubresourceType::SRV);
		camera->buffer_entitytiles_index = device->GetDescriptorIndex(&tiledLightResources.entityTiles, SubresourceType::SRV);
		camera->texture_reflection_index = device->GetDescriptorIndex(&rtReflection, SubresourceType::SRV);
		camera->texture_reflection_depth_index = device->GetDescriptorIndex(&depthBuffer_Reflection, SubresourceType::SRV);
		camera->texture_refraction_index = device->GetDescriptorIndex(&rtSceneCopy, SubresourceType::SRV);
		camera->texture_waterriples_index = device->GetDescriptorIndex(&rtWaterRipple, SubresourceType::SRV);
		camera->texture_ao_index = device->GetDescriptorIndex(&rtAO, SubresourceType::SRV);
		camera->texture_ssr_index = device->GetDescriptorIndex(&rtSSR, SubresourceType::SRV);
		camera->texture_ssgi_index = device->GetDescriptorIndex(&rtSSGI, SubresourceType::SRV);
		if (rtShadow.IsValid())
		{
			camera->shadercamera_options |= SHADERCAMERA_OPTION_USE_SHADOW_MASK;
			camera->texture_rtshadow_index = device->GetDescriptorIndex(&rtShadow, SubresourceType::SRV);
		}
		else
		{
			camera->texture_rtshadow_index = device->GetDescriptorIndex(wi::texturehelper::getWhite(), SubresourceType::SRV); // AMD descriptor branching fix
		}
		camera->texture_rtdiffuse_index = device->GetDescriptorIndex(&rtRaytracedDiffuse, SubresourceType::SRV);
		camera->texture_surfelgi_index = device->GetDescriptorIndex(&surfelGIResources.result, SubresourceType::SRV);
		camera->texture_vxgi_diffuse_index = device->GetDescriptorIndex(&vxgiResources.diffuse, SubresourceType::SRV);
		if (wi::renderer::GetVXGIReflectionsEnabled())
		{
			camera->texture_vxgi_specular_index = device->GetDescriptorIndex(&vxgiResources.specular, SubresourceType::SRV);
		}
		else
		{
			camera->texture_vxgi_specular_index = -1;
		}
		camera->texture_reprojected_depth_index = device->GetDescriptorIndex(&reprojectedDepth, SubresourceType::SRV);

		camera_reflection.canvas.init(*this);
		camera_reflection.width = (float)depthBuffer_Reflection_render.desc.width;
		camera_reflection.height = (float)depthBuffer_Reflection_render.desc.height;
		camera_reflection.scissor.left = 0;
		camera_reflection.scissor.top = 0;
		camera_reflection.scissor.right = (int)depthBuffer_Reflection_render.desc.width;
		camera_reflection.scissor.bottom = (int)depthBuffer_Reflection_render.desc.height;
		camera_reflection.sample_count = depthBuffer_Reflection_render.desc.sample_count;
		camera_reflection.shadercamera_options = SHADERCAMERA_OPTION_NONE;

		// The reflection's own fragments carry the water they are seen through,
		// the way every other draw does.
		//
		// Only while the real eye is submerged, because a mirrored segment is
		// not half wet and half dry the way it looks. Reflect across the water
		// plane and every part of it maps to a real leg on the SAME side as the
		// real eye: with the eye below, the mirrored ray's above-plane half is
		// the real surface-to-eye leg and its below-plane half the real
		// object-to-surface leg, and both are genuinely under water. With the
		// eye above, the whole mirrored ray maps to legs that are in air, and
		// fogging it would put water over a reflection of the sky.
		//
		// GetWaterFog measures against the plane, so submerged it contributes
		// exactly the object-to-surface leg and nothing more - the
		// surface-to-eye column is applied once to the whole fragment by
		// ApplyWaterFogAtSurface in oceanSurfacePS, and must not be counted
		// here as well.
		//
		// God rays deliberately left off: they sweep a screen space pattern
		// around where the sun projects, and a mirrored view has no such point.
		if (scene->weather.IsOceanEnabled() &&
			camera->Eye.y <= scene->weather.oceanParameters.waterHeight)
		{
			camera_reflection.shadercamera_options |=
				SHADERCAMERA_OPTION_WATER_FOG;
		}
		camera_reflection.texture_primitiveID_index = -1;
		camera_reflection.texture_depth_index = device->GetDescriptorIndex(&depthBuffer_Reflection, SubresourceType::SRV);
		camera_reflection.texture_velocity_index = -1;
		camera_reflection.texture_normal_roughness_index = -1;
		camera_reflection.buffer_entitytiles_index = device->GetDescriptorIndex(&tiledLightResources_planarReflection.entityTiles, SubresourceType::SRV);
		camera_reflection.texture_reflection_index = -1;
		camera_reflection.texture_reflection_depth_index = -1;
		// The opaque half of the reflection, downsampled between the two halves
		// of the reflection colour pass. Bound here rather than there because
		// the camera constant buffer carries a descriptor index, not content:
		// the pass records the copy before any draw that reads it.
		camera_reflection.texture_refraction_index = device->GetDescriptorIndex(&rtSceneCopy_Reflection, SubresourceType::SRV);
		camera_reflection.texture_waterriples_index = -1;
		camera_reflection.texture_ao_index = -1;
		camera_reflection.texture_ssr_index = -1;
		camera_reflection.texture_ssgi_index = -1;
		camera_reflection.texture_rtshadow_index = device->GetDescriptorIndex(wi::texturehelper::getWhite(), SubresourceType::SRV); // AMD descriptor branching fix
		camera_reflection.texture_rtdiffuse_index = -1;
		camera_reflection.texture_surfelgi_index = -1;
		camera_reflection.texture_vxgi_diffuse_index = -1;
		camera_reflection.texture_vxgi_specular_index = -1;
		camera_reflection.texture_reprojected_depth_index = -1;
		// The volume is built along the main camera's rays and means nothing
		// along the reflection's.
		camera_reflection.texture_volumetricfroxels_index = -1;

		video_cmd = {};
		if (getSceneUpdateEnabled() && scene->videos.GetCount() > 0)
		{
			for (size_t i = 0; i < scene->videos.GetCount(); ++i)
			{
				const wi::scene::VideoComponent& video = scene->videos[i];
				if (wi::video::IsDecodingRequired(&video.videoinstance))
				{
					video_cmd = device->BeginCommandList(QUEUE_VIDEO_DECODE);
					break;
				}
			}
			for (size_t i = 0; i < scene->videos.GetCount(); ++i)
			{
				wi::scene::VideoComponent& video = scene->videos[i];
				wi::video::DecodeVideo(&video.videoinstance, video_cmd);
			}
		}

		if (getMeshBlendEnabled() && visibility_main.IsMeshBlendVisible())
		{
			if (!meshblendResources.IsValid())
			{
				wi::renderer::CreateMeshBlendResources(meshblendResources, internalResolution);
			}
		}
		else
		{
			meshblendResources = {};
		}

		prerender_happened = true;

		RenderPath2D::PreRender();
	}

	void RenderPath3D::Render() const
	{
		if (!prerender_happened)
		{
			// Since 0.71.694: PreRender must be called before Render() because it sets up rendering resources!
			//	The proper fix is to call PreRender() yourself for a manually managed RenderPath3D
			//	But if you don't do that, as a last resort it will be called here using const_cast
			assert(0);
			const_cast<RenderPath3D*>(this)->PreRender();
		}
		prerender_happened = false;

		GraphicsDevice* device = wi::graphics::GetDevice();
		wi::jobsystem::context ctx;

		CommandList cmd_copypages;
		if (scene->terrains.GetCount() > 0)
		{
			cmd_copypages = device->BeginCommandList(QUEUE_COPY);
			wi::jobsystem::Execute(ctx, [this, cmd_copypages](wi::jobsystem::JobArgs args) {
				for (size_t i = 0; i < scene->terrains.GetCount(); ++i)
				{
					scene->terrains[i].CopyVirtualTexturePageStatusGPU(cmd_copypages);
				}
			});
		}

		// Preparing the frame:
		CommandList cmd = device->BeginCommandList();
		wi::renderer::ProcessDeferredTextureRequests(cmd); // Execute it first thing in the frame here, on main thread, to not allow other thread steal it and execute on different command list!
		CommandList cmd_prepareframe = cmd;
		wi::jobsystem::Execute(ctx, [this, cmd](wi::jobsystem::JobArgs args) {
			GraphicsDevice* device = wi::graphics::GetDevice();

			wi::renderer::BindCameraCB(
				*camera,
				camera_previous,
				camera_reflection,
				cmd
			);
			wi::renderer::UpdateRenderData(visibility_main, frameCB, cmd);

			GPUBarrier barriers[] = {
				GPUBarrier::Image(&debugUAV, debugUAV.desc.layout, ResourceState::UNORDERED_ACCESS),
				GPUBarrier::Aliasing(&rtPostprocess, &rtPrimitiveID),
			};
			device->Barrier(barriers, arraysize(barriers), cmd);

		});

		// async compute parallel with depth prepass
		cmd = device->BeginCommandList(QUEUE_COMPUTE);
		CommandList cmd_prepareframe_async = cmd;
		device->WaitCommandList(cmd, cmd_prepareframe);
		if (cmd_copypages.IsValid())
		{
			device->WaitCommandList(cmd, cmd_copypages);
		}
		wi::jobsystem::Execute(ctx, [this, cmd](wi::jobsystem::JobArgs args) {

			wi::renderer::BindCameraCB(
				*camera,
				camera_previous,
				camera_reflection,
				cmd
			);
			wi::renderer::UpdateRenderDataAsync(visibility_main, frameCB, cmd);

			if (scene->IsWetmapProcessingRequired())
			{
				wi::renderer::RefreshWetmaps(visibility_main, cmd);
			}

			if (scene->IsAccelerationStructureUpdateRequested())
			{
				wi::renderer::UpdateRaytracingAccelerationStructures(*scene, cmd);
			}

			if (wi::renderer::GetSurfelGIEnabled())
			{
				wi::renderer::SurfelGI(
					surfelGIResources,
					*scene,
					cmd
				);
			}

			if (wi::renderer::GetDDGIEnabled() && getSceneUpdateEnabled())
			{
				wi::renderer::DDGI(
					*scene,
					cmd
				);
			}

		});

		static const uint32_t drawscene_flags =
			wi::renderer::DRAWSCENE_OPAQUE |
			wi::renderer::DRAWSCENE_IMPOSTOR |
			wi::renderer::DRAWSCENE_HAIRPARTICLE |
			wi::renderer::DRAWSCENE_TESSELLATION |
			wi::renderer::DRAWSCENE_OCCLUSIONCULLING |
			wi::renderer::DRAWSCENE_MAINCAMERA
			;

		// Main camera depth prepass:
		cmd = device->BeginCommandList();
		CommandList cmd_maincamera_prepass = cmd;
		wi::jobsystem::Execute(ctx, [this, cmd](wi::jobsystem::JobArgs args) {

			GraphicsDevice* device = wi::graphics::GetDevice();

			wi::renderer::BindCameraCB(
				*camera,
				camera_previous,
				camera_reflection,
				cmd
			);

			wi::renderer::RefreshImpostors(*scene, cmd);

			if (reprojectedDepth.IsValid())
			{
				wi::renderer::ComputeReprojectedDepthPyramid(
					depthBuffer_Copy,
					rtVelocity,
					reprojectedDepth,
					cmd
				);
			}

			RenderPassImage rp[] = {
				RenderPassImage::DepthStencil(
					&depthBuffer_Main,
					RenderPassImage::LoadOp::CLEAR,
					RenderPassImage::StoreOp::STORE,
					ResourceState::DEPTHSTENCIL,
					ResourceState::DEPTHSTENCIL,
					ResourceState::DEPTHSTENCIL
				),
				RenderPassImage::RenderTarget(
					&rtPrimitiveID_render,
					RenderPassImage::LoadOp::CLEAR,
					RenderPassImage::StoreOp::STORE,
					ResourceState::SHADER_RESOURCE_COMPUTE,
					ResourceState::SHADER_RESOURCE_COMPUTE
				),
			};
			device->RenderPassBegin(rp, arraysize(rp), cmd);

			device->EventBegin("Opaque Z-prepass", cmd);
			auto range = wi::profiler::BeginRangeGPU("Z-Prepass", cmd);

			wi::graphics::Rect scissor = GetScissorInternalResolution();
			device->BindScissorRects(1, &scissor, cmd);

			Viewport vp;
			vp.width = (float)depthBuffer_Main.GetDesc().width;
			vp.height = (float)depthBuffer_Main.GetDesc().height;

			// Foreground:
			vp.min_depth = 1 - foreground_depth_range;
			vp.max_depth = 1;
			device->BindViewports(1, &vp, cmd);
			wi::renderer::DrawScene(
				visibility_main,
				RENDERPASS_PREPASS,
				cmd,
				wi::renderer::DRAWSCENE_OPAQUE |
				wi::renderer::DRAWSCENE_FOREGROUND_ONLY |
				wi::renderer::DRAWSCENE_MAINCAMERA
			);

			// Regular:
			vp.min_depth = 0;
			vp.max_depth = 1;
			device->BindViewports(1, &vp, cmd);
			wi::renderer::DrawScene(
				visibility_main,
				RENDERPASS_PREPASS,
				cmd,
				drawscene_flags
			);

			wi::profiler::EndRange(range);
			device->EventEnd(cmd);

			device->RenderPassEnd(cmd);

		});

		// Main camera compute effects:
		//	(async compute, parallel to "shadow maps" and "update textures",
		//	must finish before "main scene opaque color pass")
		cmd = device->BeginCommandList(QUEUE_COMPUTE);
		device->WaitCommandList(cmd, cmd_maincamera_prepass);
		if (video_cmd.IsValid())
		{
			device->WaitCommandList(cmd, video_cmd);
		}
		CommandList cmd_maincamera_compute_effects = cmd;
		wi::jobsystem::Execute(ctx, [this, cmd](wi::jobsystem::JobArgs args) {

			GraphicsDevice* device = wi::graphics::GetDevice();

			for (size_t i = 0; i < scene->videos.GetCount(); ++i)
			{
				wi::scene::VideoComponent& video = scene->videos[i];
				wi::video::ResolveVideoToRGB(&video.videoinstance, cmd);
			}

			wi::renderer::BindCameraCB(
				*camera,
				camera_previous,
				camera_reflection,
				cmd
			);

			wi::renderer::Visibility_Prepare(
				visibilityResources,
				rtPrimitiveID_render,
				cmd
			);

			wi::renderer::ComputeTiledLightCulling(
				tiledLightResources,
				visibility_main,
				debugUAV,
				cmd
			);

			if (
				getSSREnabled() ||
				getSSGIEnabled() ||
				getRaytracedReflectionEnabled() ||
				getRaytracedDiffuseEnabled() ||
				wi::renderer::GetScreenSpaceShadowsEnabled() ||
				wi::renderer::GetRaytracedShadowsEnabled() ||
				wi::renderer::GetVXGIEnabled()
				)
			{
				// These post effects require surface normals and/or roughness
				wi::renderer::Visibility_Surface(
					visibilityResources,
					cmd
				);
			}

			if (rtVelocity.IsValid())
			{
				wi::renderer::Visibility_Velocity(
					visibilityResources,
					rtVelocity,
					cmd
				);
			}

			if (wi::renderer::GetSurfelGIEnabled())
			{
				wi::renderer::SurfelGI_Coverage(
					surfelGIResources,
					*scene,
					depthBuffer_Copy,
					debugUAV,
					cmd
				);
			}

			RenderAO(cmd);

			if (wi::renderer::GetVariableRateShadingClassification() && device->CheckCapability(GraphicsDeviceCapability::VARIABLE_RATE_SHADING_TIER2))
			{
				wi::renderer::ComputeShadingRateClassification(
					rtShadingRate,
					debugUAV,
					cmd
				);
			}

			RenderSSR(cmd);

			RenderSSGI(cmd);

			if (wi::renderer::GetWireframeMode() == wi::renderer::WIREFRAME_DISABLED)
			{
				if (wi::renderer::GetRaytracedShadowsEnabled())
				{
					wi::renderer::Postprocess_RTShadow(
						rtshadowResources,
						*scene,
						rtShadow,
						cmd
					);
				}
				else if (wi::renderer::GetScreenSpaceShadowsEnabled()) // rtshadow or screenspace shadow, not both
				{
					wi::renderer::Postprocess_ScreenSpaceShadow(
						screenspaceshadowResources,
						rtShadow,
						cmd,
						getScreenSpaceShadowRange(),
						getScreenSpaceShadowSampleCount()
					);
				}
			}

			if (getMeshBlendEnabled() && visibility_main.IsMeshBlendVisible())
			{
				wi::renderer::PostProcess_MeshBlend_EdgeProcess(meshblendResources, cmd);
			}

		});

		// Occlusion culling:
		CommandList cmd_occlusionculling;
		if (getOcclusionCullingEnabled())
		{
			cmd = device->BeginCommandList();
			cmd_occlusionculling = cmd;
			wi::jobsystem::Execute(ctx, [this, cmd](wi::jobsystem::JobArgs args) {

				GraphicsDevice* device = wi::graphics::GetDevice();

				device->EventBegin("Occlusion Culling", cmd);
				ScopedGPUProfiling("Occlusion Culling", cmd);

				wi::renderer::BindCameraCB(
					*camera,
					camera_previous,
					camera_reflection,
					cmd
				);

				wi::renderer::OcclusionCulling_Reset(visibility_main, cmd); // must be outside renderpass!

				RenderPassImage rp[] = {
					RenderPassImage::DepthStencil(&depthBuffer_Main),
				};
				device->RenderPassBegin(rp, arraysize(rp), &scene->queryHeap, cmd);

				wi::graphics::Rect scissor = GetScissorInternalResolution();
				device->BindScissorRects(1, &scissor, cmd);

				Viewport vp;
				vp.width = (float)depthBuffer_Main.GetDesc().width;
				vp.height = (float)depthBuffer_Main.GetDesc().height;
				device->BindViewports(1, &vp, cmd);

				wi::renderer::OcclusionCulling_Render(*camera, visibility_main, cmd);

				device->RenderPassEnd(cmd);

				wi::renderer::OcclusionCulling_Resolve(visibility_main, cmd); // must be outside renderpass!

				device->EventEnd(cmd);
			});
		}

		// Planar reflections depth prepass:
		if (getReflectionsEnabled() && visibility_main.IsRequestedPlanarReflections())
		{
			cmd = device->BeginCommandList();
			wi::jobsystem::Execute(ctx, [cmd, this](wi::jobsystem::JobArgs args) {

				GraphicsDevice* device = wi::graphics::GetDevice();

				wi::renderer::BindCameraCB(
					camera_reflection,
					camera_reflection_previous,
					camera_reflection,
					cmd
				);

				device->EventBegin("Planar reflections Z-Prepass", cmd);
				auto range = wi::profiler::BeginRangeGPU("Planar Reflections Z-Prepass", cmd);

				RenderPassImage rp[] = {
					RenderPassImage::DepthStencil(
						&depthBuffer_Reflection_render,
						RenderPassImage::LoadOp::CLEAR,
						RenderPassImage::StoreOp::STORE,
						ResourceState::SHADER_RESOURCE,
						ResourceState::DEPTHSTENCIL,
						ResourceState::SHADER_RESOURCE
					)
				};
				device->RenderPassBegin(rp, arraysize(rp), cmd);

				Viewport vp;
				vp.width = (float)depthBuffer_Reflection.GetDesc().width;
				vp.height = (float)depthBuffer_Reflection.GetDesc().height;
				vp.min_depth = 0;
				vp.max_depth = 1;
				device->BindViewports(1, &vp, cmd);

				wi::renderer::DrawScene(
					visibility_reflection,
					RENDERPASS_PREPASS_DEPTHONLY,
					cmd,
					wi::renderer::DRAWSCENE_OPAQUE |
					wi::renderer::DRAWSCENE_IMPOSTOR |
					wi::renderer::DRAWSCENE_HAIRPARTICLE |
					wi::renderer::DRAWSCENE_SKIP_PLANAR_REFLECTION_OBJECTS
				);

				device->RenderPassEnd(cmd);

				if (depthBuffer_Reflection_render.desc.sample_count > 1)
				{
					wi::renderer::ResolveMSAADepthBuffer(depthBuffer_Reflection, depthBuffer_Reflection_render, cmd);
				}

				if (scene->weather.IsRealisticSky() && scene->weather.IsRealisticSkyAerialPerspective())
				{
					wi::renderer::Postprocess_AerialPerspective(
						aerialperspectiveResources_reflection,
						cmd
					);
				}

				wi::profiler::EndRange(range); // Planar Reflections
				device->EventEnd(cmd);

			});
		}

		CommandList cmd_ocean;
		if (scene->weather.IsOceanEnabled() && scene->ocean.IsValid())
		{
			// Ocean simulation can be updated async to opaque passes:
			cmd_ocean = device->BeginCommandList(QUEUE_COMPUTE);
			if (cmd_occlusionculling.IsValid())
			{
				// Ocean occlusion culling must be waited
				device->WaitCommandList(cmd_ocean, cmd_occlusionculling);
			}
			wi::renderer::UpdateOcean(visibility_main, cmd_ocean);

			// Copying to readback is done on copy queue to use DMA instead of compute warps:
			CommandList cmd_oceancopy = device->BeginCommandList(QUEUE_COPY);
			device->WaitCommandList(cmd_oceancopy, cmd_ocean);
			wi::renderer::ReadbackOcean(visibility_main, cmd_oceancopy);
		}

		// Shadow maps:
		CommandList cmd_shadowmap;
		if (getShadowsEnabled())
		{
			cmd_shadowmap = device->BeginCommandList();
			cmd = cmd_shadowmap;
			device->WaitCommandList(cmd, cmd_prepareframe_async); // shadow map waits for UpdateRenderDataAsync (particle-shadowmap interaction)
			wi::jobsystem::Execute(ctx, [this, cmd](wi::jobsystem::JobArgs args) {
				wi::renderer::DrawShadowmaps(visibility_main, cmd);
			});
		}

		if (wi::renderer::GetVXGIEnabled() && getSceneUpdateEnabled())
		{
			cmd = device->BeginCommandList();
			wi::jobsystem::Execute(ctx, [cmd, this](wi::jobsystem::JobArgs args) {
				wi::renderer::VXGI_Voxelize(visibility_main, cmd);
			});
		}

		// Updating textures:
		CommandList cmd_updatetextures;
		if (getSceneUpdateEnabled())
		{
			cmd_updatetextures = device->BeginCommandList();
			cmd = cmd_updatetextures;
			device->WaitCommandList(cmd, cmd_prepareframe_async);
			if (cmd_ocean.IsValid())
			{
				device->WaitCommandList(cmd, cmd_ocean);
			}
			wi::jobsystem::Execute(ctx, [cmd, this](wi::jobsystem::JobArgs args) {
				wi::renderer::BindCommonResources(cmd);
				wi::renderer::BindCameraCB(
					*camera,
					camera_previous,
					camera_reflection,
					cmd
				);
				wi::renderer::RefreshLightmaps(*scene, cmd);
				wi::renderer::RefreshEnvProbes(visibility_main, cmd);
				wi::renderer::PaintDecals(*scene, cmd);
			});
		}

		// Planar reflections color pass:
		if (getReflectionsEnabled() && visibility_main.IsRequestedPlanarReflections())
		{
			cmd = device->BeginCommandList();
			wi::jobsystem::Execute(ctx, [cmd, this](wi::jobsystem::JobArgs args) {

				GraphicsDevice* device = wi::graphics::GetDevice();

				wi::renderer::BindCameraCB(
					camera_reflection,
					camera_reflection_previous,
					camera_reflection,
					cmd
				);

				wi::renderer::ComputeTiledLightCulling(
					tiledLightResources_planarReflection,
					visibility_reflection,
					Texture(),
					cmd
				);

				// The splat sort is per camera and has to happen outside the
				// render pass. It writes the same buffers the main camera's
				// sort does, and reads them back in the draw below - which is
				// safe only because this command list is begun, and so
				// submitted, before the transparent pass that sorts for the
				// main camera. The env probe path already relies on the same
				// ordering.
				wi::renderer::UpdateGaussianSplatsForCamera(
					*scene, camera_reflection, cmd
				);

				device->EventBegin("Planar reflections", cmd);
				auto range = wi::profiler::BeginRangeGPU("Planar Reflections", cmd);

				// The reflection is drawn in two halves, for the same reason
				// the main view is: what a transparent refracts has to already
				// exist as a texture before the transparent is drawn. The main
				// view's own scene copy cannot stand in here, being a different
				// projection of the scene entirely, so the reflection takes its
				// own between the halves.
				const bool reflectionMSAA =
					depthBuffer_Reflection_render.desc.sample_count > 1;
				auto beginReflectionPass = [&](
					RenderPassImage::LoadOp load,
					RenderPassImage::StoreOp store
				) {
					RenderPassImage rp[3];
					uint32_t rp_count = 0;
					if (reflectionMSAA)
					{
						rp[rp_count++] = RenderPassImage::RenderTarget(
							&rtReflection_render,
							load,
							store,
							ResourceState::RENDERTARGET,
							ResourceState::RENDERTARGET
						);
						rp[rp_count++] = RenderPassImage::Resolve(&rtReflection);
					}
					else
					{
						rp[rp_count++] = RenderPassImage::RenderTarget(
							&rtReflection_render,
							load,
							store,
							ResourceState::SHADER_RESOURCE,
							ResourceState::SHADER_RESOURCE
						);
					}
					rp[rp_count++] = RenderPassImage::DepthStencil(
						&depthBuffer_Reflection_render,
						RenderPassImage::LoadOp::LOAD,
						RenderPassImage::StoreOp::STORE,
						ResourceState::SHADER_RESOURCE,
						ResourceState::DEPTHSTENCIL,
						ResourceState::SHADER_RESOURCE
					);
					device->RenderPassBegin(rp, rp_count, cmd);

					Viewport vp;
					vp.width = (float)depthBuffer_Reflection_render.GetDesc().width;
					vp.height = (float)depthBuffer_Reflection_render.GetDesc().height;
					vp.min_depth = 0;
					vp.max_depth = 1;
					device->BindViewports(1, &vp, cmd);
				};

				// Everything a transparent in the reflection can see through
				// itself. The MSAA samples have to be kept rather than dropped
				// at the resolve, since the second half loads them back.
				beginReflectionPass(
					RenderPassImage::LoadOp::CLEAR,
					RenderPassImage::StoreOp::STORE
				);

				wi::renderer::DrawScene(
					visibility_reflection,
					RENDERPASS_MAIN,
					cmd,
					wi::renderer::DRAWSCENE_OPAQUE |
					wi::renderer::DRAWSCENE_IMPOSTOR |
					wi::renderer::DRAWSCENE_HAIRPARTICLE |
					wi::renderer::DRAWSCENE_SKIP_PLANAR_REFLECTION_OBJECTS
				);
				wi::renderer::DrawSky(*scene, cmd);

				if (scene->weather.IsRealisticSky() && scene->weather.IsRealisticSkyAerialPerspective())
				{
					// Blend Aerial Perspective on top. It is the atmosphere
					// between the eye and the OPAQUE surface at each pixel,
					// resolved from the depth prepass and zero where the sky
					// shows, so it belongs to this half and only this half:
					// after it the copy below carries the haze, and nothing
					// nearer than the opaque surface is painted with the haze
					// that stands in front of it.
					device->EventBegin("Aerial Perspective Reflection Blend", cmd);
					wi::image::Params fx;
					fx.enableFullScreen();
					fx.blendFlag = BLENDMODE_PREMULTIPLIED;
					wi::image::Draw(&aerialperspectiveResources_reflection.texture_output, fx, cmd);
					device->EventEnd(cmd);
				}

				device->RenderPassEnd(cmd);

				// A linear mip chain rather than the main view's gaussian one:
				// that needs a scratch texture of its own, and at a quarter of
				// a reflection the difference is below what the roughness mip
				// lookup can show.
				device->EventBegin("Reflection Scene MIP Chain", cmd);
				wi::renderer::Postprocess_Downsample4x(
					rtReflection, rtSceneCopy_Reflection, cmd
				);
				wi::renderer::GenerateMipChain(
					rtSceneCopy_Reflection, wi::renderer::MIPGENFILTER_LINEAR, cmd
				);
				device->EventEnd(cmd);

				beginReflectionPass(
					RenderPassImage::LoadOp::LOAD,
					reflectionMSAA
						? RenderPassImage::StoreOp::DONTCARE
						: RenderPassImage::StoreOp::STORE
				);

				wi::renderer::DrawScene(
					visibility_reflection,
					RENDERPASS_MAIN,
					cmd,
					wi::renderer::DRAWSCENE_TRANSPARENT |
					wi::renderer::DRAWSCENE_SKIP_PLANAR_REFLECTION_OBJECTS
				); // separate renderscene, to be drawn after opaque and transparent sort order

				wi::renderer::DrawSoftParticles(visibility_reflection, false, cmd);
				wi::renderer::DrawGaussianSplats(*scene, camera_reflection, cmd);
				wi::renderer::DrawSpritesAndFonts(*scene, camera_reflection, false, cmd);
				// Keep the queue: the main camera draws these same trails, on a
				// job that may still be recording. RenderPath3D::Render clears
				// them once every pass is done with them.
				wi::renderer::DrawTrails(camera_reflection, cmd, false);

				device->RenderPassEnd(cmd);

				wi::profiler::EndRange(range); // Planar Reflections
				device->EventEnd(cmd);
			});
		}

		// Main camera opaque color pass:
		cmd = device->BeginCommandList();
		device->WaitCommandList(cmd, cmd_maincamera_compute_effects);
		if (cmd_ocean.IsValid())
		{
			// The wave displacement settles which side of the water each column
			// of the froxel volume is on, and the opaque shaders read it again
			// through ApplyWaterFog.
			device->WaitCommandList(cmd, cmd_ocean);
		}
		wi::jobsystem::Execute(ctx, [this, cmd](wi::jobsystem::JobArgs args) {

			GraphicsDevice* device = wi::graphics::GetDevice();
			device->EventBegin("Opaque Scene", cmd);

			wi::renderer::BindCameraCB(
				*camera,
				camera_previous,
				camera_reflection,
				cmd
			);

			// Built before anything shades, because everything that shades
			// reads it - opaque geometry included.
			volumetricFroxels.Build(
				*camera, camera_previous, camera_reflection, cmd);

			if (getRaytracedReflectionEnabled())
			{
				wi::renderer::Postprocess_RTReflection(
					rtreflectionResources,
					*scene,
					rtSSR,
					cmd,
					getRaytracedReflectionsRange(),
					getReflectionRoughnessCutoff()
				);
			}
			if (getRaytracedDiffuseEnabled())
			{
				wi::renderer::Postprocess_RTDiffuse(
					rtdiffuseResources,
					*scene,
					rtRaytracedDiffuse,
					cmd,
					getRaytracedDiffuseRange()
				);
			}
			if (wi::renderer::GetVXGIEnabled())
			{
				wi::renderer::VXGI_Resolve(
					vxgiResources,
					*scene,
					cmd
				);
			}

			if (visibility_shading_in_compute)
			{
				wi::renderer::Visibility_Shade(
					visibilityResources,
					rtMain,
					cmd
				);
			}

			Viewport vp;
			vp.width = (float)depthBuffer_Main.GetDesc().width;
			vp.height = (float)depthBuffer_Main.GetDesc().height;
			device->BindViewports(1, &vp, cmd);

			wi::graphics::Rect scissor = GetScissorInternalResolution();
			device->BindScissorRects(1, &scissor, cmd);

			if (getOutlineEnabled())
			{
				// Cut off outline source from depth:
				device->EventBegin("Outline Source", cmd);

				if (getMSAASampleCount() > 1)
				{
					RenderPassImage rp[] = {
						RenderPassImage::RenderTarget(&rtOutlineSource_MSAA, RenderPassImage::LoadOp::CLEAR, RenderPassImage::StoreOp::DONTCARE, ResourceState::RENDERTARGET, ResourceState::RENDERTARGET),
						RenderPassImage::DepthStencil(&depthBuffer_Main, RenderPassImage::LoadOp::LOAD),
						RenderPassImage::Resolve(&rtOutlineSource)
					};
					device->RenderPassBegin(rp, arraysize(rp), cmd);
				}
				else
				{
					RenderPassImage rp[] = {
						RenderPassImage::RenderTarget(&rtOutlineSource, RenderPassImage::LoadOp::CLEAR),
						RenderPassImage::DepthStencil(&depthBuffer_Main, RenderPassImage::LoadOp::LOAD)
					};
					device->RenderPassBegin(rp, arraysize(rp), cmd);
				}
				wi::image::Params params;
				params.enableFullScreen();
				params.stencilRefMode = wi::image::STENCILREFMODE_ENGINE;
				params.stencilComp = wi::image::STENCILMODE_EQUAL;
				params.stencilRef = wi::enums::STENCILREF_OUTLINE;
				wi::image::Draw(&depthBuffer_Copy, params, cmd);
				params.stencilRef = wi::enums::STENCILREF_CUSTOMSHADER_OUTLINE;
				wi::image::Draw(&depthBuffer_Copy, params, cmd);
				device->RenderPassEnd(cmd);
				device->EventEnd(cmd);
			}

			RenderPassImage rp[4] = {};
			uint32_t rp_count = 0;
			rp[rp_count++] = RenderPassImage::RenderTarget(
				&rtMain_render,
				visibility_shading_in_compute ? RenderPassImage::LoadOp::LOAD : RenderPassImage::LoadOp::CLEAR
			);
			if (device->CheckCapability(GraphicsDeviceCapability::VARIABLE_RATE_SHADING_TIER2) && rtShadingRate.IsValid())
			{
				rp[rp_count++] = RenderPassImage::ShadingRateSource(&rtShadingRate, ResourceState::UNORDERED_ACCESS, ResourceState::UNORDERED_ACCESS);
			}
			rp[rp_count++] = RenderPassImage::DepthStencil(
				&depthBuffer_Main,
				RenderPassImage::LoadOp::LOAD,
				RenderPassImage::StoreOp::STORE,
				ResourceState::DEPTHSTENCIL,
				ResourceState::DEPTHSTENCIL,
				ResourceState::DEPTHSTENCIL
			);
			device->RenderPassBegin(rp, rp_count, cmd, RenderPassFlags::ALLOW_UAV_WRITES);

			if (visibility_shading_in_compute)
			{
				// In visibility compute shading, the impostors must still be drawn using rasterization:
				wi::renderer::DrawScene(
					visibility_main,
					RENDERPASS_MAIN,
					cmd,
					wi::renderer::DRAWSCENE_IMPOSTOR
				);
			}
			else
			{
				auto range = wi::profiler::BeginRangeGPU("Opaque Scene", cmd);

				// Foreground:
				vp.min_depth = 1 - foreground_depth_range;
				vp.max_depth = 1;
				device->BindViewports(1, &vp, cmd);
				wi::renderer::DrawScene(
					visibility_main,
					RENDERPASS_MAIN,
					cmd,
					wi::renderer::DRAWSCENE_OPAQUE |
					wi::renderer::DRAWSCENE_FOREGROUND_ONLY |
					wi::renderer::DRAWSCENE_MAINCAMERA
				);

				// Regular:
				vp.min_depth = 0;
				vp.max_depth = 1;
				device->BindViewports(1, &vp, cmd);
				wi::renderer::DrawScene(
					visibility_main,
					RENDERPASS_MAIN,
					cmd,
					drawscene_flags
				);
				wi::profiler::EndRange(range); // Opaque Scene
			}

			wi::renderer::DrawSky(*scene, cmd, false); // Note: volumetric cloud sampling disabled in sky shader, instead the postprocess will be used for a high quality effect

			RenderOutline(cmd);

			device->RenderPassEnd(cmd);

			if (getMeshBlendEnabled() && visibility_main.IsMeshBlendVisible())
			{
				rp[0].loadop = RenderPassImage::LoadOp::LOAD;
				wi::renderer::PostProcess_MeshBlend_Resolve(meshblendResources, rtMain, rp, rp_count, cmd);
			}

			if (rtAO.IsValid())
			{
				device->Barrier(GPUBarrier::Aliasing(&rtAO, &rtParticleDistortion), cmd);
			}

			device->EventEnd(cmd);
		});

		CommandList cmd_weathereffect;
		if (scene->weather.IsVolumetricClouds() || (scene->weather.IsRealisticSky() && scene->weather.IsRealisticSkyAerialPerspective()))
		{
			// Weather effects will be running async to opaque scene shading pass:
			cmd_weathereffect = device->BeginCommandList(QUEUE_COMPUTE);
			cmd = cmd_weathereffect;
			if (cmd_updatetextures.IsValid())
			{
				// delay after update textures phase to balance with opaque pass
				device->WaitCommandList(cmd_weathereffect, cmd_updatetextures);
			}
			else if (cmd_shadowmap.IsValid())
			{
				// wait for shadow maps is required for shadow receiver volumetric clouds or aerial perspective
				device->WaitCommandList(cmd_weathereffect, cmd_shadowmap);
			}
			wi::jobsystem::Execute(ctx, [this, cmd](wi::jobsystem::JobArgs args) {

				wi::renderer::BindCameraCB(
					*camera,
					camera_previous,
					camera_reflection,
					cmd
				);

				if (scene->weather.IsVolumetricClouds())
				{
					wi::renderer::Postprocess_VolumetricClouds(
						volumetriccloudResources,
						cmd,
						*camera,
						camera_previous,
						camera_reflection,
						(wi::renderer::GetTemporalAAEnabled() && wi::renderer::GetWireframeMode() == wi::renderer::WIREFRAME_DISABLED) || getFSR2Enabled(),
						scene->weather.volumetricCloudsWeatherMapFirst.IsValid() ? &scene->weather.volumetricCloudsWeatherMapFirst.GetTexture() : nullptr,
						scene->weather.volumetricCloudsWeatherMapSecond.IsValid() ? &scene->weather.volumetricCloudsWeatherMapSecond.GetTexture() : nullptr
					);
				}
				if (scene->weather.IsRealisticSky() && scene->weather.IsRealisticSkyAerialPerspective())
				{
					wi::renderer::Postprocess_AerialPerspective(
						aerialperspectiveResources,
						cmd
					);
				}
			});
		}

		// Opaque color part 2 - blend weather effects or resolve MSAA:
		if (cmd_weathereffect.IsValid() || getMSAASampleCount() > 1)
		{
			cmd = device->BeginCommandList();
			if (cmd_weathereffect.IsValid())
			{
				device->WaitCommandList(cmd, cmd_weathereffect);
			}
			wi::jobsystem::Execute(ctx, [this, cmd](wi::jobsystem::JobArgs args) {
				GraphicsDevice* device = wi::graphics::GetDevice();
				device->EventBegin("Opaque Scene (weather blend / MSAA resolve)", cmd);

				wi::renderer::BindCameraCB(
					*camera,
					camera_previous,
					camera_reflection,
					cmd
				);

				Viewport vp;
				vp.width = (float)depthBuffer_Main.GetDesc().width;
				vp.height = (float)depthBuffer_Main.GetDesc().height;
				device->BindViewports(1, &vp, cmd);

				RenderPassImage rp[4] = {};
				uint32_t rp_count = 0;
				rp[rp_count++] = RenderPassImage::RenderTarget(&rtMain_render, RenderPassImage::LoadOp::LOAD);
				if (getMSAASampleCount() > 1)
				{
					rp[rp_count++] = RenderPassImage::Resolve(&rtMain);
				}
				device->RenderPassBegin(rp, rp_count, cmd);

				// Blend Aerial Perspective on top:
				if (scene->weather.IsRealisticSky() && scene->weather.IsRealisticSkyAerialPerspective())
				{
					device->EventBegin("Aerial Perspective Blend", cmd);
					wi::image::Params fx;
					fx.enableFullScreen();
					fx.blendFlag = BLENDMODE_PREMULTIPLIED;
					wi::image::Draw(&aerialperspectiveResources.texture_output, fx, cmd);
					device->EventEnd(cmd);
				}

				// Blend the volumetric clouds on top:
				if (scene->weather.IsVolumetricClouds())
				{
					wi::renderer::Postprocess_VolumetricClouds_Upsample(volumetriccloudResources, cmd);
				}

				device->RenderPassEnd(cmd);

				device->EventEnd(cmd);
			});
		}

		if (scene->terrains.GetCount() > 0)
		{
			CommandList cmd_allocation_tilerequest = device->BeginCommandList(QUEUE_COMPUTE);
			device->WaitCommandList(cmd_allocation_tilerequest, cmd); // wait for opaque scene
			wi::jobsystem::Execute(ctx, [this, cmd_allocation_tilerequest](wi::jobsystem::JobArgs args) {
				for (size_t i = 0; i < scene->terrains.GetCount(); ++i)
				{
					scene->terrains[i].AllocateVirtualTextureTileRequestsGPU(cmd_allocation_tilerequest);
				}
			});

			CommandList cmd_writeback_tilerequest = device->BeginCommandList(QUEUE_COPY);
			device->WaitCommandList(cmd_writeback_tilerequest, cmd_allocation_tilerequest);
			wi::jobsystem::Execute(ctx, [this, cmd_writeback_tilerequest](wi::jobsystem::JobArgs args) {
				for (size_t i = 0; i < scene->terrains.GetCount(); ++i)
				{
					scene->terrains[i].WritebackTileRequestsGPU(cmd_writeback_tilerequest);
				}
			});
		}

		// Transparents, post processes, etc:
		cmd = device->BeginCommandList();
		if (cmd_ocean.IsValid())
		{
			device->WaitCommandList(cmd, cmd_ocean);
		}
		wi::jobsystem::Execute(ctx, [this, cmd](wi::jobsystem::JobArgs args) {

			GraphicsDevice* device = wi::graphics::GetDevice();

			wi::renderer::BindCameraCB(
				*camera,
				camera_previous,
				camera_reflection,
				cmd
			);
			wi::renderer::BindCommonResources(cmd);

			RenderLightShafts(cmd);

			RenderTransparents(cmd);
		});

		RenderCameraComponents(ctx);

		CommandList cmd_postprocess = device->BeginCommandList();
		cmd = cmd_postprocess;
		wi::jobsystem::Execute(ctx, [this, cmd](wi::jobsystem::JobArgs args) {
			RenderPostprocessChain(cmd);

			GraphicsDevice* device = wi::graphics::GetDevice();
			device->Barrier(GPUBarrier::Image(&debugUAV, ResourceState::UNORDERED_ACCESS, debugUAV.desc.layout), cmd);
		});

		cmd = device->BeginCommandList(QUEUE_COPY);
		device->WaitCommandList(cmd, cmd_postprocess);
		wi::jobsystem::Execute(ctx, [this, cmd](wi::jobsystem::JobArgs args) {
			wi::renderer::TextureStreamingReadbackCopy(*scene, cmd);
		});

		RenderPath2D::Render();

		wi::jobsystem::Wait(ctx);

		// Every pass that draws trails leaves the queue alone, because they are
		// recorded on parallel jobs and the reflection and the main camera both
		// read it. This is the first point where they have all joined, so it is
		// the only safe place to discard them.
		wi::renderer::ClearTrails();

		first_frame = false;
	}

	void RenderPath3D::Compose(CommandList cmd) const
	{
		GraphicsDevice* device = wi::graphics::GetDevice();
		device->EventBegin("RenderPath3D::Compose", cmd);

		wi::image::Params fx;
		fx.blendFlag = BLENDMODE_OPAQUE;
		fx.quality = wi::image::QUALITY_LINEAR;
		fx.enableFullScreen();

		wi::image::Draw(&GetRenderResult3D(), fx, cmd);

		if (
			wi::renderer::GetDebugLightCulling() ||
			wi::renderer::GetVariableRateShadingClassificationDebug() ||
			wi::renderer::GetSurfelGIDebugEnabled()
			)
		{
			fx.enableFullScreen();
			fx.blendFlag = BLENDMODE_PREMULTIPLIED;
			wi::image::Draw(&debugUAV, fx, cmd);
		}

		device->EventEnd(cmd);

		RenderPath2D::Compose(cmd);
	}

	void RenderPath3D::Stop()
	{
		DeleteGPUResources();
	}

	void RenderPath3D::Start()
	{
		ResizeBuffers();
	}

	void RenderPath3D::RenderAO(CommandList cmd) const
	{
		if (rtAO.IsValid())
		{
			GetDevice()->Barrier(GPUBarrier::Aliasing(&rtParticleDistortion, &rtAO), cmd);
		}

		if (getAOEnabled())
		{
			switch (getAO())
			{
			case AO_SSAO:
				wi::renderer::Postprocess_SSAO(
					ssaoResources,
					rtAO,
					cmd,
					getAORange(),
					getAOSampleCount(),
					getAOPower()
				);
				break;
			case AO_HBAO:
				wi::renderer::Postprocess_HBAO(
					ssaoResources,
					*camera,
					rtAO,
					cmd,
					getAOPower()
				);
				break;
			case AO_MSAO:
				wi::renderer::Postprocess_MSAO(
					msaoResources,
					*camera,
					rtAO,
					cmd,
					getAOPower()
				);
				break;
			case AO_RTAO:
				wi::renderer::Postprocess_RTAO(
					rtaoResources,
					*scene,
					depthBuffer_Copy,
					rtAO,
					cmd,
					getAORange(),
					getAOPower()
				);
				break;
			case AO_DISABLED:
				break;
			}
		}
	}
	void RenderPath3D::RenderSSR(CommandList cmd) const
	{
		if (getSSREnabled() && !getRaytracedReflectionEnabled())
		{
			wi::renderer::Postprocess_SSR(
				ssrResources,
				rtSceneCopy,
				rtSSR,
				cmd,
				getReflectionRoughnessCutoff()
			);
		}
	}
	void RenderPath3D::RenderSSGI(CommandList cmd) const
	{
		if (getSSGIEnabled())
		{
			wi::renderer::Postprocess_SSGI(
				ssgiResources,
				rtSceneCopy,
				depthBuffer_Copy,
				visibilityResources.texture_normal_roughness,
				rtSSGI,
				cmd,
				getSSGIDepthRejection()
			);
		}
	}
	void RenderPath3D::RenderOutline(CommandList cmd) const
	{
		if (getOutlineEnabled())
		{
			wi::renderer::Postprocess_Outline(
				rtOutlineSource,
				cmd,
				getOutlineThreshold(),
				getOutlineThickness(),
				getOutlineColor()
			);
		}
	}
	void RenderPath3D::RenderLightShafts(CommandList cmd) const
	{
		const XMVECTOR sunDirection = XMLoadFloat3(&scene->weather.sunDirection);
		const float sunDotCamera = XMVectorGetX(XMVector3Dot(sunDirection, camera->GetAt()));

		if (getLightShaftsEnabled() && sunDotCamera > 0)
		{
			constexpr float fadeThreshold = 0.25f;

			// Calculate target fade factor based on sun-camera angle
			float targetFadeFactor = 0.0f;
			if (sunDotCamera > 0.25)
			{
				targetFadeFactor = 1.0f;
			}

			float fadeSpeed = getLightShaftsFadeSpeed();
			if (targetFadeFactor < lightShaftsFadeFactor)
			{
				// Adaptive fade-out: accelerate as we approach the cutoff threshold
				const float normalizedDistance = wi::math::saturate(sunDotCamera / fadeThreshold);
				constexpr float fadeOutMultiplier = 13.0f; // Multiplier for fast fade-out

				// When normalizedDistance is 1.0 (at threshold): slow fade (fadeSpeed)
				// When normalizedDistance is 0.0 (cutoff): very fast fade (fadeOutSpeedMax)
				fadeSpeed = wi::math::Lerp(fadeSpeed * fadeOutMultiplier, fadeSpeed, normalizedDistance);
			}

			lightShaftsFadeFactor = wi::math::Lerp(lightShaftsFadeFactor, targetFadeFactor, 1.0f - exp(-fadeSpeed * scene->dt));

			GraphicsDevice* device = wi::graphics::GetDevice();

			device->EventBegin("Light Shafts", cmd);

			const Texture* texture_fullres = nullptr;

			// Render sun stencil cutout:
			{
				if (getMSAASampleCount() > 1)
				{
					RenderPassImage rp[] = {
						RenderPassImage::RenderTarget(&rtSun[0], RenderPassImage::LoadOp::CLEAR, RenderPassImage::StoreOp::DONTCARE),
						RenderPassImage::Resolve(&rtSun_resolved),
						RenderPassImage::DepthStencil(
							&depthBuffer_Main,
							RenderPassImage::LoadOp::LOAD,
							RenderPassImage::StoreOp::STORE,
							ResourceState::DEPTHSTENCIL,
							ResourceState::DEPTHSTENCIL,
							ResourceState::DEPTHSTENCIL
						),
					};
					device->RenderPassBegin(rp, arraysize(rp), cmd);
					texture_fullres = &rtSun_resolved;
				}
				else
				{
					RenderPassImage rp[] = {
						RenderPassImage::DepthStencil(
							&depthBuffer_Main,
							RenderPassImage::LoadOp::LOAD,
							RenderPassImage::StoreOp::STORE,
							ResourceState::DEPTHSTENCIL,
							ResourceState::DEPTHSTENCIL,
							ResourceState::DEPTHSTENCIL
						),
						RenderPassImage::RenderTarget(&rtSun[0], RenderPassImage::LoadOp::CLEAR),
					};
					device->RenderPassBegin(rp, arraysize(rp), cmd);
					texture_fullres = &rtSun[0];
				}

				Viewport vp;
				vp.width = (float)depthBuffer_Main.GetDesc().width;
				vp.height = (float)depthBuffer_Main.GetDesc().height;
				device->BindViewports(1, &vp, cmd);

				wi::graphics::Rect scissor = GetScissorInternalResolution();
				device->BindScissorRects(1, &scissor, cmd);

				wi::renderer::DrawSun(cmd);

				device->RenderPassEnd(cmd);
			}

			// Radial blur on the sun:
			{
				XMVECTOR sunPos = XMVector3Project(camera->GetEye() + sunDirection * camera->zFarP, 0, 0,
					1.0f, 1.0f, 0.1f, 1.0f,
					camera->GetProjection(), camera->GetView(), XMMatrixIdentity());
				{
					// Downsample to low res first:
					wi::renderer::Postprocess_Downsample4x(*texture_fullres, rtSun[2], cmd);

					XMFLOAT2 sun;
					XMStoreFloat2(&sun, sunPos);
					wi::renderer::Postprocess_LightShafts(
						rtSun[2],
						rtSun[1],
						cmd,
						sun,
						getLightShaftsStrength()
					);
				}
			}
			device->EventEnd(cmd);
		}
	}
	void RenderPath3D::RenderSceneMIPChain(CommandList cmd) const
	{
		GraphicsDevice* device = wi::graphics::GetDevice();

		auto range = wi::profiler::BeginRangeGPU("Scene MIP Chain", cmd);
		device->EventBegin("RenderSceneMIPChain", cmd);

		{
			GPUBarrier barriers[] = {
				GPUBarrier::Aliasing(&rtPrimitiveID, &rtSceneCopy_tmp),
				GPUBarrier::Image(&rtSceneCopy_tmp, rtSceneCopy_tmp.desc.layout, ResourceState::UNORDERED_ACCESS),
			};
			device->Barrier(barriers, arraysize(barriers), cmd);
			device->ClearUAV(&rtSceneCopy_tmp, 0, cmd);
		}

		wi::renderer::Postprocess_Downsample4x(rtMain, rtSceneCopy, cmd);

		device->Barrier(GPUBarrier::Image(&rtSceneCopy_tmp, ResourceState::UNORDERED_ACCESS, rtSceneCopy_tmp.desc.layout), cmd);

		wi::renderer::MIPGEN_OPTIONS mipopt;
		mipopt.gaussian_temp = &rtSceneCopy_tmp;
		wi::renderer::GenerateMipChain(rtSceneCopy, wi::renderer::MIPGENFILTER_GAUSSIAN, cmd, mipopt);

		device->Barrier(GPUBarrier::Aliasing(&rtSceneCopy_tmp, &rtPrimitiveID), cmd);

		device->EventEnd(cmd);
		wi::profiler::EndRange(range);
	}
	void RenderPath3D::RenderTransparents(CommandList cmd) const
	{
		GraphicsDevice* device = wi::graphics::GetDevice();

		// Water ripple rendering:
		if (!scene->waterRipples.empty())
		{
			device->Barrier(GPUBarrier::Aliasing(&rtParticleDistortion, &rtWaterRipple), cmd);
			RenderPassImage rp[] = {
				RenderPassImage::RenderTarget(&rtWaterRipple, RenderPassImage::LoadOp::CLEAR),
			};
			device->RenderPassBegin(rp, arraysize(rp), cmd);

			Viewport vp;
			vp.width = (float)rtWaterRipple.GetDesc().width;
			vp.height = (float)rtWaterRipple.GetDesc().height;
			device->BindViewports(1, &vp, cmd);

			wi::renderer::DrawWaterRipples(visibility_main, cmd);
			device->RenderPassEnd(cmd);
		}

		if (getFSR2Enabled())
		{
			// Save the pre-alpha for FSR2 reactive mask:
			//	Note that rtFSR temp resource is always larger or equal to rtMain, so CopyTexture is used instead of CopyResource!
			GPUBarrier barriers[] = {
				GPUBarrier::Image(&rtMain, rtMain.desc.layout, ResourceState::COPY_SRC),
				GPUBarrier::Image(&rtFSR[1], rtFSR->desc.layout, ResourceState::COPY_DST),
			};
			device->Barrier(barriers, arraysize(barriers), cmd);
			device->CopyTexture(
				&rtFSR[1], 0, 0, 0, 0, 0,
				&rtMain, 0, 0,
				cmd
			);
			for (int i = 0; i < arraysize(barriers); ++i)
			{
				std::swap(barriers[i].image.layout_before, barriers[i].image.layout_after);
			}
			device->Barrier(barriers, arraysize(barriers), cmd);
		}

		wi::renderer::UpdateGaussianSplatsForCamera(*scene, *camera, cmd);

		wi::graphics::Rect scissor = GetScissorInternalResolution();
		device->BindScissorRects(1, &scissor, cmd);

		Viewport vp;
		vp.width = (float)depthBuffer_Main.GetDesc().width;
		vp.height = (float)depthBuffer_Main.GetDesc().height;
		vp.min_depth = 0;
		vp.max_depth = 1;
		device->BindViewports(1, &vp, cmd);

		RenderPassImage rp[3];
		uint32_t rp_count = 0;
		if (getMSAASampleCount() > 1)
		{
			rp[rp_count++] = RenderPassImage::RenderTarget(&rtMain_render, RenderPassImage::LoadOp::LOAD);
			rp[rp_count++] = RenderPassImage::Resolve(&rtMain);
			rp[rp_count++] = RenderPassImage::DepthStencil(
				&depthBuffer_Main,
				RenderPassImage::LoadOp::LOAD,
				RenderPassImage::StoreOp::STORE,
				ResourceState::DEPTHSTENCIL,
				ResourceState::DEPTHSTENCIL,
				ResourceState::DEPTHSTENCIL
			);
		}
		else
		{

			rp[rp_count++] = RenderPassImage::RenderTarget(&rtMain_render, RenderPassImage::LoadOp::LOAD);
			rp[rp_count++] = RenderPassImage::DepthStencil(
				&depthBuffer_Main,
				RenderPassImage::LoadOp::LOAD,
				RenderPassImage::StoreOp::STORE,
				ResourceState::DEPTHSTENCIL,
				ResourceState::DEPTHSTENCIL,
				ResourceState::DEPTHSTENCIL
			);
		}

		// The ocean is the only transparent that writes depth, so whatever is
		// drawn after it and lies beyond it is rejected outright, and whatever
		// is not yet in the scene copy taken just below can never be seen
		// through the water. Both halves of the transparent pass therefore have
		// to be placed relative to the ocean by which side of the water plane
		// they are on - and which side is which flips with the camera, since
		// from below the surface the submerged content is the near side.
		const bool oceanVisible =
			scene->weather.IsOceanEnabled() &&
			scene->ocean.IsValid() &&
			(!scene->ocean.IsOccluded() || !wi::renderer::GetOcclusionCullingEnabled());
		const bool eyeAboveWater =
			camera->Eye.y > scene->weather.oceanParameters.waterHeight;
		const wi::renderer::WATERSIDE farSide = !oceanVisible ? wi::renderer::WATERSIDE_ALL :
			(eyeAboveWater ? wi::renderer::WATERSIDE_SUBMERGED : wi::renderer::WATERSIDE_ABOVE);
		const wi::renderer::WATERSIDE nearSide = !oceanVisible ? wi::renderer::WATERSIDE_ALL :
			(eyeAboveWater ? wi::renderer::WATERSIDE_ABOVE : wi::renderer::WATERSIDE_SUBMERGED);

		// Restricts every following draw to one side of the water, by putting the
		// side on the camera instead of plumbing it through each draw call. The
		// shaders that split read it straight out of the camera constant buffer.
		//
		// Bracket it tightly and always restore it: the diagnostic overlays,
		// gaussian splats and the screen space effects share these passes and
		// must not be clipped by it.
		auto bindWaterSide = [&](wi::renderer::WATERSIDE side) {
			CameraComponent sidedCamera = *camera;
			if (side == wi::renderer::WATERSIDE_SUBMERGED)
			{
				sidedCamera.shadercamera_options |= SHADERCAMERA_OPTION_WATERSIDE_SUBMERGED;
			}
			else if (side == wi::renderer::WATERSIDE_ABOVE)
			{
				sidedCamera.shadercamera_options |= SHADERCAMERA_OPTION_WATERSIDE_ABOVE;
			}
			wi::renderer::BindCameraCB(
				sidedCamera,
				camera_previous,
				camera_reflection,
				cmd
			);
		};

		// Draw only the ocean first, fog and lightshafts will be blended on top:
		if (oceanVisible)
		{
			// The far side goes down before the copy, so that the ocean's
			// screen space refraction shows it through the surface instead of
			// the surface painting over it, and so the surface's depth write
			// cannot reject it outright the way it did before.
			//
			// Everything here is drawn a second time in the near pass below,
			// where it keeps the other half. That is what lets a thing spanning
			// the waterline be partly submerged rather than snapping between
			// the two sides.
			//
			// The far pass costs a scene mip chain, a render pass and a redraw
			// of every content type, and in the common case - camera above the
			// water, nothing below it - none of that draws a single pixel. So
			// sweep the scene first and find out what is actually over there.
			//
			// Every bound below is conservative: anything that might straddle
			// the surface counts as far side and is drawn in both passes, where
			// the per-fragment clip sorts it out. So this can only ever skip
			// work, never content.
			//
			// The band the surface sweeps is taken from the ocean's own bounds,
			// so there is one assumption about how far the waves displace
			// vertically rather than a second one here that could drift from
			// it.
			const wi::primitive::AABB oceanBounds = scene->ocean.GetAABB(camera->Eye);
			const float waterTop = oceanBounds.getMax().y;
			const float waterBottom = oceanBounds.getMin().y;
			auto reachesFarSide = [&](const wi::primitive::AABB& aabb) {
				return farSide == wi::renderer::WATERSIDE_SUBMERGED
					? (aabb.getMin().y <= waterTop)   // any part low enough to be under water
					: (aabb.getMax().y > waterBottom); // any part high enough to be above it
			};

			bool farSideTransparents = false;
			if (visibility_main.IsTransparentsVisible())
			{
				for (uint32_t index : visibility_main.visibleObjects)
				{
					const ObjectComponent& object = scene->objects[index];
					if (!object.IsRenderable())
						continue;
					// The far pass has no foreground draw, so a foreground
					// object is never in it - matching that here rather than
					// merely being conservative about it.
					if (object.IsForeground())
						continue;
					if ((object.GetFilterMask() & wi::enums::FILTER_TRANSPARENT) == 0)
						continue;
					if (reachesFarSide(scene->aabb_objects[index]))
					{
						farSideTransparents = true;
						break;
					}
				}
			}

			bool farSideLights = false;
			for (uint32_t index : visibility_main.visibleLights)
			{
				// The light's whole influence volume, not its visualizer gizmo,
				// and every visible light rather than only the visualized ones.
				// Both make this wider than it needs to be, which is the safe
				// direction.
				if (reachesFarSide(scene->aabb_lights[index]))
				{
					farSideLights = true;
					break;
				}
			}

			bool farSideSplats = false;
			for (size_t i = 0; i < scene->gaussian_splats.GetCount(); ++i)
			{
				if (reachesFarSide(scene->gaussian_splats[i].aabb))
				{
					farSideSplats = true;
					break;
				}
			}

			// Sprites have no world bounds to test, so any sprite at all counts
			// as far side. Fonts do have bounds, which is worth using for the
			// scenes that have text but no sprites.
			bool farSideSpritesOrFonts = scene->sprites.GetCount() > 0;
			if (!farSideSpritesOrFonts)
			{
				for (size_t i = 0; i < scene->aabb_fonts.size(); ++i)
				{
					if (reachesFarSide(scene->aabb_fonts[i]))
					{
						farSideSpritesOrFonts = true;
						break;
					}
				}
			}

			// Particles leave their emitter and travel where they like, so the
			// emitter's position bounds nothing. Any visible emitter counts.
			const bool farSideEmitters = !visibility_main.visibleEmitters.empty();

			// Trails are a per-frame queue with no bounds to test, so any trail
			// at all counts, as with sprites and emitters above.
			const bool farSideTrails = wi::renderer::AreTrailsQueued();

			// The suspended particles exist only below the surface, so they are
			// far side content exactly when the far side is the submerged one -
			// with the camera under water the far side is the air above it,
			// where every particle is clipped away. Beyond that the only
			// question is whether the field, which is centred on the camera,
			// reaches down as far as the water at all.
			const bool farSideParticles =
				getUnderwaterParticlesEnabled() &&
				farSide == wi::renderer::WATERSIDE_SUBMERGED &&
				camera->Eye.y
					- (wi::scene::environment::UnderwaterParticles(
						scene->weather.oceanParameters.waterMedium,
						getUnderwaterParticleDensity())
							.FieldSize() * 0.5f) <= waterTop;

			const bool anyFarSideContent =
				farSideTransparents ||
				farSideLights ||
				farSideSplats ||
				farSideSpritesOrFonts ||
				farSideEmitters ||
				farSideTrails ||
				farSideParticles;

			if (anyFarSideContent)
			{
				// A refractive mesh on the far side samples the scene copy, and
				// at this point it still holds the previous frame's image - the
				// only fill so far is the one taken at the end of the last
				// RenderTransparents for SSR. Refresh it, or a submerged glass
				// panel refracts a frame-old scene. The chain is built again
				// after the ocean below, for the near side, which is the price
				// of drawing the two sides in two places.
				if (farSideTransparents)
				{
					RenderSceneMIPChain(cmd);
				}

				device->EventBegin("Transparents beyond the water", cmd);
				device->RenderPassBegin(rp, rp_count, cmd);
				bindWaterSide(farSide);

				if (farSideTransparents)
				{
					// No foreground pass: a foreground object is drawn in a
					// compressed depth range right against the camera, so it is
					// always in front of the water and never on the far side.
					wi::renderer::DrawScene(
						visibility_main,
						RENDERPASS_MAIN,
						cmd,
						wi::renderer::DRAWSCENE_TRANSPARENT |
						wi::renderer::DRAWSCENE_TESSELLATION |
						wi::renderer::DRAWSCENE_OCCLUSIONCULLING |
						wi::renderer::DRAWSCENE_MAINCAMERA
					);
				}

				if (farSideLights)
				{
					wi::renderer::DrawLightVisualizers(visibility_main, cmd);
				}
				if (farSideEmitters)
				{
					wi::renderer::DrawSoftParticles(visibility_main, false, cmd);
				}
				if (farSideParticles)
				{
					wi::renderer::DrawUnderwaterParticles(
						visibility_main, getUnderwaterParticleDensity(), cmd);
				}
				if (farSideSplats)
				{
					wi::renderer::DrawGaussianSplats(*scene, *camera, cmd);
				}
				if (farSideTrails)
				{
					// Keep the queue: the near pass below draws the same trails
					// again for its own half, and consumes them then.
					wi::renderer::DrawTrails(*camera, cmd, false);
				}
				if (farSideSpritesOrFonts)
				{
					wi::renderer::DrawSpritesAndFonts(*scene, *camera, false, cmd, farSide);
				}

				bindWaterSide(wi::renderer::WATERSIDE_ALL);
				device->RenderPassEnd(cmd);
				device->EventEnd(cmd);
			}

			device->EventBegin("Copy scene tex only mip0 for ocean", cmd);
			wi::renderer::Postprocess_Downsample4x(rtMain, rtSceneCopy, cmd);
			device->EventEnd(cmd);

			// The ocean writes depth so that its own displaced grid can occlude
			// itself, but nothing drawn after it may be depth-tested against
			// the water - the transparent pass is ordered by water side
			// instead, and the editor overlays that come later are not ordered
			// at all. So it draws against a scratch copy of the opaque depth,
			// taken here: the water is still occluded by geometry in front of
			// it and by itself, and its writes leave with the copy.
			//
			// If the scratch buffer is missing the ocean falls back to the main
			// depth buffer, which is the old behaviour: correct water, but its
			// depth then erases whatever is drawn behind it afterwards.
			const bool oceanScratchDepth = depthBuffer_Ocean.IsValid();
			if (oceanScratchDepth)
			{
				device->EventBegin("Copy depth for ocean", cmd);
				GPUBarrier barriers[] = {
					GPUBarrier::Image(&depthBuffer_Main, depthBuffer_Main.desc.layout, ResourceState::COPY_SRC),
					GPUBarrier::Image(&depthBuffer_Ocean, depthBuffer_Ocean.desc.layout, ResourceState::COPY_DST),
				};
				device->Barrier(barriers, arraysize(barriers), cmd);
				device->CopyResource(&depthBuffer_Ocean, &depthBuffer_Main, cmd);
				for (int i = 0; i < arraysize(barriers); ++i)
				{
					std::swap(barriers[i].image.layout_before, barriers[i].image.layout_after);
				}
				device->Barrier(barriers, arraysize(barriers), cmd);
				device->EventEnd(cmd);
			}

			RenderPassImage rp_ocean[3];
			uint32_t rp_ocean_count = 0;
			rp_ocean[rp_ocean_count++] = RenderPassImage::RenderTarget(&rtMain_render, RenderPassImage::LoadOp::LOAD);
			if (getMSAASampleCount() > 1)
			{
				rp_ocean[rp_ocean_count++] = RenderPassImage::Resolve(&rtMain);
			}
			rp_ocean[rp_ocean_count++] = RenderPassImage::DepthStencil(
				oceanScratchDepth ? &depthBuffer_Ocean : &depthBuffer_Main,
				RenderPassImage::LoadOp::LOAD,
				// Nothing reads the scratch depth afterwards, so let a tiler
				// skip writing it back. The main buffer must still be kept.
				oceanScratchDepth ? RenderPassImage::StoreOp::DONTCARE : RenderPassImage::StoreOp::STORE,
				ResourceState::DEPTHSTENCIL,
				ResourceState::DEPTHSTENCIL,
				ResourceState::DEPTHSTENCIL
			);

			device->RenderPassBegin(rp_ocean, rp_ocean_count, cmd);

			wi::renderer::DrawScene(
				visibility_main,
				RENDERPASS_MAIN,
				cmd,
				wi::renderer::DRAWSCENE_OCEAN
			);

			device->RenderPassEnd(cmd);
		}

		if (visibility_main.IsTransparentsVisible())
		{
			RenderSceneMIPChain(cmd);
		}

		device->RenderPassBegin(rp, rp_count, cmd);

		// Transparent scene
		if (visibility_main.IsTransparentsVisible())
		{
			auto range = wi::profiler::BeginRangeGPU("Transparent Scene", cmd);
			device->EventBegin("Transparent Scene", cmd);

			// Regular: the near side only, its far half having been drawn
			// before the ocean above. Restored immediately after, because the
			// foreground pass below is always in front of the water.
			vp.min_depth = 0;
			vp.max_depth = 1;
			device->BindViewports(1, &vp, cmd);
			if (oceanVisible)
			{
				bindWaterSide(nearSide);
			}
			wi::renderer::DrawScene(
				visibility_main,
				RENDERPASS_MAIN,
				cmd,
				wi::renderer::DRAWSCENE_TRANSPARENT |
				wi::renderer::DRAWSCENE_TESSELLATION |
				wi::renderer::DRAWSCENE_OCCLUSIONCULLING |
				wi::renderer::DRAWSCENE_MAINCAMERA
			);
			if (oceanVisible)
			{
				bindWaterSide(wi::renderer::WATERSIDE_ALL);
			}

			// Foreground:
			vp.min_depth = 1 - foreground_depth_range;
			vp.max_depth = 1;
			device->BindViewports(1, &vp, cmd);
			wi::renderer::DrawScene(
				visibility_main,
				RENDERPASS_MAIN,
				cmd,
				wi::renderer::DRAWSCENE_TRANSPARENT |
				wi::renderer::DRAWSCENE_FOREGROUND_ONLY |
				wi::renderer::DRAWSCENE_MAINCAMERA
			);

			// Reset normal viewport:
			vp.min_depth = 0;
			vp.max_depth = 1;
			device->BindViewports(1, &vp, cmd);

			device->EventEnd(cmd);
			wi::profiler::EndRange(range); // Transparent Scene
		}

		wi::renderer::DrawDebugWorld(*scene, *camera, *this, cmd);

		wi::renderer::DrawWireframeOverlay(visibility_main, wi::enums::RENDERPASS_MAIN, cmd);

		// The near side only for these: the far side was drawn before the ocean
		// above, so that the water refracts it rather than rejecting it.
		//
		// The splats are sorted once per frame, for the main camera, and both
		// passes draw that same sorted set - they differ only in which half
		// each keeps, so back to front order still holds within either half.
		//
		// No pass consumes the trails, this one included: the planar reflection
		// draws them too, on a job that records in parallel with this one, so
		// clearing here would race it. Render() clears once the jobs join.
		if (oceanVisible)
		{
			bindWaterSide(nearSide);
		}
		wi::renderer::DrawLightVisualizers(visibility_main, cmd);
		wi::renderer::DrawSoftParticles(visibility_main, false, cmd);
		// Unlike the other content here the particles belong to one side of the
		// water outright, so they are drawn once rather than in both halves.
		// With the camera above the surface the near side is the air, holding
		// none of them, and the field has already gone down in the far pass
		// before the ocean - where the surface refracts it, as it should.
		if (getUnderwaterParticlesEnabled() &&
			nearSide != wi::renderer::WATERSIDE_ABOVE)
		{
			wi::renderer::DrawUnderwaterParticles(
				visibility_main, getUnderwaterParticleDensity(), cmd);
		}
		wi::renderer::DrawGaussianSplats(*scene, *camera, cmd);
		wi::renderer::DrawTrails(*camera, cmd, false);
		if (oceanVisible)
		{
			bindWaterSide(wi::renderer::WATERSIDE_ALL);
		}

		wi::renderer::DrawSpritesAndFonts(*scene, *camera, false, cmd, nearSide);

		XMVECTOR sunDirection = XMLoadFloat3(&scene->weather.sunDirection);
		if (getLightShaftsEnabled() && XMVectorGetX(XMVector3Dot(sunDirection, camera->GetAt())) > 0)
		{
			// The sun cutout these shafts are blurred from is drawn before the
			// ocean is, so the water surface is not in the depth buffer that
			// occludes it. Submerged, that paints the sun at full strength as
			// if the camera were in air: no refraction into Snell's window, no
			// attenuation by the water column, and the radial blur smears it
			// over whatever it passes - which is what lights up the underside
			// of anything floating on the surface. The underwater pass owns the
			// sun's appearance from down there.
			float lightShaftsWaterFade = 1;
			if (scene->weather.IsOceanEnabled())
			{
				// Fade across a wave's worth of height rather than snapping at
				// the still water level, so a camera bobbing at the waterline
				// does not flicker.
				const float waterHeight = scene->weather.oceanParameters.waterHeight;
				lightShaftsWaterFade =
					wi::math::SmoothStep(0.0f, 1.0f, camera->Eye.y - waterHeight);
			}

			if (lightShaftsWaterFade > 0)
			{
				device->EventBegin("Contribute LightShafts", cmd);
				wi::image::Params fx;
				fx.enableFullScreen();
				fx.blendFlag = BLENDMODE_ADDITIVE;
				fx.opacity = lightShaftsFadeFactor * lightShaftsWaterFade;
				wi::image::Draw(&rtSun[1], fx, cmd);
				device->EventEnd(cmd);
			}
		}

		if (getLensFlareEnabled())
		{
			// The window decides where a submerged eye SEES the sun, and the
			// flare has to be put where the sun looks like it is. The same
			// expression the underwater pass is given below, so the two cannot
			// disagree about how much compression there is.
			wi::renderer::DrawLensFlares(
				visibility_main,
				cmd,
				getUnderwaterSnellEnabled()
					? getUnderwaterSnellStrength()
					: 0.0f
			);
		}

		device->RenderPassEnd(cmd);

		// Distortion particles:
		{
			if (rtWaterRipple.IsValid())
			{
				device->Barrier(GPUBarrier::Aliasing(&rtWaterRipple, &rtParticleDistortion), cmd);
			}

			if (getMSAASampleCount() > 1)
			{
				RenderPassImage rp[] = {
					RenderPassImage::RenderTarget(&rtParticleDistortion_render, RenderPassImage::LoadOp::CLEAR),
					RenderPassImage::Resolve(&rtParticleDistortion),
					RenderPassImage::DepthStencil(
						&depthBuffer_Main,
						RenderPassImage::LoadOp::LOAD,
						RenderPassImage::StoreOp::STORE,
						ResourceState::DEPTHSTENCIL,
						ResourceState::DEPTHSTENCIL,
						ResourceState::DEPTHSTENCIL
					),
				};
				device->RenderPassBegin(rp, arraysize(rp), cmd);
			}
			else
			{
				RenderPassImage rp[] = {
					RenderPassImage::RenderTarget(&rtParticleDistortion, RenderPassImage::LoadOp::CLEAR),
					RenderPassImage::DepthStencil(
						&depthBuffer_Main,
						RenderPassImage::LoadOp::LOAD,
						RenderPassImage::StoreOp::STORE,
						ResourceState::DEPTHSTENCIL,
						ResourceState::DEPTHSTENCIL,
						ResourceState::DEPTHSTENCIL
					),
				};
				device->RenderPassBegin(rp, arraysize(rp), cmd);
			}

			Viewport vp;
			vp.width = (float)rtParticleDistortion.GetDesc().width;
			vp.height = (float)rtParticleDistortion.GetDesc().height;
			device->BindViewports(1, &vp, cmd);

			wi::renderer::DrawSoftParticles(visibility_main, true, cmd);
			wi::renderer::DrawSpritesAndFonts(*scene, *camera, true, cmd);

			device->RenderPassEnd(cmd);
		}

		wi::renderer::Postprocess_Downsample4x(rtMain, rtSceneCopy, cmd);
	}
	void RenderPath3D::RenderPostprocessChain(CommandList cmd) const
	{
		GraphicsDevice* device = wi::graphics::GetDevice();

		wi::renderer::BindCommonResources(cmd);
		wi::renderer::BindCameraCB(*camera, camera_previous, camera_reflection, cmd);

		const Texture* rt_first = nullptr; // not ping-ponged with read / write
		const Texture* rt_read = &rtMain;
		const Texture* rt_write = &rtPostprocess;

		// rtPostprocess aliasing transition:
		device->Barrier(GPUBarrier::Aliasing(&rtPrimitiveID, &rtPostprocess), cmd);

		// 1.) HDR post process chain
		{
			if (getFSR2Enabled() && fsr2Resources.IsValid())
			{
				wi::renderer::Postprocess_FSR2(
					fsr2Resources,
					*camera,
					rtFSR[1],
					*rt_read,
					depthBuffer_Copy,
					rtVelocity,
					rtFSR[0],
					cmd,
					scene->dt,
					getFSR2Sharpness()
				);

				// rebind these, because FSR2 binds other things to those constant buffers:
				wi::renderer::BindCameraCB(
					*camera,
					camera_previous,
					camera_reflection,
					cmd
				);
				wi::renderer::BindCommonResources(cmd);

				rt_read = &rtFSR[0];
				rt_write = &rtFSR[1];
			}
			else if (wi::renderer::GetTemporalAAEnabled() && !wi::renderer::GetTemporalAADebugEnabled() && temporalAAResources.IsValid() && wi::renderer::GetWireframeMode() == wi::renderer::WIREFRAME_DISABLED)
			{
				wi::renderer::Postprocess_TemporalAA(
					temporalAAResources,
					*rt_read,
					cmd
				);
				rt_first = temporalAAResources.GetCurrent();
			}
			if (scene->weather.IsOceanEnabled())
			{
				// Underwater magnification: refraction at the eye makes
				// submerged objects appear ~33% larger. The zoom is applied
				// inside the underwater pass and is naturally masked to below
				// the waterline by its final blend, so a partially submerged
				// view stays true scale above the surface.
				const float underwater_magnification =
					getUnderwaterMagnificationEnabled()
					? getUnderwaterMagnification()
					: 1.0f;
				const float underwater_snell =
					getUnderwaterSnellEnabled()
					? getUnderwaterSnellStrength()
					: 0.0f;
				// Ray-traced window: traces the refracted ray into the real
				// scene when the device supports it (the underwaterCS_rtapi
				// permutation is loaded then); otherwise the window shows the
				// analytic sky:
				const bool underwater_snell_rt =
					getUnderwaterSnellEnabled()
					&& getUnderwaterSnellRTEnabled()
					&& device->CheckCapability(wi::graphics::GraphicsDeviceCapability::RAYTRACING);
				wi::renderer::Postprocess_Underwater(
					rt_first == nullptr ? *rt_read : *rt_first,
					*rt_write,
					cmd,
					underwater_magnification,
					underwater_snell,
					getUnderwaterSnellFade(),
					underwater_snell_rt,
					getUnderwaterLensDistortionEnabled() ? 1.0f : 0.0f
				);

				rt_first = nullptr;
				std::swap(rt_read, rt_write);

				// Underwater vision blur: the human eye cannot focus under
				// water, so reuse the depth-of-field pipeline in underwater
				// mode. Its circle of confusion is masked to the submerged part
				// of the screen, so a partially submerged camera stays sharp
				// above the waterline. The DoF tile early-exit keeps this cheap
				// when the camera is above water.
				if (getUnderwaterBlurEnabled() && depthoffieldResources.IsValid())
				{
					const float underwater_coc = getUnderwaterBlurStrength() * 6.0f;
					wi::renderer::Postprocess_DepthOfField(
						depthoffieldResources,
						rt_first == nullptr ? *rt_read : *rt_first,
						*rt_write,
						cmd,
						underwater_coc, // coc scale
						underwater_coc, // max coc
						true // underwater mode
					);

					rt_first = nullptr;
					std::swap(rt_read, rt_write);
				}

				// A little chromatic dispersion reads as light refracting
				// through water. Masked to the submerged part of the screen.
				if (getUnderwaterChromaticAberrationEnabled() && getUnderwaterChromaticAberrationAmount() > 0)
				{
					wi::renderer::Postprocess_Chromatic_Aberration(
						rt_first == nullptr ? *rt_read : *rt_first,
						*rt_write,
						cmd,
						getUnderwaterChromaticAberrationAmount(),
						true // underwater mode (mask to submerged region)
					);

					rt_first = nullptr;
					std::swap(rt_read, rt_write);
				}
			}

			for (auto& x : custom_post_processes)
			{
				if (x.stage == CustomPostprocess::Stage::BeforeTonemap)
				{
					wi::renderer::Postprocess_Custom(
						x.computeshader,
						rt_first == nullptr ? *rt_read : *rt_first,
						*rt_write,
						cmd,
						x.params0,
						x.params1,
						x.name.c_str()
					);

					rt_first = nullptr;
					std::swap(rt_read, rt_write);
				}
			}

			if (getDepthOfFieldEnabled() && camera->aperture_size > 0.001f && getDepthOfFieldStrength() > 0.001f && depthoffieldResources.IsValid())
			{
				wi::renderer::Postprocess_DepthOfField(
					depthoffieldResources,
					rt_first == nullptr ? *rt_read : *rt_first,
					*rt_write,
					cmd,
					getDepthOfFieldStrength()
				);

				rt_first = nullptr;
				std::swap(rt_read, rt_write);
			}

			if (getMotionBlurEnabled() && getMotionBlurStrength() > 0 && motionblurResources.IsValid())
			{
				wi::renderer::Postprocess_MotionBlur(
					scene->dt,
					motionblurResources,
					rt_first == nullptr ? *rt_read : *rt_first,
					*rt_write,
					cmd,
					getMotionBlurStrength()
				);

				rt_first = nullptr;
				std::swap(rt_read, rt_write);
			}
		}

		// 2.) Tone mapping HDR -> LDR
		{
			// Bloom and eye adaption is not part of post process "chain",
			//	because they will be applied to the screen in tonemap
			if (getEyeAdaptionEnabled())
			{
				wi::renderer::ComputeLuminance(
					luminanceResources,
					rt_first == nullptr ? *rt_read : *rt_first,
					cmd,
					getEyeAdaptionRate(),
					getEyeAdaptionKey()
				);
			}
			if (getBloomEnabled())
			{
				wi::renderer::ComputeBloom(
					bloomResources,
					rt_first == nullptr ? *rt_read : *rt_first,
					cmd,
					getBloomThreshold(),
					getExposure(),
					getEyeAdaptionEnabled() ? &luminanceResources.luminance : nullptr
				);
			}

			wi::renderer::Postprocess_Tonemap(
				rt_first == nullptr ? *rt_read : *rt_first,
				*rt_write,
				cmd,
				getExposure(),
				getBrightness(),
				getContrast(),
				getSaturation(),
				getDitherEnabled(),
				getColorGradingEnabled() ? (scene->weather.colorGradingMap.IsValid() ? &scene->weather.colorGradingMap.GetTexture() : nullptr) : nullptr,
				&rtParticleDistortion,
				getEyeAdaptionEnabled() ? &luminanceResources.luminance : nullptr,
				getBloomEnabled() ? &bloomResources.texture_bloom : nullptr,
				colorspace,
				getTonemap(),
				&distortion_overlay,
				getHDRCalibration()
			);

			rt_first = nullptr;
			std::swap(rt_read, rt_write);
		}

		// 3.) LDR post process chain
		{
			for (auto& x : custom_post_processes)
			{
				if (x.stage == CustomPostprocess::Stage::AfterTonemap)
				{
					wi::renderer::Postprocess_Custom(
						x.computeshader,
						*rt_read,
						*rt_write,
						cmd,
						x.params0,
						x.params1,
						x.name.c_str()
					);

					std::swap(rt_read, rt_write);
				}
			}

			if (getSharpenFilterEnabled())
			{
				wi::renderer::Postprocess_Sharpen(*rt_read, *rt_write, cmd, getSharpenFilterAmount());

				std::swap(rt_read, rt_write);
			}

			if (getFXAAEnabled())
			{
				wi::renderer::Postprocess_FXAA(*rt_read, *rt_write, cmd);

				std::swap(rt_read, rt_write);
			}

			if (getChromaticAberrationEnabled())
			{
				wi::renderer::Postprocess_Chromatic_Aberration(*rt_read, *rt_write, cmd, getChromaticAberrationAmount());

				std::swap(rt_read, rt_write);
			}

			if (getCRTFilterEnabled())
			{
				wi::renderer::Postprocess_CRT(*rt_read, *rt_write, cmd, 0, 0, true);

				std::swap(rt_read, rt_write);
			}

			lastPostprocessRT = rt_read;

			// GUI Background blurring:
			{
				auto range = wi::profiler::BeginRangeGPU("GUI Background Blur", cmd);
				device->EventBegin("GUI Background Blur", cmd);
				bool hdrToSRGB = colorspace != ColorSpace::SRGB;
				wi::renderer::Postprocess_Downsample4x(*rt_read, rtGUIBlurredBackground[0], cmd, hdrToSRGB);
				wi::renderer::Postprocess_Downsample4x(rtGUIBlurredBackground[0], rtGUIBlurredBackground[2], cmd);
				wi::renderer::Postprocess_Blur_Gaussian(rtGUIBlurredBackground[2], rtGUIBlurredBackground[1], rtGUIBlurredBackground[2], cmd, -1, -1, true);
				device->EventEnd(cmd);
				wi::profiler::EndRange(range);
			}

			if (rtFSR[0].IsValid() && getFSREnabled())
			{
				wi::renderer::Postprocess_FSR(*rt_read, rtFSR[1], rtFSR[0], cmd, getFSRSharpness());
				lastPostprocessRT = &rtFSR[0];
			}
		}
	}

	void RenderPath3D::RenderCameraComponents(wi::jobsystem::context& ctx) const
	{
		// Render-to-texture camera components:
		for (uint32_t i = 0; i < scene->cameras.GetCount() && getSceneUpdateEnabled(); ++i)
		{
			wi::scene::CameraComponent& camera = scene->cameras[i];
			if (camera.render_to_texture.resolution.x == 0 || camera.render_to_texture.resolution.y == 0)
			{
				camera.render_to_texture = {};
				continue;
			}

			camera.render_to_texture.time_accumulator += scene->dt;
			bool should_render = (camera.render_to_texture.update_interval <= 0.0f) || 
				(camera.render_to_texture.time_accumulator >= camera.render_to_texture.update_interval);
			if (!should_render)
			{
				continue;
			}
			camera.render_to_texture.time_accumulator = 0.0f;

			GraphicsDevice* device = GetDevice();
			CommandList cmd = device->BeginCommandList();

			if (!camera.render_to_texture.rendertarget_render.IsValid() ||
				camera.render_to_texture.rendertarget_render.desc.width != camera.render_to_texture.resolution.x ||
				camera.render_to_texture.rendertarget_render.desc.height != camera.render_to_texture.resolution.y ||
				camera.render_to_texture.rendertarget_MSAA.desc.sample_count != camera.render_to_texture.sample_count
				)
			{
				TextureDesc desc;
				desc.width = camera.render_to_texture.resolution.x;
				desc.height = camera.render_to_texture.resolution.y;
				desc.format = wi::renderer::format_rendertarget_main;
				desc.bind_flags = BindFlag::RENDER_TARGET | BindFlag::SHADER_RESOURCE | BindFlag::UNORDERED_ACCESS;
				desc.mip_levels = GetMipCount(desc.width, desc.height);
				bool success = device->CreateTexture(&desc, nullptr, &camera.render_to_texture.rendertarget_render);
				assert(success);
				device->SetName(&camera.render_to_texture.rendertarget_render, "CameraComponent::rendertarget_render");
				success = device->CreateTexture(&desc, nullptr, &camera.render_to_texture.rendertarget_display);
				assert(success);
				device->SetName(&camera.render_to_texture.rendertarget_display, "CameraComponent::rendertarget_display");

				for (uint32_t i = 0; i < camera.render_to_texture.rendertarget_render.desc.mip_levels; ++i)
				{
					int subresource_index;
					subresource_index = device->CreateSubresource(&camera.render_to_texture.rendertarget_render, SubresourceType::SRV, 0, 1, i, 1);
					assert(subresource_index == i);
					subresource_index = device->CreateSubresource(&camera.render_to_texture.rendertarget_display, SubresourceType::SRV, 0, 1, i, 1);
					assert(subresource_index == i);
					subresource_index = device->CreateSubresource(&camera.render_to_texture.rendertarget_render, SubresourceType::UAV, 0, 1, i, 1);
					assert(subresource_index == i);
					subresource_index = device->CreateSubresource(&camera.render_to_texture.rendertarget_display, SubresourceType::UAV, 0, 1, i, 1);
					assert(subresource_index == i);
				}

				desc.mip_levels = 1;
				if (camera.render_to_texture.sample_count > 1)
				{
					desc.sample_count = camera.render_to_texture.sample_count;
					desc.layout = ResourceState::RENDERTARGET;
					desc.bind_flags = BindFlag::RENDER_TARGET;
					success = device->CreateTexture(&desc, nullptr, &camera.render_to_texture.rendertarget_MSAA);
					assert(success);
					device->SetName(&camera.render_to_texture.rendertarget_MSAA, "CameraComponent::rendertarget_MSAA");
				}
				else
				{
					camera.render_to_texture.rendertarget_MSAA = {};
				}

				desc.format = wi::renderer::format_depthbuffer_main;
				desc.bind_flags = BindFlag::DEPTH_STENCIL | BindFlag::SHADER_RESOURCE;
				desc.layout = ResourceState::SHADER_RESOURCE;
				success = device->CreateTexture(&desc, nullptr, &camera.render_to_texture.depthstencil);
				assert(success);
				device->SetName(&camera.render_to_texture.depthstencil, "CameraComponent::depthstencil");

				if (camera.render_to_texture.sample_count > 1)
				{
					desc.sample_count = 1;
					desc.bind_flags = BindFlag::SHADER_RESOURCE | BindFlag::UNORDERED_ACCESS;
					desc.layout = ResourceState::SHADER_RESOURCE;
					desc.format = Format::R32_FLOAT;
					success = device->CreateTexture(&desc, nullptr, &camera.render_to_texture.depthstencil_resolved);
					assert(success);
					device->SetName(&camera.render_to_texture.depthstencil_resolved, "CameraComponent::depthstencil_resolved");
				}
				else
				{
					camera.render_to_texture.depthstencil_resolved = {};
				}

				wi::renderer::TiledLightResources tiledres;
				wi::renderer::CreateTiledLightResources(tiledres, camera.render_to_texture.resolution);
				camera.render_to_texture.tileCount = tiledres.tileCount;
				camera.render_to_texture.entityTiles = tiledres.entityTiles;

				camera.render_to_texture.visibility = wi::allocator::make_shared<wi::renderer::Visibility>();
			}
			if (getSceneUpdateEnabled())
			{
				std::swap(camera.render_to_texture.rendertarget_render, camera.render_to_texture.rendertarget_display);
			}
			camera.width = (float)camera.render_to_texture.resolution.x;
			camera.height = (float)camera.render_to_texture.resolution.y;
			camera.UpdateCamera();

			if (camera.render_to_texture.depthstencil_resolved.IsValid())
			{
				camera.texture_depth_index = device->GetDescriptorIndex(&camera.render_to_texture.depthstencil_resolved, SubresourceType::SRV);
			}
			else
			{
				camera.texture_depth_index = device->GetDescriptorIndex(&camera.render_to_texture.depthstencil, SubresourceType::SRV);
			}
			camera.buffer_entitytiles_index = device->GetDescriptorIndex(&camera.render_to_texture.entityTiles, SubresourceType::SRV);


			wi::jobsystem::Execute(ctx, [this, cmd, i](wi::jobsystem::JobArgs args) {
				GraphicsDevice* device = GetDevice();
				wi::scene::CameraComponent& camera = scene->cameras[i]; // reload, not captured in lambda (alloc)
				wi::renderer::Visibility& visibility = *(wi::renderer::Visibility*)camera.render_to_texture.visibility.get();
				visibility.layerMask = getLayerMask();
				visibility.scene = scene;
				visibility.camera = &camera;
				visibility.flags = wi::renderer::Visibility::ALLOW_OBJECTS;
				visibility.flags |= wi::renderer::Visibility::ALLOW_HAIRS;
				visibility.flags |= wi::renderer::Visibility::ALLOW_LIGHTS;
				visibility.flags |= wi::renderer::Visibility::ALLOW_DECALS;
				visibility.flags |= wi::renderer::Visibility::ALLOW_ENVPROBES;
				wi::renderer::UpdateVisibility(visibility);

				ScopedGPUProfiling("Camera Entity", cmd);
				device->EventBegin("Camera Entity", cmd);
				wi::renderer::BindCommonResources(cmd);
				wi::renderer::BindCameraCB(
					camera,
					camera,
					camera,
					cmd
				);
				wi::graphics::Rect scissor;
				scissor.right = (int32_t)camera.render_to_texture.depthstencil.desc.width;
				scissor.bottom = (int32_t)camera.render_to_texture.depthstencil.desc.height;
				device->BindScissorRects(1, &scissor, cmd);
				Viewport vp;
				vp.width = (float)camera.render_to_texture.depthstencil.desc.width;
				vp.height = (float)camera.render_to_texture.depthstencil.desc.height;
				device->BindViewports(1, &vp, cmd);
				// prepass:
				{
					RenderPassImage rp[] = {
						RenderPassImage::DepthStencil(&camera.render_to_texture.depthstencil, RenderPassImage::LoadOp::CLEAR, RenderPassImage::StoreOp::STORE, camera.render_to_texture.depthstencil.desc.layout, ResourceState::DEPTHSTENCIL, ResourceState::SHADER_RESOURCE),
					};
					device->RenderPassBegin(rp, arraysize(rp), cmd);
					wi::renderer::DrawScene(
						visibility,
						RENDERPASS_PREPASS_DEPTHONLY,
						cmd,
						wi::renderer::DRAWSCENE_OPAQUE |
						wi::renderer::DRAWSCENE_IMPOSTOR |
						wi::renderer::DRAWSCENE_HAIRPARTICLE
					);
					device->RenderPassEnd(cmd);
				}
				if (camera.render_to_texture.depthstencil_resolved.IsValid())
				{
					wi::renderer::ResolveMSAADepthBuffer(camera.render_to_texture.depthstencil_resolved, camera.render_to_texture.depthstencil, cmd);
				}
				wi::renderer::TiledLightResources tiledres;
				tiledres.tileCount = camera.render_to_texture.tileCount;
				tiledres.entityTiles = camera.render_to_texture.entityTiles;
				wi::renderer::ComputeTiledLightCulling(tiledres, visibility, {}, cmd);
				// color pass:
				{
					RenderPassImage rp[3];
					uint32_t rp_count = 0;
					if (camera.render_to_texture.rendertarget_MSAA.IsValid())
					{
						rp[rp_count++] = RenderPassImage::RenderTarget(&camera.render_to_texture.rendertarget_MSAA, RenderPassImage::LoadOp::CLEAR, RenderPassImage::StoreOp::DONTCARE, ResourceState::RENDERTARGET, ResourceState::RENDERTARGET);
						rp[rp_count++] = RenderPassImage::Resolve(&camera.render_to_texture.rendertarget_render, ResourceState::SHADER_RESOURCE, ResourceState::SHADER_RESOURCE, 0);
					}
					else
					{
						rp[rp_count++] = RenderPassImage::RenderTarget(&camera.render_to_texture.rendertarget_render, RenderPassImage::LoadOp::CLEAR);
					}
					rp[rp_count++] = RenderPassImage::DepthStencil(&camera.render_to_texture.depthstencil, RenderPassImage::LoadOp::LOAD, RenderPassImage::StoreOp::DONTCARE, ResourceState::SHADER_RESOURCE, ResourceState::DEPTHSTENCIL, camera.render_to_texture.depthstencil.desc.layout);
					device->RenderPassBegin(rp, rp_count, cmd);
					wi::renderer::DrawScene(
						visibility,
						RENDERPASS_MAIN,
						cmd,
						wi::renderer::DRAWSCENE_OPAQUE |
						wi::renderer::DRAWSCENE_IMPOSTOR |
						wi::renderer::DRAWSCENE_HAIRPARTICLE
					);
					wi::renderer::DrawScene(
						visibility,
						RENDERPASS_MAIN,
						cmd,
						wi::renderer::DRAWSCENE_TRANSPARENT
					);
					wi::renderer::DrawSky(*scene, cmd);
					wi::renderer::DrawLightVisualizers(visibility, cmd);
					device->RenderPassEnd(cmd);

					if (camera.IsCRT() && getSceneUpdateEnabled())
					{
						wi::renderer::Postprocess_CRT(camera.render_to_texture.rendertarget_render, camera.render_to_texture.rendertarget_display, cmd, 0, 0, false);
						// Swap so rendertarget_render now contains the CRT-processed result for mipchain and display
						std::swap(camera.render_to_texture.rendertarget_render, camera.render_to_texture.rendertarget_display);
					}

					wi::renderer::GenerateMipChain(camera.render_to_texture.rendertarget_render, wi::renderer::MIPGENFILTER_LINEAR, cmd);
				}
				device->EventEnd(cmd);
			});
		}
	}

	void RenderPath3D::setAO(AO value)
	{
		ao = value;

		if (!rtParticleDistortion.IsValid())
			return; // ResizeBuffers hasn't been called yet

		rtAO = {};
		ssaoResources = {};
		msaoResources = {};
		rtaoResources = {};

		if (ao == AO_DISABLED)
		{
			return;
		}

		XMUINT2 internalResolution = GetInternalResolution();
		if (internalResolution.x == 0 || internalResolution.y == 0)
			return;

		TextureDesc desc;
		desc.bind_flags = BindFlag::SHADER_RESOURCE | BindFlag::UNORDERED_ACCESS | BindFlag::RENDER_TARGET; // render target binding for aliasing (in case resource heap tier < 2)
		desc.format = Format::R8_UNORM;
		desc.layout = ResourceState::SHADER_RESOURCE_COMPUTE;

		switch (ao)
		{
		case RenderPath3D::AO_SSAO:
		case RenderPath3D::AO_HBAO:
			desc.width = internalResolution.x / 2;
			desc.height = internalResolution.y / 2;
			break;
		case RenderPath3D::AO_MSAO:
			desc.width = internalResolution.x;
			desc.height = internalResolution.y;
			break;
		case RenderPath3D::AO_RTAO:
			desc.width = internalResolution.x;
			desc.height = internalResolution.y;
			break;
		default:
			break;
		}

		if (ComputeTextureMemorySizeInBytes(desc) > ComputeTextureMemorySizeInBytes(rtParticleDistortion.desc))
		{
			// There would be resource aliasing error if we proceed like this!
			//	looks like ResizeBuffers() hasn't been called yet for the current internal resolution
			//	if this happens, then ResizeBuffers() will be called next frame probably and then AO resources
			//	will be created successdully
			return;
		}

		switch (ao)
		{
		case RenderPath3D::AO_SSAO:
		case RenderPath3D::AO_HBAO:
			wi::renderer::CreateSSAOResources(ssaoResources, internalResolution);
			break;
		case RenderPath3D::AO_MSAO:
			wi::renderer::CreateMSAOResources(msaoResources, internalResolution);
			break;
		case RenderPath3D::AO_RTAO:
			wi::renderer::CreateRTAOResources(rtaoResources, internalResolution);
			break;
		default:
			break;
		}

		GraphicsDevice* device = wi::graphics::GetDevice();
		assert(ComputeTextureMemorySizeInBytes(desc) <= ComputeTextureMemorySizeInBytes(rtParticleDistortion.desc)); // aliasing check
		device->CreateTexture(&desc, nullptr, &rtAO, &rtParticleDistortion); // aliasing!
		device->SetName(&rtAO, "renderpath3D.rtAO");
	}
	void RenderPath3D::setSSREnabled(bool value)
	{
		ssrEnabled = value;

		if (value)
		{
			GraphicsDevice* device = wi::graphics::GetDevice();
			XMUINT2 internalResolution = GetInternalResolution();
			if (internalResolution.x == 0 || internalResolution.y == 0)
				return;

			TextureDesc desc;
			desc.bind_flags = BindFlag::SHADER_RESOURCE | BindFlag::UNORDERED_ACCESS;
			desc.format = Format::R16G16B16A16_FLOAT;
			desc.width = internalResolution.x;
			desc.height = internalResolution.y;
			desc.layout = ResourceState::SHADER_RESOURCE_COMPUTE;
			device->CreateTexture(&desc, nullptr, &rtSSR);
			device->SetName(&rtSSR, "renderpath3D.rtSSR");

			wi::renderer::CreateSSRResources(ssrResources, internalResolution, ssrQuality);
		}
		else
		{
			ssrResources = {};
		}
	}
	void RenderPath3D::setSSGIEnabled(bool value)
	{
		ssgiEnabled = value;

		if (value)
		{
			GraphicsDevice* device = wi::graphics::GetDevice();
			XMUINT2 internalResolution = GetInternalResolution();
			if (internalResolution.x == 0 || internalResolution.y == 0)
				return;

			TextureDesc desc;
			desc.bind_flags = BindFlag::SHADER_RESOURCE | BindFlag::UNORDERED_ACCESS;
			desc.format = Format::R16G16B16A16_FLOAT;
			desc.width = internalResolution.x;
			desc.height = internalResolution.y;
			desc.layout = ResourceState::SHADER_RESOURCE_COMPUTE;
			device->CreateTexture(&desc, nullptr, &rtSSGI);
			device->SetName(&rtSSGI, "renderpath3D.rtSSGI");

			wi::renderer::CreateSSGIResources(ssgiResources, internalResolution);
		}
		else
		{
			ssgiResources = {};
		}
	}
	void RenderPath3D::setRaytracedReflectionsEnabled(bool value)
	{
		raytracedReflectionsEnabled = value;

		if (value)
		{
			GraphicsDevice* device = wi::graphics::GetDevice();
			XMUINT2 internalResolution = GetInternalResolution();
			if (internalResolution.x == 0 || internalResolution.y == 0)
				return;

			TextureDesc desc;
			desc.bind_flags = BindFlag::SHADER_RESOURCE | BindFlag::UNORDERED_ACCESS;
			desc.format = Format::R16G16B16A16_FLOAT;
			desc.width = internalResolution.x;
			desc.height = internalResolution.y;
			device->CreateTexture(&desc, nullptr, &rtSSR);
			device->SetName(&rtSSR, "renderpath3D.rtSSR");

			wi::renderer::CreateRTReflectionResources(rtreflectionResources, internalResolution, raytracedReflectionsQuality);
		}
		else
		{
			rtreflectionResources = {};
		}
	}
	void RenderPath3D::setRaytracedDiffuseEnabled(bool value)
	{
		raytracedDiffuseEnabled = value;

		if (value)
		{
			GraphicsDevice* device = wi::graphics::GetDevice();
			XMUINT2 internalResolution = GetInternalResolution();
			if (internalResolution.x == 0 || internalResolution.y == 0)
				return;

			TextureDesc desc;
			desc.bind_flags = BindFlag::SHADER_RESOURCE | BindFlag::UNORDERED_ACCESS;
			desc.format = Format::R16G16B16A16_FLOAT;
			desc.width = internalResolution.x;
			desc.height = internalResolution.y;
			device->CreateTexture(&desc, nullptr, &rtRaytracedDiffuse);
			device->SetName(&rtRaytracedDiffuse, "renderpath3D.rtRaytracedDiffuse");

			wi::renderer::CreateRTDiffuseResources(rtdiffuseResources, internalResolution, raytracedDiffuseQuality);
		}
		else
		{
			rtRaytracedDiffuse = {};
			rtdiffuseResources = {};
		}
	}
	void RenderPath3D::setFSREnabled(bool value)
	{
		fsrEnabled = value;

		if (resolutionScale < 1.0f && fsrEnabled)
		{
			GraphicsDevice* device = wi::graphics::GetDevice();
			if (GetPhysicalWidth() == 0 || GetPhysicalHeight() == 0)
				return;

			TextureDesc desc;
			desc.bind_flags = BindFlag::SHADER_RESOURCE | BindFlag::UNORDERED_ACCESS;
			desc.format = wi::renderer::format_rendertarget_main;
			desc.width = GetPhysicalWidth();
			desc.height = GetPhysicalHeight();
			device->CreateTexture(&desc, nullptr, &rtFSR[0]);
			device->SetName(&rtFSR[0], "renderpath3D.rtFSR[0]");
			device->CreateTexture(&desc, nullptr, &rtFSR[1]);
			device->SetName(&rtFSR[1], "renderpath3D.rtFSR[1]");
		}
		else
		{
			if (!getFSR2Enabled())
			{
				// These are used both for FSR and FSR2
				rtFSR[0] = {};
				rtFSR[1] = {};
			}
		}
	}
	void RenderPath3D::setFSR2Enabled(bool value)
	{
		fsr2Enabled = value;

		if (fsr2Enabled)
		{
			GraphicsDevice* device = wi::graphics::GetDevice();
			if (GetPhysicalWidth() == 0 || GetPhysicalHeight() == 0)
				return;

			XMUINT2 displayResolution = XMUINT2(
				std::max(GetPhysicalWidth(), GetInternalResolution().x),
				std::max(GetPhysicalHeight(), GetInternalResolution().y)
			);

			wi::renderer::CreateFSR2Resources(fsr2Resources, GetInternalResolution(), displayResolution);

			TextureDesc desc;
			desc.bind_flags = BindFlag::SHADER_RESOURCE | BindFlag::UNORDERED_ACCESS;
			desc.format = wi::renderer::format_rendertarget_main;
			desc.width = displayResolution.x;
			desc.height = displayResolution.y;
			device->CreateTexture(&desc, nullptr, &rtFSR[0]);
			device->SetName(&rtFSR[0], "renderpath3D.rtFSR[0]");
			device->CreateTexture(&desc, nullptr, &rtFSR[1]);
			device->SetName(&rtFSR[1], "renderpath3D.rtFSR[1]");
		}
		else
		{
			fsr2Resources = {};
			if (!getFSREnabled())
			{
				// These are used both for FSR and FSR2
				rtFSR[0] = {};
				rtFSR[1] = {};
			}
		}

		// Depending on FSR2 is on/off, these either need to run at display or internal resolution:
		motionblurResources = {};
		depthoffieldResources = {};
	}
	void RenderPath3D::setFSR2Preset(FSR2_Preset preset)
	{
		wi::graphics::SamplerDesc desc = wi::renderer::GetSampler(wi::enums::SAMPLER_OBJECTSHADER)->GetDesc();
		switch (preset)
		{
		default:
		case FSR2_Preset::Quality:
			resolutionScale = 1.0f / 1.5f;
			desc.mip_lod_bias = -1.58f;
			break;
		case FSR2_Preset::Balanced:
			resolutionScale = 1.0f / 1.7f;
			desc.mip_lod_bias = -1.76f;
			break;
		case FSR2_Preset::Performance:
			resolutionScale = 1.0f / 2.0f;
			desc.mip_lod_bias = -2.0f;
			break;
		case FSR2_Preset::Ultra_Performance:
			resolutionScale = 1.0f / 3.0f;
			desc.mip_lod_bias = -2.58f;
			break;
		}
		wi::renderer::ModifyObjectSampler(desc);
	}
	void RenderPath3D::setMotionBlurEnabled(bool value)
	{
		motionBlurEnabled = value;
	}
	void RenderPath3D::setDepthOfFieldEnabled(bool value)
	{
		depthOfFieldEnabled = value;
	}
	void RenderPath3D::setEyeAdaptionEnabled(bool value)
	{
		eyeAdaptionEnabled = value;

		if (value)
		{
			wi::renderer::CreateLuminanceResources(luminanceResources, GetInternalResolution());
		}
		else
		{
			luminanceResources = {};
		}
	}
	void RenderPath3D::setReflectionsEnabled(bool value)
	{
		reflectionsEnabled = value;

		if (value)
		{
			GraphicsDevice* device = wi::graphics::GetDevice();
			XMUINT2 internalResolution = GetInternalResolution();
			if (internalResolution.x == 0 || internalResolution.y == 0)
				return;

			TextureDesc desc;
			desc.format = wi::renderer::format_rendertarget_main;
			desc.width = uint32_t((float)internalResolution.x * planarReflectionResolutionScale);
			desc.height = uint32_t((float)internalResolution.y * planarReflectionResolutionScale);
			desc.sample_count = planarReflectionMSAASampleCount;
			if (desc.sample_count > 1)
			{
				// Deliberately not TRANSIENT_ATTACHMENT, though nothing ever
				// samples this texture: the reflection is drawn in two render
				// passes so that its transparents have an opaque copy to
				// refract, and the second loads what the first stored. Transient
				// memory does not survive a pass boundary. rtMain_render is not
				// transient for the same reason.
				desc.bind_flags = BindFlag::RENDER_TARGET;
				desc.layout = ResourceState::RENDERTARGET;
			}
			else
			{
				desc.bind_flags = BindFlag::RENDER_TARGET | BindFlag::SHADER_RESOURCE;
			}
			device->CreateTexture(&desc, nullptr, &rtReflection_render);
			device->SetName(&rtReflection_render, "renderpath3D.rtReflection_render");

			desc.misc_flags = ResourceMiscFlag::NONE;
			desc.bind_flags = BindFlag::DEPTH_STENCIL | BindFlag::SHADER_RESOURCE;
			desc.format = wi::renderer::format_depthbuffer_main;
			desc.layout = ResourceState::SHADER_RESOURCE;
			device->CreateTexture(&desc, nullptr, &depthBuffer_Reflection_render);
			device->SetName(&depthBuffer_Reflection_render, "renderpath3D.depthBuffer_Reflection_render");

			if (desc.sample_count > 1)
			{
				desc.sample_count = 1;
				desc.format = wi::renderer::format_rendertarget_main;
				desc.bind_flags = BindFlag::RENDER_TARGET | BindFlag::SHADER_RESOURCE;
				device->CreateTexture(&desc, nullptr, &rtReflection);
				device->SetName(&rtReflection, "renderpath3D.rtReflection");

				desc.format = Format::R16_UNORM;
				desc.bind_flags = BindFlag::UNORDERED_ACCESS | BindFlag::SHADER_RESOURCE;
				device->CreateTexture(&desc, nullptr, &depthBuffer_Reflection);
				device->SetName(&depthBuffer_Reflection, "renderpath3D.depthBuffer_Reflection");
			}
			else
			{
				rtReflection = rtReflection_render;
				depthBuffer_Reflection = depthBuffer_Reflection_render;
			}

			{
				// Quartered against the reflection rather than against the
				// screen, so the refraction blurs by the same fraction of the
				// image it does in the main view however the reflection is
				// scaled.
				TextureDesc desc;
				desc.bind_flags = BindFlag::SHADER_RESOURCE | BindFlag::UNORDERED_ACCESS;
				desc.format = wi::renderer::format_rendertarget_main;
				desc.width = std::max(1u, rtReflection.desc.width / 4);
				desc.height = std::max(1u, rtReflection.desc.height / 4);
				desc.mip_levels = std::max(1u, std::min(8u, (uint32_t)std::log2(std::max(desc.width, desc.height))));
				device->CreateTextureZeroed(&desc, &rtSceneCopy_Reflection);
				device->SetName(&rtSceneCopy_Reflection, "renderpath3D.rtSceneCopy_Reflection");
				device->CreateMipgenSubresources(rtSceneCopy_Reflection);
			}

			wi::renderer::CreateTiledLightResources(tiledLightResources_planarReflection, XMUINT2(depthBuffer_Reflection.desc.width, depthBuffer_Reflection.desc.height));
		}
		else
		{
			rtReflection_render = {};
			rtReflection = {};
			rtSceneCopy_Reflection = {};
			depthBuffer_Reflection_render = {};
			depthBuffer_Reflection = {};
			tiledLightResources_planarReflection = {};
		}
	}
	void RenderPath3D::setPlanarReflectionQuality(float resolutionScale, uint32_t msaaSampleCount)
	{
		planarReflectionResolutionScale = resolutionScale;
		planarReflectionMSAASampleCount = msaaSampleCount;
		setReflectionsEnabled(getReflectionsEnabled());
	}
	void RenderPath3D::setBloomEnabled(bool value)
	{
		bloomEnabled = value;

		if (value)
		{
			const XMUINT2 internalResolution = GetInternalResolution();
			if (internalResolution.x == 0 || internalResolution.y == 0)
				return;
			wi::renderer::CreateBloomResources(bloomResources, internalResolution);
		}
		else
		{
			bloomResources = {};
		}
	}
	void RenderPath3D::setVolumeLightsEnabled(bool value)
	{
		volumeLightsEnabled = value;
	}
	void RenderPath3D::setLightShaftsEnabled(bool value)
	{
		lightShaftsEnabled = value;

		if (value)
		{
			GraphicsDevice* device = wi::graphics::GetDevice();
			XMUINT2 internalResolution = GetInternalResolution();
			if (internalResolution.x == 0 || internalResolution.y == 0)
				return;

			TextureDesc desc;
			desc.bind_flags = BindFlag::RENDER_TARGET | BindFlag::SHADER_RESOURCE;
			desc.format = wi::renderer::format_rendertarget_main;
			desc.width = internalResolution.x;
			desc.height = internalResolution.y;
			desc.sample_count = getMSAASampleCount();
			device->CreateTexture(&desc, nullptr, &rtSun[0]);
			device->SetName(&rtSun[0], "renderpath3D.rtSun[0]");

			desc.bind_flags = BindFlag::SHADER_RESOURCE | BindFlag::UNORDERED_ACCESS;
			desc.sample_count = 1;
			desc.width = internalResolution.x / 4;
			desc.height = internalResolution.y / 4;
			device->CreateTexture(&desc, nullptr, &rtSun[1]);
			device->SetName(&rtSun[1], "renderpath3D.rtSun[1]");
			device->CreateTexture(&desc, nullptr, &rtSun[2]);
			device->SetName(&rtSun[2], "renderpath3D.rtSun[2]");

			if (getMSAASampleCount() > 1)
			{
				desc.bind_flags = BindFlag::RENDER_TARGET | BindFlag::SHADER_RESOURCE;
				desc.width = internalResolution.x;
				desc.height = internalResolution.y;
				desc.sample_count = 1;
				device->CreateTexture(&desc, nullptr, &rtSun_resolved);
				device->SetName(&rtSun_resolved, "renderpath3D.rtSun_resolved");
			}
		}
		else
		{
			rtSun[0] = {};
			rtSun[1] = {};
			rtSun[2] = {};
			rtSun_resolved = {};
		}
	}
	void RenderPath3D::setOutlineEnabled(bool value)
	{
		outlineEnabled = value;

		if (value)
		{
			GraphicsDevice* device = wi::graphics::GetDevice();
			XMUINT2 internalResolution = GetInternalResolution();
			if (internalResolution.x == 0 || internalResolution.y == 0)
				return;

			TextureDesc desc;
			desc.bind_flags = BindFlag::RENDER_TARGET | BindFlag::SHADER_RESOURCE;
			desc.format = Format::R32_FLOAT;
			desc.width = internalResolution.x;
			desc.height = internalResolution.y;
			device->CreateTexture(&desc, nullptr, &rtOutlineSource);
			device->SetName(&rtOutlineSource, "renderpath3D.rtOutlineSource");

			if (getMSAASampleCount() > 1)
			{
				desc.bind_flags = BindFlag::RENDER_TARGET;
				desc.sample_count = getMSAASampleCount();
				desc.misc_flags = ResourceMiscFlag::TRANSIENT_ATTACHMENT;
				desc.layout = ResourceState::RENDERTARGET;
				device->CreateTexture(&desc, nullptr, &rtOutlineSource_MSAA);
				device->SetName(&rtOutlineSource_MSAA, "renderpath3D.rtOutlineSource_MSAA");
			}
		}
		else
		{
			rtOutlineSource = {};
			rtOutlineSource_MSAA = {};
		}
	}

	Texture RenderPath3D::CreateScreenshotWithAlphaBackground(uint8_t stencilref, wi::image::STENCILMODE stencilmode, wi::image::STENCILREFMODE stencilrefmode)
	{
		TextureDesc desc = rtMain_render.GetDesc();
		desc.format = Format::R8G8B8A8_UNORM;
		desc.bind_flags = BindFlag::RENDER_TARGET | BindFlag::SHADER_RESOURCE;
		Texture tex;
		GraphicsDevice* device = GetDevice();
		bool success = device->CreateTexture(&desc, nullptr, &tex);
		assert(success);

		Texture tex_resolved;
		if (desc.sample_count > 1)
		{
			desc.sample_count = 1;
			success = device->CreateTexture(&desc, nullptr, &tex_resolved);
			assert(success);
		}

		CommandList cmd = device->BeginCommandList();
		RenderPassImage rp[] = {
			RenderPassImage::RenderTarget(&tex, RenderPassImage::LoadOp::CLEAR),
			RenderPassImage::DepthStencil(GetDepthStencil()),
			RenderPassImage::Resolve(&tex_resolved),
		};
		device->RenderPassBegin(rp, tex_resolved.IsValid() ? 3 : 2, cmd);
		Viewport vp;
		vp.width = (float)desc.width;
		vp.height = (float)desc.height;
		device->BindViewports(1, &vp, cmd);
		wi::image::Params fx;
		fx.stencilComp = stencilmode;
		fx.stencilRef = stencilref;
		fx.stencilRefMode = stencilrefmode;
		fx.enableFullScreen();
		wi::image::Draw(GetLastPostprocessRT(), fx, cmd);
		device->RenderPassEnd(cmd);

		if (tex_resolved.IsValid())
			return tex_resolved;
		return tex;
	}

}
