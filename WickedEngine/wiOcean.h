#pragma once
#include "CommonInclude.h"
#include "wiGraphicsDevice.h"
#include "wiFFTGenerator.h"
#include "wiScene_Decl.h"
#include "wiMath.h"
#include "wiPrimitive.h"
#include "wiScene/environment/WaterMedium/WaterMedium.h"

namespace wi
{
	class Ocean
	{
	public:
		struct OceanParameters
		{
			// Must be power of 2.
			int dmap_dim = 512;
			// Typical value is 1000 ~ 2000
			float patch_length = 50.0f;

			// Adjust the time interval for simulation.
			float time_scale = 0.3f;
			// Amplitude for transverse wave. Around 1.0
			float wave_amplitude = 1000.0f;
			// Wind direction. Normalization not required.
			XMFLOAT2 wind_dir = XMFLOAT2(0.8f, 0.6f);
			// Around 100 ~ 1000
			float wind_speed = 600.0f;
			// This value damps out the waves against the wind direction.
			// Smaller value means higher wind dependency.
			float wind_dependency = 0.07f;
			// The amplitude for longitudinal wave. Must be positive.
			float choppy_scale = 1.3f;


			// Deprecated: nothing reads this. It tinted the foam blend before
			// foamColor existed, and its alpha was overwritten before use. Kept
			// because the serializer format is positional, and because the
			// alpha still carries the pre-medium water density that scenes
			// older than version 8 are migrated from. Also still bound to Lua.
			XMFLOAT4 waterColor = XMFLOAT4(0.0f, 2.0f / 255.0f, 6.0f / 255.0f, 0.6f);

			// Foam colour in rgb, and how opaque the foam is in w. The default
			// is the flat grey the foam used to be fixed at, so an untouched
			// scene looks as it did.
			XMFLOAT4 foamColor = XMFLOAT4(0.6f, 0.6f, 0.6f, 1.0f);

			// Deprecated: nothing reads this. It tinted the wave transmission
			// glow before that was derived from waterMedium below, which now
			// supplies the colour, the phase and both Fresnel interfaces.
			// Kept only so the serializer keeps its slot - the archive format
			// is positional, so dropping the field would shift every ocean
			// parameter written after it and misread every existing scene.
			XMFLOAT4 extinctionColor = XMFLOAT4(0, 0.9f, 1, 1);

			// Inherent optical properties (absorption and scattering) of the
			// water, driving turbidity and how far one can see through it.
			wi::scene::environment::WaterMedium waterMedium;

			float waterHeight = 0.0f;

			/**
			 * Resolution of the surface mesh near the viewer.
			 *
			 * Each step halves the innermost clipmap cell, so the surface can
			 * carry finer waves close in, and adds a level so that the ocean
			 * still reaches exactly as far. Detail therefore costs one ring of
			 * triangles per step and does not change coverage.
			 */
			uint32_t surfaceDetail = 4;

			/**
			 * **Unused. Retained only so existing scenes still load.**
			 *
			 * This extruded the old screen space projected grid to cope with
			 * vertices that displacement pushed outside the view. The surface
			 * is a camera-centred clipmap now, built in world space, so there
			 * is no screen rectangle to overshoot and nothing for this to do.
			 * It stays in the struct, the serializer and the Lua binding
			 * because removing it would break the scene format; the editor no
			 * longer shows it.
			 */
			float surfaceDisplacementTolerance = 2;
		};
		void Create(const OceanParameters& params);

		void UpdateDisplacementMap(wi::graphics::CommandList cmd) const;

		/**
		 * Draws the ocean into an environment probe's cubemap.
		 *
		 * @param[in] viewerPosition - World position the clipmap centres on,
		 *                             i.e. the probe being captured.
		 * @param[in] cmd - Command list to record into.
		 */
		void RenderForCubemap(const XMFLOAT3& viewerPosition, wi::graphics::CommandList cmd) const;

		/**
		 * Draws the ocean into a shadow map.
		 *
		 * @param[in] viewerPosition - World position the clipmap centres on.
		 *                             This is the **main camera**, not the
		 *                             light: the shadow map covers the stretch
		 *                             of sea the viewer can see, so the mesh
		 *                             has to be the one the viewer gets.
		 * @param[in] cmd - Command list to record into.
		 */
		void RenderForShadowmap(const XMFLOAT3& viewerPosition, wi::graphics::CommandList cmd) const;
		void Render(const wi::scene::CameraComponent& camera, wi::graphics::CommandList cmd) const;

		void CopyDisplacementMapReadback(wi::graphics::CommandList cmd) const;

		const wi::graphics::Texture* getDisplacementMap() const;
		const wi::graphics::Texture* getGradientMap() const;

		const wi::primitive::AABB GetAABB(const XMFLOAT3& camera_pos) const;

		/**
		 * Distance band over which the wave displacement flattens, in metres.
		 *
		 * Derived from the clipmap's geometry rather than authored. A level's
		 * cell size and its extent are tied together, so the cells at distance
		 * \f$d\f$ are about \f$2d/N\f$ across for \f$N\f$ cells per side.
		 * Inverting that gives the distance at which the cells grow to a
		 * chosen size, and the size that matters is the one measured against
		 * the wave patch: a patch spanned by plenty of cells is drawn
		 * faithfully, a patch spanned by two is at the Nyquist limit and
		 * cannot be drawn at all.
		 *
		 * The band therefore moves automatically with the patch size and the
		 * mesh resolution, which a pair of authored metres does not - that was
		 * the failure the previous 16-160 m band was heading for.
		 *
		 * @return `x` where flattening begins, `y` where it is complete.
		 *
		 * @note Both the surface vertex shader and everything that has to
		 *       agree with the surface as drawn read this through
		 *       `ShaderOcean::displacement_fade`; it is published once so the
		 *       two cannot drift apart.
		 */
		[[nodiscard]] XMFLOAT2 GetDisplacementFadeBand() const noexcept;

		static void Initialize();

		bool IsValid() const { return displacementMap.IsValid(); }

		// occlusion result history bitfield (32 bit->32 frame history)
		mutable uint32_t occlusionHistory = ~0u;
		mutable int occlusionQueries[wi::graphics::GraphicsDevice::GetBufferCount()];
		inline bool IsOccluded() const
		{
			// Perform a conservative occlusion test:
			// If it is visible in any frames in the history, it is determined visible in this frame
			// But if all queries failed in the history, it is occluded.
			// If it pops up for a frame after occluded, it is visible again for some frames
			return occlusionHistory == 0;
		}

		// Return the position at world space modified by the ocean displacement map
		XMFLOAT3 GetDisplacedPosition(const XMFLOAT3& worldPosition) const;

		OceanParameters params;

	protected:
		wi::graphics::Texture displacementMap;		// (RGBA32F)
		wi::graphics::Texture gradientMap;			// (RGBA16F)

		wi::graphics::Texture displacementMap_readback[wi::graphics::GraphicsDevice::GetBufferCount()];		// (RGBA32F)
		mutable bool displacement_readback_valid[arraysize(displacementMap_readback)] = {};
		mutable uint32_t displacement_readback_index = 0;
		// Number of frames the displacement map readback is still kept alive for.
		// Each GetDisplacedPosition() query refreshes this; CopyDisplacementMapReadback()
		// decrements it and skips the (expensive) GPU->CPU copy once it reaches zero, so
		// the full-texture readback is only paid for while something queries ocean height.
		mutable uint32_t displacement_readback_request = 0;

		void initHeightMap(XMFLOAT2* out_h0, float* out_omega);


		// Initial height field H(0) generated by Phillips spectrum & Gauss distribution.
		wi::graphics::GPUBuffer buffer_Float2_H0;

		// Angular frequency
		wi::graphics::GPUBuffer buffer_Float_Omega;

		// Height field H(t), choppy field Dx(t) and Dy(t) in frequency domain, updated each frame.
		wi::graphics::GPUBuffer buffer_Float2_Ht;

		// Height & choppy buffer in the space domain, corresponding to H(t), Dx(t) and Dy(t)
		wi::graphics::GPUBuffer buffer_Float_Dxyz;

	};
}
