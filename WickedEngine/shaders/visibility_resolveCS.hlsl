#include "globals.hlsli"
#include "ShaderInterop_Renderer.h"
#include "surfaceHF.hlsli"
#include "raytracingHF.hlsli"

#ifdef VISIBILITY_MSAA
Texture2DMS<uint> input_primitiveID : register(t0);
#else
Texture2D<uint> input_primitiveID : register(t0);
#endif // VISIBILITY_MSAA

groupshared uint local_bin_mask;
groupshared uint local_bin_execution_mask_0[SHADERTYPE_BIN_COUNT + 1];
groupshared uint local_bin_execution_mask_1[SHADERTYPE_BIN_COUNT + 1];

RWStructuredBuffer<ShaderTypeBin> output_bins : register(u0);
RWStructuredBuffer<VisibilityTile> output_binned_tiles : register(u1);

RWTexture2D<float> output_depth_mip0 : register(u3);
RWTexture2D<float> output_depth_mip1 : register(u4);
RWTexture2D<float> output_depth_mip2 : register(u5);
RWTexture2D<float> output_depth_mip3 : register(u6);
RWTexture2D<float> output_depth_mip4 : register(u7);

RWTexture2D<float> output_lineardepth_mip0 : register(u8);
RWTexture2D<float> output_lineardepth_mip1 : register(u9);
RWTexture2D<float> output_lineardepth_mip2 : register(u10);
RWTexture2D<float> output_lineardepth_mip3 : register(u11);
RWTexture2D<float> output_lineardepth_mip4 : register(u12);

#ifdef VISIBILITY_MSAA
RWTexture2D<uint> output_primitiveID : register(u13);
#endif // VISIBILITY_MSAA

[numthreads(VISIBILITY_BLOCKSIZE, VISIBILITY_BLOCKSIZE, 1)]
void main(uint2 Gid : SV_GroupID, uint groupIndex : SV_GroupIndex)
{
	if (groupIndex <= SHADERTYPE_BIN_COUNT)
	{
		if (groupIndex == 0)
		{
			local_bin_mask = 0;
		}
		local_bin_execution_mask_0[groupIndex] = 0;
		local_bin_execution_mask_1[groupIndex] = 0;
	}
	GroupMemoryBarrierWithGroupSync();

	const uint2 GTid = remap_lane_8x8(groupIndex);
	const uint2 pixel = Gid.xy * VISIBILITY_BLOCKSIZE + GTid.xy;
	const bool pixel_valid = (pixel.x < GetCamera().internal_resolution.x) && (pixel.y < GetCamera().internal_resolution.y);
	
	RayDesc ray = CreateCameraRay(pixel);
	
#ifdef VISIBILITY_MSAA
	uint primitiveID = input_primitiveID.Load(pixel, 0);
#else
	uint primitiveID = input_primitiveID[pixel];
#endif // VISIBILITY_MSAA

#ifdef VISIBILITY_MSAA
	output_primitiveID[pixel] = primitiveID;
#endif // VISIBILITY_MSAA

	float depth = 1; // invalid
	uint bin = SHADERTYPE_BIN_COUNT; // default to sky bin
	if (pixel_valid)
	{
		[branch]
		if (any(primitiveID))
		{
			PrimitiveID prim;
			prim.init();
			prim.unpack(primitiveID);

			Surface surface;
			surface.init();

			bool surface_loaded = surface.load(prim, ray.Origin, ray.Direction);
			if (!surface_loaded)
			{
				float fallback_depth = texture_depth[pixel];
				if (fallback_depth < 1.0f)
				{
					float fallback_lineardepth = compute_lineardepth(fallback_depth) * GetCamera().z_far_rcp;
					float depth_linear = fallback_lineardepth * GetCamera().z_range + GetCamera().z_near;
					float4 fallback_svposition = float4(float2(pixel) + 0.5f, fallback_depth, depth_linear);
					float3 fallback_positionWS = GetCamera().screen_to_world(fallback_svposition);
					surface_loaded = surface.load(prim, fallback_positionWS);
					if (surface_loaded)
					{
						depth = fallback_depth;
					}
				}
			}

			if (surface_loaded)
			{
				float4 tmp = mul(GetCamera().view_projection, float4(surface.P, 1));
				tmp.xyz /= max(0.0001, tmp.w); // max: avoid nan
				depth = saturate(tmp.z); // saturate: avoid blown up values

				bin = surface.material.GetShaderType();
			}
			else
			{
				primitiveID = 0; // mark as sky if surface data could not be resolved
				depth = 0;
			}
		}
		else
		{
			// sky:
			depth = 0;
		}
	}
	
	// Always mark this thread as active in the execution mask, even if pixel is outside viewport.
	// Consumer shaders need to know which threads exist in the group, not just which pixels are valid.
	// Invalid pixels will be skipped by pixel_valid checks in the consumer, but the execution mask
	// must represent all participating threads to avoid incorrect early returns.
	if (bin <= SHADERTYPE_BIN_COUNT)
	{
		if (groupIndex < 32)
		{
			InterlockedOr(local_bin_execution_mask_0[bin], 1u << groupIndex);
		}
		else
		{
			InterlockedOr(local_bin_execution_mask_1[bin], 1u << (groupIndex - 32u));
		}
	}

	// Accumulate which shader-type bins are used within this threadgroup.
	// Avoid Wave intrinsics for portability (Metal subgroup behavior can be unreliable),
	// so perform aggregation on the group leader after all lanes have updated the local execution masks.
	GroupMemoryBarrierWithGroupSync();
	if (groupIndex == 0)
	{
		uint tmp_mask = 0;
		for (uint b = 0; b < SHADERTYPE_BIN_COUNT + 1; ++b)
		{
			if (local_bin_execution_mask_0[b] != 0 || local_bin_execution_mask_1[b] != 0)
			{
				tmp_mask |= 1u << b;
			}
		}
		local_bin_mask = tmp_mask;
	}

	GroupMemoryBarrierWithGroupSync();
	if (groupIndex < SHADERTYPE_BIN_COUNT + 1)
	{
		if (local_bin_mask & (1u << groupIndex))
		{
			uint bin_tile_list_offset = groupIndex * GetCamera().visibility_tilecount_flat;
			uint tile_offset = 0;
			InterlockedAdd(output_bins[groupIndex].dispatchX, 1, tile_offset);

			VisibilityTile tile;
			tile.visibility_tile_id = pack_pixel(Gid.xy);
			tile.entity_flat_tile_index = flatten2D(Gid.xy / VISIBILITY_TILED_CULLING_GRANULARITY, GetCamera().entity_culling_tilecount.xy) * SHADER_ENTITY_TILE_BUCKET_COUNT;
			tile.execution_mask = uint64_t(local_bin_execution_mask_0[groupIndex]) | (uint64_t(local_bin_execution_mask_1[groupIndex]) << uint64_t(32));
			output_binned_tiles[bin_tile_list_offset + tile_offset] = tile;
		}
	}

	// Downsample depths (only write if pixel is valid to avoid out-of-bounds writes):
	if (pixel_valid)
	{
		output_depth_mip0[pixel] = depth;
		if (GTid.x % 2 == 0 && GTid.y % 2 == 0)
		{
			output_depth_mip1[pixel / 2] = depth;
		}
		if (GTid.x % 4 == 0 && GTid.y % 4 == 0)
		{
			output_depth_mip2[pixel / 4] = depth;
		}
		if (GTid.x % 8 == 0 && GTid.y % 8 == 0)
		{
			output_depth_mip3[pixel / 8] = depth;
		}
		if (GTid.x % 16 == 0 && GTid.y % 16 == 0)
		{
			output_depth_mip4[pixel / 16] = depth;
		}

		float lineardepth = compute_lineardepth(depth) * GetCamera().z_far_rcp;
		output_lineardepth_mip0[pixel] = lineardepth;
		if (GTid.x % 2 == 0 && GTid.y % 2 == 0)
		{
			output_lineardepth_mip1[pixel / 2] = lineardepth;
		}
		if (GTid.x % 4 == 0 && GTid.y % 4 == 0)
		{
			output_lineardepth_mip2[pixel / 4] = lineardepth;
		}
		if (GTid.x % 8 == 0 && GTid.y % 8 == 0)
		{
			output_lineardepth_mip3[pixel / 8] = lineardepth;
		}
		if (GTid.x % 16 == 0 && GTid.y % 16 == 0)
		{
			output_lineardepth_mip4[pixel / 16] = lineardepth;
		}
	}

}
