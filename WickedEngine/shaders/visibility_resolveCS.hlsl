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
#else
RWTexture2D<uint> output_primitiveID : register(u13);
#endif // VISIBILITY_MSAA

inline bool IntersectsSceneAABB(const RayDesc ray)
{
	float3 invdir = rcp(ray.Direction);
	float3 t0 = (GetScene().aabb_min - ray.Origin) * invdir;
	float3 t1 = (GetScene().aabb_max - ray.Origin) * invdir;
	float3 tmin3 = min(t0, t1);
	float3 tmax3 = max(t0, t1);
	float tNear = max(max(tmin3.x, tmin3.y), tmin3.z);
	float tFar = min(min(tmax3.x, tmax3.y), tmax3.z);
	return tFar >= max(tNear, 0);
}

bool TracePrimitiveFallback(RayDesc ray, uint2 pixel, uint groupIndex, inout uint primitiveID, out PrimitiveID prim)
{
	prim.init();
	if ((GetFrame().options & OPTION_BIT_PRIMITIVEID_FALLBACK) == 0)
	{
		return false;
	}
	if (GetScene().BVH_primitives < 0)
	{
		return false;
	}
	ray.Origin -= ray.Direction * 0.01f;
	ray.TMin = 0.0f;
	if (!IntersectsSceneAABB(ray))
	{
		return false;
	}

	RNG rng;
	rng.init(pixel, GetFrame().frame_count);
	RayHit hit = TraceRay_Closest(ray, ~0u, rng, groupIndex);
	if (hit.distance == FLT_MAX)
	{
		return false;
	}

	prim = hit.primitiveID;
	primitiveID = prim.pack();
	return true;
}

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

	float depth = 1; // default to far plane when nothing was resolved
	uint bin = 0;
	bool lane_has_work = false;
	if (pixel_valid)
	{
		PrimitiveID prim;
		prim.init();

		Surface surface;
		surface.init();

		bool surface_loaded = false;
		if (primitiveID != 0)
		{
			prim.unpack(primitiveID);
			surface_loaded = surface.load(prim, ray.Origin, ray.Direction);
		}

		if (!surface_loaded && TracePrimitiveFallback(ray, pixel, groupIndex, primitiveID, prim))
		{
			surface.init();
			surface_loaded = surface.load(prim, ray.Origin, ray.Direction);
		}

		if (surface_loaded)
		{
			depth = compute_inverse_lineardepth(surface.hit_depth);
			depth = saturate(depth);
			depth = max(depth, 1e-6);	// keep primary depth distinct from sky (which stays at 0)

			bin = surface.material.GetShaderType();
			lane_has_work = true;
		}
		else
		{
			// sky or nothing hit
			depth = 0;
			bin = SHADERTYPE_BIN_COUNT;
			lane_has_work = true;
			primitiveID = 0;
		}

		if (lane_has_work)
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

		#ifdef VISIBILITY_MSAA
			output_primitiveID[pixel] = primitiveID;
		#else
			if (GetFrame().options & OPTION_BIT_PRIMITIVEID_FALLBACK)
			{
				output_primitiveID[pixel] = primitiveID;
			}
		#endif // VISIBILITY_MSAA
	}

	GroupMemoryBarrierWithGroupSync();
	if (groupIndex == 0)
	{
		uint mask = 0;
		[unroll]
		for (uint i = 0; i <= SHADERTYPE_BIN_COUNT; ++i)
		{
			if ((local_bin_execution_mask_0[i] | local_bin_execution_mask_1[i]) != 0)
			{
				mask |= 1u << i;
			}
		}
		local_bin_mask = mask;
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
			tile.execution_mask[0] = local_bin_execution_mask_0[groupIndex];
			tile.execution_mask[1] = local_bin_execution_mask_1[groupIndex];
			output_binned_tiles[bin_tile_list_offset + tile_offset] = tile;
		}
	}

	// Downsample depths:
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
