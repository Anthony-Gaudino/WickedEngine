/**********************************************************************
Copyright (c) 2021 Advanced Micro Devices, Inc. All rights reserved.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
********************************************************************/

#ifndef FFX_DNSR_SHADOWS_PREPARESHADOWMASK_HLSL
#define FFX_DNSR_SHADOWS_PREPARESHADOWMASK_HLSL

#include "ffx_denoiser_shadows_util.h"

groupshared uint g_FFX_DNSR_Shadows_lane_mask_or;

uint FFX_DNSR_Shadows_GroupReduceOr(uint lane_mask, uint2 gtid)
{
#if defined(WAVE_OPS_SUPPORTED) || defined(__spirv__) || defined(__PSSL__) || defined(__XBOX_ENABLE_WAVEOPS)
    // Prefer native wave reduction when the wave covers the full threadgroup (FFX assumes 64 lanes).
    if (WaveGetLaneCount() == 64)
    {
        return WaveActiveBitOr(lane_mask);
    }
#endif // wave intrinsics available

    // Fallback path for platforms (such as Metal) where wave intrinsics either are unavailable
    // or the wave size does not span the full threadgroup.
    GroupMemoryBarrierWithGroupSync();
    if (all(gtid == 0))
    {
        g_FFX_DNSR_Shadows_lane_mask_or = 0;
    }
    GroupMemoryBarrierWithGroupSync();

    InterlockedOr(g_FFX_DNSR_Shadows_lane_mask_or, lane_mask);

    GroupMemoryBarrierWithGroupSync();
    return g_FFX_DNSR_Shadows_lane_mask_or;
}

void FFX_DNSR_Shadows_CopyResult(uint2 gtid, uint2 gid)
{
    const uint2 did = gid * uint2(8, 4) + gtid;
    const uint linear_tile_index = FFX_DNSR_Shadows_LinearTileIndex(gid, FFX_DNSR_Shadows_GetBufferDimensions().x);
    const bool hit_light = FFX_DNSR_Shadows_HitsLight(did, gtid, gid);
    const uint lane_mask = hit_light ? FFX_DNSR_Shadows_GetBitMaskFromPixelPosition(did) : 0;
    const uint merged_mask = FFX_DNSR_Shadows_GroupReduceOr(lane_mask, gtid);

    if (all(gtid == 0))
    {
        FFX_DNSR_Shadows_WriteMask(linear_tile_index, merged_mask);
    }
}
 
void FFX_DNSR_Shadows_PrepareShadowMask(uint2 gtid, uint2 gid)
{
    gid *= 4;
    uint2 tile_dimensions = (FFX_DNSR_Shadows_GetBufferDimensions() + uint2(7, 3)) / uint2(8, 4);

    for (int i = 0; i < 4; ++i)
    {
        for (int j = 0; j < 4; ++j)
        {
            uint2 tile_id = uint2(gid.x + i, gid.y + j);
            tile_id = clamp(tile_id, 0, tile_dimensions - 1);
            FFX_DNSR_Shadows_CopyResult(gtid, tile_id);
        }
    }
}

#endif
