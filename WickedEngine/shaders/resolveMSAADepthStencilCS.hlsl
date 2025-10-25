#include "globals.hlsli"
// Resolve MSAA depth buffer to a non MSAA texture

Texture2DMS<float> input : register(t0);

RWTexture2D<float> output : register(u0);

[numthreads(8, 8, 1)]
void main(uint3 dispatchThreadId : SV_DispatchThreadID)
{
	uint2 dim;
	uint sampleCount;
	input.GetDimensions(dim.x, dim.y, sampleCount);
	if (dispatchThreadId.x >= dim.x || dispatchThreadId.y >= dim.y)
	{
		return;
	}

	float resolved_depth = 1.0;
	bool has_sample = false;
	for (uint i = 0; i < sampleCount; ++i)
	{
		float sample_depth = input.Load(dispatchThreadId.xy, i).r;
		// Treat zero as "no geometry" just like the depth prepass clear value.
		if (sample_depth <= 0)
		{
			continue;
		}
		if (!has_sample)
		{
			has_sample = true;
			resolved_depth = sample_depth;
		}
		else
		{
			resolved_depth = min(resolved_depth, sample_depth);
		}
	}

	output[dispatchThreadId.xy] = has_sample ? resolved_depth : 0;
}
