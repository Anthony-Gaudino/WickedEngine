#define OBJECTSHADER_LAYOUT_COMMON
#define DISABLE_WIND
#ifndef VOXELIZATION_GEOMETRY_SHADER_ENABLED
#define OBJECTSHADER_USE_CAMERAINDEX
#define OBJECTSHADER_USE_VIEWPORTARRAYINDEX
#endif // VOXELIZATION_GEOMETRY_SHADER_ENABLED
#include "objectHF.hlsli"

struct VSOut
{
	float4 pos : SV_POSITION;
	float4 uvsets : UVSETS;
	half4 color : COLOR;
	float3 N : NORMAL;
#ifndef VOXELIZATION_GEOMETRY_SHADER_ENABLED
	float3 P : POSITION3D;
#endif // VOXELIZATION_GEOMETRY_SHADER_ENABLED
};

VSOut main(VertexInput input)
{
	VertexSurface surface;
	surface.create(GetMaterial(), input);

	VSOut Out;
	Out.pos = surface.position;
	Out.color = surface.color;
	Out.uvsets = surface.uvsets;
	Out.N = surface.normal;

#ifndef VOXELIZATION_GEOMETRY_SHADER_ENABLED
	Out.P = surface.position.xyz;

	ShaderCamera camera = GetCameraIndexed(input.GetInstancePointer().GetCameraIndex());
	Out.pos = mul(camera.view_projection, Out.pos);
#endif // VOXELIZATION_GEOMETRY_SHADER_ENABLED

	return Out;
}

