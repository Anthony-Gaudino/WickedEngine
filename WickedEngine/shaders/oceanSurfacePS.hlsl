#define DISABLE_DECALS
#define DISABLE_ENVMAPS
#define DISABLE_TRANSPARENT_SHADOWMAP
#define DISABLE_SOFT_SHADOWMAP
#define TRANSPARENT
#define WATER
#include "globals.hlsli"
#include "oceanSurfaceHF.hlsli"
#include "objectHF.hlsli"

Texture2D<float4> texture_ocean_displacementmap : register(t0);
Texture2D<float4> texture_gradientmap : register(t1);
Texture2D<float4> texture_perlin : register(t2);

[earlydepthstencil]
float4 main(PSIn input) : SV_TARGET
{
#ifdef SHADOWMAPRENDERING
	float4 color = 1;
	color.rgb += caustics(input.uv * 10);
	color.a = input.pos.z; // secondary depth
	return color;
#else

	ShaderCamera camera = GetCameraIndexed(input.cameraIndex);
	float lineardepth = camera.IsOrtho() ? ((1 - input.pos.z) * camera.z_far) : input.pos.w;
	half4 color = xOceanWaterColor;
	float2 ScreenCoord = input.pos.xy * camera.internal_resolution_rcp;

	// The surface stops short of the eye so that the waterline the underwater
	// effects draw analytically is the only thing visible there.
	// `ocean_underwater_test_position` reconstructs on the same plane this
	// discards against, and gives a 2.5 cm band - the thin water mark seen
	// when breaking the surface. Rasterizing the mesh into that metre as well
	// puts a hard geometric edge over it.
	float4 pos2D = mul(camera.view_projection, float4(input.GetPos3D(), 1));
	pos2D.xyz /= pos2D.w;
	if (pos2D.z > ocean_waterline_handoff_depth())
		discard;

	float3 V = input.GetViewVector();
	float dist = length(V);
	V /= dist;
	
	uint2 pixel = input.pos.xy;

	half4 gradient = texture_gradientmap.Sample(sampler_aniso_wrap, input.uv);

	// Waves smaller than this pixel are gone from the sample above: the mip
	// chain averaged their SLOPES, which tends to zero and leaves a mirror
	// where there should be a rough sea. What they contributed comes back as
	// roughness instead, and the variance the mip swallowed is exactly that -
	// E[|s|^2] - |E[s]|^2, the two moments stored side by side by the folding
	// pass.
	//
	// Taken from the RAW sample: gradient.rg is replaced by perlin at distance
	// just below, after which the mean and the second moment describe
	// different surfaces and their difference means nothing.
	const float2 meanSlope = gradient.rg * xOceanGridLen * 0.5;
	float slopeVariance = max(0, gradient.b - dot(meanSlope, meanSlope));

	const float g_PerlinSize = 1;
	const float2 g_UVBase = 0;
	const float3 g_PerlinOctave = float3(1.12f, 0.59f, 0.23f);
	const float g_PerlinMovement = -GetTime() * xOceanTimeScale * 0.06f;
	const float3 g_PerlinGradient = float3(1.4f, 1.6f, 2.2f);

	float2 perlin_tc = input.uv * g_PerlinSize + g_UVBase;
	float2 perlin_tc0 = perlin_tc * g_PerlinOctave.x + g_PerlinMovement;
	float2 perlin_tc1 = perlin_tc * g_PerlinOctave.y + g_PerlinMovement;
	float2 perlin_tc2 = perlin_tc * g_PerlinOctave.z + g_PerlinMovement;

	float2 perlin_0 = texture_perlin.Sample(sampler_linear_wrap, perlin_tc0).xy;
	float2 perlin_1 = texture_perlin.Sample(sampler_linear_wrap, perlin_tc1).xy;
	float2 perlin_2 = texture_perlin.Sample(sampler_linear_wrap, perlin_tc2).xy;
	
	float2 perlin = (perlin_0 * g_PerlinGradient.x + perlin_1 * g_PerlinGradient.y + perlin_2 * g_PerlinGradient.z);
	//perlin = normalize(perlin);
	perlin = perlin * 0.000005 * xOceanWaveAmplitude;

	const float gradient_fade = smoothstep(
		OCEAN_GRADIENT_FADE.x, OCEAN_GRADIENT_FADE.y, dist);
	gradient.rg = lerp(gradient.rg, (half2)perlin, (half)gradient_fade);
	//gradient.rg = perlin;
	//return float4(gradient_fade.xxx,1);
	//return float4(gradient.rg, 0, 1);


	[branch]
	if (camera.texture_waterriples_index >= 0)
	{
		gradient.rg += bindless_textures_half4[descriptor_index(camera.texture_waterriples_index)].SampleLevel(sampler_linear_clamp, ScreenCoord, 0).rg * 0.0125;
	}

	// Add some small scale detail waves to make it look less uniform:
	const uint detail_count = 3;
	half2 gradient_detail = 0;
	float detail_variance = 0;
	for (uint i = 0; i < detail_count; ++i)
	{
		const float4 detail = texture_gradientmap.Sample(sampler_aniso_wrap, input.uv * pow(2.0, half(i + 1)));
		gradient_detail += (half2)detail.rg;

		const float2 detailMean = detail.rg * xOceanGridLen * 0.5;
		detail_variance += max(0, detail.b - dot(detailMean, detailMean));
	}
	gradient_detail /= half(detail_count);
	gradient.rg += gradient_detail;

	// Each of those taps is the same map at twice the frequency of the last, so
	// they go sub-pixel long before the base does - and being added after the
	// distance fade, nothing else ever damps them. Left alone they are the
	// sharpest aliasing on the water.
	//
	// Their variance scales the way their mean does: averaging n taps divides
	// the mean by n and the variance by n squared.
	slopeVariance +=
		detail_variance / float(detail_count * detail_count);

	const float bump_strength = 0.1;
	
	float4 water_plane = camera.reflection_plane;
	const bool camera_below_water = dot(float4(camera.position, 1), water_plane) < 0; 

	Surface surface;
	surface.init();
	surface.SetReceiveShadow(true);
	surface.pixel = pixel;
	float depth = input.pos.z;
	surface.albedo = color.rgb;
	surface.f0 = 0.02;

	// Specular anti-aliasing: the waves too small to draw are folded into the
	// GGX lobe rather than left to alias. The rule is alpha'^2 = alpha^2 +
	// 2*sigma^2 for a per-axis slope variance sigma^2, and the stored moment
	// already sums both axes, so the 2 is in it. alpha is roughness SQUARED
	// here - BRDF_GetSpecular squares before calling D_GGX - hence the fourth
	// root.
	//
	// Collapses to the base roughness wherever the variance is zero, which is
	// the whole near field. Widening this also blurs the reflected sky, for
	// free: EnvironmentReflection_Global picks its mip from roughness.
	const float oceanBaseRoughness = 0.1;
	const float oceanAlphaSq =
		sqr(sqr(oceanBaseRoughness)) + slopeVariance;
	surface.roughness = (half)sqrt(sqrt(oceanAlphaSq));
	surface.P = input.GetPos3D();
	surface.N = normalize(float3(gradient.x, xOceanTexelLength * 2, gradient.y));
	surface.V = V;
	// Kept at 1 for the gate it holds open, not for the wrap diffuse it looks
	// like it buys: ApplyLighting lerps the whole diffuse term away wherever
	// the refraction shows through, which for the ocean is everywhere but the
	// foam. What this actually does is keep NdotL_sss non-zero so that
	// BACK-LIT lights survive the early-out in each light function and reach
	// the scattering term at the bottom. Setting it to 0 silently kills the
	// glow through every wave facing away from the sun.
	surface.sss = 1;
	surface.sss_inv = 1.0f / ((1 + surface.sss) * (1 + surface.sss));

	// How much water this wave puts between the sun behind it and the eye. The
	// Jacobian fold runs high where the surface compresses into a crest - a
	// thin sheet - and to zero on the flat, so it inverts into a thickness.
	surface.water_thickness = (half)lerp(
		OCEAN_SUBSURFACE_THICKNESS.x,
		OCEAN_SUBSURFACE_THICKNESS.y,
		saturate(gradient.a));
	surface.update();

	Lighting lighting;
	lighting.create(0, 0, GetAmbient(surface.N), 0);

	TiledLighting(surface, lighting, GetFlatTileIndex(pixel, camera), camera);

	if (camera_below_water)
	{
		// Total internal reflection (Snell's window): light from above only
		// refracts through the surface within the critical angle for water
		// (~48.6 deg, cos = 0.661). Past it the surface is a perfect mirror, so
		// force full reflection - surface.F = 1 both cuts the sky refraction
		// (composited as refraction * (1 - F)) and shows the reflection of the
		// scene below. The soft band around the critical angle plus the wavy
		// normal gives the shimmering window edge rather than a hard circle:
		const float VdotN = abs(dot(surface.V, surface.N));
		surface.F = 1.0 - smoothstep(0.611, 0.711, VdotN);
	}
	
	[branch]
	if (camera.texture_reflection_index >= 0)
	{
		//REFLECTION
		float4 reflectionPos = mul(camera.reflection_view_projection, float4(surface.P, 1));
		float2 reflectionUV = clipspace_to_uv(reflectionPos.xy / reflectionPos.w) + surface.N.xz * bump_strength;
		half4 reflectiveColor = bindless_textures[descriptor_index(camera.texture_reflection_index)].SampleLevel(sampler_linear_mirror, reflectionUV, 0);
		[branch]
		if(camera.texture_reflection_depth_index >=0)
		{
			float reflectiveDepth = bindless_textures[descriptor_index(camera.texture_reflection_depth_index)].SampleLevel(sampler_point_clamp, reflectionUV, 0).r;
			float3 reflectivePosition = reconstruct_position(reflectionUV, reflectiveDepth, camera.reflection_inverse_view_projection);
			float water_depth = -dot(float4(reflectivePosition, 1), water_plane);
			water_depth += texture_ocean_displacementmap.SampleLevel(sampler_linear_wrap, reflectivePosition.xz * xOceanPatchSizeRecip, 0).z; // texture contains xzy!
			// Same medium as the refraction below: what the planar reflection
			// shows has crossed water_depth of water too. Composited through
			// MakeWaterFog rather than faded towards a flat colour, so the
			// water a dying reflection leaves behind is the same water
			// everything else scatters in - lit, and carrying the sky and the
			// sun, instead of the authored base colour on its own.
			//
			// The segment runs from the reflected point up to the surface, so
			// its own direction is the mirror of the view vector. The column
			// between the surface and the eye is deliberately NOT part of it -
			// that one is applied to the whole fragment by ApplyWaterFog at the
			// end - which is why the entry depth passed here is zero.
			//
			// Clamped at zero because a reflected point can land ABOVE the
			// water plane, and a negative path would drive the transmittance
			// past one and the in-scatter negative.
			//
			// God rays deliberately off: they are a screen space pattern swept
			// around where the sun projects, and a mirrored view has no such
			// point - switching them on paints the real sun's screen position
			// onto the reflection. underwaterCS.hlsl leaves them out of
			// Snell's window for the same reason.
			const WaterFog reflectionFog = MakeWaterFog(
				MakeWaterVolumetrics(1),
				max(0, water_depth),
				-reflect(V, float3(0, 1, 0)),
				0,
				ScreenCoord,
				uv_to_clipspace(ScreenCoord),
				false
			);

			reflectiveColor.rgb = (half3)(
				reflectiveColor.rgb * reflectionFog.transmittance
				+ reflectionFog.inscatter);
		}

		// remove planar reflection at high perturbation where it gets too inaccurate
		const float3 planar_reflection_vector_flat = reflect(V, float3(0, 1, 0));
		const float3 planar_reflection_vector = reflect(V, surface.N);
		const float planar_factor = smoothstep(camera_below_water ? 0.9 : 0.95,1.0,abs(dot(planar_reflection_vector_flat, planar_reflection_vector)));
		//return float4(planar_factor.xxx,1);

		lighting.indirect.specular += reflectiveColor.rgb * surface.F * saturate(dist * 0.1) * planar_factor; // fade out very close to camera, doesn't look good
	}
	else
	{
		lighting.indirect.specular += EnvironmentReflection_Global(surface);
	}

	float water_depth = FLT_MAX;

	[branch]
	if (camera.texture_refraction_index >= 0)
	{
		// Water refraction:
		Texture2D texture_refraction = bindless_textures[descriptor_index(camera.texture_refraction_index)];
		// First sample using full perturbation:
		float2 refraction_uv = ScreenCoord.xy + surface.N.xz * bump_strength;
		float refraction_depth = find_max_depth(refraction_uv, 2, 2);
		float3 refraction_position = reconstruct_position(refraction_uv, refraction_depth);
		water_depth = -dot(float4(refraction_position, 1), water_plane);
		water_depth += texture_ocean_displacementmap.SampleLevel(sampler_linear_wrap, refraction_position.xz * xOceanPatchSizeRecip, 0).z; // texture contains xzy!
		if (camera_below_water && V.y < 0)
			water_depth = -water_depth;
		if (water_depth <= 0)
		{
			// Above water, fill holes by taking unperturbed sample:
			refraction_uv = ScreenCoord.xy;
		}
		else
		{
			// Below water, compute perturbation according to first sample water depth:
			refraction_uv = ScreenCoord.xy + surface.N.xz * bump_strength * saturate(1 - exp(-water_depth));
		}
		surface.refraction.rgb = texture_refraction.SampleLevel(sampler_linear_mirror, refraction_uv, 0).rgb;
		// Recompute depth params again with actual perturbation:
		refraction_depth = texture_depth.SampleLevel(sampler_point_clamp, refraction_uv, 0);
		refraction_position = reconstruct_position(refraction_uv, refraction_depth);
		water_depth = max(water_depth, -dot(float4(refraction_position, 1), water_plane));
		water_depth += texture_ocean_displacementmap.SampleLevel(sampler_linear_wrap, refraction_position.xz * xOceanPatchSizeRecip, 0).z; // texture contains xzy!
		if (camera_below_water && V.y < 0)
			water_depth = -water_depth;
		// The refraction is handed on WHOLE - no absorption is applied to it
		// here, and this is the point of the whole design.
		//
		// It used to be attenuated by exp(-water_path * sigmaT), with
		// water_path derived from texture_depth. That buffer is the opaque
		// prepass: transparents are not in it, and neither is the sky. So a
		// submerged particle, trail or sprite with nothing solid behind it was
		// measured against the far plane, the transmittance underflowed, and it
		// was extinguished - even though it WAS in the scene copy being sampled
		// and would otherwise have shown through. A transmissive mesh went the
		// same way, which is why it read as far too clear from above the water
		// and correct from below.
		//
		// There is no fixing that from here: one screen depth cannot describe
		// what a stack of transparents is at. So everything in the scene copy
		// now carries its own absorption instead, applied over its own path as
		// it was drawn - GetWaterFog does that for opaque geometry, for every
		// transparent, and for the sky, which is what supplies the "infinitely
		// deep water" backdrop that the far-plane depth used to fake.
		//
		// Fresnel still applies: ApplyLighting composites this as
		// refraction * (1 - F), so the surface still reflects its share away.
		surface.refraction.a = 1;
		color.a = 1;

		// The one stretch of water that no fragment fogs for itself: the crest
		// this one stands on.
		//
		// Every fragment carries the absorption of its own path, which covers
		// the refraction target whenever the target is IN the water - and since
		// the fog clips against the wave surface rather than the still plane,
		// that already includes the crest standing over the target. What it
		// cannot cover is a target that was never in the water at all. The sky
		// and anything above the surface have no water path to fog, so a crest
		// seen against them handed their radiance on undimmed and read as a
		// glass blob instead of as a body of water.
		//
		// So the test is whether the target is under water at all, and the term
		// is what is missing when it is not: the slab standing above the still
		// plane here, crossed at the angle the ray was refracted to. Adding it
		// only in that case is what keeps it from double counting against the
		// column the target already fogged.
		const float crest_height =
			max(0, surface.P.y - GetWeather().ocean.water_height);
		const float refraction_height = refraction_position.y
			- ocean_drawn_surface_height(refraction_position);

		[branch]
		if (!camera_below_water && crest_height > 0 && refraction_height >= 0)
		{
			const WaterFog crestFog = MakeWaterFog(
				MakeWaterVolumetrics(1),
				SubmergedViewPath(crest_height, V),
				V,
				0,
				ScreenCoord,
				uv_to_clipspace(ScreenCoord),
				false
			);

			surface.refraction.rgb = (half3)(
				surface.refraction.rgb * crestFog.transmittance
				+ crestFog.inscatter);
		}
	}
	
#if 1
	[branch]
	if (camera.texture_depth_index >= 0)
	{
		// FOAM:
		// How much water stands over the sea bed at this pixel, measured
		// VERTICALLY between the displaced surface and whatever lies beneath.
		//
		// This used to be the difference of the two LINEAR DEPTHS, which is a
		// distance along the view ray: for a vertical gap h at an angle theta
		// from straight down it reads h / cos(theta). So the band was narrow
		// seen along the shore and grew as the camera tilted down, spreading
		// foam over water that is not shallow at all.
		//
		// Still taken from the DISPLACED surface rather than the still plane,
		// which is what the original note here was protecting: a crest standing
		// over the sea bed really does have less water beneath it. What
		// changed is only WHERE that displaced height is read - see below.
		const float shore_device_depth = texture_depth[pixel];

		// Is there anything solid under this pixel at all? Reverse-Z clears to
		// zero, so sky and open water read as the far plane.
		//
		// Without this test the shallow-water measure below is nonsense exactly
		// where there is no sea bed: near the horizon the water surface is
		// itself almost at the far plane, so the gap between it and "the far
		// plane" is small and reads as shallow. Foam then spreads over open
		// sea, and spreads FURTHER as the camera climbs and more of the view is
		// distant water.
		//
		// The same trap as the refraction's old depth source - the opaque
		// prepass says nothing whatever about what is not in it.
		const bool shore_has_geometry = shore_device_depth > 0;

		const float3 shore_position =
			reconstruct_position(ScreenCoord, shore_device_depth);

		// Measured between the sea bed and the surface standing OVER THE SEA
		// BED, not the surface at this fragment. Those are two different
		// places: the fragment is wherever this pixel's piece of water surface
		// happens to sit, and at a grazing angle - or on a crest raised
		// between the eye and the shore - it can be a long way horizontally
		// from the ground the pixel is actually looking at. Subtracting two
		// heights sampled at different XZ is not a water column at all, and a
		// tall crest in front of dry sand reports metres of water standing
		// over ground the sea never reaches.
		//
		// Signed on purpose. Where the surface over the ground sits BELOW it
		// the ground is dry, and dry ground has no shoreline foam however much
		// water stands somewhere else along the view ray.
		const float shore_depth =
			ocean_drawn_surface_height(shore_position) - shore_position.y;
		const bool shore_is_submerged = shore_has_geometry && shore_depth > 0;

		float foam_shore =
			shore_is_submerged ? saturate(exp(-shore_depth * 2)) : 0;
		float foam_wave = pow(saturate(gradient.a), 4) * saturate(exp(-water_depth * 0.1));
		float foam_combined = saturate(foam_shore + foam_wave);
		float foam = smoothstep(0.5, 0.6, saturate(foam_combined + 0.1));

		// How much world the pixel covers. These noises are analytic, so unlike
		// the gradient map they have no mip to fall back on and the variance
		// above cannot reach them - once an octave's features are finer than
		// this, it is sampled at random and reads as per-pixel sparkle.
		//
		// MUST be computed here, not inside the branch below: that branch is
		// not uniform across a quad, and derivatives taken under divergent flow
		// are undefined.
		const float2 footprintDelta =
			max(abs(ddx(surface.P.xz)), abs(ddy(surface.P.xz)));
		const float footprint = max(footprintDelta.x, footprintDelta.y);

		// Weight per octave, dropping each as it passes Nyquist for this pixel.
		// The coarsest is never dropped - it has metre-scale features, so it
		// outlives any sensible view distance, and keeping it means the sum
		// below always has something in it.
		const float3 foam_octave_weight = float3(
			1,
			saturate(1 - footprint * 4),
			saturate(1 - footprint * 8));

		// Renormalise by the weight that survived, so losing an octave costs
		// detail rather than brightness. Every octave is the same function at a
		// different scale and so has the same mean, which is what makes
		// rescaling the sum preserve it.
		const float foam_octave_norm = 3.0 /
			(foam_octave_weight.x + foam_octave_weight.y + foam_octave_weight.z);

		// The procedural noise below is multiplied by foam_combined, so it only
		// contributes where there is actually foam (shorelines and wave crests).
		// Skip the 6 noise evaluations for the open-water majority where
		// foam_combined is ~0; the branch is spatially coherent so warps converge.
		[branch]
		if (foam_combined > 0.002)
		{
			float foam_simplex = 0;
			foam_simplex += foam_octave_weight.x * smoothstep(0, 0.8, noise_simplex_2D(surface.P.xz * 1 + GetTime()));
			foam_simplex += foam_octave_weight.y * smoothstep(0, 0.8, noise_simplex_2D(surface.P.xz * 2 + GetTime()));
			foam_simplex += foam_octave_weight.z * smoothstep(0, 0.8, noise_simplex_2D(surface.P.zx * 4 - GetTime() * 2));
			foam_simplex *= foam_octave_norm;
			float foam_voronoi = 0;
			foam_voronoi += foam_octave_weight.x * smoothstep(0.5, 0.8, noise_voronoi(surface.P.xz * 1, GetTime()).x);
			foam_voronoi += foam_octave_weight.y * smoothstep(0.5, 0.8, noise_voronoi(surface.P.xz * 2, GetTime()).x);
			foam_voronoi += foam_octave_weight.z * smoothstep(0.5, 0.8, noise_voronoi(surface.P.xz * 4, GetTime()).x);
			foam_voronoi *= foam_octave_norm;
			foam += foam_voronoi * foam_simplex * foam_combined;
		}
		foam *= 2;
		foam = saturate(foam);
		surface.albedo = lerp(surface.albedo, 0.6, foam);
		surface.refraction.a *= 1 - foam;
		// Same guard: with no sea bed under water there is no shoreline to
		// restore the refraction over.
		surface.refraction.a = saturate(surface.refraction.a +
			(shore_is_submerged ? saturate(exp(-shore_depth * 4)) : 0));
	}
#endif

	ApplyLighting(surface, lighting, color);
	
	// Blend out at distance:
	float far_fade = saturate(1 - saturate(dist / camera.z_far - 0.8) * 5.0); // fade will be on edge and inwards 20%
	color.a = far_fade;

	// What ApplyLighting just blended in from behind the surface. It reached
	// this fragment already fogged over its own, longer path - by the air when
	// it was drawn, and by the water it crossed to get out - so every applier
	// below is told to leave it alone rather than fog the same stretch twice.
	//
	// Seen from below that background is the sky; seen from above it is
	// whatever is under the water.
	const half3 background = (half3)(
		surface.refraction.rgb * (1 - surface.F) * surface.refraction.a);

	ApplyAerialPerspective(ScreenCoord, surface.P, background, color);

	ApplyFog(dist, V, background, color);

	// AtSurface: this fragment IS the interface, so what lies between it and
	// the eye is water exactly when the eye is under it. Classifying it by
	// height instead would have it fog itself in every wave trough, where the
	// displaced surface dips below the still plane.
	ApplyWaterFogAtSurface(ScreenCoord, surface.P, background, color);

	return saturateMediump(color);

#endif // SHADOWMAPRENDERING
}

