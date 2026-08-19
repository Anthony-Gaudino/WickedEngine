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

float4 main(PSIn input) : SV_TARGET
{
#ifdef SHADOWMAPRENDERING
	float4 color = 1;
	color.rgb = caustics_modulation(input.uv * WATER_CAUSTIC_TILES_PER_PATCH);
	color.a = input.pos.z; // secondary depth
	return color;
#else

	ShaderCamera camera = GetCameraIndexed(input.cameraIndex);
	float lineardepth = camera.IsOrtho() ? ((1 - input.pos.z) * camera.z_far) : input.pos.w;
	half4 color = 0;
	float2 ScreenCoord = input.pos.xy * camera.internal_resolution_rcp;

	// The surface stops short of the eye so that the waterline the underwater
	// effects draw analytically is the only thing visible there.
	// `ocean_underwater_test_position` reconstructs on the same plane this
	// discards against, and gives a 2.5 cm band - the thin water mark seen
	// when breaking the surface. Rasterizing the mesh into that metre as well
	// puts a hard geometric edge over it.
	//
	// The depth write has to follow this discard, never precede it: the surface
	// writes into the main depth buffer, so a fragment that wrote before being
	// rejected here would occlude everything drawn behind it in that metre.
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
	// Which side of the water the eye is on, and everything about how this
	// surface is shaded turns on it: the normal faces the medium the viewer is
	// in, the Fresnel switches to the dense-side curve, and past the critical
	// angle the underside becomes a mirror.
	//
	// **Measured against the DRAWN surface over the camera, never the still
	// plane.** Everything else about this surface is measured against the drawn
	// one, and the two disagree by the whole wave height: an eye sitting in a
	// trough is below the plane while plainly in the air, and a plane test
	// hands it the entire submerged treatment - flipped normal, dense-side
	// Fresnel, total internal reflection - which is what shades crests
	// differently for a camera that never got wet. An eye inside a crest is the
	// same disagreement the other way up.
	const bool camera_below_water =
		camera.position.y < ocean_drawn_surface_height(camera.position);

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

	// The wave normal is built pointing up out of the water, which is the wrong
	// side to shade from when the eye is under the surface: every dot taken
	// against it comes out negative, NdotV saturates to zero and the specular
	// lobe goes with it, so a submerged lamp casts no glint at all on the
	// underside. Face the normal into the medium the viewer is actually in.
	//
	// Costs nothing elsewhere: reflect() is even in N, so the planar reflection
	// vectors below are untouched, the critical angle test takes an abs, and
	// the refraction offset only has its direction mirrored. What does change,
	// correctly, is that the ambient now arrives from below.
	if (camera_below_water)
	{
		surface.N = -surface.N;

		// Water against air, seen from the water. Switches the direct specular
		// Fresnel to the dense-side curve, which reaches total reflection at
		// the critical angle rather than sitting at the two percent an air-side
		// f0 returns - the difference between a submerged lamp glinting off the
		// underside and it being invisible there.
		surface.internal_ior = 1.333;
	}

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
		// The reflection pass fogs its own fragments over their path up to the
		// surface, and ApplyWaterFogAtSurface below adds the column from the
		// surface to the eye. Attenuating here would count one of them twice.

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

		// Reject against the neighbourhood the SCENE COPY resolves, not the one
		// the depth buffer does. The copy is a fraction of the screen
		// resolution, so a silhouette in it is several full resolution pixels
		// wide: a UV landing merely NEAR a solid, without landing on it, still
		// reads a texel carrying that solid, and paints it into the water. Two
		// texels of reach covers the copy's own footprint and the second texel
		// bilinear filtering pulls in.
		float2 refraction_dim;
		texture_refraction.GetDimensions(refraction_dim.x, refraction_dim.y);
		const float2 refraction_texel = rcp(refraction_dim);

		float refraction_depth = find_max_depth(refraction_uv, 2, refraction_texel);
		float3 refraction_position = reconstruct_position(refraction_uv, refraction_depth);
		water_depth = -dot(float4(refraction_position, 1), water_plane);
		water_depth += texture_ocean_displacementmap.SampleLevel(sampler_linear_wrap, refraction_position.xz * xOceanPatchSizeRecip, 0).z; // texture contains xzy!
		if (camera_below_water && V.y < 0)
			water_depth = -water_depth;
		// How much water this ray crossed to reach what it refracted, which is
		// what the surface is entitled to bend.
		//
		// **Walked along the segment, never tested against the still plane.**
		// Which side of a plane the target lies on is a different question and
		// gives the wrong answer to this one in both directions: a crest stands
		// in front of things that are perfectly dry, so a target above the
		// plane is not owed zero bending; and a target below the plane can have
		// nothing but air in front of it once the wave over it has moved on.
		// Answering with the plane abandoned the refraction outright wherever
		// it read zero - a shore, a bed rising out of the water, anything
		// standing in it - and that is a whole class of surface that simply did
		// not refract.
		//
		// Zero is now a real answer rather than an artefact of the test: no
		// water was crossed, so there is nothing to bend, and the unperturbed
		// sample below is right.
		float refracted_water = TraceWaterSegment(
			GetCamera().position,
			refraction_position,
			WATER_SEGMENT_STABLE_PHASE
		).submerged;

		if (refracted_water <= 0)
		{
			// Nothing in front of it to bend the ray, so fill holes by taking
			// the unperturbed sample:
			refraction_uv = ScreenCoord.xy;
		}
		else
		{
			// Perturbation according to the water the first sample is behind:
			refraction_uv = ScreenCoord.xy + surface.N.xz * bump_strength * saturate(1 - exp(-refracted_water));
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

		// The stretch of water between this fragment and what it refracted that
		// nothing in the frame has accounted for.
		//
		// Every fragment carries the absorption of its own path, so a submerged
		// target arrives already fogged - but over the column standing over
		// ITSELF, from the surface above it down to where it is. **This ray
		// entered the water here**, at this fragment, which is somewhere else
		// entirely and, on a crest, several metres higher. The drop between the
		// two is water the ray crossed and no one measured.
		//
		// Added to what the target already charged, this comes to the whole
		// descending leg and no more: the surface height over the target
		// appears in both with opposite signs and cancels, leaving this
		// fragment's height above the target at the refracted angle. That
		// cancellation is also why the target's own side needs no test here. A
		// target ABOVE the water charged nothing, and measuring from the
		// surface over it gives the column standing above the still plane,
		// which is what a crest seen against the sky is owed.
		//
		// Zero on a flat sea, where both are the still plane, and zero looking
		// straight down a crest, where the target's own surface is that same
		// crest - which is why a crest only leaks where it has a far side to
		// see past.
		//
		// Floored at zero rather than allowed to go negative: a target lying
		// under a crest of its own counted more water than this ray crossed,
		// and a sample that has already been fogged cannot be un-fogged.
		const float uncounted_depth = max(0, surface.P.y
			- ocean_drawn_surface_height(refraction_position));

		[branch]
		if (!camera_below_water && uncounted_depth > 0)
		{
			const WaterFog waterBehind = MakeWaterFog(
				MakeWaterVolumetrics(1),
				SubmergedViewPath(uncounted_depth, V),
				V,
				0,
				ScreenCoord,
				uv_to_clipspace(ScreenCoord),
				false,
				// This branch is only taken with the camera above the surface,
				// so no froxel column holds water here: this is the only
				// description of the sun in that column, and it keeps all of it.
				1
			);

			surface.refraction.rgb = (half3)(
				surface.refraction.rgb * waterBehind.transmittance
				+ waterBehind.inscatter);
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
		const float shore_water_column =
			ocean_drawn_surface_height(shore_position) - shore_position.y;
		const bool shore_is_submerged =
			shore_has_geometry && shore_water_column > 0;

		// Whether this pixel grows shoreline foam at all: a sea bed with water
		// standing over it, and the feature left switched on.
		const bool shore_foam_present =
			shore_is_submerged && GetWeather().ocean.IsShoreFoam();

		// Reciprocal of the authored e-folding depth, which is what the
		// exponent below wants. Guarded against zero for the paths the
		// editor's slider range does not cover - a script or a hand written
		// scene can set anything.
		const float shore_foam_falloff =
			rcp(max(GetWeather().ocean.shore_foam_width, 0.001));

		// How high THIS piece of surface stands over that ground, which is the
		// other way the water here can fail to be shallow. The column above
		// establishes that the ground is under water and how deeply, but it
		// describes the surface sitting over the GROUND - and a crest is not
		// there. Standing eight metres up the face of a wave in front of a
		// shoal, this fragment is deep water by any measure that asks about
		// this fragment, while the column over the shoal goes on reporting
		// ankle depth and paints the shoal's own outline onto the wave.
		//
		// Taken as the larger of the two rather than replacing it: the column
		// is what knows the difference between shallow water and dry land, and
		// on its own this one would grow foam over a beach whenever a wave rose
		// higher than the sand.
		const float shore_depth = max(
			shore_water_column, surface.P.y - shore_position.y);

		float foam_shore = shore_foam_present
			? saturate(exp(-shore_depth * shore_foam_falloff))
				* GetWeather().ocean.shore_foam_strength
			: 0;
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
		surface.albedo = lerp(surface.albedo, (half3)xOceanFoamColor.rgb, foam);
		surface.refraction.a *= 1 - foam * (half)xOceanFoamColor.a;
		// Same guard: with no sea bed under water there is no shoreline to
		// restore the refraction over.
		//
		// Reaches half as far as the foam and is scaled by the same gain, so
		// the two cannot be set against each other: thinning the foam without
		// thinning this would brighten water that has no foam left on it.
		surface.refraction.a = saturate(surface.refraction.a +
			(shore_foam_present
				? saturate(exp(-shore_depth * shore_foam_falloff * 2))
					* GetWeather().ocean.shore_foam_strength
				: 0));
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

	ApplyVolumetricLight(
		ScreenCoord,
		surface.P,
		(half3)((1 - surface.F) * surface.refraction.a),
		color
	);

	return saturateMediump(color);

#endif // SHADOWMAPRENDERING
}

