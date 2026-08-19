/**
 * Paints one stage of the refraction lookup instead of the water.
 *
 * Each mode isolates one link in the chain that fills a pixel the scene copy
 * cannot answer for, so a single frame says which link is at fault rather than
 * leaving it to be argued from the shaded result.
 *
 * - `1` - belief, and where it was lost. **White** where the copy answered.
 *   **Red** where something nearer stands in the ray's way, **blue** where the
 *   coordinate ran off the edge of the copy, **green** where the ray projected
 *   from behind the eye. Red and blue have entirely different cures - one
 *   wants the ray traced, the other wants the copy rendered wider - so the
 *   share of each is what decides which is worth paying for.
 * - `4` - which way out of the branch chain the fragment took. **Red** took
 *   the no-water path and does not bend at all, **green** the below-water
 *   perturbation, **blue** the refraction proper. A wrong result says nothing
 *   until it says which of the three produced it.
 * - `5` - the refraction flattened to magenta, everything else shaded as
 *   normal. Anything still structured on the water is drawn by the surface's
 *   own shading - reflection, specular, foam - and not by the lookup at all.
 * - `2` - the fog's in-scatter on its own, which is the whole of what a blind
 *   pixel is filled with. Black here means the medium returns no veiling
 *   radiance, and the fill can never work while it does.
 * - `3` - the refraction after the fog has been applied. Coloured here but
 *   black on screen puts the fault after this shader's refraction entirely.
 *
 * 0 in anything but a diagnosis.
 */
#define OCEAN_REFRACTION_DEBUG 0

/**
 * Reflectance of the sea bed where it has to be invented rather than read.
 *
 * Stands in for a surface the eye cannot see, behind whatever blocked the
 * refracted ray. Deliberately dull: a guess that comes out too dark reads as
 * shadow under the object, which is roughly what belongs there, while one that
 * comes out too bright reads as a rim drawn around it.
 */
static const float UNSEEN_BED_ALBEDO = 0.2F;

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
		// Measured along this fragment's OWN view ray, which is what "how much
		// water stands in front of what it is looking through" means.
		//
		// Never a screen-space slide of it. bump_strength is a tenth of the
		// screen, so a wave normal of any strength moves the tap tens of pixels
		// sideways, and everything below is decided there instead of here - a
		// tap that lands on a dry object switches this fragment's refraction
		// off from that far away, and the unperturbed sample it then falls back
		// to is read beside that object, where the copy already holds it.
		float2 refraction_uv = ScreenCoord.xy;

		float2 refraction_dim;
		texture_refraction.GetDimensions(refraction_dim.x, refraction_dim.y);
		const float2 refraction_texel = rcp(refraction_dim);

		// Point sampled, and it has to stay that way.
		//
		// This measures how much water stands in front of the ray's target,
		// which decides both whether the surface bends at all and how far down
		// it aims. Widened to the nearest thing within a neighbourhood, it
		// reports a dry object several screen pixels before the coordinate
		// reaches it, and a target with no water in front of it switches the
		// refraction off - so a ring around everything standing out of the
		// water samples the copy at the pixel's own coordinate, which at a
		// quarter of screen resolution is exactly where that object is drawn.
		// The object then rims itself in its own colour.
		//
		// Filtering across a silhouette is a real hazard, but it belongs where
		// the colour is fetched and is handled there.
		float refraction_depth = texture_depth.SampleLevel(
			sampler_point_clamp, refraction_uv, 0);
		float3 refraction_position = reconstruct_position(refraction_uv, refraction_depth);
		water_depth = -dot(float4(refraction_position, 1), water_plane);
		water_depth += texture_ocean_displacementmap.SampleLevel(sampler_linear_wrap, refraction_position.xz * xOceanPatchSizeRecip, 0).z; // texture contains xzy!
		if (camera_below_water && V.y < 0)
			water_depth = -water_depth;
		// How much water stands in front of what this ray refracted, which is
		// what the surface is entitled to bend. Below the still plane that is
		// the target's depth, which is what this has always ramped on. Above it
		// the answer is not zero: a crest can stand in front of something that
		// is perfectly dry, and the ray cut a chord through the wave to reach
		// it. On flat water the two are the same number, since there is no
		// water above the plane to cross.
		//
		// A target with no water in front of it at all still gives zero, which
		// is what the unperturbed fallback below is for: a perturbed UV that
		// wandered onto something nearer than the water must not drag it into
		// the refraction.
		float refracted_water = water_depth;

		[branch]
		if (refracted_water <= 0 && GetCamera().IsWaterSegmentModel())
		{
			refracted_water = TraceWaterSegment(
				GetCamera().position,
				refraction_position,
				WATER_SEGMENT_STABLE_PHASE
			).submerged;
		}

		// How much of what the lookup returns can be believed, and how deep to
		// call the water where it cannot be believed at all.
		//
		// The copy is a screen buffer, so it can only answer for a point that
		// was actually drawn at the coordinate that point projects to. Where
		// something else stands there instead, the answer is NOT another part
		// of the screen: every coordinate that has scene in it is a second copy
		// of that scene, which is the whole family of duplication bugs this
		// lookup has produced. The answer is the water itself, with nothing
		// behind it, which is the one thing that cannot duplicate anything.
		float refraction_confidence = 1;
		float refraction_blind_depth = 0;

		// Which of the three reasons took the belief away, for the diagnostic
		// below. Red, green and blue in that order so one frame separates the
		// causes rather than only measuring how much was lost - they have
		// different cures, and only the count of each says which is worth
		// paying for.
		float3 refraction_blind_reason = 0;

		// Which of the three ways out of the chain below this fragment took,
		// for the diagnostic. They disagree about where to look and about
		// whether to bend at all, so a result that is wrong says nothing until
		// it says which one produced it.
		float3 refraction_branch = 0;

		if (refracted_water <= 0)
		{
			refraction_branch = float3(1, 0, 0);
			// Nothing in front of it to bend the ray, so fill holes by taking
			// the unperturbed sample:
			refraction_uv = ScreenCoord.xy;
		}
		else if (camera_below_water)
		{
			refraction_branch = float3(0, 1, 0);

			// Perturbation according to the water the first sample is behind:
			refraction_uv = ScreenCoord.xy + surface.N.xz * bump_strength * saturate(1 - exp(-refracted_water));
		}
		else
		{
			// Follow the ray the surface actually bends, rather than sliding
			// the lookup sideways by the wave's tilt. Level water refracts as
			// hard as a crest does: the bend is towards the vertical, so the
			// ray lands NEARER the eye than the straight line does, and
			// whatever stands there was drawn lower on screen.
			//
			// A tilt-driven offset is zero wherever the surface happens to be
			// level, so those pixels read the copy at their own coordinate and
			// show the scene exactly where it really is. That is what puts an
			// unrefracted image of a submerged object beside its refracted
			// one: the crest bends its share, the level water beside it does
			// not bend at all.
			//
			refraction_branch = float3(0, 0, 1);

			// Snell's law bounds the direction. Nothing entering from air runs
			// more than 48.6 degrees off vertical, so the descent rate never
			// falls below 0.661 and the path down to a given depth is at most
			// 1.51 times that depth.
			const float3 refracted_dir = refract(
				-surface.V, surface.N, 1.0 / WATER_REFRACTIVE_INDEX);
			const float descent = max(0.661, -refracted_dir.y);

			// Where there is no bed the depth buffer holds the far plane, and
			// the column measured against it runs to kilometres. Cut the ray at
			// the distance past which the water has extinguished everything
			// anyway: beyond it the sample can contribute nothing, and carrying
			// on only lands the lookup somewhere arbitrary. The range is a path
			// length and the column is a vertical drop, hence the descent rate
			// between them.
			const WaterVolumetrics medium = MakeWaterVolumetrics(1);
			float reach = medium.IsActive()
				? min(refracted_water, medium.VisibleRange() * descent)
				: refracted_water;

			// refracted_water is measured down from the still plane, so this
			// lands on the target's depth rather than on the target: the wave
			// the ray left from stands above or below that plane by its own
			// displacement. The error is one wave height against the whole
			// column, and it shifts the sample rather than dropping it.
			float3 refracted_hit =
				surface.P + refracted_dir * (reach / descent);
			float4 refracted_clip =
				mul(camera.view_projection, float4(refracted_hit, 1));
			float2 refracted_screen = clipspace_to_uv(
				refracted_clip.xy / refracted_clip.w);

			// One refinement, and it is what decides whether an object standing
			// in the water is refracted or duplicated.
			//
			// The depth above belongs to whatever lies under THIS pixel, which
			// is not what the bent ray arrives at. Aimed at the sea bed, the ray
			// sails straight past anything standing proud of it and its image is
			// flung far up the screen - so the object is drawn a second time,
			// somewhere it plainly is not, instead of being displaced a little
			// from where it is. Re-aim at the depth of whatever is drawn where
			// the first guess landed: when that really is the first thing in the
			// ray's way, the copy is holding exactly the right answer at exactly
			// that coordinate.
			//
			// **Only ever shortened.** A step that may lengthen can chase its
			// own target and settle somewhere unrelated, which is how a
			// fixed-point solve here ends up mirroring the view rather than
			// refracting it. Shortening alone is monotone and terminates.
			//
			// Measured as a depth below the plane, so the step divides by the
			// descent rate, which Snell's law bounds below at 0.661 - never by a
			// rate of recession from the eye, which collapses at a grazing view.
			// Only where the first guess landed on the screen at all. Past the
			// edge the depth tap clamps to the border and returns a pixel from
			// somewhere unrelated, so refining against it replaces a bounded
			// estimate with an arbitrary one - and the column the fog is later
			// charged for comes from the same number.
			const float2 guess_border =
				min(refracted_screen, 1 - refracted_screen);

			[branch]
			if (min(guess_border.x, guess_border.y) > 0 && refracted_clip.w > 0)
			{
				// Point sampled. The neighbourhood maximum the column measure
				// above takes would stop the ray at the nearest thing within
				// eight screen pixels of where it looked, which pulls it up
				// short all the way around every silhouette.
				const float guess_depth = texture_depth.SampleLevel(
					sampler_point_clamp, refracted_screen, 0);
				const float3 guess_position =
					reconstruct_position(refracted_screen, guess_depth);

				const float guess_water =
					-dot(float4(guess_position, 1), water_plane);

				// A target above the surface is not on this ray at all. The ray
				// descends and stays inside Snell's cone, so nothing above the
				// water can ever be in its way - a guess landing there means the
				// first aim passed BEHIND something standing out of the water,
				// not that it arrived at it.
				//
				// Taking it as a shorter target collapses the reach to nothing,
				// which stops the ray dead at the surface and reads the copy at
				// the fragment's own coordinate. Beside an object that is where
				// the object is drawn, so it rims itself in its own colour - and
				// the endpoint, now the surface, is nearer than anything the
				// probe can find, so the sample is approved on the way past.
				reach = guess_water > 0 ? min(guess_water, reach) : reach;

				refracted_hit = surface.P + refracted_dir * (reach / descent);
				refracted_clip =
					mul(camera.view_projection, float4(refracted_hit, 1));
				refracted_screen = clipspace_to_uv(
					refracted_clip.xy / refracted_clip.w);
			}

			refraction_uv = refracted_screen;

			// A projected point gives a COORDINATE, not what is drawn at it.
			// The pixel there holds whatever stood nearest along its OWN ray,
			// which is the point this ray reached only while nothing is in
			// front of it. Compare the two depths and believe the sample only
			// while they agree: an occluder reads far nearer than the point
			// the ray actually travelled to.
			//
			// Tested in view depth, not against the water plane. The plane
			// only catches occluders standing ABOVE the surface, and the ones
			// that matter most are under it - a nearer submerged object, and
			// the underwater half of the screen when the eye straddles the
			// waterline, both of which the plane test waves through.
			// Point sampled, for a different question than the column measure
			// above asks. That one asks how much water stands in front of the
			// COLOUR, where the copy's coarseness means a solid two texels away
			// still bleeds into the sample, so the nearest thing nearby is the
			// right answer. This one asks whether THIS coordinate holds the
			// point the ray reached - and answering that with the nearest thing
			// anywhere nearby condemns a ring of perfectly good pixels around
			// every silhouette, which the fill then paints.
			const float probe_depth = texture_depth.SampleLevel(
				sampler_point_clamp, refracted_screen, 0);
			const float intrusion =
				refracted_clip.w - compute_lineardepth(probe_depth);

			// Scaled by distance, because the copy's depth precision and the
			// world footprint of one of its texels both are. Floored so the
			// two smoothstep edges stay ordered - a hit landing on the eye
			// plane would otherwise be handed a zero-width range, which is
			// undefined and returns a value that poisons everything downstream.
			const float tolerance = 0.02 * max(1, refracted_clip.w);
			const float agrees =
				1 - smoothstep(tolerance, tolerance * 5, intrusion);

			// Was the ray already blocked before it got there?
			//
			// The test above cannot answer that on its own. The refinement
			// aims at the depth of whatever is drawn where its guess landed,
			// so wherever that guess fell on an occluder the endpoint is
			// DEFINED from that occluder - and comparing the two then compares
			// it against itself and always agrees. A ray passing BESIDE an
			// object settles onto it exactly as one that runs into it does,
			// which is what rims an object in its own colour.
			//
			// Self-consistency is not intersection. Walk a few points along
			// the way and require each to stand in FRONT of whatever is drawn
			// where it projects: a ray that really reaches its endpoint is
			// unobstructed the whole way, and one that merely ends somewhere
			// agreeable behind an occluder is not.
			float unblocked = 1;

			[unroll]
			for (uint walk = 1; walk <= 3; ++walk)
			{
				const float3 waypoint = surface.P
					+ refracted_dir * (reach / descent) * (walk * 0.25);
				const float4 waypoint_clip =
					mul(camera.view_projection, float4(waypoint, 1));
				const float2 waypoint_screen = clipspace_to_uv(
					waypoint_clip.xy / waypoint_clip.w);

				const float waypoint_scene = compute_lineardepth(
					texture_depth.SampleLevel(
						sampler_point_clamp, waypoint_screen, 0));

				unblocked = min(unblocked, 1 - smoothstep(
					tolerance,
					tolerance * 5,
					waypoint_clip.w - waypoint_scene));
			}

			// A hit that projects from BEHIND the eye has no coordinate at
			// all: w turns negative and the divide wraps it to the opposite
			// side of the screen, where it looks like a perfectly ordinary
			// coordinate holding something else entirely.
			const float infront = refracted_clip.w > 0 ? 1 : 0;

			// The copy simply stops at its border. Faded rather than cut, so
			// no seam runs along the edge of the screen.
			const float2 border = min(refracted_screen, 1 - refracted_screen);
			const float inside = smoothstep(0, 0.02, min(border.x, border.y));

			// With no medium there is nothing to fill a blind pixel WITH, so
			// keep the sample rather than take it out and leave black behind.
			refraction_confidence = medium.IsActive()
				? agrees * unblocked * inside * infront
				: 1;

			// Both the endpoint test and the walk are the same complaint - the
			// ray never got a clear sight of what it aimed at - so they share
			// the red channel.
			refraction_blind_reason = float3(
				1 - agrees * unblocked, 1 - infront, 1 - inside);


			// The column a blind pixel is charged for is the one its ray really
			// crossed, not the whole distance the water can be seen through.
			// Those are the same number only in water deep enough to hide its
			// own bed; over a bed a few metres down they are not close, and
			// charging the visible range there paints a patch of open-ocean
			// blue into shallow water.
			refraction_blind_depth = reach;
		}
		// Filtered, and it must stay filtered.
		//
		// The copy is a quarter of screen resolution, so at a water pixel
		// touching an object every texel within reach already mixes that object
		// with the water behind it. There is no clean tap nearby to prefer, and
		// point sampling only trades a smoothed mixture for a hard-edged one -
		// which reads as a rim drawn around the object rather than as the
		// antialiasing it is. The band is the copy's resolution and is closed by
		// raising it, not by choosing differently between contaminated samples.
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

		// Where the copy answered, this is the crest leg it did not charge for.
		// Where it did not, the sample has already been taken out above and the
		// whole visible column stands in its place, so what survives is the
		// in-scatter alone - the colour of water with nothing behind it.
		// Charging a blind pixel for the crest leg only would leave it near
		// black, with a hard edge wherever the belief ran out.
		const float fog_depth =
			lerp(refraction_blind_depth, uncounted_depth, refraction_confidence);

		// Taken for every fragment rather than only where there is a leg to
		// charge for: a blind pixel's entire colour comes from here, and a zero
		// path is a no-op anyway - transmittance 1, in-scatter nothing.
		[branch]
		if (!camera_below_water)
		{
			const WaterVolumetrics behindMedium = MakeWaterVolumetrics(1);

			const WaterFog waterBehind = MakeWaterFog(
				behindMedium,
				SubmergedViewPath(fog_depth, V),
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

			// What stands behind the water where the copy could not say.
			//
			// Not nothing: over a bed only metres down the water is nowhere
			// near thick enough to hide it, so a pixel handed an empty
			// background comes back almost black however carefully it is
			// fogged. What is missing is the bed's own colour at a point the
			// eye cannot see, and the one honest thing to say about it is that
			// it is a dull surface under the daylight that got down there -
			// which is the irradiance this fog is already built from, sun and
			// sky together, so it is lit exactly as everything else at this
			// depth is rather than by a dimmer figure of its own.
			//
			// Deliberately not sampled from anywhere on screen: every
			// coordinate that holds scene holds a second copy of that scene,
			// which is the whole family of bugs this lookup keeps producing.
			//
			// Self-limiting. It is weighed against the column actually crossed,
			// so it shows through in shallow water, where it is nearly right,
			// and is extinguished in deep water, where it is a guess.
			// Lambertian, so the daylight arriving has to be turned into the
			// radiance leaving before it is a colour at all: albedo over pi.
			// Handed straight over as irradiance it comes back about three
			// times too bright, which on a dull surface under full sun is the
			// difference between wet sand and white paper.
			//
			// Dimmed by the column standing over it as well. The fog below
			// attenuates the leg back up to the eye and this is the other one:
			// light does not arrive ten metres down as strong as it left the
			// surface, and without this the invented bed brightens with depth
			// exactly where the real one is fading out.
			const half3 unseen_bed = (half3)(
				waterBehind.downwelling
				* (UNSEEN_BED_ALBEDO / PI)
				* exp(-SubmergedViewPath(fog_depth, GetSunDirection())
					* behindMedium.sigmaT));

			surface.refraction.rgb = (half3)(
				lerp(unseen_bed, surface.refraction.rgb, refraction_confidence)
				* waterBehind.transmittance
				+ waterBehind.inscatter);

#if OCEAN_REFRACTION_DEBUG == 2
			return float4(waterBehind.inscatter, 1);
#endif // OCEAN_REFRACTION_DEBUG == 2
		}

#if OCEAN_REFRACTION_DEBUG == 1
		return float4(
			lerp(refraction_blind_reason, (float3)1, refraction_confidence), 1);
#elif OCEAN_REFRACTION_DEBUG == 3
		return float4(surface.refraction.rgb, 1);
#elif OCEAN_REFRACTION_DEBUG == 4
		return float4(refraction_branch, 1);
#endif // OCEAN_REFRACTION_DEBUG
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

#if OCEAN_REFRACTION_DEBUG == 5
	// Flattens everything the refraction lookup produced and leaves the rest of
	// the surface's shading alone. Whatever is still structured after this is
	// drawn by something that is not the lookup.
	surface.refraction.rgb = half3(1, 0, 1);
#endif // OCEAN_REFRACTION_DEBUG == 5

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

