#include "bounding_volume_hierarchy.h"
#include "disable_all_warnings.h"
#include "draw.h"
#include "image.h"
#include "ray_tracing.h"
#include "ray_differentials.h"
#include "screen.h"
#include "trackball.h"
#include "window.h"
#include "shadow.h"

// Disable compiler warnings in third-party code (which we cannot change).
DISABLE_WARNINGS_PUSH()
#include <glm/gtc/type_ptr.hpp>
#include <glm/vec2.hpp>
#include <glm/gtx/component_wise.hpp>
#include <glm/vec4.hpp>
#include <imgui.h>
DISABLE_WARNINGS_POP()
#include <chrono>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <optional>
#include <random>
#include <string>
#include <type_traits>
#ifdef USE_OPENMP
#include <omp.h>
#endif

// This is the main application. The code in here does not need to be modified.
constexpr glm::ivec2 windowResolution{ 800, 800 };
const std::filesystem::path dataPath{ DATA_DIR };
const std::filesystem::path outputPath{ OUTPUT_DIR };


struct ProgressIndicator {
	int currentProgress = 0;

	void setProgress(int numPixels) {
		float percent = (float) numPixels / (float) (windowResolution.x * windowResolution.y);
		int progress = (int)floor(percent * 10);
		progress *= 10;

		if(progress > currentProgress) {
			currentProgress = progress;
			std::cout << "Ray Tracing Progress: " << currentProgress <<"%\n";
		}
	}
};

// Textures settings - used for debugging only!
TextureFiltering textureFiltering{TextureFiltering::NearestNeighbor};
OutOfBoundsRule outOfBoundsRuleX{OutOfBoundsRule::Border}; // Border, Clamp, Repeat
OutOfBoundsRule outOfBoundsRuleY{OutOfBoundsRule::Border};
glm::vec3 textureBorderColor(0);
bool useTextures = false;

bool useBVH = false;

// the variable controls weather or not another image should be rendered in raytracing mode.
// used to put an obtained image on halt
bool rayTracing_toBeRendered = true;


enum class ViewMode {
	Rasterization = 0,
	RayTracing = 1,
	Textures = 2,
	MultipleRayDebug = 3
};


// used for debugging the textures
static glm::vec3 getFinalColorNoRayTracingJustTextures(Scene &scene, const BoundingVolumeHierarchy &bvh, Ray ray) {
	HitInfo hitInfo;



	if(bvh.intersect(ray, hitInfo, useBVH)) {

		// does not reflect rays so we only transfer the ray differentials
		transfer_ray_differentials(ray, hitInfo.normal);

		Material& mat = hitInfo.getMaterial(scene);
		if(mat.kdTexture) {
			Image& texture = mat.kdTexture.value();

			texture.setBorderColor(textureBorderColor);
			texture.setOutOfBoundsRuleX(outOfBoundsRuleX);
			texture.setOutOfBoundsRuleY(outOfBoundsRuleY);
			texture.setTextureFilteringMethod(textureFiltering);

			// verticies
			glm::vec2 textureCoordinate = hitInfo.texCoord;

			//return glm::vec3(1);

			// return the color corresponding to the texture
			float lod = computeLevelOfDetails(ray, hitInfo);
			return texture.getPixel(hitInfo.texCoord, lod);
		}

		return glm::vec3(1);
	}

	return glm::vec3(0);
}


glm::vec3 calcColor(Lighting light, Material material) {
	// difuse light
	glm::vec3 diffuse = material.kd * light.color * light.intensity * light.cosLightSurfaceAngle;
	//specular light
	glm::vec3 spec = glm::vec3(0);
	if (material.shininess > 0) {
		spec = light.color * material.ks * std::pow(light.cosLightSpecAngle, material.shininess);
	}
	return diffuse + spec;
}

static int max_reflection_level = 5;
static int sphere_light_ray_count = 10;
static int plane_light_1D_ray_count = 3;
static int glossy_ray_count = 10;
static float refraction_factor = 0.8;
// NOTE(Mathijs): separate function to make recursion easier (could also be done with lambda + std::function).
static glm::vec3 getFinalColor(Scene& scene, const BoundingVolumeHierarchy& bvh, Ray ray, int level=0)
{
	HitInfo hitInfo;

	/*For every light calulate the addition from the reflected rays,
	then add all of them*/
	if (bvh.intersect(ray, hitInfo, useBVH)) {

		tranfer_and_reflect_ray_differentials(ray, hitInfo);

		glm::vec3 color(0);

		glm::vec3 reflect = glm::reflect(glm::normalize(ray.direction), glm::normalize(hitInfo.normal));

		// if textures to be used and the material has texture, use it
		// also, avoid copying the texture over and over again
		// to achieve this, a new material is created, that does not carry over the image, only the color in that particular point

		Material& originalMaterial = hitInfo.getMaterial(scene);

		Material matForRendering;
		matForRendering.kd = originalMaterial.kd;
		matForRendering.ks = originalMaterial.ks;
		matForRendering.shininess = originalMaterial.shininess;
		matForRendering.transparency = originalMaterial.transparency;

		if (useTextures && hitInfo.is_triangle && originalMaterial.kdTexture)
		{
			Image& texture = originalMaterial.kdTexture.value();

			texture.setBorderColor(textureBorderColor);
			texture.setOutOfBoundsRuleX(outOfBoundsRuleX);
			texture.setOutOfBoundsRuleY(outOfBoundsRuleY);
			texture.setTextureFilteringMethod(textureFiltering);

			// verticies
			glm::vec2 textureCoordinate = hitInfo.texCoord;

			

			float lod = computeLevelOfDetails(ray, hitInfo); // the level of details, used in mipmapping, based on ray differentials
			matForRendering.kd = texture.getPixel(hitInfo.texCoord, lod);
		}

		// Loop over all the lights that are not in shadow
		for (const Lighting& light : getPointLights(hitInfo, reflect, scene, bvh)) {
			color += calcColor(light, matForRendering);
		}
		for (const Lighting& light : getSpherelights(hitInfo, reflect, scene, bvh, sphere_light_ray_count)) {
			color += calcColor(light, matForRendering);
		}
		for (const Lighting& light : getSpotLichts(hitInfo, reflect, scene, bvh)) {
			color += calcColor(light, matForRendering);
		}
		for (const Lighting& light : getPlaneLights(hitInfo, reflect, scene, bvh, plane_light_1D_ray_count)) {
			color += calcColor(light, matForRendering);
		}

		if (level >= max_reflection_level) {
			return color;
		}

		if (matForRendering.transparency == 1.0f)
		{ // not tranparent
			// calculate reflected light
			if (matForRendering.ks.x > 0 || matForRendering.ks.y > 0 || matForRendering.ks.z > 0)
			{

				glm::vec3 reflectColor = glm::vec3(0);
				// the ( + 0.01f * reflect ) is to prevent a surface reflecting itself
				Ray refRay = { hitInfo.hitPoint + 0.01f * reflect, reflect };
				// Calculate the color that is reflected
				reflectColor += matForRendering.ks * getFinalColor(scene, bvh, refRay, level + 1);

				// prevent divide by 0 and allow for prefect reflections
				if (matForRendering.shininess != 0)
				{

					//hitInfo.material.shininess = shininess;

					// a vector that is not in line whith reflect
					glm::vec3 notr = reflect;
					if (reflect.x != 0) {
						notr.y = -reflect.x;
						notr.x = reflect.y;
					}
					else {
						notr.y = -reflect.z;
						notr.z = reflect.y;
					}
					// right angle between reflection ray and pr1 and pr2
					glm::vec3 pr1 = glm::cross(reflect, notr);
					glm::vec3 pr2 = glm::cross(reflect, pr1);

					// how mutch should the rayes diviate (cos(angle)^shininess > 0.5)
					float d = std::pow(0.5f, -1 / (float)matForRendering.shininess) * std::sqrt(1 - std::pow(0.5, 2 / (float)matForRendering.shininess));

					//srand(6437376437);
					for (int i = 1; i < glossy_ray_count; i++) {
						// find random direction
						glm::vec3 shineDir;
						float a, b;
						int loopcount = 0; // prevent infinent loop
						do {
							do {
								a = rand() / (float)RAND_MAX;
								b = rand() / (float)RAND_MAX;
							} while (a == 0 && b == 0 && a * a + b * b < 1);
							a = (2 * a - 1) * d;
							b = (2 * b - 1) * d;
							shineDir = glm::normalize(reflect + a * pr1 + b * pr2);
							loopcount++;
						} while (glm::dot(shineDir, hitInfo.normal) <= 0 && loopcount<glossy_ray_count/4);
						if (glm::dot(shineDir, hitInfo.normal) > 0) {
							// cast the ray
							Ray shineRay = { hitInfo.hitPoint + 0.01f * shineDir, shineDir };
							glm::vec3 ctmp = getFinalColor(scene, bvh, shineRay, level + 1) * std::max(std::pow(glm::dot(reflect, shineDir), matForRendering.shininess), 0.0f);
							//std::cout << "r_r: " << ctmp.r << std::endl;
							reflectColor += ctmp;
						}
					}
					color += matForRendering.ks * reflectColor / (float)glossy_ray_count;
				}
				else {
					color += matForRendering.ks * reflectColor;
				}
			}
		}
		else { // transparent
			
			// sources
			// https://www.cs.drexel.edu/~david/Classes/Papers/p343-whitted.pdf
			// https://en.wikipedia.org/wiki/Snell%27s_law#Vector_form
			// https://www.youtube.com/watch?v=iKNSPETJNgo&list=PLujxSBD-JXgnGmsn7gEyN28P1DnRZG7qi&index=5
			// https://stackoverflow.com/questions/29758545/how-to-find-refraction-vector-from-incoming-vector-and-surface-normal : cosI = abs(dot(normal, incident))


			// calculate the refraction direction
			glm::vec3 l = glm::normalize(ray.direction);
			glm::vec3 n = glm::normalize(hitInfo.normal);

			float r = refraction_factor; // n1/n2 : indexes of refraction
			float c = std::abs(glm::dot(l, n));

			    // vector form of snell's law
			glm::vec3 refract = r*l + ( r*c - std::sqrt(1 - r*r*(1 - c*c )) )*n; // std::sqrt might return NaN
			refract = glm::normalize(refract);
			

			// calculate how how much light should reflect and how much should refract
			float &R0 = matForRendering.transparency; //R0static; // probebility of reflection at 0 rad/deg to normal

			// simplified fersnel equasion
			float reflectionChance = R0 + (1-R0)*(std::pow(1-c,5));
			float refractionChance = 1 - reflectionChance;

			// get the color the reflected and refrected rays see
			color += reflectionChance * getFinalColor(scene, bvh, { hitInfo.hitPoint + 0.01f * reflect, reflect }, level + 1);
			if (r * r * (1 - c * c) <= 1.0f) { // check if total internal reflection occures
				color += refractionChance * getFinalColor(scene, bvh, { hitInfo.hitPoint + 0.01f * refract, refract }, level + 1);
			}
		}
		drawRay(ray, color);
		return color;
	}
	// if the ray did not hit anything
	else { 
		// Draw a red debug ray if the ray missed.
		drawRay(ray, glm::vec3(1.0f, 0.0f, 0.0f));
		// Set the color of the pixel to black if the ray misses.
		return glm::vec3(0.0f);
	}
}


static void setOpenGLMatrices(const Trackball& camera);
static void renderOpenGL(const Scene& scene, const Trackball& camera, int selectedLight);

// Subdivide the pixel according to the given sampleSize
// Returns a vector of the centers of every one of the subdivisions
std::vector<glm::vec2> getPixelRays(glm::vec2 pixelCenter, int sampleSize) {

	float offsetX = (1.0f / windowResolution.x) * (1.0f / (glm::sqrt(sampleSize) * 2));
	float offsetY = (1.0f / windowResolution.y) * (1.0f / (glm::sqrt(sampleSize) * 2));

	std::vector<glm::vec2> origins;
	std::array<glm::vec2, 4> quadrantSigns = { glm::vec2(-1.0f,1.0f), //top-left - +
											glm::vec2(1.0f,1.0f),		//top-right	+ +
											glm::vec2(-1.0f,-1.0f),	//bottom-left - -
											glm::vec2(1.0f,-1.0f) };	//bottom-right + -
	
	// max moves needed to move the pixelCenter to the center of the subdivision
	// x4 --> 1 move by offset for furthest
	// x16 --> 5 moves by offset for furthest
	// x64 --> 7 moves by offset for furthest
	int moves = glm::sqrt(sampleSize) - 1; 

	for (int i = 0; i < 4; i++) { // move over quadrants
		for (int x = 1; x <= moves; x = x + 2) { // move horizontally in pixel
			for (int y = 1; y <= moves; y = y + 2) { // move vertically in pixel
				origins.push_back(glm::vec2(pixelCenter.x + (offsetX * quadrantSigns[i].x * x), 
											pixelCenter.y + (offsetY * quadrantSigns[i].y * y)));
			}
		}
	}
	return origins; // the origins generated by subdividing the original pixel
}


// This is the main rendering function. You are free to change this function in any way (including the function signature).
// if texture debugging is set to true, the method will call the getFinalColorNoRayTracingJustTextures, instead of getFinalColor.
static void renderRayTracing(Scene& scene, const Trackball& camera, const BoundingVolumeHierarchy& bvh, Screen& screen,
								bool textureDebugging = false, bool anti_aliasing = false, bool multipleRays = false, int sampleSize = 4) {
	ProgressIndicator indicator;
	int pixelCount = 0;
	#ifdef USE_OPENMP
	#pragma omp parallel for
	#endif
	for (int y = 0; y < windowResolution.y; y++) {
		for (int x = 0; x != windowResolution.x; x++) {
			// NOTE: (-1, -1) at the bottom left of the screen, (+1, +1) at the top right of the screen.
			const glm::vec2 normalizedPixelPos{
				float(x) / windowResolution.x * 2.0f - 1.0f,
				float(y) / windowResolution.y * 2.0f - 1.0f
			};
			const Ray cameraRay = camera.generateRay(normalizedPixelPos);
			if (textureDebugging)
				screen.setPixel(x, y, getFinalColorNoRayTracingJustTextures(scene, bvh, cameraRay));
			else {
				if (anti_aliasing) {
					//get a quarter pixel's worth of offset depending on the window resolution
					float offsetX = 1.0f / windowResolution.x * 0.25f;
					float offsetY = 1.0f / windowResolution.y * 0.25f;

					std::array <glm::vec2, 4> offsets; // calculate quadrant offsets for each pixel
					offsets[0] = (glm::vec2(normalizedPixelPos.x - offsetX, normalizedPixelPos.y + offsetY)); // - + top-left quadrant
					offsets[1] = (glm::vec2(normalizedPixelPos.x + offsetX, normalizedPixelPos.y + offsetY)); // + + top-right quadrant
					offsets[2] = (glm::vec2(normalizedPixelPos.x - offsetX, normalizedPixelPos.y - offsetY)); // - - bottom-left quadrant
					offsets[3] = (glm::vec2(normalizedPixelPos.x + offsetX, normalizedPixelPos.y - offsetY)); // + - bottom-right quadrant

					glm::vec3 avgColor(0);
					for (int i = 0; i < 4; i++) { //get color values for each pixel offset
						Ray ray = camera.generateRay(offsets[i]);
						avgColor += getFinalColor(scene, bvh, ray);
					}
					avgColor *= 0.25; // divide by the number of quadrants to get the average color for that pixel
					screen.setPixel(x, y, avgColor);
				}
				else if (multipleRays) {
					std::vector<glm::vec2> rayOrigins = getPixelRays(normalizedPixelPos, sampleSize);
					glm::vec3 avgColor(0);
					for (auto& rayOrigin : rayOrigins) {
						Ray ray = camera.generateRay(rayOrigin);
						avgColor += getFinalColor(scene, bvh, ray);
					}
					avgColor = avgColor * (float)(1.0f / sampleSize);
					screen.setPixel(x, y, avgColor);
				}
				else {
					screen.setPixel(x, y, getFinalColor(scene, bvh, cameraRay)); //anti-aliasing off 
				}
			}
		
			pixelCount++;
			indicator.setProgress(pixelCount);
		}
	}

	// applies bloom filter or gamma effect
	screen.postprocessImage();

}

int main(int argc, char** argv)
{
    Trackball::printHelp();
    std::cout << "\n Press the [R] key on your keyboard to create a ray towards the mouse cursor" << std::endl
			  << "\n Press the [A] key on your keyboard to switch to OpenGL rendering mode" << std::endl
			  << "\n Prsss the [B] key on your keyboard to switch to RayTracing mode" << std::endl
			  << "\n Press the [C] key on your keyboard to switch to RayTracint Only Texture mode" << std::endl
              << std::endl;

    Window window { "Final Project - Part 2", windowResolution, OpenGLVersion::GL2 };
    Screen screen { windowResolution };
    Trackball camera { &window, glm::radians(50.0f), 3.0f };
    camera.setCamera(glm::vec3(0.0f, 0.0f, 0.0f), glm::radians(glm::vec3(20.0f, 20.0f, 0.0f)), 3.0f);
	ViewMode viewMode{ViewMode::Rasterization};

	SceneType sceneType { SceneType::SingleTriangle };
    std::optional<Ray> optDebugRay;
    Scene scene = loadScene(sceneType, dataPath);
    BoundingVolumeHierarchy bvh { &scene };

    int bvhDebugLevel = 0;
    bool bvhShowLeafNodes {false}; // controls whether or not the leaf nodes should be highlighted in red
    bool debugBVH { false };

	// Bloom settings
	bool liveDebug = false;
	FilteringOption bloom_filtering_option = FilteringOption::None;
	Kernel kernel = Kernel::BoxKernel;
	int filterSize = 5; // 1 to 30
	bool gammaCorrection = false;
	float gammaValue = 2.2; // 1 to 5
	float exposure = 0.5;	// 0 to 1
	float sigma = 2.0f;
	int kernel_num_repetitions = 1;

	// Anti-Alising toggle
	bool anti_aliasing = false;
	
	// Multiple ray casting settings
	bool multipleRays = false;
	int sampleSize = 4;
	std::vector<Ray> debugRays;

    window.registerKeyCallback([&](int key, int /* scancode */, int action, int /* mods */) {
        if (action == GLFW_PRESS) {
            switch (key) {
            case GLFW_KEY_R: {
                // Shoot a ray. Produce a ray from camera to the far plane.
                const auto tmp = window.getNormalizedCursorPos();
                optDebugRay = camera.generateRay(tmp * 2.0f - 1.0f);
				if (multipleRays) { //if multiple rays is turned on
					debugRays.clear();
					std::vector<glm::vec2> origins = getPixelRays(tmp * 2.0f - 1.0f, sampleSize);
					for (auto& origin : origins) {
						debugRays.push_back(camera.generateRay(origin));
					}
					viewMode = ViewMode::MultipleRayDebug;
				}
				else {
					viewMode = ViewMode::Rasterization;
				}
            } break;
            case GLFW_KEY_ESCAPE: {
                window.close();
            } break;
			case GLFW_KEY_A: {
				viewMode = ViewMode::Rasterization;
				break;
			}
			case GLFW_KEY_B: {
				viewMode = ViewMode::RayTracing;
				rayTracing_toBeRendered = true;
				break;
			}
			case GLFW_KEY_C: {
				viewMode = ViewMode::Textures;
				rayTracing_toBeRendered = true;
				break;
			}
			
            };
        }
    });

    int selectedLight { 0 };
    while (!window.shouldClose()) {
        window.updateInput();

        // === Setup the UI ===
        ImGui::Begin("Final Project - Part 2");
        {
			constexpr std::array items{"SingleTriangle", "Bookshelf", "Cube", "Cornell Box (with mirror)", "Cornell Box (spherical light and mirror)", "Cornell Box (plane light and mirror)", "Monkey",
				"Teapot", "Dragon", "Spheres", "Chess", "Custom", "AndreasScene", "CatalinScene", "MikeScene", "MikeScene2"};
			if (ImGui::Combo("Scenes", reinterpret_cast<int *>(&sceneType), items.data(), int(items.size())))
			{
				optDebugRay.reset();
                scene = loadScene(sceneType, dataPath);
                bvh = BoundingVolumeHierarchy(&scene);
                if (optDebugRay) {
                    HitInfo dummy {};
                    bvh.intersect(*optDebugRay, dummy, useBVH);
                }
			}
		}
        {
            constexpr std::array items { "Rasterization", "Ray Traced", "Textures" , "Multiple Ray Debug"};
            ImGui::Combo("View mode", reinterpret_cast<int*>(&viewMode), items.data(), int(items.size()));
        }

		ImGui::Checkbox("Use BVH", &useBVH);

        if (ImGui::Button("Render to file")) {
            {
                using clock = std::chrono::high_resolution_clock;
                const auto start = clock::now();
                renderRayTracing(scene, camera, bvh, screen, false, anti_aliasing, multipleRays, sampleSize);
                const auto end = clock::now();
                std::cout << "Time to render image: " << std::chrono::duration<float, std::milli>(end - start).count() << " milliseconds" << std::endl;
            }
            screen.writeBitmapToFile(outputPath / "render.bmp");
        }

		if (ImGui::TreeNode("Ray count settings")) {
			ImGui::SliderInt("Max reflections", &max_reflection_level, 1, 15);

			ImGui::SliderInt("Sphere light rays", &sphere_light_ray_count, 1, 40);
			ImGui::SliderInt("Plane light rays(1D)", &plane_light_1D_ray_count, 1, 6);
			
			ImGui::SliderInt("Glossy reflection ray count", &glossy_ray_count, 1, 40);

			ImGui::SliderFloat("Refraction", &refraction_factor, 0, 2);
			ImGui::TreePop();
		}
		
		// Gamma correction 
		ImGui::Checkbox("Gamma Correction", &gammaCorrection);
		if(gammaCorrection)
			ImGui::SliderFloat("Gamme Value", &gammaValue, 1, 5);


		screen.enableGammaCorrection(gammaCorrection);
		screen.setGammaValue(gammaValue);

		//Anti-Aliasing
		ImGui::Checkbox("Anti-Aliasing", &anti_aliasing);

		//Mutliple Ray casting per pixel
		ImGui::Checkbox("Cast multiple rays", &multipleRays);
		if (multipleRays) {
			ImGui::SliderInt("Sample size", &sampleSize, 4, 64);
			if (sampleSize < 16) {
				sampleSize = 4;
			}
			else if(sampleSize >= 16 && sampleSize < 64){
				sampleSize = 16;
			}
			else {
				sampleSize = 64;
			}
		}

		if(ImGui::TreeNode("Bloom Filtering Settings")) {
			
			ImGui::Checkbox("Live Debug", &liveDebug);

			constexpr std::array filtering_modes{"None", "Bloom", "Bloom With Reinhard Tone Mapping", "Bloom With Exposure Based Tone Mapping", "Show Only Bright Areas (Debugging)", "Show Only Bright Areas and Apply Kernel (Debugging)"};
			ImGui::Combo("Bloom Effect Options", reinterpret_cast<int *>(&bloom_filtering_option), filtering_modes.data(), int(filtering_modes.size()));

			constexpr std::array kernel_modes{"Box Kernel", "Gaussian Kernel"};
			ImGui::Combo("Kernel", reinterpret_cast<int *>(&kernel), kernel_modes.data(), int(kernel_modes.size()));

			ImGui::SliderInt("FilterSize", &filterSize, 1, 300);
			ImGui::SliderInt("No. Kernel Repetitions", &kernel_num_repetitions, 1, 10);
			if(kernel == Kernel::GaussianKernel) 
				ImGui::SliderFloat("Sigma", &sigma, 0.1f, 5.0f);


			if(bloom_filtering_option == FilteringOption::BloomWithExposureHdr)
				ImGui::SliderFloat("Exposure Level", &exposure, 0, 1);

			ImGui::TreePop();

			// setup the bloom filtering settings
			screen.setBloomFilterLive(liveDebug);
			screen.setBloomFilter(bloom_filtering_option);
			screen.setKernel(kernel);
			screen.setFilterSize(filterSize);
			screen.setKernelNumRepetitions(kernel_num_repetitions);
			screen.setExposure(exposure);
			screen.setSigma(sigma);
		}

		if (ImGui::TreeNode("Texture Settings"))
		{


			ImGui::Checkbox("Use Textures", &useTextures);


			constexpr std::array filtering_modes{"Nearest Neighbor", "Bilinear", "Nearest Level Mipmapping", "Bilinear Mipmapping", "Trilinear Mipmapping"};
			ImGui::Combo("Filtering Mode", reinterpret_cast<int *>(&textureFiltering), filtering_modes.data(), int(filtering_modes.size()));

			constexpr std::array oob_x{"Border", "Clamping", "Repeat"};
			ImGui::Combo("Out Of Bounds X", reinterpret_cast<int *>(&outOfBoundsRuleX), oob_x.data(), int(oob_x.size()));

			constexpr std::array oob_y{"Border", "Clamping", "Repeat"};
			ImGui::Combo("Out Of Bounds Y", reinterpret_cast<int *>(&outOfBoundsRuleY), oob_y.data(), int(oob_y.size()));

			ImGui::ColorEdit3("Border color", glm::value_ptr(textureBorderColor));

			ImGui::TreePop();
		}

		ImGui::Spacing();
        ImGui::Separator();
        ImGui::Text("Debugging");
        if (viewMode == ViewMode::Rasterization) {
            ImGui::Checkbox("Draw BVH", &debugBVH);
            if (debugBVH) {
                ImGui::SliderInt("BVH Level", &bvhDebugLevel, 0, bvh.numLevels() - 1);
                ImGui::Checkbox("Show BVH leaf nodes", &bvhShowLeafNodes);
            }

	
            
        } 

        ImGui::Spacing();
        ImGui::Separator();
        ImGui::Text("Lights");
        if (!scene.pointLights.empty() || !scene.sphericalLight.empty() || !scene.spotLight.empty() || !scene.planeLight.empty()) {
            {
                std::vector<std::string> options;
                for (size_t i = 0; i < scene.pointLights.size(); i++) {
                    options.push_back("Point Light " + std::to_string(i + 1));
                }
                for (size_t i = 0; i < scene.sphericalLight.size(); i++) {
                    options.push_back("Spherical Light " + std::to_string(i + 1));
                }
				for (size_t i = 0; i < scene.spotLight.size(); i++) {
					options.push_back("Spot Light " + std::to_string(i + 1));
				}
				for (size_t i = 0; i < scene.planeLight.size(); i++) {
					options.push_back("Plane Light " + std::to_string(i + 1));
				}

                std::vector<const char*> optionsPointers;
                std::transform(std::begin(options), std::end(options), std::back_inserter(optionsPointers),
                    [](const auto& str) { return str.c_str(); });

                ImGui::Combo("Selected light", &selectedLight, optionsPointers.data(), static_cast<int>(optionsPointers.size()));
            }

            {
                const auto showLightOptions = [](auto& light) {
                    ImGui::DragFloat3("Light position", glm::value_ptr(light.position), 0.01f, -3.0f, 3.0f);
                    ImGui::ColorEdit3("Light color", glm::value_ptr(light.color));
                    if constexpr (std::is_same_v<std::decay_t<decltype(light)>, SphericalLight>) {
                        ImGui::DragFloat("Light radius", &light.radius, 0.01f, 0.01f, 0.5f);
                    }
					if constexpr (std::is_same_v<std::decay_t<decltype(light)>, SpotLight>) {
						ImGui::DragFloat("Max light angle", &light.angle, 0.5f, 1.0f, 20.0f);
						ImGui::DragFloat3("Light direction", glm::value_ptr(light.direction), 0.01f, -1.0f, 1.0f);
					}
					if constexpr (std::is_same_v<std::decay_t<decltype(light)>, PlaneLight>) {
						ImGui::DragFloat3("Width", glm::value_ptr(light.width), 0.01f, -3.0f, 3.0f);
						ImGui::DragFloat3("Height", glm::value_ptr(light.height), 0.01f, -3.0f, 3.0f);
					}
                };
                if (selectedLight < static_cast<int>(scene.pointLights.size())) {
                    // Draw a big yellow sphere and then the small light sphere on top.
                    showLightOptions(scene.pointLights[selectedLight]);
                }
				else if (selectedLight < static_cast<int>(scene.pointLights.size()+scene.sphericalLight.size())) {
                    // Draw a big yellow sphere and then the smaller light sphere on top.
                    showLightOptions(scene.sphericalLight[selectedLight - scene.pointLights.size()]);
                }
				else if (selectedLight < static_cast<int>(scene.pointLights.size() + scene.sphericalLight.size() + scene.spotLight.size() )) {
					showLightOptions(scene.spotLight[selectedLight - (scene.pointLights.size() + scene.sphericalLight.size())]);
				}
				else {
					showLightOptions(scene.planeLight[selectedLight - (scene.pointLights.size() + scene.sphericalLight.size() + scene.spotLight.size())]);
				}
            }
        }
		if (ImGui::TreeNode("Add/remove light")) {
			if (ImGui::Button("Add point light")) {
				scene.pointLights.push_back(PointLight{ glm::vec3(0.0f), glm::vec3(1.0f) });
				selectedLight = int(scene.pointLights.size() - 1);
			}
			if (ImGui::Button("Add spherical light")) {
				scene.sphericalLight.push_back(SphericalLight{ glm::vec3(0.0f), 0.1f, glm::vec3(1.0f) });
				selectedLight = int(scene.pointLights.size() + scene.sphericalLight.size() - 1);
			}
			if (ImGui::Button("Add spot light")) {
				scene.spotLight.push_back(SpotLight{ glm::vec3(0.0f), glm::vec3(1), 10, glm::vec3(1.0f) });
				selectedLight = int(scene.pointLights.size() + scene.sphericalLight.size() + scene.spotLight.size() - 1);
			}
			if (ImGui::Button("Add plane light")) {
				scene.planeLight.push_back(PlaneLight{ glm::vec3(0.0f), glm::vec3(0.2f,0,0), glm::vec3(0,0,0.2f), glm::vec3(1.0f) });
				selectedLight = int(scene.pointLights.size() + scene.sphericalLight.size() + scene.spotLight.size() + scene.planeLight.size() - 1);
			}
			if (ImGui::Button("Remove selected light")) {
				if (selectedLight < static_cast<int>(scene.pointLights.size())) {
					scene.pointLights.erase(std::begin(scene.pointLights) + selectedLight);
				}
				else if (selectedLight < static_cast<int>(scene.pointLights.size() + scene.sphericalLight.size())) {
					scene.sphericalLight.erase(std::begin(scene.sphericalLight) + (selectedLight - scene.pointLights.size()));
				}
				else if (selectedLight < static_cast<int>(scene.pointLights.size() + scene.sphericalLight.size() + scene.spotLight.size())) {
					scene.spotLight.erase(std::begin(scene.spotLight) + (selectedLight - (scene.pointLights.size() + scene.sphericalLight.size())));
				}
				else {
					scene.planeLight.erase(std::begin(scene.planeLight) + (selectedLight - (scene.pointLights.size() + scene.sphericalLight.size() + scene.spotLight.size())));
				}
				selectedLight = 0;
			}
			ImGui::TreePop();
		}

        // Clear screen.
        glClearDepth(1.0f);
        glClearColor(0.0, 0.0, 0.0, 0.0);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	

		// Draw either using OpenGL (rasterization) or the ray tracing function.
        switch (viewMode) {
        case ViewMode::Rasterization: {
            glPushAttrib(GL_ALL_ATTRIB_BITS);
            renderOpenGL(scene, camera, selectedLight);
            if (optDebugRay) {
                // Call getFinalColor for the debug ray. Ignore the result but tell the function that it should
                // draw the rays instead.
                enableDrawRay = true;
                (void)getFinalColor(scene, bvh, *optDebugRay);
                enableDrawRay = false;
            }

			glPopAttrib();
        } break;
        case ViewMode::RayTracing: {
			if(rayTracing_toBeRendered) {
				screen.clear(glm::vec3(0.0f));
				renderRayTracing(scene, camera, bvh, screen, false, anti_aliasing, multipleRays, sampleSize);
				rayTracing_toBeRendered = false;
			}
			screen.setPixel(0, 0, glm::vec3(1.0f));
			screen.draw(); // Takes the image generated using ray tracing and outputs it to the screen using OpenGL.
		} break;
		case ViewMode::Textures: {
			// render in textureDebugging mode
			if(rayTracing_toBeRendered) {
				screen.clear(glm::vec3(0.0f));
				renderRayTracing(scene, camera, bvh, screen, true);
				rayTracing_toBeRendered = false;
			}
			screen.setPixel(0, 0, glm::vec3(1.0f));
			screen.draw(); // Takes the image generated using ray tracing and outputs it to the screen using OpenGL.
		} break;
		case ViewMode::MultipleRayDebug: {
			glPushAttrib(GL_ALL_ATTRIB_BITS);
			renderOpenGL(scene, camera, selectedLight);
			if (optDebugRay) {
				// Call getFinalColor for the debug ray. Ignore the result but tell the function that it should
				// draw the rays instead.
				for (auto& ray : debugRays) {
					enableDrawRay = true;
					(void)getFinalColor(scene, bvh, ray);
					enableDrawRay = false;
				}
			}

			glPopAttrib();
		} break;
        default:
            break;
        };

        if (debugBVH) {
            glPushAttrib(GL_ALL_ATTRIB_BITS);
            setOpenGLMatrices(camera);
            glDisable(GL_LIGHTING);
            glEnable(GL_DEPTH_TEST);

            // Enable alpha blending. More info at:
            // https://learnopengl.com/Advanced-OpenGL/Blending
            glEnable(GL_BLEND);
            glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
            bvh.debugDraw(bvhDebugLevel, bvhShowLeafNodes);

			glPopAttrib();
        }


		ImGui::End();
        window.swapBuffers();
    }

    return 0; // execution never reaches this point
}

static void setOpenGLMatrices(const Trackball& camera)
{
	// Load view matrix.
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	const glm::mat4 viewMatrix = camera.viewMatrix();
	glMultMatrixf(glm::value_ptr(viewMatrix));

	// Load projection matrix.
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	const glm::mat4 projectionMatrix = camera.projectionMatrix();
	glMultMatrixf(glm::value_ptr(projectionMatrix));
}

static void renderOpenGL(const Scene& scene, const Trackball& camera, int selectedLight)
{
	// Normals will be normalized in the graphics pipeline.
	glEnable(GL_NORMALIZE);
	// Activate rendering modes.
	glEnable(GL_DEPTH_TEST);
	// Draw front and back facing triangles filled.
	glPolygonMode(GL_FRONT, GL_FILL);
	glPolygonMode(GL_BACK, GL_FILL);
	// Interpolate vertex colors over the triangles.
	glShadeModel(GL_SMOOTH);
	setOpenGLMatrices(camera);

	glDisable(GL_LIGHTING);
	// Render point lights as very small dots
	for (const auto& light : scene.pointLights)
		drawSphere(light.position, 0.01f, light.color);
	for (const auto& light : scene.sphericalLight)
		drawSphere(light.position, light.radius, light.color);

	if (!scene.pointLights.empty() || !scene.sphericalLight.empty() || !scene.spotLight.empty() || !scene.planeLight.empty()) {
		if (selectedLight < static_cast<int>(scene.pointLights.size())) {
			// Draw a big yellow sphere and then the small light sphere on top.
			const auto& light = scene.pointLights[selectedLight];
			drawSphere(light.position, 0.05f, glm::vec3(1, 1, 0));
			glDisable(GL_DEPTH_TEST);
			drawSphere(light.position, 0.01f, light.color);
			glEnable(GL_DEPTH_TEST);
		}
		else if(selectedLight < static_cast<int>(scene.pointLights.size() + scene.sphericalLight.size())) {
			// Draw a big yellow sphere and then the smaller light sphere on top.
			const auto& light = scene.sphericalLight[selectedLight - scene.pointLights.size()];
			drawSphere(light.position, light.radius + 0.01f, glm::vec3(1, 1, 0));
			glDisable(GL_DEPTH_TEST);
			drawSphere(light.position, light.radius, light.color);
			glEnable(GL_DEPTH_TEST);
		}
		else if (selectedLight < static_cast<int>(scene.pointLights.size() + scene.sphericalLight.size() + scene.spotLight.size())) {
			// Draw a big yellow sphere and then the small light sphere on top.
			const auto& light = scene.spotLight[selectedLight - (scene.pointLights.size() + scene.sphericalLight.size())];
			drawSphere(light.position, 0.05f, glm::vec3(1, 1, 0));
			glDisable(GL_DEPTH_TEST);
			drawSphere(light.position, 0.01f, light.color);
			glEnable(GL_DEPTH_TEST);
		}
		else {
			const auto& light = scene.planeLight[selectedLight - (scene.pointLights.size() + scene.sphericalLight.size() + scene.spotLight.size())];
			drawPlane(light.position, light.width, light.height, glm::vec3(1, 1, 0));
			glDisable(GL_DEPTH_TEST);
			drawSphere(light.position, 0.01f, light.color);
			glEnable(GL_DEPTH_TEST);
		}
	}

	// Activate the light in the legacy OpenGL mode.
	glEnable(GL_LIGHTING);

	int i = 0;
	const auto enableLight = [&](const auto& light) {
		glEnable(GL_LIGHT0 + i);
		const glm::vec4 position4{ light.position, 1 };
		glLightfv(GL_LIGHT0 + i, GL_POSITION, glm::value_ptr(position4));
		const glm::vec4 color4{ glm::clamp(light.color, 0.0f, 1.0f), 1.0f };
		const glm::vec4 zero4{ 0.0f, 0.0f, 0.0f, 1.0f };
		glLightfv(GL_LIGHT0 + i, GL_AMBIENT, glm::value_ptr(zero4));
		glLightfv(GL_LIGHT0 + i, GL_DIFFUSE, glm::value_ptr(color4));
		glLightfv(GL_LIGHT0 + i, GL_SPECULAR, glm::value_ptr(zero4));
		// NOTE: quadratic attenuation doesn't work like you think it would in legacy OpenGL.
		// The distance is not in world space but in NDC space!
		glLightf(GL_LIGHT0 + i, GL_CONSTANT_ATTENUATION, 1.0f);
		glLightf(GL_LIGHT0 + i, GL_LINEAR_ATTENUATION, 0.0f);
		glLightf(GL_LIGHT0 + i, GL_QUADRATIC_ATTENUATION, 0.0f);
		i++;
	};
	for (const auto& light : scene.pointLights)
		enableLight(light);
	for (const auto& light : scene.sphericalLight)
		enableLight(light);
	for (const auto& light : scene.spotLight)
		enableLight(light);
	for (const auto& light : scene.planeLight)
		enableLight(light);

	// Draw the scene and the ray (if any).
	drawScene(scene);

	// Draw a colored sphere at the location at which the trackball is looking/rotating around.
	glDisable(GL_LIGHTING);
	drawSphere(camera.lookAt(), 0.01f, glm::vec3(0.2f, 0.2f, 1.0f));
}
