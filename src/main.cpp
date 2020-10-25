#include "bounding_volume_hierarchy.h"
#include "disable_all_warnings.h"
#include "draw.h"
#include "image.h"
#include "ray_tracing.h"
#include "screen.h"
#include "trackball.h"
#include "window.h"
#include "shadow.h"

// Disable compiler warnings in third-party code (which we cannot change).
DISABLE_WARNINGS_PUSH()
#include <glm/gtc/type_ptr.hpp>
#include <glm/vec2.hpp>
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

// Textures settings - used for debugging only!
TextureFiltering textureFiltering{TextureFiltering::NearestNeighbor};
OutOfBoundsRule outOfBoundsRuleX{OutOfBoundsRule::Border}; // Border, Clamp, Repeat
OutOfBoundsRule outOfBoundsRuleY{OutOfBoundsRule::Border};
glm::vec3 textureBorderColor(0);


enum class ViewMode {
	Rasterization = 0,
	RayTracing = 1,
	Textures = 2
};


bool checkShadow(HitInfo hitInfo, PointLight light, const BoundingVolumeHierarchy& bvh) {
	Ray testRay;
	HitInfo old = hitInfo;
	testRay.origin = old.hitPoint;
	glm::vec3 direction = light.position - old.hitPoint;
	testRay.direction = glm::normalize(direction);
	if (bvh.intersect(testRay, hitInfo)) {
		glm::vec3 t = (hitInfo.hitPoint - testRay.origin) / direction;
		if (glm::all(glm::greaterThan(t, glm::vec3(0))) && glm::all(glm::lessThan(t, glm::vec3(1)))) {
			drawRay(testRay, glm::vec3(1, 0, 0));
			return false;
		}
	}
	drawRay(testRay, glm::vec3(1, 1, 0));
	return true;
}

glm::vec3 calcSpecular(int level, const BoundingVolumeHierarchy& bvh, PointLight light, Ray ray, HitInfo hitInfo, glm::vec3 cameraPos) {
	glm::vec3 lightDir = glm::normalize(light.position - hitInfo.hitPoint);
	glm::vec3 resColor(0);
	//Check if the intersecting surface has a non black Ks value and that we haven't passed the relfected ray count
	if (glm::all(glm::notEqual(hitInfo.material.ks, glm::vec3(0))) && level > 0) {
		if (checkShadow(hitInfo, light, bvh)) {
			glm::vec3 reflectedLight = 2 * glm::dot(lightDir, hitInfo.normal) * hitInfo.normal - lightDir;
			glm::vec3 viewDir = glm::normalize(cameraPos - hitInfo.hitPoint);

			resColor = hitInfo.material.ks * light.color * glm::pow(glm::dot(glm::normalize(reflectedLight), viewDir), hitInfo.material.shininess);

			//Keep the old hitInfo in case the reflected ray has no further intersections
			HitInfo old = hitInfo;

			Ray reflectedRay;
			reflectedRay.origin = old.hitPoint;
			reflectedRay.direction = 2 * glm::dot(old.normal, viewDir) * old.normal - viewDir;

			Ray normalRay;
			normalRay.origin = old.hitPoint;
			normalRay.direction = old.normal;

			Ray lightRay;
			lightRay.origin = old.hitPoint;
			lightRay.direction = lightDir;

			if (bvh.intersect(reflectedRay, hitInfo)) {
				//std::cout << "Refl Hit" << std::endl;
				drawRay(reflectedRay, glm::vec3(1));
				//drawRay(normalRay, glm::vec3(1, 0, 0));
				level--;
				resColor += calcSpecular(level, bvh, light, reflectedRay, hitInfo, cameraPos);
			}
			else {
				//Restore old hitInfo because there was no hit
				hitInfo = old;
				//std::cout << "No Hit" << std::endl;
				drawRay(reflectedRay, glm::vec3(1));
				//Stop the recursion
				level = 0;
				return resColor;
			}
		}
		//Point in shadow skip computations
		else {
			return resColor;
		}
	}
	//Black specularity return 0 vector back
	else {
		
		return glm::vec3(0);
	}

}





// used for debugging the textures
static glm::vec3 getFinalColorNoRayTracingJustTextures(const Scene &scene, const BoundingVolumeHierarchy &bvh, Ray ray) {
	HitInfo hitInfo;



	if(bvh.intersect(ray, hitInfo)) {
		Material mat = hitInfo.material;
		if(mat.kdTexture) {
			Image texture = mat.kdTexture.value();

			texture.setBorderColor(textureBorderColor);
			texture.setOutOfBoundsRuleX(outOfBoundsRuleX);
			texture.setOutOfBoundsRuleY(outOfBoundsRuleY);
			texture.setTextureFilteringMethod(textureFiltering);

			// verticies
			glm::vec2 textureCoordinate = hitInfo.texCoord;

			//return glm::vec3(1);

			// return the color corresponding to the texture
			return texture.getPixel(textureCoordinate);
		}

		return glm::vec3(1);
	}

	return glm::vec3(0);
}


glm::vec3 calcColor(Lighting light, Material material) {
	// difuse light
	glm::vec3 diffuse = material.kd * light.color * light.intensity * light.cosLightSurfaceAngle;
	//specular light
	glm::vec3 spec = light.color * material.ks * std::pow(light.cosLightSpecAngle, material.shininess);
	return diffuse + spec;
}

static int max_reflection_level = 5;
static int sphere_light_ray_count = 10;
static int plane_light_1D_ray_count = 3;
static int glossy_ray_count = 10;
static float refraction_factor = 0.8;
// NOTE(Mathijs): separate function to make recursion easier (could also be done with lambda + std::function).
static glm::vec3 getFinalColor(const Scene& scene, const BoundingVolumeHierarchy& bvh, Ray ray, int level=0)
{
	HitInfo hitInfo;

	/*For every light calulate the addition from the reflected rays,
	then add all of them*/
	if (bvh.intersect(ray, hitInfo)) {
		drawRay(ray, glm::vec3(1));
		glm::vec3 color(0);
		hitInfo.normal = glm::normalize(hitInfo.normal);
		glm::vec3 reflect = glm::reflect(glm::normalize(ray.direction), hitInfo.normal);
		
		// Loop over all the lights that are not in shadow
		for (const Lighting& light : getPointLights(hitInfo, reflect, scene, bvh)) {
			color += calcColor(light, hitInfo.material);
		}
		for (const Lighting& light : getSpherelights(hitInfo, reflect, scene, bvh, sphere_light_ray_count)) {
			color += calcColor(light, hitInfo.material);
		}
		for (const Lighting& light : getSpotLichts(hitInfo, reflect, scene, bvh)) {
			color += calcColor(light, hitInfo.material);
		}
		for (const Lighting& light : getPlaneLights(hitInfo, reflect, scene, bvh, plane_light_1D_ray_count)) {
			color += calcColor(light, hitInfo.material);
		}

		if (level >= max_reflection_level) {
			return color;
		}
		
		if (hitInfo.material.transparency == 1.0f) { // not tranparent
			// calculate reflected light
			if (hitInfo.material.ks.x > 0 || hitInfo.material.ks.y > 0 || hitInfo.material.ks.z > 0) {

				glm::vec3 reflectColor = glm::vec3(0);
				// the ( + 0.01f * reflect ) is to prevent a surface reflecting itself
				Ray refRay = { hitInfo.hitPoint + 0.01f * reflect, reflect };
				// Calculate the color that is reflected
				reflectColor += hitInfo.material.ks * getFinalColor(scene, bvh, refRay, level + 1);
				
				// prevent divide by 0 and allow for prefect reflections
				if (hitInfo.material.shininess != 0) {

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
					float d = std::pow(0.5f, -1 / (float)hitInfo.material.shininess) * std::sqrt(1-std::pow(0.5, 2/(float)hitInfo.material.shininess));

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
							glm::vec3 ctmp = getFinalColor(scene, bvh, shineRay, level + 1) *std::max(std::pow(glm::dot(reflect, shineDir), hitInfo.material.shininess), 0.0f);
							//std::cout << "r_r: " << ctmp.r << std::endl;
							reflectColor += ctmp;
						}
					}
					//std::cout << "r_final: " << (glossColor / (float)glossrays).r << std::endl;
					color += hitInfo.material.ks * reflectColor / (float)glossy_ray_count;

				}
				else {
					color += hitInfo.material.ks * reflectColor;
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
			glm::vec3& n = hitInfo.normal;

			float r = refraction_factor; // n1/n2 : indexes of refraction
			float c = std::abs(glm::dot(l, n));

			    // vector form of snell's law
			glm::vec3 refract = r*l + ( r*c - std::sqrt(1 - r*r*(1 - c*c )) )*n; // std::sqrt might return NaN
			refract = glm::normalize(refract);
			

			// calculate how how much light should reflect and how much should refract
			float& R0 = hitInfo.material.transparency;//R0static; // probebility of reflection at 0 rad/deg to normal
			    
				// simplified fersnel equasion
			float reflectionChance = R0 + (1-R0)*(std::pow(1-c,5));
			float refractionChance = 1 - reflectionChance;
			

			// get the color the reflected and refrected rays see
			color += reflectionChance * getFinalColor(scene, bvh, { hitInfo.hitPoint + 0.01f * reflect, reflect }, level + 1);
			if (r * r * (1 - c * c) <= 1.0f) { // check if total internal reflection occures
				color += refractionChance * getFinalColor(scene, bvh, { hitInfo.hitPoint + 0.01f * refract, refract }, level + 1);
			}

		}
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



// This is the main rendering function. You are free to change this function in any way (including the function signature).
// if texture debugging is set to true, the method will call the getFinalColorNoRayTracingJustTextures, instead of getFinalColor.
static void renderRayTracing(const Scene& scene, const Trackball& camera, const BoundingVolumeHierarchy& bvh, Screen& screen, bool textureDebugging = false)
{
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
			else
				screen.setPixel(x, y, getFinalColor(scene, bvh, cameraRay));
		}
	}
}

int main(int argc, char** argv)
{
    Trackball::printHelp();
    std::cout << "\n Press the [R] key on your keyboard to create a ray towards the mouse cursor" << std::endl
              << std::endl;

    Window window { "Final Project - Part 2", windowResolution, OpenGLVersion::GL2 };
    Screen screen { windowResolution };
    Trackball camera { &window, glm::radians(50.0f), 3.0f };
    camera.setCamera(glm::vec3(0.0f, 0.0f, 0.0f), glm::radians(glm::vec3(20.0f, 20.0f, 0.0f)), 3.0f);

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

	ViewMode viewMode { ViewMode::Rasterization };

	bool drawWhenInTextureMode {false}; // controls whether the object is drawn or not when in texture mode. Used to reduce lag

	

    window.registerKeyCallback([&](int key, int /* scancode */, int action, int /* mods */) {
        if (action == GLFW_PRESS) {
            switch (key) {
            case GLFW_KEY_R: {
                // Shoot a ray. Produce a ray from camera to the far plane.
                const auto tmp = window.getNormalizedCursorPos();
                optDebugRay = camera.generateRay(tmp * 2.0f - 1.0f);
                viewMode = ViewMode::Rasterization;
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
				break;
			}
			case GLFW_KEY_C: {
				viewMode = ViewMode::Textures;
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
            constexpr std::array items { "SingleTriangle", "Cube", "Cornell Box (with mirror)", "Cornell Box (spherical light and mirror)", "Cornell Box (plane light and mirror)", "Monkey", "Teapot", "Dragon", /* "AABBs",*/ "Spheres", /*"Mixed",*/ "Custom" };
            if (ImGui::Combo("Scenes", reinterpret_cast<int*>(&sceneType), items.data(), int(items.size()))) {
                optDebugRay.reset();
                scene = loadScene(sceneType, dataPath);
                bvh = BoundingVolumeHierarchy(&scene);
                if (optDebugRay) {
                    HitInfo dummy {};
                    bvh.intersect(*optDebugRay, dummy);
                }
			}
		}
        {
            constexpr std::array items { "Rasterization", "Ray Traced", "Textures" };
            ImGui::Combo("View mode", reinterpret_cast<int*>(&viewMode), items.data(), int(items.size()));
        }
        if (ImGui::Button("Render to file")) {
            {
                using clock = std::chrono::high_resolution_clock;
                const auto start = clock::now();
                renderRayTracing(scene, camera, bvh, screen);
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

        ImGui::Spacing();
        ImGui::Separator();
        ImGui::Text("Debugging");
        if (viewMode == ViewMode::Rasterization) {
            ImGui::Checkbox("Draw BVH", &debugBVH);
            if (debugBVH) {
                ImGui::SliderInt("BVH Level", &bvhDebugLevel, 0, bvh.numLevels() - 1);
                ImGui::Checkbox("Show BVH leaf nodes", &bvhShowLeafNodes);
            }
            
        } else if(viewMode == ViewMode::Textures) {
			{
				ImGui::Checkbox("Draw", &drawWhenInTextureMode);

				constexpr std::array filtering_modes{"Nearest Neighbor", "Bilinear", "Nearest Level Mipmapping", "Bilinear Mipmapping", "Trilinear Mipmapping"};
				ImGui::Combo("Filteting Mode", reinterpret_cast<int *>(&textureFiltering), filtering_modes.data(), int(filtering_modes.size()));

				constexpr std::array oob_x{"Border", "Clamping", "Repeat"};
				ImGui::Combo("Out Of Bounds X", reinterpret_cast<int *>(&outOfBoundsRuleX), oob_x.data(), int(oob_x.size()));

				constexpr std::array oob_y{"Border", "Clamping", "Repeat"};
				ImGui::Combo("Out Of Bounds Y", reinterpret_cast<int *>(&outOfBoundsRuleY), oob_y.data(), int(oob_y.size()));

				ImGui::ColorEdit3("Border color", glm::value_ptr(textureBorderColor));
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

			

            screen.clear(glm::vec3(0.0f));
            renderRayTracing(scene, camera, bvh, screen);
            screen.setPixel(0, 0, glm::vec3(1.0f));
            screen.draw(); // Takes the image generated using ray tracing and outputs it to the screen using OpenGL.
        } break;
		case ViewMode::Textures: {
			screen.clear(glm::vec3(0.0f));
			// render in textureDebugging mode
			if(drawWhenInTextureMode) renderRayTracing(scene, camera, bvh, screen, true);
			screen.setPixel(0, 0, glm::vec3(1.0f));
			screen.draw(); // Takes the image generated using ray tracing and outputs it to the screen using OpenGL.
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
