

#pragma once
#include <vector>
#include "scene.h"
#include "bounding_volume_hierarchy.h"

#define SHADOW_ERROR_OFFSET 0.0005f

//#define USE_MONTE_CARLO_INTEGRATION

bool isInHardShadow(const glm::vec3 point, const Scene& scene, const BoundingVolumeHierarchy& bvh);

struct Lighting {
	glm::vec3 color;
	glm::vec3 position;
	float intensity;
	float cosLightSurfaceAngle;
	float cosLightSpecAngle;
};


std::vector<Lighting> getPointLights(const HitInfo& point, const glm::vec3 reflectdir, Scene& scene, const BoundingVolumeHierarchy& bvh);

struct SoftShadowSphere {
	SphericalLight light;
	float intensity;
};
std::vector<Lighting> getSpherelights(const HitInfo& point, const glm::vec3& reflect, Scene& scene, const BoundingVolumeHierarchy& bvh, int rayCount = 7);

std::vector<Lighting> getSpotLichts(const HitInfo& point, const glm::vec3& reflectdir, Scene& scene, const BoundingVolumeHierarchy& bvh);

struct SoftShadowPlane
{
	PlaneLight light;
	float intensity;
	float cosAngle;
};
std::vector<Lighting> getPlaneLights(const HitInfo& point, const glm::vec3& reflectdir, Scene& scene, const BoundingVolumeHierarchy& bvh, int rayCount1D = 3);

