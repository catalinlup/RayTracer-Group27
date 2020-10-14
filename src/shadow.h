#include <vector>
#include "scene.h"
#include "bounding_volume_hierarchy.h"

#define SHADOW_ERROR_OFFSET 0.0005f

bool isInHardShadow(const glm::vec3 point, const Scene& scene, const BoundingVolumeHierarchy& bvh);

std::vector<PointLight> getVisablePointLights(const glm::vec3 point, const Scene& scene, const BoundingVolumeHierarchy& bvh);