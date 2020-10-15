#include "shadow.h"

#include <vector>
#include "scene.h"
#include "bounding_volume_hierarchy.h"

#include <glm/glm.hpp>
#include "draw.h"

// Check if the given point is in a hard shadow.
// When there are multyple light scources the function only returns true if the point can not see any point light.
bool isInHardShadow(const glm::vec3 point, const Scene& scene, const BoundingVolumeHierarchy& bvh) {

	// check if the point can see any light
	for (const PointLight& light : scene.pointLights) {

		// calculate the ray
		glm::vec3 o = point;
		glm::vec3 d = light.position - point;
		float distance = glm::length(d);
		d = glm::normalize(d);
		o += SHADOW_ERROR_OFFSET * d;
		Ray ray = { o, d };

		// see if it hit somthing
		HitInfo hitInfo;
		bool hit = bvh.intersect(ray, hitInfo);

		// the ray hit nothing or the hit was only after the light
		if (!hit || (ray.t > distance - 2 * SHADOW_ERROR_OFFSET)) {
			drawRay(ray, glm::vec3(0.5f, 1.0f, 0.5f));
			// the point is illuminated by a light so we return false
			return false;
		}
		else {
			drawRay(ray, glm::vec3(0.8f, 0.2f, 0.0f));
			return true;
		}
	}
	// no light has been found that illuminates the point
}

// Returns a list of PointLights that illuminate the given point.
std::vector<PointLight> getVisablePointLights(const glm::vec3 point, const Scene& scene, const BoundingVolumeHierarchy& bvh) {

	std::vector<PointLight> visibleLights;

	// check if the point can see any light
	for (const PointLight& light : scene.pointLights) {

		// calculate the ray
		glm::vec3 o = point;
		glm::vec3 d = light.position - point;
		float distance = glm::length(d);
		d = glm::normalize(d);
		o += SHADOW_ERROR_OFFSET * d;
		Ray ray = { o, d };

		// see if it hit somthing
		HitInfo hitInfo;
		bool hit = bvh.intersect(ray, hitInfo);

		// the ray hit nothing or the hit was only after the light
		if (!hit || (ray.t > distance - 2 * SHADOW_ERROR_OFFSET)) {
			// debug ray when the point is illuminated
			drawRay(ray, glm::vec3(0.5f, 1.0f, 0.5f));

			// store the light scources that we can see
			visibleLights.push_back(light);
		}
		else {
			// debug ray when the point is in shadow
			drawRay(ray, glm::vec3(0.8f, 0.2f, 0.0f));
		}
	}
	return visibleLights;
}

bool isHitFromLight(const glm::vec3 point, const PointLight& light, const BoundingVolumeHierarchy& bvh) {

	// calculate the ray
	glm::vec3 o = point;
	glm::vec3 d = light.position - point;
	float distance = glm::length(d);
	d = glm::normalize(d);
	o += SHADOW_ERROR_OFFSET * d;
	Ray ray = { o, d };

	// see if it hit somthing
	HitInfo hitInfo;
	bool hit = bvh.intersect(ray, hitInfo);

	// the ray hit nothing or the hit was only after the light
	if (!hit || (ray.t > distance - 2 * SHADOW_ERROR_OFFSET)) {
		drawRay(ray, glm::vec3(0.5f, 1.0f, 0.5f));
		// the point is illuminated by a light so we return true
		return true;
	}
	else {
		drawRay(ray, glm::vec3(0.8f, 0.2f, 0.0f));
		return false;
	}
}