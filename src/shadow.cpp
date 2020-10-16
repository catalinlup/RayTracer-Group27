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
		}
	}
	// no light has been found that illuminates the point
	return true;
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
		drawRay(ray, glm::vec3(0.5f));

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

// a 3d rotation matrix
glm::mat3 rotatetionMatrix(const float angle, const glm::vec3 axis) {
	glm::mat3 C = { {0, axis.z, -axis.y}, {-axis.z, 0, axis.x}, {axis.y, -axis.x, 0} };
	return glm::mat3(1) + C * std::sin(angle) + C * C * (1 - std::cos(angle));
}
// returns a list of visable light scoures and how mutch of that lightscource is visable
std::vector<SoftShadow> getSoftlights(const glm::vec3 point, const Scene& scene, const BoundingVolumeHierarchy& bvh, int rayCount) {
	std::vector<SoftShadow> lights;
	for (const SphericalLight& light : scene.sphericalLight) {

		// calculate the starting position of the rays to the light
		glm::vec3 o = point;
		glm::vec3 d = light.position - point;
		float distance = glm::length(d);
		d = glm::normalize(d);
		o += SHADOW_ERROR_OFFSET * d;

		// ray for the center of the light
		Ray ray = { o, d };

		// the number of rays that hit the light
		int hits = 0;

		// see if the ray to the center hit somthing
		HitInfo hitInfo;
		bool hit = bvh.intersect(ray, hitInfo);

		// the ray hit nothing or the hit was only after the light
		if (!hit || (ray.t > distance - 2 * SHADOW_ERROR_OFFSET)) {
			hits++;
		}


		// thare are rayCount-1 rayes that go to the edge of the circle
		glm::mat3 rotate = rotatetionMatrix(2 * 3.14159365358979f / (rayCount - 1), d); // rotaion matrix around the vector d

		// a matrix that is not in line whith d
		glm::vec3 notd = d;
		if (d.x != 0) {
			notd.y = -d.x;
			notd.x = d.y;
		}
		else {
			notd.y = -d.z;
			notd.z = d.y;
		}

		// perpendicular to the ray t the center of the light whith the same length as the radius of the light
		glm::vec3 perp = glm::normalize(glm::cross(d, notd)) * light.radius; 

		// loop over all the rays that go the the edge of the circle
		for (int i = 0; i < rayCount - 1; i++) {
			Ray edgeRay = { o,d + perp };
			bool hit = bvh.intersect(edgeRay, hitInfo);
			if (!hit || (edgeRay.t > distance - 2 * SHADOW_ERROR_OFFSET)) {
				hits++;
			}
			drawRay(edgeRay, { 1.0f, 1.0f, 0.0f });
			// rotate the perpenicular vector to be able to calculate the next position on the light where a ray must go
			perp = rotate * perp;
		}

		// save the visable lights
		if (hits > 0) {
			lights.push_back({ light, hits / (float)rayCount });

		}
	}
	return lights;
}