#include "shadow.h"

#include <vector>
#include "scene.h"
#include "bounding_volume_hierarchy.h"

#include <glm/glm.hpp>
#include "draw.h"

// check if one point can see an other
bool cansee(const glm::vec3& p1, const glm::vec3& p2, const Scene& scene, const BoundingVolumeHierarchy& bvh) {
	glm::vec3 o = p1;
	glm::vec3 d = p2 - p1;
	float distance = glm::length(d);
	d = glm::normalize(d);
	o += SHADOW_ERROR_OFFSET * d;
	Ray ray = { o, d };

	HitInfo hitInfo;
	bool hit = bvh.intersect(ray, hitInfo, false);

	if (!hit || (ray.t > distance - 2 * SHADOW_ERROR_OFFSET)) {
		drawRay(ray, glm::vec3(0.5f, 1.0f, 0.5f));
		return true;
	}
	else {
		drawRay(ray, glm::vec3(0.8f, 0.2f, 0.0f));
		return false;
	}
}



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
		bool hit = bvh.intersect(ray, hitInfo, false);

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
std::vector<Lighting> getPointLights(const HitInfo& point, const glm::vec3 reflectdir, const Scene& scene, const BoundingVolumeHierarchy& bvh) {

	std::vector<Lighting> visibleLights;

	// check if the point can see any light
	for (const PointLight& light : scene.pointLights) {

		// the ray hit nothing or the hit was only after the light
		if (cansee(point.hitPoint, light.position, scene, bvh)) {
			// store the light scources that we can see
			Lighting l;
			l.color = light.color;
			l.position = light.position;
			l.intensity = 1;
			l.cosLightSurfaceAngle = std::max(0.0f, glm::dot(point.normal, glm::normalize(light.position - point.hitPoint)));
			l.cosLightSpecAngle = std::max(0.0f, glm::dot(glm::normalize(reflectdir), glm::normalize(light.position - point.hitPoint)));
			visibleLights.push_back(l);
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
std::vector<Lighting> getSpherelights(const HitInfo& point, const glm::vec3& reflectdir, const Scene& scene, const BoundingVolumeHierarchy& bvh, int rayCount) {
	std::vector<Lighting> lights;
	for (const SphericalLight& light : scene.sphericalLight) {

		// calculate the starting position of the rays to the light
		glm::vec3 o = point.hitPoint;
		glm::vec3 d = light.position - point.hitPoint;
		float distance = glm::length(d);
		d = glm::normalize(d);
		o += SHADOW_ERROR_OFFSET * d;

		// ray for the center of the light
		Ray ray = { o, d };

		// the number of rays that hit the light
		int hits = 0;

		// see if the ray to the center hit somthing
		HitInfo hitInfo;
		bool hit = bvh.intersect(ray, hitInfo, false);

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
			bool hit = bvh.intersect(edgeRay, hitInfo, false);
			if (!hit || (edgeRay.t > distance - 2 * SHADOW_ERROR_OFFSET)) {
				hits++;
			}
			drawRay(edgeRay, { 1.0f, 1.0f, 0.0f });
			// rotate the perpenicular vector to be able to calculate the next position on the light where a ray must go
			perp = rotate * perp;
		}

		// save the visable lights
		if (hits > 0) {
			// setting all the properites of this light
			Lighting l;
			l.color = light.color;
			l.intensity = hits / (float)rayCount;
			l.position = light.position;
			l.cosLightSurfaceAngle = std::max(0.0f, glm::dot(point.normal, glm::normalize(light.position - point.hitPoint)));
			l.cosLightSpecAngle = std::max(0.0f, glm::dot(glm::normalize(reflectdir), glm::normalize(light.position - point.hitPoint)));
			lights.push_back(l);

		}
	}
	return lights;
}

// returns a list of all spotlights that iluminate a point
std::vector<Lighting> getSpotLichts(const HitInfo& point, const glm::vec3& reflectdir, const Scene& scene, const BoundingVolumeHierarchy& bvh) {
	std::vector<Lighting> lights;

	for (const SpotLight& light : scene.spotLight) {

		// if the point is in the beam of the spot light
		if (glm::dot(glm::normalize(light.direction), glm::normalize(point.hitPoint - light.position)) > std::cos(glm::radians(light.angle))) {

			// check if the path is not obstructed
			if (cansee(point.hitPoint, light.position, scene, bvh)) {
				// setting all the properites of this light
				Lighting l;
				l.color = light.color;
				l.position = light.position;
				l.intensity = 1;
				l.cosLightSurfaceAngle = std::max(0.0f, glm::dot(point.normal, glm::normalize(light.position - point.hitPoint)));
				l.cosLightSpecAngle = std::max(0.0f, glm::dot(glm::normalize(reflectdir), glm::normalize(light.position - point.hitPoint)));
				lights.push_back(l);
			}
		}
	}

	return lights;
}


std::vector<Lighting> getPlaneLights(const HitInfo& point, const glm::vec3& reflectdir, const Scene& scene, const BoundingVolumeHierarchy& bvh, int rayCount1D) {
	std::vector<Lighting> lights;

	for (const PlaneLight& light : scene.planeLight) {
		float hit = 0;
		float maxCosAngle = 0; // angle for specularity so that the specular highlight will be square
		glm::vec3 dx = (1.0f / (rayCount1D-1)) * (light.width); // split the width and hight so that i can cast multype rays
		glm::vec3 dy = (1.0f / (rayCount1D-1)) * (light.height);
		glm::vec3 py = light.position;
		glm::vec3 normal = glm::normalize(glm::cross(light.width, light.height));

		// optimization so only points in front of the light get checked whith raytacing
		if (glm::dot(glm::normalize(point.hitPoint - (light.position+0.5f*(light.width+light.height))), normal) > 0) {
			// loop over points on the light
			for (int i = 0; i < rayCount1D; i++) {
				glm::vec3 px = py;
				for (int j = 0; j < rayCount1D; j++) {
					// if light ray not obstructed
					if (cansee(point.hitPoint, px, scene, bvh)) {
						// becouse of this(the hit calculation) you dont later need to multyply by the cosine of the angle(normal and light dir)
						hit += std::max(glm::dot(glm::normalize(point.hitPoint - px), normal), 0.0f) / glm::length(point.hitPoint - px);
						// find the smalledst angle to the light, so max cosine
						maxCosAngle = std::max(maxCosAngle, glm::dot(glm::normalize(reflectdir), glm::normalize(px - point.hitPoint)));
					}
					px += dx;
				}
				py += dy;
			}
		}
		if (hit > 0) {
			// setting all the properites of this light
			Lighting l;
			l.color = light.color;
			l.position = light.position;
			l.intensity = hit / (float)(rayCount1D * rayCount1D); // this inclues the cosLightSurfaceAngle
			l.cosLightSurfaceAngle = 1; 
			l.cosLightSpecAngle = maxCosAngle;
			lights.push_back(l);
		}
	}

	return lights;
}