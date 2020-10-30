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

bool cansee(const glm::vec3& p1, const glm::vec3& p2, float& intensity, Scene& scene, const BoundingVolumeHierarchy& bvh) {
	glm::vec3 o = p1;
	glm::vec3 d = p2 - p1;
	float distance = glm::length(d);
	d = glm::normalize(d);
	o += SHADOW_ERROR_OFFSET * d;
	Ray ray = { o, d };
	HitInfo hitInfo;

	while (distance > SHADOW_ERROR_OFFSET) {
		bool hit = bvh.intersect(ray, hitInfo, true);
		Material& matForRendering = hitInfo.getMaterial(scene);
		if (!hit || (ray.t > distance - 2 * SHADOW_ERROR_OFFSET)) {
			drawRay(ray, glm::vec3(0.5f, 1.0f, 0.5f));
			return true;
		}
		else if (matForRendering.transparency!=1.0f) {
			drawRay(ray, glm::vec3(0.3f, 0.3f, 8.0f));
			distance -= ray.t;
			ray.t = std::numeric_limits<float>::max();
			ray.origin = hitInfo.hitPoint + SHADOW_ERROR_OFFSET * ray.direction;

			float c = std::abs(glm::dot(ray.direction, hitInfo.normal));
			// calculate how how much light should reflect and how much should refract
			float& R0 = matForRendering.transparency;//R0static; // probebility of reflection at 0 rad/deg to normal
				// simplified fersnel equasion : refraction
			intensity *= 1 - (R0 + (1 - R0) * (std::pow(1 - c, 5)));

		}
		else {
			drawRay(ray, glm::vec3(0.8f, 0.2f, 0.0f));
			return false;
		}
	}
	return true;
	

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
std::vector<Lighting> getPointLights(const HitInfo& point, const glm::vec3 reflectdir, Scene& scene, const BoundingVolumeHierarchy& bvh) {

	std::vector<Lighting> visibleLights;

	// check if the point can see any light
	for (const PointLight& light : scene.pointLights) {
		float intensity = 1.0f;


		//std::cout << "norm" << point.normal.x << ", " << point.normal.y << ", " << point.normal.z << " len " << glm::length(point.normal) << std::endl;
		//drawRay({ point.hitPoint, point.normal, 1 }, {0,0,1});

		// the ray hit nothing or the hit was only after the light
		if (cansee(point.hitPoint, light.position, intensity, scene, bvh)) {
			// store the light scources that we can see
			Lighting l;
			l.color = light.color;
			l.position = light.position;
			l.intensity = intensity;
			l.cosLightSurfaceAngle = glm::abs((glm::dot(glm::normalize(point.normal), glm::normalize(light.position - point.hitPoint))));
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
std::vector<Lighting> getSpherelights(const HitInfo& point, const glm::vec3& reflectdir, Scene& scene, const BoundingVolumeHierarchy& bvh, int rayCount) {
	std::vector<Lighting> lights;
	for (const SphericalLight& light : scene.sphericalLight) {

		float intensity = 1; // used for shadow of transparent objects
		float intensitySum = 1;
		int hits = 0; // the number of rays that hit the light

		// check if the center of the light is visable
		if (cansee(point.hitPoint, light.position, intensitySum, scene, bvh)) {
			hits++;
		}

		// calculate direction to the center of the light ray
		glm::vec3 d = light.position - point.hitPoint;
		float distance = glm::length(d);
		d = glm::normalize(d);

		// a vector that is not in line whith d
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
#ifdef USE_MONTE_CARLO_INTEGRATION
		glm::vec3 perp2 = glm::normalize(glm::cross(d, perp)) * light.radius;
		float a = 1;
		float b = 1;
		for (int i = 0; i < rayCount - 1; i++) {
			intensity = 1;
			// get a random point on the light source
			do {
				a = (rand() / (float)RAND_MAX) * 2 - 1;
				b = (rand() / (float)RAND_MAX) * 2 - 1;
			} while (a * a + b * b > 1);
			// see if the random point is visable
			if (cansee(point.hitPoint, light.position + a * perp + b * perp2, intensity, scene, bvh)) {
				hits++;
				intensitySum += intensity;
			}
		}
#else
		
		// m: how meny rigns, n how meny rays per ring
		int m = std::max(1, (int)(rayCount/std::round(std::sqrt(2 * 3.14159365358979f * rayCount))));
		int n = (rayCount-1) / m;
		rayCount = m * n + 1;

		// thare are n rayes that go to the edge of the circle
		glm::mat3 rotate = rotatetionMatrix(2 * 3.14159365358979f /n, d); // rotation matrix around the vector d

		// loop over all the rays that go the the edge of the circle
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < m; j++) {
				intensity = 1;
				if (cansee(point.hitPoint, light.position + ((m-j)/(float)m)*perp, intensity, scene, bvh)) {
					hits++;
					intensitySum += intensity;
				}
				// rotate the perpenicular vector to be able to calculate the next position on the light where a ray must go
			}
			perp = rotate * perp;
		}
#endif
		

		// save the visable lights
		if (hits > 0) {
			// setting all the properites of this light
			Lighting l;
			l.color = light.color;
			l.intensity = intensitySum / (float)(rayCount);
			l.position = light.position;
			l.cosLightSurfaceAngle = std::abs(glm::dot(glm::normalize(point.normal), glm::normalize(light.position - point.hitPoint)));
			l.cosLightSpecAngle = std::max(0.0f, glm::dot(glm::normalize(reflectdir), glm::normalize(light.position - point.hitPoint)));
			lights.push_back(l);

		}
	}
	return lights;
}

// returns a list of all spotlights that iluminate a point
std::vector<Lighting> getSpotLichts(const HitInfo& point, const glm::vec3& reflectdir, Scene& scene, const BoundingVolumeHierarchy& bvh) {
	std::vector<Lighting> lights;

	for (const SpotLight& light : scene.spotLight) {

		// if the point is in the beam of the spot light
		if (glm::dot(glm::normalize(light.direction), glm::normalize(point.hitPoint - light.position)) > std::cos(glm::radians(light.angle))) {
			float intensity = 1;
			// check if the path is not obstructed
			if (cansee(point.hitPoint, light.position, intensity, scene, bvh)) {
				// setting all the properites of this light
				Lighting l;
				l.color = light.color;
				l.position = light.position;
				l.intensity = intensity;
				l.cosLightSurfaceAngle = std::abs(glm::dot(glm::normalize(point.normal), glm::normalize(light.position - point.hitPoint)));
				l.cosLightSpecAngle = std::max(0.0f, glm::dot(glm::normalize(reflectdir), glm::normalize(light.position - point.hitPoint)));
				lights.push_back(l);
			}
		}
	}

	return lights;
}


std::vector<Lighting> getPlaneLights(const HitInfo& point, const glm::vec3& reflectdir, Scene& scene, const BoundingVolumeHierarchy& bvh, int rayCount1D) {
	std::vector<Lighting> lights;

	for (const PlaneLight& light : scene.planeLight) {
		float hit = 0;
		int hitCount = 0;
		float maxCosAngle = 0; // angle for specularity so that the specular highlight will be square
		float intensitySum = 0; // used for shadow of transparent objects
		float intensity = 1;
		glm::vec3 dx = (1.0f / (rayCount1D-1)) * (light.width); // split the width and hight so that i can cast multype rays
		glm::vec3 dy = (1.0f / (rayCount1D-1)) * (light.height);
		glm::vec3 py = light.position;
		glm::vec3 normal = glm::normalize(glm::cross(light.width, light.height));

		// optimization so only points in front of the light get checked whith raytacing
		if (glm::dot(glm::normalize(point.hitPoint - (light.position+0.5f*(light.width+light.height))), normal) > 0) {
#ifdef USE_MONTE_CARLO_INTEGRATION
			for (int i = 0; i < rayCount1D * rayCount1D; i++) {
				float a = (rand() / (float)RAND_MAX);
				float b = (rand() / (float)RAND_MAX);
				intensity = 1;
				// to random point on the light source to sample
				glm::vec3 px = light.position + a * light.width + b * light.height;
				if (cansee(point.hitPoint, px, intensity, scene, bvh)) {
					intensitySum += intensity;
					// becouse of this(the hit calculation) you dont later need to multyply by the cosine of the angle(normal and light dir)
					hit += std::max(glm::dot(glm::normalize(point.hitPoint - px), normal), 0.0f) / glm::length(point.hitPoint - px);
					hitCount++;
					// find the smalledst angle to the light, so max cosine
					maxCosAngle = std::max(maxCosAngle, glm::dot(glm::normalize(reflectdir), glm::normalize(px - point.hitPoint)));
				}
			}
#else
			// loop over points on the light
			for (int i = 0; i < rayCount1D; i++) {
				glm::vec3 px = py;
				for (int j = 0; j < rayCount1D; j++) {
					// if light ray not obstructed
					intensity = 1;
					if (cansee(point.hitPoint, px, intensity, scene, bvh)) {
						intensitySum += intensity;
						// becouse of this(the hit calculation) you dont later need to multyply by the cosine of the angle(normal and light dir)
						hit += std::max(glm::dot(glm::normalize(point.hitPoint - px), normal), 0.0f) / glm::length(point.hitPoint - px);
						hitCount++;
						// find the smalledst angle to the light, so max cosine
						maxCosAngle = std::max(maxCosAngle, glm::dot(glm::normalize(reflectdir), glm::normalize(px - point.hitPoint)));
					}
					px += dx;
				}
				py += dy;
			}
#endif
		}
		if (hit > 0) {
			// setting all the properites of this light
			Lighting l;
			l.color = light.color;
			l.position = light.position;
			l.intensity = (intensitySum/hitCount) * hit / (float)(rayCount1D * rayCount1D); // this inclues the cosLightSurfaceAngle
			l.cosLightSurfaceAngle = 1; 
			l.cosLightSpecAngle = maxCosAngle;
			lights.push_back(l);
		}
	}

	return lights;
}