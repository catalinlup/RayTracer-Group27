#include "ray_tracing.h"
#include "disable_all_warnings.h"
// Suppress warnings in third-party code.
DISABLE_WARNINGS_PUSH()
#include <glm/geometric.hpp>
#include <glm/gtx/component_wise.hpp>
#include <glm/vector_relational.hpp>
DISABLE_WARNINGS_POP()
#include <cmath>
#include <iostream>
#include <limits>

bool pointInTriangle(const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& v2, const glm::vec3& n, const glm::vec3& p)
{
    //Check sings of cross products for each side of the triangle
    bool s0 = glm::dot(glm::cross((p - v0), (v2 - v0)), n) >= 0;
    bool s1 = glm::dot(glm::cross((p - v2), (v1 - v2)), n) >= 0;
    bool s2 = glm::dot(glm::cross((p - v1), (v0 - v1)), n) >= 0;

    //If all positive then point is inside the triangle
    if (s0 && s1 && s2) {
        return true;
    }
    //If all negative then point is inside the triangle
    else if (!s0 && !s1 && !s2) {
        return true;
    }
    //Point is outisde the triangle
    return false;
}

bool intersectRayWithPlane(const Plane& plane, Ray& ray)
{
    float nDotd = glm::dot(glm::normalize(ray.direction), plane.normal);
    //Check that ray is not parallel to the plane of the triangle
    if (nDotd == 0) {
        return false;
    }
    float t = (plane.D - glm::dot(ray.origin, plane.normal)) / nDotd;
    if (t >= 0) {
        //Check if found t is actually closer to the origin than the current t stored in ray.t
        if (t < ray.t) {
            ray.t = t;
            return true;
        }
        else {
            return false;
        }
    }
    else {
        return false;
    }
}

Plane trianglePlane(const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& v2)
{
    Plane plane;
    glm::vec3 n = glm::normalize(glm::cross((v0 - v2), (v1 - v2)));
    float distance = glm::dot(n, v0);
    plane.D = distance;
    plane.normal = n;

    return plane;
}

/// Input: the three vertices of the triangle
/// Output: if intersects then modify the hit parameter ray.t and return true, otherwise return false
bool intersectRayWithTriangle(const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& v2, Ray& ray, HitInfo& hitInfo)
{
    //Store old t
    float tOld = ray.t;
    Plane plane = trianglePlane(v0, v1, v2); //get the plane
    if (intersectRayWithPlane(plane, ray)) { //does ray intersect plane
        glm::vec3 point = ray.origin + ray.direction * ray.t;
        if (pointInTriangle(v0, v1, v2, plane.normal, point)) { //does point lie in triangle
            return true; //return found point
        }
        else {
            //Restore old t if no hit
            ray.t = tOld;
        }
    }
    return false;
}

/// Input: a sphere with the following attributes: sphere.radius, sphere.center
/// Output: if intersects then modify the hit parameter ray.t and return true, otherwise return false
bool intersectRayWithShape(const Sphere& sphere, Ray& ray, HitInfo& hitInfo)
{
    //Translate the ray by the sphere center. Considering sphere to be centered around (0,0,0)
    glm::vec3 rayMoved = ray.origin - sphere.center;

    float A = glm::pow(ray.direction.x, 2) + glm::pow(ray.direction.y, 2) + glm::pow(ray.direction.z, 2);
    float B = 2 * (ray.direction.x * rayMoved.x + ray.direction.y * rayMoved.y + ray.direction.z * rayMoved.z);
    float C = glm::pow(rayMoved.x, 2) + glm::pow(rayMoved.y, 2) + glm::pow(rayMoved.z, 2) - glm::pow(sphere.radius, 2);

    float disc = glm::pow(B, 2) - 4 * A * C;

    if (disc >= 0) {
        float t0 = (-B + glm::sqrt(disc)) / (2 * A);
        float t1 = (-B - glm::sqrt(disc)) / (2 * A);
        ray.t = glm::min(t0, t1);
        return true;
    }
    return false;
}

/// Input: an axis-aligned bounding box with the following parameters: minimum coordinates box.lower and maximum coordinates box.upper
/// Output: if intersects then modify the hit parameter ray.t and return true, otherwise return false
bool intersectRayWithShape(const AxisAlignedBox& box, Ray& ray)
{
    glm::vec3 dir = glm::normalize(ray.direction);

    float txmin = (box.lower.x - ray.origin.x) / dir.x;
    float txmax = (box.upper.x - ray.origin.x) / dir.x;

    float tymin = (box.lower.y - ray.origin.y) / dir.y;
    float tymax = (box.upper.y - ray.origin.y) / dir.y;

    float tzmin = (box.lower.z - ray.origin.z) / dir.z;
    float tzmax = (box.upper.z - ray.origin.z) / dir.z;

    if (glm::all(glm::greaterThan(ray.origin, box.lower))) {
        if (glm::all(glm::lessThan(ray.origin, box.upper))) {
            std::vector<float> tValues = { txmin, txmax, tymin, tymax, tzmin, tzmax };
            float min = FLT_MAX;
            for (float v : tValues) {
                if (v > 0 && v < min) {
                    min = v;
                }
            }
            ray.t = min;
            return true;
        }
    }

    //Get min values for tin
    float tinx = std::min(txmin, txmax);
    float tiny = std::min(tymin, tymax);
    float tinz = std::min(tzmin, tzmax);
    //Get max values for tout
    float toutx = std::max(txmin, txmax);
    float touty = std::max(tymin, tymax);
    float toutz = std::max(tzmin, tzmax);

    //Tin gets the the value after all tins have been crossed [min]
    //Tout gets the first value for tout that is found [max]
    float tin = glm::max(glm::max(tinx, tiny), tinz);
    float tout = glm::min(glm::min(toutx, touty), toutz);

    //No hit
    if (tin > tout || tout < 0) {
        return false;
    }
    //Hit
    ray.t = tin;
    return true;
}
