#include "ray_tracing.h"
#include "draw.h"
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

// returns true if the floats are equal, otherwise false
bool isEqual(float x, float y)
{
    return glm::abs(x - y) < 1e-4;
}

// returns true if the floa is zero, otherwise false
bool isZero(float x)
{
    return isEqual(x, 0);
}

// returns true if x >= y, otherwise false
bool greaterThanOrEqual(float x, float y)
{
    return isEqual(x, y) || x >= y;
}

// returns true if x <= y, otherwise false
bool smallerThanOrEqual(float x, float y)
{
    return greaterThanOrEqual(y, x);
}
float getSign(const glm::vec3 &v, const glm::vec3 &n)
{
    return glm::dot(glm::normalize(v), n);
}

bool pointInTriangle(const glm::vec3 &v0, const glm::vec3 &v1, const glm::vec3 &v2, const glm::vec3 &n, const glm::vec3 &p)
{
    //Check sings of cross products for each side of the triangle
    bool s0 = glm::dot(glm::cross((p - v0), (v2 - v0)), n) >= 0;
    bool s1 = glm::dot(glm::cross((p - v2), (v1 - v2)), n) >= 0;
    bool s2 = glm::dot(glm::cross((p - v1), (v0 - v1)), n) >= 0;

    //If all positive then point is inside the triangle
    if (s0 && s1 && s2)
    {
        return true;
    }
    //If all negative then point is inside the triangle
    else if (!s0 && !s1 && !s2)
    {
        return true;
    }
    //Point is outisde the triangle
    return false;
}

bool intersectRayWithPlane(const Plane &plane, Ray &ray)
{
    float nDotd = glm::dot(glm::normalize(ray.direction), plane.normal);
    //Check that ray is not parallel to the plane of the triangle
    if (nDotd == 0)
    {
        return false;
    }
    float t = (plane.D - glm::dot(ray.origin, plane.normal)) / nDotd;
    if (t >= 0)
    {
        //Check if found t is actually closer to the origin than the current t stored in ray.t
        if (t < ray.t)
        {
            ray.t = t;
            return true;
        }
        else
        {
            return false;
        }
    }
    else
    {
        return false;
    }
}

Plane trianglePlane(const glm::vec3 &v0, const glm::vec3 &v1, const glm::vec3 &v2)
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
bool intersectRayWithTriangle(const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& v2, Ray& ray, HitInfo& hitInfo, int material_index)
{
    //Store old t
    float tOld = ray.t;
    Plane plane = trianglePlane(v0, v1, v2); //get the plane
    if (intersectRayWithPlane(plane, ray))
    { //does ray intersect plane
        glm::vec3 point = ray.origin + ray.direction * ray.t;
        if (pointInTriangle(v0, v1, v2, plane.normal, point))
        { //does point lie in triangle
            hitInfo.normal = plane.normal;
            hitInfo.hitPoint = point;
            hitInfo.material_index = material_index;


            return true; //return found point
        }
        else
        {
            //Restore old t if no hit
            ray.t = tOld;
        }
    }
    return false;
}

/// Input: the three vertices of the triangle
/// Output: if intersects then modify the hit parameter ray.t and return true, otherwise return false.
/// In addition to the method 'intersectRayWithTriangle' it also interpolates the normals and the texture coordinates of the vertices
bool intersectRayWithTriangleWithInterpolation(const Vertex& v0, const Vertex& v1, const Vertex& v2, Ray& ray, HitInfo& hitInfo, int material_index) {
    
    float old_t = ray.t;
    
    // if the ray intersects the triangle, interpolate the normals and the texture coordinates
    if(intersectRayWithTriangle(v0.p, v1.p, v2.p, ray, hitInfo, material_index)) {

        // if the object was actually hit, update the hitInfo with the corresponding triangle

        if(ray.t < old_t) {
            hitInfo.is_triangle = true;
            hitInfo.intersected_triangle = std::array<Vertex, 3> {v0, v1, v2};
        }

        glm::vec3 barCoords;
        // should always find barycentric coordinates, since the hitPoint is inside the triangle
        barycentricCoordinates(v0.p, v1.p, v2.p, hitInfo.hitPoint, barCoords);

        // INTERPOLATION OF THE NORMALS
        glm::vec3 facenormal = hitInfo.normal;
        //drawRay({ v0.p, v0.n, 0.1f }, glm::vec3(0, 0, 1));
        //drawRay({ v1.p, v1.n, 0.1f }, glm::vec3(0, 0, 1));
        //drawRay({ v2.p, v2.n, 0.1f }, glm::vec3(0, 0, 1));
        std::array<glm::vec3, 3> normals {v0.n, v1.n, v2.n};
        hitInfo.normal = interpolateProperty(barCoords, normals);
        if (glm::dot(hitInfo.normal, facenormal) < 0) {
            hitInfo.normal = -hitInfo.normal;
        }
        //drawRay({ hitInfo.hitPoint, hitInfo.normal, 0.2f }, glm::vec3(0, 1, 1));
        // update the 



        // INTERPOLATION FOR THE TEXTURE COORDINATES
        // Warning: Even if the triangle does not have a texture, the texture coordinates of the provided vertices will still be interpolated
        std::array<glm::vec2, 3> textCoords {v0.texCoord, v1.texCoord, v2.texCoord};
        hitInfo.texCoord = interpolateProperty(barCoords, textCoords);

        //std::cout << barCoords.x << " " << barCoords.y << std::endl;

        return true;
        
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
        if (t0 < 0) { t0 = t1; }
        if (t1 < 0) { t1 = t0; }
        if (glm::min(t0, t1) > 0 && glm::min(t0, t1)<ray.t) {
            ray.t = glm::min(t0, t1);
            hitInfo.hitPoint = ray.origin + ray.t * ray.direction;
            hitInfo.normal = glm::normalize(hitInfo.hitPoint - sphere.center);
            hitInfo.sphere_material = sphere.material;
            hitInfo.is_triangle = false;
            hitInfo.intersected_sphere = sphere;
            return true;
        }
    }
    return false;
}

/// Input: an axis-aligned bounding box with the following parameters: minimum coordinates box.lower and maximum coordinates box.upper
/// Output: if intersects then modify the hit parameter ray.t and return true, otherwise return false
bool intersectRayWithShape(const AxisAlignedBox& box, Ray& ray)
{
    if (glm::all(glm::equal(box.lower, glm::vec3(std::numeric_limits<float>::max()))) && glm::all(glm::equal(box.upper, glm::vec3(-std::numeric_limits<float>::max())))) {
        return false;
    }
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
            float min = std::numeric_limits<float>::max();
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

// computes the area of the parallelogram defined by the verticies v0, v1 and v2
// the fourth  vertex (v3) is not specified explicitly, but conceptually you
// can think of it as being placed such that v0, v1, v2, v3 form a parallelogram
float parallelogramArea(const glm::vec3 &v0, const glm::vec3 &v1, const glm::vec3 &v2) {
    return glm::length(glm::cross(v1 - v0, v2 - v0));
}

// If the point p is in the triangle formed by the vertices v0, v1, v2, the function sets the coords argument to the barycentric
// coordinates of p, in terms of v0, v1, v2 and returns true
// else returns false
bool barycentricCoordinates(const glm::vec3 &v0, const glm::vec3 &v1, const glm::vec3 &v2, const glm::vec3 &p, glm::vec3 &coords) {
    Plane trPlane = trianglePlane(v0, v1, v2);


    // if the point is not in the triangles plane
    if(!isZero(glm::dot(trPlane.normal, p - v0)))
        return false;

   
    
    // if the point is not in the triangle, return false
    if(!pointInTriangle(v0, v1, v2, trPlane.normal, p))
        return false;

    // compute the barycentric coordinates
    float totalArea = parallelogramArea(v0, v1, v2);

    //the shape is not a triangle, so the barycentric coordinates do not make sense
    if(isZero(totalArea))
        return false;

    float c0 = parallelogramArea(p, v1, v2) / totalArea;
    float c1 = parallelogramArea(v0, p, v2) / totalArea;
    float c2 = parallelogramArea(v0, v1, p) / totalArea;

    coords.x = c0;
    coords.y = c1;
    coords.z = c2;
    
    

    return true;
}

// interpolates a property of type P, over the provided barycentric coordinates
// barCoords - the barycentric coordinates
// properties - an array of size 3 specifying the property at each vertex of the barycentric coordinate system
template <typename P>
P interpolateProperty(const glm::vec3 &barCoords, const std::array<P, 3> &properties) {
    return properties[0] * barCoords.x + properties[1] * barCoords.y + properties[2] * barCoords.z;
}

