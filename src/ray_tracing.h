#pragma once
#include "scene.h"
#include <optional>

struct HitInfo {
    glm::vec3 normal;
    glm::vec3 hitPoint;
    Material material;

    // the interpolated texture coordinate, based on the 3 vertices defining the plane
    glm::vec2 texCoord;
};

bool intersectRayWithPlane(const Plane& plane, Ray& ray);

// Returns true if the point p is inside the triangle spanned by v0, v1, v2 with normal n.
bool pointInTriangle(const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& v2, const glm::vec3& n, const glm::vec3& p);

Plane trianglePlane(const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& v2);

bool intersectRayWithTriangle(const glm::vec3 &v0, const glm::vec3 &v1, const glm::vec3 &v2, Ray &ray, HitInfo &hitInfo, Material &m);

/// Input: the three vertices of the triangle
/// Output: if intersects then modify the hit parameter ray.t and return true, otherwise return false.
/// In addition to the method 'intersectRayWithTriangle' it also interpolates the normals and the texture coordinates of the vertices
bool intersectRayWithTriangleWithInterpolation(const Vertex &v0, const Vertex &v1, const Vertex &v2, Ray &ray, HitInfo &hitInfo, Material &m);

bool intersectRayWithShape(const Sphere& sphere, Ray& ray, HitInfo& hitInfo);
bool intersectRayWithShape(const AxisAlignedBox& box, Ray& ray);

// If the point p is in the triangle formed by the vertices v0, v1, v2, sets the coords argument to the barycentric
// coordinates of p, in terms of v0, v1, v2 and returns true
// else returns false
bool barycentricCoordinates(const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& v2, const glm::vec3& p, glm::vec3& coords);


// interpolates a property of type P, over the provided barycentric coordinates
// barCoords - the barycentric coordinates
// properties - an array of size 3 specifying the property at each vertex of the barycentric coordinate system
template <typename P>
P interpolateProperty(const glm::vec3& barCoords, const std::array<P, 3>& properties);