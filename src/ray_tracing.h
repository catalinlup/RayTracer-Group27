
#pragma once
#include "scene.h"
#include <optional>

struct HitInfo {
    glm::vec3 normal;
    glm::vec3 hitPoint;
    int material_index; // stores the index of the mesh in the scene that contains this material

    Material sphere_material; // spheres don't have textures, so it is fine to store the material

    // the interpolated texture coordinate, based on the 3 vertices defining the plane
    glm::vec2 texCoord;

    // the nature of the intersected object is necessary for the ray differentials update
    bool is_triangle = false; // true if the intersected object is a triangle, false if it is a sphere
    std::array<Vertex, 3> intersected_triangle; // the intersected triangle if the ray intersected a triangle
    Sphere intersected_sphere; // the intersected sphere if the intersected object is a sphere.

    // returns a reference to the material
    // fast, but the material is lost once the object gets out of scope
    Material& getMaterial(Scene& scene) {
        if(is_triangle)
            return scene.meshes[material_index].material;
        return sphere_material;
    }


    // the material is not lost, but it is very slow
    Material getMaterialCopy(Scene& scene) {
        if(is_triangle)
            return scene.meshes[material_index].material;
        return sphere_material;
    }

};

bool intersectRayWithPlane(const Plane& plane, Ray& ray);

// Returns true if the point p is inside the triangle spanned by v0, v1, v2 with normal n.
bool pointInTriangle(const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& v2, const glm::vec3& n, const glm::vec3& p);

Plane trianglePlane(const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& v2);

bool intersectRayWithTriangle(const glm::vec3 &v0, const glm::vec3 &v1, const glm::vec3 &v2, Ray &ray, HitInfo &hitInfo, int material_index);

/// Input: the three vertices of the triangle
/// Output: if intersects then modify the hit parameter ray.t and return true, otherwise return false.
/// In addition to the method 'intersectRayWithTriangle' it also interpolates the normals and the texture coordinates of the vertices
bool intersectRayWithTriangleWithInterpolation(const Vertex &v0, const Vertex &v1, const Vertex &v2, Ray &ray, HitInfo &hitInfo, int material_index);

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

