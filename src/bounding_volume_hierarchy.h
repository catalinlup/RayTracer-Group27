#pragma once
#include "ray_tracing.h"
#include "scene.h"
#include <array>
#include <gsl-lite/gsl-lite.hpp>
#include <iostream>
#include <queue>
#include <map>

//NOTE:
// is_triangle = true - means index of triangle object
// is_triangle = false - means index of sphere object

struct Node {
    bool isLeaf = false;
    std::vector<int> children;
    std::vector<bool> is_triangle;
    AxisAlignedBox AABB;
    
};

class BoundingVolumeHierarchy {
public:
    BoundingVolumeHierarchy(Scene* pScene);

    // Use this function to visualize your BVH. This can be useful for debugging.
    void debugDraw(int level);
    int numLevels() const;

    // Return true if something is hit, returns false otherwise.
    // Only find hits if they are closer than t stored in the ray and the intersection
    // is on the correct side of the origin (the new t >= 0).
    bool intersect(Ray &ray, HitInfo &hitInfo, bool useBVH) const;

private:

    void constructBVH();

    int getRootIndex() const;

    void createNodeAndUpdateStats(std::vector<int> &object_indices, std::vector<bool> &is_triangle, int level);
    Node createNodeFromObjects(std::vector<int> &object_indices, std::vector<bool> &is_triangle, int level) const;

    AxisAlignedBox createAabbFromObjects(std::vector<int>& object_indices, std::vector<bool>& is_triangle) const;


    void sortObjects(std::vector<int>& object_indices, std::vector<bool>& is_triangle, int level) const;

    float getSortingAttribute(int object_index, bool is_triangle, int level) const;
    float getSortingAttributeTriangle(int triangle_index, int level) const;
    float getSortingAttributeSphere(int sphere_index, int level) const;

    bool validateObjectVectorPair(std::vector<int> &object_indices, std::vector<bool> &is_triangle) const;
    bool validateSphereIndex(int index) const;
    bool validateTriangleIndex(int index) const;

    void loadObjectsFromScene();

    bool intersectNode(int node_index, Ray& ray) const;

    bool intersectObject(int object_index, bool is_triangle, Ray& ray, HitInfo& hitInfo) const;

    bool intersectBVH(int node_index, Ray& ray, HitInfo& hitInfo) const;

  

    int max_level = 4; // inclusive (indexed from 0)
    int max_level_achieved = 0; // inclusive (indexed from 0)

    Scene *m_pScene;

    std::vector<Node> nodes;
    std::map<int, std::vector<int>> node_indices_per_level;
    
    std::map<int, std::vector<int>> object_indices_per_node;
    std::map<int, std::vector<bool>> is_triangle_stats_per_node;

    std::vector<Sphere> spheres;
    std::vector<std::array<Vertex, 3>> triangles;
    std::vector<Material> triangle_materials;
};
