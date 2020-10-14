#pragma once
#include "ray_tracing.h"
#include "scene.h"
#include <map>
#include <climits>
#include <algorithm>
#include <array>
#include <exception>
#include <utility>
#include <queue>
#include <iostream>

// It would have made a lot more sense to have a BvhSphere class and a BvhTriangle class extend from BvhObject,
// but unfortunetely, polymorphism cannot be achieved without pointers :(

// Encapsulates an atomic object of Bounding Volume Hierarchy.
// Can be either a triangle or a sphere, but has the potential to be extended to other shapes.

class BvhObject
{

    static unsigned long long id;

public:
    BvhObject();

    // constructs an object from a triangle
    BvhObject(glm::vec3 &v0, glm::vec3 &v1, glm::vec3 &v2);

    // constructs an object of a sphere
    BvhObject(Sphere &sphere);

    // returns the AABB corresponding to the object
    AxisAlignedBox &getAABB();

    // returns a string with the object type, either 'triangle' or 'sphere'
    std::string getType();

    // returns the object's unique id
    unsigned long long getId();

    std::array<glm::vec3, 3> getTriangle();
    Sphere getSphere();

    // after calling the this function, the ids of the further declared objects will start again at 0.
    static void resetIdGenerator();

private:
    // default constructor, called by all of the public constructor

    // the type of the object, either triangle or sphere
    std::string _type;

    const std::string _SPHERE_TYPE = "sphere";
    const std::string _TRIANGLE_TYPE = "triangle";

    // unique id for the object.
    unsigned long long _id;

    AxisAlignedBox _boundingBox;

    // use 'type' to distiguish between triangle and sphere
    std::array<glm::vec3, 3> _triangle; // stores the triangle, in case this object wraps around a triangle
    std::array<Sphere, 1> _sphere;      
};

struct Node
{

    // ctr
    // takes a vector with all the BvhObjects in it's subtree. Used to construct the AABB
    // isLeaf should be true if this is a leaf node, false otherwise.
    Node(std::vector<BvhObject> &containingBvhObjectss, bool isLeaf);

    Node();

    // add child Node, only for non-leaves
    void addChild(Node& other);

    // add child BvhObject, only for leaves
    void addChild(BvhObject &object);

    // returns true if the node is a leaf, otherwise false
    bool isLeaf();

    // returns the id of the node
    unsigned long long getId();

    // return the bounding box corresponding to this node
    AxisAlignedBox getBoundingBox();

    // returns the ids of the children of the node.
    std::vector<unsigned long long> getChildren();




    // after calling the this function, the ids of the further declared nodes will start again at 0.
    static void resetIdGenerator();

private:
    void updateAABB(AxisAlignedBox& other);

    static unsigned long long id;

    bool _isLeaf;

    std::vector<unsigned long long> _children; // contains the ids of the children of this node. These ids are either node ids, or bvhObject ids, in the case of leaf nodes

    unsigned long long _id;

    AxisAlignedBox _boundingBox;
};




class BoundingVolumeHierarchy {
public:
    BoundingVolumeHierarchy(Scene* pScene);

    // Use this function to visualize your BVH. This can be useful for debugging.
    void debugDraw(int level, bool showLeafNodes);
    int numLevels() const;

    // Return true if something is hit, returns false otherwise.
    // Only find hits if they are closer than t stored in the ray and the intersection
    // is on the correct side of the origin (the new t >= 0).
    bool intersect(Ray& ray, HitInfo& hitInfo) const;

    // loads the data from the scene as BvhObjects
    void addBvhObjectsFromScene();

    // bvhObjectIds - the ids of the objects that are port of this node's subtree
    // level - the level of this node
    // isLeaf - true if the node should be a leaf node, false otherwise
    unsigned long long createNode(std::vector<unsigned long long> &bvhObjectIds, unsigned int level);

    private :

        // prints a textual representation of the hierarchy tree. Useful for debugging
        void printHierarchy();

        unsigned int _max_level = 7;

        unsigned int _max_level_achieved = 0;

        // the root node
        unsigned _root_id;

        // performs a split of the provided objects, based on the median value
        // split_axis = 0 - for x_axis; 1 - for y_axis; 2 - for z_axis
        std::array<std::vector<unsigned long long>, 2> medianSplit(std::vector<unsigned long long> bvhObjectIds, unsigned int split_axis);

        Scene *m_pScene;

        std::vector<unsigned long long> leafNodes;                                     // stores the ids of all of the leaf nodes
        std::map<unsigned long long, Node> nodes;                                      // map storing all of the nodes in the BVH indexed by their id.
        std::map<unsigned long long, BvhObject> bvhObjects;                            // map storing all BvhObjects indexed by their id
        std::map<unsigned long long, std::vector<unsigned long long>> objects_at_node; // map storing the ids of all the objects corresponding to one Node
        std::map<unsigned int, std::vector<unsigned long long>> nodes_at_level;        // stores the ids of all of the nodes located at a certain level

        // Utility functions to get all keys and all values from a map
        template <typename K, typename V>
        std::vector<K> getAllMapKeys(std::map<K, V> mp);

        template <typename K, typename V>
        std::vector<V> getAllMapValues(std::map<K, V> mp);

        template <typename K, typename V>
        std::vector<V> getMapValuesOfKeys(std::map<K, V> mp, std::vector<K> keys);
    };
