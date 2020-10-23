#include "bounding_volume_hierarchy.h"
#include "draw.h"


BoundingVolumeHierarchy::BoundingVolumeHierarchy(Scene* pScene)
    : m_pScene(pScene)
{
    constructBVH();
}

// Use this function to visualize your BVH. This can be useful for debugging. Use the functions in
// draw.h to draw the various shapes. We have extended the AABB draw functions to support wireframe
// mode, arbitrary colors and transparency.
void BoundingVolumeHierarchy::debugDraw(int level)
{

    // Draw the AABB as a transparent green box.
    //AxisAlignedBox aabb{ glm::vec3(-0.05f), glm::vec3(0.05f, 1.05f, 1.05f) };
    //drawShape(aabb, DrawMode::Filled, glm::vec3(0.0f, 1.0f, 0.0f), 0.2f);

    // Draw the AABB as a (white) wireframe box.
    // AxisAlignedBox aabb { glm::vec3(-0.05f), glm::vec3(0.05f, 1.05f, 1.05f) };
    // //drawAABB(aabb, DrawMode::Wireframe);
    // drawAABB(aabb, DrawMode::Filled, glm::vec3(0.05f, 1.0f, 0.05f), 0.1f);

    for(const auto& index : node_indices_per_level[level]) {
        Node n = nodes[index];
        drawAABB(n.AABB, DrawMode::Wireframe, glm::vec3(1.0f, 0.0f, 0.0f));
    }
}

int BoundingVolumeHierarchy::numLevels() const
{
    return max_level_achieved + 1;
}

// Return true if something is hit, returns false otherwise. Only find hits if they are closer than t stored
// in the ray and if the intersection is on the correct side of the origin (the new t >= 0). Replace the code
// by a bounding volume hierarchy acceleration structure as described in the assignment. You can change any
// file you like, including bounding_volume_hierarchy.h .
bool BoundingVolumeHierarchy::intersect(Ray& ray, HitInfo& hitInfo) const
{
    // bool hit = false;
    // // Intersect with all triangles of all meshes.
    // for (const auto& mesh : m_pScene->meshes) {
    //     for (const auto& tri : mesh.triangles) {
    //         const auto v0 = mesh.vertices[tri[0]];
    //         const auto v1 = mesh.vertices[tri[1]];
    //         const auto v2 = mesh.vertices[tri[2]];
    //         if (intersectRayWithTriangle(v0.p, v1.p, v2.p, ray, hitInfo)) {
    //             hitInfo.material = mesh.material;
    //             hit = true;
    //         }
    //     }
    // }
    // // Intersect with spheres.
    // for (const auto& sphere : m_pScene->spheres)
    //     hit |= intersectRayWithShape(sphere, ray, hitInfo);
    // return hit;

    if(getRootIndex() < 0)
        return false;

    return intersectBVH(getRootIndex(), ray, hitInfo);
}

void BoundingVolumeHierarchy::loadObjectsFromScene() 
{
    for (const auto& mesh : m_pScene->meshes) {
        for (const auto &tri : mesh.triangles)
        {
            const auto v0 = mesh.vertices[tri[0]];
            const auto v1 = mesh.vertices[tri[1]];
            const auto v2 = mesh.vertices[tri[2]];
            std::array<Vertex, 3> triangle {v0, v1, v2};
            triangles.push_back(triangle);
            triangle_materials.push_back(mesh.material);
        }
    }
    // Intersect with spheres.
    for (const auto& sphere : m_pScene->spheres) {
        spheres.push_back(sphere);
    }
}

int BoundingVolumeHierarchy::getRootIndex() const {
    if(nodes.size() == 0)
        return -1;

    return 0;
}

void BoundingVolumeHierarchy::constructBVH() {
    std::queue<std::pair<int, int>> node_indices_level_queue; // queue storing a node index and it's level

    // load the objects from the scene
    loadObjectsFromScene();


    // create a vector of all indices and the is_triangle vectors
    std::vector<int> object_indices;
    std::vector<bool> is_triangle;
    for(int i = 0; i < triangles.size(); i++) {
        object_indices.push_back(i);
        is_triangle.push_back(true);
    }

    for(int i = 0; i < spheres.size(); i++) {
        object_indices.push_back(i);
        is_triangle.push_back(false);
    }

    // cannot create BVH from 0 objects
    if(object_indices.size() == 0) {
        max_level_achieved = -1;
    }

    int root_level = 0;
    max_level_achieved = 0;


    // create the root node
   createNodeAndUpdateStats(object_indices, is_triangle, root_level);



   node_indices_level_queue.push(std::make_pair(nodes.size() - 1, root_level));

   // generate the tree using BFS

   while (!node_indices_level_queue.empty())
   {
       // get the node from the queue
       auto node_index_level = node_indices_level_queue.front();
       int node_index = node_index_level.first;
       int level = node_index_level.second;
       max_level_achieved = std::max(max_level_achieved, level);
       node_indices_level_queue.pop();


       auto node_object_indices = object_indices_per_node[node_index];
       auto node_is_triangle = is_triangle_stats_per_node[node_index];

       // if the node is a leaf, add the objects as children and and stop splitting
       if (nodes.at(node_index).isLeaf)
       {
           nodes[node_index].children = node_object_indices;
           nodes[node_index].is_triangle = node_is_triangle;

           continue;
       }

       // we are now at the next level
       level++;

       // create the left and the right nodes

       // sort the objects based on the relevant axis
       sortObjects(node_object_indices, node_is_triangle, level);


       // split the objects into 2
       std::vector<int> node_object_indices_left, node_object_indices_right;
       std::vector<bool> node_is_triangle_left, node_is_triangle_right;

       // the left vectors
       for (int i = 0; i < (node_object_indices.size() + 1) / 2; i++)
       {
           node_object_indices_left.push_back(node_object_indices[i]);
           node_is_triangle_left.push_back(node_is_triangle[i]);
       }

       // the right vectors
       for (int i = (node_object_indices.size() + 1) / 2; i < node_object_indices.size(); i++)
       {
           node_object_indices_right.push_back(node_object_indices[i]);
           node_is_triangle_right.push_back(node_is_triangle[i]);
       }


       // if the size is smaller than 0, there is no point in creating a node
       if (node_object_indices_left.size() > 0)
       {
           // create the left and the right vectors and add them to the queue
           // left
           createNodeAndUpdateStats(node_object_indices_left, node_is_triangle_left, level);
           node_indices_level_queue.push(std::make_pair(nodes.size() - 1, level)); // add the left node to the queue
           // add the child to the parent
           nodes[node_index].children.push_back(nodes.size() - 1);
        }


        if(node_object_indices_right.size() > 0) {
            // right node
            createNodeAndUpdateStats(node_object_indices_right, node_is_triangle_right, level);
            node_indices_level_queue.push(std::make_pair(nodes.size() - 1, level)); // add the right node to the queue

            nodes[node_index].children.push_back(nodes.size() - 1);
        }

    }
}

void BoundingVolumeHierarchy::createNodeAndUpdateStats(std::vector<int> &object_indices, std::vector<bool> &is_triangle, int level) {

    if(!validateObjectVectorPair(object_indices, is_triangle))
        return;


    nodes.push_back(createNodeFromObjects(object_indices, is_triangle, level));


    object_indices_per_node.insert(std::make_pair(nodes.size() - 1, object_indices));


    is_triangle_stats_per_node.insert(std::make_pair(nodes.size() - 1, is_triangle));

    node_indices_per_level[level].push_back(nodes.size() - 1);

}

Node BoundingVolumeHierarchy::createNodeFromObjects(std::vector<int> &object_indices, std::vector<bool> &is_triangle, int level) const {

    if(!validateObjectVectorPair(object_indices, is_triangle))
        return Node();


    if(object_indices.size() <= 1 || level >= max_level) {
        // leaf node
        Node n;
        n.isLeaf = true;
        n.AABB = createAabbFromObjects(object_indices, is_triangle);


        return n;
    }


    // normal node
    Node n;
    n.isLeaf = false;
    n.AABB = createAabbFromObjects(object_indices, is_triangle);


    return n;

}


AxisAlignedBox BoundingVolumeHierarchy::createAabbFromObjects(std::vector<int> &object_indices, std::vector<bool> &is_triangle) const {
    if(!validateObjectVectorPair(object_indices, is_triangle))
        return AxisAlignedBox {glm::vec3(0), glm::vec3(0)};
    
    glm::vec3 p_min = glm::vec3(std::numeric_limits<float>::max());
    glm::vec3 p_max = glm::vec3(-std::numeric_limits<float>::max());


    for(int i = 0; i < object_indices.size(); i++) {


        if(is_triangle[i]) {
            // triangle case
            std::array<Vertex, 3> triangle = triangles.at(object_indices[i]);

            p_min = glm::min(glm::min(p_min, triangle[0].p), glm::min(triangle[1].p, triangle[2].p));
            p_max = glm::max(glm::max(p_max, triangle[0].p), glm::max(triangle[1].p, triangle[2].p));
        }
        else {
            // sphere case
            Sphere sph = spheres.at(object_indices[i]);
            glm::vec3 sph_min = sph.center - glm::vec3(sph.radius);
            glm::vec3 sph_max = sph.center + glm::vec3(sph.radius);

            p_min = glm::min(p_min, glm::min(sph_min, sph_max));
            p_max = glm::max(p_max, glm::max(sph_min, sph_max));

        }
    }


    return AxisAlignedBox {p_min, p_max};
}

void BoundingVolumeHierarchy::sortObjects(std::vector<int> &object_indices, std::vector<bool> &is_triangle, int level) const {
    if(!validateObjectVectorPair(object_indices, is_triangle))
        return;

    std::vector<std::pair<float, int>> attribute_index;

    for(int i = 0; i < object_indices.size(); i++) {
        attribute_index.push_back(std::make_pair(getSortingAttribute(object_indices[i], is_triangle[i], level), i));
    }

    std::sort(attribute_index.begin(), attribute_index.end());

    std::vector<int> object_indices_copy(object_indices);
    std::vector<bool> is_triangle_copy(is_triangle);

    for(int i = 0; i < object_indices_copy.size(); i++) {
        object_indices[i] = object_indices_copy[attribute_index[i].second];
        is_triangle[i] = is_triangle_copy[attribute_index[i].second];
    }
    
}

float BoundingVolumeHierarchy::getSortingAttribute(int object_index, bool is_triangle, int level) const {
    if(is_triangle)
        return getSortingAttributeTriangle(object_index, level);
    
    return getSortingAttributeSphere(object_index, level);
}

float BoundingVolumeHierarchy::getSortingAttributeTriangle(int triangle_index, int level) const {
    if (!validateTriangleIndex(triangle_index))
        return 0.0f;

    int attr = level % 3;

    std::array<Vertex, 3> triangle = triangles.at(triangle_index);

    if(attr == 0) {
        return (triangle[0].p.x + triangle[1].p.x + triangle[2].p.x) / 3;
    }
    else if(attr == 1) {
        return (triangle[0].p.y + triangle[1].p.y + triangle[2].p.y) / 3;
    }

    return (triangle[0].p.z + triangle[1].p.z + triangle[2].p.z) / 3;
}


float BoundingVolumeHierarchy::getSortingAttributeSphere(int sphere_index, int level) const {

    if(!validateSphereIndex(sphere_index))
        return 0.0f;

    int attr = level % 3;

    Sphere sphere = spheres.at(sphere_index);

    // return either the x, y or z coordinate of the center of the sphere based on the level

    if(attr == 0) {
        return sphere.center.x;
    }
    else if(attr == 1) {
        return sphere.center.y;
    }

    return sphere.center.z;
}

bool BoundingVolumeHierarchy::validateObjectVectorPair(std::vector<int> &object_indices, std::vector<bool> &is_triangle) const {
    if (object_indices.size() != is_triangle.size()) {
        std::cerr << "Invalid pair!" << std::endl;
        return false;
    }

    return true;
}

bool BoundingVolumeHierarchy::validateSphereIndex(int index) const {
    if (index < 0 || index >= spheres.size()) {
        std::cerr << "Invalid sphere index " << index << std::endl;
        return false;
    }

    return true;       
}

bool BoundingVolumeHierarchy::validateTriangleIndex(int index) const {
    if (index < 0 || index >= triangles.size()) {
        std::cerr << "Invalid triangle index " << index << std::endl;
        return false;
    }

    return true;
}

bool BoundingVolumeHierarchy::intersectNode(int node_index, Ray &ray) const{
    float t_old = ray.t;
    bool intersect = intersectRayWithShape(nodes.at(node_index).AABB, ray);
    ray.t = t_old;

    return intersect;
}

bool BoundingVolumeHierarchy::intersectObject(int object_index, bool is_triangle, Ray &ray, HitInfo &hitInfo) const {
    if(is_triangle) {
        auto triangle = triangles.at(object_index);
        auto material = triangle_materials.at(object_index);
        return intersectRayWithTriangleWithInterpolation(triangle[0], triangle[1], triangle[2], ray, hitInfo, material);
    }

    // is sphere
    auto sphere = spheres.at(object_index);
    return intersectRayWithShape(sphere, ray, hitInfo);
}

bool BoundingVolumeHierarchy::intersectBVH(int node_index, Ray &ray, HitInfo &hitInfo) const {

    // if the current node is not intersected, return false
    if(!intersectNode(node_index, ray))
        return false;

    Node n = nodes.at(node_index);

    if(n.isLeaf) {
        // do the intersection with the child shapes and return
        bool interesected = false;

        for(int i = 0; i < n.children.size(); i++) {
            int object_index = n.children.at(i);
            bool is_triangle = n.is_triangle.at(i);

            interesected |= intersectObject(object_index, is_triangle, ray, hitInfo);
        }



        return interesected;
    }

    // if it is not a leaf, do the intersection with the nodes
    bool intersected = false;
    for(int i = 0; i < n.children.size(); i++) {
        int node_index = n.children.at(i);

        intersected |= intersectBVH(node_index, ray, hitInfo);

    }

    return intersected;
}