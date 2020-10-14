#include "bounding_volume_hierarchy.h"
#include "draw.h"


// utility function, used for debugging
void printAABB(AxisAlignedBox aabb) {
    std::cout << "Lower " << aabb.lower.x << " " << aabb.lower.y << " " << aabb.lower.z << std::endl;
    std::cout << "Upper " << aabb.upper.x << " " << aabb.upper.y << " " << aabb.upper.z << std::endl;
    std::cout << std::endl;
}

unsigned long long Node::id = 0;
unsigned long long BvhObject::id = 0;

BoundingVolumeHierarchy::BoundingVolumeHierarchy(Scene* pScene)
    : m_pScene(pScene)
{
    // reset the node and BvhObject reference ids;

    Node::resetIdGenerator();
    BvhObject::resetIdGenerator();

    // as an example of how to iterate over all meshes in the scene, look at the intersect method below

    // load the data from the scene
    addBvhObjectsFromScene();

    // BVH Tree Creation

    // get all the bvhObjectIds
    std::vector<unsigned long long> bvhObjectIds = getAllMapKeys(bvhObjects);


    // create the root node at level 0
    _root_id = createNode(bvhObjectIds, 0);


    // create a queue for a BFS traversal of the tree
    std::queue<std::pair<unsigned long long, unsigned int>> node_level_queue;

    // insert the node _root_id and it's level inside the queue
    node_level_queue.push(std::make_pair(_root_id, 0));

    // do a BFS traversal of the tree
    while(node_level_queue.size() > 0) {
        std::pair<unsigned long long, unsigned int> node_level = node_level_queue.front();
        node_level_queue.pop();
        
        unsigned long long node_id = node_level.first;
        unsigned int level = node_level.second;
        Node& node = nodes[node_id];


        _max_level_achieved = std::max(_max_level_achieved, level);



        // if the node is a leaf, we should stop splitting
        if(node.isLeaf())
            continue;
        
        // otherwise, split the objects of the node into 2 groups
        std::vector<unsigned long long> objectsAtNode = objects_at_node[node_id];
        std::array<std::vector<unsigned long long>, 2> first_second = medianSplit(objectsAtNode, level % 3);

        // create 2 new nodes based on the split, and add them as children to the current node as well as in the queue, for further processing
        unsigned long long node_left_id = createNode(first_second[0], level + 1);
        unsigned long long node_right_id = createNode(first_second[1], level + 1);

        node.addChild(nodes[node_left_id]);
        node.addChild(nodes[node_right_id]);

        // add them to the queue
        node_level_queue.push(std::make_pair(node_left_id, level + 1));
        node_level_queue.push(std::make_pair(node_right_id, level + 1));
    }

    // for debugging only
    //printHierarchy();

    
}

unsigned long long BoundingVolumeHierarchy::createNode(std::vector<unsigned long long>& bvhObjectIds, unsigned int level) {

    bool isLeaf = (bvhObjectIds.size() <= 1 || _max_level == level) ? true : false;

    std::vector<BvhObject> objects = getMapValuesOfKeys(bvhObjects, bvhObjectIds);

    

    Node n = Node(objects, isLeaf);


    // if the node is a leaf node, add references to the objects
    if(isLeaf) {
        for(int i = 0; i < objects.size(); i++) {
            n.addChild(objects[i]);
        }
        // add the node to the leaf nodes list
        leafNodes.push_back(n.getId());
    }

    // add node to the nodes map
    nodes.insert(std::make_pair(n.getId(), n));
    // add node to nodes at level map
    nodes_at_level[level].push_back(n.getId());
    // asscociate the objects to the node
    objects_at_node.insert(std::make_pair(n.getId(), bvhObjectIds));

    return n.getId();

}

// Use this function to visualize your BVH. This can be useful for debugging. Use the functions in
// draw.h to draw the various shapes. We have extended the AABB draw functions to support wireframe
// mode, arbitrary colors and transparency.
void BoundingVolumeHierarchy::debugDraw(int level, bool showLeafNodes)
{

    // // Draw the AABB as a transparent green box.
    // //AxisAlignedBox aabb{ glm::vec3(-0.05f), glm::vec3(0.05f, 1.05f, 1.05f) };
    // //drawShape(aabb, DrawMode::Filled, glm::vec3(0.0f, 1.0f, 0.0f), 0.2f);

    // // Draw the AABB as a (white) wireframe box.
    // AxisAlignedBox aabb { glm::vec3(-0.05f), glm::vec3(0.05f, 1.05f, 1.05f) };
    // //drawAABB(aabb, DrawMode::Wireframe);
    // drawAABB(aabb, DrawMode::Filled, glm::vec3(0.05f, 1.0f, 0.05f), 0.1);

    unsigned int level_to_extract = std::min((unsigned int)level, _max_level);

    std::vector<unsigned long long> nodes_to_draw_ids = nodes_at_level[level_to_extract];
    std::vector<Node> nodes_to_draw = getMapValuesOfKeys(nodes, nodes_to_draw_ids);

    for(int i = 0; i < nodes_to_draw.size(); i++) {
        Node &n = nodes_to_draw[i];
        // draw the AABB in green
        if(!showLeafNodes)
            drawAABB(n.getBoundingBox(), DrawMode::Filled, glm::vec3(0.05f, 1.0f, 0.05f), 0.1);
    }

    // if the option is active, draw all of the leaf nodes in red.
    if(showLeafNodes) {
        // get all of the leaf nodes
        std::vector<Node> leafNodesToDraw = getMapValuesOfKeys(nodes, leafNodes);

        //std::cout << leafNodesToDraw.size() << std::endl;

        for (int i = 0; i < leafNodesToDraw.size(); i++) {
            Node &n = leafNodesToDraw[i];
            // draw the leaf node AABB in red.
            drawAABB(n.getBoundingBox(), DrawMode::Filled, glm::vec3(1.0f, 0.05f, 0.05f), 0.5);
        }   
    }
}


// loads the data from the scene as BvhObjects
void BoundingVolumeHierarchy::addBvhObjectsFromScene() {
    // create BvhObjects for all elements of the scene (spheres and triangles), and add them to the bvhObjects map

    // triangles

    for (const auto &mesh : m_pScene->meshes)
    {
        for (const auto &tri : mesh.triangles)
        {
            glm::vec3 v0 = mesh.vertices[tri[0]].p;
            glm::vec3 v1 = mesh.vertices[tri[1]].p;
            glm::vec3 v2 = mesh.vertices[tri[2]].p;

            BvhObject obj(v0, v1, v2);

            bvhObjects.insert(std::make_pair(obj.getId(), obj));
        }
    }

    // spheres

    for (int i = 0; i < m_pScene->spheres.size(); i++)
    {
        Sphere sphere = m_pScene->spheres[i];

        BvhObject obj = BvhObject(sphere);

        bvhObjects.insert(std::make_pair(obj.getId(), obj));
    }
}

    // performs a split of the provided objects, based on the median value
    // split_axis = 0 - for x_axis; 1 - for y_axis; 2 - for z_axis
std::array<std::vector<unsigned long long>, 2> BoundingVolumeHierarchy::medianSplit(std::vector<unsigned long long> bvhObjectIds, unsigned int split_axis)
{
    if(split_axis > 2)
        throw std::invalid_argument("Invalid axis");

    std::vector<std::pair<float, unsigned long long>> axis_values;

    // extract the coordinates of the correct axis
    for(auto id : bvhObjectIds) {
        BvhObject obj = bvhObjects[id];
        AxisAlignedBox aabb = obj.getAABB();

        float weight = aabb.lower.x;

        if(split_axis == 1)
            weight = aabb.lower.y;
        else
            weight = aabb.lower.z;

        axis_values.push_back(std::make_pair(weight, id));
    }

    // the half vector
    std::vector<unsigned long long> firstHalf;
    std::vector<unsigned long long> secondHalf;

    // Old method, less efficient
    // // sort the coordinate vector
    // std::sort(axis_values.begin(), axis_values.end());

    // // split in 2 halves
    
    // for(int i = 0; i < axis_values.size() / 2; i++)
    //     firstHalf.push_back(axis_values[i].second);
    
    // for(int i = axis_values.size() / 2 ; i < axis_values.size(); i++)
    //     secondHalf.push_back(axis_values[i].second);

    // New method, more efficient: Use quickselect (n_th element) to find the median point

    int medianPosition = axis_values.size() / 2;
    std::nth_element(axis_values.begin(), axis_values.begin() + medianPosition + 1, axis_values.end());
    // get the median value:
    std::pair<float, unsigned long long> medianElement = axis_values[medianPosition];

    // split the data in 2 halves, one smaller or equal to the median, the other larger

    for(int i = 0; i < axis_values.size(); i++) {
        std::pair<float, unsigned long long> element = axis_values[i];

        // if it is smaller, put it in the first half
        if(element.first <= medianElement.first) {
            firstHalf.push_back(element.second);
        }
        // else, put it in the second half
        else {
            secondHalf.push_back(element.second);
        }
    }


    std::array<std::vector<unsigned long long>, 2> arr;
    arr[0] = firstHalf;
    arr[1] = secondHalf;

    return arr;
}

int BoundingVolumeHierarchy::numLevels() const
{
    return _max_level_achieved + 1;
}

// Return true if something is hit, returns false otherwise. Only find hits if they are closer than t stored
// in the ray and if the intersection is on the correct side of the origin (the new t >= 0). Replace the code
// by a bounding volume hierarchy acceleration structure as described in the assignment. You can change any
// file you like, including bounding_volume_hierarchy.h .
bool BoundingVolumeHierarchy::intersect(Ray& ray, HitInfo& hitInfo) const
{
    bool hit = false;
    // Intersect with all triangles of all meshes.
    for (const auto& mesh : m_pScene->meshes) {
        for (const auto& tri : mesh.triangles) {
            const auto v0 = mesh.vertices[tri[0]];
            const auto v1 = mesh.vertices[tri[1]];
            const auto v2 = mesh.vertices[tri[2]];
            if (intersectRayWithTriangle(v0.p, v1.p, v2.p, ray, hitInfo)) {
                hitInfo.material = mesh.material;
                hit = true;
            }
        }
    }
    // Intersect with spheres.
    for (const auto& sphere : m_pScene->spheres)
        hit |= intersectRayWithShape(sphere, ray, hitInfo);
    return hit;
}

// prints a textual representation of the hierarchy
void BoundingVolumeHierarchy::printHierarchy() {
    std::cout << "View per level" << std::endl;
    std::cout << "----\n\n" << std::endl;

    for(int i = 0; i <= _max_level_achieved; i++) {
        std::vector<unsigned long long> ids = nodes_at_level[i];

        std::cout << "Level " << i << std:: endl;

        for(const auto& id : ids) {
            std::cout << "Node " << id << std::endl;
            std::vector<unsigned long long> objs = objects_at_node[id];
            for(const auto& obj_id : objs) {
                std::cout << obj_id << " ";
            }
            std::cout << std::endl;
        }

        std::cout << std::endl;
        std::cout << std::endl;
    }

}

// BvhObject Implementation


//ctr default
BvhObject::BvhObject() {
    

}


// ctr triangle
BvhObject::BvhObject(glm::vec3& v0, glm::vec3& v1, glm::vec3& v2) {

    _id = id++;

    // construct AABB for the triangle
    _boundingBox.lower.x = std::min(std::min(v0.x, v1.x), v2.x);
    _boundingBox.lower.y = std::min(std::min(v0.y, v1.y), v2.y);
    _boundingBox.lower.z = std::min(std::min(v0.z, v1.z), v2.z);

    _boundingBox.upper.x = std::max(std::max(v0.x, v1.x), v2.x);
    _boundingBox.upper.y = std::max(std::max(v0.y, v1.y), v2.y);
    _boundingBox.upper.z = std::max(std::max(v0.z, v1.z), v2.z);

    // set the type of the object to triangle
    _type = _TRIANGLE_TYPE;

    // store reference to the triangle
    _triangle[0] = v0;
    _triangle[1] = v1;
    _triangle[2] = v2;
}

// cre sphere
BvhObject::BvhObject(Sphere& sphere) {

    _id = id++;



    // construct the AABB for the sphere
    _boundingBox.lower.x = sphere.center.x - sphere.radius;
    _boundingBox.lower.y = sphere.center.y - sphere.radius;
    _boundingBox.lower.z = sphere.center.z - sphere.radius;

    _boundingBox.upper.x = sphere.center.x + sphere.radius;
    _boundingBox.upper.y = sphere.center.y + sphere.radius;
    _boundingBox.upper.z = sphere.center.z + sphere.radius;

    

    // set the type
    _type = _SPHERE_TYPE;

    // store reference to the sphere
    _sphere[0] = sphere;
}

std::array<glm::vec3, 3>& BvhObject::getTriangle() {
    if (_type != _TRIANGLE_TYPE)
        throw std::exception(); 

    return _triangle;
}

Sphere& BvhObject::getSphere() {
    if (_type != _SPHERE_TYPE)
        throw std::exception();

    return _sphere[0];
}

AxisAlignedBox& BvhObject::getAABB() {
    return _boundingBox;
}

std::string BvhObject::getType() {
    return _type;
}

unsigned long long BvhObject::getId() {

    return _id;
}

// after calling the this function, the ids of the further declared objects will start again at 0.
void BvhObject::resetIdGenerator() {
    BvhObject::id = 0;
}

// Node Implementation


Node::Node() {

}

Node::Node(std::vector<BvhObject>& containingBvhObjectss, bool isLeaf) {

    // set the id
    _id = id++;

    // set the isLeaf

    _isLeaf = isLeaf;

    // initialize the bounding box
    _boundingBox.lower.x = std::numeric_limits<float>::max();
    _boundingBox.lower.y = std::numeric_limits<float>::max();
    _boundingBox.lower.z = std::numeric_limits<float>::max();

    _boundingBox.upper.x = -std::numeric_limits<float>::max();
    _boundingBox.upper.y = -std::numeric_limits<float>::max();
    _boundingBox.upper.z = -std::numeric_limits<float>::max();


    

    for (auto object : containingBvhObjectss) {
        updateAABB(object.getAABB());
    }
}

bool Node::isLeaf() {
    return _isLeaf;
}

void Node::addChild(Node& other) {
    if(isLeaf()) {
        throw std::exception();
    }

    _children.push_back(other.getId());
    
}

void Node::addChild(BvhObject& object) {
    if (!isLeaf()){
        throw std::exception();
    }

    _children.push_back(object.getId());
}

void Node::updateAABB(AxisAlignedBox& other) {

    // the AABB corresponding to the node should be the smallest one that contains all other AABBs
    // hence the minimum point should be the minimum of all the children's minimum points
    // the maximum point should be the maximum of all the children's maximum points.
    _boundingBox.lower.x = std::min(other.lower.x, _boundingBox.lower.x);
    _boundingBox.lower.y = std::min(other.lower.y, _boundingBox.lower.y);
    _boundingBox.lower.z = std::min(other.lower.z, _boundingBox.lower.z);

    _boundingBox.upper.x = std::max(other.upper.x, _boundingBox.upper.x);
    _boundingBox.upper.y = std::max(other.upper.y, _boundingBox.upper.y);
    _boundingBox.upper.z = std::max(other.upper.z, _boundingBox.upper.z);
}


unsigned long long Node::getId() {
    return _id;
}

AxisAlignedBox Node::getBoundingBox() {
    return _boundingBox;
}

// after calling the this function, the ids of the further declared ndoes will start again at 0.
void Node::resetIdGenerator()
{
    Node::id = 0;
}

// Utility functions to get all keys and all values from a map
template <typename K, typename V>
std::vector<K> BoundingVolumeHierarchy::getAllMapKeys(std::map<K, V> mp){

    std::vector<K> keys;

    for(auto const& it : mp) {
        keys.push_back(it.first);
    }

    return keys;
}

template <typename K, typename V>
std::vector<V> BoundingVolumeHierarchy::getAllMapValues(std::map<K, V> mp) {

    std::vector<V> values;

    for(const auto& it : mp) {
        values.push_back(it->second);
    }

    return values;
}

template <typename K, typename V> 
std::vector<V> BoundingVolumeHierarchy::getMapValuesOfKeys(std::map<K, V> mp, std::vector<K> keys) {
    std::vector<V> values;

    for(const auto& key : keys) {
        values.push_back(mp[key]);
    }

    return values;
}