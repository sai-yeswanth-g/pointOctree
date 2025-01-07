#include "pointOctree.h"
#include <chrono>
#include <numeric>
#include <memory>

pLine::pLine()
{
    // Default constructor for pLine
    // Default line is unit vector along x-axis
    p0.resize(3, 0.0);
    p1.resize(3, 0.0);
    dir.resize(3, 0.0);
    p1[0] = 1.0;
    dir[0] = 1.0;
}

pLine::pLine(vector<double> &_p0, vector<double> &p1_dir, int isP1orDir)
{
    // pLine constructor with p0 and p1 or dir
    // if isP1orDir==0, then p1_dir is p1
    // if isP1orDir==1, then p1_dir is dir
    p0 = _p0;
    if (isP1orDir == 0)
    {
        p1 = p1_dir;
        getDir();
    }
    else if (isP1orDir == 1)
    {
        dir = p1_dir;
        getP1();
    }
}

pLine::~pLine() {}

void pLine::getDir()
{
    // Get unit vector defining direction of pLine
    vector<double> p0p1(3);
    double dmag = 0.0;
    for (unsigned int i = 0; i < 3; i++)
    {
        p0p1[i] = p1[i] - p0[i];
        dmag += pow(p0p1[i], 2.0);
    }
    dmag = sqrt(dmag);
    dir = p0p1;
    for (vector<double>::iterator it = dir.begin(); it != dir.end(); ++it)
        *it /= dmag;
}

void pLine::getP1()
{
    // Get a point on the pLine, p1, located 1.0 units away from the origin, p0
    vector<double> p1(3);
    for (unsigned int i = 0; i < 3; i++)
        p1[i] = p0[i] + dir[i];
}

// ------------------------------------------------------

pointOctNode::pointOctNode()
{
    // Default octNode constructor
    level = 0;
    size = 1.0;
    position.resize(3, 0.0);
    getLowUppVerts();
    data.reserve(MAX_OCTNODE_OBJECTS);
}

pointOctNode::pointOctNode(int _level, vector<double> _position, double _size)
{
    level = _level;
    position = _position;
    size = _size;
    getLowUppVerts();
    data.reserve(MAX_OCTNODE_OBJECTS);
}

// octNode destructor
pointOctNode::~pointOctNode()
{
    // cout << "Calling destructor for pointOctNode " << nid << endl;
}

bool pointOctNode::isLeafNode()
{
    // Checks if pointOctNode is a leaf node by counting the number of branches. A
    // leaf node has no branches
    return branches.size() == 0;
}

void pointOctNode::getLowUppVerts()
{
    // Get coordinates of the lower and upper vertices of the pointOctNode
    low.resize(3);
    upp.resize(3);
    double halfSize = size / 2.0;
    for (int i = 0; i < 3; i++)
    {
        low[i] = position[i] - halfSize;
        upp[i] = position[i] + halfSize;
    }
}

void pointOctNode::addPoint(int _indx) { data.push_back(_indx); }

size_t pointOctNode::numPoints() { return data.size(); }

void pointOctNode::addNode(int _level, vector<double> _position, double _size)
{
    branches.push_back(pointOctNode(_level, _position, _size));
}

bool pointOctNode::containsPoint(const vector<double> point)
{
    return (low[0] <= point[0] && point[0] <= upp[0]) &&
           (low[1] <= point[1] && point[1] <= upp[1]) &&
           (low[2] <= point[2] && point[2] <= upp[2]);
}

bool pointOctNode::containsSphere(const vector<double> point, double radius)
{
    auto pos = position;
    vector<double> clipped {
        clamp(point[0],low[0],upp[0]),
        clamp(point[1],low[1],upp[1]),
        clamp(point[2],low[2],upp[2])
    };
    vector<double> dist = vectSubtract(point,clipped);
    bool contains = dotProduct(dist,dist) <= radius * radius;

    return contains;
}

bool pointOctNode::sphereRayIntersect(pLine &ray)
{
    // Quick test for determining if a ray is *likely* to intersect a given node

    // Radius of sphere that contains node
    double radius = distBetweenPoints(low, position);

    // Project centre of sphere (node.position) onto ray
    vector<double> oc = vectSubtract(position, ray.p0);
    double s = dotProduct(oc, ray.dir);
    vector<double> projpnt = vectAdd(ray.p0, ray.dir, s);
    double dist = distBetweenPoints(projpnt, position);

    // If distance between spherical centre and projected point is
    // less than the radius of the sphere, then an intersection is
    // *possible*
    return (dist <= radius);
}

bool pointOctNode::intersectsBox(const vector<double> &boxMin, const vector<double> &boxMax) const
{

    return ((upp[0] >= boxMin[0]) && (low[0] <= boxMax[0])) &&
           ((upp[1] >= boxMin[1]) && (low[1] <= boxMax[1])) &&
           ((upp[2] >= boxMin[2]) && (low[2] <= boxMax[2]));
}

bool pointOctNode::boxRayIntersect(pLine &ray)
{
    // An accurate test for determining if a ray will intersect a given node.
    // Tests for intersections between the ray and all 6 faces of the node.

    vector<double> p;
    double D, sDen, sNum, s, tol = 1.0e-06;
    int i, j;

    for (int faceId = 1; faceId <= 6; faceId++)
    {
        // Get D (distance of plane to origin) and N (face normal) of node face
        vector<double> N(3, 0.0);
        switch (faceId)
        {
        case 1:
        {
            D = -low[0];
            N[0] = -1.0;
            break;
        } // -x face
        case 2:
        {
            D = -low[1];
            N[1] = -1.0;
            break;
        } // -y face
        case 3:
        {
            D = -low[2];
            N[2] = -1.0;
            break;
        } // -z face
        case 4:
        {
            D = upp[0];
            N[0] = 1.0;
            break;
        } // +x face
        case 5:
        {
            D = upp[1];
            N[1] = 1.0;
            break;
        } // +y face
        case 6:
        {
            D = upp[2];
            N[2] = 1.0;
        } // +z face
        }

        // Get intersection point between face plane and ray. If no intersection is
        // possible (i.e. the normal of the face is perp. to the line) then skip face
        sDen = dotProduct(ray.dir, N);
        if (fabs(sDen) > tol)
        {

            // Find intersection point p
            sNum = D - dotProduct(ray.p0, N);
            s = sNum / sDen;
            p = vectAdd(ray.p0, ray.dir, s);

            // Check if intersection point is within bounds of face. If so, then
            // return true. If not, then skip face
            if (faceId == 1 || faceId == 4)
            {
                i = 1;
                j = 2;
            } // -x,+x
            else if (faceId == 2 || faceId == 5)
            {
                i = 0;
                j = 2;
            } // -y,+y
            else if (faceId == 3 || faceId == 6)
            {
                i = 0;
                j = 1;
            } // -z,+z
            if ((p[i] >= low[i] && p[i] <= upp[i]) && (p[j] >= low[j] && p[j] <= upp[j]))
            {
                return true;
            }
        }
    }
    return false;
}

pointOctree::pointOctree(int _max_depth)
{
    if (_max_depth <= MAX_OCTREE_LEVELS)
    {
        max_depth = _max_depth;
    }
    else
    {
        max_depth = MAX_OCTREE_LEVELS;
    }
}

void pointOctree::create(vector<vector<double>> &_points)
{
    points = &_points;
    int _offsets[][3] = {{-1, -1, -1}, {+1, -1, -1}, {-1, +1, -1}, {+1, +1, -1}, {-1, -1, +1}, {+1, -1, +1}, {-1, +1, +1}, {+1, +1, +1}};

    for (int i = 0; i < 8; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            branchOffsets[i][j] = _offsets[i][j];
        }
    }

    vector<double> position = getPositionRoot();
    double size = getSizeRoot();
    root = pointOctNode(1, position, size);
    insertPoints();
}

size_t pointOctree::numPoints() { return points->size(); }

void pointOctree::insertPoint(pointOctNode &node, vector<double> &point, int32_t _idx)
{
    if (node.isLeafNode())
    {
        if (node.containsPoint(point))
        {
            if (node.numPoints() < node.MAX_OCTNODE_OBJECTS)
            {
                // Add point if the node has capacity
                node.addPoint(_idx);
            }
            else
            {
                // Add point to node and split if max capacity reached
                node.addPoint(_idx);
                if (node.level < max_depth)
                {
                    splitNodeAndReallocate(node);
                }
            }
        }
    }
    else
    {
        auto sign = vectSubtract(point, node.position);
        int branchIndex = 0;

        for (size_t i = 0; i < point.size(); ++i)
        {
            if (point[i] > node.position[i]) // Assume `center` is the node's center
            {
                branchIndex |= (1 << i); // Set bit corresponding to dimension
            }
        }

        insertPoint(node.branches[branchIndex], point, _idx);
    }
}

void pointOctree::insertPoints()
{
    auto start = std::chrono::high_resolution_clock::now();
    auto function_start = std::chrono::high_resolution_clock::now();
    auto n_points = numPoints();
    for (uint32_t i = 0; i < n_points; i++)
    {

        insertPoint(root, (*points)[i], i);
        if (i % 1000000 == 0)
        {
            std::chrono::duration<double> dt = std::chrono::high_resolution_clock::now() - start;
            cout << "creating octree : " << i * 100 / (float)numPoints() << " % " << "interval time : " << dt.count() << "s \r";
            start = std::chrono::high_resolution_clock::now();
        }
    }
    cout << endl;
    std::chrono::duration<double> dt = std::chrono::high_resolution_clock::now() - function_start;
    cout << "total time taken : " << dt.count() << "s" << endl;
}

vector<double> pointOctree::getPositionRoot()
{

    // Get low and upp
    vector<double> low, upp, position(3);
    auto points_ref = (*points);
    low = points_ref[0];
    upp = points_ref[0];
    for (size_t i = 1; i < points_ref.size(); i++)
    {
        for (int j = 0; j < 3; j++)
        {
            if (points_ref[i][j] < low[j])
            {
                low[j] = points_ref[i][j];
            }
            if (points_ref[i][j] > upp[j])
            {
                upp[j] = points_ref[i][j];
            }
        }
    }
    // Center of node is average of low and upp
    for (int i = 0; i < 3; i++)
    {
        position[i] = 0.5 * (low[i] + upp[i]);
    }
    return position;
}

double pointOctree::getSizeRoot()
{

    // Get low and upp
    vector<double> low, upp, range;
    auto points_ref = (*points);
    low = points_ref[0];
    upp = points_ref[0];
    for (size_t i = 1; i < points_ref.size(); i++)
    {
        for (int j = 0; j < 3; j++)
        {
            if (points_ref[i][j] < low[j])
            {
                low[j] = points_ref[i][j];
            }
            if (points_ref[i][j] > upp[j])
            {
                upp[j] = points_ref[i][j];
            }
        }
    }
    // Range is the size of the node in each coord direction
    range = vectSubtract(upp, low);
    double size = range[0];
    for (int i = 1; i < 3; i++)
    {
        if (range[i] > size)
        {
            size = range[i];
        }
    }
    size *= 1.05;
    // Scale up size of node by 5%
    return size;
}

void pointOctree::splitNodeAndReallocate(pointOctNode &node)
{
    // Split node into 8 branches
    vector<double> position(3);
    for (int i = 0; i < node.NUM_BRANCHES_OCTNODE; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            position[j] = node.position[j] + 0.25 * node.size * branchOffsets[i][j];
        }
        node.addNode(node.level + 1, position, 0.5 * node.size);
    }

    // Reallocate date from node to branches
    for (int i = 0; i < node.NUM_BRANCHES_OCTNODE; i++)
    {
        for (int j = 0; j < node.numPoints(); j++)
        {
            int indx = node.data[j];
            if (node.branches[i].containsPoint((*points)[indx])) //.isInNode(node.branches[i])
            {
                if (node.branches[i].numPoints() < node.MAX_OCTNODE_OBJECTS)
                {
                    node.branches[i].addPoint(indx);
                }
                else
                {
                    // cout << "splitting cloud " << endl;
                    splitNodeAndReallocate(node.branches[i]);
                }
            }
        }
    }
    node.data.resize(0);
}

pointOctree::~pointOctree()
{
    // cout << "Destroying the pointOctree" << endl;
}

void pointOctree::getPointsToCheck(pointOctNode &node, pLine &ray, set<int> &intTestPoints)
{
    // Utility function for getListPolysToCheck. Finds all OctNodes hit by a given ray
    // and returns a list of the objects contained within
    if (node.sphereRayIntersect(ray))
    {
        if (node.boxRayIntersect(ray))
        {
            if (node.isLeafNode())
            {
                for (int i = 0; i < node.numPoints(); i++)
                {
                    intTestPoints.insert(node.data[i]);
                }
            }
            else
            {
                for (int i = 0; i < node.NUM_BRANCHES_OCTNODE; i++)
                {
                    getPointsToCheck(node.branches[i], ray, intTestPoints);
                }
            }
        }
    }
}

vector<pointOctNode *> pointOctree::getSortedNodesToCheck(pLine &ray)
{
    // Finds all the nodes that intersect with given ray. Uses the nodes "position"
    // to sort the nodes by distance from the ray origin (in ascending order).
    // Nodes that are closest to the ray origin will be checked first for poly
    // intersections
    vector<pair<pointOctNode *, double>> nodeList;
    getNodesToCheck(root, ray, nodeList);
    sort(nodeList.begin(), nodeList.end(), sortNodes);
    vector<pointOctNode *> nodeListSorted;
    vector<pair<pointOctNode *, double>>::iterator it;
    for (it = nodeList.begin(); it != nodeList.end(); it++)
    {
        nodeListSorted.push_back((*it).first);
    }
    return nodeListSorted;
}

void pointOctree::getNodesToCheck(pointOctNode &node, pLine &ray, vector<pair<pointOctNode *, double>> &nodeList)
{
    // Utility function for getSortedNodesToCheck
    // Finds all the nodes that intersect with given ray. Projects the node "position" (node
    // centre) onto the ray to facilitate sorting of the nodes by distance from ray origin
    if (node.sphereRayIntersect(ray))
    {
        if (node.boxRayIntersect(ray))
        {
            if (node.isLeafNode())
            {
                // Project node "position" on to ray and find distance from ray origin
                vector<double> oc = vectSubtract(node.position, ray.p0);
                double s = dotProduct(oc, ray.dir);
                // Add node and distance to list
                nodeList.push_back(std::make_pair(&node, s));
            }
            else
            {
                for (int i = 0; i < node.NUM_BRANCHES_OCTNODE; i++)
                {
                    getNodesToCheck(node.branches[i], ray, nodeList);
                }
            }
        }
    }
}
void getNodesIntersectingBox(pointOctNode &node, vector<double> boxMin, vector<double> boxMax, vector<pointOctNode *> &nodes)
{
    // Check if the node's bounding box intersects with the input box
    bool intersects = node.intersectsBox(boxMin, boxMax);

    if (intersects)
    {
        // If the node is a leaf and has data, add it to the list
        if (node.isLeafNode())
        {
            if (!node.data.empty())
                nodes.push_back(&node);
        }
        else
        {
            // If the node is not a leaf, recursively check its branches
            for (int i = 0; i < 8; i++)
            {
                auto &branchNode = node.branches[i];
                getNodesIntersectingBox(branchNode, boxMin, boxMax, nodes);
            }
        }
    }
}

void getNodesIntersectingSphere(pointOctNode &node, vector<double> origin, double radius, vector<pointOctNode *> &nodes)
{
    bool contains = node.containsSphere(origin, radius);

    if (contains)
    {
        if (node.isLeafNode())
        {
            if (!node.data.empty())
                nodes.push_back(&node);
        }

        else
        {
            for (int i = 0; i < 8; i++)
            {
                auto &branchNode = node.branches[i];
                getNodesIntersectingSphere(branchNode, origin, radius, nodes);
            }
        }
    }
}
void pointOctree::getPointIndicesIntersectingSphere(vector<double> origin, double radius, vector<size_t> &indices)
{
    vector<pointOctNode *> nodes;

    getNodesIntersectingSphere(root, origin, radius, nodes);

    for (auto node : nodes)
    {
        for (auto _p_idx : node->data)
        {
            auto point = (*points)[_p_idx];
            auto distvec = vectSubtract(origin, point);
            if (dotProduct(distvec, distvec) <= radius * radius)
            {
                indices.push_back(_p_idx);
            }
        }
    }
}

void pointOctree::getPointIndicesIntersectingBox(vector<double> origin, double size, vector<size_t> &indices)
{
    vector<double> halfsize{size / 2.0, size / 2.0, size / 2.0};
    vector<double> min = vectSubtract(origin, halfsize);
    vector<double> max = vectAdd(origin, halfsize);
    vector<pointOctNode *> nodes;
    getNodesIntersectingBox(root, min, max, nodes);
    for (auto node : nodes)
    {
        for (auto _p_idx : node->data)
        {
            auto point = (*points)[_p_idx];
            if (node->containsPoint(point))
            {
                indices.push_back(_p_idx);
            }
        }
    }
}

void pointOctree::getPointsIntersectingSphere(vector<double> origin, double size, vector<vector<double>> &outPoints)
{
    vector<size_t> indices;
    getPointIndicesIntersectingSphere(origin, size, indices);
    outPoints.reserve(indices.size());
    for (auto idx : indices)
    {
        outPoints.push_back((*points)[idx]);
    }
}

void pointOctree::getPointsIntersectingBox(vector<double> origin, double radius, vector<vector<double>> &outPoints)
{
    vector<size_t> indices;
    getPointIndicesIntersectingBox(origin, radius, indices);
    outPoints.reserve(indices.size());
    for (auto idx : indices)
    {
        outPoints.push_back((*points)[idx]);
    }
}

void pointOctree::sphereCast(pLine &ray, double radius, double max_distance, HitInfo &hitInfo)
{
    set<int> point_indices;
    auto nodes = getSortedNodesToCheck(ray);
    int node_idx = -1;
    int point_idx = -1;
    vector<double> hitLocation;
    double d = INFINITY;
    double closest_distance = INFINITY;
    hitInfo.hit = false;
    for (size_t n_idx = 0; n_idx < nodes.size(); n_idx++)
    {
        auto node = nodes[n_idx];
        if (node->data.size() != 0)
        {
            double min_dist_r = max_distance;
            for (size_t p_idx = 0; p_idx < node->data.size(); p_idx++)
            {
                int idx = node->data[p_idx];
                auto p = (*points)[idx];
                auto OP = vectSubtract(p, ray.p0);
                double dist_r = dotProduct(OP, ray.dir);
                auto pl_vec = crossProduct(OP, ray.dir);
                double dist_p_sq = dotProduct(pl_vec, pl_vec);

                if (dist_p_sq <= radius * radius)
                {
                    double d_minus = sqrtf(radius * radius - dist_p_sq);
                    double hit_distance = dist_r - d_minus;
                    if (dist_r < 0)
                        continue;
                    if (hit_distance <= max_distance && hit_distance < closest_distance)
                    {
                        closest_distance = hit_distance;
                        point_idx = idx;
                        node_idx = n_idx;
                        d = hit_distance;
                        auto lvec = vectScalarMult(ray.dir, hit_distance);
                        hitLocation = vectAdd(ray.p0, lvec);
                    }
                }
            }
            if (node_idx != -1)
            {
                hitInfo.d = d;
                hitInfo.hit = true;
                hitInfo.hitPoint = (*points)[point_idx];
                hitInfo.pointIdx = point_idx;
                hitInfo.hitLocation = hitLocation;
                break;
            }
            else
            {
                hitInfo.hit = false;
            }
        }
    }
}

void pointOctree::batchSphereCast(vector<pLine> rays, double radius, double max_distance, vector<HitInfo> &hitInfos)
{

#pragma omp parallel for
    for (size_t i = 0; i < rays.size(); i++)
    {
        sphereCast(rays[i], radius, max_distance, hitInfos[i]);
    }
}
void pointOctree::getNodes(vector<pointOctNode *> &nodes, int node_type)
{
    /*
    node_type :
        0 - non empty leaf nodes
        1 - all leaf nodes
    */
    vector<pointOctNode *> stack;
    stack.push_back(&root);
    while (stack.size() != 0)
    {
        auto current = stack.back();
        stack.pop_back();
        if (current->isLeafNode())
        {
            if (node_type == 0)
            {
                if (current->numPoints() != 0)
                    nodes.push_back(current);
            }
            else
            {
                nodes.push_back(current);
            }
        }
        else
        {

            for (int i = 0; i < 8; i++)
            {
                stack.push_back(&current->branches[i]);
            }
        }
    }
}

uint64_t pointOctree::getMortonCodeForPoint(vector<double> &point, int32_t level)
{
    if (!root.containsPoint(point))
        return -1; // Use max uint32_t to indicate an invalid code.

    uint64_t code = 1; // Initialize Morton code as 0.

    // Normalize the point to the [0, 1] range within the root bounds.
    auto numerator = vectSubtract(point, root.low);
    double scale = 1.0;
    auto normalized_p = vectScalarMult(numerator, 1.0 / root.size);
    // normalized_p = vectScalarMult(normalized_p, scale);
    // Compute Morton code by interleaving bits.
    for (int i = 0; i < level; i++)
    {
        for (int j = 0; j < 3; j++) // Assuming a 3D point.
        {
            code <<= 1;
            // code += (int)floor(ldexp(normalized_p[j] / scale, i));
            code |= (int)floor(ldexp(normalized_p[j], i));
            // normalized_p[j] = fmod(normalized_p[j], ldexp(scale, -i));
            normalized_p[j] = fmod(normalized_p[j], ldexp(1.0, -i));
        }
    }
    return code;
}

void pointOctree::createMortonHashMap(pointOctNode *node)
{
    // TODO store median depth level

    if (node->isLeafNode())
    {
        auto code = getMortonCodeForPoint(node->position, node->level);
        mortonHashMap[code] = node;
    }
    else
    {
        for (int i = 0; i < 8; i++)
        {
            createMortonHashMap(&(node->branches[i]));
        }
    }
}

double dotProduct(vector<double> &v1, vector<double> &v2)
{
    // Calculates dot product v1.v2
    double dp = 0.0;
    for (unsigned int i = 0; i < 3; i++)
        dp += v1[i] * v2[i];
    return dp;
}

double distBetweenPoints(vector<double> &p1, vector<double> &p2)
{
    // Calculate the distance between points p1 and p2, |p1-p2|
    double sum = 0.0;
    for (unsigned int i = 0; i < 3; i++)
        sum += pow((p1[i] - p2[i]), 2.0);
    sum = sqrt(sum);
    return sum;
}

vector<double> crossProduct(vector<double> &v1, vector<double> &v2)
{
    // Calculates cross product v1xv2
    vector<double> cp(3);
    cp[0] = v1[1] * v2[2] - v1[2] * v2[1];
    cp[1] = v1[2] * v2[0] - v1[0] * v2[2];
    cp[2] = v1[0] * v2[1] - v1[1] * v2[0];
    return cp;
}

vector<double> vectAdd(vector<double> &a, vector<double> &b)
{
    // Vector addition, c=a+b
    return vectAdd(a, b, 1.0);
}

vector<double> vectAdd(vector<double> &a, vector<double> &b, double sf)
{
    // Vector addition and scaling, c=a+sf*b
    vector<double> c(a.size());
    for (int i = 0; i < a.size(); i++)
        c[i] = a[i] + sf * b[i];
    return c;
}
vector<double> vectSubtract(const vector<double> &a, const vector<double> &b)
{
    // Vector subtraction, c=a-b
    vector<double> c(a.size());
    for (int i = 0; i < a.size(); i++)
        c[i] = a[i] - b[i];
    return c;
}

vector<double> vectScalarMult(vector<double> &a, double b)
{
    vector<double> r{a[0] * b, a[1] * b, a[2] * b};
    return r;
}

bool sortNodes(const pair<pointOctNode *, double> &i, const pair<pointOctNode *, double> &j)
{
    // Function used to sort a vector of pointOctNode,double pairs by the value of
    // the double. The double will typically represent distance from the ray
    // origin in a ray-node intersection test
    return i.second < j.second;
}

// ------------------------------------------------------
