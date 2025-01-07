#include <iostream>
#include <vector>
#include <set>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string>
#include <sstream>
#include <algorithm>
#include <utility>
#include <unordered_map>
#include <stdint.h>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;

typedef struct HitInfo
{
    int pointIdx;
    double d;
    vector<double> hitPoint;
    vector<double> hitLocation;
    bool hit;
    HitInfo()
    {
        pointIdx = 0;
        hitPoint.resize(3, 0.0);
        d = 0.0;
        hit = false;
    }
    HitInfo(int _pointIdx, vector<double> _hitPoint, vector<double> _hitLocation, double _d, bool _hit)
    {
        pointIdx = _pointIdx;
        hitPoint = _hitPoint;
        hitLocation = _hitLocation;
        d = _d;
        hit = _hit;
    }
    bool operator<(const HitInfo &hitInfo) const
    {
        return (d < hitInfo.d);
    }
} HitInfo;

class pLine
{
public:
    vector<double> p0, p1, dir;
    pLine();
    pLine(vector<double> &lp0, vector<double> &p1_dir, int isp1orDir);
    ~pLine();
    void getDir();
    void getP1();
};

class pointOctNode
{
public:
    static const int MAX_OCTNODE_OBJECTS = 50;
    static const int NUM_BRANCHES_OCTNODE = 8;
    double size;
    int level;
    vector<double> position;
    vector<pointOctNode> branches;
    vector<int> data;
    vector<double> low, upp;
    pointOctNode();
    pointOctNode(int _level, vector<double> _position, double _size);
    ~pointOctNode();
    bool isLeafNode();
    void addNode(int _level, vector<double> _position, double _size);
    void getLowUppVerts();
    bool containsPoint(vector<double> point);
    bool containsSphere(const vector<double> point, double radius);
    bool boxRayIntersect(pLine &ray);
    bool sphereRayIntersect(pLine &ray);

    bool intersectsBox(const vector<double> &boxMin, const vector<double> &boxMax) const;

    size_t numPoints();
    void addPoint(int _indx);
};

class pointOctree
{
public:
    static const int MAX_OCTREE_LEVELS = 16;
    int branchOffsets[8][3];
    int max_depth;
    pointOctNode root;
    vector<vector<double>> *points;
    unordered_map<uint32_t, pointOctNode *> mortonHashMap;
    pointOctree(int max_depth);
    void create(vector<vector<double>> &_points);

    ~pointOctree();
    double getSizeRoot();
    size_t numPoints();
    vector<double> getPositionRoot();
    vector<int> findRayIntersects(vector<pLine> &rayList);
    vector<int> findRayIntersectsSorted(vector<pLine> &rayList);
    vector<pointOctNode *> getSortedNodesToCheck(pLine &ray);
    void insertPoint(pointOctNode &node, vector<double> &point, int32_t _idx);
    void insertPoints();

    void splitNodeAndReallocate(pointOctNode &node);
    void getPointsToCheck(pointOctNode &node, pLine &ray, set<int> &intTestPoints);
    void getNodesToCheck(pointOctNode &node, pLine &ray, vector<pair<pointOctNode *, double>> &nodeList);

    // Queries

    void getPointIndicesIntersectingSphere(vector<double> origin, double radius, vector<size_t> &indices);
    void getPointIndicesIntersectingBox(vector<double> origin, double size, vector<size_t> &indices);

    void getPointsIntersectingSphere(vector<double> origin, double radius, vector<vector<double>> &outPoints);
    void getPointsIntersectingBox(vector<double> origin, double radius, vector<vector<double>> &outPoints);

    void sphereCast(pLine &ray, double radius, double max_distance, HitInfo &hitInfo);
    void batchSphereCast(vector<pLine> rays, double radius, double max_distance, vector<HitInfo> &hitInfos);
    void getNodes(vector<pointOctNode *> &nodes, int node_type);
    uint64_t getMortonCodeForPoint(vector<double> &point, int32_t level);
    void createMortonHashMap(pointOctNode *node);
};

bool sortNodes(const pair<pointOctNode *, double> &i, const pair<pointOctNode *, double> &j);
double distBetweenPoints(vector<double> &p1, vector<double> &p2);
vector<double> vectAdd(vector<double> &a, vector<double> &b);
vector<double> vectAdd(vector<double> &a, vector<double> &b, double sf);
vector<double> vectSubtract(vector<double> &a, vector<double> &b);
vector<double> vectScalarMult(vector<double> &a, double b);
double dotProduct(vector<double> &v1, vector<double> &v2);
vector<double> crossProduct(vector<double> &v1, vector<double> &v2);