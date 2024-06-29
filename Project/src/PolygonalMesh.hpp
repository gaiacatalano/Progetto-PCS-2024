#ifndef PolygonalMesh_H
#define PolygonalMesh_H

#include <iostream>
#include <vector>
#include <Eigen/Eigen>

using namespace std;
using namespace Eigen;

namespace PolygonalMeshLibrary{


struct PolygonalMesh{

    unsigned int numberCell0D = 0;
    vector<unsigned int> idVertices = {};
    vector<Vector3d> coordVertices = {};

    unsigned int numberCell1D = 0;
    vector<unsigned int> idEdges = {};
    vector<Vector2i> extremitiesEdges = {};
    vector<bool> active_edge = {};
    vector<Vector2i> nearPolygons = {};
    vector<Vector2i> newedge = {};


    unsigned int numberCell2D = 0;
    vector<unsigned int> idPolygon = {};
    vector<vector<unsigned int>> verticesPolygons = {};
    vector<vector<unsigned int>> edgesPolygons = {};
    vector<bool> active_polygon = {};

};

struct plm{
    vector<PolygonalMesh> meshes = {};
};

}
#endif
