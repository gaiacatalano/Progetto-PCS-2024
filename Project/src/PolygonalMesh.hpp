#ifndef PolygonalMesh_H
#define PolygonalMesh_H

#include <iostream>
#include <vector>
#include <Eigen/Eigen>

using namespace std;
using namespace Eigen;

namespace PolygonalMeshLibrary{

struct PolygonalMesh{

    unsigned int numVertices;
    vector<unsigned int> idVertices;
    vector<Vector3d> coordVertices; //ci metto dentro un vector per farci operazioni matematiche

    unsigned int numEdges;
    vector<unsigned int> idEdges;
    vector<array<unsigned int, 2>> extremitiesEdges; //sono id, non serve fare operazioni matematiche

    unsigned int numPolygons;
    list<vector<unsigned int>> verticesPolygons; //non so a priori quanti sono e non mi serve accedere per indice (non hanno id)
    list<vector<unsigned int>> edgesPolygons;

};

struct plm{
    vector<PolygonalMesh> meshes = {};
};

}
#endif
