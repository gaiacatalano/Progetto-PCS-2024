#ifndef PolygonalMesh_H
#define PolygonalMesh_H

#include <iostream>
#include <vector>
#include <Eigen/Eigen>

using namespace std;
using namespace Eigen;

namespace PolygonalMeshLibrary{


struct PolygonalMesh{

    unsigned int numberCell0D = 0;                          // numero di celle 0D
    vector<unsigned int> idVertices = {};                   // vettore degli id delle celle 0D
    vector<Vector3d> coordVertices = {};                    // vettore di vettori di coordinate delle celle 0D

    unsigned int numberCell1D = 0;                          // numero di celle 1D
    vector<unsigned int> idEdges = {};                      // vettore degli id delle celle 1D
    vector<Vector2i> extremitiesEdges = {};                 // vettore contentne i vettori di interi con gli id dei vertici del lato
    vector<bool> active_edge = {};                          // vettore di booleani, true se il lato è attivo
    vector<Vector2i> nearPolygons = {};                     // vettore contenenti i vettori di interi con gli id dei poligoni adiacenti
    vector<Vector2i> newedge = {};                          // vettore contenenti i vettori di interi con gli id dei nuovi lati


    unsigned int numberCell2D = 0;                          // numero di celle 2D
    vector<unsigned int> idPolygon = {};                    // vettore degli id delle celle 2D
    vector<vector<unsigned int>> verticesPolygons = {};     // vettore contenenti i vettori di interi con gli id dei vertici del poligono
    vector<vector<unsigned int>> edgesPolygons = {};        // vettore contenenti i vettori di interi con gli id dei lati del poligono
    vector<bool> active_polygon = {};                       // vettore di booleani, true se il poligono è attivo

};

struct Meshes {
    vector<PolygonalMesh> meshes = {};                      // vettore contenente tutte le polygonal mesh
};

}
#endif
