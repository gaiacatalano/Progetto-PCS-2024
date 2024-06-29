#ifndef DFN_HPP
#define DFN_HPP

#endif // DFN_HPP

#pragma once
#include <vector>
#include "Eigen/Eigen"

using namespace std;
using namespace Eigen;

namespace DiscreteFractureNetworkLibrary{

struct Fracture
{
    unsigned int id = 0;                            // id frattura
    unsigned int verticesNumber = 0;                // numero di vertici
    MatrixXd verticesCoordinates = {};              // coordinate dei vertici
    Vector3d barycenter = {};                       // baricentro (coordinate)
    double radius = 0;                              // raggio = distanza baricentro - vertice più lontano
    Vector3d normal = {};                           // normale al piano contenente la frattura
    Vector4d plane = {};                            // piano contenente la frattura
    vector<unsigned int> passingTraces = {};        // lista contenente id delle tracce passanti
    vector<unsigned int> notPassingTraces = {};     // lista contenente id delle tracce non passanti
};

struct Trace
{
    unsigned int idTrace;                             // id traccia
    array <Vector3d, 2> extremitiesCoordinates = {};  // array delle coordinate dei punti estremi della traccia
    array <unsigned int, 2> fracturesIds = {};        // array degli id delle due fratture che formano la traccia
    double length = 0;                                // lunghezza della traccia (distanza estremità)
    array<bool,2> tips = {};                          // false se passante, true se non passante (prima e seconda frattura)
    array <Vector3d, 2> lineTrace = {};               // retta contenente la traccia (punto e direzione)
};

struct DFN
{
    unsigned int fractureNumber = 0;      // numero di fratture
    vector<Fracture> fractures = {};      // vettore di fratture
    vector<Trace> traces = {};            // vettore di tracce
};

}
