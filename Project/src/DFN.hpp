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
    int id = 0;                                     // id frattura
    unsigned int verticesNumber = 0;                // numero di vertici
    MatrixXd verticesCoordinates = {};              // coordinate dei vertici
    Vector3d barycenter = {};                       // baricentro (coordinate)
    double radius = 0;                              // raggio = distanza baricentro - vertice più lontano
    Vector3d normal = {};                           // normale al piano contenente la frattura
    Vector4d plane = {};                            // piano contenente la frattura
    vector<unsigned int> passingTraces = {};        // lista contenente id delle tracce passanti
    vector<unsigned int> notPassingTraces = {};     // lista contenente id delle tracce non passanti




public:
    void setId(int newId);
    void setVerticesNumber(unsigned int newVerticesNumber);
    void setVerticesCoordinates(const MatrixXd &newVerticesCoordinates);
    void setBarycenter(const Vector3d &newBarycenter);
    void setRadius(double newRadius);
    void setNormal(const Vector3d &newNormal);
    void setPlane(const Vector4d &newPlane);
    void setPassingTraces(const vector<unsigned int> &newPassingTraces);
    void setNotPassingTraces(const vector<unsigned int> &newNotPassingTraces);
};

struct Trace
{
    unsigned int idTrace;                             // id traccia
    array <Vector3d, 2> extremitiesCoordinates = {};  // array delle coordinate dei punti estremi della traccia
    array <int, 2> fracturesIds = {};                 // array degli id delle due fratture che formano la traccia
    double length = 0;                                // lunghezza della traccia (distanza estremità)
    array<bool,2> Tips = {};                          // false se passante, true se non passante (prima e seconda frattura)
    array <Vector3d, 2> lineTrace = {};               // retta contenente la traccia (punto e direzione)



// Costruisco una serie di metodi:

public:                                                 // La parola chiave public è un modificatore di accesso che rende i membri di una classe o struttura accessibili da qualsiasi parte del programma
    unsigned int getIdTrace() const;
    void setIdTrace(unsigned int newIdTrace);
    array<Vector3d, 2> getExtremitiesCoordinates() const;
    void setExtremitiesCoordinates(const array<Vector3d, 2> &newExtremitiesCoordinates);
    array<int, 2> getFracturesIds() const;
    void setFracturesIds(const array<int, 2> &newFracturesIds);
    double getLength() const;
    void setLength(double newLength);
    array<bool, 2> getTips() const;
    void setTips(const array<bool, 2> &newTips);
    array<Vector3d, 2> getLineTrace() const;
    void setLineTrace(const array<Vector3d, 2> &newLineTrace);
};

struct DFN
{
    unsigned int fractureNumber = 0;      // numero di fratture
    vector<Fracture> Fractures = {};      // vettore di fratture
    vector<Trace> Traces = {};            // vettore di tracce
    MatrixXi Intersections = {};          // matrice di intersezione tra le fratture
};


/* Costruisco i metodi nel seguente modo:
 * inline: suggerisce al compilatore di sostituire la chiamata a questa funzione con il corpo della funzione stessa.
 *         questo può ridurre il sovraccarico della chiamata alla funzione e migliorare le prestazioni.
 * viene assegnato al parametro verticesNumber di Trace il numero newVerticesNumber, per esempio
*/

inline void Fracture::setVerticesNumber(unsigned int newVerticesNumber)
{
    verticesNumber = newVerticesNumber;
}

inline void Fracture::setVerticesCoordinates(const MatrixXd &newVerticesCoordinates)
{
    verticesCoordinates = newVerticesCoordinates;
}

inline void Fracture::setBarycenter(const Vector3d &newBarycenter)
{
    barycenter = newBarycenter;
}

inline void Fracture::setRadius(double newRadius)
{
    radius = newRadius;
}

inline void Fracture::setNormal(const Vector3d &newNormal)
{
    normal = newNormal;
}

inline void Fracture::setPlane(const Vector4d &newPlane)
{
    plane = newPlane;
}

inline void Fracture::setPassingTraces(const vector<unsigned int> &newPassingTraces)
{
    passingTraces = newPassingTraces;
}

inline void Fracture::setNotPassingTraces(const vector<unsigned int> &newNotPassingTraces)
{
    notPassingTraces = newNotPassingTraces;
}

inline void Fracture::setId(int newId)
{
    id = newId;
}

inline array<Vector3d, 2> Trace::getExtremitiesCoordinates() const   // estrarre le estremità di una traccia
{
    return extremitiesCoordinates;
}

inline void Trace::setExtremitiesCoordinates(const array<Vector3d, 2> &newExtremitiesCoordinates)
{
    extremitiesCoordinates = newExtremitiesCoordinates;
}

inline array<int, 2> Trace::getFracturesIds() const
{
    return fracturesIds;
}

inline void Trace::setFracturesIds(const array<int, 2> &newFracturesIds)
{
    fracturesIds = newFracturesIds;
}

inline double Trace::getLength() const
{
    return length;
}

inline void Trace::setLength(double newLength)
{
    length = newLength;
}

inline array<bool, 2> Trace::getTips() const
{
    return Tips;
}

inline void Trace::setTips(const array<bool, 2> &newTips)
{
    Tips = newTips;
}

inline array<Vector3d, 2> Trace::getLineTrace() const
{
    return lineTrace;
}

inline void Trace::setLineTrace(const array<Vector3d, 2> &newLineTrace)
{
    lineTrace = newLineTrace;
}

inline unsigned int Trace::getIdTrace() const
{
    return idTrace;
}

inline void Trace::setIdTrace(unsigned int newIdTrace)
{
    idTrace = newIdTrace;
}



}
