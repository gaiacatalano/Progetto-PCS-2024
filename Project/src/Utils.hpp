#pragma once
#include <iostream>
#include "DFN.hpp"
#include "PolygonalMesh.hpp"

using namespace std;
using namespace DiscreteFractureNetworkLibrary;

namespace DiscreteFractureNetworkLibrary{

// lettura da file, creazione delle strutture e salvataggio dati
bool ImportFractures(const string& filePath, DFN& dfn, double tol);

// distanza tra due punti
double PointsDistance(const Vector3d p1, const Vector3d p2);

// metodo che retituisce true se 2 poligoni sono paralleli
bool Parallel(Fracture& f1, Fracture& f2, double tol);

// metodo che restituisce true se due sfere si intersecano
bool IntersectionSphere(Fracture& f1, Fracture& f2, double tol);

array<Vector3d, 2> LineIntersection(Fracture &f1, Fracture &f2);

pair<bool,Vector2d> InterFractureLine(Fracture &f1, array<Vector3d, 2> &r, double tol);

void FindTraces(vector<Fracture> &fractures, double tol, DFN &dfn);

void PrintGlobalResults (const string& fileName, vector<Trace>& traces);

void PrintLocalResults (const string& fileName, vector<Fracture>& fractures, const vector<Trace>& traces);

bool CompareTraceLength(const unsigned int& id1, const unsigned int& id2, const std::vector<Trace>& traces);

void SortTracesByLength(vector<unsigned int>& vecIdTraces, const vector<Trace>& traces);

}

namespace PolygonalMeshLibrary{

int PositionVert(const Vector3d& point, array<Vector3d,2> line, double tol);

void CutAndSave(PolygonalMesh &mesh, unsigned int &polygonId, array<Vector3d, 2> &line, array<unsigned int, 2>& vertIdsHelp, array<unsigned int, 2>& pointsIntersIds, double tol);

void CutFracture(PolygonalMesh &Mesh, DFN &dfn, unsigned int &polygonId, vector<unsigned int>& traces, double tol);

void Turn(PolygonalMesh &mesh, unsigned int &polygonId, unsigned int &edgeId);

void CorrectMesh(PolygonalMesh& mesh);

void CreateMesh(vector<Fracture> &fractures, double tol, DFN &dfn, Meshes &meshesV);

void PrintMeshes (const string& fileName, Meshes &meshesV);

}
