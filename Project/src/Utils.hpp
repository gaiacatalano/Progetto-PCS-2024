#pragma once
#include "DFN.hpp"
#include "PolygonalMesh.hpp"

using namespace std;
using namespace DiscreteFractureNetworkLibrary;

namespace DiscreteFractureNetworkLibrary{

bool ImportFractures(const string& filePath, DFN& dfn, const double tol);

double PointsDistance(const Vector3d p1, const Vector3d p2);

bool Parallel(const Fracture& f1, const Fracture& f2, const double tol);

bool IntersectionSphere(const Fracture& f1, const Fracture& f2, const double tol);

array<Vector3d, 2> LineIntersection(const Fracture &f1, const Fracture &f2);

pair<bool,Vector2d> InterFractureLine(const Fracture &f1, const array<Vector3d, 2> &line, const double tol);

void FindTraces(vector<Fracture> &fractures, DFN &dfn, const double tol);

void PrintGlobalResults (const string& fileName, const vector<Trace>& traces);

void PrintLocalResults (const string& fileName, vector<Fracture>& fractures, vector<Trace>& traces, const double tol);

bool CompareTraceLength(const unsigned int& id1, const unsigned int& id2, const vector<Trace>& traces, const double tol);

void SortTracesByLength(vector<unsigned int>& vecIdTraces, vector<Trace>& traces, const double tol);

}

namespace PolygonalMeshLibrary{

int PositionVert(const Vector3d& point, const array<Vector3d,2> line, const double tol);

void CutAndSave(PolygonalMesh &mesh, const unsigned int &polygonId, const array<Vector3d, 2> &line, array<unsigned int, 2> &vertIdsHelp, array<unsigned int, 2> &pointsIntersIds, const double tol);

void CutFracture(PolygonalMesh &Mesh, const DFN &dfn, unsigned int &polygonId, const double tol);

void CorrectMesh(PolygonalMesh& mesh);

void CreateMesh(const DFN &dfn, Meshes &meshesVector, const double tol);

void PrintMeshes (const string& fileName, Meshes &meshesVector);

}
