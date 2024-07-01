#ifndef __TEST_HPP
#define __TEST_HPP

#include <gtest/gtest.h>
#include "Utils.hpp"
#include "Eigen/Eigen"
#include <cmath>
#include <iostream>
#include "PolygonalMesh.hpp"

using namespace Eigen;
using namespace std;
using namespace PolygonalMeshLibrary;
using namespace DiscreteFractureNetworkLibrary;


double tol=10*numeric_limits<double>::epsilon();



TEST(IMPORTTEST, TestImportFractures) {
    string path = "DFN/FR8_test.txt";
    string patherr = "DFN/FR_errato.txt";
    DFN dfn;
    double tol=10*numeric_limits<double>::epsilon();

    bool b = ImportFractures(path, dfn, tol);
    bool c = ImportFractures(patherr, dfn, tol);
    ASSERT_TRUE(b);
    ASSERT_FALSE(c);
    EXPECT_EQ(dfn.fractureNumber, 8);
    EXPECT_EQ(dfn.fractures[0].verticesNumber,4);
    EXPECT_EQ(dfn.fractures[1].verticesNumber,4);
    EXPECT_EQ(dfn.fractures[2].verticesNumber,4);
    EXPECT_EQ(dfn.fractures[3].verticesNumber,4);
    EXPECT_EQ(dfn.fractures[4].verticesNumber,4);
    EXPECT_EQ(dfn.fractures[5].verticesNumber,3);
    EXPECT_EQ(dfn.fractures[6].verticesNumber,3);
    EXPECT_EQ(dfn.fractures[7].verticesNumber,4);


    Matrix<double, 3, 4> vm0 {
        {0, 4, 4, 0}, {-2, -2, 2, 2}, {0, 0, 0, 0}
    };

    Matrix<double, 3, 4> vm1 {
        {0, 4, 4, 0}, {-2, -2, 2, 2}, {2, 2, 2, 2}
    };

    Matrix<double, 3, 4> vm2 {
        {-8, -8, -8, -8}, {-6, -4, -4, -6}, {0, 0, 1, 1}
    };

    Matrix<double, 3, 4> vm3 {
        {4, 4, -2, -2}, {-2, 2, 2, -2}, {0, 0, 3, 3}
    };

    Matrix<double, 3, 4> vm4 {
        {0, 0, 0, 0}, {4, 6, 6, 4}, {2, 2, 3, 3}
    };

    Matrix<double, 3, 3> vm5 {
        {4, 4, 7}, {-2, 2, 0}, {2, 2, 2}
    };

    Matrix<double, 3, 3> vm6 {
        {-2, -2, 0}, {-2, 2, 0}, {3, 3, 2}
    };

    Matrix<double, 3, 4> vm7 {
        {0, 4, 4, 0}, {-2, -2, 0, 0}, {0, 0, 0, 0}
    };

    EXPECT_LT((vm0-dfn.fractures[0].verticesCoordinates).norm(), tol);
    EXPECT_LT((vm1-dfn.fractures[1].verticesCoordinates).norm(), tol);
    EXPECT_LT((vm2-dfn.fractures[2].verticesCoordinates).norm(), tol);
    EXPECT_LT((vm3-dfn.fractures[3].verticesCoordinates).norm(), tol);
    EXPECT_LT((vm4-dfn.fractures[4].verticesCoordinates).norm(), tol);
    EXPECT_LT((vm5-dfn.fractures[5].verticesCoordinates).norm(), tol);
    EXPECT_LT((vm6-dfn.fractures[6].verticesCoordinates).norm(), tol);
    EXPECT_LT((vm7-dfn.fractures[7].verticesCoordinates).norm(), tol);

    Vector3d b0(2, 0, 0);
    Vector3d b1(2, 0, 2);
    Vector3d b2(-8, -5, 0.5);
    Vector3d b3(1, 0, 1.5);
    Vector3d b4(0, 5, 2.5);
    Vector3d b5(5, 0, 2);
    Vector3d b6(-1.33333333333333333333, 0, 2.666666666666666666666);
    Vector3d b7(2, -1, 0);


    EXPECT_LT((b0-dfn.fractures[0].barycenter).norm(), tol);
    EXPECT_LT((b1-dfn.fractures[1].barycenter).norm(), tol);
    EXPECT_LT((b2-dfn.fractures[2].barycenter).norm(), tol);
    EXPECT_LT((b3-dfn.fractures[3].barycenter).norm(), tol);
    EXPECT_LT((b4-dfn.fractures[4].barycenter).norm(), tol);
    EXPECT_LT((b5-dfn.fractures[5].barycenter).norm(), tol);
    EXPECT_LT((b6-dfn.fractures[6].barycenter).norm(), tol);
    EXPECT_LT((b7-dfn.fractures[7].barycenter).norm(), tol);

    Vector3d n0(0, 0, 1);
    Vector3d n1(0, 0, 1);
    Vector3d n2(1, 0, 0);
    Vector3d n3(-1, 0, -2);
    Vector3d n4(1, 0, 0);
    Vector3d n5(0, 0, 1);
    Vector3d n6(-1, 0, -2);
    Vector3d n7(0, 0, 1);

    EXPECT_LT((dfn.fractures[0].normal.cross(n0)).norm(), tol);
    EXPECT_LT((dfn.fractures[1].normal.cross(n1)).norm(), tol);
    EXPECT_LT((dfn.fractures[2].normal.cross(n2)).norm(), tol);
    EXPECT_LT((dfn.fractures[3].normal.cross(n3)).norm(), tol);
    EXPECT_LT((dfn.fractures[4].normal.cross(n4)).norm(), tol);
    EXPECT_LT((dfn.fractures[5].normal.cross(n5)).norm(), tol);
    EXPECT_LT((dfn.fractures[6].normal.cross(n6)).norm(), tol);
    EXPECT_LT((dfn.fractures[7].normal.cross(n7)).norm(), tol);



}

// DFN dfn;
// ImportFractures("DFN/FR8_test.txt", dfn, tol);
// Fracture f0 = dfn.Fractures[0];
// Fracture f1 = dfn.Fractures[1];
// Fracture f2 = dfn.Fractures[2];
// Fracture f3 = dfn.Fractures[3];
// Fracture f4 = dfn.Fractures[4];
// Fracture f5 = dfn.Fractures[5];
// Fracture f6 = dfn.Fractures[6];
// Fracture f7 = dfn.Fractures[7];

TEST(DISTANCETEST, TestPointsDistance){

    Vector3d p1(1, 0, 0);
    Vector3d p2(7, 0, 0);
    Vector3d p3(7+ tol/10, 0, 0);
    double d12 = PointsDistance(p1, p2);
    double d23 = PointsDistance(p2, p3);
    EXPECT_EQ(d12, 6);
    EXPECT_LT(d23, tol);

}

TEST(DISTANCETEST, TestParallel) {
    DFN dfn;
    ImportFractures("DFN/FR8_test.txt", dfn, tol);
    Fracture f0 = dfn.fractures[0];
    Fracture f1 = dfn.fractures[1];
    Fracture f3 = dfn.fractures[3];
    EXPECT_TRUE(Parallel(f0, f1, tol));
    EXPECT_FALSE(Parallel(f0, f3, tol));
}

TEST(DISTANCETEST, TestIntersectionSphere) {
    DFN dfn;
    ImportFractures("DFN/FR8_test.txt", dfn, tol);
    Fracture f0 = dfn.fractures[0];
    Fracture f1 = dfn.fractures[1];
    Fracture f2 = dfn.fractures[2];
    Fracture f4 = dfn.fractures[4];
    EXPECT_TRUE(IntersectionSphere(f2, f4, tol));
    EXPECT_FALSE(IntersectionSphere(f0, f1, tol));
}

Vector3d v(0,1,0);  // vettore lungo EH
Vector3d null(0, 0, 0);

TEST(TRACESTEST, TestLineIntersection)
{
    DFN dfn;
    ImportFractures("DFN/FR8_test.txt", dfn, tol);
    Fracture f1 = dfn.fractures[1];
    Fracture f3 = dfn.fractures[3];
    Vector3d prod = LineIntersection(f1, f3)[1].cross(v);
    EXPECT_LT(prod.norm(), tol);
}

TEST(TRACESTEST, TestInterFractureLine)
{
    DFN dfn;
    ImportFractures("DFN/FR8_test.txt", dfn, tol);
    Fracture f3 = dfn.fractures[3];
    Fracture f5 = dfn.fractures[5];
    Fracture f6 = dfn.fractures[6];
    Vector3d E(0, -2, 2);
    array<Vector3d, 2> r = {E, v};  // retta EH
    pair<bool, Vector2d> a = InterFractureLine(f6, r, tol);
    pair<bool, Vector2d> b = InterFractureLine(f3, r, tol);
    pair<bool, Vector2d> c = InterFractureLine(f5, r, tol);

    EXPECT_EQ(a.first, false);  // un solo punto non costituisce intersezione
    EXPECT_EQ(b.first, true);
    EXPECT_EQ(c.first, false);

}

TEST(TRACESTEST, TestFindTraces)
{
    DFN dfn;
    ImportFractures("DFN/FR8_test.txt", dfn, tol);
    FindTraces(dfn.fractures,tol,dfn);
    Fracture f0 = dfn.fractures[0];
    Fracture f1 = dfn.fractures[1];
    Fracture f2 = dfn.fractures[2];
    Fracture f3 = dfn.fractures[3];
    Fracture f4 = dfn.fractures[4];
    Fracture f5 = dfn.fractures[5];
    Fracture f6 = dfn.fractures[6];
    Fracture f7 = dfn.fractures[7];


    // controllo che f4 non abbia tracce (quindi che non ne abbia anche con f6)
    EXPECT_TRUE(f4.passingTraces.empty());
    EXPECT_TRUE(f4.notPassingTraces.empty());

    // controllo che f6 non abbia tracce (così verifico che un punto non costituisca una traccia e che le sovrapposzizioni non siano considerate)
    EXPECT_TRUE(f6.passingTraces.empty());
    EXPECT_TRUE(f6.notPassingTraces.empty());

    // controllo che l'unica traccia non passante di f3 sia quella con f7

    EXPECT_EQ(f3.notPassingTraces[0], f7.passingTraces[0]);

    // controllo che f3 abbia due tracce passanti (una con f0 e una con f1)
    EXPECT_EQ(f3.passingTraces.size(), 2);

}

TEST(TRACESTEST, TestSort){
    vector<unsigned int> ids = {0, 1};
    Trace a;
    a.idTrace = 0;
    a.length = 2;
    Trace b;
    b.idTrace = 1;
    b.length = 1;

    vector<Trace> trs = {a, b};

    SortTracesByLength(ids, trs);
    vector<unsigned int> w = {0, 1};
    EXPECT_TRUE(w == ids);

}

#endif
