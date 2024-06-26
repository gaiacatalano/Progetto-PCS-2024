#include <iostream>
#include "./src/Utils.hpp"
#include "./src/DFN.hpp"
#include <iomanip>

using namespace std;
using namespace Eigen;
using namespace DiscreteFractureNetworkLibrary;
using namespace PolygonalMeshLibrary;

int main()
{
    DFN dfn;
    plm plm;
    string path = "DFN/FR10_data.txt";
    double tol=10*numeric_limits<double>::epsilon();

    //cout << "Insert file path: ";
    //cin >> path;

    if(!ImportFractures(path, dfn, tol))
    {
        return 1;
    }
    else
    {
        cout << "Il numero di fratture del DFN e': " << dfn.fractureNumber << endl;
        for(unsigned int i=0; i<dfn.fractureNumber; i++){
            Fracture fract = dfn.Fractures[i];
            cout << "L'Id della frattura e': " << fract.id << " e ha "<< fract.verticesNumber << " vertici" << endl;
            cout << "Le coordinate dei vertici sono: " << endl;
            for(unsigned int j=0; j<3; j++){
                cout << "[ " ;
                for(unsigned int k=0; k<fract.verticesNumber; k++){
                    cout << fixed << setprecision(16) << fract.verticesCoordinates(j,k) << " ";
                }
                cout <<  " ]" << endl;
            }
            // cout << "Il baricentro della frattura e':" << fract.barycenter << endl;
            // cout << "La normale della frattura e':" << fract.normal << endl;
            // cout << "Il raggio della frattura e':" << fract.radius << endl;
            // cout << "Il piano della frattura e':" << fract.plane << endl;
        }
    }

    FindTraces(dfn.Fractures,tol,dfn);
    cout << "Il numero di tracce e': " << dfn.Traces.size() << endl;

    printGlobalResults("results.txt", dfn.Traces);
    printLocalResults("lresults.txt",dfn.Fractures,dfn.Traces);

    CreateMesh(dfn.Fractures,tol,dfn,plm);
    for(int i=0; i<plm.meshes.size(); i++){
        cout << "Numero di celle 2D: " << plm.meshes[i].numPolygons << endl;
    }

    //tryOutput("try.txt",plm);

    return 0;
}
