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
    Meshes plm;

    string path;
    double tol = 10*numeric_limits<double>::epsilon();

    cout << "Insert file path (ex: DFN/FR10_data.txt): ";
    cin >> path;

    if(!ImportFractures(path, dfn, tol))
    {
        return 1;
    }
    else
    {
        cout << "Il numero di fratture del DFN e': " << dfn.fractureNumber << endl;
        for(unsigned int i=0; i<dfn.fractureNumber; i++){
            Fracture fract = dfn.fractures[i];
            cout << "L'Id della frattura e': " << fract.id << " e ha "<< fract.verticesNumber << " vertici" << endl;
            cout << "Le coordinate dei vertici sono: " << endl;
            for(unsigned int j=0; j<3; j++){
                cout << "[ " ;
                for(unsigned int k=0; k<fract.verticesNumber; k++){
                    cout << fixed << setprecision(16) << fract.verticesCoordinates(j,k) << " ";
                }
                cout <<  " ]" << endl;
            }
        }
    }

    FindTraces(dfn.fractures,tol,dfn);
    cout << "Il numero di tracce e': " << dfn.traces.size() << endl;

    PrintGlobalResults("results.txt", dfn.traces);
    PrintLocalResults("lresults.txt",dfn.fractures,dfn.traces);

    CreateMesh(dfn.fractures,tol,dfn,plm);

    for(unsigned int i=0; i<plm.meshes.size(); i++){
        cout << "Numero di celle 2D: " << plm.meshes[i].numberCell2D << endl;
    }

    PrintMeshes("meshes.txt",plm);

    return 0;
}
