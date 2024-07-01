#include "Utils.hpp"
#include "DFN.hpp"
#include <iostream>
#include <fstream>
#include <numeric>
#include <sstream>
#include <algorithm>


using namespace std;
using namespace Eigen;

namespace DiscreteFractureNetworkLibrary
{

// funzione che legge il file e salva le informazioni iniziali sulle fratture
bool ImportFractures(const string& filePath, DFN& dfn, double tol) {
    ifstream file(filePath);

    if (file.fail()) {
        cerr << "Error in opening the file" << endl;
        return false;
    }

    string line;
    while (getline(file, line)) {
        if (line[0] != '#') {
            break;
        }
    }

    if (line.empty() || line[0] == '#') {
        cerr << "Invalid file format" << endl;
        return false;
    }

    istringstream convertFrac(line);
    convertFrac >> dfn.fractureNumber;

<<<<<<< Updated upstream
    // ridefinisco la dimensione del vettore di fratture in dfn
=======
>>>>>>> Stashed changes
    dfn.fractures.resize(dfn.fractureNumber);

    for (unsigned int i = 0; i < dfn.fractureNumber; i++) {
        Fracture fra;
        getline(file, line);
        getline(file, line, ';');

        istringstream convertId(line);
        // salvo l'id nella frattura
        unsigned int idFra;
        convertId >> idFra;
        fra.id = idFra;

        getline(file, line);
        istringstream convertVert(line);
        unsigned int numVert;
        convertVert >> numVert;
        fra.verticesNumber = numVert;

        MatrixXd VerMatrix(3, numVert);

        getline(file, line);

        for (unsigned int j = 0; j < 3; j++) {
            getline(file, line);
            replace(line.begin(), line.end(), ';', ' ');
            istringstream convertLine(line);

            for (unsigned int k = 0; k < numVert; k++) {
                convertLine >> VerMatrix(j, k);
            }
        }

        fra.verticesCoordinates = VerMatrix;

        Vector3d barycenter = VerMatrix.rowwise().mean();
        fra.barycenter = barycenter;

        double radius = 0;
        for(unsigned int i = 0; i<fra.verticesNumber; i++)
        {
            double d = (barycenter - fra.verticesCoordinates.col(i)).norm();
            if(d - radius > tol){
                radius = d;
            }
        }
        fra.radius = radius;

        Vector3d v1 = fra.verticesCoordinates.col(0);
        Vector3d v2 = fra.verticesCoordinates.col(1);
        Vector3d v3 = fra.verticesCoordinates.col(2);
        Vector3d d1 = v1-v2;
        Vector3d d2 = v3-v2;
        Vector3d normal = d1.cross(d2);
        fra.normal = normal;

        double d = (fra.normal[0]*fra.barycenter[0] + fra.normal[1]*fra.barycenter[1] + fra.normal[2]*fra.barycenter[2]);
        Vector4d plane = {fra.normal[0], fra.normal[1], fra.normal[2], d};
        fra.plane = plane;

<<<<<<< Updated upstream
        // salvo la frattura nel vettore delle fratture
=======
>>>>>>> Stashed changes
        dfn.fractures[i] = fra;
    }

    file.close();
    return true;
}

// funzione che calcola la distanza tra due punti in norma euclidea
double PointsDistance(const Vector3d p1, const Vector3d p2)
{
    double distance = sqrt((p1-p2).transpose() * (p1-p2));
    return distance;
}

// funzione che verifica se due fratture si trovano su piani paralleli
bool Parallel(Fracture &f1, Fracture &f2, double tol)
{
    Vector3d n1 = f1.normal;
    Vector3d n2 = f2.normal;
    double tol2 = max(10*numeric_limits<double>::epsilon(), tol*tol);
    if((n1.cross(n2)).norm() < tol2){
        return true;
    }
    return false;
}

// funzione che calcola la più piccola sfera contiene la frattura
bool IntersectionSphere(Fracture &f1,  Fracture &f2, double tol)
{
    double r1 = f1.radius;
    double r2 = f2.radius;
    Vector3d b1 = f1.barycenter;
    Vector3d b2 = f2.barycenter;
    double distanceBarycenters = PointsDistance(b1,b2);
    if(distanceBarycenters - (r1 + r2) > tol){
        return true;
    }
    return false;
}

// funzione che trova la retta d'intersezione tra due piani
array<Vector3d, 2> LineIntersection(Fracture &f1, Fracture &f2) {

    // Vettore direzione della retta (prodotto vettoriale dei normali dei due piani)
    Vector3d n1 = f1.normal;
    Vector3d n2 = f2.normal;
    Vector4d p1 = f1.plane;
    Vector4d p2 = f2.plane;
    Vector3d direzione = n1.cross(n2);


    // Troviamo un punto sulla retta di intersezione. Metto a sistema i due piani e l'equazione della retta eguagliata ad un valore (0 per semplicità)
    Matrix3d A;
    A.row(0) = n1;
    A.row(1) = n2;
    A.row(2) = direzione;

    Vector3d b = {p1[3], p2[3], 0};

    // risolvo il sistema lineare e trovo il punto che uso per trovare la retta in forma parametrica
    Vector3d punto = A.fullPivLu().solve(b);

    array<Vector3d, 2> retta = {punto, direzione};

    return retta;
}

// funzione che trova le intersezioni di una frattura con una retta
pair<bool,Vector2d> InterFractureLine(Fracture &f1, array<Vector3d, 2> &r, double tol){
    double tol2 = max(10*numeric_limits<double>::epsilon(), tol*tol);
    Vector3d punto = r[0];
    Vector3d direzioneRetta = r[1];
    double den = 1/(direzioneRetta.dot(direzioneRetta));
    MatrixXd vertices = f1.verticesCoordinates;

    unsigned int count = 0;
    Vector2d beta;  // creo il vettore delle ascisse curvilinee
    bool con = false;


    for (unsigned int i = 0; i < f1.verticesCoordinates.cols(); i++) {
        // Calcolo degli indici dei vertici del segmento di linea
        unsigned int j = (i + 1) % f1.verticesCoordinates.cols();
        Vector3d ver1 = f1.verticesCoordinates.col(i);
        Vector3d ver2 = f1.verticesCoordinates.col(j);
        Vector3d direzioneLato = ver2 - ver1;

        Vector3d prodotto = direzioneLato.cross(direzioneRetta);

        if (prodotto.norm() > tol2){   // se non sono paralleli
            double alpha = ((punto-ver1).cross(direzioneRetta)).dot(prodotto)/(prodotto.dot(prodotto));

            if (alpha >= -tol && alpha < 1-tol){
                double betaTemp = ((ver1-punto).cross(direzioneLato)).dot(-prodotto)/(prodotto.dot(prodotto));
                beta[count] = betaTemp;  // aggiungo betaTemporaneo al vettore dei beta
                count++;

                if(count==2){
                    sort(beta.begin(), beta.end());   // se il poligono interseca la retta due volte, (può farlo 0 o 1 o 2), ordino le acsisse, in ordine crescente
                    con = true;
                    break;
                }
            }

        }


        else{
            double alpha = ((ver1-punto).dot(direzioneRetta))*den;
            Vector3d proiezione = punto + alpha*direzioneRetta;
            if ((proiezione-ver1).norm()<tol){
                double betaTemp = ((ver1-punto).dot(direzioneRetta))*den;

                beta[count] = betaTemp;
                count++;
                con = true;

                if(count==2){
                    sort(beta.begin(), beta.end());
                    con = true;
                    break;
                }
            }
        }

    }
    return {con, beta};
    // quindi per ogni frattura ottengo due info:
    // con: true se c'è intersezione con la retta
    // beta che contiene le ascisse dell'intersezione
}

// funzione che trova le tracce e ne salva i dati
void FindTraces(vector<Fracture> &fractures, double tol, DFN &dfn)
{
    unsigned int count = 0;
    double tol2 = max(10*numeric_limits<double>::epsilon(), tol*tol);

    for (unsigned int i=0; i<fractures.size(); i++){
        for(unsigned int j=i+1; j<fractures.size(); j++){
            Fracture &f1 = fractures[i];
            Fracture &f2 = fractures[j];
            array<Vector3d, 2> estremiTraccia;

            // controlli ed esclusioni
            if (Parallel(f1,f2,tol)){
                continue;
            }
            if (IntersectionSphere(f1,f2,tol)){
                continue;
            }

            // funzione che trova retta di intersezione dei piani contenenti le fratture f1 e f2
            array<Vector3d, 2> rettaf1xf2 = LineIntersection(f1,f2);

            if (rettaf1xf2[1].norm() < tol2){
                continue;
            }

            // funzione che trova i punti di intersezione delle due fratture con la retta
            auto [con1, betaf1] = InterFractureLine(f1, rettaf1xf2, tol);
            auto [con2, betaf2] = InterFractureLine(f2, rettaf1xf2, tol);

            if (!con1 || !con2) {  // se almeno una frattura non interseca la retta, allora le due fratture non si intersecano
                continue;
            }

            // punti in cui risp. la frattura 1 e la 2 intersecano la retta
            array<Vector3d, 2> puntiInterf1 = {rettaf1xf2[0] + (betaf1[0]*rettaf1xf2[1]), rettaf1xf2[0] + (betaf1[1]*rettaf1xf2[1])};
            array<Vector3d, 2> puntiInterf2 = {rettaf1xf2[0] + (betaf2[0]*rettaf1xf2[1]), rettaf1xf2[0] + (betaf2[1]*rettaf1xf2[1])};


            if (betaf1[1]<betaf2[0]+tol || betaf2[1]<betaf1[0]+tol){ // caso in cui sono disgiunte lungo la retta
                continue;
            }
            else { // c'è intersezione
                if (betaf1[0]<betaf2[0]+tol){
                    estremiTraccia[0] = rettaf1xf2[0] + (betaf2[0]*rettaf1xf2[1]);
                }
                else{
                    estremiTraccia[0] = rettaf1xf2[0] + (betaf1[0]*rettaf1xf2[1]);
                }

                if (betaf1[1]<betaf2[1]+tol){
                    estremiTraccia[1] = rettaf1xf2[0] + (betaf1[1]*rettaf1xf2[1]);
                }
                else{
                    estremiTraccia[1] = rettaf1xf2[0] + (betaf2[1]*rettaf1xf2[1]);
                }
            }

            array<bool,2> tips = {false, false};
            if (estremiTraccia != puntiInterf1){
                tips[0] = true;
            }
            if (estremiTraccia != puntiInterf2){
                tips[1] = true;
            }

            // calcolo la lunghezza ed escludo il caso in cui l'intersezione consista di un unico punto
            double length = PointsDistance(estremiTraccia[0],estremiTraccia[1]);
            if (length > tol){
                // creo la traccia
                Trace tr;
                tr.idTrace = count;
                tr.extremitiesCoordinates = estremiTraccia;
                tr.fracturesIds = {f1.id,f2.id};
                tr.length = length;
                tr.tips = tips;
                tr.lineTrace = rettaf1xf2;

                if(tips[0]){
                    f1.notPassingTraces.push_back(tr.idTrace);
                }
                else{
                    f1.passingTraces.push_back(tr.idTrace);

                }
                if(tips[1]){
                    f2.notPassingTraces.push_back(tr.idTrace);

                }
                else{
                    f2.passingTraces.push_back(tr.idTrace);
                }

                dfn.traces.push_back(tr);
                count++;
            }

        }
    }
}


// funzione che stampa le informazioni sulle tracce
void PrintGlobalResults (const string& fileName, vector<Trace>& traces){
    ofstream ofstr(fileName); // se il file non esiste, lo crea
    ofstr << "# Number of Traces" << endl;
    ofstr << traces.size() << endl;
    ofstr << "# TraceId; FractureId1; FractureId2; X1; Y1; Z1; X2; Y2; Z2" << endl;
    for (Trace& tr:traces){
        ofstr << tr.idTrace << "; " << tr.fracturesIds[0] << "; " << tr.fracturesIds[1] << "; " << (tr.extremitiesCoordinates[0])[0] << "; "
              << (tr.extremitiesCoordinates[0])[1] << "; " << (tr.extremitiesCoordinates[0])[2] << "; " << (tr.extremitiesCoordinates[1])[0] << "; "
              << (tr.extremitiesCoordinates[1])[1] << "; " << (tr.extremitiesCoordinates[1])[2] << endl;
    }
    ofstr.close();
}

// funzione che stampa le informazioni sulle fratture e sulle tracce corrispondenti
void PrintLocalResults (const string& fileName, vector<Fracture>& fractures, const vector<Trace>& traces){

    ofstream ofstr(fileName);
    bool firstTime;
    for (Fracture& fr:fractures){
        firstTime=true;
        ofstr << "# FractureId; NumTraces" << endl;
        ofstr << fr.id << "; " << (fr.passingTraces.size())+(fr.notPassingTraces.size()) << endl;

        SortTracesByLength(fr.passingTraces, traces);
        SortTracesByLength(fr.notPassingTraces, traces);

        for (unsigned int idTrace:fr.passingTraces){
            if (firstTime){
                ofstr << "# TraceId; Tips; Length" << endl;
                firstTime=false;
            }
            ofstr << idTrace << "; " << "false; " << traces[idTrace].length << endl;

        }
        for (unsigned int idTrace:fr.notPassingTraces){
            if (firstTime){
                ofstr << "# TraceId; Tips; Length" << endl;
                firstTime=false;
            }
            ofstr << idTrace << "; " << "true; " << traces[idTrace].length << endl;
        }
    }
    ofstr.close();
}


// funzione che confronta le tracce per lunghezza (decrescente)
bool CompareTraceLength(const unsigned int& id1, const unsigned int& id2, const vector<Trace>& traces) {
    return traces[id1].length > traces[id2].length;
}

// funzione che ordina le tracce
void SortTracesByLength(vector<unsigned int>& vecIdTraces, const vector<Trace>& traces) {
    sort(vecIdTraces.begin(), vecIdTraces.end(), [&traces](const unsigned int& id1, const unsigned int& id2) {
        return CompareTraceLength(id1, id2, traces);
    });
}

}


namespace PolygonalMeshLibrary
{

// funzione che determina la posizione di un punto rispetto ad una retta del piano,
int PositionVert(const Vector3d& point, array<Vector3d,2> retta, double tol) {
    double tol2 = max(10*numeric_limits<double>::epsilon(), tol*tol);
    Vector3d relPoint = point - retta[0];
    Vector3d crossProduct = relPoint.cross(retta[1]);
    double dotProduct = crossProduct.dot(Vector3d(0, 0, 1));
    if (dotProduct > tol2) {
        return 1;
    }
    if (dotProduct < tol2){
        return -1;
    }
    return 0;
}

// funzione che esegue il taglio del poligono salvando le informazioni successive al taglio
void CutAndSave(PolygonalMesh &mesh, unsigned int &polygonId, array<Vector3d, 2> &r, array<unsigned int, 2>& vertIdsHelp, array<unsigned int, 2>& interIds, double tol){

    double tol2 = max(10*numeric_limits<double>::epsilon(), tol*tol);
    Vector3d punto = r[0];
    Vector3d direzioneRetta = r[1];
    double den = 1/(direzioneRetta.dot(direzioneRetta));

    unsigned int count = 0;
    Vector2d beta;
    Vector2i edgeOffIds;  // id dei lati da spegnere
    array<Vector2i,2> edgeSpVert = {Vector2i(-1, -1), Vector2i(-1,-1)}; // id dei vertici dei lati da spegnere

    vector<unsigned int> vertIds = mesh.verticesPolygons[polygonId];

    for (unsigned int i = 0; i < vertIds.size(); i++) {

        // calcolo degli indici dei vertici del segmento di linea
        unsigned int j = (i + 1) % vertIds.size();
        Vector3d ver1 = mesh.coordVertices[vertIds[i]];
        Vector3d ver2 = mesh.coordVertices[vertIds[j]];
        Vector3d direzioneLato = ver2 - ver1;

        Vector3d prodotto = direzioneLato.cross(direzioneRetta);

        if (prodotto.norm() > tol2){
            double alpha = ((punto-ver1).cross(direzioneRetta)).dot(prodotto)/(prodotto.dot(prodotto));

            if (alpha >= -tol && alpha < 1-tol){
                double betaTemp = ((ver1-punto).cross(direzioneLato)).dot(-prodotto)/(prodotto.dot(prodotto));
                beta[count] = betaTemp;

                Vector2i latoSpVert = {vertIds[i],vertIds[j]};
                vertIdsHelp[count] = vertIds[i];
                edgeSpVert[count] = latoSpVert;
                unsigned int latoSpId = distance(mesh.extremitiesEdges.begin(),find(mesh.extremitiesEdges.begin(), mesh.extremitiesEdges.end(), latoSpVert));
                edgeOffIds[count] = latoSpId;

                count++;

                if(count==2){
                    break;
                }
            }
        }
        else{
            double alpha = ((ver1-punto).dot(direzioneRetta))*den;
            Vector3d proiezione = punto + alpha*direzioneRetta;
            if ((proiezione-ver1).norm()<tol){
                double betaTemp = ((ver1-punto).dot(direzioneRetta))*den;
                beta[count] = betaTemp;

                Vector2i latoSpVert = {vertIds[i],vertIds[j]};
                vertIdsHelp[count] = vertIds[i];
                edgeSpVert[count] = latoSpVert;
                unsigned int latoSpId = distance(mesh.extremitiesEdges.begin(),find(mesh.extremitiesEdges.begin(), mesh.extremitiesEdges.end(), latoSpVert));
                edgeOffIds[count] = latoSpId;

                count++;

                if(count==2){
                    break;
                }
            }
        }
    }

    // trovo e salvo le intersezioni
    array<Vector3d, 2> puntiInterFraTra = {r[0] + (beta[0]*r[1]), r[0] + (beta[1]*r[1])};

    if (edgeSpVert[0][0]!=-1 && edgeSpVert[0][1]!=-1 && edgeSpVert[1][0]!=-1 && edgeSpVert[1][1]!=-1){
        // se le intersezioni non coincidono con nessun vertice
        if (PointsDistance(puntiInterFraTra[0], mesh.coordVertices[edgeSpVert[0][0]])>tol &&
            PointsDistance(puntiInterFraTra[0], mesh.coordVertices[edgeSpVert[0][1]])>tol
            && PointsDistance(puntiInterFraTra[1], mesh.coordVertices[edgeSpVert[1][0]])>tol
            && PointsDistance(puntiInterFraTra[1], mesh.coordVertices[edgeSpVert[1][1]])>tol){

            const unsigned int numVertices = mesh.idVertices.size();

            mesh.idVertices.resize(numVertices+2);
            mesh.coordVertices.resize(numVertices+2);

            for (unsigned int in = 0; in < 2; in++){

                interIds[in] = numVertices+in;
                mesh.idVertices[interIds[in]] = interIds[in];
                mesh.coordVertices[interIds[in]] = puntiInterFraTra[in];
            }

            // nuovi lati
            array<unsigned int, 4> edgeNewIds;

            const unsigned int numEdges = mesh.idEdges.size();

            mesh.idEdges.resize(numEdges+4);
            mesh.extremitiesEdges.resize(numEdges+4);
            mesh.active_edge.resize(numEdges+4);
            mesh.nearPolygons.resize(numEdges+4);
            mesh.newedge.resize(numEdges+4);

            for (unsigned int ne = 0; ne < 4; ne++){

                edgeNewIds[ne] = numEdges+ne;
                mesh.idEdges[edgeNewIds[ne]] = edgeNewIds[ne];
                mesh.active_edge[edgeNewIds[ne]] = true;
                mesh.nearPolygons[edgeNewIds[ne]] = {-1,-1};
                mesh.newedge[edgeNewIds[ne]] = {-1,-1};
            }

            mesh.extremitiesEdges[numEdges] = {edgeSpVert[0][0],interIds[0]};
            mesh.extremitiesEdges[numEdges+1] = {interIds[0],edgeSpVert[0][1]};
            mesh.extremitiesEdges[numEdges+2] = {edgeSpVert[1][0],interIds[1]};
            mesh.extremitiesEdges[numEdges+3] = {interIds[1],edgeSpVert[1][1]};

            // spengo i lati e aggiungo i pezzi che li formano
            mesh.active_edge[edgeOffIds[0]] = false;
            if (mesh.nearPolygons[edgeOffIds[0]][0]!=-1 && mesh.active_polygon[mesh.nearPolygons[edgeOffIds[0]][0]]==true)
            {
                auto it_edge = find(mesh.edgesPolygons[mesh.nearPolygons[edgeOffIds[0]][0]].begin(), mesh.edgesPolygons[mesh.nearPolygons[edgeOffIds[0]][0]].end(), edgeOffIds[0]);
                const unsigned int local_position = it_edge - mesh.edgesPolygons[mesh.nearPolygons[edgeOffIds[0]][0]].begin();

                mesh.edgesPolygons[mesh.nearPolygons[edgeOffIds[0]][0]][local_position] = edgeNewIds[0];
                mesh.edgesPolygons[mesh.nearPolygons[edgeOffIds[0]][0]].insert(it_edge + 1, edgeNewIds[1]);

                mesh.nearPolygons[edgeNewIds[0]][0] = mesh.nearPolygons[edgeOffIds[0]][0];
                mesh.nearPolygons[edgeNewIds[1]][0] = mesh.nearPolygons[edgeOffIds[0]][0];

                auto it_vertex = mesh.verticesPolygons[mesh.nearPolygons[edgeOffIds[0]][0]].begin() + local_position + 1;
                mesh.verticesPolygons[mesh.nearPolygons[edgeOffIds[0]][0]].insert(it_vertex, interIds[0]);
            }
            if (mesh.nearPolygons[edgeOffIds[0]][1]!=-1 && mesh.active_polygon[mesh.nearPolygons[edgeOffIds[0]][1]]==true)
            {
                auto it_edge = find(mesh.edgesPolygons[mesh.nearPolygons[edgeOffIds[0]][1]].begin(), mesh.edgesPolygons[mesh.nearPolygons[edgeOffIds[0]][1]].end(), edgeOffIds[0]);
                const unsigned int local_position = it_edge - mesh.edgesPolygons[mesh.nearPolygons[edgeOffIds[0]][1]].begin();

                mesh.edgesPolygons[mesh.nearPolygons[edgeOffIds[0]][1]][local_position] = edgeNewIds[1];
                mesh.edgesPolygons[mesh.nearPolygons[edgeOffIds[0]][1]].insert(it_edge + 1, edgeNewIds[0]);

                mesh.nearPolygons[edgeNewIds[0]][0] = mesh.nearPolygons[edgeOffIds[0]][1];
                mesh.nearPolygons[edgeNewIds[1]][0] = mesh.nearPolygons[edgeOffIds[0]][1];

                auto it_vertex = mesh.verticesPolygons[mesh.nearPolygons[edgeOffIds[0]][1]].begin() + local_position + 1;
                mesh.verticesPolygons[mesh.nearPolygons[edgeOffIds[0]][1]].insert(it_vertex, interIds[0]);
            }
            mesh.newedge[edgeOffIds[0]] = {edgeNewIds[0], edgeNewIds[1]};

            mesh.active_edge[edgeOffIds[1]] = false;
            if (mesh.nearPolygons[edgeOffIds[1]][0]!=-1 && mesh.active_polygon[mesh.nearPolygons[edgeOffIds[1]][0]]==true)
            {
                auto it_edge = find(mesh.edgesPolygons[mesh.nearPolygons[edgeOffIds[1]][0]].begin(), mesh.edgesPolygons[mesh.nearPolygons[edgeOffIds[1]][0]].end(), edgeOffIds[0]);
                const unsigned int local_position = it_edge - mesh.edgesPolygons[mesh.nearPolygons[edgeOffIds[1]][0]].begin();

                mesh.edgesPolygons[mesh.nearPolygons[edgeOffIds[1]][0]][local_position] = edgeNewIds[2];
                mesh.edgesPolygons[mesh.nearPolygons[edgeOffIds[1]][0]].insert(it_edge + 1, edgeNewIds[3]);

                mesh.nearPolygons[edgeNewIds[2]][0] = mesh.nearPolygons[edgeOffIds[1]][0];
                mesh.nearPolygons[edgeNewIds[3]][0] = mesh.nearPolygons[edgeOffIds[1]][0];

                auto it_vertex = mesh.verticesPolygons[mesh.nearPolygons[edgeOffIds[1]][0]].begin() + local_position + 1;
                mesh.verticesPolygons[mesh.nearPolygons[edgeOffIds[1]][0]].insert(it_vertex, interIds[1]);
            }
            if (mesh.nearPolygons[edgeOffIds[1]][1]!=-1 && mesh.active_polygon[mesh.nearPolygons[edgeOffIds[1]][1]]==true)
            {
                auto it_edge = find(mesh.edgesPolygons[mesh.nearPolygons[edgeOffIds[1]][1]].begin(), mesh.edgesPolygons[mesh.nearPolygons[edgeOffIds[1]][1]].end(), edgeOffIds[1]);
                const unsigned int local_position = it_edge - mesh.edgesPolygons[mesh.nearPolygons[edgeOffIds[1]][1]].begin();

                mesh.edgesPolygons[mesh.nearPolygons[edgeOffIds[1]][1]][local_position] = edgeNewIds[3];
                mesh.edgesPolygons[mesh.nearPolygons[edgeOffIds[1]][1]].insert(it_edge + 1, edgeNewIds[2]);

                mesh.nearPolygons[edgeNewIds[2]][0] = mesh.nearPolygons[edgeOffIds[1]][1];
                mesh.nearPolygons[edgeNewIds[3]][0] = mesh.nearPolygons[edgeOffIds[1]][1];

                auto it_vertex = mesh.verticesPolygons[mesh.nearPolygons[edgeOffIds[1]][1]].begin() + local_position + 1;
                mesh.verticesPolygons[mesh.nearPolygons[edgeOffIds[1]][1]].insert(it_vertex, interIds[1]);
            }

            mesh.newedge[edgeOffIds[1]] = {edgeNewIds[2], edgeNewIds[3]};
        }

        // se coincidono con il primo vertice e non con il secondo
        else if ((PointsDistance(puntiInterFraTra[0], mesh.coordVertices[edgeSpVert[0][0]])<tol || PointsDistance(puntiInterFraTra[0], mesh.coordVertices[edgeSpVert[0][1]])<tol)
                 && (PointsDistance(puntiInterFraTra[1], mesh.coordVertices[edgeSpVert[1][0]])>tol && PointsDistance(puntiInterFraTra[1], mesh.coordVertices[edgeSpVert[1][1]])>tol)){

            if (PointsDistance(puntiInterFraTra[0], mesh.coordVertices[edgeSpVert[0][0]])<tol){
                interIds[0] = edgeSpVert[0][0];
            }
            else{
                interIds[0] = edgeSpVert[0][1];
            }

            const unsigned int numVertices = mesh.idVertices.size();

            mesh.idVertices.resize(numVertices+1);
            mesh.coordVertices.resize(numVertices+1);

            interIds[1] = numVertices;
            mesh.idVertices[numVertices] = numVertices;
            mesh.coordVertices[numVertices] = puntiInterFraTra[1];

            // nuovi lati
            array<unsigned int, 2> edgeNewIds;

            const unsigned int numEdges = mesh.idEdges.size();

            mesh.idEdges.resize(numEdges+2);
            mesh.extremitiesEdges.resize(numEdges+2);
            mesh.active_edge.resize(numEdges+2);
            mesh.nearPolygons.resize(numEdges+2);
            mesh.newedge.resize(numEdges+2);

            for (unsigned int ne = 0; ne < 2; ne++){

                edgeNewIds[ne] = numEdges+ne;
                mesh.idEdges[edgeNewIds[ne]] = edgeNewIds[ne];
                mesh.active_edge[edgeNewIds[ne]] = true;
                mesh.nearPolygons[edgeNewIds[ne]] = {-1,-1};
                mesh.newedge[edgeNewIds[ne]] = {-1,-1};
            }

            mesh.extremitiesEdges[numEdges] = {edgeSpVert[1][0],interIds[1]};
            mesh.extremitiesEdges[numEdges+1] = {interIds[1],edgeSpVert[1][1]};

            // spengo i lati e aggiungo i pezzi che li formano

            mesh.active_edge[edgeOffIds[1]] = false;
            if (mesh.nearPolygons[edgeOffIds[1]][0]!=-1 && mesh.active_polygon[mesh.nearPolygons[edgeOffIds[1]][0]]==true){

                auto it_edge = find(mesh.edgesPolygons[mesh.nearPolygons[edgeOffIds[1]][0]].begin(), mesh.edgesPolygons[mesh.nearPolygons[edgeOffIds[1]][0]].end(), edgeOffIds[0]);
                const unsigned int local_position = it_edge - mesh.edgesPolygons[mesh.nearPolygons[edgeOffIds[1]][0]].begin();

                mesh.edgesPolygons[mesh.nearPolygons[edgeOffIds[1]][0]][local_position] = edgeNewIds[0];
                mesh.edgesPolygons[mesh.nearPolygons[edgeOffIds[1]][0]].insert(it_edge + 1, edgeNewIds[1]);

                mesh.nearPolygons[edgeNewIds[0]][0] = mesh.nearPolygons[edgeOffIds[1]][0];
                mesh.nearPolygons[edgeNewIds[1]][0] = mesh.nearPolygons[edgeOffIds[1]][0];

                auto it_vertex = mesh.verticesPolygons[mesh.nearPolygons[edgeOffIds[1]][0]].begin() + local_position + 1;
                mesh.verticesPolygons[mesh.nearPolygons[edgeOffIds[1]][0]].insert(it_vertex, interIds[1]);
            }
            if (mesh.nearPolygons[edgeOffIds[1]][1]!=-1 && mesh.active_polygon[mesh.nearPolygons[edgeOffIds[1]][1]]==true){

                auto it_edge = find(mesh.edgesPolygons[mesh.nearPolygons[edgeOffIds[1]][1]].begin(), mesh.edgesPolygons[mesh.nearPolygons[edgeOffIds[1]][1]].end(), edgeOffIds[1]);
                const unsigned int local_position = it_edge - mesh.edgesPolygons[mesh.nearPolygons[edgeOffIds[1]][1]].begin();

                mesh.edgesPolygons[mesh.nearPolygons[edgeOffIds[1]][1]][local_position] = edgeNewIds[1];
                mesh.edgesPolygons[mesh.nearPolygons[edgeOffIds[1]][1]].insert(it_edge + 1, edgeNewIds[0]);

                mesh.nearPolygons[edgeNewIds[0]][0] = mesh.nearPolygons[edgeOffIds[1]][1];
                mesh.nearPolygons[edgeNewIds[1]][0] = mesh.nearPolygons[edgeOffIds[1]][1];

                auto it_vertex = mesh.verticesPolygons[mesh.nearPolygons[edgeOffIds[1]][1]].begin() + local_position + 1;
                mesh.verticesPolygons[mesh.nearPolygons[edgeOffIds[1]][1]].insert(it_vertex, interIds[1]);
            }
            mesh.newedge[edgeOffIds[1]] = {edgeNewIds[0],edgeNewIds[1]};
        }

        // se coincide il secondo e non il primo
        else if ((PointsDistance(puntiInterFraTra[0], mesh.coordVertices[edgeSpVert[0][0]])>tol && PointsDistance(puntiInterFraTra[0], mesh.coordVertices[edgeSpVert[0][1]])>tol)
                 && (PointsDistance(puntiInterFraTra[1], mesh.coordVertices[edgeSpVert[1][0]])<tol || PointsDistance(puntiInterFraTra[1], mesh.coordVertices[edgeSpVert[1][1]])<tol)){

            if (PointsDistance(puntiInterFraTra[1], mesh.coordVertices[edgeSpVert[1][0]])<tol){
                interIds[1] = edgeSpVert[1][0];
            }
            else{
                interIds[1] = edgeSpVert[1][1];
            }

            const unsigned int numVertices = mesh.idVertices.size();

            mesh.idVertices.resize(numVertices+1);
            mesh.coordVertices.resize(numVertices+1);

            interIds[0] = numVertices;
            mesh.idVertices[numVertices] = numVertices;
            mesh.coordVertices[numVertices] = puntiInterFraTra[0];

            // nuovi lati
            array<unsigned int, 2> edgeNewIds;

            const unsigned int numEdges = mesh.idEdges.size();

            mesh.idEdges.resize(numEdges+2);
            mesh.extremitiesEdges.resize(numEdges+2);
            mesh.active_edge.resize(numEdges+2);
            mesh.nearPolygons.resize(numEdges+2);
            mesh.newedge.resize(numEdges+2);

            for (unsigned int ne = 0; ne < 2; ne++){

                edgeNewIds[ne] = numEdges+ne;
                mesh.idEdges[edgeNewIds[ne]] = edgeNewIds[ne];
                mesh.active_edge[edgeNewIds[ne]] = true;
                mesh.nearPolygons[edgeNewIds[ne]] = {-1,-1};
                mesh.newedge[edgeNewIds[ne]] = {-1,-1};
            }

            mesh.extremitiesEdges[numEdges] = {edgeSpVert[0][0],interIds[0]};
            mesh.extremitiesEdges[numEdges+1] = {interIds[0],edgeSpVert[0][1]};

            // spengo i lati e aggiungo i pezzi che li formano

            mesh.active_edge[edgeOffIds[0]] = false;
            if (mesh.nearPolygons[edgeOffIds[0]][0]!=-1 && mesh.active_polygon[mesh.nearPolygons[edgeOffIds[0]][0]]==true){

                auto it_edge = find(mesh.edgesPolygons[mesh.nearPolygons[edgeOffIds[0]][0]].begin(), mesh.edgesPolygons[mesh.nearPolygons[edgeOffIds[0]][0]].end(), edgeOffIds[0]);
                const unsigned int local_position = it_edge - mesh.edgesPolygons[mesh.nearPolygons[edgeOffIds[0]][0]].begin();

                mesh.edgesPolygons[mesh.nearPolygons[edgeOffIds[0]][0]][local_position] = edgeNewIds[0];
                mesh.edgesPolygons[mesh.nearPolygons[edgeOffIds[0]][0]].insert(it_edge + 1, edgeNewIds[1]);

                mesh.nearPolygons[edgeNewIds[0]][0] = mesh.nearPolygons[edgeOffIds[0]][0];
                mesh.nearPolygons[edgeNewIds[1]][0] = mesh.nearPolygons[edgeOffIds[0]][0];

                auto it_vertex = mesh.verticesPolygons[mesh.nearPolygons[edgeOffIds[0]][0]].begin() + local_position + 1;
                mesh.verticesPolygons[mesh.nearPolygons[edgeOffIds[0]][0]].insert(it_vertex, interIds[0]);
            }
            if (mesh.nearPolygons[edgeOffIds[0]][1]!=-1 && mesh.active_polygon[mesh.nearPolygons[edgeOffIds[0]][1]]==true){

                auto it_edge = find(mesh.edgesPolygons[mesh.nearPolygons[edgeOffIds[0]][1]].begin(), mesh.edgesPolygons[mesh.nearPolygons[edgeOffIds[0]][1]].end(), edgeOffIds[0]);
                const unsigned int local_position = it_edge - mesh.edgesPolygons[mesh.nearPolygons[edgeOffIds[0]][1]].begin();

                mesh.edgesPolygons[mesh.nearPolygons[edgeOffIds[0]][1]][local_position] = edgeNewIds[1];
                mesh.edgesPolygons[mesh.nearPolygons[edgeOffIds[0]][1]].insert(it_edge + 1, edgeNewIds[0]);

                mesh.nearPolygons[edgeNewIds[0]][0] = mesh.nearPolygons[edgeOffIds[0]][1];
                mesh.nearPolygons[edgeNewIds[1]][0] = mesh.nearPolygons[edgeOffIds[0]][1];

                auto it_vertex = mesh.verticesPolygons[mesh.nearPolygons[edgeOffIds[0]][1]].begin() + local_position + 1;
                mesh.verticesPolygons[mesh.nearPolygons[edgeOffIds[0]][1]].insert(it_vertex, interIds[0]);
            }
            mesh.newedge[edgeOffIds[0]] = {edgeNewIds[0],edgeNewIds[1]};
        }

        // se entrambe le inter coincidono coi vertici del lato che "tagliano"
        else {

            if (PointsDistance(puntiInterFraTra[0], mesh.coordVertices[edgeSpVert[0][0]])<tol){
                interIds[0] = edgeSpVert[0][0];
            }
            else{
                interIds[0] = edgeSpVert[0][1];
            }

            if (PointsDistance(puntiInterFraTra[1], mesh.coordVertices[edgeSpVert[1][0]])<tol){
                interIds[1] = edgeSpVert[1][0];
            }
            else{
                interIds[1] = edgeSpVert[1][1];
            }

        }



    }

    //salvo il lato-traccia
    const unsigned int numEdges = mesh.idEdges.size();

    mesh.idEdges.resize(numEdges+1);
    mesh.extremitiesEdges.resize(numEdges+1);
    mesh.active_edge.resize(numEdges+1);
    mesh.nearPolygons.resize(numEdges+1);
    mesh.newedge.resize(numEdges+1);

    mesh.idEdges[numEdges] = numEdges;
    mesh.extremitiesEdges[numEdges] = {interIds[0],interIds[1]};
    mesh.active_edge[numEdges] = true;
    mesh.nearPolygons[numEdges] = {-1,-1};
    mesh.newedge[numEdges] = {-1,-1};


}

// funzione che verifica che ci siano tracce che tagliano la frattura e crea le nuove celle 2D
void CutFracture(PolygonalMesh &mesh, DFN &dfn, unsigned int &polygonId, vector<unsigned int>& traces, double tol){

    if (!traces.empty()){

        // 1) Utilizzare la funzione InterFractureLine per trovare i punti di intersezione della
        //    traccia (passante o non passannte con l'insieme di vertici in considerazione)
        // 2) Dividere i vertici a seconda che si trovino a destra o a sinistra della traccia in due sottoinsiemi
        // 3) Dividere le tracce rimanenti
        // 4) Richiamare la funzione stessa con le nuove sottofratture

        mesh.active_polygon[polygonId] = false;

        unsigned int traId = traces.front();
        Trace tra = dfn.traces[traId];

        array<Vector3d, 2> rettaTraccia = tra.lineTrace;

        array<unsigned int, 2> interIds;
        array<unsigned int, 2> vertIdsHelp;
        CutAndSave(mesh, polygonId, rettaTraccia, vertIdsHelp, interIds, tol);

        traces.erase(traces.begin());

        // divido i vertici nelle due sotto figure
        vector<unsigned int> vertIdsSub1;
        vector<unsigned int> vertIdsSub2;

        bool where = true;
        vector<unsigned int> vertIds = mesh.verticesPolygons[polygonId];

        for (unsigned int i = 0; i < vertIds.size(); i++) {

            if (where){
                if (find(vertIdsSub1.begin(), vertIdsSub1.end(), vertIds[i]) == vertIdsSub1.end()){

                    vertIdsSub1.push_back(vertIds[i]);

                    if (vertIds[i] - vertIdsHelp[0] < tol){
                        where = false;
                        if (PointsDistance(mesh.coordVertices[interIds[0]], mesh.coordVertices[vertIds[i]])>tol){
                            vertIdsSub1.push_back(interIds[0]);
                        }
                        vertIdsSub1.push_back(interIds[1]);
                    }
                }
            }
            else if (!where){
                if (find(vertIdsSub2.begin(), vertIdsSub2.end(), vertIds[i]) == vertIdsSub2.end()){

                    vertIdsSub2.push_back(vertIds[i]);

                    if (vertIds[i] - vertIdsHelp[1] < tol){
                        where = true;
                        if (PointsDistance(mesh.coordVertices[interIds[1]], mesh.coordVertices[vertIds[i]])>tol){
                            vertIdsSub2.push_back(interIds[1]);
                        }
                        vertIdsSub2.push_back(interIds[0]);
                    }
                }
            }
        }

        // salvo le due nuove sottofratture

        const unsigned int numPolygons = mesh.idPolygon.size();

        mesh.idPolygon.resize(numPolygons+2);
        mesh.verticesPolygons.resize(numPolygons+2);
        mesh.edgesPolygons.resize(numPolygons+2);
        mesh.active_polygon.resize(numPolygons+2);

        // salvo C1

        vector<unsigned int> edgeIdsSub1;
        edgeIdsSub1.resize(vertIdsSub1.size());

        for (unsigned int e=0; e<vertIdsSub1.size(); e++){
            unsigned int j = (e + 1) % vertIdsSub1.size();

            Vector2i lato = {vertIdsSub1[e],vertIdsSub1[j]};
            Vector2i lato2 = {vertIdsSub1[j],vertIdsSub1[e]};

            if (find(mesh.extremitiesEdges.begin(), mesh.extremitiesEdges.end(), lato) != mesh.extremitiesEdges.end()){
                edgeIdsSub1[e] = distance(mesh.extremitiesEdges.begin(),find(mesh.extremitiesEdges.begin(), mesh.extremitiesEdges.end(), lato));

            }
            else if (find(mesh.extremitiesEdges.begin(), mesh.extremitiesEdges.end(), lato2) != mesh.extremitiesEdges.end()){
                edgeIdsSub1[e] = distance(mesh.extremitiesEdges.begin(),find(mesh.extremitiesEdges.begin(), mesh.extremitiesEdges.end(), lato2));
            }


            if (mesh.nearPolygons[edgeIdsSub1[e]][0]!=-1 && mesh.nearPolygons[edgeIdsSub1[e]][1]!=-1){
                if (mesh.nearPolygons[edgeIdsSub1[e]][0] == polygonId){
                    mesh.nearPolygons[edgeIdsSub1[e]][0] = numPolygons;
                }
                else if (mesh.nearPolygons[edgeIdsSub1[e]][1] == polygonId){
                    mesh.nearPolygons[edgeIdsSub1[e]][1] = numPolygons;
                }
            }
            else if (mesh.nearPolygons[edgeIdsSub1[e]][0]!=-1 && mesh.nearPolygons[edgeIdsSub1[e]][1]==-1){
                mesh.nearPolygons[edgeIdsSub1[e]][1] = numPolygons;
            }
        }


        mesh.idPolygon[numPolygons] = numPolygons;
        mesh.verticesPolygons[numPolygons] = vertIdsSub1;
        mesh.edgesPolygons[numPolygons] = edgeIdsSub1;
        mesh.active_polygon[numPolygons] = true;

        // salvo C2

        vector<unsigned int> edgeIdsSub2;
        edgeIdsSub2.resize(vertIdsSub2.size());


        for (unsigned int e=0; e<vertIdsSub2.size(); e++){
            unsigned int j = (e + 1) % vertIdsSub2.size();

            Vector2i lato = {vertIdsSub2[e],vertIdsSub2[j]};
            Vector2i lato2 = {vertIdsSub2[j],vertIdsSub2[e]};

            if (find(mesh.extremitiesEdges.begin(), mesh.extremitiesEdges.end(), lato) != mesh.extremitiesEdges.end()){
                edgeIdsSub2[e] = distance(mesh.extremitiesEdges.begin(),find(mesh.extremitiesEdges.begin(), mesh.extremitiesEdges.end(), lato));

            }
            else if (find(mesh.extremitiesEdges.begin(), mesh.extremitiesEdges.end(), lato2) != mesh.extremitiesEdges.end()){
                edgeIdsSub2[e] = distance(mesh.extremitiesEdges.begin(),find(mesh.extremitiesEdges.begin(), mesh.extremitiesEdges.end(), lato2));
            }


            if (mesh.nearPolygons[edgeIdsSub2[e]][0]!=-1 && mesh.nearPolygons[edgeIdsSub2[e]][1]!=-1){
                if (mesh.nearPolygons[edgeIdsSub2[e]][0] == polygonId){
                    mesh.nearPolygons[edgeIdsSub2[e]][0] = numPolygons+1;
                }
                else if (mesh.nearPolygons[edgeIdsSub2[e]][1] == polygonId){
                    mesh.nearPolygons[edgeIdsSub2[e]][1] = numPolygons;
                }
            }
            else if (mesh.nearPolygons[edgeIdsSub2[e]][0]!=-1 && mesh.nearPolygons[edgeIdsSub2[e]][1]==-1){
                mesh.nearPolygons[edgeIdsSub2[e]][1] = numPolygons+1;
            }
        }

        mesh.idPolygon[numPolygons+1] = numPolygons+1;
        mesh.verticesPolygons[numPolygons+1] = vertIdsSub2;
        mesh.edgesPolygons[numPolygons+1] = edgeIdsSub2;
        mesh.active_polygon[numPolygons+1] = true;

        unsigned int polygonId1 = mesh.idPolygon[numPolygons];
        unsigned int polygonId2 = mesh.idPolygon[numPolygons+1];
        mesh.nearPolygons[mesh.idEdges[mesh.idEdges.size()-1]] = {polygonId1,polygonId2};

        // divido le tracce
        vector<unsigned int> tracesSub1;
        vector<unsigned int> tracesSub2;

        for (unsigned int& traId : traces){

            Trace tra = dfn.traces[traId];

            array<Vector3d, 2> estr = tra.extremitiesCoordinates;

            int w1 = PositionVert(estr[0],rettaTraccia, tol);
            int w2 = PositionVert(estr[1],rettaTraccia, tol);

            if (w1>tol && w2>tol){
                tracesSub1.push_back(traId);
            }
            else if (w1<tol && w2<tol){
                tracesSub2.push_back(traId);
            }
            else {
                tracesSub1.push_back(traId);
                tracesSub2.push_back(traId);
            }
        }

        // richiamo la funzione per le sottofratture fino alla fine delle tracce
        CutFracture(mesh, dfn, polygonId1, tracesSub1, tol);
        CutFracture(mesh, dfn, polygonId2, tracesSub2, tol);
    }

}

// funzione che si occupa di inserire i dati corretti nella mesh
void CorrectMesh(PolygonalMesh& mesh){

    mesh.numberCell0D = mesh.idVertices.size();

    for (unsigned int i=0; i<mesh.idEdges.size(); i++){
        if (mesh.active_edge[i]==true){
            mesh.numberCell1D++;
        }
    }

    for (unsigned int i=0; i<mesh.idPolygon.size(); i++){
        if (mesh.active_polygon[i]==true){
            mesh.numberCell2D++;
        }
    }
}

// funzione che crea la mesh
void CreateMesh(vector<Fracture>& fractures, double tol, DFN &dfn, Meshes &meshesV){

    for (unsigned int i=0; i<fractures.size(); i++)
    {
        Fracture fra = fractures[i];
        vector<unsigned int> allidT(fra.passingTraces.size() + fra.notPassingTraces.size());

        for(unsigned int t = 0; t < fra.passingTraces.size(); t++){
            allidT[t] = fra.passingTraces[t];
        }

        for(unsigned int t = 0; t < fra.notPassingTraces.size(); t++){
            allidT[t + fra.passingTraces.size()] = fra.notPassingTraces[t];
        }

        PolygonalMesh mesh;

        const unsigned int numVertices = fra.verticesCoordinates.cols();

        mesh.idVertices.resize(numVertices);
        mesh.coordVertices.resize(numVertices);

        mesh.idEdges.resize(numVertices);
        mesh.active_edge.resize(numVertices);
        mesh.extremitiesEdges.resize(numVertices);
        mesh.nearPolygons.resize(numVertices);
        mesh.newedge.resize(numVertices);

        for (unsigned int v=0; v < numVertices; v++){

            mesh.idVertices[v] = v;
            mesh.coordVertices[v] = fra.verticesCoordinates.col(v);

            mesh.idEdges[v] = v;
            mesh.extremitiesEdges[v] = {v, (v + 1) % numVertices};
            mesh.active_edge[v] = true;
            mesh.nearPolygons[v] = {-1,-1};
            mesh.newedge[v] = {-1,-1};

        }

        mesh.idPolygon.resize(1);
        mesh.idPolygon[0] = 0;

        mesh.verticesPolygons.resize(1);
        mesh.verticesPolygons[0].resize(numVertices);
        iota(mesh.verticesPolygons[0].begin(),  mesh.verticesPolygons[0].end(), 0);

        mesh.edgesPolygons.resize(1);

        mesh.edgesPolygons[0].resize(numVertices);
        iota(mesh.edgesPolygons[0].begin(),  mesh.edgesPolygons[0].end(), 0);

        mesh.active_polygon.resize(1);
        mesh.active_polygon[0] = true;


        CutFracture(mesh, dfn, mesh.idPolygon[0], allidT, tol);

        CorrectMesh(mesh);

        meshesV.meshes.push_back(mesh);


    }
}

// funzione di stampa delle mesh
void PrintMeshes (const string& fileName, Meshes &meshesV){ // primo file di ouput, con le informazioni sulle tracce
    ofstream ofstr(fileName); // se il file non esiste, lo crea
    ofstr << "# Number of Meshes" << endl;
    ofstr << meshesV.meshes.size() << endl;
    for (PolygonalMesh& meh: meshesV.meshes){

        ofstr << "# Number of Vertices" << endl;
        ofstr << meh.numberCell0D << endl;
        ofstr << "# IdVertici; X1; Y1; Z1" << endl;
        for (unsigned int i=0;i<meh.numberCell0D;i++)
            ofstr << meh.idVertices[i] << "; " << meh.coordVertices[i][0] << "; " << meh.coordVertices[i][1]
                  << "; " << meh.coordVertices[i][2] << endl;

        ofstr << "# Number of Edges" << endl;

        unsigned int numEdge = 0;
        for (unsigned int i=0;i<meh.numberCell1D;i++){
            if (meh.active_edge[i]==true){
                numEdge++;
            }
        }

        ofstr << numEdge << endl;
        ofstr << "# IdEdges; IdVertices1, IdVertices2" << endl;
        for (unsigned int i=0;i<meh.numberCell1D;i++){
            if (meh.active_edge[i]==true){
                ofstr << meh.idEdges[i] << "; " << meh.extremitiesEdges[i][0] << "; " << meh.extremitiesEdges[i][1] << endl;
            }
        }

        ofstr << "# Number of Polygons" << endl;
        unsigned int numPoly = 0;
        for (unsigned int i=0;i<meh.numberCell2D;i++){
            if (meh.active_polygon[i]==true){
                numPoly++;
            }
        }
        ofstr << numPoly << endl;
    }
    ofstr.close();
}

}
