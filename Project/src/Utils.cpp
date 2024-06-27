#include "Utils.hpp"
#include "DFN.hpp"
#include <iostream>
#include <fstream>
#include <sstream>


using namespace std;
using namespace Eigen;

namespace DiscreteFractureNetworkLibrary
{

bool ImportFractures(const string& filePath, DFN& dfn, double tol) {
    ifstream file(filePath);

    if (file.fail()) {
        cerr << "Error in opening the file" << endl;
        return false;
    }

    string line;
    while (getline(file, line)) { // leggi finché trovi una linea che non è un commento
        if (line[0] != '#') {
            break;
        }
    }

    if (line.empty() || line[0] == '#') {
        cerr << "Invalid file format" << endl;
        return false;
    }

    istringstream convertFrac(line); // salvo il numero di fratture
    convertFrac >> dfn.fractureNumber;

    // ridefinisco la dimensione del vettore di fratture in dfn
    dfn.Fractures.resize(dfn.fractureNumber);

    for (unsigned int i = 0; i < dfn.fractureNumber; i++) {
        Fracture fra;
        getline(file, line);             // linea commenti
        getline(file, line, ';');

        istringstream convertId(line);   // salvo l'id
        // salvo l'id nella frattura
        unsigned int idFra;
        convertId >> idFra;
        fra.setId(idFra);

        getline(file, line);
        istringstream convertVert(line);   // numero di vertici
        // salvo il numero di vertici nella frattura
        unsigned int numVert;
        convertVert >> numVert;
        fra.setVerticesNumber(numVert);

        // trovo la matrice delle coordinate dei vertici e la salvo
        MatrixXd VerMatrix(3, numVert);

        getline(file, line); // linea commenti

        for (unsigned int j = 0; j < 3; j++) {
            getline(file, line);
            replace(line.begin(), line.end(), ';', ' ');
            istringstream convertLine(line);

            for (unsigned int k = 0; k < numVert; k++) {
                convertLine >> VerMatrix(j, k);
            }
        }

        fra.setVerticesCoordinates(VerMatrix);

        // calcolo il baricentro della frattura e lo salvo
        Vector3d barycenter = VerMatrix.rowwise().mean();
        fra.setBarycenter(barycenter);

        // calcolo e salvo il raggio
        double radius = 0;
        for(unsigned int i = 0; i<fra.verticesNumber; i++)
        {
            double d = (barycenter - fra.verticesCoordinates.col(i)).norm();
            if(d - radius > tol){
                radius = d;
            }
        }
        fra.setRadius(radius);

        // calcolo e salvo la normale
        Vector3d v1 = fra.verticesCoordinates.col(0);
        Vector3d v2 = fra.verticesCoordinates.col(1);
        Vector3d v3 = fra.verticesCoordinates.col(2);
        Vector3d d1 = v1-v2;
        Vector3d d2 = v3-v2;
        Vector3d normal = d1.cross(d2);
        fra.setNormal(normal);

        // calcolo e salvo il piano
        double d = (fra.normal[0]*fra.barycenter[0] + fra.normal[1]*fra.barycenter[1] + fra.normal[2]*fra.barycenter[2]);
        Vector4d plane = {fra.normal[0], fra.normal[1], fra.normal[2], d};
        fra.setPlane(plane);

        // salvo la frattura nel vettore delle fratture
        dfn.Fractures[i] = fra;
    }

    file.close();
    return true;
}

double PointsDistance(const Vector3d p1, const Vector3d p2)
{
    double distance = sqrt((p1-p2).transpose() * (p1-p2));   // distanza fra due punti in norma euclidea
    return distance;
}

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

// retta di intersezione tra due piani
array<Vector3d, 2> LineIntersection(Fracture &f1, Fracture &f2) {

    // Vettore direzione della retta (prodotto vettoriale dei normali dei due piani)
    Vector3d n1 = f1.normal;
    Vector3d n2 = f2.normal;
    Vector4d p1 = f1.plane;
    Vector4d p2 = f2.plane;
    Vector3d direzione = n1.cross(n2);


    // Troviamo un punto sulla retta di intersezione. Metto a sistema i due piani e l'equazione della retta eguagliata ad un valore a caso (0 per semplicità)
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

// intersezioni di una frattura con una retta ( che poi sarà la retta di intersez dei due piani per esempio)
pair<bool,Vector2d> InterFractureLine(Fracture &f1, array<Vector3d, 2> &r, double tol){
    double tol2 = max(10*numeric_limits<double>::epsilon(), tol*tol);  // numeric_limits<double>::epsilon() calcola la precisione di macchina per il tipo double
    Vector3d punto = r[0];
    Vector3d direzioneRetta = r[1];
    double den = 1/(direzioneRetta.dot(direzioneRetta));
    MatrixXd vertices = f1.verticesCoordinates;

    unsigned int count = 0;
    Vector2d beta;  // creo il vettore delle ascisse curvilinee
    bool con = false;


    for (unsigned int i = 0; i < f1.verticesCoordinates.cols(); i++) {   // f1.verticesCoordinates.cols() è il numero di vertici di f1
        // Calcolo degli indici dei vertici del segmento di linea
        unsigned int j = (i + 1) % f1.verticesCoordinates.cols();  // (i+1) MODULO numeroDiVertici: fa lo stesso di i+1 fintanto che i+1< numeroDiVertici. Quando i+1=numeroDiVertici, j=0 e quindi connetto l'ultimo vertice con il primo
        Vector3d ver1 = f1.verticesCoordinates.col(i);
        Vector3d ver2 = f1.verticesCoordinates.col(j);
        Vector3d direzioneLato = ver2 - ver1;

        Vector3d prodotto = direzioneLato.cross(direzioneRetta);  // prodotto scalare tra lato del poligono e retta

        if (prodotto.norm() > tol2){   // se non sono paralleli
            double alpha = ((punto-ver1).cross(direzioneRetta)).dot(prodotto)/(prodotto.dot(prodotto));  // VEDI (!)

            if (alpha >= -tol && alpha < 1-tol){
                double betaTemp = ((ver1-punto).cross(direzioneLato)).dot(-prodotto)/(prodotto.dot(prodotto));  // VEDI (!!); BetaTemp = beta temporaneo
                beta[count] = betaTemp;  // aggiungo betaTemporaneo al vettore dei beta
                count++;

                if(count==2){
                    sort(beta.begin(), beta.end());   // se il poligono interseca la retta due volte, (può farlo 0 o 1 o 2), ordino le acsisse, in ordine crescente
                    con = true;
                    break;
                }
            }

        }

         /* (!) SPIEGAZIONE ALPHA
          * La retta r si parametrizza come punto + t * direzioneRetta, con t parametro.
          * Il lato si parametrizza come ver1 + alpha * direzioneLato, con alpha parametro --> se alpha sta tra 0 e 1, l'intersezione
          * tra la retta che contiene il lato e la retta r avviene dentro al lato (quello che vogliamo)
          * trovo alpha ponendo l'intersezione: punto + t * direzioneRetta = ver1 + alpha * direzioneLato e prendendo
          * da entrambi i lati il prodotto vettoriale con direzioneRetta per eliminare t
          *
          * (!!) SPIEGAZIONE BETA
          * La retta r si parametrizza come punto + betaTemp * direzioneRetta, con betaTemp parametro -->
          * --> beta rappresenta l'ascissa che indica la posizione dell'intersezione tra r e lato lungo r
          * Il lato si parametrizza come ver1 + alpha * direzioneLato, con alpha parametro
          * trovo betaTemp ponendo l'intersezione: punto + betTemp * direzioneRetta = ver1 + alpha * direzioneLato e prendendo
          * da entrambi i lati il prodotto vettoriale con direzioneLato per eliminare alpha
         */

        // NORMALIZZARE VETTORE DIREZIONE RETTA ALL'INIZIO COSì DA FARE MENO DIVISIONI
        // vd libro Calafiore - Optimization Models
        else{  // VEDI (!!!)
            // const double b0 = (ver1[0] - punto[0])/ direzioneRetta[0];
            // const double b1 = (ver1[1] - punto[1])/ direzioneRetta[1];
            //const double b2 = (ver1[2] - punto[2])/ direzioneRetta[2];
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

        /* (!!!) SPIEGAZIONE CONTROLLO PARALLELISMO (else):
         * b0, b1, e b2 sono i parametri ai quali le coordinate x, y, e z del vertice ver1 stanno lungo r (definita da punto e direzioneRetta).
         * I parametri sono calcolati proiettando ver1 su r in ogni direzione della base (x, y, z).
         * La condizione "if (abs(b0-b1) < tol2 && abs(b0 - b2) < tol2)" controlla se b0, b1, e b2 sono approssimativamente uguali.
         * Se lo sono, signifca che ver1 giace su r (fai disegnino e lo vedi).
         * betaTemp ha lo stesso signifcato che ha in (!!).
         * "((ver1 - punto).dot(direzioneRetta))" proietta il vettore (ver1 - punto) sulla direzione direzioneRetta;
         * ; il denominatore di betaTemp serve per normalizzare.
         *
         * Quindi, riassumento l'else controlla il caso in cui un lato ed il segmento siano paralleli e in particolare
         * sovrapposti. Salva il vertice in beta con il metodo dell'ascissa coerente con il passo (!!).
        */


    }
    return {con, beta};
    // quindi per ogni frattura ottengo due info:
    // * con: true se c'è intersezione con la retta
    // * beta che contiene le ascisse dell'intersezione, in caso
}

void FindTraces(vector<Fracture> &fractures, double tol, DFN &dfn)
{
    unsigned int count = 0;
    double tol2 = max(10*numeric_limits<double>::epsilon(), tol*tol);

    for (unsigned int i=0; i<fractures.size(); i++){ //devo controllare ogni coppia di fratture possibile
        for(unsigned int j=i+1; j<fractures.size(); j++){
            Fracture &f1 = fractures[i];
            Fracture &f2 = fractures[j];
            array<Vector3d, 2> estremiTraccia;

            // controlli ed esclusioni
            if (Parallel(f1,f2,tol)){
                continue;  // fa andare alla prossima iterazione del for interno: esclude la coppia i,j
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


            if (betaf1[1]<betaf2[0]+tol || betaf2[1]<betaf1[0]+tol){ // caso in cui sono disgiunte lungo la retta (viene "prima" f1 o f2)
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

            array<bool,2> tips = {false, false}; // se resta così è passante per entrambi
            if (estremiTraccia != puntiInterf1){
                tips[0] = true;  // in questo caso è non passante
            }
            if (estremiTraccia != puntiInterf2){
                tips[1] = true;
            }

            // calcolo la lunghezza ed escludo il caso in cui l'intersezione consista di un unico punto
            double length = PointsDistance(estremiTraccia[0],estremiTraccia[1]);
            if (length > tol){
                // creo la traccia
                Trace tr;
                tr.setIdTrace(count);
                tr.setExtremitiesCoordinates(estremiTraccia);
                tr.setFracturesIds({f1.id,f2.id});
                tr.setLength(length);
                tr.setTips(tips);
                tr.setLineTrace(rettaf1xf2);

                if(tips[0]){
                    f1.notPassingTraces.push_back(tr.getIdTrace());
                }
                else{
                    f1.passingTraces.push_back(tr.getIdTrace());

                }
                if(tips[1]){
                    f2.notPassingTraces.push_back(tr.getIdTrace());

                }
                else{
                    f2.passingTraces.push_back(tr.getIdTrace());
                }

                dfn.Traces.push_back(tr);
                count++;
            }

        }
    }
}

void printGlobalResults (const string& fileName, vector<Trace>& traces){ // primo file di ouput, con le informazioni sulle tracce
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

void printLocalResults (const string& fileName, vector<Fracture>& fractures, const vector<Trace>& traces)
{
    // secondo file di ouput, con le informazioni sulle fratture e sulle tracce corrispondenti
    ofstream ofstr(fileName);
    bool firstTime;
    for (Fracture& fr:fractures){
        firstTime=true;
        ofstr << "# FractureId; NumTraces" << endl;
        ofstr << fr.id << "; " << (fr.passingTraces.size())+(fr.notPassingTraces.size()) << endl;
        sortTracesByLength(fr.passingTraces, traces); // ordino le tracce passanti per lunghezza decrescente
        sortTracesByLength(fr.notPassingTraces, traces); // ordino le tracce non passanti per lunghezza decrescente
        for (unsigned int idTrace:fr.passingTraces){
            if (firstTime){
                ofstr << "# TraceId; Tips; Length" << endl;
                firstTime=false;
            }
            ofstr << idTrace << "; " << "false; " << traces[idTrace].length << endl; //sto stampando prima tutte quelle passanti, quindi avranno tutte tips=false

        }
        for (unsigned int idTrace:fr.notPassingTraces){
            if (firstTime){
                ofstr << "# TraceId; Tips; Length" << endl;
                firstTime=false;
            }
            ofstr << idTrace << "; " << "true; " << traces[idTrace].length << endl; // sto stampando tutte quelle non passanti, quindi tips = true
        }
    }
    ofstr.close();
}

bool compareTraceLength(const unsigned int& id1, const unsigned int& id2, const std::vector<Trace>& traces) {
    return traces[id1].length > traces[id2].length; // > per ordine decrescente
}

void sortTracesByLength(vector<unsigned int>& vecIdTraces, const vector<Trace>& traces) {
    sort(vecIdTraces.begin(), vecIdTraces.end(), [&traces](const unsigned int& id1, const unsigned int& id2) {
        return compareTraceLength(id1, id2, traces);
    });
}
}

namespace PolygonalMeshLibrary
{

int PositionVert(const Vector3d& point, array<Vector3d,2> retta) {
    Vector3d relPoint = point - retta[0];
    Vector3d crossProduct = relPoint.cross(retta[1]);
    double dotProduct = crossProduct.dot(Vector3d(0, 0, 1));
    if (dotProduct > 0) {
        return 1;
    }
    if (dotProduct < 0){
        return -1;
    }
    return 0;
}

pair<Vector2d, Vector2i> InterFractureLine2(PolygonalMesh &mesh, vector<unsigned int> &vertIds, array<Vector3d, 2> &r, double tol){
    double tol2 = max(10*numeric_limits<double>::epsilon(), tol*tol);
    Vector3d punto = r[0];
    Vector3d direzioneRetta = r[1];

    unsigned int count = 0;
    Vector2d beta;
    Vector2i idVertEdgesTrace;

    for (unsigned int i = 0; i < vertIds.size(); i++) {
        // Calcolo degli indici dei vertici del segmento di linea
        unsigned int j = (i + 1) % vertIds.size();
        Vector3d ver1 = mesh.coordVertices[vertIds[i]];
        Vector3d ver2 = mesh.coordVertices[vertIds[j]];
        Vector3d direzioneLato = ver2 - ver1;

        Vector3d prodotto = direzioneLato.cross(direzioneRetta);

        if (prodotto.norm() > tol2){
            double alpha = ((punto-ver1).cross(direzioneRetta)).dot(prodotto)/(prodotto.dot(prodotto));

            if (alpha > 0 && alpha < 1){
                double betaTemp = ((ver1-punto).cross(direzioneLato)).dot(-prodotto)/(prodotto.dot(prodotto));
                beta[count] = betaTemp;
                idVertEdgesTrace[count] = i;
                count++;

                if(count==2){
                    break;
                }
            }
        }
        else{
            const double b0 = (ver1[0] - punto[0])/ direzioneRetta[0];
            const double b1 = (ver1[1] - punto[1])/ direzioneRetta[1];
            const double b2 = (ver1[2] - punto[2])/ direzioneRetta[2];
            if (abs(b0-b1) < tol2 && abs(b0 - b2) < tol2){
                double betaTemp = ((ver1-punto).dot(direzioneRetta))/direzioneRetta.dot(direzioneRetta);

                beta[count] = betaTemp;
                idVertEdgesTrace[count] = i;
                count++;

                if(count==2){
                    break;
                }
            }
        }
    }
    return {beta, idVertEdgesTrace};
}

void cutFracture(PolygonalMesh &mesh, Fracture& fra, vector<Vector3d>& vertCoor, vector<Trace>& traces, unsigned int countIdV, unsigned int countIdE, unsigned int countIdP, double tol){
    // unsigned int countIdV = mesh.numVertices;
    // unsigned int countIdE = mesh.numEdges;
    // unsigned int countIdP = mesh.numPolygons;

    vector<unsigned int> vertIds;
    vertIds.resize(vertCoor.size());

    for (unsigned int i=0; i<vertCoor.size(); i++){

        if (find(mesh.coordVertices.begin(), mesh.coordVertices.end(), vertCoor[i]) != mesh.coordVertices.end()){
            vertIds[i] = distance(mesh.coordVertices.begin(), find(mesh.coordVertices.begin(), mesh.coordVertices.end(), vertCoor[i]));
        }
        else{
            vertIds[i] = countIdV;
            countIdV++;
            mesh.idVertices.push_back(vertIds[i]);
            mesh.coordVertices.push_back(vertCoor[i]);
        }
    }

    vector<unsigned int> edgeIds;
    edgeIds.resize(vertIds.size());

    for (unsigned int i=0; i<vertIds.size(); i++){
        unsigned int j = (i + 1) % vertIds.size();

        Vector2i lato = {vertIds[i],vertIds[j]};
        Vector2i lato2 = {vertIds[j],vertIds[i]};

        if (find(mesh.extremitiesEdges.begin(), mesh.extremitiesEdges.end(), lato) != mesh.extremitiesEdges.end()){
            edgeIds[i] = distance(mesh.extremitiesEdges.begin(),find(mesh.extremitiesEdges.begin(), mesh.extremitiesEdges.end(), lato));
            mesh.nearPolygons[i][1] = countIdP;
        }
        else if (find(mesh.extremitiesEdges.begin(), mesh.extremitiesEdges.end(), lato2) != mesh.extremitiesEdges.end()){
            edgeIds[i] = distance(mesh.extremitiesEdges.begin(),find(mesh.extremitiesEdges.begin(), mesh.extremitiesEdges.end(), lato2));
            mesh.nearPolygons[i][1] = countIdP;
        }
        else{
            edgeIds[i] = countIdE;
            countIdE++;
            mesh.idEdges.push_back(edgeIds[i]);
            mesh.extremitiesEdges.push_back(lato);
            mesh.active_edge.push_back(true);
            Vector2i vicini = {countIdP,-1};
            mesh.nearPolygons.push_back(vicini);
            mesh.newedge.push_back({-1,-1});
        }
    }

    mesh.idPolygon.push_back(countIdP);
    countIdP++;
    mesh.verticesPolygons.push_back(vertIds);
    mesh.edgesPolygons.push_back(edgeIds);

    if (traces.empty()){
        mesh.active_polygon.push_back(true);
    }
    else {
        // 1) Utilizzare la funzione InterFractureLine per trovare i punti di intersezione della
        //    traccia (passante o non passannte con l'insieme di vertici in considerazione)
        // 2) Dividere i vertici a seconda che si trovino a destra o a sinistra della traccia in due sottoinsiemi
        // 3) Dividere le tracce rimanenti
        // 4) Richiamare la funzione stessa con le nuove sottofratture

        mesh.active_polygon.push_back(false);

        Trace &tra = traces.front();
        traces.erase(traces.begin());

        array<Vector3d, 2> rettaTraccia = tra.lineTrace;

        // divido i vertici nelle due sotto figure
        vector<Vector3d> vertCoorSub1;
        vector<Vector3d> vertCoorSub2;

        auto [beta, idVertEdgesTrace] = InterFractureLine2(mesh, vertIds, rettaTraccia, tol);
        array<Vector3d, 2> puntiInterFraTra = {rettaTraccia[0] + (beta[0]*rettaTraccia[1]), rettaTraccia[0] + (beta[1]*rettaTraccia[1])};

        Vector2i idInterVert = {};
        for (unsigned int i = 0; i < puntiInterFraTra.size(); i++) {
            idInterVert[i] = countIdV;
            mesh.idVertices.push_back(countIdV);
            countIdV++;
            mesh.coordVertices.push_back(puntiInterFraTra[i]);
        }

        for (unsigned int i = 0; i < idVertEdgesTrace.size(); i++) {
            unsigned int j = (idVertEdgesTrace[i] + 1) % vertIds.size();
            Vector2i lato = {vertIds[idVertEdgesTrace[i]],vertIds[j]};
            for (unsigned int ed = 0; ed < mesh.idEdges.size(); ed++) {
                if (mesh.extremitiesEdges[ed]==lato){

                    Vector2i idInterEdg = {};

                    Vector2i latoNew1 = {vertIds[idVertEdgesTrace[i]],idInterVert[i]};
                    idInterEdg[0] = countIdE;
                    mesh.idEdges.push_back(countIdE);
                    countIdE++;
                    mesh.extremitiesEdges.push_back(latoNew1);
                    mesh.active_edge.push_back(true);
                    mesh.newedge.push_back({-1,-1});

                    Vector2i latoNew2 = {idInterVert[i],vertIds[j]};
                    idInterEdg[1] = countIdE;
                    mesh.idEdges.push_back(countIdE);
                    countIdE++;
                    mesh.extremitiesEdges.push_back(latoNew2);
                    mesh.active_edge.push_back(true);
                    mesh.newedge.push_back({-1,-1});


                    mesh.active_edge[ed] = false;
                    mesh.newedge[ed] = idInterEdg;

                    for (unsigned int p = 0; p < 2; p++) {
                        int idP = mesh.nearPolygons[ed][i];
                        if (idP != -1){
                            vector<unsigned int> v = mesh.edgesPolygons[idP];
                            v.erase(find(v.begin(),v.end(),ed));
                            mesh.verticesPolygons[idP].push_back(idInterVert[0]);
                            mesh.verticesPolygons[idP].push_back(idInterVert[1]);
                            v.push_back(idInterEdg[0]);
                            v.push_back(idInterEdg[1]);
                        }
                    }
                }
            }
        }

        bool where = true;

        for (unsigned int i = 0; i < vertIds.size(); i++) {
            Vector3d ver1 = mesh.coordVertices[i];

            if (where){
                vertCoorSub1.push_back(ver1);

                if (i - idVertEdgesTrace[0] == 0){
                    where = false;
                    vertCoorSub1.push_back(puntiInterFraTra[0]);
                    vertCoorSub1.push_back(puntiInterFraTra[1]);
                }
            }
            else if (!where){
                vertCoorSub2.push_back(ver1);

                if (i - idVertEdgesTrace[1] == 0){
                    where = true;
                    vertCoorSub2.push_back(puntiInterFraTra[1]);
                    vertCoorSub2.push_back(puntiInterFraTra[0]);
                }
            }


            /*
            if (((PointsDistance(ver1, vertEdgesTrace[0][0]) < tol && PointsDistance(ver2, vertEdgesTrace[0][1]) < tol)
                 || (PointsDistance(ver2, vertEdgesTrace[0][1]) < tol && PointsDistance(ver1, vertEdgesTrace[0][0]) < tol))){
                if (where){
                    vertCoorSub1.push_back(ver1);
                    vertCoorSub1.push_back(puntiInterFraTra[0]);
                    where = false;
                }
                else if (!where){
                    vertCoorSub2.push_back(ver1);
                    vertCoorSub2.push_back(puntiInterFraTra[0]);
                    where = true;
                }
            }
            else if (((PointsDistance(ver1, vertEdgesTrace[1][0]) < tol && PointsDistance(ver2, vertEdgesTrace[1][1]) < tol)
                      || (PointsDistance(ver2, vertEdgesTrace[1][1]) < tol && PointsDistance(ver1, vertEdgesTrace[1][0]) < tol))){
                if (where){
                    vertCoorSub1.push_back(ver1);
                    vertCoorSub1.push_back(puntiInterFraTra[1]);
                    where = false;
                }
                else if (!where){
                    vertCoorSub2.push_back(ver1);
                    vertCoorSub2.push_back(puntiInterFraTra[1]);
                    where = true;
                }
            }
            else{
                if (where){
                    vertCoorSub1.push_back(ver1);
                }
                else if (!where){
                    vertCoorSub2.push_back(ver1);
                }
            }
            */
        }

        /*
        bool precedente = PositionVert(vertCoor[0],rettaTraccia) > 0;
        for (Vector3d& vert : vertCoor){

            int where = PositionVert(vert,rettaTraccia);

            if (where > tol){
                if (!precedente){
                    vertCoorSub1.push_back(puntiInterFraTra[0]);
                    vertCoorSub1.push_back(puntiInterFraTra[1]);
                    precedente = true;
                    }
                vertCoorSub1.push_back(vert);
            }
            else if (where < tol) {
                if (precedente) {
                    vertCoorSub2.push_back(puntiInterFraTra[1]);
                    vertCoorSub2.push_back(puntiInterFraTra[0]);
                    precedente = false;
                }
                vertCoorSub2.push_back(vert);
            }
        }
*/
        // divido le tracce
        vector<Trace> tracesSub1;
        vector<Trace> tracesSub2;

        for (Trace& tra : traces){
            array<Vector3d, 2> estr = tra.extremitiesCoordinates;

            int w1 = PositionVert(estr[0],rettaTraccia);
            int w2 = PositionVert(estr[1],rettaTraccia);

            if (w1>tol && w2>tol){
                tracesSub1.push_back(tra);
            }
            if (w1<tol && w2<tol){
                tracesSub2.push_back(tra);
            }
            else {
                tracesSub1.push_back(tra);
                tracesSub2.push_back(tra);
            }
        }

        // creo le due sotto fratture
        Fracture fraSub1;
        Fracture fraSub2;

        MatrixXd VerMatrixSub1(3, vertCoorSub1.size());
        for (unsigned int j = 0; j < vertCoorSub1.size(); j++) {
            for (unsigned int k = 0; k < 3; k++) {
                VerMatrixSub1(k,j) = vertCoorSub1[j][k];
            }
        }
        fraSub1.setVerticesCoordinates(VerMatrixSub1);

        MatrixXd VerMatrixSub2(3, vertCoorSub2.size());
        for (unsigned int j = 0; j < vertCoorSub2.size(); j++) {
            for (unsigned int k = 0; k < 3; k++) {
                VerMatrixSub2(k,j) = vertCoorSub2[j][k];
            }
        }
        fraSub2.setVerticesCoordinates(VerMatrixSub2);

        // richiamo la funzione per le sottofratture fino alla fine delle tracce
        cutFracture(mesh, fraSub1, vertCoorSub1, tracesSub1, countIdV, countIdE, countIdP, tol);
        cutFracture(mesh, fraSub2, vertCoorSub2, tracesSub2, countIdV, countIdE, countIdP, tol);

    }

}

void correctMesh(PolygonalMesh& mesh){

    mesh.NumberCell0D = mesh.idVertices.size();

    for (unsigned int i=0; i<mesh.idEdges.size(); i++){
        if (mesh.active_edge[i]==true){
            mesh.NumberCell1D++;
        }
    }

    for (unsigned int i=0; i<mesh.idPolygon.size(); i++){
        if (mesh.active_polygon[i]==true){
            mesh.NumberCell2D++;
        }
    }

}

void CreateMesh(vector<Fracture>& fractures, double tol, DFN &dfn, plm &plm){

    for (unsigned int i=0; i<fractures.size(); i++){
        Fracture fra = fractures[i];
        vector<unsigned int> allidT;
        allidT.resize(fra.passingTraces.size() + fra.notPassingTraces.size());

        // Usa la funzione merge per fondere i due vettori in uno
        merge(fra.passingTraces.begin(), fra.passingTraces.end(), fra.notPassingTraces.begin(), fra.notPassingTraces.end(), allidT.begin());

        vector<Trace> allTraces;
        for (unsigned int t=0; t<allidT.size(); t++){
            for (unsigned int p=0; p<dfn.Traces.size(); p++){
                if (dfn.Traces[p].idTrace==t){
                    allTraces.push_back(dfn.Traces[p]);
                }
            }
        }

        vector<Vector3d> vertCoord;
        for (unsigned int c=0; c<fra.verticesCoordinates.cols(); c++){
            vertCoord.push_back(fra.verticesCoordinates.col(c));
        }

        /*
        for (unsigned int k=0; k<allTraces.size(); k++){
            cout << dfn.Traces[allTraces[k]].Tips[0] << endl;
            cout << dfn.Traces[allTraces[k]].length << endl;
        }*/


        PolygonalMesh Mesh;

        unsigned int countIdV = 0;
        unsigned int countIdE = 0;
        unsigned int countIdP = 0;

        cutFracture(Mesh, fra, vertCoord, allTraces, countIdV, countIdE, countIdP, tol);

        correctMesh(Mesh);

        plm.meshes.push_back(Mesh);


    }
}

void tryOutput (const string& fileName, plm &plm){ // primo file di ouput, con le informazioni sulle tracce
    ofstream ofstr(fileName); // se il file non esiste, lo crea
    ofstr << "# Number of Meshes" << endl;
    ofstr << plm.meshes.size() << endl;
    for (PolygonalMesh& meh: plm.meshes){
        ofstr << "# Number of Vertices" << endl;
        ofstr << meh.NumberCell0D << endl;
        ofstr << "# IdVertici; X1; Y1; Z1" << endl;
        for (unsigned int i=0;i<meh.NumberCell0D;i++)
            ofstr << meh.idVertices[i] << "; " << meh.coordVertices[i][0] << "; " << meh.coordVertices[i][1]
                  << "; " << meh.coordVertices[i][2] << endl;

        ofstr << "# Number of Edges" << endl;
        ofstr << meh.NumberCell1D << endl;
        ofstr << "# IdEdges; IdVertices1, IdVertices2" << endl;
        for (unsigned int i=0;i<meh.NumberCell1D;i++)
            ofstr << meh.idEdges[i] << "; " << meh.extremitiesEdges[i][0] << "; " << meh.extremitiesEdges[i][1] << endl;

        ofstr << "# Number of Polygons" << endl;
        ofstr << meh.NumberCell2D << endl;
    }
    ofstr.close();
}

}





