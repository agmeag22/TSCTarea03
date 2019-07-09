#include <iostream>
#include "math_tools.h"
#include "classes.h"
#include "tools.h"
#include "sel.h"

int main(int argc, char *argv[])
{
    char filename[150];
    strcpy(filename,"");

    vector<Matrix> localKs;
    vector<Vector> localbs;
    vector<Vector> Ts;
    Matrix K;
    Vector b;

    cout << "IMPLEMENTACI"<<char(224)<<"N DEL M"<<char(144)<<"TODO DE LOS ELEMENTOS FINITOS\n"
         << "\t- ECUACIONES DE NAVIER-STOKES\n" << "\t- 2 DIMENSIONES\n"
         << "\t- FUNCIONES DE FORMA LINEALES\n" << "\t- PESOS DE GALERKIN\n"
         << "*********************************************************************************\n\n";

    mesh m;
    leerMallayCondiciones(m,filename);
    cout << "Datos obtenidos correctamente\n********************\n";

    crearSistemasLocales(m,localKs,localbs);
    showKs(localKs); showbs(localbs);
    cout << "******************************\n";

    zeroes(K,4*m.getSize(NODES));
    zeroes(b,4*m.getSize(NODES));
    ensamblaje(m,localKs,localbs,K,b);
    showMatrix(K); showVector(b);
    cout << "******************************\n";
    //cout << K.size() << " - "<<K.at(0).size()<<"\n";
    //cout << b.size() <<"\n";

    applyDirichlet(m,K,b);
    showMatrix(K); showVector(b);
    cout << "******************************\n";
    //cout << K.size() << " - "<<K.at(0).size()<<"\n";
    //cout << b.size() <<"\n";

    //zeroes(T,b.size());
    calculate(K,b,Ts,3*m.getSize(NODES)-m.getSize(DIRICHLET));
    showbs(Ts);

    //cout << "La respuesta es: \n";
    //showVector(T);

    writeResults(m,Ts,filename);

    return 0;
}
