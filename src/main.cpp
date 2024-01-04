#include <ilcplex/ilocplex.h>
#include <fstream>
#include <iostream>
#include <vector>
#include "main.h"

ILOSTLBEGIN


struct Vertice {
    IloInt i;
    IloInt j;
    IloNum d;
    IloNum D;
};


void getData(char filename[], IloInt& n, IloInt& s, IloInt& t, IloNum& S,
    IloNum& d1, IloNum& d2, IloNumArray& p, IloNumArray& ph, std::vector<Vertice>& Mat) {

    char readChar;
    int readInt;

    ifstream file(filename);

    if (!file) {
        std::cerr << "Error when opening file" << std::endl;
        exit(1);
    }

    file >> readChar >> readChar >> n;
    file >> readChar >> readChar >> s;
    file >> readChar >> readChar >> t;
    file >> readChar >> readChar >> S;
    file >> readChar >> readChar >> readChar >> d1;
    file >> readChar >> readChar >> readChar >> d2;
    file >> readChar >> readChar;
    for (int i=0; i<n; i++) {
        file >> readChar >> readInt;
        p.add(readInt);
    }
    // ] ph =
    file >> readChar >> readChar >> readChar >> readChar;
    for (int i=0; i<n; i++) {
        file >> readChar >> readInt;
        ph.add(readInt);
    }
    // ] Mat = [
    file >> readChar >> readChar >> readChar >> readChar >> readChar >> readChar;

    while (readChar != ']') {
        Vertice v;
        file >> v.i;
        file >> v.j;
        file >> v.d;
        file >> v.D;
        file >> readChar; // either ';' or ']'
        Mat.push_back(v);
    }
    file.close();
}


int main(int argc, char **argv) {

    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <filename>" << std::endl;
        exit(1);
    }
    char* filename = argv[1];
    std::cout << filename << std::endl;

    // Data reading
    IloEnv env;

    IloInt n;
    IloInt s;
    IloInt t;
    IloNum S;
    IloNum d1;
    IloNum d2;
    IloNumArray p=IloNumArray(env);
    IloNumArray ph=IloNumArray(env);
    std::vector<Vertice> Mat;

    getData(filename, n, s, t, S, d1, d2, p, ph, Mat);
    std::cout << "n = " << n << std::endl;

    return 0;
}