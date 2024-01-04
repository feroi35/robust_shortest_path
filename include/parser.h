// parser.h
#ifndef PARSER_H
#define PARSER_H

#include <fstream>
#include <iostream>
#include <vector>
#include <limits>

// Magic tricks to have CPLEX behave well:
#ifndef IL_STD
#define IL_STD
#endif
#include <cstring>
#include <ilcplex/ilocplex.h>
ILOSTLBEGIN
// End magic tricks

struct Arc {
    IloInt i;
    IloInt j;
    IloNum d;
    IloNum D;
};

const double undefinedValue = std::numeric_limits<double>::quiet_NaN();


struct Instance {
    char* name;
    IloInt n;
    IloInt s;
    IloInt t;
    IloNum S;
    IloNum d1;
    IloNum d2;
    IloNumArray p;
    IloNumArray ph;
    std::vector<Arc> mat;

    std::vector<std::vector<double>> d;
    std::vector<std::vector<double>> D;

    std::vector<IloInt> sol;

    Instance(IloEnv env, char filename[]);
    void display() const;
    // void exportSol(char filename[]) const;
};


#endif
