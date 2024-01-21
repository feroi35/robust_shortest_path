// parser.h
#ifndef PARSER_H
#define PARSER_H

#include <vector>
#include <ilcplex/ilocplex.h>

ILOSTLBEGIN

struct Arc {
    IloInt i;
    IloInt j;
    IloNum d;
    IloNum D;
};


struct Instance {
    std::string name;
    IloInt n; // nombre de noeuds
    IloInt n_arc; // nombre d'arcs
    IloInt s; // ville d'entrée
    IloInt t; // ville d'entrée
    IloNum S; // contrainte poids de villes
    IloNum d1; // max incertitude poids arcs
    IloNum d2; // max incertitude poids villes
    IloNumArray p; // poids des villes
    IloNumArray ph; // incertitudes poids des villes
    IloNumArray d_vec; // durée de trajet des arcs
    IloNumArray D_vec; // incertitude durée de trajet des arcs
    std::vector<Arc> mat;

    std::vector<IloInt> sol; // liste des villes visitées dans l'ordre

    Instance(IloEnv env, char filename[]);
    ~Instance(){};
    void display() const;

    double compute_static_score() const { return compute_static_score(sol); }
    double compute_robust_score(IloEnv env, const unsigned int& time_limit=60, const int& verbose=0) const { return compute_robust_score(env, sol, time_limit, verbose); }
    double compute_static_constraint() const { return compute_static_constraint(sol); }
    double compute_robust_constraint(IloEnv env, const unsigned int& time_limit=60, const int& verbose=0) const { return compute_robust_constraint(env, sol, time_limit, verbose); }

    double compute_static_score(const std::vector<IloInt>& sol) const;
    double compute_robust_score(IloEnv env, const std::vector<IloInt>& sol, const unsigned int& time_limit =60, const int& verbose=0) const;
    double compute_static_constraint(const std::vector<IloInt>& sol) const;
    double compute_robust_constraint(IloEnv env, const std::vector<IloInt>& sol, const unsigned int& time_limit=60, const int& verbose=0) const;

    void exportSol(char filename[]) const;
};


#endif
