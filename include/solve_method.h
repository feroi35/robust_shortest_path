// solve_method.h
#ifndef SOLVE_METHOD_H
#define SOLVE_METHOD_H

#include <chrono>
#include <vector>
#include <string>
class IloEnv;
class IloCplex;
class IloExpr;
class IloModel;
class IloNumVar;
class IloBoolVarArray;
class IloNumArray;
class Instance;

namespace Subproblems {
    double solve_objective_subproblem(IloEnv& env, const IloNumArray& xValues, const IloBoolVarArray& x, const Instance& inst, IloExpr& expr);
    double solve_constraint_subproblem(IloEnv& env, const IloNumArray& yValues, const IloBoolVarArray& y, const Instance& inst, IloExpr& expr);
}

struct SolveMethod {
    std::string method_name;
    int nodesExplored = 0;
    int nCallBacks = 0; // Branch and cut and plans coupants
    double callBacksTimeSpan = 0.0; // Temps passé à résoudre des sous-problèmes

    void add_static_constraints(IloEnv& env, IloModel& model, IloBoolVarArray& x, IloBoolVarArray& y, IloNumVar z, Instance& inst);
    void retrieveCplexSolution(const IloCplex& cplex, const IloNumArray& xValues, Instance& inst);
    void parametrizeCplex(IloCplex& cplex, const unsigned int& time_limit, const int& verbose);

    virtual void solve(IloEnv& env, Instance& inst, const unsigned int& time_limit, const int& verbose) = 0;
    
    void solve_and_display(IloEnv& env, Instance& inst, const unsigned int& time_limit, const int& verbose);
};

#endif