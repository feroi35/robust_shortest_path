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
    SolveMethod() {};

    std::string method_name;
    double infBound = 0.0;
    int nodesExplored = 0;
    int nCallBacks = 0;
    double callBacksTimeSpan = 0.0;

    void add_static_constraints(IloEnv& env, IloModel& model, IloBoolVarArray& x, IloBoolVarArray& y, IloNumVar z, const Instance& inst) const;
    void retrieveCplexSolution(const IloCplex& cplex, const IloNumArray& xValues, Instance& inst);
    void parametrizeCplex(IloCplex& cplex, const unsigned int& time_limit, const int& verbose) const;

    virtual void solve(IloEnv& env, Instance& inst, const unsigned int& time_limit, const int& verbose) = 0;
    
    void solve_and_display(IloEnv& env, Instance& inst, const unsigned int& time_limit, const int& verbose);
};

#endif