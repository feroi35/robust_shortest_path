#include "static_solve.h"
#include "parser.h"


void StaticMethod::solve(IloEnv& env, Instance& inst, const unsigned int& time_limit, const int& verbose) {
    IloModel model(env);

    // Variables
    IloBoolVarArray x(env, inst.n_arc);
    IloBoolVarArray y(env, inst.n);
    IloNumVar z(env, 0.0, IloInfinity, "z");

    // Objective
    IloObjective obj(env, z, IloObjective::Minimize);
    model.add(obj);

    // Constraints
    add_static_constraints(env, model, x, y, z, inst);

    // Solve
    IloCplex cplex(model);
    parametrizeCplex(cplex, time_limit, verbose);
    cplex.solve();

    if (cplex.getStatus() == IloAlgorithm::Infeasible) {
        throw std::domain_error("Infeasible " + method_name + " model for instance " + inst.name);
    } else if (cplex.getStatus() == IloAlgorithm::Unknown) {
        throw std::domain_error("No solution found (yet) for method " + method_name + " with instance " + inst.name + ". Maybe not enough time");
    }
    // Retrieve solution
    IloNumArray xValues(env);
    cplex.getValues(xValues, x);
    retrieveCplexSolution(cplex, xValues, inst);
    xValues.end();

    nodesExplored = cplex.getNnodes();
    infBound = cplex.getBestObjValue();

    if (abs(inst.compute_static_score() - cplex.getObjValue()) > TOL) {
        throw std::domain_error("Not the same static objective value for instance " + inst.name);
    }
}
