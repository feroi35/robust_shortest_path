#include "dualized_formulation.h"
#include "parser.h"


void DualizedMethod::solve(IloEnv& env, Instance& inst, const unsigned int& time_limit, const int& verbose) {
    IloModel model(env);

    // Variables
    IloBoolVarArray x(env, inst.n_arc);
    IloNumVarArray lambda(env, inst.n_arc, 0.0, IloInfinity);
    IloBoolVarArray y(env, inst.n);
    IloNumVarArray beta(env, inst.n, 0.0, IloInfinity);
    IloNumVar alpha(env, 0.0, IloInfinity);
    IloNumVar eta(env, 0.0, IloInfinity);
    IloNumVar z(env, 0.0, IloInfinity);

    // Objective
    IloObjective obj(env, z, IloObjective::Minimize);
    model.add(obj);

    // Constraints
    add_static_constraints(env, model, x, y, z, inst);
    model.add(z >= inst.d1*eta + IloScalProd(x, inst.d_vec) + IloScalProd(lambda, inst.D_vec));
    for (unsigned int a = 0; a < inst.n_arc; ++a) {
        model.add(eta + lambda[a] >= inst.d_vec[a]*x[a]);
    }
    model.add(inst.d2*alpha + 2*IloSum(beta) + IloScalProd(y, inst.p) <= inst.S);
    for (unsigned int i=0; i<inst.n; i++) {
        model.add(alpha + beta[i] >= inst.ph[i]*y[i]);
    }

    // Solve
    IloCplex cplex(model);
    parametrizeCplex(cplex, time_limit, verbose);

    cplex.solve();

    // Retrieve solution
    IloNumArray xValues(env);
    cplex.getValues(xValues, x);
    retrieveCplexSolution(cplex, xValues, inst);
    xValues.end();

    nodesExplored = cplex.getNnodes();
    infBound = cplex.getBestObjValue();
}
