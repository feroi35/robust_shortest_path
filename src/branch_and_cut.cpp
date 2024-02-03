#include "branch_and_cut.h"
#include "parser.h"


ILOLAZYCONSTRAINTCALLBACK7(myCallBack, const IloBoolVarArray&, x, const IloBoolVarArray&, y,
        const IloNumVar&, z, const Instance&, inst, const unsigned int&, verbose, int&, nIteration,
        double&, spentTime) {
    nIteration++;
    if (verbose > 0) std::cout << "Lazy constraint callback, it " << nIteration << std::endl;
    std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();

    IloEnv env = getEnv();
    IloNumArray xValues(env);
    IloNumArray yValues(env);
    getValues(xValues, x);
    getValues(yValues, y);
    double zValue = getValue(z);

    IloExpr exprObjective(env);
    double robust_score = Subproblems::solve_objective_subproblem(env, xValues, x, inst, exprObjective);
    if (robust_score > zValue + TOL) {
        add(exprObjective <= z);
    }
    exprObjective.end();

    IloExpr exprConstraint(env);
    double robust_constraint = Subproblems::solve_constraint_subproblem(env, yValues, y, inst, exprConstraint);
    if (robust_constraint > inst.S + TOL) {
        add(exprConstraint <= inst.S);
    }
    exprConstraint.end();

    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::chrono::microseconds duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    double time_spent = static_cast<double>(duration.count()) / 1e6;
    if (verbose > 0) std::cout << "time: " << time_spent << std::endl;
    spentTime += time_spent;

    xValues.end();
    yValues.end();
}


void BranchAndCutMethod::solve(IloEnv& env, Instance& inst, const unsigned int& time_limit, const int& verbose) {
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

    IloCplex cplex(model);
    parametrizeCplex(cplex, time_limit, verbose);

    cplex.use(myCallBack(env,x,y,z,inst,verbose,nCallBacks,callBacksTimeSpan));
    cplex.solve();

    // Retrieve solution
    IloNumArray xValues(env);
    cplex.getValues(xValues, x);
    retrieveCplexSolution(cplex, xValues, inst);
    xValues.end();
}
