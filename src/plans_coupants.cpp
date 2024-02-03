#include "plans_coupants.h"
#include "parser.h"


void PlansCoupantsMethod::solve(IloEnv& env, Instance& inst, const unsigned int& time_limit, const int& verbose) {
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

    // Stock the current best solution
    // If at some point, an admissible solution is found but its score does not match, it will be excluded
    // And the final solution could be not admissible (because of the time limit)
    IloNumArray best_xValues(env);
    IloNumArray best_yValues(env);
    double best_score = 1e9;
    double best_bound = 0.0;
    IloNumArray xValues(env); // to retrieve the values
    IloNumArray yValues(env);

    std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
    double spentTime = 0.0;

    IloCplex cplex;
    cplex = IloCplex(model);
    parametrizeCplex(cplex, time_limit, verbose);

    bool optimality = false;
    bool new_best_sol_found = false;
    while (!optimality && spentTime < time_limit) {
        if (verbose > 0) {
            std::cout << "Plans coupants: iteration " << nCallBacks << std::endl;
            std::cout << "Time: " << spentTime << std::endl;
        }

        // Warm start
        if (new_best_sol_found) {
            if (verbose > 0) std::cout << "Warm start" << std::endl;
            IloNumVarArray startVarX(env);
            for (unsigned int a = 0; a < inst.n_arc; ++a) {
                startVarX.add(x[a]);
            }
            cplex.addMIPStart(startVarX, best_xValues);
            startVarX.end();

            IloNumVarArray startVarY(env);
            for (unsigned int i = 0; i < inst.n; ++i) {
                startVarY.add(y[i]);
            }
            cplex.addMIPStart(startVarY, best_yValues);
            startVarY.end();

            IloNumVarArray startVarZ(env);
            IloNumArray zArray(env, 1, best_score);
            startVarZ.add(z);
            zArray[0] = best_score;
            cplex.addMIPStart(startVarZ, zArray);
            startVarZ.end();
            zArray.end();
        }
        new_best_sol_found = false;

        cplex.solve();

        if (cplex.getStatus() == IloAlgorithm::Infeasible) {
            throw std::domain_error("Infeasible plans_coupants model for instance " + inst.name + " at iteration " + std::to_string(nCallBacks));
        } else if (cplex.getStatus() == IloAlgorithm::Unknown) {
            throw std::domain_error("No solution found for instance " + inst.name  + " at iteration " + std::to_string(nCallBacks) + ". Maybe not enough time");
        }

        cplex.getValues(xValues, x); // cplex.getValues does not work for IloBoolArray...
        cplex.getValues(yValues, y);
        double current_best_value = cplex.getObjValue();
        double current_best_bound = cplex.getBestObjValue();

        // Update best bound
        if (best_bound < current_best_bound) {
            best_bound = current_best_bound;
            model.add(z >= best_bound);
        }

        // Solve subproblems
        std::chrono::steady_clock::time_point startSubProblemsResolution = std::chrono::steady_clock::now();
        IloExpr exprObjective(env);
        double robust_objective = Subproblems::solve_objective_subproblem(env, xValues, x, inst, exprObjective);
        bool violated_objective = robust_objective > current_best_value + TOL;
        if (violated_objective) {
            model.add(exprObjective <= z);
        }
        exprObjective.end();

        IloExpr exprConstraint(env);
        double robust_constraint = Subproblems::solve_constraint_subproblem(env, yValues, y, inst, exprConstraint);
        bool violated_constraint = robust_constraint > inst.S + TOL;
        if (violated_constraint) {
            model.add(exprConstraint <= inst.S);
        }
        exprConstraint.end();
        callBacksTimeSpan += static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - startSubProblemsResolution).count()) / 1e6;

        // Update best solution
        if (!violated_constraint && robust_objective < best_score) {
            best_score = robust_objective;
            best_xValues.end();
            best_yValues.end();
            best_xValues = xValues;
            best_yValues = yValues;
            new_best_sol_found = true;
        }
        optimality = !violated_objective && !violated_constraint;
        nCallBacks++;
        nodesExplored += cplex.getNnodes();
        spentTime = static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - start).count()) / 1e6;
    }

    // Retrieve solution
    retrieveCplexSolution(cplex, best_xValues, inst);
    best_xValues.end();
    best_yValues.end();
    if (xValues.getSize() > 0) {
        // It the best solution is obtained at the last iteration, there can be segfault
        xValues.end();
        yValues.end();
    }
    infBound = cplex.getBestObjValue();
}
