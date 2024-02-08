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
    for (unsigned int a = 0; a < inst.n_arc; ++a) {
        best_xValues.add(0.0);
    }
    for (unsigned int i = 0; i < inst.n; ++i) {
        best_yValues.add(0.0);
    }
    double best_score = 1e9;
    double less_violated_constraint = 1e9;
    double best_violated_score = 1e9;
    // Even if the solution is not admissible, it is still useful to return the less violated solution found
    bool admissible_solution_found = false;
    double best_bound = 0.0;
    IloNumArray xValues(env);
    IloNumArray yValues(env);

    std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
    double spentTime = 0.0;

    IloCplex cplex;
    cplex = IloCplex(model);
    parametrizeCplex(cplex, time_limit, verbose);
    // The IloCplex object remains the same, but the model is different
    // It allows to keep warm starts.

    bool optimality = false;
    bool new_best_sol_found = false;
    unsigned int nWarmStarts = 0;
    while (!optimality && spentTime < time_limit) {
        if (verbose > 0) {
            std::cout << "Plans coupants: iteration " << nCallBacks << std::endl;
            std::cout << "Time: " << spentTime << std::endl;
        }

        if (warmStart && new_best_sol_found) {
            nWarmStarts++;
            if (verbose > 0) {
                std::cout << "Manual added warm start " << nWarmStarts << std::endl;
                std::cout << "Real warm starts: " << cplex.getNMIPStarts() << std::endl;
            }
            IloNumVarArray allStartVar(env);
            IloNumArray allStartVal(env);
            for (unsigned int a = 0; a < inst.n_arc; ++a) {
                allStartVar.add(x[a]);
                allStartVal.add(best_xValues[a]);
            }
            for (unsigned int i = 0; i < inst.n; ++i) {
                allStartVar.add(y[i]);
                allStartVal.add(best_yValues[i]);
            }
            allStartVar.add(z);
            allStartVal.add(best_score);
            cplex.addMIPStart(allStartVar, allStartVal);
            allStartVar.end();
            allStartVal.end();
            new_best_sol_found = false;
        }

        if (verbose > 1) std::cout << "\nSolving..." << std::endl;
        cplex.setParam(IloCplex::Param::TimeLimit, (int) std::max(120., (double) time_limit-spentTime));
        cplex.solve();
        if (verbose > 1) std::cout << "Solved!" << std::endl;
        nodesExplored += cplex.getNnodes();
        // cplexCheckStatus(cplex, inst);

        if (cplex.getStatus() == IloAlgorithm::Unknown) {
            break; 
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

        // Breaking symmetries
        std::vector<IloInt> sol;
        unsigned int current_node = inst.s-1;
        while (current_node != inst.t-1) {
            sol.push_back(current_node+1);
            if (yValues[current_node] < 1 - TOL) {
                std::cerr << "No equivalence between xValues and yValues" << std::endl;
                throw std::domain_error("A node is not reached");
            }
            for (unsigned int a=0; a<inst.n_arc; a++) {
                if (inst.mat[a].tail == current_node+1 && xValues[a] >= 0.5) {
                    current_node = inst.mat[a].head-1;
                    break;
                }
            }
            if (current_node == sol[sol.size()-1]-1) {
                throw std::domain_error("Using arc that does not exist for instance " + inst.name);
            }
        }
        sol.push_back(inst.t);
        if (yValues[inst.t-1] < 1 - TOL) {
            std::cerr << "t is not reached" << std::endl;
            throw std::domain_error("You can't get into t");
        }

        for (unsigned int i=0; i<sol.size()-2;i++){
            if (inst.pair_nodes[sol[i]-1][sol[i+2]-1]) {
                continue; // Already explored possible subpath
            }
            inst.pair_nodes[sol[i]-1][sol[i+2]-1] = true;
            int node_i = (int) sol[i]-1;
            int node_j = (int) sol[i+2]-1;
            std::vector<std::vector<int>> to_forbid = arcs_to_forbid(inst, node_i, node_j);
            for (unsigned int k=0; k<to_forbid.size(); k++) {
                if (to_forbid[k][2] == to_forbid[k][3]) {
                    std::cerr << "Arc to forbid loop on arc " << to_forbid[k][3] << " for node " << to_forbid[k][0] << "and node " << to_forbid[k][1] << endl;
                }
                model.add(x[to_forbid[k][2]] + x[to_forbid[k][3]] <= 1);
            }
        }

        callBacksTimeSpan += static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - startSubProblemsResolution).count()) / 1e6;

        // Update best solution
        if ((!violated_constraint && robust_objective < best_score)
                || (!admissible_solution_found && robust_constraint < less_violated_constraint)
                || (!admissible_solution_found && robust_constraint == less_violated_constraint && robust_objective < best_violated_score)) {
            admissible_solution_found = admissible_solution_found || !violated_constraint;
            if (admissible_solution_found) {
                best_score = robust_objective;
            } else {
                best_violated_score = robust_objective;
                less_violated_constraint = robust_constraint;
            }

            for (unsigned int a = 0; a < inst.n_arc; ++a) {
                best_xValues[a] = xValues[a];
            }
            for (unsigned int i = 0; i < inst.n; ++i) {
                best_yValues[i] = yValues[i];
            }
            new_best_sol_found = true;
        }
        optimality = !violated_objective && !violated_constraint;
        nCallBacks++;
        spentTime = static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - start).count()) / 1e6;
    }

    // Retrieve solution
    retrieveCplexSolution(best_xValues, inst);
    best_xValues.end();
    best_yValues.end();
    xValues.end();
    yValues.end();
    infBound = best_bound;
}
