#include <chrono>
#include "plans_coupants.h"
#include "parser.h"


double solve_objective_subproblem(IloEnv env, IloCplex cplex, IloModel model, const IloNumArray& xValues,
        const IloBoolVarArray& x, const IloNumVar& z, const double current_best_value, const Instance& inst) {
    // Compute the robust objective (solving a continuous knapsack)
    // And add the violated constraint (if any)
    // It is done in the same function because it allows to be more efficient
    std::vector<IloNum> uncertainties; // D
    std::vector<IloNum> weights; // d
    std::vector<IloInt> idx_edges;
    double static_score = 0.0;
    unsigned int n_edges = 0;

    for (unsigned int a=0; a<inst.n_arc; a++) {
        if (xValues[a] > 1e-3) {
            uncertainties.push_back(inst.mat[a].D);
            weights.push_back(inst.mat[a].d);
            idx_edges.push_back(a);
            static_score += inst.mat[a].d;
            n_edges++;
        }
    }
    std::vector<size_t> argsorted_weights = argsort(weights);
    double robust_attack = 0.0;
    double used_budget = 0.0;
    int idx = n_edges-1;
    IloExpr expr(env);
    while (used_budget < inst.d1 && idx >= 0) {
        unsigned int arc_idx = argsorted_weights[idx];
        float delta1_i = std::min(inst.d1 - used_budget, uncertainties[arc_idx]);
        used_budget += delta1_i;
        robust_attack += delta1_i * weights[arc_idx];

        IloInt real_arc = idx_edges[arc_idx];
        expr += inst.mat[real_arc].d * delta1_i * x[real_arc];
        idx--;
    }
    expr += IloScalProd(x, inst.d_vec);
    double robust_score = static_score + robust_attack;

    // Add the violated constraint
    if (robust_score > current_best_value + 1e-3) {
        model.add(expr <= z);
    }
    expr.end();
    return robust_score;
}


double solve_constraint_subproblem(IloEnv env, IloCplex cplex, IloModel model, const IloNumArray& yValues,
        const IloBoolVarArray& y, const Instance& inst) {
    // Compute the robust constraint (solving a continuous knapsack)
    // and add the violated constraint (if any)
    // It is done in the same function because it allows to be more efficient
    std::vector<IloNum> weight; // ph
    std::vector<IloInt> idx_nodes;
    double static_constraint = 0.0;
    unsigned int n_nodes = 0;

    for (unsigned int i=0; i<inst.n; i++) {
        if (yValues[i] > 1e-3) {
            weight.push_back(inst.ph[i]);
            idx_nodes.push_back(i);
            static_constraint += inst.p[i];
            n_nodes++;
        }
    }

    std::vector<size_t> argsorted_weight = argsort(weight);
    double robust_attack2 = 0.0;
    double used_budget2 = 0.0;
    int idx2 = n_nodes-1;
    IloExpr expr(env);
    while (used_budget2 < inst.d2 && idx2 >= 0) {
        unsigned int node_idx2 = argsorted_weight[idx2];
        float delta2_i = std::min(inst.d2 - used_budget2, weight[node_idx2]);
        used_budget2 += delta2_i;
        robust_attack2 += delta2_i * weight[node_idx2];

        IloInt real_node = idx_nodes[node_idx2];
        expr += inst.ph[real_node] * delta2_i * y[real_node];
        idx2--;
    }
    expr += IloScalProd(y, inst.p);
    double robust_constraint = static_constraint + robust_attack2;
    // Add the violated constraint
    if (robust_constraint > inst.S + 1e-3) {
        model.add(expr <= inst.S);
    }
    expr.end();
    return robust_constraint;
}


void plans_coupants_solve(IloEnv env, Instance& inst, const unsigned int& time_limit, const int& verbose) {
    IloModel model(env);

    ///////////////////////////
    // Static model
    // Variables
    IloBoolVarArray x(env, inst.n_arc);
    IloBoolVarArray y(env, inst.n);
    IloNumVar z(env, 0.0, IloInfinity, "z");

    // Objective
    IloObjective obj(env, z, IloObjective::Minimize);
    model.add(obj);

    // Constraints
    // Flow conservation
    for (unsigned int i=0; i<inst.n; i++) {
        IloExpr out_arcs_i(env);
        IloExpr in_arcs_i(env);
        bool empty_out_arcs_i = true;
        bool empty_in_arcs_i = true;
        for (unsigned int a=0; a<inst.n_arc; a++) {
            if (inst.mat[a].tail == i+1) {
                out_arcs_i += x[a];
                empty_out_arcs_i = false;
            }
            if (inst.mat[a].head == i+1) {
                in_arcs_i += x[a];
                empty_in_arcs_i = false;
            }
        }
        if (!empty_out_arcs_i) {
            if (i != inst.t-1) {
                model.add(out_arcs_i == y[i]);
            } else {
                model.add(out_arcs_i == 0);
                // you can't get out of t
            }
        }
        if (!empty_in_arcs_i) {
            if (i != inst.s-1) {
                model.add(in_arcs_i == y[i]);
            } else {
                model.add(in_arcs_i == 0);
                // you can't get into s
            }
        }
        out_arcs_i.end();
        in_arcs_i.end();
    }
    model.add(y[inst.s-1] == 1);
    model.add(y[inst.t-1] == 1);
    model.add(IloScalProd(y, inst.p) <= inst.S);
    model.add(z >= IloScalProd(x, inst.d_vec));
    ///////////////////////////

    ///////////////////////////
    // Stock the current best solution
    // Indeed, if at some point, an admissible solution is found but its score does not match, it will be excluded
    // And the final solution could be not admissible (because of the time limit)
    IloNumArray best_xValues(env);
    IloNumArray best_yValues(env);
    double best_score = 1e6;
    double best_bound = 0.0;
    ///////////////////////////
    IloNumArray xValues(env); // to retrieve the values
    IloNumArray yValues(env);

    std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
    float time_spent = 0.0;

    IloCplex cplex;
    cplex = IloCplex(model); // Je repars du même modèle à chaque fois auquel j'ajoute des contraintes au fur et à mesure
    cplex.setParam(IloCplex::Param::TimeLimit, time_limit);

    bool optimality = false;
    bool new_best_sol_found = false;
    unsigned int iteration = 0;
    while (!optimality && time_spent < time_limit) {
        if (verbose > 0) {
            std::cout << "Plans coupants: iteration " << iteration << std::endl;
            std::cout << "Time spent: " << time_spent << std::endl;
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
        if (verbose < 2) cplex.setOut(env.getNullStream());
        cplex.solve();

        if (cplex.getStatus() == IloAlgorithm::Infeasible) {
            std::cout << inst.name << "," << "plans_coupants,,,,,,,,," << std::endl;
            throw std::domain_error("Infeasible plans_coupants model for instance " + inst.name + " at iteration " + std::to_string(iteration));
        } else if (cplex.getStatus() == IloAlgorithm::Unknown) {
            std::cout << inst.name << "," << "plans_coupants,,,,,,,,," << std::endl;
            throw std::domain_error("No solution found for instance " + inst.name  + " at iteration " + std::to_string(iteration) + ". Maybe not enough time");
        }

        // Retrieve the values. Done here because y cant be retrieve after a new constraint has been added
        cplex.getValues(xValues, x); // cplex.getValues does not work for IloBoolArray...
        cplex.getValues(yValues, y);
        double current_best_value = cplex.getObjValue();
        double current_best_bound = cplex.getBestObjValue();
        if (best_bound < current_best_bound) {
            best_bound = current_best_bound;
            model.add(z >= best_bound); // Add the bound
        }

        double robust_objective = solve_objective_subproblem(env, cplex, model, xValues, x, z, current_best_value, inst);
        double robust_constraint = solve_constraint_subproblem(env, cplex, model, yValues, y, inst);
        bool violated_objective = robust_objective > current_best_value + 1e-3;
        bool violated_constraint = robust_constraint > inst.S + 1e-3;

        if (!violated_constraint && robust_objective < best_score) {
            best_score = robust_objective;
            best_xValues = xValues;
            best_yValues = yValues;
            new_best_sol_found = true;
        }
        bool violated = violated_objective || violated_constraint;
        optimality = !violated;
        iteration++;
        time_spent = static_cast<float>(std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - start).count()) / 1e6;
    }
    if (!inst.sol.empty()) {
        std::cerr << "Warning: solution vector not empty for instance " << inst.name << std::endl;
        inst.sol.clear();
    }
    std::string path_str = "[";
    unsigned int current_node = inst.s-1;
    while (current_node != inst.t-1) {
        path_str += std::to_string(current_node+1) + ";";
        inst.sol.push_back(current_node+1);
        for (unsigned int a=0; a<inst.n_arc; a++) {
            if (inst.mat[a].tail == current_node+1 && best_xValues[a] == 1) {
                // Différence par rapport aux autres méthodes, best_xValues au lieu de model.getValue(x[a])
                current_node = inst.mat[a].head-1;
                break;
            }
        }
        if (current_node == inst.sol[inst.sol.size()-1]-1) {
            std::cout << inst.name << "," << "branch_and_cut,,,,,,,,," << std::endl;
            throw std::domain_error("Using arc that does not exist for instance " + inst.name);
        }
    }
    inst.sol.push_back(inst.t);
    path_str += std::to_string(inst.t) + "]";

    float robust_constraint = inst.compute_robust_constraint(env);
    float robust_score = (robust_constraint < inst.S + 1e-3) ? inst.compute_robust_score(env) : 1e6;
    std::cout << inst.name << ","
        << "plans_coupants,"
        << robust_score << ","
        << best_bound << ","
        << time_spent << ","
        << cplex.getNnodes() << ","
        << robust_constraint << ","
        << inst.compute_static_score() << ","
        << inst.compute_static_constraint() << ","
        << inst.S << ","
        << path_str << std::endl;
}
