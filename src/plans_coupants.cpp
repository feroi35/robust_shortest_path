#include <chrono>
#include "plans_coupants.h"
#include "parser.h"


bool solve_subproblems(IloEnv env, IloCplex cplex, IloModel model, IloBoolVarArray x, IloBoolVarArray y, IloNumVar z,
        Instance inst, unsigned int verbose) {
    std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
    bool violated = false;
    // Retrieve the values. Done here because y cant be retrieve after a new constraint has been added
    // But it is possible to do so in CallBack...
    IloNumArray xValues(env);
    cplex.getValues(xValues, x); // Fait gagner un facteur 10,000 en temps par rapport à getValue
    IloNumArray yValues(env);
    cplex.getValues(yValues, y); // Fait gagner un facteur 10,000 en temps par rapport à getValue

    ////////////////////////
    // Compute the robust objective (solving a continuous knapsack)
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
    if (robust_score > cplex.getValue(z) + 1e-3) {
        model.add(expr <= z);
        violated = true;
    }
    expr.end();
    ////////////////////////

    ////////////////////////
    // Compute the robust constraint
    std::vector<IloNum> weights2; // ph
    std::vector<IloInt> idx_nodes;
    double static_constraint = 0.0;
    unsigned int n_nodes = 0;

    for (unsigned int i=0; i<inst.n; i++) {
        if (yValues[i] > 1e-3) {
            weights2.push_back(inst.ph[i]);
            idx_nodes.push_back(i);
            static_constraint += inst.p[i];
            n_nodes++;
        }
    }

    std::vector<size_t> argsorted_weights2 = argsort(weights2);
    double robust_attack2 = 0.0;
    double used_budget2 = 0.0;
    int idx2 = n_nodes-1;
    IloExpr expr2(env);
    while (used_budget2 < inst.d2 && idx2 >= 0) {
        unsigned int node_idx2 = argsorted_weights2[idx2];
        float delta2_i = std::min(inst.d2 - used_budget2, weights2[node_idx2]);
        used_budget2 += delta2_i;
        robust_attack2 += delta2_i * weights2[node_idx2];

        IloInt real_node = idx_nodes[node_idx2];
        expr2 += inst.ph[real_node] * delta2_i * y[real_node];
        idx2--;
    }
    expr2 += IloScalProd(y, inst.p);
    double robust_constraint = static_constraint + robust_attack2;
    // Add the violated constraint
    if (robust_constraint > inst.S + 1e-3) {
        model.add(expr2 <= inst.S);
        violated = true;
    }
    expr2.end();
    ////////////////////////

    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::chrono::microseconds duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    if (verbose > 0) std::cout << "Resolution sous problèmes: " << static_cast<double>(duration.count()) / 1e6 << std::endl;
    return violated;
}


void plans_coupants_solve(IloEnv env, Instance& inst, const unsigned int& time_limit, const int& verbose) {
    IloModel model(env);

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

    std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
    float time_spent = 0.0;

    IloCplex cplex;
    bool optimality = false;
    unsigned int count = 0;
    while (!optimality) {
        if (verbose > 0) {
            std::cout << "Plans coupants: iteration " << std::to_string(count) << std::endl;
            std::cout << "Time spent: " << time_spent << std::endl; 
        }

        cplex = IloCplex(model); // Je repars du même modèle à chaque fois auquel j'ajoute des contraintes au fur et à mesure
        cplex.setParam(IloCplex::Param::TimeLimit, time_limit);
        if (verbose < 2) cplex.setOut(env.getNullStream());
        cplex.solve();
        count++;

        if (cplex.getStatus() == IloAlgorithm::Infeasible) {
            std::cout << inst.name << "," << "plans_coupants,,,,,,,,," << std::endl;
            throw std::domain_error("Infeasible plans_coupants model for instance " + inst.name + " at iteration " + std::to_string(count));
        } else if (cplex.getStatus() == IloAlgorithm::Unknown) {
            std::cout << inst.name << "," << "plans_coupants,,,,,,,,," << std::endl;
            throw std::domain_error("No solution found for instance " + inst.name  + " at iteration " + std::to_string(count) + ". Maybe not enough time");
        }

        std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
        time_spent = static_cast<float>(std::chrono::duration_cast<std::chrono::microseconds>(end - start).count()) / 1e6;
        if (time_spent > time_limit) {
            // si on ne break pas ici et qu'on ajoute des contraintes, on peut ne peut plus accéder aux valeurs des variables
            // car on peut avoir un problème infaisable
            // Remarque: on peut perdre une solution réalisable contre une qui a un meilleur score mais n'est pas réalisable...
            break;
        }

        bool violated = false;
        violated = solve_subproblems(env, cplex, model, x, y, z, inst, verbose);
        optimality = !violated;
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
            if (inst.mat[a].tail == current_node+1 && cplex.getValue(x[a]) == 1) {
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
    float robust_score = (robust_constraint < inst.S + 1e-3) ? inst.compute_robust_score(env) : 1000000;
    std::cout << inst.name << ","
        << "plans_coupants,"
        << robust_score << ","
        << cplex.getBestObjValue() << ","
        << time_spent << ","
        << cplex.getNnodes() << ","
        << robust_constraint << ","
        << inst.compute_static_score() << ","
        << inst.compute_static_constraint() << ","
        << inst.S << ","
        << path_str << std::endl;
}
