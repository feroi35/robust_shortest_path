#include <fstream>
#include "solve_method.h"
#include "parser.h"


double Subproblems::solve_objective_subproblem(IloEnv& env, const IloNumArray& xValues,
            const IloBoolVarArray& x, const Instance& inst, IloExpr& expr) {
    // Compute the robust objective (solving a continuous knapsack) and modify the expression
    // It is done in the same function because it allows to be more efficient
    std::vector<IloNum> uncertainties; // D
    std::vector<IloNum> weights; // d
    std::vector<IloInt> idx_edges;
    double static_score = 0.0;
    unsigned int n_edges = 0;

    for (unsigned int a=0; a<inst.n_arc; a++) {
        if (xValues[a] > 1 - TOL) {
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
    return robust_score;
}


double Subproblems::solve_constraint_subproblem(IloEnv& env, const IloNumArray& yValues,
        const IloBoolVarArray& y, const Instance& inst, IloExpr& expr) {
    // Compute the robust constraint (solving a continuous knapsack) and modify the expression
    // It is done in the same function because it allows to be more efficient
    std::vector<IloNum> weight; // ph
    std::vector<IloInt> idx_nodes;
    double static_constraint = 0.0;
    unsigned int n_nodes = 0;

    for (unsigned int i=0; i<inst.n; i++) {
        if (yValues[i] > TOL) {
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
    return robust_constraint;
}


std::string get_path(std::vector<IloInt> solution, const Instance& inst){
    std::string path_str = "[";
    unsigned int current_node = inst.sol[0];
    unsigned int next_node;
    if (current_node != inst.s) {
        throw std::domain_error("First node of solution is not s for instance " + inst.name);
    }
    path_str += std::to_string(current_node) + ";" ;
    for (unsigned int i=1; i<inst.sol.size(); i++) {
        next_node = inst.sol[i];
        if (inst.d[current_node-1][next_node-1] !=  inst.d[current_node-1][next_node-1]) {
            throw std::domain_error("Using arc that does not exist for instance " + inst.name);
        }
        current_node = next_node;
        path_str += std::to_string(current_node);
        path_str += (i == inst.sol.size()-1) ? "]" : ";";
    }
    if (current_node != inst.t) {
        throw std::domain_error("Last node of solution is not t for instance " + inst.name);
    }
    return path_str;
}


void SolveMethod::add_static_constraints(IloEnv& env, IloModel& model, IloBoolVarArray& x, IloBoolVarArray& y, IloNumVar z, const Instance& inst) const {
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
                // you can't get out of t
                model.add(out_arcs_i == 0);
            }
        }
        if (!empty_in_arcs_i) {
            if (i != inst.s-1) {
                model.add(in_arcs_i == y[i]);
            } else {
                // you can't get into s
                model.add(in_arcs_i == 0);
            }
        }
        out_arcs_i.end();
        in_arcs_i.end();
    }
    model.add(y[inst.s-1] == 1);
    model.add(y[inst.t-1] == 1);
    model.add(IloScalProd(y, inst.p) <= inst.S);
    model.add(z >= IloScalProd(x, inst.d_vec));
}


void SolveMethod::retrieveCplexSolution(const IloCplex& cplex, const IloNumArray& xValues, Instance& inst) {
    // Makes some checks and updates inst.sol, nodesExplored
    if (cplex.getStatus() == IloAlgorithm::Infeasible) {
        throw std::domain_error("Infeasible " + method_name + " model for instance " + inst.name);
    } else if (cplex.getStatus() == IloAlgorithm::Unknown) {
        throw std::domain_error("No solution found for method " + method_name + " with instance " + inst.name + ". Maybe not enough time");
    }

    if (!inst.sol.empty()) {
        std::cerr << "Warning: solution vector not empty for instance " << inst.name << std::endl;
        inst.sol.clear();
    }
    unsigned int current_node = inst.s-1;
    while (current_node != inst.t-1) {
        inst.sol.push_back(current_node+1);
        for (unsigned int a=0; a<inst.n_arc; a++) {
            if (inst.mat[a].tail == current_node+1 && xValues[a] == 1) {
                current_node = inst.mat[a].head-1;
                break;
            }
        }
        if (current_node == inst.sol[inst.sol.size()-1]-1) {
            throw std::domain_error("Using arc that does not exist for instance " + inst.name);
        }
    }
    inst.sol.push_back(inst.t);
}


void SolveMethod::parametrizeCplex(IloCplex& cplex, const unsigned int& time_limit, const int& verbose) const {
    cplex.setParam(IloCplex::Param::TimeLimit, time_limit);
    if (verbose < 2) cplex.setOut(cplex.getEnv().getNullStream());
}


void SolveMethod::solve_and_display(IloEnv& env, Instance& inst, const unsigned int& time_limit, const int& verbose) {
    try {
        std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
        solve(env, inst, time_limit, verbose);
        std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
        std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::microseconds>(end - start);

        double static_cost = inst.compute_static_score();
        double static_constraint = inst.compute_static_constraint();
        double robust_constraint = inst.compute_robust_constraint_milp(env);
        bool admissibility = robust_constraint < inst.S + TOL;
        double robust_score = admissibility ? inst.compute_robust_score_milp(env) : 1e9;
        std::string path_str = get_path(inst.sol, inst);

        // Display results
        std::cout << inst.name << ","
            << method_name << ","
            << robust_score << ","
            << infBound << ","
            << static_cast<double>(time_span.count()) << ","
            << nodesExplored << ","
            << robust_constraint << ","
            << static_cost << ","
            << static_constraint << ","
            << inst.S << ","
            << path_str << ","
            << admissibility << ","
            << nCallBacks << ","
            << callBacksTimeSpan << std::endl;
    }
    catch (IloException& e) {
        std::cout << inst.name << "," << method_name << "," << ",,,,,,,," << std::endl;
        throw e;
    }
    catch (std::domain_error& e) {
        std::cout << inst.name << "," << method_name << "," << ",,,,,,,," << std::endl;
        throw e;
    }
    catch (std::exception& e) {
        std::cout << inst.name << "," << method_name << "," << ",,,,,,,," << std::endl;
        throw e;
    }
}
