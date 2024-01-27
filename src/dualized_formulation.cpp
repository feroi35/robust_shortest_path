#include "dualized_formulation.h"
#include "parser.h"


void dualized_solve(IloEnv env, Instance& inst, const unsigned int& time_limit, const int& verbose) {
    IloModel model(env);

    // Variables
    IloBoolVarArray x(env, inst.n_arc);
    IloNumVarArray lambda(env, inst.n_arc, 0.0, IloInfinity);
    IloBoolVarArray y(env, inst.n);
    IloNumVarArray beta(env, inst.n, 0.0, IloInfinity);
    IloNumVar alpha(env, 0.0, IloInfinity);
    IloNumVar eta(env, 0.0, IloInfinity);

    // Objective
    IloObjective obj(env, inst.d1*eta + IloScalProd(x, inst.d_vec) + IloScalProd(lambda, inst.D_vec), IloObjective::Minimize);
    model.add(obj);

    // Constraints
    // Flow conservation
    for (unsigned int i=0; i<inst.n; i++) {
        IloExpr out_arcs_i(env);
        IloExpr in_arcs_i(env);
        for (unsigned int a=0; a<inst.n_arc; a++) {
            if (inst.mat[a].tail == i+1)
                out_arcs_i += x[a];
            if (inst.mat[a].head == i+1)
                in_arcs_i += x[a];
        }
        if (i != inst.t-1) {
            model.add(out_arcs_i == y[i]);
        } else {
            model.add(out_arcs_i == 0);
            // you can't get out of t
        }
        if (i != inst.s-1) {
            model.add(in_arcs_i == y[i]);
        } else {
            model.add(in_arcs_i == 0);
            // you can't get into s
        }
        out_arcs_i.end();
        in_arcs_i.end();
    }
    model.add(y[inst.s-1] == 1);
    model.add(y[inst.t-1] == 1);

    for (unsigned int a = 0; a < inst.n_arc; ++a) {
        model.add(eta + lambda[a] >= inst.d_vec[a]*x[a]);
    }
    model.add(inst.d2*alpha + 2*IloSum(beta) + IloScalProd(y, inst.p) <= inst.S);
    for (unsigned int i=0; i<inst.n; i++) {
        model.add(alpha + beta[i] >= inst.ph[i]*y[i]);
    }

    // Solve
    IloCplex cplex(model);
    cplex.setParam(IloCplex::Param::TimeLimit, time_limit);
    if (verbose < 2) cplex.setOut(env.getNullStream());

    std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
    cplex.solve();
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::chrono::microseconds duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);

    if (cplex.getStatus() == IloAlgorithm::Infeasible) {
        std::cout << inst.name << "," << "dualized,,,,,,,,," << std::endl;
        throw std::domain_error("Infeasible dualized model for instance " + inst.name);
    } else if (cplex.getStatus() == IloAlgorithm::Unknown) {
        std::cout << inst.name << "," << "dualized,,,,,,,,," << std::endl;
        throw std::domain_error("No solution found for instance " + inst.name + ". Maybe not enough time");
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
            std::cout << inst.name << "," << "dualized,,,,,,,,," << std::endl;
            throw std::domain_error("Using arc that does not exist for instance " + inst.name);
        }
    }
    inst.sol.push_back(inst.t);
    path_str += std::to_string(inst.t) + "]";

    double robust_cost = inst.compute_robust_score(env);
    if (abs(robust_cost - cplex.getObjValue()) > 1e-3) {
        std::cout << inst.name << "," << "dualized,,,,,,,,," << std::endl;
        throw std::domain_error("Robust cost and CPLEX cost do not match for instance " + inst.name);
    }

    std::cout << inst.name << ","
        << "dualized,"
        << cplex.getObjValue() << ","
        << cplex.getBestObjValue() << ","
        << static_cast<double>(duration.count()) / 1e6 << ","
        << cplex.getNnodes() << ","
        << inst.compute_robust_constraint(env) << ","
        << inst.compute_static_score() << ","
        << inst.compute_static_constraint() << ","
        << inst.S << ","
        << path_str << std::endl;
}
