#include "branch_and_cut.h"
#include "parser.h"


ILOLAZYCONSTRAINTCALLBACK7(myCallBack, const IloBoolVarArray&, x, const IloBoolVarArray&, y,
        const IloNumVar&, z, Instance&, inst, const unsigned int&, verbose, int&, nIteration,
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

    try {
        IloExpr exprObjective(env);
        double robust_score = Subproblems::solve_objective_subproblem(env, xValues, x, inst, exprObjective);
        if (robust_score > zValue + TOL) {
            add(exprObjective <= z);
        }
        exprObjective.end();
    } catch (std::exception& e) {
        std::cerr << "Error in lazy constraint callback for objective at iter : " << nIteration
             << " " << e.what() << std::endl;
        throw e;
    }

    try {
        IloExpr exprConstraint(env);
        double robust_constraint = Subproblems::solve_constraint_subproblem(env, yValues, y, inst, exprConstraint);
        if (robust_constraint > inst.S + TOL) {
            add(exprConstraint <= inst.S);
        }
        exprConstraint.end();
    } catch (std::exception& e) {
        std::cerr << "Error in lazy constraint callback for constraint at iter : " << nIteration
             << " " << e.what() << std::endl;
        throw e;
    }

    try {
        // To brake symmetry, we forbid paths that are equivalent to the current one
        // if current solution takes a->b->d and there is a->c->d possible with the same cost
        // we forbid a->c->d
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
                add(x[to_forbid[k][2]] + x[to_forbid[k][3]] <= 1);
            }
        }
    } catch (std::exception& e) {
        std::cerr << "Error in lazy constraint callback breaking symmetries: " << e.what() << std::endl;
        throw e;
    }

    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::chrono::microseconds duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    double spentTime_ = static_cast<double>(duration.count()) / 1e6;
    if (verbose > 0) std::cout << "Time: " << spentTime_ << "s" << std::endl;
    spentTime += spentTime_;

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

    cplexCheckStatus(cplex, inst);
    // Retrieve solution
    IloNumArray xValues(env);
    cplex.getValues(xValues, x);
    retrieveCplexSolution(xValues, inst);
    xValues.end();

    nodesExplored = cplex.getNnodes();
    infBound = cplex.getBestObjValue();
}
