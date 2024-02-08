#include "branch_and_cut.h"
#include "parser.h"


struct VarCallBack {
    // Because CallBack can not take more than 7 arguments...
    VarCallBack(IloBoolVarArray& x_, IloBoolVarArray& y_, IloNumVar& z_) {
        x = x_;
        y = y_;
        z = z_;
    }
    IloBoolVarArray x;
    IloBoolVarArray y;
    IloNumVar z;
};

struct BestSolFound {
    BestSolFound(IloNumArray& best_xValues_, IloNumArray& best_yValues_, double& best_score_,
            double& less_violated_constraint_, double& best_violated_score_, double best_bound_,
            bool& admissible_solution_found_) {
        best_xValues = best_xValues_;
        best_yValues = best_yValues_;
        best_score = best_score_;
        less_violated_constraint = less_violated_constraint_;
        best_violated_score = best_violated_score_;
        best_bound = best_bound_;
        admissible_solution_found = admissible_solution_found_;
    }
    IloNumArray best_xValues;
    IloNumArray best_yValues;
    double best_score;
    double less_violated_constraint;
    double best_violated_score;
    double best_bound;
    bool admissible_solution_found;
};


ILOLAZYCONSTRAINTCALLBACK7(myCallBack, const VarCallBack&, varCallback, Instance&, inst,
        BestSolFound&, bestSolFound, const unsigned int&, verbose, int&, nIteration, double&, spentTime,
        const IloCplex&, cplex) {
    nIteration++;
    if (verbose > 0) std::cout << "Lazy constraint callback, it " << nIteration << std::endl;
    std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();

    IloEnv env = getEnv();
    IloNumArray xValues(env);
    IloNumArray yValues(env);

    getValues(xValues, varCallback.x);
    getValues(yValues, varCallback.y);
    double zValue = getValue(varCallback.z);

    IloExpr exprObjective(env);
    double robust_objective = Subproblems::solve_objective_subproblem(env, xValues, varCallback.x, inst, exprObjective);
    if (robust_objective > zValue + TOL) {
        add(exprObjective <= varCallback.z);
    }
    exprObjective.end();

    IloExpr exprConstraint(env);
    double robust_constraint = Subproblems::solve_constraint_subproblem(env, yValues, varCallback.y, inst, exprConstraint);
    bool violated_constraint = robust_constraint > inst.S + TOL;
    if (violated_constraint) {
        add(exprConstraint <= inst.S);
    }
    exprConstraint.end();


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
            add(varCallback.x[to_forbid[k][2]] + varCallback.x[to_forbid[k][3]] <= 1);
        }
    }

    // Update best solution found
    double current_bound = cplex.getBestObjValue();
    if (bestSolFound.best_bound < current_bound) {
        bestSolFound.best_bound = current_bound;
    }
    if ((!violated_constraint && robust_objective < bestSolFound.best_score)
            || (!bestSolFound.admissible_solution_found && robust_constraint < bestSolFound.less_violated_constraint)
            || (!bestSolFound.admissible_solution_found && robust_constraint == bestSolFound.less_violated_constraint && robust_objective < bestSolFound.best_violated_score)) {
        bestSolFound.admissible_solution_found = bestSolFound.admissible_solution_found || !violated_constraint;
        if (bestSolFound.admissible_solution_found) {
            bestSolFound.best_score = robust_objective;
        } else {
            bestSolFound.best_violated_score = robust_objective;
            bestSolFound.less_violated_constraint = robust_constraint;
        }

        for (unsigned int a = 0; a < inst.n_arc; ++a) {
            bestSolFound.best_xValues[a] = xValues[a];
        }
        for (unsigned int i = 0; i < inst.n; ++i) {
            bestSolFound.best_yValues[i] = yValues[i];
        }
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
    // To retrieve manually the bestSolution (even when it is not admissible because of the robust constraint)
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
    bool admissible_solution_found = false;
    double best_bound = 0.0;
    BestSolFound bestSolFound(best_xValues, best_yValues, best_score, less_violated_constraint,
        best_violated_score, best_bound, admissible_solution_found);
    //

    IloModel model(env);

    // Variables
    IloBoolVarArray x(env, inst.n_arc);
    IloBoolVarArray y(env, inst.n);
    IloNumVar z(env, 0.0, IloInfinity, "z");

    VarCallBack varCallback(x,y,z);

    // Objective
    IloObjective obj(env, z, IloObjective::Minimize);
    model.add(obj);

    // Constraints
    add_static_constraints(env, model, x, y, z, inst);

    IloCplex cplex(model);
    parametrizeCplex(cplex, time_limit, verbose);
    cplex.use(myCallBack(env,varCallback,inst,bestSolFound,verbose,nCallBacks,callBacksTimeSpan,cplex));
    cplex.solve();

    if (cplex.getStatus() == IloAlgorithm::Infeasible) {
        throw std::domain_error("Infeasible " + method_name + " model for instance " + inst.name);
    }
    // Retrieve solution
    retrieveCplexSolution(bestSolFound.best_xValues, inst);

    nodesExplored = cplex.getNnodes();
    infBound = cplex.getBestObjValue();

    best_xValues.end();
    best_yValues.end();
}
