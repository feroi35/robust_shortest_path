#include "dualized_formulation.h"
#include <chrono>

void dualized_solve(IloEnv env, Instance& inst, const unsigned int& time_limit, const int& verbose) {

    IloModel model(env);

    IloArray<IloBoolVarArray> x(env, inst.n);
    IloArray<IloNumVarArray> lambda(env, inst.n);
    IloBoolVarArray y(env, inst.n);
    IloNumVarArray beta(env, inst.n);

    std::stringstream name;

    name << "alpha";
    IloNumVar alpha(env, 0.0, IloInfinity, name.str().c_str());
    name.str("");
    name << "eta";
    IloNumVar eta(env, 0.0, IloInfinity, name.str().c_str());
    name.str("");

    // Create variables x
    for(unsigned int i = 0; i < inst.n; ++i) {
        x[i] = IloBoolVarArray(env, inst.n);
        lambda[i] = IloNumVarArray(env, inst.n);
        for(unsigned int j = 0; j < inst.n; ++j) {
        name << "x_" << i << j;
        x[i][j] = IloBoolVar(env, name.str().c_str());
        name.str("");
        name << "lambda_" << i << j;
        lambda[i][j] = IloNumVar(env, 0.0, IloInfinity, name.str().c_str());
        name.str("");
        }
    }
    // Create variables y
    for(unsigned int i = 0; i < inst.n; ++i) {
        name << "y_" << i;
        y[i] = IloBoolVar(env, name.str().c_str());
        name.str("");
        name << "beta_" << i;
        beta[i] = IloNumVar(env, 0.0, IloInfinity, name.str().c_str());
        name.str("");
    }

    IloExpr expression_obj(env);

    expression_obj += inst.d1*eta;
    for(unsigned int k = 0; k < inst.mat.size(); k++) {
        Arc v = inst.mat[k];
        expression_obj += v.d*x[v.i-1][v.j-1] + v.D*lambda[v.i-1][v.j-1];
    }
    // Fixer les x[i][j] = 0 si on n'a pas d'arc entre i et j ?
    // format de l'instance à revoir, c'est pour voir si ça tourne
    IloObjective obj(env, expression_obj, IloObjective::Minimize);
    model.add(obj);
    expression_obj.end();

    // Constraints

    for(unsigned int i=0; i<inst.n; i++) {
        for(unsigned int j=0; j<inst.n; j++) {
        model.add(eta + lambda[i][j] >= inst.d[i][j]*x[i][j]);
        }
    }

    IloExpr expression_cstr(env);
    for(unsigned int i=0; i<inst.n; i++) {
        expression_cstr += 2*beta[i] + inst.p[i]*y[i];
    }
    model.add(inst.d2*alpha + expression_cstr <= inst.S);
    expression_cstr.end();

    for(unsigned int i=0; i<inst.n; i++){
        model.add(alpha + beta[i] >= inst.ph[i]*y[i]);
    }

    for (unsigned int i=0; i<inst.n; i++) {
        IloExpr expression_cstr(env);
        for (unsigned int j=0; j<inst.n; j++) {
        if (inst.d[i][j] == inst.d[i][j]) {
            // not nan
            expression_cstr += x[i][j];
        }
        }
        if (i != inst.t-1) {
        model.add(expression_cstr == y[i]);
        } else {
        model.add(expression_cstr == 0);
        // On ne peut pas sortir de t
        }
        expression_cstr.end();
    }

    for (unsigned int j=0; j<inst.n; j++) {
        IloExpr expression_cstr(env);
        for (unsigned int i=0; i<inst.n; i++) {
        if (inst.d[i][j] == inst.d[i][j])
            expression_cstr += x[i][j];
        }
        if (j != inst.s-1) {
        model.add(expression_cstr == y[j]);
        } else {
        model.add(expression_cstr == 0);
        // On ne peut pas entrer dans s
        }
        expression_cstr.end();
    }

    model.add(y[inst.s-1] == 1);
    model.add(y[inst.t-1] == 1);

    IloCplex cplex(model);
    cplex.setParam(IloCplex::Param::TimeLimit, time_limit);
   if(verbose <2)
        cplex.setOut(env.getNullStream());

    std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
    cplex.solve();
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::chrono::microseconds duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);

    if (cplex.getStatus() == IloAlgorithm::Infeasible)
        cout << "No dualized Solution for file " << inst.name << endl;
    else{

        if(verbose >= 1){
            std::cout << "objective: " << cplex.getObjValue() << std::endl;
            std::cout << "time: " << static_cast<double>(duration.count()) / 1e6 << std::endl;
            std::cout << "path: ";
        }

        assert(inst.sol.empty());
        unsigned int current_node = inst.s-1;
        while (current_node != inst.t-1) {
            if(verbose >=1){
                std::cout << current_node+1 << " ";
            }
            inst.sol.push_back(current_node+1);
            for (unsigned int j=0; j<inst.n; j++) {
                if ((inst.d[current_node][j] == inst.d[current_node][j])
                    && (cplex.getValue(x[current_node][j]) == 1)) {
                current_node = j;
                break;
                }
            }
        }
        if(verbose >=1)
            std::cout << inst.t << std::endl;
        inst.sol.push_back(inst.t);

        if(verbose >=1){
            std::cout << "size sol: " << inst.sol.size() << std::endl;
            for (unsigned int i=0; i<inst.sol.size(); i++) {
                std::cout << inst.sol[i] << " ";
            }
            std::cout << std::endl;
        }

        // A homogénéiser, selon si on crée un csv ou pas
        // en faire un paramètre?
        std::cout << inst.name << ", " << cplex.getObjValue() << ", " 
            << static_cast<double>(duration.count()) / 1e6 <<
            ", " << cplex.getNnodes() << ", " << cplex.getBestObjValue() << std::endl; 
    }
}