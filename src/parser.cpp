#include <fstream>
#include "parser.h"


Instance::Instance(IloEnv env, char filename[]) {
    name = filename;
    p = IloNumArray(env);
    ph = IloNumArray(env);
    d_vec = IloNumArray(env);
    D_vec = IloNumArray(env);
    mat = std::vector<Arc>();

    char readChar;
    int readInt;
    std::ifstream file(filename);
    if (!file) {
        std::cerr << "Error when opening file " << name <<  std::endl;
        exit(1);
    }

    file >> readChar >> readChar >> n;
    file >> readChar >> readChar >> s;
    file >> readChar >> readChar >> t;
    file >> readChar >> readChar >> S;
    file >> readChar >> readChar >> readChar >> d1;
    file >> readChar >> readChar >> readChar >> d2;
    file >> readChar >> readChar;
    for (unsigned int k=0; k<n; k++) {
        file >> readChar >> readInt;
        p.add(readInt);
    }
    // ] ph =
    file >> readChar >> readChar >> readChar >> readChar;
    for (unsigned int k=0; k<n; k++) {
        file >> readChar >> readInt;
        ph.add(readInt);
    }
    // ] mat = [
    file >> readChar >> readChar >> readChar >> readChar >> readChar >> readChar;

    // std::map<Index,double> d_;
    // std::map<Index,double> D_;
    std::vector<std::vector<double>> d_(n, std::vector<double>(n, undefinedValue));
    std::vector<std::vector<double>> D_(n, std::vector<double>(n, undefinedValue));

    while (readChar != ']') {
        Arc v;
        file >> v.tail;
        file >> v.head;
        file >> v.d;
        file >> v.D;
        file >> readChar; // either ';' or ']'
        // d_[Index(v.i-1,v.j-1)] = v.d;
        // D_[Index(v.i-1,v.j-1)] = v.D;
        d_[v.i-1][v.j-1] = v.d;
        D_[v.i-1][v.j-1] = v.D;
        mat.push_back(v);
        d_vec.add(v.d);
        D_vec.add(v.D);
    }
    file.close();
    d = d_;
    D = D_;
    n_arc = mat.size();

    std::vector<std::vector<int>>* neighbors = new std::vector<std::vector<int>>(n);
    std::vector<std::vector<int>>* reverse_neighbors = new std::vector<std::vector<int>>(n);
    for (unsigned int a=0; a<n_arc; a++) {
        neighbors->at(mat[a].i-1).push_back(mat[a].j-1);
        reverse_neighbors->at(mat[a].j-1).push_back(mat[a].i-1);
    }
    neighbors_list = neighbors;
    reverse_neighbors_list = reverse_neighbors;
}


void Instance::display() const {
    std::cout << "n = " << n << std::endl;
    std::cout << "s = " << s << std::endl;
    std::cout << "t = " << t << std::endl;
    std::cout << "S = " << S << std::endl;
    std::cout << "d1 = " << d1 << std::endl;
    std::cout << "d2 = " << d2 << std::endl;
    std::cout << "p = [";
    for (unsigned int i=0; i<n; i++) {
        std::cout << p[i] << " ";
    }
    std::cout << "]" << std::endl;
    std::cout << "ph = [";
    for (unsigned int i=0; i<n; i++) {
        std::cout << ph[i] << " ";
    }
    std::cout << "]" << std::endl;

    std::cout << "mat = [" << std::endl;
    for (unsigned int i=0; i<mat.size(); i++) {
        std::cout << "(" << mat[i].tail << ", " << mat[i].head << ", " << mat[i].d << ", "
            << mat[i].D << ") " << std::endl;
    }
    std::cout << "]" << std::endl;
}


double Instance::compute_static_score(const std::vector<IloInt>& solution) const {
    double static_score = 0.0;
    IloInt current_node = solution[0];
    IloInt next_node = solution[1];
    assert(current_node == s);
    for (unsigned int i=1; i<solution.size(); i++) {
        next_node = solution[i];
        for (unsigned int a=0; a<n_arc; a++) {
            if (mat[a].tail == current_node && mat[a].head == next_node) {
                static_score += mat[a].d;
                break;
            }
        }
        current_node = next_node;
    }
    assert(current_node == t);
    return static_score;
}


double Instance::compute_static_constraint(const std::vector<IloInt>& solution) const {
    double static_constraint = 0.0;
    for (unsigned int i=0; i < solution.size(); i++) {
        static_constraint += p[solution[i]-1];
    }
    return static_constraint;
}


double Instance::compute_robust_score(IloEnv env, const std::vector<IloInt>& solution,
        const unsigned int& time_limit, const int& verbose) const {
    unsigned int n_var = solution.size()-1; // one variable for each arc
    IloModel model(env);

    // Variables
    IloNumVarArray delta1(env, n_var);
    unsigned int start_node;
    unsigned int end_node;
    for (unsigned int k=0; k < n_var; k++) {
        start_node = solution[k];
        end_node = solution[k+1];
        for (unsigned int a = 0; a < n_arc; a++) {
            if (mat[a].tail == start_node && mat[a].head == end_node) {
                delta1[k] = IloNumVar(env, 0.0, D_vec[a]);
                break;
            }
        }
    }

    // Objective
    IloExpr expression_obj(env);
    for (unsigned int k=0; k < n_var; k++) {
        start_node = solution[k];
        end_node = solution[k+1];
        for (unsigned int a = 0; a < n_arc; a++) {
            if (mat[a].tail == start_node && mat[a].head == end_node) {
                expression_obj += mat[a].d * (1 + delta1[k]);
                break;
            }
        }
    }
    IloObjective obj(env, expression_obj, IloObjective::Maximize);
    model.add(obj);
    expression_obj.end();

    // Constraints
    model.add(IloSum(delta1) <= d1);

    IloCplex cplex(model);
    cplex.setParam(IloCplex::Param::TimeLimit, time_limit);
    if (verbose < 2) cplex.setOut(env.getNullStream());

    std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
    cplex.solve();
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::chrono::microseconds duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);

    if (cplex.getStatus() == IloAlgorithm::Infeasible) {
        std::cerr << "No Solution" << std::endl;
        throw std::domain_error("No solution in robust objective problem");
    }
    return cplex.getObjValue();
}


double Instance::compute_robust_score(IloEnv env, const std::vector<IloInt>& solution, std::vector<IloNum>& delta1_vec,
        const unsigned int& time_limit, const int& verbose) const {
    // How to write it cleanly?

    unsigned int n_var = solution.size()-1; // one variable for each arc
    IloModel model(env);

    // Variables
    IloNumVarArray delta1(env, n_var);
    unsigned int start_node;
    unsigned int end_node;
    for (unsigned int k=0; k < n_var; k++) {
        start_node = solution[k];
        end_node = solution[k+1];
        for (unsigned int a = 0; a < n_arc; a++) {
            if (mat[a].tail == start_node && mat[a].head == end_node) {
                delta1[k] = IloNumVar(env, 0.0, D_vec[a]);
                break;
            }
        }
    }

    // Objective
    IloExpr expression_obj(env);
    for (unsigned int k=0; k < n_var; k++) {
        start_node = solution[k];
        end_node = solution[k+1];
        for (unsigned int a = 0; a < n_arc; a++) {
            if (mat[a].tail == start_node && mat[a].head == end_node) {
                expression_obj += mat[a].d * (1 + delta1[k]);
                break;
            }
        }
    }
    IloObjective obj(env, expression_obj, IloObjective::Maximize);
    model.add(obj);
    expression_obj.end();

    // Constraints
    model.add(IloSum(delta1) <= d1);

    IloCplex cplex(model);
    cplex.setParam(IloCplex::Param::TimeLimit, time_limit);
    if (verbose < 2) cplex.setOut(env.getNullStream());

    std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
    cplex.solve();
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::chrono::microseconds duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);

    if (cplex.getStatus() == IloAlgorithm::Infeasible) {
        std::cerr << "No Solution" << std::endl;
        throw std::domain_error("No solution in robust objective problem");
    }

    // Get the values of delta1
    delta1_vec = std::vector<IloNum>(n, 0.0);
    for (unsigned int k=0; k < n_var; k++) {
        IloInt node = solution[k];
        delta1_vec[node-1] = cplex.getValue(delta1[node-1]);
    }
    return cplex.getObjValue();
}



double Instance::compute_robust_constraint(IloEnv env, const std::vector<IloInt>& solution,
        const unsigned int& time_limit,const int& verbose) const {
    unsigned int n_sol = solution.size();
    IloModel model(env);

    // Variables
    IloNumVarArray delta2(env, n_sol, 0.0, 2.0);

    // Objective
    IloExpr expression_obj(env);
    unsigned int node;
    for (unsigned int k = 0; k < n_sol; k++) {
        node = solution[k];
        expression_obj += p[node-1] + ph[node-1]*delta2[k];
    }
    IloObjective obj(env, expression_obj, IloObjective::Maximize);
    model.add(obj);
    expression_obj.end();

    // Constraints
    model.add(IloSum(delta2) <= d2);

    IloCplex cplex(model);
    cplex.setParam(IloCplex::Param::TimeLimit, time_limit);
    if (verbose < 2) cplex.setOut(env.getNullStream());

    std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
    cplex.solve();
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::chrono::microseconds duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);

    if (cplex.getStatus() == IloAlgorithm::Infeasible){
        std::cerr << "No Solution for robust constraint with instance " << name << std::endl;
        throw std::domain_error("No solution in robust constraint problem");
    }
    if (verbose >= 1) {
        std::cout << "robust constraint: " << cplex.getObjValue() << std::endl;
        std::cout << "time: " << static_cast<double>(duration.count()) / 1e6 << std::endl;
    }
    return cplex.getObjValue();
}
