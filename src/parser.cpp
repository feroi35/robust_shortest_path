#include "parser.h"
#include <chrono>


Instance::Instance(IloEnv env, char filename[]) {
    // I'm not sure why the instance depends of the environment
    // and I think it make more sense not to do so
    // But that's what the professor did in the TP
    name = filename;
    p = IloNumArray(env);
    ph = IloNumArray(env);
    mat = std::vector<Arc>();
    
    char readChar;
    int readInt;
    std::ifstream file(filename);

    if (!file) {
        std::cerr << "Error when opening file" << std::endl;
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

    std::vector<std::vector<double>> d_(n, std::vector<double>(n, undefinedValue));
    std::vector<std::vector<double>> D_(n, std::vector<double>(n, undefinedValue));
    // débattable, c'est pour avoir une première version qui tourne

    while (readChar != ']') {
        Arc v;
        file >> v.i;
        file >> v.j;
        file >> v.d;
        file >> v.D;
        file >> readChar; // either ';' or ']'
        mat.push_back(v);
        d_[v.i-1][v.j-1] = v.d;
        D_[v.i-1][v.j-1] = v.D;
    }
    file.close();
    d = d_;
    D = D_;
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
        std::cout << "(" << mat[i].i << ", " << mat[i].j << ", " << mat[i].d << ", " 
            << mat[i].D << ") " << std::endl;
    }
    std::cout << "]" << std::endl;
    std::cout << "d = [" << std::endl;
    for (unsigned int i=0; i<d.size(); i++) {
        std::cout << "[";
        for (unsigned int j=0; j<d[i].size(); j++) {
            std::cout << d[i][j] << " ";
        }
        std::cout << "]" << std::endl;
    }
}

double Instance::compute_static_score(std::vector<IloInt> solution) const {
    double static_score = 0.0;
    IloInt current_node = solution[0];
    assert(current_node == s); 
    for (unsigned int i=1; i<solution.size(); i++) {
        IloInt next_node = solution[i];
        static_score += d[current_node-1][next_node-1];
        current_node = next_node;
    }
    assert(current_node == t);
    return static_score;
}

double Instance::compute_robust_score(IloEnv env, std::vector<IloInt> solution, unsigned int time_limit) const {
    
    IloModel model(env);

    IloArray<IloBoolVarArray> delta1(env, n);
    std::stringstream name;

      // Create variables delta1
    for(unsigned int i = 0; i < n; ++i) {
        delta1[i] = IloBoolVarArray(env, n);
        for(unsigned int j = 0; j < n; ++j) {
        name << "delta1_" << i << j;
        delta1[i][j] = IloBoolVar(env, name.str().c_str());
        name.str("");
        }
    }

    IloExpr expression_obj(env);
    unsigned int start_node;
    unsigned int end_node;

    for(unsigned int k = 0; k < solution.size()-1; k++) {
        start_node = solution[k];
        end_node = solution[k+1];
        expression_obj += d[start_node-1][end_node-1]*(1+delta1[start_node-1][end_node-1]);
    }

    IloObjective obj(env, expression_obj, IloObjective::Maximize);
    model.add(obj);
    expression_obj.end();

    // Constraints
    IloExpr expression_cstr(env);
    for (unsigned int i=0; i<n; i++) {
        for (unsigned int j=0; j<n; j++) {
            expression_cstr += delta1[i][j];
            model.add(delta1[i][j] <= D[i][j]);
        }
    }
    model.add(expression_cstr <= d1);
    expression_cstr.end();

    IloCplex cplex(model);
    cplex.setParam(IloCplex::Param::TimeLimit, time_limit);
    // cplex.setOut(env.getNullStream());

    std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
    cplex.solve();
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::chrono::microseconds duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);

    if (cplex.getStatus() == IloAlgorithm::Infeasible){
        cout << "No Solution" << endl;
        throw std::domain_error("No solution in robust objective");
        return 0;
    }
    else{
        std::cout << "robust objective: " << cplex.getObjValue() << std::endl;
        std::cout << "time: " << static_cast<double>(duration.count()) / 1e6 << std::endl;
        return cplex.getObjValue();        
    }
}

double Instance::compute_robust_constraint(IloEnv env, std::vector<IloInt> solution, unsigned int time_limit) const {
    
    IloModel model(env);

    IloBoolVarArray delta2(env, n);
    std::stringstream name;

      // Create variables delta2
    for(unsigned int i = 0; i < n; ++i) {
        name << "delta2_" << i;
        delta2[i] = IloBoolVar(env, name.str().c_str());
        name.str("");
    }

    IloExpr expression_obj(env);
    unsigned int node;

    for(unsigned int k = 0; k < solution.size()-1; k++) {
        node = solution[k];
        expression_obj += p[node-1] + ph[node-1]*delta2[node-1];
    }

    IloObjective obj(env, expression_obj, IloObjective::Maximize);
    model.add(obj);
    expression_obj.end();

    // Constraints
    IloExpr expression_cstr(env);
    for (unsigned int i=0; i<n; i++) {
        expression_cstr += delta2[i];
        model.add(delta2[i] <= 2);
    }
    model.add(expression_cstr <= d2);
    expression_cstr.end();

    IloCplex cplex(model);
    cplex.setParam(IloCplex::Param::TimeLimit, time_limit);
    // cplex.setOut(env.getNullStream());

    std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
    cplex.solve();
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::chrono::microseconds duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);

    if (cplex.getStatus() == IloAlgorithm::Infeasible){
        cout << "No Solution" << endl;
        throw std::domain_error("No solution in robust objective");
        return 0;
    }
    else{
        std::cout << "robust constraint: " << cplex.getObjValue() << std::endl;
        std::cout << "time: " << static_cast<double>(duration.count()) / 1e6 << std::endl;
        return cplex.getObjValue();        
    }
}
