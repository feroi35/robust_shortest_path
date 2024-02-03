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

    std::vector<std::vector<float>> d_(n, std::vector<float>(n, undefinedValue));
    std::vector<std::vector<float>> D_(n, std::vector<float>(n, undefinedValue));
    std::vector<std::vector<int>> map_mat_(n, std::vector<int>(n, -1));
    int index = 0;
    while (readChar != ']') {
        Arc v;
        file >> v.tail;
        file >> v.head;
        file >> v.d;
        file >> v.D;
        file >> readChar; // either ';' or ']'
        d_[v.tail-1][v.head-1] = v.d;
        D_[v.tail-1][v.head-1] = v.D;
        mat.push_back(v);
        d_vec.add(v.d);
        D_vec.add(v.D);
        map_mat_[v.tail-1][v.head-1] = index;
        index++;
    }
    file.close();
    d = d_;
    D = D_;
    map_mat = map_mat_;
    n_arc = mat.size();

    std::vector<std::vector<int>> neighbors(n, std::vector<int>());
    std::vector<std::vector<int>> reverse_neighbors(n, std::vector<int>());
    for (unsigned int a=0; a<n_arc; a++) {
        neighbors[mat[a].tail-1].push_back(mat[a].head-1);
        reverse_neighbors[mat[a].head-1].push_back(mat[a].tail-1);
    }
    neighbors_list = neighbors;
    reverse_neighbors_list = reverse_neighbors;
}


Instance::Instance(const Instance& instan){
    name = instan.name;
    n = instan.n;
    n_arc = instan.n_arc;
    s = instan.s;
    t = instan.t;
    S = instan.S;
    d1 = instan.d1;
    d2 = instan.d2;
    p = instan.p;
    ph = instan.ph;
    d_vec = instan.d_vec;
    D_vec = instan.D_vec;
    d = instan.d;
    D = instan.D;
    mat = instan.mat;
    neighbors_list = instan.neighbors_list;
    reverse_neighbors_list = instan.reverse_neighbors_list;
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
    if (solution.empty()) {
        throw std::domain_error("Empty solution");
    }
    double static_score = 0.0;
    IloInt current_node = solution[0];
    IloInt next_node;
    if (current_node != s) {
        throw std::domain_error("First node of solution is not s");
    }
    for (unsigned int k=0; k < solution.size()-1; k++) {
        next_node = solution[k+1];
        if (d[current_node-1][next_node-1] != d[current_node-1][next_node-1]) {
            throw std::domain_error("Nan Value: no arc between " + std::to_string(current_node) + " and " + std::to_string(next_node));
        }
        static_score += d[current_node-1][next_node-1];
        current_node = next_node;
    }
    if (current_node != t) {
        throw std::domain_error("Last node of solution is not t");
    }
    return static_score;
}


double Instance::compute_robust_score_milp(IloEnv env, const std::vector<IloInt>& solution,
        const unsigned int& time_limit, const unsigned int& verbose) const {
    std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
    if (solution.empty()) {
        throw std::domain_error("Empty solution");
    }
    unsigned int n_var = solution.size()-1; // one variable for each arc
    if (n_var == 0) {
        if (solution[0] == s && solution[0] == t) {
            return 0.0;
        } else {
            throw std::domain_error("Solution of size 1 is not s=t");
        }
    }

    IloModel model(env);

    IloNumVarArray delta1(env, n_var); // Variables delta1
    IloExpr expression_obj(env); // Expression of the objective

    IloInt current_node = solution[0];
    IloInt next_node;
    if (current_node != s) {
        throw std::domain_error("First node of solution is not s");
    }
    for (unsigned int k=0; k < n_var; k++) {
        next_node = solution[k+1];
        if (D[current_node-1][next_node-1] != D[current_node-1][next_node-1] || d[current_node-1][next_node-1] != d[current_node-1][next_node-1]) {
            throw std::domain_error("No arc between " + std::to_string(current_node) + " and " + std::to_string(next_node));
        }
        delta1[k] = IloNumVar(env, 0.0, D[current_node-1][next_node-1]);
        expression_obj += d[current_node-1][next_node-1] * (1 + delta1[k]);
        current_node = next_node;
    }
    if (current_node != t) {
        throw std::domain_error("Last node of solution is not t");
    }
    IloObjective obj(env, expression_obj, IloObjective::Maximize);
    model.add(obj);
    expression_obj.end();

    model.add(IloSum(delta1) <= d1);

    IloCplex cplex(model);
    cplex.setParam(IloCplex::Param::TimeLimit, time_limit);
    if (verbose < 3) cplex.setOut(env.getNullStream());

    cplex.solve();

    if (cplex.getStatus() == IloAlgorithm::Infeasible) {
        throw std::domain_error("No solution in robust objective problem");
    } else if (cplex.getStatus() == IloAlgorithm::Unknown) {
        throw std::domain_error("No solution found for instance " + name + ". Maybe not enough time");
    } else if (cplex.getStatus() != IloAlgorithm::Optimal) {
        throw std::domain_error("Solution of subproblem not optimal for instance " + name);
    }
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::chrono::microseconds duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    if (verbose >= 2) {
        std::cout << "(Time: " << duration.count() << " microseconds) ";
    }
    return cplex.getObjValue();
}


double Instance::compute_robust_score_knapsack(const std::vector<IloInt>& solution, const unsigned int& verbose) const {
    std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
    if (solution.empty()) {
        throw std::domain_error("Empty solution");
    }
    unsigned int n_edges = solution.size()-1;
    if (n_edges == 0) {
        if (solution[0] == s && solution[0] == t) {
            return 0.0;
        } else {
            throw std::domain_error("Solution of size 1 is not s=t");
        }
    }

    std::vector<IloNum> weights(n_edges, 0.0);
    std::vector<IloNum> uncertainties(n_edges, 0.0);

    unsigned int current_node = solution[0];
    unsigned int next_node;
    if (current_node != s) {
        throw std::domain_error("First node of solution is not s");
    }
    double static_score = 0.0;
    for (unsigned int k=1; k < n_edges+1; k++) {
        next_node = solution[k];
        weights[k-1] = d[current_node-1][next_node-1];
        uncertainties[k-1] = D[current_node-1][next_node-1];
        if (weights[k-1] != weights[k-1] || uncertainties[k-1] != uncertainties[k-1]) {
            throw std::domain_error("Arc between " + std::to_string(current_node) + ", " + std::to_string(next_node) + " does not exist");
        }
        static_score += d[current_node-1][next_node-1];
        current_node = next_node;

    }
    if (current_node != t) {
        throw std::domain_error("Last node of solution is not t");
    }

    std::vector<size_t> argsorted_weights = argsort(weights);
    double robust_attack = 0.0;
    double used_budget = 0.0;
    int idx = n_edges-1;
    while (used_budget < d1 && idx >= 0) {
        unsigned int arc_idx = argsorted_weights[idx];
        float delta1_i = std::min(d1 - used_budget, uncertainties[arc_idx]);
        used_budget += delta1_i;
        robust_attack += delta1_i * weights[arc_idx];
        idx--;
    }
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::chrono::microseconds duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    if (verbose >= 2) {
        std::cout << "(Time: " << duration.count() << " microseconds) ";
    }
    return static_score + robust_attack;
}


double Instance::compute_static_constraint(const std::vector<IloInt>& solution) const {
    double static_constraint = 0.0;
    for (unsigned int i=0; i < solution.size(); i++) {
        static_constraint += p[solution[i]-1];
    }
    return static_constraint;
}


double Instance::compute_robust_constraint_milp(IloEnv env, const std::vector<IloInt>& solution,
        const unsigned int& time_limit,const unsigned int& verbose) const {
    std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
    if (solution.empty()) {
        throw std::domain_error("Empty solution");
    }
    unsigned int n_sol = solution.size();

    IloModel model(env);

    IloNumVarArray delta2(env, n_sol, 0.0, 2.0);
    IloExpr expression_obj(env);
    IloInt node = solution[0];
    if (node != s) {
        throw std::domain_error("First node of solution is not s");
    }
    for (unsigned int k = 0; k < n_sol; k++) {
        node = solution[k];
        expression_obj += p[node-1] + ph[node-1]*delta2[k];
    }
    if (node != t) {
        throw std::domain_error("Last node of solution is not t");
    }
    IloObjective obj(env, expression_obj, IloObjective::Maximize);
    model.add(obj);
    expression_obj.end();

    model.add(IloSum(delta2) <= d2);

    IloCplex cplex(model);
    cplex.setParam(IloCplex::Param::TimeLimit, time_limit);
    if (verbose < 3) cplex.setOut(env.getNullStream());

    cplex.solve();

    if (cplex.getStatus() == IloAlgorithm::Infeasible){
        throw std::domain_error("No solution in robust constraint problem");
    } else if (cplex.getStatus() == IloAlgorithm::Unknown) {
        throw std::domain_error("No solution found for instance " + name + ". Maybe not enough time");
    } else if (cplex.getStatus() != IloAlgorithm::Optimal) {
        throw std::domain_error("Solution of subproblem not optimal for instance " + name);
    }
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::chrono::microseconds duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    if (verbose >= 2) {
        std::cout << "(Time: " << duration.count() << " microseconds) ";
    }
    return cplex.getObjValue();
}


double Instance::compute_robust_constraint_knapsack(const std::vector<IloInt>& solution, const unsigned int& verbose) const {
    std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
    if (solution.empty()) {
        throw std::domain_error("Empty solution");
    }
    unsigned int n_cities = solution.size();
    std::vector<IloNum> uncertain_weights(n_cities, 0.0);

    unsigned int current_node = solution[0];
    if (current_node != s) {
        throw std::domain_error("First node of solution is not s");
    }
    double static_constraint = 0.0;
    for (unsigned int k=0; k < n_cities; k++) {
        current_node = solution[k];
        uncertain_weights[k] = ph[current_node-1];
        static_constraint += p[current_node-1];
    }
    if (current_node != t) {
        throw std::domain_error("Last node of solution is not t");
    }

    std::vector<size_t> argsorted_weights = argsort(uncertain_weights);
    double robust_constraint = 0.0;
    double used_budget = 0.0;
    int idx = n_cities-1;
    while (used_budget < d2 && idx >= 0) {
        unsigned int arc_idx = argsorted_weights[idx];
        float delta2_i = std::min(d2 - used_budget, 2.0);
        used_budget += delta2_i;
        robust_constraint += delta2_i * uncertain_weights[arc_idx];
        idx--;
    }
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::chrono::microseconds duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    if (verbose >= 2) {
        std::cout << "(Time: " << duration.count() << " microseconds) ";
    }
    return static_constraint + robust_constraint;
}



std::vector<std::vector<int>> arcs_to_forbid(const Instance& inst, const int& i, const int& j){
    std::vector<tuple<Arc2,Arc2,float,float>> sub_paths;

    for (unsigned int l = 0; l<inst.neighbors_list[i].size(); ++l) {
        int k = inst.neighbors_list[i][l];
        if (std::count(inst.reverse_neighbors_list[j].begin(), inst.reverse_neighbors_list[j].end(), k)){
            Arc2 arc1 = Arc2(k,i,inst.d[i][k],inst.D[i][k]);
            Arc2 arc2 = Arc2(j,k,inst.d[k][j],inst.D[k][j]);
            sub_paths.push_back(std::make_tuple(arc1,arc2,inst.p[k],inst.ph[k]));
        }
    }

    std::vector<std::vector<int>> to_forbid;
    std::vector<bool> already_forbade(sub_paths.size(),false);

    for(unsigned int l=0; l<sub_paths.size(); l++){
        if (already_forbade[l]) continue;
        bool found_better = false;
        for(unsigned int m=0; m<sub_paths.size(); m++){
            if (m==l || already_forbade[m]) continue;
            if (std::get<2>(sub_paths[m])==std::get<2>(sub_paths[l])){
                if (std::get<3>(sub_paths[m])==std::get<3>(sub_paths[l])){
                    bool same_ds = (std::get<0>(sub_paths[m]).d == std::get<0>(sub_paths[l]).d) && (std::get<1>(sub_paths[m]).d == std::get<1>(sub_paths[l]).d) && (std::get<0>(sub_paths[m]).D == std::get<0>(sub_paths[l]).D) && (std::get<1>(sub_paths[m]).D == std::get<1>(sub_paths[l]).D);
                    bool inverted_ds = (std::get<0>(sub_paths[m]).d == std::get<1>(sub_paths[l]).d) && (std::get<1>(sub_paths[m]).d == std::get<0>(sub_paths[l]).d) && (std::get<0>(sub_paths[m]).D == std::get<1>(sub_paths[l]).D) && (std::get<1>(sub_paths[m]).D == std::get<0>(sub_paths[l]).D);
                    if (same_ds || inverted_ds){
                            if(std::get<0>(sub_paths[m]).head > std::get<0>(sub_paths[l]).head){
                                already_forbade[m] = true;
                                int k = std::get<0>(sub_paths[m]).head;
                                int index_first_arc = inst.map_mat[i][k];
                                int index_second_arc = inst.map_mat[k][j];
                                std::vector<int> to_forbid_l{i, j, index_first_arc, index_second_arc};
                                to_forbid.push_back(to_forbid_l);
                            }
                            else{
                                found_better = true;
                            }
                        }
                    }
                }
            }
        if (found_better==true){
            already_forbade[l] = true;
            int k = std::get<0>(sub_paths[l]).head;
            int index_first_arc = inst.map_mat[i][k];
            int index_second_arc = inst.map_mat[k][j];
            std::vector<int> to_forbid_l{i, j, index_first_arc, index_second_arc};
            to_forbid.push_back(to_forbid_l);
        }
    }
    return to_forbid;
}

std::vector<std::vector<int>> arcs_to_forbid(const Instance& inst){

    std::vector<std::vector<int>> to_forbid;
    for (unsigned int i=0; i<inst.n; i++){
        for (unsigned int j=0; j<inst.n; j++){
            if (i != j){
                std::vector<std::vector<int>> to_forbid_ij = arcs_to_forbid(inst,i,j);
                to_forbid.insert(to_forbid.end(), to_forbid_ij.begin(), to_forbid_ij.end());
            }
        }
    }
    return to_forbid;
}