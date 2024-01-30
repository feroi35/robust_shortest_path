#include "parser.h"
#include "heuristics.h"


bool cmp_index_val(const std::tuple<int, float>& a, const std::tuple<int, float>& b) {
    return std::get<1>(a) > std::get<1>(b);
}

bool cmp_dij(const std::tuple<float, float>& a, const std::tuple<float, float>& b) {
    return std::get<0>(a) > std::get<0>(b);
}

float node_info::sum_knapsack_capa() const{
    float sum = 0;
    for(auto it = knapsack_dij.begin(); it != knapsack_dij.end(); ++it){
        sum += std::get<1>(*it);
    }
    return sum;
}

void node_info::compute_dist_robust(const float d1){
    float sum = 0;
    float weight = 0;
    for(auto it = knapsack_dij.begin(); it != knapsack_dij.end(); ++it){
        if (weight + std::get<1>(*it) <= d1){
            sum += std::get<0>(*it)*std::get<1>(*it);
            weight += std::get<1>(*it);
        }
        else{
            sum += std::get<0>(*it)*(d1-weight);
            break;
        }
    }
    robust_dist = sum;
}

void node_info::compute_nodes_robust(const float d2){
    float sum = 0;
    float weight = 0;
    for(auto it = knapsack_phi.begin(); it != knapsack_phi.end(); ++it){
        if (weight + 2 <= d2){
            sum += 2*(*it);
            weight += 2;
        }
        else{
            sum += (*it)*(d2-weight);
            break;
        }
    }
    robust_dist_nodes = sum;
}


void node_info::redo_knapsack_dij(const Instance& inst, const int& new_node){

    float new_dij = inst.d[parent][new_node];
    float new_Dij = inst.D[parent][new_node];

    knapsack_dij.push_back(std::make_tuple(new_dij, new_Dij));
    std::sort(knapsack_dij.begin(), knapsack_dij.end(), cmp_dij);

    float sum_knapsack = sum_knapsack_capa();

    while(sum_knapsack-std::get<1>(knapsack_dij.back()) >= inst.d1){
        std::tuple<float,float> last = knapsack_dij.back();
        knapsack_dij.pop_back();
        sum_knapsack -= std::get<1>(last);
    }

    compute_dist_robust(inst.d1);
}

void node_info::redo_knapsack_phi(const Instance& inst, const int& new_node){
    float new_phi = inst.ph[new_node];

    knapsack_phi.push_back(new_phi);
    std::sort(knapsack_phi.begin(), knapsack_phi.end(), std::greater<float>());

    float sum_knapsack = 2*knapsack_phi.size();

    while(sum_knapsack-2 >= inst.d2){
        knapsack_phi.pop_back();
        sum_knapsack -= 2;
    }

    compute_nodes_robust(inst.d2);
}

double compute_static_score(const Instance& inst, const std::vector<int>& sol){
    double score = 0;
    for (unsigned int i=0; i<sol.size()-1; i++) {
        score += inst.d[sol[i]-1][sol[i+1]-1];
    }
    return score;
}

double compute_robust_score(const Instance& inst, const std::vector<int>& sol){
    double score = 0;
    std::vector<std::tuple<float,float>> knapsack_dij;
    for (unsigned int i=0; i<sol.size()-1; i++) {
        score += inst.d[sol[i]-1][sol[i+1]-1];
        float dij = inst.d[sol[i]-1][sol[i+1]-1];
        float Dij = inst.D[sol[i]-1][sol[i+1]-1];
        knapsack_dij.push_back(std::make_tuple(dij, Dij));
    }
    std::sort(knapsack_dij.begin(), knapsack_dij.end(), cmp_dij);
    float sum_knapsack = 0;
    for(auto it = knapsack_dij.begin(); it != knapsack_dij.end(); ++it){
        if (sum_knapsack + std::get<1>(*it) <= inst.d1){
            score += std::get<0>(*it)*std::get<1>(*it);
            sum_knapsack += std::get<1>(*it);
        }
        else{
            score += std::get<0>(*it)*(inst.d1-sum_knapsack);
            break;
        }
    }
    return score;
}

double compute_static_constraint(const Instance& inst, const std::vector<int>& sol){
    double score = 0;
    for (unsigned int i=0; i<sol.size(); i++) {
        score += inst.p[sol[i]-1];
    }
    return score;
}

double compute_robust_constraint(const Instance& inst, const std::vector<int>& sol){
    double score = 0;
    std::vector<float> knapsack_phi;
    for (unsigned int i=0; i<sol.size(); i++) {
        score += inst.p[sol[i]-1];
        float phi = inst.ph[sol[i]-1];
        knapsack_phi.push_back(phi);
    }
    std::sort(knapsack_phi.begin(), knapsack_phi.end(), std::greater<float>());
    float sum_knapsack = 0;
    for(auto it = knapsack_phi.begin(); it != knapsack_phi.end(); ++it){
        if (sum_knapsack + 2 <= inst.d2){
            score += 2*(*it);
            sum_knapsack += 2;
        }
        else{
            score += (*it)*(inst.d2-sum_knapsack);
            break;
        }
    }
    return score;
}


Solution_info::Solution_info(const Instance& inst, const std::vector<int>& solu, const double& K_, const std::chrono::steady_clock::time_point& start){
    sol = solu;
    K = K_;
    time_obtained = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - start);
    static_score = compute_static_score(inst,solu);
    robust_score = compute_robust_score(inst, solu);
    static_constraint = compute_static_constraint(inst,solu);
    robust_constraint = compute_robust_constraint(inst, solu);
}



Heuristic_instance::Heuristic_instance(const Instance& inst){
    inf_dist = std::vector<float>(inst.n, pow(10, 8));
    inf_dist_nodes = std::vector<float>(inst.n, pow(10, 8));
}


std::vector<float> Heuristic_instance::backward_dijkstra_distance(const Instance& inst) const{
    std::vector<float> dist(inst.n, pow(10, 8));
    std::vector<bool> visited(inst.n, false);
    std::vector<std::tuple<int, float>> to_visit;
    to_visit.emplace_back(inst.t-1, 0);
    std::make_heap(to_visit.begin(), to_visit.end(), cmp_index_val);
    dist[inst.t-1] = 0;
    visited[inst.t-1] = true;
    while (!to_visit.empty()) {
        pop_heap(to_visit.begin(), to_visit.end(), cmp_index_val);
        int i = std::get<0>(to_visit.back());
        to_visit.pop_back();
        for (unsigned int k = 0; k<inst.reverse_neighbors_list[i].size(); ++k) {
            int j = inst.reverse_neighbors_list[i][k];
            if (!visited[j]) {
                if (dist[j] == pow(10, 8)) {
                    dist[j] = dist[i] + inst.d[j][i];
                    to_visit.emplace_back(j, dist[j]);
                    std::push_heap(to_visit.begin(), to_visit.end(), cmp_index_val);
                }
                else if (dist[j] > dist[i] + inst.d[j][i]) {
                    dist[j] = dist[i] + inst.d[j][i];
                    // replace value of tuple of first value j
                    for (auto it = to_visit.begin(); it != to_visit.end(); ++it) {
                        if (std::get<0>(*it) == j) {
                            std::get<1>(*it) = dist[j];
                            std::push_heap(to_visit.begin(), it+1, cmp_index_val);
                            break;
                        }
                    }
                }
            }
        }
        visited[i] = true;
    }
    return dist;
}

std::vector<float> Heuristic_instance::backward_dijkstra_nodes(const Instance& inst) const{
    // renvoie le poids non robuste minimale  du chemin restant pour arriver a t
    // en particulier le poids restant ne prends pas en compte le poids du noeud courant
    std::vector<float> dist_nodes(inst.n, pow(10, 8));
    std::vector<bool> visited(inst.n, false);
    std::vector<std::tuple<int, float>> to_visit;
    to_visit.emplace_back(inst.t-1, 0);
    std::make_heap(to_visit.begin(), to_visit.end(), cmp_index_val);
    dist_nodes[inst.t-1] = 0;
    visited[inst.t-1] = true;
    while (!to_visit.empty()) {
        pop_heap(to_visit.begin(), to_visit.end(), cmp_index_val);
        int i = std::get<0>(to_visit.back());
        to_visit.pop_back();
        for (unsigned int k = 0; k<inst.reverse_neighbors_list[i].size(); ++k) {
            int j = inst.reverse_neighbors_list[i][k];
            if (!visited[j]) {
                if (dist_nodes[j] == pow(10, 8)) {
                    dist_nodes[j] = dist_nodes[i] + inst.p[i];
                    to_visit.emplace_back(j, dist_nodes[j]);
                    std::push_heap(to_visit.begin(), to_visit.end(), cmp_index_val);
                }
                else if (dist_nodes[j] > dist_nodes[i] + inst.p[i]) {
                    dist_nodes[j] = dist_nodes[i] + inst.p[i];
                    // replace value of tuple of first value j
                    for (auto it = to_visit.begin(); it != to_visit.end(); ++it) {
                        if (std::get<0>(*it) == j) {
                            std::get<1>(*it) = dist_nodes[j];
                            std::push_heap(to_visit.begin(), it+1, cmp_index_val);
                            break;
                        }
                    }
                }
            }
        }
        visited[i] = true;
    }
    return dist_nodes;
}

std::vector<int> Heuristic_instance::retrieve_feaible_sol(const Instance& inst) const{
    std::vector<float> dist_nodes(inst.n, pow(10, 8));
    std::vector<int> predecessors(inst.n, -1);
    std::vector<bool> visited(inst.n, false);
    std::vector<std::tuple<int, float>> to_visit;
    to_visit.emplace_back(inst.s-1, inst.p[inst.s-1] + 2*inst.ph[inst.s-1]);
    std::make_heap(to_visit.begin(), to_visit.end(), cmp_index_val);
    dist_nodes[inst.s-1] = inst.p[inst.s-1] + 2*inst.ph[inst.s-1];
    visited[inst.s-1] = true;
    while (!to_visit.empty()) {
        pop_heap(to_visit.begin(), to_visit.end(), cmp_index_val);
        int i = std::get<0>(to_visit.back());
        to_visit.pop_back();
        for (unsigned int k = 0; k<inst.neighbors_list[i].size(); ++k) {
            int j = inst.neighbors_list[i][k];
            if (!visited[j]) {
                if (dist_nodes[j] == pow(10, 8)) {
                    dist_nodes[j] = dist_nodes[i] + inst.p[i] + 2*inst.ph[i];
                    to_visit.emplace_back(j, dist_nodes[j]);
                    std::push_heap(to_visit.begin(), to_visit.end(), cmp_index_val);
                    predecessors[j]=i;
                }
                else if (dist_nodes[j] > dist_nodes[i] + inst.p[i] + 2*inst.ph[i]) {
                    dist_nodes[j] = dist_nodes[i] + inst.p[i] + 2*inst.ph[i];
                    predecessors[j]=i;
                    // replace value of tuple of first value j
                    for (auto it = to_visit.begin(); it != to_visit.end(); ++it) {
                        if (std::get<0>(*it) == j) {
                            std::get<1>(*it) = dist_nodes[j];
                            std::push_heap(to_visit.begin(), it+1, cmp_index_val);
                            break;
                        }
                    }
                }
            }
        }
        visited[i] = true;
    }
    std::vector<int> sol;
    int current_node = inst.t-1;
    sol.push_back(current_node+1);
    while(current_node != inst.s-1){
        current_node = predecessors[current_node];
        sol.push_back(current_node+1);
    }
    return sol;
}


std::vector<IloInt> Heuristic_instance::retrieve_feaible_sol_2(const Instance& inst, IloEnv env, const int& verbose) const{

    IloModel model(env);

    // Variables
    IloBoolVarArray x(env, inst.n_arc);
    IloBoolVarArray y(env, inst.n);
    IloNumVarArray beta(env, inst.n, 0.0, IloInfinity);
    IloNumVar alpha(env, 0.0, IloInfinity);

    // Objective
    IloObjective obj(env, IloScalProd(x, inst.d_vec), IloObjective::Minimize);
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


    model.add(inst.d2*alpha + 2*IloSum(beta) + IloScalProd(y, inst.p) <= inst.S);
    for (unsigned int i=0; i<inst.n; i++) {
        model.add(alpha + beta[i] >= inst.ph[i]*y[i]);
    }

    // Solve
    IloCplex cplex(model);
    cplex.setParam(IloCplex::Param::TimeLimit, 300);
    if (verbose < 2) cplex.setOut(env.getNullStream());

    std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
    cplex.solve();
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::chrono::microseconds duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);

    if (cplex.getStatus() == IloAlgorithm::Infeasible) {
        std::cout << inst.name << "," << "heuristics,,,,,,,,," << std::endl;
        throw std::domain_error("Infeasible heuristic retrieval model for instance " + inst.name);
    } else if (cplex.getStatus() == IloAlgorithm::Unknown) {
        std::cout << inst.name << "," << "heursitics,,,,,,,,," << std::endl;
        throw std::domain_error("No solution found for instance " + inst.name + ". Maybe not enough time");
    }

    std::vector<IloInt> sol;
    unsigned int current_node = inst.s-1;
    while (current_node != inst.t-1) {
        sol.push_back(current_node+1);
        for (unsigned int a=0; a<inst.n_arc; a++) {
            if (inst.mat[a].tail == current_node+1 && cplex.getValue(x[a]) == 1) {
                current_node = inst.mat[a].head-1;
                break;
            }
        }
        if (current_node == sol[sol.size()-1]-1) {
            std::cout << inst.name << "," << "retrieval,,,,,,,,," << std::endl;
            throw std::domain_error("Using arc that does not exist for instance " + inst.name);
        }
    }
    sol.push_back(inst.t);

    return sol;

}



std::vector<int> Heuristic_instance::astar_solve(const Instance& inst, const double& K, const int& verbose) const{
    std::vector<int> sol;
    std::vector<float>* dist = new std::vector<float>(inst.n, pow(10, 8));
    std::vector<float>* dist_star = new std::vector<float>(inst.n, pow(10, 8));
    std::vector<node_info>* node_infos = new std::vector<node_info>(inst.n);
    std::vector<bool> visited(inst.n, false);
    std::vector<std::tuple<int, float>> to_visit;
    (*node_infos)[inst.s-1] = node_info(inst.s-1, -1, 0, 0, inst.p[inst.s-1], 0, std::vector<std::tuple<float,float>>(), std::vector<float>());
    (*node_infos)[inst.s-1].redo_knapsack_phi(inst, inst.s-1);
    (*dist)[inst.s-1] = (*node_infos)[inst.s-1].compute_tot_dist(K);
    (*dist_star)[inst.s-1] = (*dist)[inst.s-1] + (inf_dist[inst.s-1] + K*inf_dist_nodes[inst.s-1]);
    to_visit.emplace_back(inst.s-1, (*dist_star)[inst.s-1]);
    std::make_heap(to_visit.begin(), to_visit.end(), cmp_index_val);
    visited[inst.s-1] = true;
    while (!to_visit.empty() && !visited[inst.t-1]) {
        pop_heap(to_visit.begin(), to_visit.end(), cmp_index_val);
        int i = std::get<0>(to_visit.back());
        to_visit.pop_back();
        for (unsigned int k = 0; k<inst.neighbors_list[i].size(); ++k) {
            int j = inst.neighbors_list[i][k];
            if (!visited[j]) {
                if ((*dist)[j] == pow(10, 8)) {
                    (*node_infos)[j] = node_info(j, i, (*node_infos)[i].dist+inst.d[i][j], (*node_infos)[i].robust_dist, (*node_infos)[i].dist_nodes+inst.p[j], (*node_infos)[i].robust_dist_nodes, (*node_infos)[i].knapsack_dij, (*node_infos)[i].knapsack_phi);
                    if(!(*node_infos)[j].is_knapsack_dij_full(inst.d1) || inst.d[i][j] > std::get<0>((*node_infos)[j].knapsack_dij.back())){
                        (*node_infos)[j].redo_knapsack_dij(inst, j);
                    }
                    if(!(*node_infos)[j].is_knapsack_phi_full(inst.d2) || inst.p[i] > (*node_infos)[j].knapsack_phi.back()){
                        (*node_infos)[j].redo_knapsack_phi(inst, j);
                    }

                    (*dist)[j] = (*node_infos)[j].compute_tot_dist(K);
                    (*dist_star)[j] = (*dist)[j] + (inf_dist[j] + K*inf_dist_nodes[j]);

                    to_visit.emplace_back(j, (*dist_star)[j]);
                    std::push_heap(to_visit.begin(), to_visit.end(), cmp_index_val);
                }
                else{
                    node_info new_node_info =  node_info(j, i, (*node_infos)[i].dist+inst.d[i][j], (*node_infos)[i].robust_dist, (*node_infos)[i].dist_nodes+inst.p[i], (*node_infos)[i].robust_dist_nodes, (*node_infos)[i].knapsack_dij, (*node_infos)[i].knapsack_phi);
                    if(!(*node_infos)[j].is_knapsack_dij_full(inst.d1) || inst.d[i][j] > std::get<0>((*node_infos)[j].knapsack_dij.back())){
                        (*node_infos)[j].redo_knapsack_dij(inst, j);
                    }
                    if(!(*node_infos)[j].is_knapsack_phi_full(inst.d2) || inst.p[i] > (*node_infos)[j].knapsack_phi.back()){
                        (*node_infos)[j].redo_knapsack_phi(inst, j);
                    }
                    float new_dist = new_node_info.compute_tot_dist(K);
                    if (new_dist < (*dist)[j]){
                        (*node_infos)[j] = new_node_info;
                        (*dist)[j] = new_dist;
                        (*dist_star)[j] = (*dist)[j] + (inf_dist[j] + K*inf_dist_nodes[j]);
                        // replace value of tuple of first value j
                        for (auto it = to_visit.begin(); it != to_visit.end(); ++it) {
                            if (std::get<0>(*it) == j) {
                                std::get<1>(*it) = (*dist_star)[j];
                                std::push_heap(to_visit.begin(), it+1, cmp_index_val);
                                break;
                            }
                        }
                    }
                }
            }
        }
        visited[i] = true;
    }

    int current_node = inst.t-1;
    sol.push_back(current_node+1);
    while(current_node != inst.s-1){
        current_node = (*node_infos)[current_node].parent;
        sol.push_back(current_node+1);
    }
    std::reverse(sol.begin(), sol.end());

    return sol;
}

void Heuristic_instance::complete_astar_solve(const Instance& inst, IloEnv env, const double& precision_K, const int& max_iter, const float& max_duration, const int& verbose) const{

    std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();

    std::vector<int> sol_0 = astar_solve(inst, 0., verbose);
    Solution_info sol_inf(inst, sol_0, 0, start);

    std::vector<int> sol_1 = astar_solve(inst, 1., verbose);
    Solution_info sol_sup(inst, sol_1, 1., start);

    double sup_K = 1.;
    double inf_K = 0.;
    int counter_iter_init = 0;

    while(sol_sup.robust_constraint > inst.S && counter_iter_init < 20){
        sup_K *= 2;
        std::vector<int> sol_K = astar_solve(inst, sup_K, verbose);
        Solution_info new_sol_sup(inst, sol_K, sup_K, start);
        sol_sup = new_sol_sup;
        counter_iter_init += 1;
    }



    if(counter_iter_init == 20){
        if (verbose>0){
            std::cout << "initialisation of K took too much iter" << counter_iter_init << "/" << 20 << std::endl;
        }

        std::vector<int> retrieved_feasible_sol = retrieve_feaible_sol(inst);

        if(compute_robust_constraint(inst, retrieved_feasible_sol) > inst.S){
            retrieved_feasible_sol.clear();
            std::vector<IloInt> retrieved_feasible_sol_2 = retrieve_feaible_sol_2(inst, env, verbose);
            for (unsigned int i=0; i<retrieved_feasible_sol_2.size(); i++) {
                retrieved_feasible_sol.push_back(retrieved_feasible_sol_2[i]);
            }
        }
        Solution_info retrieved_sol(inst, retrieved_feasible_sol, -1, start);

        std::string path_str = "[";
        for (auto it = retrieved_sol.sol.begin(); it != retrieved_sol.sol.end(); ++it) {
            path_str += std::to_string(*it) + ";";
        }
        path_str += std::to_string(inst.t) + "]";
        std::chrono::microseconds duration = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - start);

        std::cout << inst.name << ","
        << "heuristics,"
        << retrieved_sol.robust_score << ","
        << std::to_string(0) << ","
        << static_cast<double>(duration.count()) / 1e6 << ","
        << std::to_string(0) << ","
        << retrieved_sol.robust_constraint << ","
        << retrieved_sol.static_score<< ","
        << retrieved_sol.static_constraint << ","
        << inst.S << ","
        << path_str << std::endl;
    }

    else{
        std::chrono::microseconds duration = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - start);
        int iter = 0;

        while ((static_cast<double>(duration.count()) / 1e6) < max_duration && sup_K-inf_K > precision_K && iter<max_iter){
            if(sol_inf.robust_score == sol_sup.robust_score){break;}
            if(verbose>1){
                cout << "iter = " << iter << ", duration in dichotomie = " << (static_cast<double>(duration.count()) / 1e6) << ", inf_K = " << std::to_string(inf_K) << ", sup_K =" << std::to_string(sup_K) << endl;
            }
            double new_K = (sup_K+inf_K)/2;
            std::vector<int> new_solu_K = astar_solve(inst, new_K, verbose);
            Solution_info new_sol(inst, new_solu_K, new_K,start);
            iter+=1;
            if(new_sol.robust_constraint > inst.S){
                sol_inf = new_sol;
                inf_K = new_K;
            }
            else{
                sol_sup = new_sol;
                sup_K = new_K;
            }
            duration = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - start);
        }

        if(verbose>0){
            std::cout << inst.name << std::endl;
            std::cout << "borne inf = " << sol_inf.robust_score << ", weight = " << sol_inf.robust_constraint << ", for S = " << inst.S << ", found in "<< static_cast<double>(sol_inf.time_obtained.count()) / 1e6 << ", for K = "<< std::to_string(inf_K) << std::endl;
            std::cout << "with path : ";
            for (auto it = sol_inf.sol.begin(); it != sol_inf.sol.end(); ++it) {
                std::cout << *it << " ";
            }
            std::cout << " " << std::endl;
            std::cout << "borne sup = " << sol_sup.robust_score << ", weight = " << sol_sup.robust_constraint << ", for S = " << inst.S << ", found in "<< static_cast<double>(sol_sup.time_obtained.count()) / 1e6 << ", for K = "<< std::to_string(sup_K) <<std::endl;
            std::cout << "with path : ";
            for (auto it = sol_sup.sol.begin(); it != sol_sup.sol.end(); ++it) {
                std::cout << *it << " ";
            }
            std::cout << " " << std::endl;
        }

        if (verbose >0){
            if ((static_cast<double>(duration.count()) / 1e6) > max_duration){
                std::cout << "ended because of duration after " << (static_cast<double>(duration.count()) / 1e6) << " s " << std::endl;
            }
            if (iter == max_iter){
                std::cout << "ended because of iter number " << iter << "/" << max_iter << std::endl;
            }
            if(sup_K-inf_K <= precision_K){
                std::cout << "ended because of K gap with inf_K = " << inf_K << ", sup_K = " << sup_K << ", and precision_K = "<< precision_K << std::endl;
            }
        }

        if(verbose>1){
            std::vector<int> sol = astar_solve(inst, inf_K, verbose);
            std::cout << "Solution  for K = " << inf_K << std::endl;
            for (auto it = sol.begin(); it != sol.end(); ++it) {
                std::cout << *it << " ";
            }
            cout << " " << endl;
            std::cout << "static objective = " << compute_static_score(inst, sol) << std::endl;
            std::cout << "robust objective = " << compute_robust_score(inst, sol) << std::endl;
            std::cout << "static constraint = " << compute_static_constraint(inst, sol) << std::endl;
            std::cout << "robust constraint  = " << compute_robust_constraint(inst, sol) << " for S = " << inst.S << std::endl;

            sol = astar_solve(inst, sup_K, verbose);
            std::cout << "Solution  for K = " << sup_K << std::endl;
            for (auto it = sol.begin(); it != sol.end(); ++it) {
                std::cout << *it << " ";
            }
            cout << " " << endl;
            std::cout << "static objective = " << compute_static_score(inst, sol) << std::endl;
            std::cout << "robust objective = " << compute_robust_score(inst, sol) << std::endl;
            std::cout << "static constraint = " << compute_static_constraint(inst, sol) << std::endl;
            std::cout << "robust constraint  = " << compute_robust_constraint(inst, sol) << " for S = " << inst.S << std::endl;
        }

        std::string path_str = "[";
        for (auto it = sol_sup.sol.begin(); it != sol_sup.sol.end(); ++it) {
            path_str += std::to_string(*it) + ";";
        }
        path_str += std::to_string(inst.t) + "]";

        std::cout << inst.name << ","
        << "heuristics,"
        << sol_sup.robust_score << ","
        << sol_inf.robust_score << ","
        << static_cast<double>(duration.count()) / 1e6 << ","
        << std::to_string(0) << ","
        << sol_sup.robust_constraint << ","
        << sol_sup.static_score<< ","
        << sol_sup.static_constraint << ","
        << inst.S << ","
        << path_str << std::endl;
    }
}

