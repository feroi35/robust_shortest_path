#include "parser.h"
#include "heuristics.h"


// #include <limits>

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

std::vector<int> Heuristic_instance::astar_solve(const Instance& inst, const float& K, const int& verbose) const{
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

