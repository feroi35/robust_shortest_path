#include "parser.h"
#include "heuristics.h"


// #include <limits>

bool cmp_index_val(const std::tuple<int, double>& a, const std::tuple<int, double>& b) {
    return std::get<1>(a) > std::get<1>(b);
}

bool cmp_dij(const std::tuple<double, double>& a, const std::tuple<double, double>& b) {
    return std::get<0>(a) > std::get<0>(b);
}

double node_info::sum_knapsack_capa() const{
    double sum = 0;
    for(auto it = knapsack_dij.begin(); it != knapsack_dij.end(); ++it){
        sum += std::get<1>(*it);
    }
    return sum;
}

void node_info::compute_dist_robust(const double d1){
    double sum = 0;
    double weight = 0;
    for(auto it = knapsack_dij.begin(); it != knapsack_dij.end(); ++it){
        if (weight + std::get<1>(*it) <= d1){
            sum += std::get<0>(*it)*(1+std::get<1>(*it));
            weight += std::get<1>(*it);
        }
        else{
            sum += std::get<0>(*it)*(1+ (d1-weight));
            break;
        }
    }
    robust_dist = sum;
}

void node_info::compute_nodes_robust(const double d2){
    double sum = 0;
    double weight = 0;
    for(auto it = knapsack_phi.begin(); it != knapsack_phi.end(); ++it){
        if (weight + 2 <= d2){
            sum += *it;
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
    
    double new_dij = inst.d[parent][new_node];
    double new_Dij = inst.D[parent][new_node];

    knapsack_dij.push_back(std::make_tuple(new_dij, new_Dij));
    std::sort(knapsack_dij.begin(), knapsack_dij.end(), cmp_dij);

    double sum_knapsack = sum_knapsack_capa();

    while(sum_knapsack-std::get<1>(knapsack_dij.back()) >= inst.d1){
        std::tuple<double,double> last = knapsack_dij.back();
        knapsack_dij.pop_back();
        sum_knapsack -= std::get<1>(last);
    }

    compute_dist_robust(inst.d1);
}

void node_info::redo_knapsack_phi(const Instance& inst, const int& new_node){
    double new_phi = inst.ph[new_node];

    knapsack_phi.push_back(new_phi);
    std::sort(knapsack_phi.begin(), knapsack_phi.end());

    double sum_knapsack = 2*knapsack_phi.size();

    while(sum_knapsack-2 >= inst.d2){
        knapsack_phi.pop_back();
        sum_knapsack -= 2;
    }

    compute_nodes_robust(inst.d2);
}


std::vector<double>* backward_dijkstra_distance(const Instance& inst) {
    std::vector<double>* dist = new std::vector<double>(inst.n, pow(10, 8));
    std::vector<bool> visited(inst.n, false);
    std::vector<std::tuple<int, double>> to_visit;
    to_visit.emplace_back(inst.t-1, 0);
    std::make_heap(to_visit.begin(), to_visit.end(), cmp_index_val);
    (*dist)[inst.t-1] = 0;
    visited[inst.t-1] = true;
    while (!to_visit.empty()) {
        pop_heap(to_visit.begin(), to_visit.end(), cmp_index_val);
        int i = std::get<0>(to_visit.back());
        to_visit.pop_back();
        for (int j : (*inst.reverse_neighbors_list)[i]) {
            if (!visited[j]) {
                if ((*dist)[j] == pow(10, 8)) {
                    (*dist)[j] = (*dist)[i] + inst.d[j][i];
                    to_visit.emplace_back(j, (*dist)[j]);
                    std::push_heap(to_visit.begin(), to_visit.end(), cmp_index_val);
                }
                else if ((*dist)[j] > (*dist)[i] + inst.d[j][i]) {
                    (*dist)[j] = (*dist)[i] + inst.d[j][i];
                    // replace value of tuple of first value j
                    for (auto it = to_visit.begin(); it != to_visit.end(); ++it) {
                        if (std::get<0>(*it) == j) {
                            std::get<1>(*it) = (*dist)[j];
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

std::vector<double>* backward_dijkstra_nodes(const Instance& inst){
    // renvoie le poids non robuste minimale  du chemin restant pour arriver a t
    // en particulier le poids restant ne prends pas en compte le poids du noeud courant
    std::vector<double>* dist_nodes = new std::vector<double>(inst.n, pow(10, 8));
    std::vector<bool> visited(inst.n, false);
    std::vector<std::tuple<int, double>> to_visit;
    to_visit.emplace_back(inst.t-1, 0);
    std::make_heap(to_visit.begin(), to_visit.end(), cmp_index_val);
    (*dist_nodes)[inst.t-1] = 0;
    visited[inst.t-1] = true;
    while (!to_visit.empty()) {
        pop_heap(to_visit.begin(), to_visit.end(), cmp_index_val);
        int i = std::get<0>(to_visit.back());
        to_visit.pop_back();
        for (int j : (*inst.reverse_neighbors_list)[i]) {
            if (!visited[j]) {
                if ((*dist_nodes)[j] == pow(10, 8)) {
                    (*dist_nodes)[j] = (*dist_nodes)[i] + inst.p[i];
                    to_visit.emplace_back(j, (*dist_nodes)[j]);
                    std::push_heap(to_visit.begin(), to_visit.end(), cmp_index_val);
                }
                else if ((*dist_nodes)[j] > (*dist_nodes)[i] + inst.p[i]) {
                    (*dist_nodes)[j] = (*dist_nodes)[i] + inst.p[i];
                    // replace value of tuple of first value j
                    for (auto it = to_visit.begin(); it != to_visit.end(); ++it) {
                        if (std::get<0>(*it) == j) {
                            std::get<1>(*it) = (*dist_nodes)[j];
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

std::vector<int> astar_solve(const Instance& inst, const std::vector<double>& inf_dist, const std::vector<double>& inf_dist_nodes, const int& verbose) {
    std::vector<int> sol;
    // std::vector<double>* dist = new std::vector<double>(inst.n, pow(10, 8));
    // std::vector<double>* dist_star = new std::vector<double>(inst.n, pow(10, 8));
    // std::vector<bool> visited(inst.n, false);
    // std::vector<std::tuple<int, double>> to_visit;
    // to_visit.emplace_back(inst.s-1, 0);
    // std::make_heap(to_visit.begin(), to_visit.end(), cmp_index_val);
    // visited[inst.s-1] = true;
    // while (!to_visit.empty()) {
    //     pop_heap(to_visit.begin(), to_visit.end(), cmp_index_val);
    //     int i = std::get<0>(to_visit.back());
    //     to_visit.pop_back();
    //     sol.push_back(i+1);
    //     if (i == inst.t-1) {
    //         break;
    //     }
    //     for (int j : (*inst.neighbors_list)[i]) {
    //         if (!visited[j]) {
    //             if (inf_dist[j] == pow(10, 8)) {
    //                 inf_dist[j] = inf_dist[i] + inst.d[i][j];
    //                 to_visit.emplace_back(j, inf_dist[j] + inf_dist_nodes[j]);
    //                 std::push_heap(to_visit.begin(), to_visit.end(), cmp_index_val);
    //             }
    //             else if (inf_dist[j] > inf_dist[i] + inst.d[i][j]) {
    //                 inf_dist[j] = inf_dist[i] + inst.d[i][j];
    //                 // replace value of tuple of first value j
    //                 for (auto it = to_visit.begin(); it != to_visit.end(); ++it) {
    //                     if (std::get<0>(*it) == j) {
    //                         std::get<1>(*it) = inf_dist[j] + inf_dist_nodes[j];
    //                         std::push_heap(to_visit.begin(), it+1, cmp_index_val);
    //                         break;
    //                     }
    //                 }
    //             }
    //         }
    //     }
    //     visited[i] = true;
    // }
    return sol;
}