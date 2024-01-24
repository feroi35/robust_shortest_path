// heuristics.h
#ifndef HEURISTICS_H
#define HEURISTICS_H

#include <vector>
#include <tuple>
#include <algorithm>
#include <numeric>

class Instance; // forward declaration

bool cmp_index_val(const std::tuple<int, double>& a, const std::tuple<int, double>& b);

bool cmp_dij(const std::tuple<double, double>& a, const std::tuple<double, double>& b);

class node_info{
    public:
        int index;
        int parent;   
        double dist;
        double robust_dist;
        double dist_nodes;
        double robust_dist_nodes;
        std::vector<std::tuple<double,double>> knapsack_dij; // (dij, Dij), trié du plus grand dij au plus petit
        std::vector<double> knapsack_phi;

        node_info(int index, int parent, double dist, double robust_dist, double dist_nodes, double robust_dist_nodes, std::vector<std::tuple<double,double>> knapsack_dij, std::vector<double> knapsack_phi): index(index), parent(parent), dist(dist), robust_dist(robust_dist), dist_nodes(dist_nodes), robust_dist_nodes(robust_dist_nodes), knapsack_dij(knapsack_dij), knapsack_phi(knapsack_phi) {};
        
        double sum_knapsack_capa() const;
        void compute_dist_robust(const double d1);
        void compute_nodes_robust(const double d2);
        
        void redo_knapsack_dij(const Instance& inst, const int& new_node);
        void redo_knapsack_phi(const Instance& inst, const int& new_node);
};


std::vector<double>* backward_dijkstra_distance(const Instance& inst);
std::vector<double>* backward_dijkstra_nodes(const Instance& inst);

std::vector<int> astar_solve(const Instance& inst, const std::vector<double>& inf_dist, const std::vector<double>& inf_dist_nodes, const int& verbose=0);

#endif