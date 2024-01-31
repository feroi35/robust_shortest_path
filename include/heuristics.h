// heuristics.h
#ifndef HEURISTICS_H
#define HEURISTICS_H

#include <vector>
#include <tuple>
#include <algorithm>
#include <numeric>

class Instance; // forward declaration
class IloEnv;

bool cmp_index_val(const std::tuple<int, float>& a, const std::tuple<int, float>& b);
bool cmp_dij(const std::tuple<float, float>& a, const std::tuple<float, float>& b);

class node_info{
    public:
        int index;
        int parent;
        float dist;
        float robust_dist;
        float dist_nodes;
        float robust_dist_nodes;
        std::vector<std::tuple<float,float>> knapsack_dij; // (dij, Dij), trié du plus grand dij au plus petit
        std::vector<float> knapsack_phi;

        node_info(): index(-1), parent(-1), dist(0), robust_dist(0), dist_nodes(0), robust_dist_nodes(0), knapsack_dij(), knapsack_phi() {};
        node_info(int index, int parent, float dist, float robust_dist, float dist_nodes, float robust_dist_nodes, std::vector<std::tuple<float,float>> knapsack_dij, std::vector<float> knapsack_phi): index(index), parent(parent), dist(dist), robust_dist(robust_dist), dist_nodes(dist_nodes), robust_dist_nodes(robust_dist_nodes), knapsack_dij(knapsack_dij), knapsack_phi(knapsack_phi) {};

        float sum_knapsack_capa() const;
        void compute_dist_robust(const float d1);
        void compute_nodes_robust(const float d2);

        bool is_knapsack_dij_full(const float& d1) const { return sum_knapsack_capa() >= d1; };
        bool is_knapsack_phi_full(const float& d2) const { return 2*knapsack_phi.size() >= d2; };

        float compute_tot_dist(const double K) const { return dist + robust_dist + K*dist_nodes + K*robust_dist_nodes; };

        void redo_knapsack_dij(const Instance& inst, const int& new_node);
        void redo_knapsack_phi(const Instance& inst, const int& new_node);
};


class Solution_info{
    public:
        std::vector<int> sol;
        double K;
        float static_score;
        float robust_score;
        float static_constraint;
        float robust_constraint;
        std::chrono::microseconds time_obtained;

        Solution_info(const Instance& inst, const std::vector<int>& solu, const double& K_, const std::chrono::steady_clock::time_point& start);
        ~Solution_info(){};
};


double compute_static_score(const Instance& inst, const std::vector<int>& sol);
double compute_robust_score(const Instance& inst, const std::vector<int>& sol);
double compute_static_constraint(const Instance& inst, const std::vector<int>& sol);
double compute_robust_constraint(const Instance& inst, const std::vector<int>& sol);


class Heuristic_instance{
    public:
        std::vector<float> inf_dist;
        std::vector<float> inf_dist_nodes;

        // Heuristic_instance(): inst(), inf_dist(), inf_dist_nodes() {};
        Heuristic_instance(const Instance& inst);
        ~Heuristic_instance(){};

        std::vector<float> backward_dijkstra_distance(const Instance& inst) const;
        std::vector<float> backward_dijkstra_nodes(const Instance& inst) const;

        std::vector<int> retrieve_feaible_sol(const Instance& inst) const;
        std::vector<IloInt> retrieve_feaible_sol_2(const Instance& inst, IloEnv env, const int& verbose) const;

        std::vector<int> astar_solve(const Instance& inst, const double& K, const int& verbose=0) const;
        void complete_astar_solve(Instance& inst,IloEnv env, const double& precision_K, const int& max_iter, const float& max_duration, const int& verbose=0) const;
};

#endif