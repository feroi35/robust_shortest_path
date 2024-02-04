// heuristics.h
#ifndef HEURISTICS_H
#define HEURISTICS_H

#include "solve_method.h"
#include <vector>
#include <tuple>
#include <algorithm>
#include <numeric>


class Instance;
class IloEnv;

struct HeuristicMethod : public SolveMethod {
    HeuristicMethod(const double& precision_K_ = 1e-6, const int& max_iter_ = 2000,
            const float& max_duration_ = 20.0): precision_K(precision_K_), max_iter(max_iter_), max_duration(max_duration_) {
        method_name = "heuristic";
    };
    double precision_K;
    int max_iter;
    float max_duration;
    std::vector<float> inf_dist;
    std::vector<float> inf_dist_nodes;

    std::vector<float> backward_dijkstra_distance(const Instance& inst) const;
    std::vector<float> backward_dijkstra_nodes(const Instance& inst) const;
    std::vector<IloInt> retrieve_feasible_sol(const Instance& inst) const;
    std::vector<IloInt> retrieve_feasible_sol_2(const Instance& inst, IloEnv& env, const int& verbose);
    std::vector<IloInt> astar_solve(const Instance& inst, const double& K, const int& verbose=0) const;
    void complete_astar_solve(Instance& inst, IloEnv& env, const double& precision_K, const int& max_iter, const float& max_duration, const int& verbose=0);

    void solve(IloEnv& env, Instance& inst, const unsigned int& time_limit, const int& verbose) override;
};


struct NodesInfo{
    int index;
    int parent;
    float dist;
    float robust_dist;
    float dist_nodes;
    float robust_dist_nodes;
    std::vector<std::tuple<float,float>> knapsack_dij; // (dij, Dij), trié du plus grand dij au plus petit
    std::vector<float> knapsack_phi;

    NodesInfo(): index(-1), parent(-1), dist(0), robust_dist(0), dist_nodes(0), robust_dist_nodes(0),
        knapsack_dij(), knapsack_phi() {};
    NodesInfo(int index, int parent, float dist, float robust_dist, float dist_nodes, float robust_dist_nodes,
        std::vector<std::tuple<float,float>> knapsack_dij, std::vector<float> knapsack_phi): index(index),
        parent(parent), dist(dist), robust_dist(robust_dist), dist_nodes(dist_nodes), robust_dist_nodes(robust_dist_nodes),
        knapsack_dij(knapsack_dij), knapsack_phi(knapsack_phi) {};

    float sum_knapsack_capa() const;
    void compute_dist_robust(const float d1);
    void compute_nodes_robust(const float d2);

    bool is_knapsack_dij_full(const float& d1) const { return sum_knapsack_capa() >= d1; };
    bool is_knapsack_phi_full(const float& d2) const { return 2*knapsack_phi.size() >= d2; };

    float compute_tot_dist(const double K) const { return dist + robust_dist + K*(dist_nodes + robust_dist_nodes); };

    void redo_knapsack_dij(const Instance& inst, const int& new_node);
    void redo_knapsack_phi(const Instance& inst, const int& new_node);
};


struct SolutionInfo{
    std::vector<IloInt> sol;
    double K;
    float static_score;
    float robust_score;
    float static_constraint;
    float robust_constraint;
    std::chrono::microseconds time_obtained;

    SolutionInfo(const Instance& inst, const std::vector<IloInt>& solu, const double& K_, const std::chrono::steady_clock::time_point& start);
};

#endif