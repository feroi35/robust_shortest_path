#include "parser.h"
#include "heuristics.h"


bool cmp_index_val(const std::tuple<int, float>& a, const std::tuple<int, float>& b) {return std::get<1>(a) > std::get<1>(b);}
bool cmp_dij(const std::tuple<float, float>& a, const std::tuple<float, float>& b) {return std::get<0>(a) > std::get<0>(b);}


void HeuristicMethod::solve(IloEnv& env, Instance& inst, const unsigned int& time_limit, const int& verbose) {
    inf_dist = backward_dijkstra_distance(inst);
    inf_dist_nodes = backward_dijkstra_nodes(inst);
    if (verbose>0) std::cout <<"Heuristic initialization done" << std::endl;
    complete_astar_solve(inst, env, precision_K, max_iter, max_duration, verbose);
}


std::vector<float> HeuristicMethod::backward_dijkstra_distance(const Instance& inst) const {
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


std::vector<float> HeuristicMethod::backward_dijkstra_nodes(const Instance& inst) const{
    // renvoie le poids non robuste minimal du chemin restant pour arriver a t
    // en particulier le poids restant ne prend pas en compte le poids du noeud courant
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


std::vector<IloInt> HeuristicMethod::retrieve_feasible_sol(const Instance& inst) const{
    // Qu'est ce que ça fait ?? 
    // Utilité de garder cette méthode ?
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
    std::vector<IloInt> sol;
    int current_node = inst.t-1;
    sol.push_back(current_node+1);
    while(current_node != inst.s-1) {
        current_node = predecessors[current_node];
        sol.push_back(current_node+1);
    }
    std::reverse(sol.begin(), sol.end());
    return sol;
}


std::vector<IloInt> HeuristicMethod::retrieve_feasible_sol_2(const Instance& inst, IloEnv& env, const int& verbose) {
    // Looks a lot like dualized method but with a relaxed objective
    // There is no more proof of optimality but it gives a feasible solution
    IloModel model(env);

    // Variables
    IloBoolVarArray x(env, inst.n_arc);
    IloBoolVarArray y(env, inst.n);
    IloNumVarArray beta(env, inst.n, 0.0, IloInfinity);
    IloNumVar alpha(env, 0.0, IloInfinity);
    IloNumVar z(env, 0.0, IloInfinity);

    // Objective
    IloObjective obj(env, z, IloObjective::Minimize);
    model.add(obj);

    // Constraints
    add_static_constraints(env, model, x, y, z, inst);
    model.add(inst.d2*alpha + 2*IloSum(beta) + IloScalProd(y, inst.p) <= inst.S);
    for (unsigned int i=0; i<inst.n; i++) {
        model.add(alpha + beta[i] >= inst.ph[i]*y[i]);
    }

    // Solve
    IloCplex cplex(model);
    parametrizeCplex(cplex, 300, verbose);
    cplex.solve();

    if (cplex.getStatus() == IloAlgorithm::Infeasible) {
        throw std::domain_error("Infeasible heuristic retrieval model for instance " + inst.name);
    } else if (cplex.getStatus() == IloAlgorithm::Unknown) {
        throw std::domain_error("No solution found for heuristic retrieval model for instance "
            + inst.name + ". Maybe not enough time");
    }

    std::vector<IloInt> sol;
    unsigned int current_node = inst.s-1;

    IloNumArray xValues(env);
    cplex.getValues(xValues, x);
    while (current_node != inst.t-1) {
        sol.push_back(current_node+1);
        for (unsigned int a=0; a<inst.n_arc; a++) {
            if (inst.mat[a].tail == current_node+1 && xValues[a] >= 1 - TOL) {
                current_node = inst.mat[a].head-1;
                break;
            }
        }
        if (current_node == sol[sol.size()-1]-1) {
            throw std::domain_error("Using arc that does not exist for instance " + inst.name);
        }
    }
    sol.push_back(inst.t);
    infBound = cplex.getObjValue();
    return sol;
}


std::vector<IloInt> HeuristicMethod::astar_solve(const Instance& inst, const double& K, const int& verbose) const {
    std::vector<IloInt> sol;
    std::vector<bool> visited(inst.n, false);
    std::vector<std::tuple<int, float>> to_visit;
    std::vector<float>* dist = new std::vector<float>(inst.n, pow(10, 8));
    std::vector<float>* dist_star = new std::vector<float>(inst.n, pow(10, 8));
    std::vector<NodesInfo>* nodesInfos = new std::vector<NodesInfo>(inst.n);

    (*nodesInfos)[inst.s-1] = NodesInfo(inst.s-1, -1, 0, 0, inst.p[inst.s-1], 0,
        std::vector<std::tuple<float,float>>(),std::vector<float>());
    (*nodesInfos)[inst.s-1].redo_knapsack_phi(inst, inst.s-1);
    (*dist)[inst.s-1] = (*nodesInfos)[inst.s-1].compute_tot_dist(K);
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
                    (*nodesInfos)[j] = NodesInfo(j, i, (*nodesInfos)[i].dist+inst.d[i][j], (*nodesInfos)[i].robust_dist,
                                                (*nodesInfos)[i].dist_nodes+inst.p[j], (*nodesInfos)[i].robust_dist_nodes,
                                                (*nodesInfos)[i].knapsack_dij, (*nodesInfos)[i].knapsack_phi);
                    if (!(*nodesInfos)[j].is_knapsack_dij_full(inst.d1) ||
                            inst.d[i][j] > std::get<0>((*nodesInfos)[j].knapsack_dij.back())) {
                        (*nodesInfos)[j].redo_knapsack_dij(inst, j);
                    }
                    if(!(*nodesInfos)[j].is_knapsack_phi_full(inst.d2) || inst.p[i] > (*nodesInfos)[j].knapsack_phi.back()) {
                        (*nodesInfos)[j].redo_knapsack_phi(inst, j);
                    }

                    (*dist)[j] = (*nodesInfos)[j].compute_tot_dist(K);
                    (*dist_star)[j] = (*dist)[j] + (inf_dist[j] + K*inf_dist_nodes[j]);

                    to_visit.emplace_back(j, (*dist_star)[j]);
                    std::push_heap(to_visit.begin(), to_visit.end(), cmp_index_val);
                }
                else {
                    NodesInfo new_NodesInfo =  NodesInfo(j, i, (*nodesInfos)[i].dist+inst.d[i][j],
                                                        (*nodesInfos)[i].robust_dist, (*nodesInfos)[i].dist_nodes+inst.p[i],
                                                        (*nodesInfos)[i].robust_dist_nodes, (*nodesInfos)[i].knapsack_dij,
                                                        (*nodesInfos)[i].knapsack_phi);
                    if (!(*nodesInfos)[j].is_knapsack_dij_full(inst.d1) ||
                            inst.d[i][j] > std::get<0>((*nodesInfos)[j].knapsack_dij.back())) {
                        (*nodesInfos)[j].redo_knapsack_dij(inst, j);
                    }
                    if (!(*nodesInfos)[j].is_knapsack_phi_full(inst.d2) || inst.p[i] > (*nodesInfos)[j].knapsack_phi.back()) {
                        (*nodesInfos)[j].redo_knapsack_phi(inst, j);
                    }
                    float new_dist = new_NodesInfo.compute_tot_dist(K);
                    if (new_dist < (*dist)[j]) {
                        (*nodesInfos)[j] = new_NodesInfo;
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
    while(current_node != inst.s-1) {
        current_node = (*nodesInfos)[current_node].parent;
        sol.push_back(current_node+1);
    }
    std::reverse(sol.begin(), sol.end());

    // Free memory
    delete dist;
    delete dist_star;
    delete nodesInfos;
    return sol;
}


void HeuristicMethod::complete_astar_solve(Instance& inst, IloEnv& env, const double& precision_K,
        const int& max_iter, const float& max_duration, const int& verbose) {
    std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();

    // Looking for a feasible solution to start the dichotomy
    double inf_K = 0.;
    double sup_K = 1.;
    std::vector<IloInt> sol_inf_ = astar_solve(inst, inf_K, verbose);
    std::vector<IloInt> sol_sup_ = astar_solve(inst, sup_K, verbose);
    SolutionInfo sol_inf(inst, sol_inf_, inf_K, start);
    SolutionInfo sol_sup(inst, sol_sup_, sup_K, start);
    int counter_iter_init = 0;
    while (sol_sup.robust_constraint > inst.S && counter_iter_init < 20) {
        // Increasing the weight of the robust constraint until a feasible solution is found
        sup_K *= 2;
        std::vector<IloInt> sol_K = astar_solve(inst, sup_K, verbose);
        SolutionInfo new_sol_sup(inst, sol_K, sup_K, start);
        sol_sup = new_sol_sup;
        counter_iter_init += 1;
    }

    if (counter_iter_init == 20) {
        // The solution is not feasible
        if (verbose > 0) {
            std::cout << "Initialization of K took too much iterations" << std::endl;
            std::cout << "Trying first other method to get admissible solution" << std::endl;
        }
        std::vector<IloInt> retrieved_feasible_sol = retrieve_feasible_sol(inst);
        if (inst.compute_robust_constraint_knapsack(retrieved_feasible_sol) > inst.S) {
            // The solution is still not feasible
            if (verbose > 0) {
                std::cout << "Trying second other method to get admissible solution" << std::endl;
            }
            retrieved_feasible_sol.clear();
            std::vector<IloInt> retrieved_feasible_sol_2 = retrieve_feasible_sol_2(inst, env, verbose);
            for (unsigned int i=0; i<retrieved_feasible_sol_2.size(); i++) {
                retrieved_feasible_sol.push_back(retrieved_feasible_sol_2[i]);
            }
        }

        SolutionInfo retrieved_sol(inst, retrieved_feasible_sol, -1, start);
        if (!inst.sol.empty()) {
            std::cerr << "Warning: solution vector not empty for instance " << inst.name << std::endl;
            inst.sol.clear();
        }
        for (auto it = retrieved_sol.sol.begin(); it != retrieved_sol.sol.end(); ++it) {
            inst.sol.push_back(*it);
        }
        return;
    }

    // A feasible solution has been found, now we can start the dichotomy
    std::chrono::microseconds duration =  std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - start);
    int iter = 0;
    while ((static_cast<double>(duration.count()) / 1e6) < max_duration && sup_K-inf_K > precision_K
            && iter<max_iter && sol_inf.robust_score!=sol_sup.robust_score) {
        if (sol_inf.robust_score == sol_sup.robust_score)  break;
        if (verbose>1) {
            std::cout << "iter = " << iter << ", Time spent in dichotomy = " << (static_cast<double>(duration.count()) / 1e6)
            << ", inf_K = " << std::to_string(inf_K) << ", sup_K =" << std::to_string(sup_K) << endl;
        }

        double new_K = (sup_K+inf_K)/2;
        std::vector<IloInt> new_solu_K = astar_solve(inst, new_K, verbose);
        SolutionInfo new_sol(inst, new_solu_K, new_K,start);
        iter+=1;
        if (new_sol.robust_constraint > inst.S) {
            sol_inf = new_sol;
            inf_K = new_K;
        } else {
            sol_sup = new_sol;
            sup_K = new_K;
        }
        duration = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - start);
    }

    if (verbose>0) {
        std::cout << inst.name << std::endl;
        std::cout << "borne inf = " << sol_inf.robust_score << ", weight = " << sol_inf.robust_constraint << ", for S = " << inst.S << ", found in "<< static_cast<double>(sol_inf.time_obtained.count()) / 1e6 << ", for K = "<< std::to_string(inf_K) << std::endl;
        std::cout << "with path : ";
        for (auto it = sol_inf.sol.begin(); it != sol_inf.sol.end(); ++it) {
            std::cout << *it << " ";
        }
        std::cout << std::endl;
        std::cout << "borne sup = " << sol_sup.robust_score << ", weight = " << sol_sup.robust_constraint << ", for S = " << inst.S << ", found in "<< static_cast<double>(sol_sup.time_obtained.count()) / 1e6 << ", for K = "<< std::to_string(sup_K) <<std::endl;
        std::cout << "with path : ";
        for (auto it = sol_sup.sol.begin(); it != sol_sup.sol.end(); ++it) {
            std::cout << *it << " ";
        }
        std::cout << std::endl;

        if ((static_cast<double>(duration.count()) / 1e6) > max_duration) {
            std::cout << "ended because of duration after " << (static_cast<double>(duration.count()) / 1e6) << " s " << std::endl;
        }
        if (iter == max_iter) {
            std::cout << "ended because of iter number " << iter << "/" << max_iter << std::endl;
        }
        if (sup_K-inf_K <= precision_K) {
            std::cout << "ended because of K gap with inf_K = " << inf_K << ", sup_K = " << sup_K << ", and precision_K = "<< precision_K << std::endl;
        }
    }

    if (verbose>1) {
        std::vector<IloInt> sol = astar_solve(inst, inf_K, verbose);
        std::cout << "Solution  for K = " << inf_K << std::endl;
        for (auto it = sol.begin(); it != sol.end(); ++it) {
            std::cout << *it << " ";
        }
        cout << " " << endl;
        std::cout << "static objective = " << inst.compute_static_score(sol) << std::endl;
        std::cout << "robust objective = " << inst.compute_robust_score_knapsack(sol) << std::endl;
        std::cout << "static constraint = " << inst.compute_static_constraint(sol) << std::endl;
        std::cout << "robust constraint  = " << inst.compute_robust_constraint_knapsack(sol) << " for S = " << inst.S << std::endl;

        sol = astar_solve(inst, sup_K, verbose);
        std::cout << "Solution  for K = " << sup_K << std::endl;
        for (auto it = sol.begin(); it != sol.end(); ++it) {
            std::cout << *it << " ";
        }
        std::cout << std::endl;
        std::cout << "static objective = " << inst.compute_static_score(sol) << std::endl;
        std::cout << "robust objective = " << inst.compute_robust_score_knapsack(sol) << std::endl;
        std::cout << "static constraint = " << inst.compute_static_constraint(sol) << std::endl;
        std::cout << "robust constraint  = " << inst.compute_robust_constraint_knapsack(sol) << " for S = " << inst.S << std::endl;
    }

    if (!inst.sol.empty()) {
        std::cerr << "Warning: solution vector not empty for instance " << inst.name << std::endl;
        inst.sol.clear();
    }
    for (auto it = sol_sup.sol.begin(); it != sol_sup.sol.end(); ++it) {
        inst.sol.push_back(*it);
    }
}


float NodesInfo::sum_knapsack_capa() const {
    float sum = 0;
    for (auto it = knapsack_dij.begin(); it != knapsack_dij.end(); ++it) {
        sum += std::get<1>(*it);
    }
    return sum;
}


void NodesInfo::compute_dist_robust(const float d1) {
    float sum = 0;
    float weight = 0;
    for (auto it = knapsack_dij.begin(); it != knapsack_dij.end(); ++it) {
        if (weight + std::get<1>(*it) <= d1) {
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


void NodesInfo::compute_nodes_robust(const float d2) {
    float sum = 0;
    float weight = 0;
    for (auto it = knapsack_phi.begin(); it != knapsack_phi.end(); ++it) {
        if (weight + 2 <= d2) {
            sum += 2*(*it);
            weight += 2;
        } else {
            sum += (*it)*(d2-weight);
            break;
        }
    }
    robust_dist_nodes = sum;
}


void NodesInfo::redo_knapsack_dij(const Instance& inst, const int& new_node) {
    float new_dij = inst.d[parent][new_node];
    float new_Dij = inst.D[parent][new_node];

    knapsack_dij.push_back(std::make_tuple(new_dij, new_Dij));
    std::sort(knapsack_dij.begin(), knapsack_dij.end(), cmp_dij);
    float sum_knapsack = sum_knapsack_capa();
    while(sum_knapsack-std::get<1>(knapsack_dij.back()) >= inst.d1) {
        std::tuple<float,float> last = knapsack_dij.back();
        knapsack_dij.pop_back();
        sum_knapsack -= std::get<1>(last);
    }
    compute_dist_robust(inst.d1);
}


void NodesInfo::redo_knapsack_phi(const Instance& inst, const int& new_node) {
    float new_phi = inst.ph[new_node];
    knapsack_phi.push_back(new_phi);
    std::sort(knapsack_phi.begin(), knapsack_phi.end(), std::greater<float>());
    float sum_knapsack = 2*knapsack_phi.size();
    while (sum_knapsack-2 >= inst.d2) {
        knapsack_phi.pop_back();
        sum_knapsack -= 2;
    }
    compute_nodes_robust(inst.d2);
}


SolutionInfo::SolutionInfo(const Instance& inst, const std::vector<IloInt>& solu, const double& K_, const std::chrono::steady_clock::time_point& start) {
    sol = solu;
    K = K_;
    time_obtained = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - start);
    static_score = inst.compute_static_score(sol);
    robust_score = inst.compute_robust_score_knapsack(sol);
    static_constraint = inst.compute_static_constraint(sol);
    robust_constraint =  inst.compute_robust_constraint_knapsack(sol);
}
