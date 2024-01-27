#include "main.h"
#include "parser.h"
#include "static_solve.h"
#include "dualized_formulation.h"
#include "heuristics.h"
#include "branch_and_cut.h"
#include "plans_coupants.h"

#include <typeinfo> // for debug purpose


int main(int argc, char **argv) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <filename>" << "method" << "(time_limit)" << "(verbose)" << std::endl;
        exit(1);
    }
    char* filename = argv[1];
    char* method = argv[2];
    unsigned int time_limit = 60;
    if (argc > 3) {
        time_limit = atoi(argv[3]);
    }
    int verbose = 0;
    if (argc > 4) {
        verbose = atoi(argv[4]);
    }

    IloEnv env;
    Instance instance(env, filename);
    if (verbose > 3) instance.display();

    try {
        if(strcmp(method, "heuristics") == 0){

            Heuristic_instance heur(instance);
            heur.inf_dist = heur.backward_dijkstra_distance(instance);
            heur.inf_dist_nodes = heur.backward_dijkstra_nodes(instance);

            if(verbose>0){
                std::cout <<"initialisation heur faite" << std::endl;
            }
            
            double precision_K = 0.00001;
            heur.complete_astar_solve(instance, env, precision_K, 2000, 20, verbose);
            
            return 0;
        }
        else if(strcmp(method, "static") == 0) {
            static_solve(env, instance, time_limit, verbose);
        } else if (strcmp(method, "dualized") == 0) {
            dualized_solve(env, instance, time_limit, verbose);
        } else if (strcmp(method, "branch_and_cut") == 0) {
            branch_and_cut_solve(env, instance, time_limit, verbose);
        } else if (strcmp(method, "plans_coupants") == 0) {
            plans_coupants_solve(env, instance, time_limit, verbose);
        } else {
            std::cerr << "Method not recognized: " << method << std::endl;
            exit(1);
        }
        if (verbose > 0) {
            std::cout << std::endl;
            std::cout << "static objective = " << instance.compute_static_score() << std::endl;
            std::cout << "robust objective = " << instance.compute_robust_score(env) << std::endl;
            std::cout << "robust objective knapsack = " << instance.compute_robust_score_bis() << std::endl;
            std::cout << "static constraint = " << instance.compute_static_constraint() << std::endl;
            std::cout << "robust constraint  = " << instance.compute_robust_constraint(env) << std::endl;
            std::cout << "S = " << instance.S << std::endl;
        }
    } catch (IloException& e) {
        std::cerr << "Ilo exception caught: " << e << std::endl;
    } catch (std::domain_error& e) {
        std::cerr << "Domain error caught: " << e.what() << std::endl;
    } catch (...) {
        std::cerr << "Unknown exception caught" << std::endl;
    }

    env.end();
    return 0;
}
