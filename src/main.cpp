#include "main.h"
#include "parser.h"
#include "static_solve.h"
#include "dualized_formulation.h"
#include "heuristics.h"
#include "branch_and_cut.h"
#include "plans_coupants.h"


int main(int argc, char **argv) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <filename>" << "method_name" << "(time_limit)" << "(verbose)" << std::endl;
        exit(1);
    }
    char* filename = argv[1];
    char* method_name = argv[2];
    unsigned int time_limit = 60;
    if (argc > 3) {
        time_limit = atoi(argv[3]);
    }
    unsigned int verbose = 0;
    if (argc > 4) {
        verbose = atoi(argv[4]);
    }

    IloEnv env;
    Instance instance(env, filename);
    if (verbose > 3) instance.display();

    std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();

    try {
        if (strcmp(method_name, "heuristic") == 0) {
            double precision_K = 1e-5;
            int max_iter = 2000;
            float max_duration = 20.0;
            HeuristicMethod method(precision_K, max_iter, max_duration);
            method.solve_and_display(env, instance, time_limit, verbose);
        }
        else if(strcmp(method_name, "static") == 0) {
            StaticMethod method;
            method.solve_and_display(env, instance, time_limit, verbose);
        } else if (strcmp(method_name, "dualized") == 0) {
            DualizedMethod method(false);
            method.solve_and_display(env, instance, time_limit, verbose);
        } else if (strcmp(method_name, "branch_and_cut") == 0) {
            BranchAndCutMethod method;
            method.solve_and_display(env, instance, time_limit, verbose);
        } else if (strcmp(method_name, "plans_coupants") == 0) {
            PlansCoupantsMethod method;
            method.solve_and_display(env, instance, time_limit, verbose);
        } else {
            std::cerr << "Method_name not recognized: " << method_name << std::endl;
            exit(2);
        }

        if (verbose > 0) {
            std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
            std::chrono::microseconds duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);

            std::cout << std::endl;
            std::cout << "static objective = " << instance.compute_static_score() << std::endl;
            std::cout << "robust objective MILP = " << instance.compute_robust_score_milp(env, time_limit, verbose) << std::endl;
            std::cout << "robust objective knapsack = " << instance.compute_robust_score_knapsack(verbose) << std::endl;
            std::cout << "static constraint = " << instance.compute_static_constraint() << std::endl;
            std::cout << "robust constraint MILP = " << instance.compute_robust_constraint_milp(env, time_limit, verbose) << std::endl;
            std::cout << "robust constraint knapsack = " << instance.compute_robust_constraint_knapsack(verbose) << std::endl;
            std::cout << "S = " << instance.S << std::endl;
            std::cout << "Time taken = " << static_cast<double>(duration.count()) / 1e6 << std::endl;
        }
    } catch (IloException& e) {
        std::cerr << "Ilo exception caught: " << e << std::endl;
    } catch (std::domain_error& e) {
        std::cerr << "Domain error caught: " << e.what() << std::endl;
    } catch (std::exception& e) {
        std::cerr << "Other exception caught: " << e.what() << std::endl;
    }

    env.end();
    return 0;
}
