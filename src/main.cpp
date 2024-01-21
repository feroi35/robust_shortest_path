#include <cstdlib>
#include "main.h"
#include "parser.h"
#include "static_solve.h"
#include "dualized_formulation.h"


int main(int argc, char **argv) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <filename>" << "method" << "(verbose)" << std::endl;
        exit(1);
    }
    char* filename = argv[1];
    char* method = argv[2];
    int verbose = 0;
    if (argc > 3) {
        verbose = atoi(argv[3]);
    }
    unsigned int time_limit = 300;

    IloEnv env;
    Instance instance(env, filename);
    if (verbose > 3)
        instance.display();

    try {
        if(strcmp(method, "static") == 0) {
            static_solve(env, instance, time_limit, verbose);
        } else if (strcmp(method, "dualized") == 0){
            dualized_solve(env, instance, time_limit, verbose);
        } else {
            std::cerr << "Method not recognized: " << method << std::endl;
            std::cerr << "Method should be either 'static' or 'dualized'" << std::endl;
            exit(1);
        }
        if (verbose > 0) {
            std::cout << std::endl;
            std::cout << "static objective = " << instance.compute_static_score() << std::endl;
            std::cout << "robust objective = " << instance.compute_robust_score(env) << std::endl;
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
