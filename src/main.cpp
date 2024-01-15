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
    unsigned int time_limit = 120;
    unsigned int time_limit_robust_obj = 60;
    unsigned int time_limit_robust_cstr = 60;

    IloEnv env;
    Instance instance(env, filename);
    // instance.display();

    if(strcmp(method, "static") == 0) {
        static_solve(env, instance, time_limit, verbose);
    } else if (strcmp(method, "dualized") == 0){
        dualized_solve(env, instance, time_limit, verbose);
    } else {
        std::cerr << "Method not recognized: " << method << std::endl;
        std::cerr << "Method should be either 'static' or 'dualized'" << std::endl;
        exit(1);
    }

    if(verbose > 0){
        std::cout << std::endl;
        double static_obj = instance.compute_static_score(verbose);
        std::cout << "static objective = " << static_obj << std::endl;
        double robust_obj = instance.compute_robust_score(env, instance.sol, time_limit_robust_obj, verbose-2);
        std::cout << "robust objective = " << robust_obj << std::endl;
        double robust_cstr = instance.compute_robust_constraint(env, instance.sol, time_limit_robust_cstr, verbose-2);
        std::cout << "robust constraint  = " << robust_cstr << " with S = " << instance.S << std::endl;
    }

    env.end(); 
    return 0;
}
