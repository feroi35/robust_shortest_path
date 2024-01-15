#include "main.h"
#include "parser.h"
#include "static_solve.h"
#include "dualized_formulation.h"
#include <cstdlib>

// #include <string>


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
    unsigned int time_limit = 60;
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
        std::cout << "Method not recognized" << method << std::endl;
        std::cout << "Method should be either 'static' or 'dualized'" << std::endl;
    }

    if(verbose > 0){
        double obj; 
        obj = instance.compute_static_score(verbose);
        std::cout << "obj static = " << obj << std::endl;
        double robust_associated_obj;
        robust_associated_obj = instance.compute_robust_score(env, instance.sol, time_limit_robust_obj, verbose);
        std::cout << "obj robust = " << robust_associated_obj << std::endl;
        double robust_associated_cstr;
        robust_associated_cstr = instance.compute_robust_constraint(env, instance.sol, time_limit_robust_cstr, verbose);
        std::cout << "cstr robust = " << robust_associated_cstr << " with S = " << instance.S << std::endl;
    }

    env.end(); 
    return 0;
}
