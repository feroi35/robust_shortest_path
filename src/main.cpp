#include "main.h"
#include "parser.h"
#include "static_solve.h"
#include "dualized_formulation.h"

// #include <string>


int main(int argc, char **argv) {

    if (argc < 4) {
        std::cerr << "Usage: " << argv[0] << " <filename>" << std::endl;
        exit(1);
    }
    char* filename = argv[1];
    unsigned int time_limit = 60;
    char* method = argv[2];

    // char* verb = argv[3];
    // int verbose = std::stoi(verb);

    int verbose = *argv[3] - '0';


    unsigned int time_limit_robust_obj = 60;
    unsigned int time_limit_robust_cstr = 60;


    if(strcmp(method, "static") == 0) {
        IloEnv env;
        Instance instance(env, filename);
        // instance.display();
        static_solve(env, instance, time_limit, verbose);

        if(verbose >= 1){
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
    }

    if(strcmp(method, "dualized") == 0){
        IloEnv env;
        Instance instance(env, filename);
        // instance.display();
        dualized_solve(env, instance, time_limit, verbose);

        if(verbose >= 1){
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
    }

    return 0;
}
