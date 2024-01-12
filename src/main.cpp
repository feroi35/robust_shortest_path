#include "main.h"
#include "parser.h"
#include "static_solve.h"


int main(int argc, char **argv) {

    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <filename>" << std::endl;
        exit(1);
    }
    char* filename = argv[1];
    unsigned int time_limit = 60;
    char* method = argv[2];

    unsigned int time_limit_robust_obj = 60;
    unsigned int time_limit_robust_cstr = 60;


    if(strcmp(method, "static") == 0) {
        // mettre une disjonction sur la valeur de methode
        IloEnv env;
        Instance instance(env, filename);
        // instance.display();
        static_solve(env, instance, time_limit);
        double obj; 
        obj = instance.compute_static_score();
        std::cout << "obj static = " << obj << std::endl;

        double robust_associated_obj;
        robust_associated_obj = instance.compute_robust_score(env, instance.sol, time_limit_robust_obj);
        std::cout << "obj robust = " << robust_associated_obj << std::endl;
        double robust_associated_cstr;
        robust_associated_cstr = instance.compute_robust_constraint(env, instance.sol, time_limit_robust_cstr);
        std::cout << "cstr robust = " << robust_associated_cstr << " with S = " << instance.S << std::endl;

        env.end(); 

    }
    // // mettre une disjonction sur la valeur de methode
    // IloEnv env;
    // Instance instance(env, filename);
    // // instance.display();
    // static_solve(env, instance, time_limit);
    // double obj; 
    // obj = instance.compute_static_score();
    // std::cout << "obj static = " << obj << std::endl;

    // double robust_associated_obj;
    // robust_associated_obj = instance.compute_robust_score(env, instance.sol, time_limit_robust_obj);
    // std::cout << "obj robust = " << robust_associated_obj << std::endl;
    // double robust_associated_cstr;
    // robust_associated_cstr = instance.compute_robust_constraint(env, instance.sol, time_limit_robust_cstr);
    // std::cout << "cstr robust = " << robust_associated_cstr << " with S = " << instance.S << std::endl;

    // env.end();    
    return 0;
}
