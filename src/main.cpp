#include "main.h"
#include "parser.h"
#include "static_solve.h"


int main(int argc, char **argv) {

    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <filename>" << std::endl;
        exit(1);
    }
    char* filename = argv[1];
    unsigned int time_limit = 10;

    IloEnv env;
    Instance instance(env, filename);
    // instance.display();
    static_solve(env, instance, time_limit);

    env.end();
    return 0;
}
