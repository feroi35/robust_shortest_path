// static_solve.h
#ifndef STATIC_SOLVE_H
#define STATIC_SOLVE_H

#include "solve_method.h"
class Instance;
class IloEnv;

struct StaticMethod : public SolveMethod {
    StaticMethod() {method_name = "static";}
    void solve(IloEnv& env, Instance& inst, const unsigned int& time_limit, const int& verbose) override;
};

#endif