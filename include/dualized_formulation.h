// dualized_formulation.h
#ifndef DUALIZED_FORMULATION_H
#define DUALIZED_FORMULATION_H

#include "solve_method.h"
class Instance;
class IloEnv;


struct DualizedMethod : public SolveMethod {
    DualizedMethod() {method_name = "dualized";}
    void solve(IloEnv& env, Instance& inst, const unsigned int& time_limit, const bool& reduce_symetry, const int& verbose) override;
};

#endif