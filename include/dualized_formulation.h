// dualized_formulation.h
#ifndef DUALIZED_FORMULATION_H
#define DUALIZED_FORMULATION_H

#include "solve_method.h"
class Instance;
class IloEnv;


struct DualizedMethod : public SolveMethod {
    bool reduce_symetry;
    
    DualizedMethod(const bool& reduce_sym) {method_name = "dualized"; reduce_symetry = reduce_sym;};
    void solve(IloEnv& env, Instance& inst, const unsigned int& time_limit, const int& verbose) override;
};

#endif