// branch_and_cut.h
#ifndef BRANCH_AND_CUT_H
#define BRANCH_AND_CUT_H

#include "solve_method.h"
class Instance;
class IloEnv;

struct BranchAndCutMethod : public SolveMethod {
    BranchAndCutMethod() {method_name = "branch_and_cut";}
    void solve(IloEnv& env, Instance& inst, const unsigned int& time_limit, const int& verbose) override;
};

#endif