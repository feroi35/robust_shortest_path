// plans_coupants.h
#ifndef PLANS_COUPANTS_H
#define PLANS_COUPANTS_H

#include "solve_method.h"
class Instance;
class IloEnv;


struct PlansCoupantsMethod : public SolveMethod {
    PlansCoupantsMethod() {method_name = "plans_coupants";}
    void solve(IloEnv& env, Instance& inst, const unsigned int& time_limit, const int& verbose) override;
};

#endif