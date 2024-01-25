// plans_coupants.h
#ifndef PLANS_COUPANTS_H
#define PLANS_COUPANTS_H

class Instance; // forward declaration
class IloEnv;

void plans_coupants_solve(IloEnv env, Instance& inst, const unsigned int& time_limit=300, const int& verbose=0);

#endif