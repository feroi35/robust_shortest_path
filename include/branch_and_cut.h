// branch_and_cut.h
#ifndef BRANCH_AND_CUT_H
#define BRANCH_AND_CUT_H

class Instance; // forward declaration
class IloEnv;

void branch_and_cut_solve(IloEnv env, Instance& inst, const unsigned int& time_limit=10, const int& verbose=0);

#endif