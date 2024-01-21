// dualized_formulation.h
#ifndef DUALIZED_FORMULATION_H
#define DUALIZED_FORMULATION_H

class Instance; // forward declaration
class IloEnv;

void dualized_solve(IloEnv env, Instance& inst, const unsigned int& time_limit=10, const int& verbose=0);

#endif