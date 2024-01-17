// dualized_formulation.h
#ifndef DUALIZED_FORMULATION_H
#define DUALIZED_FORMULATION_H

#include "parser.h"

void dualized_solve(IloEnv env, Instance& inst, const unsigned int& time_limit=10, const int& verbose=0);

void dualized_solve_2(IloEnv env, Instance& inst, const unsigned int& time_limit=10, const int& verbose=0);


#endif