// static_solve.h
#ifndef STATIC_SOLVE_H
#define STATIC_SOLVE_H

#include "parser.h"

void static_solve(IloEnv env, Instance& inst, const unsigned int& time_limit=60, const int& verbose=0);

#endif