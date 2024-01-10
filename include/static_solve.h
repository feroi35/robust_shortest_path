// static_solve.h
#ifndef STATIC_SOLVE_H
#define STATIC_SOLVE_H

#include "parser.h"

void static_solve(IloEnv env, Instance& inst, unsigned int time_limit=10);

#endif