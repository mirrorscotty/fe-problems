#ifndef LIN_MAXWELL_H
#define LIN_MAXWELL_H

#include "matrix.h"
#include "fe-solver.h"

struct fe1d;

double ResSolid(struct fe1d *, matrix *, Elem1D *, double, int, int);
double ResDtSolid(struct fe1d *, matrix *, Elem1D *, double, int, int);
double ResFSolid(struct fe1d *, matrix *, Elem1D *, double, int, int);

#endif

