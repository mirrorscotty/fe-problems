#ifndef KINEMATICS_H
#define KINEMATICS_H

#include "fe-solver.h"

double ResDtSolid_dvde(struct fe1d *, matrix *, Elem1D *, double, int, int);
double ResSolid_dvdv(struct fe1d *, matrix *, Elem1D *, double, int, int);

double ResSolid_dudu(struct fe1d *, matrix *, Elem1D *, double, int, int);
double ResFSolid_u(struct fe1d *, matrix *, Elem1D *, double, int, int);
double ResSolid_dude(struct fe1d *, matrix *, Elem1D *, double, int, int);

#endif

