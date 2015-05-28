#ifndef DEFORMATION_H
#define DEFORMATION_H

#include "fe-solver.h"

double DeformationGrad(struct fe1d*, double, double);
double DeformGradPc(struct fe1d*, double, double);
double DeformGradBeta(struct fe1d*, double, double);
double FindPoisson(struct fe1d*, double, double);
double Porosity(struct fe1d*, double, int);

double EffPorePress(double, double);

#endif

