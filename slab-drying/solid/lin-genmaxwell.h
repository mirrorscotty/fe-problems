#ifndef LIN_GENMAXWELL_H
#define LIN_GENMAXWELL_H

double ResSolid_dTde(struct fe1d *, matrix *, Elem1D *, double, int, int);
double ResSolid_dTdq1(struct fe1d *, matrix *, Elem1D *, double, int, int);
double ResSolid_dTdq2(struct fe1d *, matrix *, Elem1D *, double, int, int);
double ResSolid_zero(struct fe1d *, matrix *, Elem1D *, double, int, int);
double ResSolid_dP1dq1(struct fe1d *, matrix *, Elem1D *, double, int, int);
double ResSolid_dP2dq2(struct fe1d *, matrix *, Elem1D *, double, int, int);

double ResDtSolid_dP1de(struct fe1d *, matrix *, Elem1D *, double, int, int);
double ResDtSolid_dP1dq1(struct fe1d *, matrix *, Elem1D *, double, int, int);
double ResDtSolid_dP2de(struct fe1d *, matrix *, Elem1D *, double, int, int);
double ResDtSolid_dP2dq2(struct fe1d *, matrix *, Elem1D *, double, int, int);

double ResFSolid_T(struct fe1d *, matrix *, Elem1D *, double, int, int);

#endif

