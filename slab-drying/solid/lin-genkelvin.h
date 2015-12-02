#ifndef LIN_GENMAXWELL_H
#define LIN_GENMAXWELL_H

double ResSolid_dTde(struct fe1d *, matrix *, Elem1D *, double, int, int);
double ResSolid_dTdr1(struct fe1d *, matrix *, Elem1D *, double, int, int);
double ResSolid_dTdr2(struct fe1d *, matrix *, Elem1D *, double, int, int);
double ResSolid_zero(struct fe1d *, matrix *, Elem1D *, double, int, int);
double ResSolid_dP1dr1(struct fe1d *, matrix *, Elem1D *, double, int, int);
double ResSolid_dP2dr2(struct fe1d *, matrix *, Elem1D *, double, int, int);

double ResDtSolid_dP1de(struct fe1d *, matrix *, Elem1D *, double, int, int);
double ResDtSolid_dP1dr1(struct fe1d *, matrix *, Elem1D *, double, int, int);
double ResDtSolid_dP2de(struct fe1d *, matrix *, Elem1D *, double, int, int);
double ResDtSolid_dP2dr2(struct fe1d *, matrix *, Elem1D *, double, int, int);

double ResFSolid_T(struct fe1d *, matrix *, Elem1D *, double, int, int);
double ResFSolid_P1(struct fe1d *, matrix *, Elem1D *, double, int, int);
double ResFSolid_P2(struct fe1d *, matrix *, Elem1D *, double, int, int);

double ResSolid_dSdS(struct fe1d *, matrix *, Elem1D *, double, int, int);
double ResSolid_dTdS(struct fe1d *, matrix *, Elem1D *, double, int, int);
double ResSolid_dP1dS(struct fe1d *, matrix *, Elem1D *, double, int, int);
double ResSolid_dP2dS(struct fe1d *, matrix *, Elem1D *, double, int, int);
double ResFSolid_S(struct fe1d *, matrix *, Elem1D *, double, int, int);
#endif

