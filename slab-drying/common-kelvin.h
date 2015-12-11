#ifndef COMMON_H
#define COMMON_H

#define NEQ 3

#define CVAR 0 /* Concentration */
#define STVAR 1 /* Total strain */
#define SUVAR 2 /* Displacement */

//#define TAMB 310 // K
//#define TAMB 500 // K
//#define TINIT 273 //K
//#define TINIT 313.15 //K
//#define TINIT 353.15 //K
#ifndef TINIT
#define TINIT 333 //K
#endif
#ifndef HCONV
#define HCONV 50
#endif

#ifndef CAMB
//#define CAMB 0.0305371 // kg/kg db
//#define CAMB 0.160262 // kg/kg db
//#define CAMB 0.140196
//#define CAMB 0.074695
//#define CAMB .06
//#define CAMB 0.0861867 // kg/kg db (RH=.65, T=80C)
//#define CAMB 0.138 // kg/kg db (RH=.65, T=40C)
//#define CAMB 0.0975028 // kg/kg db (RH=.7, T=80C)
//#define CAMB 0.15049 // kg/kg db (RH=.7, T=40C)
#endif
#ifndef CINIT
#define CINIT .33 // kg/kg db
#endif
#ifndef KC_CONV
#define KC_CONV 5
#endif

#ifndef THICKNESS
#define THICKNESS 1e-3
#endif

#ifndef POISSON
#define POISSON .37
#endif

#define SHRINKAGE_DIFFEQ

//#define DIFF(X, T) DiffCh10((X), (T))
#define DIFF(X, T) DiffLitchfield((X), (T))
//#define DIFF(X, T) DiffAvg((X), (T))

#define CREEP(TIME, T, X, P, DERIV)  CreepLaura2((TIME), (T), (X), (P), (DERIV))

#define BETA .6

double DiffAvg(double, double);

matrix* CreateElementMatrix(struct fe1d *, Elem1D *, matrix *);
matrix* CreateDTimeMatrix(struct fe1d *, Elem1D *, matrix *);
matrix* CreateElementLoad(struct fe1d *, Elem1D *, matrix *);
int IsOnRightBoundary(struct fe1d *, int);
int IsOnLeftBoundary(struct fe1d *, int);
void ApplyAllBCs(struct fe1d *);
//double DeformationGrad(struct fe1d *, double, double);

#endif

