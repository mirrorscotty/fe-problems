/**
 * @file lin-genkelvin.c
 *
 * The viscoelasticity equations are taken from The Finite Element Method for
 * Solid and Structural Mechanics (Seventh Edition) by Zienkiewicz, Taylor, and
 * Fox.
 * Viscoelasticity differential equations for generalized Kelvin model
 * Primary equation
 * \f[
 * \epsilon = J_0\sigma + \sum_{m=1}^M J_m r^m
 * \f]
 * Equation for each Kelvin element included in the model
 * \f[
 * \frac{\partial}{\partial t}r^m + \frac{1}{\tau_m}r^m = \sigma
 * \f]
 *
 * Weak forms
 * \f[
 * R_i^T = \int_{x_i}^{x_{i+1}} \left\{J_0\sigma + \sum_{m=1}^M J_m r_i^m\phi_i - \epsilon_i\phi_i\right\}\phi_j dx = 0
 * \f]
 * \f[
 * R_i^{P_m} = \int_{x_i}^{x_{i+1}} \left\{\frac{\partial r_i^m}{\partial t}\phi_i + \frac{1}{\tau_m} r_i^m\phi_i - \sigma\right\}\phi_j dx = 0
 * \f]
 */

#include "matrix.h"
#include "fe-solver.h"
#include "material-data.h"
#include "../common-kelvin.h"
#include "deformation.h"
#include <math.h>
#include <stdlib.h>

#define STRESS0(T) EffPorePress(CINIT/Camb, (T))
#define STRESS(X, T, e) \
    ( EffPorePress((X), (T))/STRESS0((T)) * .0612 )

/* Stress relaxation parameters from Rozzi */
/*
#define EA(M, T) \
    1e6 * (68.18*(1/(1+exp(((M)*100-250.92*exp(-0.0091*(T)))/2.19))+0.078)) / STRESS0(T)
#define E1(M, T) \
    1e6 * (20.26*exp(-0.0802*((M)*100+0.0474*(T)-14.238))) / STRESS0(T)
#define E2(M, T) \
    1e6 * (2.484 + 6.576/(1+exp(((M)*100-19.36)/0.848))) / STRESS0(T)
#define LAMBDA1 scaleTime(p->chardiff, 7)
#define LAMBDA2 scaleTime(p->chardiff, 110)
*/

/* T=333, M=.05 */
/*
#define J0(M, T) 1.081284e-08 * STRESS0(T)
#define J1(M, T) 1.616248e-09 * STRESS0(T)
#define J2(M, T) 1.614146e-09 * STRESS0(T)
#define TAU1(M, T) scaleTime(p->chardiff, 8.057304e+00)
#define TAU2(M, T) scaleTime(p->chardiff, 1.241419e+02)
*/

/* T=333, M=.4 */
/*
#define J0(M, T) 1.036386e-07 * STRESS0(T)
#define J1(M, T) 2.268417e-08 * STRESS0(T)
#define J2(M, T) 6.040908e-08 * STRESS0(T)
#define TAU1(M, T) scaleTime(p->chardiff, 8.575243e+00)
#define TAU2(M, T) scaleTime(p->chardiff, 1.618177e+02)
*/
#define CREEPFILE "output.csv"
#define J0(M, T) CreepLookupJ0(CREEPFILE, T, M) * STRESS0(T)
#define J1(M, T) CreepLookupJ1(CREEPFILE, T, M) * STRESS0(T)
#define J2(M, T) CreepLookupJ2(CREEPFILE, T, M) * STRESS0(T)
#define TAU1(M, T) scaleTime(p->chardiff, CreepLookupTau1(CREEPFILE, T, M))
#define TAU2(M, T) scaleTime(p->chardiff, CreepLookupTau2(CREEPFILE, T, M))

extern double Camb;

/**
 * Derivative of the main differential equation with respect to strain
 */
double ResSolid_dTde(struct fe1d *p, matrix *guess, Elem1D *elem,
                    double x, int f1, int f2)
{
    double value;
    basis *b;

    b = p->b;

    value = -1 * b->phi[f1](x) * b->phi[f2](x) / IMap1D(p, elem, x);
    return value;
}

/**
 * Derivative of the main differential equation with respect to the partial
 * stress in the first Maxwell element.
 */
double ResSolid_dTdr1(struct fe1d *p, matrix *guess, Elem1D *elem,
                      double x, int f1, int f2)
{
    double T = TINIT,
           C = 0,
           j1 = 0,
           Ci, value;
    int i;
    basis *b;
    b = p->b;

    solution *s;
    s = CreateSolution(p->t, p->dt, guess);

    for(i=0; i<b->n; i++) {
        Ci = EvalSoln1D(p, CVAR, elem, s, valV(elem->points, i));
        Ci = uscaleTemp(p->chardiff, Ci);

        C += Ci * b->phi[i](x);
        j1 += J1(Ci, T) * b->phi[i](x);
    }

    value = j1 * b->phi[f1](x) * b->phi[f2](x) / IMap1D(p, elem, x);
    free(s);
    return value;
}

/**
 * Derivative of the main differential equation with respect to the partial
 * stress in the second Maxwell element.
 */
double ResSolid_dTdr2(struct fe1d *p, matrix *guess, Elem1D *elem,
                      double x, int f1, int f2)
{
    double T = TINIT,
           C = 0,
           j2 = 0,
           Ci, value;
    int i;
    basis *b;
    b = p->b;

    solution *s;
    s = CreateSolution(p->t, p->dt, guess);

    for(i=0; i<b->n; i++) {
        Ci = EvalSoln1D(p, CVAR, elem, s, valV(elem->points, i));
        Ci = uscaleTemp(p->chardiff, Ci);

        C += Ci * b->phi[i](x);
        j2 += J2(Ci, T) * b->phi[i](x);
    }

    value = j2 * b->phi[f1](x) * b->phi[f2](x) / IMap1D(p, elem, x);
    free(s);
    return value;
}

double ResSolid_zero(struct fe1d *p, matrix *guess, Elem1D *elem,
                     double x, int f1, int f2)
{
    return 0;
}

/**
 * Derivative of the stress equation for the first maxwell element  with
 * respect to the partial stress in that element.
 */
double ResSolid_dP1dr1(struct fe1d *p, matrix *guess, Elem1D *elem,
                      double x, int f1, int f2)
{
    double T = TINIT,
           C = 0,
           tau1 = 0,
           Ci;
    int i;
    basis *b;
    b = p->b;

    solution *s;
    s = CreateSolution(p->t, p->dt, guess);

    for(i=0; i<b->n; i++) {
        Ci = EvalSoln1D(p, CVAR, elem, s, valV(elem->points, i));
        Ci = uscaleTemp(p->chardiff, Ci);

        C += Ci * b->phi[i](x);
        tau1 += TAU1(Ci, T) * b->phi[i](x);
    }
    free(s);
    return 1/tau1
        * p->b->phi[f1](x) * p->b->phi[f2](x) / IMap1D(p, elem, x);
}

/**
 * Derivative of the stress equation for the second maxwell element  with
 * respect to the partial stress in that element.
 */
double ResSolid_dP2dr2(struct fe1d *p, matrix *guess, Elem1D *elem,
                      double x, int f1, int f2)
{
    double T = TINIT,
           C = 0,
           tau2 = 0,
           Ci;
    int i;
    basis *b;
    b = p->b;

    solution *s;
    s = CreateSolution(p->t, p->dt, guess);

    for(i=0; i<b->n; i++) {
        Ci = EvalSoln1D(p, CVAR, elem, s, valV(elem->points, i));
        Ci = uscaleTemp(p->chardiff, Ci);

        C += Ci * b->phi[i](x);
        tau2 += TAU2(Ci, T) * b->phi[i](x);
    }
    free(s);
    return 1/tau2
        * p->b->phi[f1](x) * p->b->phi[f2](x) / IMap1D(p, elem, x);
}

double ResDtSolid_dP1dr1(struct fe1d *p, matrix *guess, Elem1D *elem,
                      double x, int f1, int f2)
{
    return p->b->phi[f1](x) * p->b->phi[f2](x) / IMap1D(p, elem, x);
}

double ResDtSolid_dP2dr2(struct fe1d *p, matrix *guess, Elem1D *elem,
                      double x, int f1, int f2)
{
    return p->b->phi[f1](x) * p->b->phi[f2](x) / IMap1D(p, elem, x);
}


double ResFSolid_T(struct fe1d *p, matrix *guess, Elem1D *elem,
                   double x, int f1, int a)
{
    double T = TINIT,
           C = 0,
           sigma = 0,
           epsilon = 0,
           j0 = 0,
           Ci, ei;
    int i;
    solution *s;
    basis *b;

    b = p->b;
    s = CreateSolution(p->t, p->dt, guess);

    for(i=0; i<b->n; i++) {
        Ci = EvalSoln1D(p, CVAR, elem, s, valV(elem->points, i));
        Ci = uscaleTemp(p->chardiff, Ci);
        ei = EvalSoln1D(p, STVAR, elem, s, valV(elem->points, i));
        epsilon += ei * b->phi[i](x);
        C += Ci * b->phi[i](x);
        j0 += J0(Ci, T) * b->phi[i](x);
        sigma += STRESS(Ci, T, epsilon) * b->phi[i](x);
    }
    free(s);
    return sigma*j0 / IMap1D(p, elem, x);
}

double ResFSolid_P1(struct fe1d *p, matrix *guess, Elem1D *elem,
                   double x, int f1, int a)
{
    double T = TINIT,
           C = 0,
           sigma = 0,
           epsilon = 0,
           Ci, ei;
    int i;
    solution *s;
    basis *b;

    b = p->b;
    s = CreateSolution(p->t, p->dt, guess);

    for(i=0; i<b->n; i++) {
        Ci = EvalSoln1D(p, CVAR, elem, s, valV(elem->points, i));
        Ci = uscaleTemp(p->chardiff, Ci);
        ei = EvalSoln1D(p, STVAR, elem, s, valV(elem->points, i));
        epsilon += ei * b->phi[i](x);
        C += Ci * b->phi[i](x);
        sigma += STRESS(Ci, T, epsilon) * b->phi[i](x);
    }
    free(s);
    return -1*sigma / IMap1D(p, elem, x);
}

double ResFSolid_P2(struct fe1d *p, matrix *guess, Elem1D *elem,
                   double x, int f1, int a)
{
    double T = TINIT,
           C = 0,
           sigma = 0,
           epsilon = 0,
           Ci, ei;
    int i;
    solution *s;
    basis *b;

    b = p->b;
    s = CreateSolution(p->t, p->dt, guess);

    for(i=0; i<b->n; i++) {
        Ci = EvalSoln1D(p, CVAR, elem, s, valV(elem->points, i));
        Ci = uscaleTemp(p->chardiff, Ci);
        ei = EvalSoln1D(p, STVAR, elem, s, valV(elem->points, i));

        epsilon += ei * b->phi[i](x);
        C += Ci * b->phi[i](x);

        sigma += STRESS(Ci, T, epsilon) * b->phi[i](x);
    }
    free(s);
    return -1*sigma / IMap1D(p, elem, x);
}

