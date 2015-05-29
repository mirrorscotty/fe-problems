/**
 * @file lin-genmaxwell.c
 *
 * Viscoelasticity differential equations for generalized Maxwell model
 * Primary equation
 * \f[
 * \sigma = E_0\epsilon + \sum_{m=1}^M E_m q^m
 * \f]
 * Equation for each Maxwell element included in the model
 * \f[
 * \frac{\partial}{\partial t}q^m + \frac{1}{\lambda_m}q^m = \frac{\partial}{\partial t}\epsilon
 * \f]
 *
 * Weak forms
 * \f[
 * R_i^T = \int_{x_i}^{x_{i+1}} \left\{D_0\epsilon_i\phi_i + \sum_{m=1}^M D_m q_i^m\phi_i - \sigma\right\}\phi_j dx = 0
 * \f]
 * \f[
 * R_i^{P_m} = \int_{x_i}^{x_{i+1}} \left\{\frac{\partial \epsilon_i}{\partial t}\phi_i - \frac{\partial q_i^m}{\partial t}\phi_i -\frac{1}{\lambda_m} q_i^m\phi_i\right\}\phi_j dx = 0
 * \f]
 */

#include "matrix.h"
#include "fe-solver.h"
#include "material-data.h"
#include "common.h"
#include "deformation.h"
#include <math.h>
#include <stdlib.h>

#define STRESS0(T) EffPorePress(CINIT/CAMB, (T))
#define STRESS(X, T) (EffPorePress((X), (T)) / STRESS0(T))

#define EA(M, T) \
    1e6 * (68.18*(1/(1+exp(((M)-250.92*exp(-0.0091*(T)))/2.19))+0.078)) / STRESS0(T)
#define E1(M, T) \
    1e6 * (20.26*exp(-0.0802*((M)+0.0474*(T)-14.238))) / STRESS0(T)
#define E2(M, T) \
    1e6 * (2.484 + 6.576/(1+exp(((M)-19.36)/0.848))) / STRESS0(T)
#define LAMBDA1 scaleTime(p->chardiff, 7)
#define LAMBDA2 scaleTime(p->chardiff, 110)

/**
 * Derivative of the main differential equation with respect to strain
 */
double ResSolid_dTde(struct fe1d *p, matrix *guess, Elem1D *elem,
                    double x, int f1, int f2)
{
    double T = TINIT,
           C = 0,
           Ci, Ea = 0, value;
    int i;
    basis *b;
    b = p->b;

    solution *s;
    s = CreateSolution(p->t, p->dt, guess);

    for(i=0; i<b->n; i++) {
        Ci = EvalSoln1D(p, CVAR, elem, s, valV(elem->points, i));
        Ci = uscaleTemp(p->chardiff, Ci);
        
        C += Ci * b->phi[i](x);
        Ea += EA(Ci, T) * b->phi[i](x);
    }

    value = Ea * b->phi[f1](x) * b->phi[f2](x) / IMap1D(p, elem, x);
    free(s);
    return value;
}

/**
 * Derivative of the main differential equation with respect to the partial
 * stress in the first Maxwell element.
 */
double ResSolid_dTdq1(struct fe1d *p, matrix *guess, Elem1D *elem,
                      double x, int f1, int f2)
{
    double T = TINIT,
           C = 0,
           e1 = 0,
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
        e1 += E1(Ci, T) * b->phi[i](x);
    }

    value = e1 * b->phi[f1](x) * b->phi[f2](x) / IMap1D(p, elem, x);
    free(s);
    return value;
}

/**
 * Derivative of the main differential equation with respect to the partial
 * stress in the second Maxwell element.
 */
double ResSolid_dTdq2(struct fe1d *p, matrix *guess, Elem1D *elem,
                      double x, int f1, int f2)
{
    double T = TINIT,
           C = 0,
           e2 = 0,
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
        e2 += E2(Ci, T) * b->phi[i](x);
    }

    value = e2 * b->phi[f1](x) * b->phi[f2](x) / IMap1D(p, elem, x);
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
double ResSolid_dP1dq1(struct fe1d *p, matrix *guess, Elem1D *elem,
                      double x, int f1, int f2)
{
    return 1/LAMBDA1
        * p->b->phi[f1](x) * p->b->phi[f2](x) / IMap1D(p, elem, x);
}

/**
 * Derivative of the stress equation for the second maxwell element  with
 * respect to the partial stress in that element.
 */
double ResSolid_dP2dq2(struct fe1d *p, matrix *guess, Elem1D *elem,
                      double x, int f1, int f2)
{
    return 1/LAMBDA2
        * p->b->phi[f1](x) * p->b->phi[f2](x) / IMap1D(p, elem, x);
}

double ResDtSolid_dP1de(struct fe1d *p, matrix *guess, Elem1D *elem,
                      double x, int f1, int f2)
{
    return -1 * p->b->phi[f1](x) * p->b->phi[f2](x) / IMap1D(p, elem, x);
}

double ResDtSolid_dP1dq1(struct fe1d *p, matrix *guess, Elem1D *elem,
                      double x, int f1, int f2)
{
    return p->b->phi[f1](x) * p->b->phi[f2](x) / IMap1D(p, elem, x);
}

double ResDtSolid_dP2de(struct fe1d *p, matrix *guess, Elem1D *elem,
                      double x, int f1, int f2)
{
    return -1 * p->b->phi[f1](x) * p->b->phi[f2](x) / IMap1D(p, elem, x);
}

double ResDtSolid_dP2dq2(struct fe1d *p, matrix *guess, Elem1D *elem,
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
           Ci;
    int i;
    solution *s;
    basis *b;

    b = p->b;
    s = CreateSolution(p->t, p->dt, guess);

    for(i=0; i<b->n; i++) {
        Ci = EvalSoln1D(p, CVAR, elem, s, valV(elem->points, i));
        Ci = uscaleTemp(p->chardiff, Ci);
        C += Ci;
        sigma = STRESS(Ci, T) * b->phi[i](x);
    }
    free(s);
    return -1*sigma;
}

