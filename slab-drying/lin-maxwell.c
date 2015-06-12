/**
 * @file lin-maxwell.c
 * Infintesimal strain constituative model for the three parameter standard
 * linear solid model.
 * \f[
 * \sigma + \lambda\dot{\sigma} = E_\infty\epsilon + \lambda(E_\infty+E)\dot{\epsilon}
 * \f]
 *
 *  Finite element formulation:
 *  \f[
 *  R_S^i = \int_0^L\left\{
 *      \lambda(E_a+E)\frac{\partial\epsilon_i}{\partial t}\phi_i
 *      + E_a\epsilon_i\phi_i - \sigma_i\phi_i
 *      + \lambda\frac{\partial \sigma}{\partial t}\phi_i
 *  \right\}\phi_j dx = 0
 *  \f]
 */

#include <stdlib.h>
#include <math.h>

#include "fe-solver.h"
#include "matrix.h"
#include "common.h"
#include "deformation.h"

#define EA(X, T) RozziEa((X), (T))
#define E1(X, T) RozziEavg((X), (T))
#define LAMBDA(X, T) RozziLambda((X), (T))
#define STRESS(X, T) EffPorePress((X), (T))

double RozziEa(double M, double T)
{
    double Ea = 68.18*(1/(1+exp((M-250.92*exp(-0.0091*T))/2.19))+0.078);
    return Ea * 1e6;
}

double RozziEavg(double M, double T)
{
    double E1 = 20.26*exp(-0.0802*(M+0.0474*T-14.238)),
           E2 = 2.484 + 6.576/(1+exp((M-19.36)/0.848));
    return (E1+E2)/2 * 1e6;
}

double RozziLambda(double M, double T)
{
    double l1 = 7,
           l2 = 110;
    return (l1+l2)/2;
}

double ResSolid(struct fe1d *p, matrix *guess, Elem1D *elem,
                double x, int f1, int f2)
{
    double T = TINIT,
           C = 0,
           Ci, Ea, value;
    int i;
    basis *b;
    b = p->b;

    solution *s;
    s = CreateSolution(p->t, p->dt, guess);

    for(i=0; i<b->n; i++) {
        Ci = EvalSoln1D(p, CVAR, elem, s, valV(elem->points, i));
        Ci = uscaleTemp(p->chardiff, Ci);

        C += C * b->phi[i](x);
    }

    Ea = EA(C, T);
    value = Ea
            * b->phi[f1](x)*b->phi[f2](x) / IMap1D(p, elem, x);
    free(s);
    return value;
}

double ResDtSolid(struct fe1d *p, matrix *guess, Elem1D *elem,
                  double x, int f1, int f2)
{
    double T = TINIT,
           C = 0,
           Ci, Ea, E, lambda, sigma, term1;
    int i;
    basis *b;
    b = p->b;

    solution *s;
    s = CreateSolution(p->t, p->dt, guess);

    for(i=0; i<b->n; i++) {
        Ci = EvalSoln1D(p, CVAR, elem, s, valV(elem->points, i));
        Ci = uscaleTemp(p->chardiff, Ci);

        C += Ci;
    }

    Ea += EA(C, T);
    E += E1(C, T);
    lambda += LAMBDA(C, T);
    sigma += STRESS(C, T);

    term1 = lambda*(Ea+E) * b->phi[f1](x)*b->phi[f2](x) / IMap1D(p, elem, x);
    free(s);
    return term1;
}

double ResFSolid(struct fe1d *p, matrix *guess, Elem1D *elem, double x, int f1, int a)
{
    double T = TINIT,
           C = 0,
           Cp = 0,
           Ci, Cpi, lambda, sigma, sigmap, DsigmaDt;
    int i;
    solution *s, *sp;
    basis *b;

    b = p->b;
    s = CreateSolution(p->t, p->dt, guess);
    sp = FetchSolution(p, p->t-1);

    for(i=0; i<b->n; i++) {
        Ci = EvalSoln1D(p, CVAR, elem, s, valV(elem->points, i));
        Ci = uscaleTemp(p->chardiff, Ci);

        if(p->t > 1) {
            Cpi = EvalSoln1D(p, CVAR, elem, sp, valV(elem->points, i));
            Cpi = uscaleTemp(p->chardiff, Cpi);
        } else {
            Cpi = Ci;
        }

        C += Ci;
        Cp += Cpi;
    }

    lambda += LAMBDA(C, T);
    sigma += STRESS(C, T);
    sigmap += STRESS(Cp, T);

    DsigmaDt = (sigma-sigmap)/p->dt;
    free(s);
    return sigma - lambda*DsigmaDt;
}

