/**
 * @file kinematics.c
 *
 * Differential equation to calculate velocity
 * \f[
 * \underline{\underline{\dot{F}}} = \frac{\partial}{\partial t}[1-\epsilon]
 *  = \frac{\partial}{\partial\underline{X}}[ \underline{V}(\underline{X}, t)]
 *  = \frac{\partial}{\partial\underline{x}}[ \underline{v}(\underline{x}, t)]
 * \f]
 * Weak Form:
 * \f[
 * R_i^{v} = \frac{\partial\epsilon}{\partial t}
 * + v_i\frac{\partial\phi_i}{\partial x}
 * \f]
 *
 * Displacement:
 * \f[
 * \nabla\underline{u} = \underline{\underline{F}}-\underline{\underline{I}}
 * \f]
 * \f[
 * R_i^u = u_i \frac{\partial\phi_i}{\partial x}\phi_j = \epsilon_i\phi_i\phi_j
 * \f]
 */

#include "kinematics.h"
#include "fe-solver.h"
#include "../common-kelvin.h"

double ResDtSolid_dvde(struct fe1d *p, matrix *guess, Elem1D *elem,
                       double x, int f1, int f2)
{
    return p->b->phi[f1](x) * p->b->phi[f2](x) / IMap1D(p, elem, x);
}

double ResSolid_dvdv(struct fe1d *p, matrix *guess, Elem1D *elem,
                     double x, int f1, int f2)
{
    return p->b->phi[f1](x) * p->b->dphi[f2](x);
}

double ResSolid_dudu(struct fe1d *p, matrix *guess, Elem1D *elem,
                     double x, int f1, int f2)
{
    return p->b->phi[f1](x) * p->b->dphi[f2](x);// * IMap1D(p, elem, x);
}

double ResSolid_dude(struct fe1d *p, matrix *guess, Elem1D *elem,
                     double x, int f1, int f2)
{
    return p->b->phi[f1](x) * p->b->phi[f2](x) / IMap1D(p, elem, x);
}
/*
double ResFSolid_u(struct fe1d *p, matrix *guess, Elem1D *elem,
                   double x, int f1, int f2)
{
    int i;
    double ei = 0, e = 0;
    basis *b;
    solution *s;

    b = p->b;

    s = CreateSolution(p->t, p->dt, guess);
    for(i=0; i<b->n; i++) {
        ei = EvalSoln1D(p, STVAR, elem, s, valV(elem->points, i));
        e += ei * b->phi[i](x);
    }
    free(s);

    e = .05;
    return e * b->phi[f1](x) / IMap1D(p, elem, x);
}
*/
