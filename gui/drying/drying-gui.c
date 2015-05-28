#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "matrix.h"
#include "basis.h"
#include "integrate.h"
#include "mesh1d.h"
#include "isoparam.h"
#include "finite-element1d.h"
#include "auxsoln.h"

#include "material-data/choi-okos/choi-okos.h"

#include "drying-gui.h"

extern double EaA, EaB, AA, AB;
double h = 5;

/* The following two functions are for heat conduction. */
/* Creates the Jacobian and helps solve for the current time step */
double ResT(struct fe1d *p, matrix *guess, Elem1D *elem, double x, int f1, int f2)
{
    double value = 0;
    basis *b;
    b = p->b;

    /* This is just to fetch the value of T */
    double T;
    solution *s;
    s = FetchSolution(p, p->t-1);
    if(!s)
        T = p->charvals.Tc;
    else
        T = EvalSoln1D(p, 0, elem, s, x);
   
    /* Normally, this should be multiplied by the thermal diffusivity. However,
     * we will assume that alpha is constant and, because of dimensionless
     * groups, this is all done later. */
    value  = b->dphi[f1](x) * b->dphi[f2](x);
    value *= IMap1D(p, elem, x);
    value *= (alpha(T)/p->charvals.alpha);
    value *= IMapCyl1D(p, elem, x);

    //printf("alpha_c = %g, alpha = %g, ", p->charvals.alpha, alpha(T));
    //printf("alpha_c/alpha = %g\n", p->charvals.alpha/alpha(T));

    return value;
}

/* Calculate the coefficient matrix for the time derivative unknowns */
double ResTDt(struct fe1d *p, matrix *guess, Elem1D *elem, double x, int f1, int f2)
{
    double value;
    basis *b;
    b = p->b;

    value = b->phi[f1](x) * b->phi[f2](x) * 1/IMap1D(p, elem, x);
    value *= IMapCyl1D(p, elem, x);
//    value *= -1;

    return value;
}

/* The following two functions are for diffusion. */
/* Creates the Jacobian and helps solve for the current time step */
double ResM(struct fe1d *p, matrix *guess, Elem1D *elem, double x, int f1, int f2)
{
    double value = 0;
    basis *b;
    b = p->b;

    /* This is just to fetch the value of T */
    double T;
    solution *s;
    s = FetchSolution(p, p->t-1);
    if(!s)
        T = p->charvals.Tc;
    else
        T = EvalSoln1D(p, 0, elem, s, x);
   
    /* Normally, this should be multiplied by the diffusion coef However,
     * we will assume that D is constant and, because of dimensionless
     * groups, this is all done later. */
    value  = b->dphi[f1](x) * b->dphi[f2](x);
    value *= IMap1D(p, elem, x);
    value *= (alpha(T)/p->charvals.alpha);
    value *= IMapCyl1D(p, elem, x);

    //printf("alpha_c = %g, alpha = %g, ", p->charvals.alpha, alpha(T));
    //printf("alpha_c/alpha = %g\n", p->charvals.alpha/alpha(T));

    return value;
}

/* Calculate the coefficient matrix for the time derivative unknowns */
double ResMDt(struct fe1d *p, matrix *guess, Elem1D *elem, double x, int f1, int f2)
{
    double value;
    basis *b;
    b = p->b;

    value = b->phi[f1](x) * b->phi[f2](x) * 1/IMap1D(p, elem, x);
    value *= IMapCyl1D(p, elem, x);
//    value *= -1;

    return value;
}

double ResZero(struct fe1d *p, matrix *guess, Elem1D *elem, double x, int f1, int f2)
{
    return 0.0;
}

matrix* CreateElementMatrix(struct fe1d *p, Elem1D *elem, matrix *guess)
{
    basis *b;
    b = p->b;
    
    int v = p->nvars;
    
    int i, j;
    double value = 0;
    matrix *m;
    
    m = CreateMatrix(b->n*v, b->n*v);
    
    for(i=0; i<b->n*v; i+=v) {
        for(j=0; j<b->n*v; j+=v) {
            value = quad1d3generic(p, guess, elem, &ResT, i/v, j/v);
            setval(m, value, i, j);

            value = quad1d3generic(p, guess, elem, &ResZero, i/v, j/v);
            setval(m, value, i+1, j);

            value = quad1d3generic(p, guess, elem, &ResZero, i/v, j/v);
            setval(m, value, i, j+1);

            value = quad1d3generic(p, guess, elem, &ResM, i/v, j/v);
            setval(m, value, i+1, j+1);
        }
    }

    return m;
}

/* Create the coefficient matrix for the time derivatives */
matrix* CreateDTimeMatrix(struct fe1d *p, Elem1D *elem, matrix *guess) {
    basis *b;
    b = p->b;
    
    int v = p->nvars;
    
    int i, j;
    double value = 0;
    matrix *m;
    
    m = CreateMatrix(b->n*v, b->n*v);
    
    for(i=0; i<b->n*v; i+=v) {
        for(j=0; j<b->n*v; j+=v) {
            value = quad1d3generic(p, guess, elem, &ResTDt, i/v, j/v);
            setval(m, value, i, j);

            value = quad1d3generic(p, guess, elem, &ResZero, i/v, j/v);
            setval(m, value, i+1, j);

            value = quad1d3generic(p, guess, elem, &ResZero, i/v, j/v);
            setval(m, value, i, j+1);

            value = quad1d3generic(p, guess, elem, &ResMDt, i/v, j/v);
            setval(m, value, i+1, j+1);
        }
    }

    return m;
}

/* Create the load vector... thing */
matrix* CreateElementLoad(struct fe1d *p, Elem1D *elem, matrix *guess) {
    basis *b;
    b = p->b;
    
    int v = p->nvars;
    
    matrix *m;
    
    m = CreateMatrix(b->n*v, 1);

    return m;
}

/* Returns true if the specified row (node) is on the right-most boundary of the
 * domain. */
int IsOnRightBoundary(struct fe1d *p, int row)
{
    double width = p->mesh->x2 - p->mesh->x1;
    double x = valV(p->mesh->nodes, row/p->nvars);
    
    if(fabs(x - width) < 1e-5)
        return 1;
    else
        return 0;
}

/* Same as above, only for the left boundary */
int IsOnLeftBoundary(struct fe1d *p, int row)
{
    double x = valV(p->mesh->nodes, row/p->nvars);
  
    if(fabs(x) < 1e-5)
        return 1;
    else
        return 0;
}

double ExternalTemp(struct fe1d *p, int row)
{
    return scaleTemp(p->charvals, T_ext(uscaleTime(p->charvals, p->t*p->dt)));
}


/* The way this function is implemented is probably not mathematically accurate.
 * Strictly speaking, for the implicit solver, the temperature fetched should be
 * the temperature at the next time step, not the one that we already have the
 * solution for. The way it is now should be good enough (tm).*/
double ConvBC(struct fe1d *p, int row)
{
    double Tinf = scaleTemp(p->charvals, T_ext(uscaleTime(p->charvals, p->dt*p->t)));
    double T;
    double Bi = BiotNumber(p->charvals);

    /* Fetch the value of T from the previous solution. If this is being
     * applied at the initial condition, then simply return 0.*/
    solution *s;
    s = FetchSolution(p, p->t-1);
    if(s) {
        T = val(s->val, row, 0);
        if(T==0)
            return 0;
        return -Bi*(T-Tinf);
    } else {
        return 0;
    }
}

/* Diffusion boundary conditions */
double ExternalMoisture(struct fe1d *p, int row)
{
    return scaleTemp(p->chardiff, T_ext(uscaleTime(p->chardiff, p->t*p->dt)));
}

void ApplyAllBCs(struct fe1d *p)
{
    double Bi = BiotNumber(p->charvals); 
    
    /* BC at x=L
     * This approximates any Biot number larger than 100 as Bi->infty. This is a
     * good approximation for this problem since it results in the outside of
     * the can reaching the external temperature incredibly quickly (<1sec). */
    if(Bi<100.00)
        ApplyNaturalBC1D(p, 0, &IsOnRightBoundary, &ConvBC);
    else
        ApplyEssentialBC1D(p, 0, &IsOnRightBoundary, &ExternalTemp);

    ApplyEssentialBC1D(p, 1, &IsOnRightBoundary, &ExternalMoisture);

    return;
}

