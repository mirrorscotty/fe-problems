#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "matrix.h"
#include "fe-solver.h"

#include "material-data.h"

#include "heat-transfer.h"
#include "diffusion.h"
#include "common-mod.h"
#include "lin-genmaxwell.h"
#include "deformation.h"

extern choi_okos *comp_global;

/**
 * Calculate an average diffusivity for the entire drying process. The
 * temperature is taken to be the initial temperature, and the average between
 * the initial and final moisture content is used.
 * @param Xdb Moisture content (Not used)
 * @param T Temperature (Not used)
 * @returns Diffusivity [m^2/s]
 */
double DiffAvg(double Xdb, double T)
{
    return DiffCh10((CINIT+CAMB)/2, TINIT);
}

/**
 * This takes the function that calculates the residual and integrates it
 * to form the element-level matrix used in solving the problem. This whole
 * function is given to the finite element solver and run repeatedly to generate
 * the global matrix given to the matrix solver.
 * @param p Finite element problem structure
 * @param elem Element to calculate the matrix for
 * @param guess Estimated value of the solution at the current time step
 *
 * @returns An n*v x n*v matrix containing the integrated values of the
 *      residual, where "n" is the number of basis functions used. For a
 *      linear basis, n = 2. The value for "v" is the number of PDEs being
 *      solved simultaneously.
 */
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
            value = quad1d3generic(p, guess, elem, &ResMass, i/v, j/v);
            setval(m, value, i+CVAR, j+CVAR);

            value = quad1d3generic(p, guess, elem, &ResSolid_dTde, i/v, j/v);
            setval(m, value, i+STVAR, j+STVAR);
            value = quad1d3generic(p, guess, elem, &ResSolid_dTdq1, i/v, j/v);
            setval(m, value, i+STVAR, j+SP1VAR);
            value = quad1d3generic(p, guess, elem, &ResSolid_dTdq2, i/v, j/v);
            setval(m, value, i+STVAR, j+SP2VAR);

            value = quad1d3generic(p, guess, elem, &ResSolid_zero, i/v, j/v);
            setval(m, value, i+SP1VAR, j+STVAR);
            value = quad1d3generic(p, guess, elem, &ResSolid_dP1dq1, i/v, j/v);
            setval(m, value, i+SP1VAR, j+SP1VAR);
            value = quad1d3generic(p, guess, elem, &ResSolid_zero, i/v, j/v);
            setval(m, value, i+SP1VAR, j+SP2VAR);

            value = quad1d3generic(p, guess, elem, &ResSolid_zero, i/v, j/v);
            setval(m, value, i+SP2VAR, j+STVAR);
            value = quad1d3generic(p, guess, elem, &ResSolid_zero, i/v, j/v);
            setval(m, value, i+SP2VAR, j+SP1VAR);
            value = quad1d3generic(p, guess, elem, &ResSolid_dP2dq2, i/v, j/v);
            setval(m, value, i+SP2VAR, j+SP2VAR);
        }
    }

    return m;
}

/**
 * This takes the function that calculates the element-level coefficient matrix
 * multiplying \f$\frac{\partial T_i}{\partial t}\f$. As with the
 * CreateElementMatrix function, this is passed to the FEM solver and run
 * repeatedly to generate the global matrix.
 * @param p Finite element problem structure
 * @param elem Element to calculate the matrix for
 * @param guess Estimated value of the solution at the current time step
 *
 * @returns An n*v x n*v matrix containing the integrated values of the
 *      portion of the residual multiplying
 *      \f$\frac{\partial T_i}{\partial t}\f$, where "n" is the number of basis
 *      functions used. For a linear basis, n = 2. The value for "v" is the
 *      number of PDEs being solved simultaneously.
 * @see CreateElementMatrix
 */
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
            value = quad1d3generic(p, guess, elem, &ResDtMass, i/v, j/v);
            setval(m, value, i+CVAR, j+CVAR);

            value = quad1d3generic(p, guess, elem, &ResSolid_zero, i/v, j/v);
            setval(m, value, i+STVAR, j+STVAR);
            value = quad1d3generic(p, guess, elem, &ResSolid_zero, i/v, j/v);
            setval(m, value, i+STVAR, j+SP1VAR);
            value = quad1d3generic(p, guess, elem, &ResSolid_zero, i/v, j/v);
            setval(m, value, i+STVAR, j+SP2VAR);

            value = quad1d3generic(p, guess, elem, &ResDtSolid_dP1de, i/v, j/v);
            setval(m, value, i+SP1VAR, j+STVAR);
            value = quad1d3generic(p, guess, elem, &ResDtSolid_dP1dq1, i/v, j/v);
            setval(m, value, i+SP1VAR, j+SP1VAR);
            value = quad1d3generic(p, guess, elem, &ResSolid_zero, i/v, j/v);
            setval(m, value, i+SP1VAR, j+SP2VAR);

            value = quad1d3generic(p, guess, elem, &ResDtSolid_dP2de, i/v, j/v);
            setval(m, value, i+SP2VAR, j+STVAR);
            value = quad1d3generic(p, guess, elem, &ResSolid_zero, i/v, j/v);
            setval(m, value, i+SP2VAR, j+SP1VAR);
            value = quad1d3generic(p, guess, elem, &ResDtSolid_dP2dq2, i/v, j/v);
            setval(m, value, i+SP2VAR, j+SP2VAR);
        }
    }

    return m;
}

/**
 * Creates an integrated, element-level load vector (source term + boundary
 * conditions). The boundary conditions are added separately, and the source
 * term for this equation is equal to zero, so this simply returns a column
 * matrix of zeros of the appropriate size.
 */
matrix* CreateElementLoad(struct fe1d *p, Elem1D *elem, matrix *guess) {
    basis *b;
    double value;
    int v = p->nvars, i;
    matrix *m;

    b = p->b;

    m = CreateMatrix(b->n*v, 1);

    for(i=0; i<b->n*v; i+=v) {
        value = quad1d3generic(p, guess, elem, &ResFSolid_T, i/v, 0);
            printf("Value = %g\n", value);
        setval(m, value, i+STVAR, 0);
    }

    return m;
}

/**
 * Returns true if the specified row (node) is on the right-most boundary of the
 * domain.
 * @param p Finite element problem structure
 * @param row Row to check to see if it's on the boundary. Since this is a 1D
 *      problem, the row number corresponds to the node number.
 * @returns True if the node is on the right-hand boundary and false otherwise.
 */
int IsOnRightBoundary(struct fe1d *p, int row)
{
    if(row == len(p->mesh->nodes)-1)
        return 1;
    else
        return 0;
}

/**
 * Returns true if the specified row (node) is on the left-most boundary of the
 * domain.
 * @param p Finite element problem structure
 * @param row Row to check to see if it's on the boundary. Since this is a 1D
 *      problem, the row number corresponds to the node number.
 * @returns True if the node is on the left-hand boundary and false otherwise.
 */
int IsOnLeftBoundary(struct fe1d *p, int row)
{
    if(row == 0)
        return 1;
    else
        return 0;
}

/**
 * Take the tests functions that define where the boundaries are and the
 * functions for boundary conditions and apply them to the problem. This is
 * given to the FEM solver and run each time the global matrices are
 * recalculated.
 * @param p Finite element problem structure
 */
void ApplyAllBCs(struct fe1d *p)
{
#ifdef TVAR
    double Bi = BiotNumber(p->charvals);
#endif
#ifdef CVAR
    double Bim = BiotNumber(p->chardiff);
#endif

    /* BC at x=L:
     * This approximates any Biot number larger than 100 as Bi->infty. This is a
     * good approximation for this problem since it results in the outside of
     * the can reaching the external temperature incredibly quickly (<1sec). */
#ifdef TVAR
    if(Bi<100.00)
        ApplyNaturalBC1D(p, TVAR, &IsOnRightBoundary, &ConvBCHeat);
    else
        ApplyEssentialBC1D(p, TVAR, &IsOnRightBoundary, &ExternalTemp);
#endif

    /* Do the same for the mass transfer boundary condition. */
#ifdef CVAR
    if(Bim<100.00)
        ApplyNaturalBC1D(p, CVAR, &IsOnRightBoundary, &ConvBCMass);
    else
        ApplyEssentialBC1D(p, CVAR, &IsOnRightBoundary, &ExternalConc);
#endif

#ifdef VAP_MODEL
    if(Bim<100.00)
        ApplyNaturalBC1D(p, PVAR, &IsOnRightBoundary, &ConvBCVap);
    else
        ApplyEssentialBC1D(p, PVAR, &IsOnRightBoundary, &ExternalConc);
#endif
    return;
}

