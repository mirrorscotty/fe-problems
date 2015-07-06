/**
 * @file output.c
 * Spit out simulation results into a CSV file!
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "fe-solver.h"
#include "material-data.h"
#include "common-kelvin.h"

extern choi_okos *comp_global;

/**
 * @brief Output bunches of data for a single node
 *
 * Spits out the time, temperature, etc. for one node n a CSV file. This needs
 * to be changed, depending on whether the freezing model or the sterilization
 * model is built.
 * @param p The problem with the data in it
 * @param row The number of the row (node) to output data for
 * @param filename The name of the file to spit stuff out into
 */
void CSVOutFixedNodeDiff(struct fe1d *p, int row, char *filename)
{
    int i;
    FILE *fp;
    solution *s, *sp;
    double C, X, dC, dx, dt, M, J, rhot;
    choi_okos *comp_dry, *comp_wet;

    row = row * p->nvars;

    fp = fopen(filename, "w+");

    if(!fp) {
        fprintf(stderr, "Unable to open file for writing.");
        return;
    }

    comp_dry = CreateChoiOkos(PASTACOMP);

    /* Print out the column headers */
    //fprintf(fp, "Time,Position,Temperature,Density,HeatCapacity,ThermalConductivity\n");
    fprintf(fp, "Time,DimensionlessTime,Position,WaterFlux,MassFlux\n");

    /* Print out the values */
    for(i=0; i<p->t; i++) {
        s = FetchSolution(p, i);
        if(i<1)
            sp = s;
        else
            sp = FetchSolution(p, i-1);

        C = uscaleTemp(p->chardiff, val(s->val, row+CVAR, 0));
        //dC = EvalDSoln1DG(p, 0, s, valV(p->mesh->orig->nodes, row), 0);
        dC = uscaleTemp(p->chardiff, val(s->val, row+CVAR, 0))
                - uscaleTemp(p->chardiff, val(s->val, row-p->nvars+CVAR, 0));
        dx = uscaleLength(p->chardiff, valV(s->mesh->nodes, row/p->nvars))
                - uscaleLength(p->chardiff, valV(s->mesh->nodes, row/p->nvars-1));
        J = -1*DiffCh10(C, TINIT)*dC/dx;
        comp_wet = AddDryBasis(comp_dry, C);
        rhot = rho(comp_wet, TINIT);
        DestroyChoiOkos(comp_wet);
        dt = uscaleTime(p->chardiff, CurrentTime(p, i))
                - uscaleTime(p->chardiff, CurrentTime(p, i-1));
        M = rhot*dx/dt;

        X = valV(s->mesh->nodes, row);

        fprintf(fp, "%g,%g,%g,%g,%g\n",
                uscaleTime(p->chardiff, CurrentTime(p,i)), CurrentTime(p,i), X, J, M);
    }
    fprintf(fp, "\n");

    DestroyChoiOkos(comp_dry);
    fclose(fp);

    return;
}

/**
 * @brief Output bunches of data for a single node
 *
 * Spits out the time, temperature, etc. for one node n a CSV file. This needs
 * to be changed, depending on whether the freezing model or the sterilization
 * model is built.
 * @param p The problem with the data in it
 * @param row The number of the row (node) to output data for
 * @param filename The name of the file to spit stuff out into
 */
void CSVOutAvg(struct fe1d *p, int var, char *filename)
{
    int i;
    FILE *fp;
    double C, u;
    solution *s;

    fp = fopen(filename, "w+");

    if(!fp) {
        fprintf(stderr, "Unable to open file for writing.");
        return;
    }

    /* Print out the column headers */
    fprintf(fp, "Time,Concentration,Displacement\n");

    /* Print out the values */
    for(i=0; i<p->t; i++) {
        s = FetchSolution(p, i);
        //u = valV(s->mesh->nodes, len(s->mesh->nodes)-1);
#ifdef SUVAR
        u = EvalDSoln1DG(p, SUVAR, s, 1.0, 0);
#else
        u = 0;
#endif

        C = AvgSoln1DG(p, i, var);

        fprintf(fp, "%g,%g,%g\n",
                uscaleTime(p->chardiff, CurrentTime(p, i)),
                uscaleTemp(p->chardiff, C),
                u);
    }
    fprintf(fp, "\n");

    fclose(fp);

    return;
}

void CSVOutProfiles(struct fe1d *p, int n, char *filename)
{
    int dt = floor(p->t/n),
        i, j;
    matrix *data, *tmp1, *tmp2;
    FILE *fp;
    char *header, *hdrtmp;
    solution *s;

    /* Get the first solution... */
    s = FetchSolution(p, 0);
    /* and make the first two columns of the matrix with it. */
    tmp1 = CatColVector(1, s->mesh->nodes);
    data = AugmentMatrix(tmp1, s->val);
    DestroyMatrix(tmp1);

    /* Assume we need an average of 15 characters per profile for a header */
    header = (char*) calloc(sizeof(char), n*50);
    sprintf(header, ",t=%d", 0);

    for(i=1; i<n; i++) {
        /* Get the next solution */
        s = FetchSolution(p, i*dt);
        /* Add a column for the x-coordinates */
        tmp1 = CatColVector(1, s->mesh->nodes);
        tmp2 = AugmentMatrix(data, tmp1);
        DestroyMatrix(tmp1);
        DestroyMatrix(data);
        /* And then for the concentrations */
        data = AugmentMatrix(tmp2, s->val);
        DestroyMatrix(tmp2);

        /* Then do the header */
        hdrtmp = (char*) calloc(sizeof(char), 50);
        sprintf(hdrtmp, ",,t=%g", uscaleTime(p->chardiff, i*dt));
        strcat(header, hdrtmp);
        free(hdrtmp);
    }

    strcat(header, "\n");

    fp = fopen(filename, "w+");
    fprintf(fp, header);
    free(header);
    for(i=0; i<nRows(data); i++) {
        for(j=0; j<nCols(data)-1; j++) {
            if(j%2)
                fprintf(fp, "%g,", uscaleTemp(p->chardiff, val(data, i, j)));
            else
                fprintf(fp, "%g,", val(data, i, j));
        }
        fprintf(fp, "%g\n",
                uscaleTemp(p->chardiff, val(data, i, nCols(data)-1)));
    }
    fclose(fp);
    DestroyMatrix(data);

    return;
}

