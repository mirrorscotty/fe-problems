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
#include "common.h"

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
    int i, n;
    FILE *fp;
    solution *s;
    double C, X, dC, J, rhot;
    Mesh1D *mesh;
    choi_okos *comp_dry, *comp_wet;

    fp = fopen(filename, "w+");

    if(!fp) {
        fprintf(stderr, "Unable to open file for writing.");
        return;
    }

    comp_dry = CreateChoiOkos(PASTACOMP);

    /* Print out the column headers */
    //fprintf(fp, "Time,Position,Temperature,Density,HeatCapacity,ThermalConductivity\n");
    fprintf(fp, "Time,Position,Concentration,WaterFlux,Density\n");

    /* Print out the values */
    for(i=0; i<p->maxsteps; i++) {
        s = FetchSolution(p, i);

        C = uscaleTemp(p->chardiff, val(s->val, row, 0));
        dC = EvalDSoln1DG(p, 0, s, valV(p->mesh->orig->nodes, row), 0);
        dC = uscaleTemp(p->chardiff, dC);
        J = -1*dC*DiffCh10(C, TINIT);
        comp_wet = AddDryBasis(comp_dry, C);
        rhot = rho(comp_wet, TINIT);
        DestroyChoiOkos(comp_wet);

        mesh = p->mesh;
        for(n=0; n<p->maxsteps-i-1; n++)
            mesh = mesh->prev;
        X = valV(mesh->nodes, row);
        
        fprintf(fp, "%g,%g,%g,%g,%g\n", i*p->dt, X, C, J, rhot);
        //fprintf(fp, "%g,%g,%g\n", rho(comp_global, T), Cp(comp_global, T), k(comp_global, T));
    }
    fprintf(fp, "\n");

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
    for(i=0; i<p->maxsteps; i++) {
        s = FetchSolution(p, i);
        u = valV(s->mesh->nodes, len(s->mesh->nodes)-1);

        C = AvgSoln1DG(p, i, var);

        fprintf(fp, "%g,%g,%g\n",
                uscaleTime(p->chardiff, i*p->dt),
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
    matrix *main, *tmp1, *tmp2;
    FILE *fp;
    char *header, *hdrtmp;
    solution *s;

    /* Get the first solution... */
    s = FetchSolution(p, 0);
    /* and make the first two columns of the matrix with it. */
    tmp1 = CatColVector(1, s->mesh->nodes);
    main = AugmentMatrix(tmp1, s->val);
    DestroyMatrix(tmp1);

    /* Assume we need an average of 15 characters per profile for a header */
    header = (char*) calloc(sizeof(char), n*15);
    sprintf(header, ",t=%d", 0);
    
    for(i=1; i<n; i++) {
        /* Get the next solution */
        s = FetchSolution(p, i*dt);
        /* Add a column for the x-coordinates */
        tmp1 = CatColVector(1, s->mesh->nodes);
        tmp2 = AugmentMatrix(main, tmp1);
        DestroyMatrix(tmp1);
        DestroyMatrix(main);
        /* And then for the concentrations */
        main = AugmentMatrix(tmp2, s->val);
        DestroyMatrix(tmp2);
        
        /* Then do the header */
        hdrtmp = (char*) calloc(sizeof(char), 15);
        sprintf(hdrtmp, ",,t=%g", uscaleTime(p->chardiff, i*dt));
        strcat(header, hdrtmp);
    }

    strcat(header, "\n");

    fp = fopen(filename, "w+");
    fprintf(fp, header);
    for(i=0; i<nRows(main); i++) {
        for(j=0; j<nCols(main)-1; j++) {
            if(j%2)
                fprintf(fp, "%g,", uscaleTemp(p->chardiff, val(main, i, j)));
            else
                fprintf(fp, "%g,", val(main, i, j));
        }
        fprintf(fp, "%g\n",
                uscaleTemp(p->chardiff, val(main, i, nCols(main)-1)));
    }
    fclose(fp);

    return;
}

