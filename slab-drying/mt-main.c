/**
 * @file mt-main.c
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "matrix.h"
#include "fe-solver.h"
#include "material-data.h"

#include "common.h"
#include "deformation.h"

choi_okos *comp_global;

int main(int argc, char *argv[])
{
    Mesh1D *mesh;
    basis *b;
    matrix *IC_mass;
    struct fe1d* problem;
    scaling_ht scale_mass;

    solution *s;
    FILE *FPnu;
    double C, x, dt = 1e-3;
    int i, tfinal;
    char *outfile;

    comp_global = CreateChoiOkos(0, 0, 0, 1, 0, 0, 0);
    scale_mass = SetupScaling(DIFF(CINIT, TINIT), CINIT, CAMB, THICKNESS, DIFF(CINIT, TINIT), KC_CONV*1e10);
    //scale_mass = SetupScaling(1, CINIT, CAMB, 1, 1, HCONV);

    /* Make a linear 1D basis */
    b = MakeLinBasis(1);

    /* Create a uniform mesh */
    mesh = GenerateUniformMesh1D(b, 0.0, scaleLength(scale_mass, THICKNESS), 10);

    //tfinal = floor(scaleTime(scale_mass, 72000)/.01);
    tfinal = 3*10*floor(scaleTime(scale_mass, 7200)/dt);
    //tfinal = floor(scaleTime(scale_mass, 1080000/3)/.01);
    printf("tf = %d\n", tfinal);
    problem = CreateFE1D(b, mesh,
                         &CreateDTimeMatrix,
                         &CreateElementMatrix,
                         &CreateElementLoad,
                         &ApplyAllBCs,
                         tfinal);
#ifdef SHRINKAGE_DIFFEQ
    problem->nvars = 4; /* Number of simultaneous PDEs to solve */
#else
    problem->nvars = 1;
#endif

    problem->dt = dt; /* Dimensionless time step size */
    problem->chardiff = scale_mass;

    /* Set the initial temperature */
    IC_mass = GenerateInitCondConst(problem, CVAR, scaleTemp(problem->chardiff, CINIT));

    /* Apply the initial condition to the problem and set up the transient
     * solver. */
    FE1DTransInit(problem, IC_mass);
    FE1DInitAuxSolns(problem, 1);

    FPnu = fopen("poisson.csv", "w");
    fprintf(FPnu, "t,x,Xdb,nu\n");
    while(problem->t<problem->maxsteps) {
        //printf("Step %d of %d\r", problem->t, problem->maxsteps);
        printf("Step %d of %d\r", problem->t, problem->maxsteps);
        fflush(stdout);
        NLinSolve1DTransImp(problem, NULL);

#ifdef SHRINKAGE_DIFFEQ
        //problem->t -= 1;
        s = FetchSolution(problem, problem->t-1);
        if(problem->t > 5) {
            problem->dt = StepSize(problem, PredictSolnO2(problem), s->val);
            printf("New dt = %g\n", StepSize(problem, PredictSolnO2(problem), s->val));
        } else {
            printf("\n");
        }
        //problem->t += 1;
#endif

#ifdef SHRINKAGE
        if(problem->t-1 > 0) {
            problem->mesh =
                MoveMeshF(problem, problem->mesh->orig,
                          problem->t-1, &DeformGradPc);

            //PrintVector(problem->mesh->nodes);
            for(i=5; i<6; i++) {
            //for(i=0; i<len(mesh->nodes); i++) {
                x = valV(mesh->nodes, i);
                s = FetchSolution(problem, problem->t-1);
                C = EvalSoln1DG(problem, CVAR, s, x, 0);
                fprintf(FPnu, "%g,%g,%g,%g\n",
                        uscaleTime(problem->chardiff, problem->t-1),
                        uscaleLength(problem->chardiff, x),
                        uscaleTemp(problem->chardiff, C),
                        //FindPoisson(problem, x, problem->t-1));
                        Porosity(problem, x, problem->t-1));
            }
        }
#endif
    }
    fclose(FPnu);

    PrintScalingValues(problem->chardiff);

    //PrintSolution(problem, 1);
    printf("Solution at t = %g (s):\n", uscaleTime(problem->chardiff, problem->t*problem->dt));
    PrintSolution(problem, problem->t-1);
/*
    CSVOutFixedNode2(problem, 0, "output00.csv");
    CSVOutFixedNode2(problem, 1, "output01.csv");
    CSVOutFixedNode2(problem, 2, "output02.csv");
    CSVOutFixedNode2(problem, 3, "output03.csv");
    CSVOutFixedNode2(problem, 4, "output04.csv");
    CSVOutFixedNode2(problem, 5, "output05.csv");
    CSVOutFixedNode2(problem, 6, "output06.csv");
    CSVOutFixedNode2(problem, 7, "output07.csv");
    CSVOutFixedNode2(problem, 8, "output08.csv");
    CSVOutFixedNode2(problem, 9, "output09.csv");
    CSVOutFixedNode2(problem, 10, "output10.csv");
*/
    outfile = (char*) calloc(sizeof(char), 30);

    sprintf(outfile, "OutAvg-%gK.csv", (double) TINIT);
    CSVOutAvg(problem, CVAR, outfile);
    CSVOutFixedNodeDiff(problem, 10, "OutDx10.csv");
    sprintf(outfile, "OutProfile-%gK.csv", (double) TINIT);
    CSVOutProfiles(problem, 15, outfile);

    free(outfile);

    PrintVector(problem->mesh->orig->nodes);
    PrintVector(problem->mesh->nodes);

    /* Clean up */
    DestroyFE1D(problem);
    DestroyChoiOkos(comp_global);

    return 0;
}

