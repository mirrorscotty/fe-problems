#include "drying-gui.h"

fe1d* setupDomain()
{
    fe1d *problem;

    int NNodes = 10;
    double Deltat = 0.01;
    int NTimeSteps = 1000;

    double Tc = 300;
    double Lc = .1;
    double h = 300;

    double Mc = 1;

    Mesh1D *mesh;
    basis *b;
    matrix *IC;

    scaling_ht temp;
    scaling_ht diff;

    temp = SetupScaling(1e-5, Tc, Lc, 1e-2, h);
    diff = SetupScaling(1e-8, Mc, Lc, 1, 1);

    b = MakeLinBasis(1);

    mesh = GenerateUniformMesh1D(b, 0.0, 1.0, NNodes-1);

    problem = CreateFE1D(b, mesh,
                         &CreateDTimeMatrix,
                         &CreateElementMatrix,
                         &CreateElementLoad,
                         &ApplyAllBCs,
                         NTimeSteps);

    problem->charvals = temp;
    problem->chardiff = diff;

    problem->nvars = 2;

    problem->dt = scaleTime(charvals, Deltat);

    IC = GenerateInitCondConst(problem, 0, 1);
    IC = GenerateInitCondConst(problem, 1, 1);



