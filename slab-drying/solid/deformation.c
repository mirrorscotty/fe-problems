#include "fe-solver.h"
#include "../common.h"
#include "deformation.h"
#include "material-data.h"
#include "matrix.h"

extern choi_okos *comp_global;

double CreepGina(double t, double T, double X, double P, int deriv)
{
    double J;
    burgerse *b;
    b = CreateBurgersE();
    if(deriv)
        J = DBurgersECreep(b, t, T, X, P);
    else
        J = BurgersECreep(b, t, T, X, P);
    DestroyBurgersE(b);
    return J;
}

double CreepGinaBulk(double t, double T, double X, double P, int deriv)
{
    double nu = .45, /* Poisson ratio */
           J, /* Creep compliance */
           B; /* Bulk compliance */
    burgerse *b;
    b = CreateBurgersE();
    if(deriv)
        J = DBurgersECreep(b, t, T, X, P);
    else
        J = BurgersECreep(b, t, T, X, P);
    DestroyBurgersE(b);

    B = 6*(.5-nu)*J;
    return B;
}

double CreepLaura(double t, double T, double X, double P, int deriv)
{
    double J;
    if(deriv)
        J = DMaxwellCreepConverted(t, T, X);
    else
        J = MaxwellCreepConverted(t, T, X);
    return J;
}

double CreepLaura2(double t, double T, double X, double P, int deriv)
{
    if(t<0.01)
        t=0.01;
    double J;
    if(deriv)
        J = DLLauraCreep(t, T, X, P);
    else
        J = LLauraCreep(t, T, X, P);
    return J;
}

double CreepZhu(double t, double T, double X, double P, int deriv)
{
    double J, Bc = 8.2e-14;
    maxwell *m;
    m = CreateMaxwellZhu();
    if(deriv)
        J = DMaxwellCreep(m, t, T, X);
    else
        J = MaxwellCreep(m, t, T, X);
    DestroyMaxwell(m);
    return J*1e0*Bc;
}

double CreepCummings(double t, double T, double X, double P, int deriv)
{
    double J;
    if(t<1e-2)
        t = 1e-2;
    if(deriv)
        J = DLCummingsCreep(t, T, X, P);
    else
        J = LCummingsCreep(t, T, X, P);
    return J;
}

double CreepLauraL(double t, double T, double X, double P, int deriv)
{
    double J;
    if(deriv)
        J = DLLauraCreep(t, T, X, P);
    else
        J = LLauraCreep(t, T, X, P);
    return J;
}

/**
 * Calculate the value of the deformation gradient at the specified time and
 * spatial coordinates.
 * \f[ \underline{\underline{F}}(\underline{X}, T)\f]
 * This value is calculated from the density Choi-Okos equations, and supplied
 * to the moving mesh function to recalculate the nodal values at each time
 * step.
 * @param p Finite element problem structure
 * @param X Spatial coordinate (global)
 * @param t Time step number
 * @returns Calculated value for the deformation gradient
 */
double DeformationGrad(struct fe1d *p, int i, double X, double t)
{
    solution *s0, *sn;
    double rho0, rhon;
    double T0 = TINIT, Tn = TINIT;

#ifdef CVAR
    double C0, Cn;
    choi_okos *cowet0, *cowetn;
#endif

    s0 = FetchSolution(p, 0);
    sn = FetchSolution(p, t);

#ifdef TVAR
    Tn = uscaleTemp(p->charvals, EvalSoln1DG(p, TVAR, sn, X, 0));
    T0 = uscaleTemp(p->charvals, EvalSoln1DG(p, TVAR, s0, X, 0));
#endif
#ifdef CVAR
    Cn = uscaleTemp(p->chardiff, EvalSoln1DG(p, CVAR, sn, X, 0));
    C0 = uscaleTemp(p->chardiff, EvalSoln1DG(p, CVAR, s0, X, 0));

    cowet0 = AddDryBasis(comp_global, C0);
    cowetn = AddDryBasis(comp_global, Cn);

    rho0 = rho(cowet0, T0);
    rhon = rho(cowetn, Tn);

    DestroyChoiOkos(cowet0);
    DestroyChoiOkos(cowetn);
#else
    rhon = rho(comp_global, Tn);
    rho0 = rho(comp_global, T0);
#endif

    return rho0/rhon;
}

/**
 * Calculate the effective pore pressure based on the initial conditions,
 * temperature, and moisture content.
 */
double EffPorePressDefault(double X, double T)
{
    return pore_press(X, T);
}

/**
 * Calculate the effective pore pressure based on the initial conditions,
 * temperature, and moisture content.
 */
double EffPorePress(double X, double T)
{
    double P = pore_press(X, T),
           P0 = pore_press(CINIT, T),
           Pnet = 0;
    Pnet = P-P0;
    //if(Pnet < 0)
        //Pnet = 0;

    return Pnet;
}

/**
 * Calculate the porosity of the sample at point X (material coordinates) using
 * the assumption that the pores were initially filled with water. Here, the
 * porosity is taken to be the porosity of the solid matrix itself, without
 * considering the volume of either water or air that may be filling the pores.
 * @param p Finite element problem definition
 * @param X Material coordinate [-]
 * @param t Time [-]
 * @returns Porosity [-]
 */
double Porosity(struct fe1d *p, double X, int t)
{
    solution *s;
    double vf, vs, e, rhos, rhow, T=TINIT, phi;
    choi_okos *co;

    /* Calculate the strain at the desired point */
    s = FetchSolution(p, t);
    if(t<1)
        e = 0;
    else
        e = (EvalSoln1DG(p, -1, s, X, 0)-EvalSoln1DG(p, -1, s, X, 1))
            /EvalSoln1DG(p, -1, s, X, 0);


    /* Find the water and pasta densities based on the Choi-Okos equations. */
    co = CreateChoiOkos(WATERCOMP);
    rhow = rho(co, T);
    DestroyChoiOkos(co);
    co = CreateChoiOkos(PASTACOMP);
    rhos = rho(co, T);
    DestroyChoiOkos(co);

    vs = 1/(1+(rhos*CINIT/rhow)); /* Solid volume */
    vf = 1-e; /* Total volume */
    phi = (vf-vs)/(vf); /* Porosity */

    return phi*0.1;
}


/**
 * Unoptimized time integration algorithm to calculate part of the strain.
 */
double _StrainPc(struct fe1d *p, int z, double X, double t)
{
    solution *s;
    double T = TINIT,
           /* TODO: Fix this so that the time step can change */
           tf = uscaleTime(p->chardiff, t*p->dt),
           dt, ti, Xdb, P, e=0;
    int i;

    for(i=1; i<p->t; i++) {
        s = FetchSolution(p, i);
        dt = s->dt;
        /* Time */
        ti = uscaleTime(p->chardiff, i*dt);
        /* Moisture content */
        Xdb = uscaleTemp(p->chardiff, EvalSoln1DG(p, CVAR, s, X, 0));
        /* Calculate pore pressure */
        P = EffPorePress(Xdb, T) * Porosity(p, X, i*dt);
        /* Evaluate the creep function and add the current strain to the total
         */
        e += CREEP(tf-ti, T, Xdb, -1*P, 1) * P * s->dt;
    }
    return e;
}

double _StrainPcOpt(struct fe1d *p, int i, double X, double t)
{
    solution *s, *sa;
    double T = TINIT,
           tf = uscaleTime(p->chardiff, t*p->dt),
           ti, Xdb, P, h;

    sa = 0;

    if(t==1) {
        h=0;
    } else {
        sa = FetchAuxSoln(p, 0, t-1);

        /* Internal strain variable value */
        h = val(sa->val, i, 0);
    }

    s = FetchSolution(p, t);
    /* Time */
    ti = uscaleTime(p->chardiff, t);
    /* Moisture content */
    Xdb = uscaleTemp(p->chardiff, EvalSoln1DG(p, CVAR, s, X, 0));
    /* Calculate pore pressure */
    P = EffPorePress(Xdb, T) * Porosity(p, X, t);
    /* Evaluate the creep function and add the current strain to the total
     */
    h += CREEP(tf-ti, T, Xdb, -1*P, 1) * P * p->dt;

    sa = FetchAuxSoln(p, 0, t);
    if(!sa) {
        InitAuxStep(p, 0, t);
        sa = FetchAuxSoln(p, 0, t);
    }
    setval(sa->val, h, i, 0);

    return h;
}

/**
 * Calculate the value of the strain at the specified time and spatial
 * coordinates.
 * The value here is based on the pore pressure and viscoelastic strain from
 * the specified models.
 * @param p Finite element problem structure
 * @param i Node number
 * @param X Spatial coordinate (global)
 * @param t Time step number
 * @returns Calculated value for the deformation gradient
 */
double StrainPc(struct fe1d *p, int i, double X, double t)
{
#ifndef CVAR
    return 1;
#endif

    solution *s;
    double T = TINIT, /* Initial temperature */
           Xdb, /* Moisture Content */
           tf = uscaleTime(p->chardiff, t*p->dt), /* Final time [s] */
           e = 0, /* Strain [-] */
           P; /* Pore pressure [Pa] */

    /* Integrate the creep function from t=0 to t=tf. */
    //e = _StrainPc(p, X, t);
    e = _StrainPc(p, i, X, t);

    /* Because we used integration by parts, add in the rest of the integration
     * formula. Because we're using inverse Laplace transforms, use a value
     * that's slightly above zero to prevent numerical errors. */
    s = FetchSolution(p, t);
    Xdb = uscaleTemp(p->chardiff, EvalSoln1DG(p, CVAR, s, X, 0));
    P = EffPorePress(Xdb, T) * Porosity(p, X, t);
    e += -1*CREEP(tf, T, Xdb, -1*P, 0) * P;
    s = FetchSolution(p, 0);
    Xdb = uscaleTemp(p->chardiff, EvalSoln1DG(p, CVAR, s, X, 0));
    P = EffPorePress(Xdb, T) * Porosity(p, X, 0);
    e += CREEP(0.01, T, Xdb, -1*P, 0) * P;

    if(e<0)
        return 0;
    else
        return e;
}

/**
 * The creep function is converted to bulk compliance using the formula in
 * Bazang 1975 (assuming constant poisson ratio).
 * @param p Finite element problem definition
 * @param i Node number
 * @param X Lagrangian x-coordinate
 * @param t Time (dimensionless)
 *
 * @returns Deformation Gradient
 */
double DeformGradPc(struct fe1d *p, int i, double X, double t)
{
    return 1 - StrainPc(p, i, X, t)*6*(.5-POISSON);
}

/**
 * Deformation gradient based on the hygroscopic expansion coefficient given in Cummings et al. 1993.
 * \f[
 * \frac{\Delta L}{L_0} = \beta\Delta X_{db}
 * \f]
 * where \f$\beta\f$ is the hygroscopic expansion coefficient.
 */
double DeformGradBeta(struct fe1d *p, int i, double X, double t)
{
#ifndef CVAR
    return 1;
#endif

    solution *s;
    double e = 0,
           beta = BETA,
           C = 0;

    s = FetchSolution(p, t-1);
    C = EvalSoln1DG(p, CVAR, s, X, 0);

    e = beta * (C-CINIT);

    return 1+e;
}

double FindPoisson(struct fe1d *p, int i, double X, double t)
{
    double e = StrainPc(p, i, X, t),
           F = DeformGradBeta(p, i, X, t),
           nu;

    nu = (3*e + F-1)/(6*e);
    //if(nu<0)
        //printf("e = %g, F = %g\n", e, F);

    return nu;
}

