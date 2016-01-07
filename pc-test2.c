#include "slab-drying/solid/deformation.h"
#include "matrix.h"
#include "material-data.h"

choi_okos *comp_global;

int main(int argc, char *argv[])
{
    matrix *out;
    vector *Xdb,
           *Pc,
           *PcJ0,
           *PcJ,
           *PcG0;
    double T = 313;
    int i, n = 300;

    Xdb = linspaceV(0.05, .4, n);
    Pc = CreateVector(n);
    PcJ0 = CreateVector(n);
    PcJ = CreateVector(n);
    PcG0 = CreateVector(n);
    for(i=0; i<n; i++) {
        setvalV(Pc, i, EffPorePress(valV(Xdb, i), T));
        setvalV(PcJ0, i, EffPorePressExp(valV(Xdb, i), T));
        setvalV(PcJ, i, EffPorePressExpJ(valV(Xdb, i), T));
        setvalV(PcG0, i, EffPorePressExpG0(valV(Xdb, i), T));
    }

    out = CatColVector(5, Xdb, Pc, PcJ, PcJ0, PcG0);
    mtxprntfilehdr(out, "stress-test.csv", "Xdb, Pc, PcJ, PcJ0, PcG0\n");

    return 0;
}

