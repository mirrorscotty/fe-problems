CC=gcc
CFLAGS=-lm -I. -Imatrix -Imaterial-data -Ife-solver -Wall -g3 -O2
VPATH=other slab-drying problems/viscoelasticity gui solver/mesh solver/ode solver/integration matrix material-data scaling solver output

all: diffusion

doc:
	doxygen DoxyFile
	make -C doc/latex
	cp doc/latex/refman.pdf doc/Reference.pdf

2dlaplace: 2dlaplace.o fe-solver.a material-data.a matrix.a
	$(CC) -o $@ $^ $(CFLAGS)

ce675p1: ce675p1.o fe-solver.a material-data.a matrix.a
	$(CC) -o $@ $^ $(CFLAGS)

spheroid: spheroid.o fe-solver.a material-data.a matrix.a
	$(CC) -o $@ $^ $(CFLAGS)

ce675p2: ce675p2.o fe-solver.a material-data.a matrix.a
	$(CC) -o $@ $^ $(CFLAGS)

heat-explicit: heat-explicit.o fe-solver.a material-data.a matrix.a
	$(CC) -o $@ $^ $(CFLAGS)

heat-cyl: heat-cyl.o heat-gui.o fe-solver.a material-data.a matrix.a
	$(CC) -o $@ $^ $(CFLAGS)

heat-transfer: heat-transfer.o ht-main.o common.o fe-solver.a material-data.a matrix.a
	$(CC) -o $@ $^ $(CFLAGS)

ht-mt: diffusion.o heat-transfer.o main.o common.o fe-solver.a material-data.a matrix.a
	$(CC) -o $@ $^ $(CFLAGS)

diffusion: diffusion.o deformation.o mt-main.o common.o output.o fe-solver.a material-data.a matrix.a
	$(CC) -o $@ $^ $(CFLAGS)

diffusion-mod: diffusion.o deformation.o lin-genmaxwell.o mt-main.o common-mod.o fe-solver.a material-data.a matrix.a
	$(CC) -o $@ $^ $(CFLAGS)

diffusion-kelvin: diffusion.o deformation.o lin-genkelvin.o mt-main.o common-kelvin.o fe-solver.a material-data.a matrix.a
	$(CC) -o $@ $^ $(CFLAGS)

diffusion.o: diffusion.h common.h
deformation.o: deformation.h common.h
lin-maxwell.o: lin-maxwell.h
lin-genmaxwell.o: lin-genmaxwell.h
lin-genkelvin.o: lin-genkelvin.h
mt-main.o: common.h
common.o: common.h
common-mod.o: common-mod.h
common-kelvin.o: common-kelvin.h

clean:
#rm -rf spheroid 2dlaplace ce675p1 ce675p2 heat-explicit heat-cyl meshtest
	rm -rf diffusion-kelvin
	rm -rf *.o *.a
	$(MAKE) -C matrix clean
	$(MAKE) -C material-data clean
	$(MAKE) -C fe-solver clean
	rm -rf doc

freezing.o: material-data/freezing/freezing.c material-data/freezing/freezing.h
	$(CC) -c material-data/freezing/freezing.c $(CFLAGS)

heat-implicit.o: problems/heat-implicit.c
	$(CC) -c problems/heat-implicit.c $(CFLAGS)

heat-gui.o: heating/heat-gui.c heating/heat-gui.h
	$(CC) -c gui/heating/heat-gui.c $(CFLAGS)

fe-solver.a:
	$(MAKE) -C fe-solver fe-solver.a
	cp fe-solver/fe-solver.a .

matrix.a:
	$(MAKE) -C matrix matrix.a
	cp matrix/matrix.a .

material-data.a:
	$(MAKE) -C material-data material-data.a
	cp material-data/material-data.a .

