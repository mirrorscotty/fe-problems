CC=gcc
CFLAGS=-lm -I. -Imatrix -Imaterial-data -Ife-solver -Wall -g3 -O2
SRC=$(wildcard other/*.c) \
    $(wildcard slab-drying/*.c) \
    $(wildcard gui/*.c)

all: diffusion

doc:
	doxygen DoxyFile
	make -C doc/latex
	cp doc/latex/refman.pdf doc/Reference.pdf

2dlaplace: other/2dlaplace.o fe-solver.a material-data.a matrix.a
	$(CC) -o $@ $^ $(CFLAGS)

ce675p1: other/ce675p1.o fe-solver.a material-data.a matrix.a
	$(CC) -o $@ $^ $(CFLAGS)

spheroid: other/spheroid.o fe-solver.a material-data.a matrix.a
	$(CC) -o $@ $^ $(CFLAGS)

ce675p2: other/ce675p2.o fe-solver.a material-data.a matrix.a
	$(CC) -o $@ $^ $(CFLAGS)

heat-explicit: other/heat-explicit.o fe-solver.a material-data.a matrix.a
	$(CC) -o $@ $^ $(CFLAGS)

heat-cyl: other/heat-cyl.o heat-gui.o fe-solver.a material-data.a matrix.a
	$(CC) -o $@ $^ $(CFLAGS)

heat-transfer: slab-drying/heat-transfer.o slab-drying/ht-main.o slab-drying/common.o fe-solver.a material-data.a matrix.a
	$(CC) -o $@ $^ $(CFLAGS)

ht-mt: slab-drying/diffusion.o slab-drying/heat-transfer.o slab-drying/main.o slab-drying/common.o fe-solver.a material-data.a matrix.a
	$(CC) -o $@ $^ $(CFLAGS)

diffusion: slab-drying/diffusion.o slab-drying/deformation.o slab-drying/mt-main.o slab-drying/common.o slab-drying/output.o fe-solver.a material-data.a matrix.a
	$(CC) -o $@ $^ $(CFLAGS)

diffusion-mod: slab-drying/diffusion.o slab-drying/deformation.o slab-drying/lin-genmaxwell.o slab-drying/mt-main.o slab-drying/common-mod.o fe-solver.a material-data.a matrix.a
	$(CC) -o $@ $^ $(CFLAGS)

diffusion-kelvin: slab-drying/diffusion.o slab-drying/deformation.o slab-drying/lin-genkelvin.o slab-drying/mt-main.o slab-drying/common-kelvin.o fe-solver.a material-data.a matrix.a
	$(CC) -o $@ $^ $(CFLAGS)

clean:
#rm -rf spheroid 2dlaplace ce675p1 ce675p2 heat-explicit heat-cyl meshtest
	rm -rf diffusion-kelvin diffusion
	rm -rf $(SRC:.c=.o)
	rm -rf $(SRC:.c=.d)
	rm -rf *.a
	$(MAKE) -C matrix clean
	$(MAKE) -C material-data clean
	$(MAKE) -C fe-solver clean
	rm -rf doc

fe-solver.a:
	$(MAKE) -C fe-solver fe-solver.a
	cp fe-solver/fe-solver.a .

matrix.a:
	$(MAKE) -C matrix matrix.a
	cp matrix/matrix.a .

material-data.a:
	$(MAKE) -C material-data material-data.a
	cp material-data/material-data.a .

%.o: %.c
	$(CC) -c $(CFLAGS) $*.c -o $*.o
	$(CC) -MM $(CFLAGS) $*.c > $*.d
	@mv -f $*.d $*.d.tmp
	@sed -e 's|.*:|$*.o:|' < $*.d.tmp > $*.d
	@sed -e 's/.*://' -e 's/\\$$//' < $*.d.tmp | fmt -1 | \
	  sed -e 's/^ *//' -e 's/$$/:/' >> $*.d
	@rm -f $*.d.tmp

-include $(OBJ:.o=.d)

