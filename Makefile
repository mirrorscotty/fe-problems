CC=gcc
CFLAGS=-lm -I. -Islab-drying -Imatrix -Imaterial-data -Ife-solver -Wall -g3 -O2 -lblas -llapack
SRC=$(wildcard other/*.c) \
    $(wildcard slab-drying/*.c) \
    $(wildcard slab-drying/mass/*.c) \
    $(wildcard slab-drying/heat/*.c) \
    $(wildcard slab-drying/solid/*.c) \
    $(wildcard gui/*.c)

all: diffusion-kelvin

doc:
	doxygen DoxyFile
	make -C doc/latex
	cp doc/latex/refman.pdf doc/Reference.pdf

diffusion-kelvin: slab-drying/mass/diffusion.o slab-drying/solid/deformation.o slab-drying/solid/lin-genkelvin.o slab-drying/solid/kinematics.o slab-drying/mt-main.o slab-drying/output.o slab-drying/common-kelvin.o fe-solver/fe-solver.a material-data/material-data.a matrix/matrix.a
	$(CC) -o $@ $^ $(CFLAGS)

pc-test2: slab-drying/solid/deformation.o pc-test2.o fe-solver/fe-solver.a material-data/material-data.a matrix/matrix.a
	$(CC) -o $@ $^ $(CFLAGS)

clean:
#rm -rf spheroid 2dlaplace ce675p1 ce675p2 heat-explicit heat-cyl meshtest
	rm -rf diffusion-kelvin
	rm -rf $(SRC:.c=.o)
	rm -rf $(SRC:.c=.d)
	rm -rf *.a
	$(MAKE) -C matrix clean
	$(MAKE) -C material-data clean
	$(MAKE) -C fe-solver clean
	rm -rf doc

mostlyclean:
	rm -rf diffusion-kelvin
	rm -rf $(SRC:.c=.o)
	rm -rf $(SRC:.c=.d)

fe-solver/fe-solver.a: force_build
	$(MAKE) -C fe-solver fe-solver.a

matrix/matrix.a: force_build
	$(MAKE) -C matrix matrix.a

material-data/material-data.a: force_build
	$(MAKE) -C material-data material-data.a

force_build:
	true

%.o: %.c
	$(CC) -c $(CFLAGS) $*.c -o $*.o
	$(CC) -MM $(CFLAGS) $*.c > $*.d
	@mv -f $*.d $*.d.tmp
	@sed -e 's|.*:|$*.o:|' < $*.d.tmp > $*.d
	@sed -e 's/.*://' -e 's/\\$$//' < $*.d.tmp | fmt -1 | \
	  sed -e 's/^ *//' -e 's/$$/:/' >> $*.d
	@rm -f $*.d.tmp

-include $(OBJ:.o=.d)

