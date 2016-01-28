FE Problems
===========
* Author: Alex Griessman (<alex.griessman@gmail.com>)
* Repository: https://github.com/mirrorscotty/fe-problems/

This repository contains a set of PDEs to be solved using the fe-solver library.

Contents
--------
* `gui` - Graphical interface for the a simple heat and mass transfer problem written in C++
* `other` - Assorted (old and untested) problems
* `pasta-cwp` - Non-functional start of a model based on Zhu (2012)
* `slab-drying` - Pasta drying model for slab geometry

Building
--------
Type `make` to compile the pasta drying model.

Usage
-----
This program requires a table of creep compliance data at the correct temperature to be saved in a file called `output.csv` in the same directory that the program is run from. This file can be generated from the `creep-table` program in the regression library located [here](https://github.com/mirrorscotty/regression).

Dependencies
------------
Compilation requires GCC and GNU Make. Building documentation requires Doxygen and LaTeX to generate PDF output. This code also requires the matrix library found [here](https://github.com/mirrorscotty/matrix), the material data library from [here](https://github.com/mirrorscotty/material-data), and the finite element solver library from [here](https://github.com/mirrorscotty/fe-solver). Additionally, the program won't run without a set of creep data generated using the regression code [here](https://github.com/mirrorscotty/regression).
