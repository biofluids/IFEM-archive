This is the Version 2 of IFEM program.
July 2003.

new feature:
- removed the artificial viscosity effect sigma_(ij,j)^f of the fluid in the solid domain.
- quadrature integration used in solid, r_stang.f, is modified from the qudrature points per direction, instead, 
it is now using the qudrature definitions from the fluid solver.




usage:
compile&link: make
run: mpirun -np 4 xts-mpi<xts.in(input)

-----------------------------------------
changes made from original:

modified:
all c functions.
title=>lower case and '_'

files with included file "fortran.h":
removed _fcd, _fcdtocp, etc..

in ewdmem.c:
comment out the tracebk() and _memerror.
