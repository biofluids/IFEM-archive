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
