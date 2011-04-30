#!/bin/bash 
echo 'Starting job'
date

mpirun -mode VN -cwd /gpfs/small/CFSI/home/CFSIwngc/IFEM.PointSet.Nonuniform/test_nu_V3/ ./IFEM

echo 'Job completed'
date
