#!/bin/bash 
echo 'Starting job'
date

mpirun -mode VN -cwd /gpfs/small/CFSI/home/CFSIwngc/IFEM.PointSet.Dam.V3/bubble_dam_long/ ./IFEM

echo 'Job completed'
date
