#!/bin/bash 
echo 'Starting job'
date

mpirun -mode VN -cwd /gpfs/small/CFSI/home/CFSIwngc/IFEM.PointSet.CA.V2/ca_long_sm/ ./IFEM

echo 'Job completed'
date
