#!/bin/bash 
echo 'Starting job'
date

mpirun -mode VN -cwd /gpfs/small/CFSI/home/CFSIwngc/IFEM.PointSet.CA.V2/test_ca_den/ ./IFEM

echo 'Job completed'
date
