#!/bin/bash 
echo 'Starting job'
date

mpirun -mode VN -cwd /gpfs/small/CFSI/home/CFSIwngc/IFEM.PointSet.Merge.V2/droplet_test1/ ./IFEM

echo 'Job completed'
date
