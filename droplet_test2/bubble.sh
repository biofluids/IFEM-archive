#!/bin/bash 
echo 'Starting job'
date

mpirun -mode VN -cwd /gpfs/small/CFSI/home/CFSIwngc/IFEM.PointSet.Merge.V2/droplet_test2/ ./IFEM

echo 'Job completed'
date
