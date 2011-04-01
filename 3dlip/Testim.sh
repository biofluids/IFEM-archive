#!/bin/bash 
echo 'Starting job'
date

mpirun -mode VN -cwd /gpfs/small/CFSI/home/CFSIwx12/IFEM.imrkpm.3D.clean/3dlip/ ./IFEM

echo 'Job completed'
date
