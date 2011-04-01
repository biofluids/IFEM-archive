#!/bin/bash 
echo 'Starting job'
date

mpirun -mode VN -cwd /gpfs/small/CFSI/home/CFSIwx12/IFEM.imrkpm.3D.clean/Liz_3d_test_CVF/ ./IFEM

echo 'Job completed'
date
