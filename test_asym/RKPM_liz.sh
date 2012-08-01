#!/bin/bash 
echo 'Starting job'
date

mpirun -mode VN -cwd /gpfs/small/CFSI/home/CFSIwx12/IFEM.linels.com/vf_liz_den_soft/ ./IFEM

echo 'Job completed'
date
