#!/bin/bash 
echo 'Starting job'
date

mpirun -mode VN -cwd /gpfs/small/CFSI/home/CFSIwx12/IFEM.linels.com.V5.3part/vf_half_male_3p/ ./IFEM

echo 'Job completed'
date
