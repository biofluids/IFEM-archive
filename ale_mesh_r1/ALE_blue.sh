#!/bin/bash 
echo 'Starting job'
date

mpirun -mode VN -cwd /gpfs/small/CFSI/home/CFSIwx12/ALE.pa.3/ale_mesh_r1/ ./IFEM

echo 'Job completed'
date
