#!/bin/bash 
echo 'Starting job'
date

mpirun -mode VN -cwd /gpfs/sb/home/CFSI/CFSIwngc/code/IFEM.RT/RT ./IFEM

echo 'Job completed'
date
