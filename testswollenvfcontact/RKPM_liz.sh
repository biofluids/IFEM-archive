#!/bin/bash 
echo 'Starting job'
date

mpirun -mode VN -cwd /gpfs/small/CFSI/home/CFSIynjb/IFEM.linels.com/vf20120731a/ ./IFEM

echo 'Job completed'
date
