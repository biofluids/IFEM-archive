#!/bin/bash 
echo 'Starting job'
date

mpirun -mode VN -cwd /gpfs/small/CFSI/home/CFSIweda/ALE.pa.3.dyn/movingnaca00095aoatest4/ ./IFEM

echo 'Job completed'
date
