#!/bin/bash 
echo 'Starting job'
date

mpirun -mode VN -cwd /gpfs/small/CFSI/home/CFSIwngc/IFEM.SLIP.Droplet/droplet2_v2 ./IFEM

echo 'Job completed'
date
