#!/bin/bash 
echo 'Starting job'
date

mpirun -mode VN -cwd /gpfs/small/CFSI/home/CFSIwngc/IFEM.SLIP.Droplet/droplet1 ./IFEM

echo 'Job completed'
date
