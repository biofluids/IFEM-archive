#!/bin/bash 
echo 'Starting job'
date

mpirun -mode VN -cwd /gpfs/small/CFSI/home/CFSIwx12/IFEM.pa.4.im_tao.ccni/thom_air_non/ ./IFEM

echo 'Job completed'
date
