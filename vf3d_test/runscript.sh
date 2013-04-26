#!/bin/bash
#SBATCH --job-name=3dVFfulltest
#SBATCH -t 23:30:00
#SBATCH -n 1024

date
srun ./IFEM
date
