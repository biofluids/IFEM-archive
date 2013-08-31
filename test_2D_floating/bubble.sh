#!/bin/bash
#SBATCH --job-name=floating
#SBATCH -t 12:00:00
#SBATCH -n 512

date
srun ./IFEM   
date

