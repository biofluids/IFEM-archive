#!/bin/bash
srun hostname -s > /tmp/host.$SLURM_JOB_ID
if [ "x$SLURM_NPROCS" = "x" ]
then
  if [ "x$SLURM_NTASKS_PER_NODE" = "x" ]
  then
    SLURM_NTASKS_PER_NODE=1
    fi
    SLURM_NPROCS=`expr $SLURM_JOB_NUM_NODES \* $SLURM_NTASKS_PER_NODE`
fi

MPIRUN=/opt/openmpi/1.2.6-V00/gcc/bin/mpirun

$MPIRUN -hostfile /tmp/host.$SLURM_JOB_ID -np $SLURM_NPROCS /gpfs/small/CFSI/home/CFSIwx12/IFEM.DIAG.pa.vocal.opterons/test_air_uni/IFEM

rm /tmp/host.$SLURM_JOB_ID
~                             
