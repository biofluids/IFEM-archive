#!/bin/bash

HOSTFILE=/tmp/hosts.$SLURM_JOB_ID


srun hostname -s > $HOSTFILE

if [ -z "$SLURM_NPROCS" ] ; then
  if [ -z "$SLURM_NTASKS_PER_NODE" ] ; then
    SLURM_NTASKS_PER_NODE=1
  fi
  SLURM_NPROCS=$(( $SLURM_JOB_NUM_NODES * $SLURM_NTASKS_PER_NODE ))
fi

echo 'Starting job'
date

mpirun -machinefile $HOSTFILE -np $SLURM_NPROCS /fasttmp/jyang/IFEM.linels.pscom.bloodvessel/bloodvessel130222/IFEM

cat /tmp/hosts.$SLURM_JOB_ID
rm /tmp/hosts.$SLURM_JOB_ID

echo 'Job completed'
date

