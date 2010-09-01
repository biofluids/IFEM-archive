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

/usr/local/mpich/latest/ch_p4/bin/mpirun -machinefile $HOSTFILE -np $SLURM_NPROCS /borg/xwang/IFEM.pa.new.2/test_panew/IFEM

cat /tmp/hosts.$SLURM_JOB_ID
rm /tmp/hosts.$SLURM_JOB_ID

echo 'Job completed'
date

