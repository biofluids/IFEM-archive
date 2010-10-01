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


srun --mpi=none --nodes $SLURM_NPROCS /gpfs/small/CFSI/home/CFSIwx12/IFEM.fluid.pi.pa.3/chantest_water_givens_pa/IFEM

rm /tmp/host.$SLURM_JOB_ID
~                             
