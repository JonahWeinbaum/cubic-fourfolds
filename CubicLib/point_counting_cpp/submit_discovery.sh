#!/bin/bash -l

# Name of the job
#SBATCH --job-name=count_bigq

# Number of cores
#SBATCH --ntasks-per-node=1

# Array jobs.  This example will create 25 jobs
#SBATCH --array=1-25

# Walltime (job duration)
#SBATCH --time=00:15:00

# Email notifications
#SBATCH --mail-type=BEGIN,END,FAIL

./count_bigq_$SLURM_ARRAY_TASK_ID

sleep 300
hostname -s