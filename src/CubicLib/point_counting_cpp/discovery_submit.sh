#!/bin/bash -l

# Name of the job
#SBATCH --job-name=point_counting

# Name of the partition / queue
#SBATCH --partition=standard

# Name of the cluster account
# How long should I job run for
#SBATCH --time=01:00:00

# Number of CPU cores, in this case 4 cores
#SBATCH --ntasks=4

# Number of compute nodes to use, in this case 2
#SBATCH --nodes=2

# Name of the output files to be created. If not specified the outputs will be joined
#SBATCH --output=output/%x.%j.out
#SBATCH --error=output/%x.%j.err

# The code you want to run in your job
g++ -std=c++11 -mpclmul -O4 -DCOEFFSFILE='"coefficient_headers/coeffs1.h"' -DN=11 count_bigq.cpp -o a1.out
./a1.out
