#!/bin/bash --login

#SBATCH --nodes=1
#SBATCH --ntasks=28
#SBATCH --cpus-per-task=1
#SBATCH --time=1-00:00:00
#SBATCH --account=y95
#SBATCH --export=NONE

module load singularity/3.3.0

export OMP_NUM_THREADS=1

srun --export=all -N 1 -n 28 -d 1 singularity exec ./panann-plus.sif run_pasa.sh
