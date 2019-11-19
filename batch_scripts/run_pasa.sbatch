#!/bin/bash --login

#SBATCH --array=0,1,2,3,4
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --partition=workq
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --account=fl3
#SBATCH --export=NONE

module load singularity/3.3.0
module load samtools

export OMP_NUM_THREADS=1

srun --export=all -N 1 -n 24 -d 1 singularity exec ./panann-plus.sif batch_scripts/run_pasa.sh
