#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --array=1-10
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH -J=felix_test
#SBATCH -o slurmouts/slurmout_%a.out
#SBATCH -e slurmouts/slurmout_%a.errarray

number=$(printf "%d" ${SLURM_ARRAY_TASK_ID})

./comparaison_lotsiz 20 4 2000 100 ${number}

