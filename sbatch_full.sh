#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --array=1-5
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH -J=felix_test
#SBATCH -o slurmouts/slurmout_%a.out
#SBATCH -e slurmouts/slurmout_%a.errarray

./modele_entier inst_${SLURM_ARRAY_TASK_ID} param1.txt

