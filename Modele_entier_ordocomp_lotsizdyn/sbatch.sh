#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --array=1-288
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH -J=inst_test
#SBATCH -o slurmouts_sansrelax/slurmout_%a_12.out
#SBATCH -e slurmouts_sansrelax/slurmout_%a_12.errarray
#SBATCH --exclude=trencavel.10g,balibaste.10g,nestorius.10g


./modele_entier inst_${SLURM_ARRAY_TASK_ID} param12.txt

