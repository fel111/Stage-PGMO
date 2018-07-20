#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --array=1-288
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH -J=inst_test
#SBATCH -o slurmouts_nomaxboucle/slurmout_%a_1.out
#SBATCH -e slurmouts_nomaxboucle/slurmout_%a_1.errarray
#SBATCH --exclude=trencavel.10g,balibaste.10g,nestorius.10g


./modele_entier inst_${SLURM_ARRAY_TASK_ID} param1.txt

