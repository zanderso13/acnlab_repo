#!/bin/bash
#source /home/cca7133/.bash_profile
#SBATCH -A p30352
#SBATCH --mem=1000
#SBATCH --job-name=MID_stats
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -p normal
#SBATCH -t 00:10:00

#change name of .fsf file#

/projects/p30352/Programs/fsl/bin/feat /projects/p30352/MID_data/Analysis/template/design_slurm_test.fsf
