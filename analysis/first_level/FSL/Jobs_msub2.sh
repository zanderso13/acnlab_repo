#!/bin/bash
#source /home/cca7133/.bash_profile
#MSUB -A p30352
#MSUB -l mem=1gb
#MSUB -N MID_stats
#MSUB -l nodes=1:ppn=1
#MSUB -q normal


/projects/p30352/Programs/fsl/bin/feat /projects/p30352/MID_data/Analysis/template/design081818.fsf
