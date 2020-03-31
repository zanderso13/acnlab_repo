#!/bin/bash

#SBATCH -A p30954
#SBATCH -p short
#SBATCH -t 4:00:00
#SBATCH --mem=8G
#SBATCH -J fsl_first_level_test

module purge fsl
export FSLDIR='/home/zaz3744/projects/software/fsl'

. /home/zaz3744/projects/software/fsl/etc/fslconf/fsl.sh

subdir=$1

cd /home/zaz3744/projects/${subdir}/analysis/MID_Run1.feat

/home/zaz3744/projects/software/fsl/bin/feat Anticipation.fsf
