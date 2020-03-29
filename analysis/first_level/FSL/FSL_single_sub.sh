#!/bin/bash

#SBATCH -A p30954
#SBATCH -p short
#SBATCH -t 1:00:00
#SBATCH --mem=8G
#SBATCH -J fsl_first_level_single_sub

module purge fsl
export FSLDIR='/home/zaz3744/projects/software/fsl'

. /home/zaz3744/projects/software/fsl/etc/fslconf/fsl.sh

subdir=${1}

cd /home/zaz3744/ACNlab/projects/MID_FSL_contrasts/MID_data/${subdir}
pwd

/home/zaz3744/projects/software/fsl/bin/feat /home/zaz3744/repo/acnlab_repo/analysis/first_level/FSL/GainvNoGain.fsf
