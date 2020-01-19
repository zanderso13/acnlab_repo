#!/usr/bin/bash

#SBATCH -A p30954
#SBATCH -p normal
#SBATCH -t 24:00:00
#SBATCH --mem=64G
#SBATCH -J fmriprep_test

module purge
module load singularity/latest
echo "modules loaded" 
cd /projects/b1108
pwd
echo "Begin Preprocessing"

singularity run --cleanenv /projects/b1108/software/singularity_images/fmriprep-1.5.4.simg /projects/b1108/data/BrainMAPD/ /projects/b1108/projects/BrainMAPD_func_conn/ participant --participant-label 10001 --fs-license-file ~/freesurfer_license/license.txt --fs-no-reconall -w /projects/b1108/projects/BrainMAPD_func_conn/work
