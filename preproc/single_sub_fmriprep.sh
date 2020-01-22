#!/usr/bin/bash

#SBATCH -A p30954
#SBATCH -p normal
#SBATCH -t 48:00:00
#SBATCH --mem=64G
#SBATCH -J fmriprep_test

module purge
module load singularity/latest
echo "modules loaded" 
cd /projects/b1108
pwd
echo "beginning preprocessing"

singularity run --cleanenv /projects/b1108/software/singularity_images/fmriprep-1.5.4.simg /projects/b1108/data/BrainMAPD/ /projects/b1108/projects/BrainMAPD_func_conn/ participant --participant-label 10002 --fs-license-file ~/freesurfer_license/license.txt -w /projects/b1108/projects/BrainMAPD_func_conn/work
