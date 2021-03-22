#!/usr/bin/bash

#SBATCH -A p30954
#SBATCH -p normal
#SBATCH -t 12:00:00
#SBATCH --mem=64G
#SBATCH -J qsiprep_single_sub

module purge
module load singularity/latest
echo "modules loaded"

echo "beginning preprocessing"

singularity run --cleanenv -B /projects/b1108:/projects/b1108/ /projects/b1108/software/singularity_images/qsiprep-0.11.0.sif /projects/b1108/data/MWMH /projects/b1108/projects/MWMH_project participant --participant-label ${1} --fs-license-file /projects/b1108/software/freesurfer_license/license.txt --skip-bids-validation --output-resolution 222 -w /projects/b1108/projects/MWMH_project/work
