#!/usr/bin/bash

#SBATCH -A p30954
#SBATCH -p normal
#SBATCH -t 24:00:00
#SBATCH --mem=64G
#SBATCH -J fmriprep_single_sub

module purge
module load singularity/latest
echo "modules loaded" 
cd /projects/p30954

echo "beginning preprocessing"

singularity run --cleanenv -B /projects/p30954:/projects/p30954/ singularity_images/fmriprep-20.2.0.simg --skip-bids-validation /projects/p30954/Schizconnect_final/bids/COBRE/ /projects/p30954/Schizconnect_final/preproc/ participant --participant-label ${1} -t rest --fs-license-file /home/zaz3744/projects/freesurfer_license/license.txt --fs-no-reconall --clean-workdir -w /home/zaz3744/projects/Schizconnect_final/preproc/work
# change above to --clean-workdir so that the work directory is cleared and you save space
