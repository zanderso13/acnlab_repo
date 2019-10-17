#!/usr/bin/bash
#SBATCH --account p30954               # Allocation
#SBATCH --partition=short                # Queue
#SBATCH --time 00:00:30             # Walltime/duration of the job
#SBATCH -N 1                    # Number of Nodes
#SBATCH -n 4
#SBATCH --mem=30G
#SBATCH --job-name="wrapper test"       # Name of job

# load necessary modules 
# module load python/2.7.5
# module load freesurfer/6.0

# freesurfer setup
# source $FREESURFER_HOME/SetUpFreeSurfer.sh

# set your new subject directory where everyones data lives

# change to that directory I think this is a good place to eventually save the aparc table here

echo $basedir
echo $SUBJECTS_DIR

# recon-all -all -i $SUBJECTS_DIR/mprage.nii -s stats
