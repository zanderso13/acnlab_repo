#!/usr/bin/bash
#SBATCH --account p30954               # Allocation
#SBATCH --partition=normal                # Queue
#SBATCH --time 12:00:00             # Walltime/duration of the job
#SBATCH -N 1                    # Number of Nodes
#SBATCH -n 4
#SBATCH --mem=30G
#SBATCH --job-name="freesurfer 1 subject"       # Name of job

# include a sanity check and output the directories youre working with
# export SUBJECTS_DIR=/projects/p30954/struct_data/rerun_freesurfer/10034_MR1

echo $basedir
echo $SUBJECTS_DIR

# call freesurfer

recon-all -all -i $SUBJECTS_DIR/mprage.nii -s stats
