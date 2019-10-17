#!/usr/bin/bash
# Above is the shebang, you have to start every script with this
# Below is all the information about your slurm job or basically what kind of computer do you want to use. This is an example of some of the parameters you might use

#SBATCH --account p30954               # Allocation
#SBATCH --partition=short                # Queue
#SBATCH --time 04:00:00             # Walltime/duration of the job
#SBATCH -N 1                    # Number of Nodes
#SBATCH -n 4
#SBATCH --mem=30G
#SBATCH --job-name="create stats tables"       # Name of job

# This is how you make a comment! The terminal will ignore this and you can take notes about whatever is happening next
# Slurm works with modules, so this is how you access software on the cluster

module purge 
# start fresh by getting rid of all your modules
module load python/2.7.5 
# you need this version of python for freeview
module load freesurfer/6.0 
# here's the most recent version of freesurfer

# set environment variables for freesurfer. You're familiar with this one!
source $FREESURFER_HOME/SetUpFreeSurfer.sh

# this is how you create a variable! This is setting up a base directory to refer to freesurfer filenames down the line
export basedir=/projects/p30954/struct_data/session1/

# Beginning a for loop!
for d in $basedir/*/; do
# Now I want to set a subject specific directory
export SUBJECTS_DIR=$d

# if statements check for a condition to see if it exists. If it does, it will do what you say if not it will move to a new statement
if [ -f $SUBJECTS_DIR/stats/mri/*2009* ]; then
	# echo is a useful bash command. It just prints whatever you place after it. Here I echo the subject directory as a sanity check to make sure bash and I are on the same page as to who we are running
	echo $d
	# this is a freesurfer specific command. Google what it does!
	aparcstats2table --subjects stats --parc aparc.a2009s --hemi lh --meas volume --tablefile lh_vol_a2009s.txt
	aparcstats2table --subjects stats --parc aparc.a2009s --hemi rh --meas volume --tablefile rh_vol_a2009s.txt
	# mv is also a command you'll use a lot. It can rename a file or move it around. 
	mv lh_vol_a2009s.txt rh_vol_a2009s.txt $d
fi # You have to tell the if statement that you're finished. In Bash this is how you do it
done # you also have to tell the for loop that you're finished and this is how you do it.
