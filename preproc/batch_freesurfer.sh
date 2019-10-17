#!/usr/bin/bash


export basedir=/projects/p30954/struct_data/rerun_freesurfer/
cd $basedir

# load necessary modules
module load python/2.7.5
module load freesurfer/6.0

# freesurfer setup
source $FREESURFER_HOME/SetUpFreeSurfer.sh

for d in */; do
	if [ -d "$d/stats" ]; then
		continue
	else
		export SUBJECTS_DIR=$basedir$d
		sbatch /home/zaz3744/repo/acnlab_scripts/preproc/single_sub_freesurfer_job.sh
	fi
done
