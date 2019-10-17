#!/usr/bin/bash

# In this script I'm going to be moving files around so that we can download
# NIFTI's from NUNDA, run this script, and have it be ready for bidskit after
# that. 

# USER SPECIFICATION
# Set your project directory here!

export raw_directory='/projects/p30954/func_data'

# This is where your sourcedata directory will live
# Where do you want to put your sourcedata folder?

export project_directory='/projects/p30954'

# We have lots of longitudinal data! Put in how many sessions we need to account for
# BrainMAPD has four sessions for instance

export sessions=2

# END USER SPECIFICATION

# start by checking for folders needed for bidskit conversion
# create them if they don't already exist

if [ -d "$project_directory/sourcedata" ]; then
	echo "sourcedata directory exists!"
else
	echo "Creating sourcedata folder"
	mkdir $project_directory/sourcedata
fi


for sess in `seq 1 $sessions`;do
        echo "Beginning Session $sess"

	if [ -d "$project_directory/sourcedata/$sess" ]; then
		echo "sourcedata sess $sess directory exists!"
	else
		echo "Creating sourcedata folder"
		mkdir $project_directory/sourcedata/$sess
	fi

	echo "Finding NIFTI files for session $sess...."
	echo "Transferring!"
	find $raw_directory/*MR$sess/ -iname *.nii -exec ln -s {} $project_directory/sourcedata/$sess/ \;
done
