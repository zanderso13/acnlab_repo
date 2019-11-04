#!/bin/bash

# In this script I'm going to be moving files around so that we can download
# NIFTI's from NUNDA, run this script, and have it be ready for bidskit after
# that. 

# USER SPECIFICATION
# Set your project directory here!

export raw_directory='/Users/anncarroll/BIDS/data/MWMH_dicoms'

# This is where your sourcedata directory will live
# Where do you want to put your sourcedata folder?

export project_directory='/Users/anncarroll/BIDS/data/MWMH_BIDS'

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

# Next, I'm running into some folder name issues. So I'm going to sort all of the initial data folders by session first
# This will make it easier for me to reference/find them later


for sort in `seq 1 $sessions`;do
	if [ -d "$raw_directory/$sort" ]; then
		echo "Raw data has already been sorted!"
	else
		echo "Sorting raw directory by session"
		# This will actually move subject folders into session folders
		# So we kind of begin the restructuring here
		mkdir $raw_directory/$sort
		mv $raw_directory/*CV$sort $raw_directory/$sort
	fi
	echo "renaming subject folders"
	cd $raw_directory/$sort
	rename -s CV$sort '' *
done

# Need to rename my folders so they don't have session numbers in them.
# I swear this makes sense...

for sess in `seq 1 $sessions`;do
	echo "Starting session $sess"
	cd $raw_directory/$sess
	for sub in */;do
        	echo "Beginning Subject $sub"

		if [ -d "$project_directory/sourcedata/$sub/$sess" ]; then
			echo "sourcedata sub $sub sess $sess directory exists!"
		else
			echo "Creating subject folder"
			mkdir $project_directory/sourcedata/$sub
			echo "Creating session folder $sess"
			mkdir $project_directory/sourcedata/$sub/$sess
		fi
	
				
		echo "Finding NIFTI files for session $sess...."
		echo "Transferring!"
		find $raw_directory/$sess/$sub/ -iname *.dcm -exec ln -s {} $project_directory/sourcedata/$sub/$sess/ \;
	done
done

# do a final renaming of subject folders to remove the _MR
# This syntax is going to be janky but will work
# OLD CODE
# cd $project_directory/sourcedata
# for d in *CV;do
# mv "$d" "${d::(-3)}"
# done



