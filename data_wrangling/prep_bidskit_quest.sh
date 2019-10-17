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


# END USER SPECIFICATION

# start by checking for folders needed for bidskit conversion
# create them if they don't already exist

if [ -d "$project_directory/sourcedata" ]; then
	echo "sourcedata directory exists!"
else
	echo "Creating sourcedata folder"
	mkdir $project_directory/sourcedata
fi


if [ -d "$project_directory/sourcedata/1" ]; then
	echo "Session 1 project directory exists!"
else
	echo "Creating a session 1 project directory"
	mkdir $project_directory/sourcedata/1
fi

if [ -d "$project_directory/sourcedata/2" ]; then
	echo "Session 2 project directory exists!"
else
	echo "Creating a session 2 project directory"
	mkdir $project_directory/sourcedata/2
fi

echo "Finding NIFTI files for session 1...."
echo "Transferring!"
find $raw_directory/*MR1/ -iname *.nii -exec ln -s {} $project_directory/sourcedata/1/ \;

echo "Finding NIFTI files for session 2...."
echo "Transferring!"
find $raw_directory/*MR2/ -iname *.nii -exec ln -s {} $project_directory/sourcedata/2/ \;
