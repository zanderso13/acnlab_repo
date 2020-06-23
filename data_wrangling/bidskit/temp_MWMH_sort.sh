#!/bin/bash

# In this script I'm going to be moving files around so that we can download
# NIFTI's from NUNDA, run this script, and have it be ready for bidskit after
# that.

# USER SPECIFICATION
# Set your project directory here!

export raw_directory='/projects/b1108/data/to_format/MWMH_dicoms'

# This is where your sourcedata directory will live
# Where do you want to put your sourcedata folder?

export project_directory='/projects/b1108/data/to_format/MWMH_bidskit'

# We have lots of longitudinal data! Put in how many sessions we need to accoun$
# BrainMAPD has four sessions for instance

export sessions=2

# END USER SPECIFICATION

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

