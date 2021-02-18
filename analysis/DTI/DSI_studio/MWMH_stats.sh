#!/bin/bash

sub=$1         #subject txt file with list of subjects

projdir='/Volumes/ZachExternal/ACNlab/MWMH/dti/ses-1'   # path for folder with subject data
                # Output directory for all generated .trk files
fa_threshold=0            # QA threshold
thread_count=4               # Set low for personal computer, but should be able to max out at 24 threads on the cluster
seed_ori=primary 	       # Seed orientation 
seed_pos=voxel	       # Seed position

seed_count=2000000

conn_val=qa,count,ncount

#connectivity value- use qa only if track goes end to end over length of seed
#fullfib='20130228_F046Y_130228155338.src.gz.odf8.f5rec.bal.reg1i2.qsdr.1.25.2mm.jac.map.R70.fib.gz  # Full file name and extension for subject'

fibfile=$(ls ${projdir}/${fibdir}/${sub})
subid=$(echo ${sub}|cut -c1-11)
echo $fibfile

#step 1: generate .trk files for each subject

/Applications/dsi_studio.app/Contents/MacOS/dsi_studio --action=trk --source=${fibfile} --seed_count=${seed_count} --track_id=Uncinate_Fasciculus_L --export=stat,qa,dti_fa --output=${subid}_left_uncinate_fasciculus_output.tt.gz

/Applications/dsi_studio.app/Contents/MacOS/dsi_studio --action=trk --source=${fibfile} --seed_count=${seed_count} --track_id=Uncinate_Fasciculus_R --export=stat,qa,dti_fa --output=${subid}_right_uncinate_fasciculus_output.tt.gz

mv *output* /Volumes/ZachExternal/ACNlab/MWMH/dti/results/Uncinate
