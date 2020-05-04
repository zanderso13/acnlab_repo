#!/bin/bash

#SBATCH -A p30954
#SBATCH -p short
#SBATCH -t 00:10:00
#SBATCH --mem=32G
#SBATCH -J dti_create_trk

fibfile=$1         #subject txt file with list of subjects

trkdir=tracks                # Output directory for all generated .trk files
fa_threshold=0            # QA threshold
thread_count=4               # Set low for personal computer, but should be able to max out at 24 threads on the cluster
seed_ori=primary 	       # Seed orientation 
seed_pos=voxel	       # Seed position
fiber_count=500               # FROM KEVIN’S SCRIPT - N/A Number of streamlines to generate; not using
step_size=0.05              # FROM KEVIN’S SCRIPT - N/A Step size
turning_angle=0           # FROM KEVIN’S SCRIPT - N/A Turning angle cutoff threshold
smoothing=0                      # Degree of smoothing
min_length=0                         # Minimum fiber length in mm
max_length=100                        # Maximum fiber length in mm
turn_thresh=0
seed_count=2000000

conn_val=qa,count,ncount

#connectivity value- use qa only if track goes end to end over length of seed
#fullfib='20130228_F046Y_130228155338.src.gz.odf8.f5rec.bal.reg1i2.qsdr.1.25.2mm.jac.map.R70.fib.gz  # Full file name and extension for subject'

echo $fibfile

#step 1: generate .trk files for each subject

singularity exec /home/zaz3744/ACNlab/software/singularity_images/dsi-studio-docker.sif /dsistudio/dsi_studio_64/dsi_studio --action=trk --source=${fibfile} --method=0 --fa_threshold=${fa_threshold} --turning_angle=${turn_thresh} --initial_dir=0 --seed_plan=1 --interpolation=0 --step_size=${step_size} --smoothing=${smoothing} --min_length=${min_length} --max_length=${max_length}} --fiber_count=${fiber_count} --seed_count=${seed_count} --seed=FreeSurferSeg:Right-Accumbens-area --end=FreeSurferSeg:Right-Amygdala  --export=stat,qa,dti_fa --output=${projdir}/$1-tracks.trk.gz

