#!/bin/bash

#SBATCH -A p30954
#SBATCH -p short
#SBATCH -t 00:10:00
#SBATCH --mem=30G

matlab -nodisplay -nosplash -nodesktop -r "addpath(genpath('~/repo')); single_sub_smooth(10172, 1,2,0); quit"

