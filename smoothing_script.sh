#!/bin/bash

#SBATCH -A p30954
#SBATCH -p short
#SBATCH -t 00:10:00
#SBATCH --mem=30G

matlab -nodisplay -nosplash -nodesktop -r "addpath(genpath('~/repo')); single_sub_smooth(21684, 2,1,0); quit"

