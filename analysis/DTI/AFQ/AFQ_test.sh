#!/usr/bin/bash

#SBATCH -A p30954
#SBATCH -p short
#SBATCH -t 00:10:00
#SBATCH --mem=30G

module add matlab/r2016a

matlab -nosplash -nodesktop -r "addpath(genpath('~/repo')); AFQ_single_sub; quit"

