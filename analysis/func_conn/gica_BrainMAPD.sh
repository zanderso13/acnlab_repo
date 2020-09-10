#!/bin/bash

#SBATCH -A p30954
#SBATCH -p normal
#SBATCH -t 12:00:00
#SBATCH --mem=64G

matlab -nodisplay -nosplash -nodesktop -r "addpath(genpath('~/repo')); addpath(genpath('~/ACNlab/software/gift')); gica_cmd --data smooth_ICA_sub_list.txt --o '/projects/b1108/projects/BrainMAPD_func_conn/ICA'; quit"

