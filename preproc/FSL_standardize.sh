#!/bin/bash

#SBATCH -A p30954
#SBATCH -p normal
#SBATCH -t 24:00:00
#SBATCH --mem=32G
#SBATCH -J fsl_standardization

# module load fsl

. ${FSLDIR}/etc/fslconf/fsl.sh

Usage() {
        echo ‘Usage: FSL_loop <subject list>’
        exit 0
}

[ ‘$1’ = ‘’ ] && Usage

# Subject list for loop
subs=`cat ${1}`          # make sure txt file is full fib file names

projdir='/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD/MID_FSL_contrasts/anticipation'

cd ${projdir}

for sub in ${subs}
do

file=${sub}'/cope1.nii.gz'

echo $file	
fnirt --in=${file} --ref=MNI152_T1_2mm --verbose --iout=${file}_standard

echo ‘done for subject: ‘ ${sub}
done

