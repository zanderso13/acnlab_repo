#!/bin/bash

sub_ids=`cat ${1}`

starting_id=$2

for sub in ${sub_ids}
do

sed -i.bak "s/${starting_id}/${sub}/g" /projects/b1108/projects/MID_FSL_contrasts/zach_scripts/Anticipation.fsf

sbatch /projects/b1108/projects/MID_FSL_contrasts/zach_scripts/FSL_single_sub.sh ${sub}

echo $sub
echo $starting_id

starting_id=${sub}

done
