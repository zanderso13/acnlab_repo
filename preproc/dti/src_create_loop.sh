#!/bin/bash


# DSI studio is very touchy about paths. You need to have everything running from a directory that's on the absolute path
# of the target file you're trying to convert.


# This input is a txt file with all the dwi file names and their absolute paths. Loop through this

subs=`cat ${1}`
echo "Just before loop"
for sub in ${subs}
do
ibvec=${sub/nii.gz/bvecs}
ibval=${sub/nii.gz/bvals}
echo $ibvec
echo $ibval
sbatch /projects/b1108/single_sub_create_src_quest.sh ${sub} ${ibvec} ${ibval}
echo "after call to single sub script"
echo ‘done for subject: ‘ ${sub}
done
