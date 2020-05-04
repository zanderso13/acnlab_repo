#!/bin/bash

# DSI studio is very touchy about paths. You need to have everything running from a directory that's on the absolute path
# of the target file you're trying to convert.

# Input should be a txt file with filenames and absolute paths to those files

subs=`cat ${1}`

for sub in ${subs}
do

sbatch /projects/b1108/single_sub_recon_fib_quest.sh ${sub}

done

