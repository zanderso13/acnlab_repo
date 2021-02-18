#!/bin/bash

# Subject list for loop
subs=`cat ${1}`          # make sure txt file is full fib file names

projdir=‘/Users/zaz3744/Documents/current_projects/dti_BrainMAPD/src’

cd ${projdir}

for sub in ${subs}
do
	
./MWMH_uncinate_stats.sh ${sub}

echo ‘done for subject: ‘ ${sub}
done

