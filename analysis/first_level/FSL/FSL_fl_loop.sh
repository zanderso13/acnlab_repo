basedir='/home/zaz3744/ACNlab/projects/MID_FSL_contrasts/MID_data'

subs=`cat ${1}`

for sub in ${subs}
do

/home/zaz3744/repo/acnlab_repo/analysis/first_level/FSL/FSL_single_sub.sh ${sub}

done
