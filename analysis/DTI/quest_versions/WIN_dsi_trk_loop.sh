#!/bin/bash

# Subject list for loop
subs=`cat ${1}`          # make sure txt file is full fib file names

for sub in ${subs}
do

# can only run these one at a time because the output files have the same names. Kinda a bummer but submitting 1500 jobs all at once would also probably be frowned upon
	
sbatch BrainMAPD_dsi_trk_L-NAcc_Amyg ${sub}

#sbatch BrainMAPD_dsi_trk_L-Amyg_mOFC ${sub}

#sbatch BrainMAPD_dsi_trk_L-NAcc_mOFC ${sub}

#sbatch BrainMAPD_dsi_trk_R-NAcc_Amyg ${sub}

#sbatch BrainMAPD_dsi_trk_R-Amyg_mOFC ${sub}

#sbatch BrainMAPD_dsi_trk_R-NAcc_mOFC ${sub}

echo ‘done for subject: ‘ ${sub}
done

