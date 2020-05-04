#!/bin/bash

# Reconstruction Parameters                                                                                           
param0="1.25"
record_odf="1"
reg_method="1"
scheme_balance="1"
check_btable="1"
voxel_res="2"     
thread="4"


/Applications/dsi_studio.app/Contents/MacOS/dsi_studio --action=rec --thread=${thread} --source=${1} --method=7 --param0=${param0} --odf_order=8 --param1=${voxel_res} --output_jac=1 --output_map=1 --record_odf=1 --reg_method=${reg_method}


