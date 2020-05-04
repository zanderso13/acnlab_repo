#!/bin/bash                                                                                                                             

# Subject list                                                                                                                          
subs=`cat /Users/zaz3744/Documents/current_projects/dti_BrainMAPD/src/src_files.txt`                                                                                                    

# Directory information                                                                                                                 
projdir="/Users/zaz3744/Documents/current_projects/dti_BrainMAPD/src"            # Primary project directory                                                                     

# Loop through tractography           


                                                                                                  
for sub in $subs
do
	echo "Running Subject " ${sub}
	/Users/zaz3744/Documents/repo/acnlab_repo/data_wrangling/recon_WIN ${projdir}/${sub}
	echo "Done Subject " ${sub}		
done
                                                                                                                      

                                                                                                                      
                                                                                             

                                                                                                               
                                                                   

                                                                                                           
                                                                                                   
                                                                                                      


