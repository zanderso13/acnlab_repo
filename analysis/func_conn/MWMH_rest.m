% Need to write a script to apply the 264 ROI set to MWMH rest data. This
% will also act as a template for any other study I want to apply the 264
% masks that are now in my repo. 

% 1. Visualize masks
% 2. Load in residual files. 
% 3. Write those residuals out into a single nii file per sub to make
% loading easier in the future
% 4. Apply all 264 masks to time course data for each subject
% 5. Will result in another set of time course data for each ROI, run
% corrcoef on that to get corr matrix
% 6. Fishers r to z
% 7. Save resulting files .mat