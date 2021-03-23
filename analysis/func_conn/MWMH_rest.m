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

% First, a quick thing for Robin. Need to check how many volumes are
% dropped per sub and apply various motion cutoffs.

% 1.
maskdir = '/home/zaz3744/repo/acnlab_repo/masks';

load(fullfile(maskdir,'Power_atlas_obj.mat'))

%% 2. 

datadir = '/projects/b1108/projects/MWMH_project/first_levels';
testdir = 'sub-MWMH102/*/*/rest/'

fnames = filenames(fullfile(datadir,testdir,'Res*.nii'));

dat = fmri_data(fnames);




