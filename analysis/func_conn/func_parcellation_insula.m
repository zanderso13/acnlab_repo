% This script is meant to perform a functional parcellation of the insula
% using functional conncetivity and k-means clustering. We're replicating
% methods found at this link:
% https://academic.oup.com/cercor/article/21/7/1498/331510 
% Steps involved in this script are as follows: 1. Identify functional
% files that we want to use to localize functionally distinct parcels of
% insula. 2. Load those images in 3. Load in an appropriate atlas to use as
% an insula mask 4. apply the mask and perform whole brain voxelwise func
% conn 5. k-means clustering function to isolate functionally distinct
% portions of insula.

% Hypotheses: Stay tuned but they're going to have a developmental flair.
% Publication above found 3 clusters optimally describe the functional
% parcellation (ant, post, middle). The goal here is to track these
% functional parcellations as they develop and also to observe differences
% in those with depression. A separate publication:
% https://www.sciencedirect.com/science/article/pii/S1878929313000455 noted
% decreased anterior insula volume mapped onto decreased self reported
% impulsivity. It's possible similar decreases in volume might predict
% improvements in other clinical measures. Since we are functionally
% defining our ROI's, we would be predicting reduced cluster size in ant
% insula would correspond with reduced depressive symtpom experience

% find all the func files
fnames = fullfile(filenames('/projects/p30954/func_data/*/SCANS/15/NIFTI/*.nii'));
% Loop through all the func files
for sub=1:length(fnames)
    % load func files into data obj
    disp('Working on subject:',fnames{sub})
    fmri_data(fnames{sub})
    