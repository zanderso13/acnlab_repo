% This is going to be a neat little script that won't do much other than
% make it so I won't have to type all this in everytime I want to check how
% my first levels are looking. Also it will remind me that nothing is
% smoothed
% con_0001 - loss anticipation
% con_0002 - gain anticipation
% con_0003 - loss consumption
% con_0004 - gain consumption

% SPM directories
fldir = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD/func_conn/first_levels/first_level_output';
ses_str = 'ses-2';
run_str = 'run-2';
task_str = 'MID';

% FSL directories
% fldir = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD/MID_FSL_contrasts/anticipation';
% ses_str = '*';
% run_str = 'analysis/MID_Run1.feat';
% task_str = 'stats';

smoothing_kernel = 6;
data = fmri_data(filenames(fullfile(fldir,'*',ses_str,run_str,task_str,'con_0001.nii')));

data = preprocess(data, 'smooth', smoothing_kernel);

thresh_data = threshold(ttest(data),.001,'unc');

montage(thresh_data)
