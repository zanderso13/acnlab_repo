% This is going to be a neat little script that won't do much other than
% make it so I won't have to type all this in everytime I want to check how
% my first levels are looking. Also it will remind me that nothing is
% smoothed

fldir = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD/func_conn/first_levels';
ses_str = 'ses-2';
run_str = 'run-1';
task_str = 'MID';

smoothing_kernel = 4;
data = fmri_data(filenames(fullfile(fldir,'*',ses_str,run_str,task_str,'con_0001.nii')));

data = preprocess(data, 'smooth', smoothing_kernel);

thresh_data = threshold(ttest(data),.05,'unc');

montage(thresh_data)
