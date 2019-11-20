% Alright time to write my first func conn script for the lab! Going to use
% the PPI approach Robin likes. This is going to mean the following series
% of steps

% 1. Load in functional data. This data has been preprocessed on NUNDA so
% I'm not entirely sure how crazy I am about that but for the sake of
% getting initial results and prepping a poster I think it's ok.
% 2. Plot FD metrics. Gotta make sure exclusions have already been taken
% into account
% 3. Plot averaged time course across all subs. Should roughly look like a
% task
% 4. Load in task/timing data. Visualize that. Should also basically look
% like a task
% 5. Parse func data into separate tasks
% 6. Load in mask, and extract time series data for that mask across all
% these different windows
% 7. Correlate masked out data in each time window

% set up directories
data_dir = '/home/zaz3744/projects/BrainMAPD';
repo_dir = '~/repo/';
spm_dir = '~/spm12';

%% Only run once
addpath(genpath(repo_dir))
addpath(genpath(spm_dir))
%% Alright let's get into the meat of it now
fnames = filenames(fullfile(data_dir, '*/func/MID1.nii'));
data = fmri_data(fnames)