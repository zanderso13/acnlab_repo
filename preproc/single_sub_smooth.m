
function smooth_single_sub(PID,ses,run,overwrite)
%% var set up
if nargin==0 % defaults just for testing
    PID = 10001;  
    overwrite = 1;
    ses = 1;
    run = 1;
end

% directories
% first is where your stats files will be output to
directories{1} = '/projects/b1108/projects/BrainMAPD_func_conn/first_levels/first_level_output';
% next is where the preprocessed data is
directories{2} = '/projects/b1108/projects/BrainMAPD_func_conn/fmriprep';
% where the raw data lives (raw meaning before preprocessing)
directories{3} = '/projects/b1108/data/BrainMAPD';
% the timing files for modelling (onsets, durations, names)
directories{4} = '/projects/b1108/projects/BrainMAPD_func_conn/timing_files';
% where framewise displacement files will be saved
directories{5} = '/projects/b1108/projects/BrainMAPD_func_conn/first_levels/additional_files';

fl_dir = directories{1};
preproc_dir = directories{2};
raw_dir = directories{3};
timing_dir = directories{4};
save_dir = directories{5};

if nargin==1
    overwrite = 1;
end  

PID = strcat('sub-',num2str(PID));
ndummies=0;
% FL directory for saving 1st level results: beta images, SPM.mat, etc.
% in{1} = {fullfile(fl_dir, PID, strcat('ses-',num2str(ses)), strcat('run-', num2str(run)), 'MID')};

rundir = fullfile(preproc_dir, PID, strcat('ses-', num2str(ses)), 'func');
in{1} = cellstr(spm_select('ExtFPList', rundir, strcat('.*task-MID_run-0',num2str(run),'_space-MNI152NLin2009cAsym_desc-preproc_bold.nii'), ndummies+1:9999));

if isempty(in{1}{1})
    warning('No preprocd functional found')
    return
end

jobfile = {'/home/zaz3744/repo/acnlab_repo/preproc/smooth_template.m'};
jobs = 'smooth_template.m';


spm('defaults', 'FMRI');
spm_jobman('run', jobs, in{:});


end
