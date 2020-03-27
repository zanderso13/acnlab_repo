% run first level models for sponpain task for a subject
% the URSI and whether to overwrite existing first levels.  Overwrite = 0
% (default) or 1. before running this file, must run
% read_timings_make_onsets.m

function run_subject_firstlevel_MID(PID, ses, run, directories, overwrite)


%% var set up
if nargin==0 % defaults just for testing
    PID = 20341;  
    overwrite = 1;
    ses = 2;
    run = 2;
    % directories
    fl_dir = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD/func_conn/first_levels';
    preproc_dir = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD/func_conn/fmriprep';
    raw_dir = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD/func_conn/raw';
    timing_dir = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD/func_conn/final_timing_files';
    save_dir = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD/func_conn/first_levels';
end

fl_dir = directories{1};
preproc_dir = directories{2};
raw_dir = directories{3};
timing_dir = directories{4};
save_dir = directories{5};

if nargin==1
    overwrite = 1;
end  

PID = strcat('sub-',num2str(PID));

fprintf(['Preparing 1st level model for MID task for ' PID ' / ' ses], ['Overwrite = ' num2str(overwrite)]);


ndummies = 2;
TR = 2;

%% Model for MID task

% FL directory for saving 1st level results: beta images, SPM.mat, etc.
in{1} = {fullfile(fl_dir, PID, 'ses-2', strcat('run-', num2str(run)), 'MID')};

% preproc images
rundir = fullfile(preproc_dir, PID, 'ses-2', 'func');
in{2} = cellstr(spm_select('ExtFPList', rundir, strcat('.*task-MID_run-0',num2str(run),'_space-MNI152NLin2009cAsym_desc-preproc_bold.nii'), ndummies+1:9999));

if isempty(in{2}{1})
    warning('No preprocd functional found')
    return
end

% onset files
sub_str = PID(5:9);
in{3} = filenames(fullfile(timing_dir, strcat('run-', num2str(run)), strcat(sub_str,'.mat')));

if isempty(in{3})
    warning('No modeling found (behav data might be missing)')
    return
end

%% nuisance covs

% fmriprep output
confound_fname = filenames(fullfile(preproc_dir, PID, 'ses-2', 'func', '*task-MID*confounds*.tsv'));

% find the raw image file, for spike detection
rawrun = filenames(fullfile(raw_dir, PID, 'ses-2', 'func', strcat('*task-MID_run-0',num2str(run),'_bold.nii*')));

% get nuis covs
[Rfull, Rselected, n_spike_regs, FD] = make_nuisance_from_fmriprep_output(confound_fname{run}, rawrun, TR);%, 4);
save(fullfile(save_dir, strcat(sub_str, '_ses', num2str(ses), '_run', num2str(run), '.mat')), 'FD')

% choose which matrix to use
R = Rselected;

% its now possible that some of the spike regs are all zero, b/c the spikes
% were discarded in the step above. find all-zero regs and drop
R(:, ~any(table2array(R))) = [];
R = R(ndummies+1:end, :); %discard dummy vols

% put in SPM format: matrix called 'R', and 'names'
names = R.Properties.VariableNames;
R = table2array(R);

confoundFile = fullfile(fileparts(in{3}{1}),'MID_confounds.mat');
save(confoundFile,'R', 'names');
in{4} = {confoundFile};

% checks
if any(cellfun( @(x) isempty(x{1}), in))
    in
    error('Some input to the model is missing')
end


% check for SPM.mat and overwrite if needed
skip = 0;
if exist(fullfile(in{1}{1},'SPM.mat'),'file')
    if overwrite
        fprintf('\n\nWARNING: EXISTING SPM.MATAND BETA FILES WILL BE OVERWRITTEN\n%s\n\n',fullfile(in{1}{1},'SPM.mat'));
        rmdir(in{1}{1},'s');
    else
        fprintf('\n\nFirst levels already exist, wont ovewrite: %s\n\n',fullfile(in{1}{1},'SPM.mat'));
        skip=1;
    end
end


if ~skip
    % make dir for beta and contrast files
    if ~isdir(in{1}{1}), mkdir(in{1}{1}); end


    % run spm FL estimation
    cwd = pwd;
    job = 'MID_SPM_template.m';
    %%
    spm('defaults', 'FMRI')
    spm_jobman('serial',job,'',in{:});

    cd(cwd);
end

end