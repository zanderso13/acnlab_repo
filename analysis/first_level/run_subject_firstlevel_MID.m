% run first level models for sponpain task for a subject
% the URSI and whether to overwrite existing first levels.  Overwrite = 0
% (default) or 1. before running this file, must run
% read_timings_make_onsets.m

function run_subject_firstlevel_MID(PID, ses, overwrite)


%% var set up
if nargin==0 % defaults just for testing
    URSI = 'sub-M80303741';    
    overwrite = 1;
    ses = 1;
end

if nargin==1
    overwrite = 1;
end


    
% subject id    
[subID, ~, ~, sub_str, ses_str] = get_subjID_from_URSI(URSI);    



fprintf(['Preparing 1st level model for SPONPAIN task for ' URSI ' / ' subID], ['Overwrite = ' num2str(overwrite)]);

% directories
fl_dir = '/home/zaz3744/ACNlab/projects/BrainMAPD_func_conn/first_levels';
preproc_dir = '/home/zaz3744/ACNlab/projects/BrainMAPD_func_conn/fmriprep';
raw_dir = '/home/zaz3744/ACNlab/data/BrainMAPD';
save_dir = '/home/zaz3744/repo/acnlab_repo/analysis/first_level';

ndummies = 2;
TR = 2;

%% Model for MID task

% FL directory for saving 1st level results: beta images, SPM.mat, etc.
in{1} = {fullfile(fl_dir, sub_str, ses_str(ses,:), 'MID')};

% preproc images
rundir = fullfile(preproc_dir, URSI, ses_str(ses,[1:4 6]), 'func');
in{2} = cellstr(spm_select('ExtFPList', rundir, '.*task-MID_space-MNI152NLin2009cAsym_desc-preproc_bold.nii', ndummies+1:9999));

if isempty(in{2}{1})
    warning('No preprocd functional found')
    return
end

% onset files
in{3} = filenames(fullfile(fl_dir,sub_str, ses_str(ses,:), 'spm_modeling', 'sponpain_ratings_model.mat'));

if isempty(in{3})
    warning('No modeling found (behav data might be missing)')
    return
end

%% nuisance covs

% fmriprep output
confound_fname = filenames(fullfile(preproc_dir, URSI, ses_str(ses,[1:4 6]), 'func', '*task-sponpain*confounds*.tsv'), 'char');

% find the raw image file, for spike detection
rawrun = filenames(fullfile(raw_dir, URSI, ses_str(ses,[1:4 6]), 'func', '*task-sponpain_bold*.nii*'));

% get nuis covs
[Rfull, Rselected, n_spike_regs, FD] = make_nuisance_covs_from_fmriprep_output_ZA(confound_fname, rawrun, TR);%, 4);
save(fullfile(save_dir, strcat(sub_str, '_ses', num2str(ses), '.mat')), 'n_spike_regs', 'FD')

% choose which matrix to use
R = Rselected;

% its now possible that some of the spike regs are all zero, b/c the spikes
% were discarded in the step above. find all-zero regs and drop
R(:, ~any(table2array(R))) = [];
R = R(ndummies+1:end, :); %discard dummy vols

% put in SPM format: matrix called 'R', and 'names'
names = R.Properties.VariableNames;
R = table2array(R);

confoundFile = fullfile(fileparts(in{3}{1}),'sponpain_confounds.mat');
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
    job = 'subject_fl_spm_batch_template_sponpain.m';
    %%
    spm('defaults', 'FMRI')
    spm_jobman('serial',job,'',in{:});

    cd(cwd);
end

end