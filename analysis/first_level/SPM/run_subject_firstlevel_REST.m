% run first level models for REST task for a subject
% the URSI and whether to overwrite existing first levels.  Overwrite = 0
% (default) or 1. before running this file, must run
% read_timings_make_onsets.m

    
 
function run_subject_firstlevel_REST(ids, ses, run, overwrite)
%% var set up
if nargin==0 % defaults just for testing    
    ids = 21684;
    overwrite = 1;
    ses = 2;
    run = 1;
end


% directories
% first is where your stats files will be output to
directories{1} = '/projects/b1108/projects/BrainMAPD_func_conn/first_levels/sarah_first_level_output';
% next is where the preprocessed data is
directories{2} = '/projects/b1108/projects/BrainMAPD_func_conn/fmriprep';
% where the raw data lives (raw meaning before preprocessing)
directories{3} = '/projects/b1108/data/BrainMAPD';
% the timing files for modelling (onsets, durations, names)
directories{4} = '/projects/b1108/projects/BrainMAPD_func_conn/first_levels/sarah_first_level_output/timing_files';
% where framewise displacement files will be saved
directories{5} = '/projects/b1108/projects/BrainMAPD_func_conn/first_levels/sarah_first_level_output/framewise_displacement';

fl_dir = directories{1};
preproc_dir = directories{2};
raw_dir = directories{3};
save_dir = directories{5};

if nargin==1
    overwrite = 1;
end 

ids = strcat('sub-',num2str(ids));

fprintf(['Preparing 1st level model for REST task for ' ids ' / ' ses], ['Overwrite = ' num2str(overwrite)]);


ndummies = 2;
TR = .555;

%% Model for REST task

% FL directory for saving 1st level results: beta images, SPM.mat, etc.
in{1} = {fullfile(fl_dir, ids, 'ses-2', strcat('run-', num2str(run)), 'REST')};

% preproc images
rundir = fullfile(preproc_dir, ids, 'ses-2', 'func');
in{2} = cellstr(spm_select('ExtFPList', rundir, strcat('.*task-REST_run-0',num2str(run),'_space-MNI152NLin2009cAsym_desc-preproc_bold.nii'), ndummies+1:9999));
    
if isempty(in{2}{1})
    warning('No preprocd functional found')
    return
end

%% nuisance covs

% fmriprep output
confound_fname = filenames(fullfile(preproc_dir, ids, 'ses-2', 'func', '*task-REST*confounds*.tsv'));

% find the raw image file, for spike detection
rawrun = filenames(fullfile(raw_dir, ids, 'ses-2', 'func', strcat('*task-REST_run-0',num2str(run),'_bold.nii*')));
cd(fullfile(preproc_dir, ids));

% get nuis covs
sub_str = ids(5:9);

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

save(fullfile(fl_dir,'REST_confounds.mat'),'R','names');
confoundFile = fullfile(fl_dir,'REST_confounds.mat');
in{3} = {confoundFile};

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
    job = 'REST_SPM_template.m';
    %%
    spm('defaults', 'FMRI')
    spm_jobman('serial',job,'',in{:});

    cd(cwd);
end

end


