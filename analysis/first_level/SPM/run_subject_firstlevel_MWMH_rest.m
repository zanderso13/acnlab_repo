% run first level models for REST task for a subject
% the URSI and whether to overwrite existing first levels.  Overwrite = 0
% (default) or 1. before running this file, must run
% read_timings_make_onsets.m

    
 
function run_subject_firstlevel_REST(PID, ses, run, overwrite,directories)
%% var set up
if nargin==0 % defaults just for testing 
    % Define some paths
    basedir = '/Users/zaz3744/Documents/current_projects/ACNlab/MWMH/rest';
    % directories
    % first is where your stats files will be output to
    directories{1} = fullfile(basedir,'first_levels');
    % next is where the preprocessed data is
    directories{2} = fullfile(basedir,'fmriprep');
    % where the raw data lives (raw meaning before preprocessing)
    % directories{3} = '/projects/b1108/data/BrainMAPD';
    directories{3} = '/Users/zaz3744/Documents/current_projects/ACNlab/MWMH/rest/raw';
    % where framewise displacement files will be saved
    directories{4} = fullfile(basedir,'first_levels/FD');
    PID = 113;
    overwrite = 1;
    ses = 1;
    run = 1;
end

fl_dir = directories{1};
preproc_dir = directories{2};
raw_dir = directories{3};
save_dir = directories{4};

if nargin==1
    overwrite = 1;
end 

PID = strcat('sub-MWMH',num2str(PID));

fprintf(['Preparing 1st level model for REST task for ' PID ' / ' ses], ['Overwrite = ' num2str(overwrite)]);


ndummies = 8;
TR = .555;

%% Model for REST task

% FL directory for saving 1st level results: beta images, SPM.mat, etc.
in{1} = {fullfile(fl_dir, PID, strcat('ses-',num2str(ses)), strcat('run-', num2str(run)), 'rest')};

% preproc images
rundir = fullfile(preproc_dir, PID, strcat('ses-',num2str(ses)), 'func');
in{2} = cellstr(spm_select('ExtFPList', rundir, strcat('.*task-rest_run-',num2str(run),'_space-MNI152NLin2009cAsym_desc-preproc_bold.nii'), ndummies+1:9999));
    
if isempty(in{2}{1})
    warning('No preprocd functional found')
    return
end

%% nuisance covs

% fmriprep output
confound_fname = filenames(fullfile(preproc_dir, strcat(PID), strcat('ses-',num2str(ses)), 'func', '*task-rest*confounds*.tsv'));

% find the raw image file, for spike detection
rawrun = filenames(fullfile(raw_dir, PID, strcat('ses-',num2str(ses)), 'func', strcat('*task-rest_run-0',num2str(run),'_bold.nii*')));
cd(fullfile(preproc_dir, PID));

% get nuis covs

[Rfull, Rselected, n_spike_regs, FD] = make_nuisance_from_fmriprep_output_restMWMH(confound_fname{run}, rawrun, TR, 4);
save(fullfile(save_dir, strcat(PID, '_ses', num2str(ses), '_run', num2str(run), '.mat')), 'FD')


% choose which matrix to use
R = Rselected;

% its now possible that some of the spike regs are all zero, b/c the spikes
% were discarded in the step above. find all-zero regs and drop
R(:, ~any(table2array(R))) = [];
R = R(ndummies+1:end, :); %discard dummy vols

% put in SPM format: matrix called 'R', and 'names'
names = R.Properties.VariableNames;
R = table2array(R);

save(fullfile(preproc_dir,'rest_confounds.mat'),'R','names');
confoundFile = fullfile(preproc_dir,'rest_confounds.mat');
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
    job = 'MWMH_rest_SPM_template.m';
    %%
    spm('defaults', 'FMRI')
    spm_jobman('serial',job,'',in{:});

    cd(cwd);
end

end


