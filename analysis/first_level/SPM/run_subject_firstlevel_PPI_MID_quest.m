% run first level models for sponpain task for a subject
% the URSI and whether to overwrite existing first levels.  Overwrite = 0
% (default) or 1. before running this file, must run
% read_timings_make_onsets.m

function run_subject_firstlevel_PPI_MID(PID, ses, run, overwrite)


%% var set up
if nargin==0 % defaults just for testing
    PID = '10004';  
    overwrite = 1;
    ses = 2;
    run = 1;
end


con_string = 'ant'; %'ant' or 'con'
job = 'MID_PPI_template.m'; %
fl_dir = '/projects/b1108/projects/BrainMAPD_func_analysis/PPI/ppi_fldir/anticipation/HOAmyg';
preproc_dir = '/projects/b1108/projects/BrainMAPD_func_analysis/fmriprep';
nuisance_ppi_dir = '/projects/b1108/projects/BrainMAPD_func_analysis/PPI/mdl_dir/anticipation/HOAmyg';
timing_dir = '/projects/b1108/projects/BrainMAPD_func_analysis/timing_files/final/run1';
save_dir = '/projects/b1108/projects/BrainMAPD_func_analysis/PPI/additional_files';
PID = num2str(PID);
if nargin==1
    overwrite = 1;
end  

%fprintf(['Preparing 1st level model for MID task for ' PID ' / ' num2str(ses)], ['Overwrite = ' num2str(overwrite)]);


ndummies = 2;
TR = 2.05;

%% Model for MID task

% FL directory for saving 1st level results: beta images, SPM.mat, etc.
in{1} = {fullfile(fl_dir, PID, strcat('ses-',num2str(ses)), strcat('run-', num2str(run)), 'MID')};

% preproc images
rundir = fullfile(preproc_dir, strcat('sub-',PID,'/ses-',num2str(ses),'/func'));
in{2} = cellstr(spm_select('ExtFPList', rundir, strcat('ssub-',PID,'.*task-MID_run-0',num2str(run),'_space-MNI152NLin2009cAsym_desc-preproc_bold.nii'), ndummies+1:9999));

if isempty(in{2}{1})
    warning('No preprocd functional found')
    return
end

% onset files
in{3} = filenames(fullfile(timing_dir, strcat(PID,'*',con_string,'*.mat')));

if isempty(in{3})
    warning('No modeling found (behav data might be missing)')
    return
end

%% nuisance covs

% regressors = filenames(fullfile(mdl_dir,strcat(mask_string, '*', PID, '*.mat')));
regressors = filenames(fullfile(nuisance_ppi_dir,strcat('*', PID, '*.mat')));
in{4} = regressors;

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
    % job = 'MID_PPI_template.m';
    %%
    spm('defaults', 'FMRI')
    spm_jobman('serial',job,'',in{:});

    cd(cwd);
end

end
