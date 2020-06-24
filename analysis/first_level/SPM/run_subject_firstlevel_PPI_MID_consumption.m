% run first level models for sponpain task for a subject
% the URSI and whether to overwrite existing first levels.  Overwrite = 0
% (default) or 1. before running this file, must run
% read_timings_make_onsets.m

function run_subject_firstlevel_MID_consumption(PID, ses, run, mask_string, overwrite)


%% var set up
if nargin==0 % defaults just for testing
    PID = 20695;  
    overwrite = 1;
    ses = 2;
    mask_string = 'HO_VMPFC';
    run = 1;
end

% directories
directories{1} = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD/func_conn/first_levels/first_level_output';
% next is where the preprocessed data is
directories{2} = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD/func_conn/ICA/MID_data';
% the timing files for modelling (onsets, durations, names)
directories{3} = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD/func_conn/final_timing_files/run-1/anticipation/spm_all_vs_0_timing';
% where your extra covariates are including PPI regressors
% directories{4} = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD/func_conn/PPI/mdl_dir/consumption';
% rerunning original first levels with temporal derivatives because of
% rapid event related design (Wager et al, 2005)
directories{4} = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD/MID_contrasts_w_temp_der/nuisance_regressors';

fl_dir = directories{1};
preproc_dir = directories{2};
timing_dir = directories{3};
mdl_dir = directories{4};

if nargin==1
    overwrite = 1;
end  

fprintf(['Preparing 1st level model for MID task for ' PID ' / ' ses], ['Overwrite = ' num2str(overwrite)]);


ndummies = 2;
TR = 2;

%% Model for MID task

% FL directory for saving 1st level results: beta images, SPM.mat, etc.
in{1} = {fullfile(fl_dir, PID, strcat('ses-',num2str(ses)), strcat('run-', num2str(run)), 'MID')};

% preproc images
rundir = fullfile(preproc_dir);
in{2} = cellstr(spm_select('ExtFPList', rundir, strcat('ssub-',PID,'.*task-MID_run-0',num2str(run),'_space-MNI152NLin2009cAsym_desc-preproc_bold.nii'), ndummies+1:9999));

if isempty(in{2}{1})
    warning('No preprocd functional found')
    return
end

% onset files
in{3} = filenames(fullfile(timing_dir, strcat('consumption_',PID,'.mat')));

if isempty(in{3})
    warning('No modeling found (behav data might be missing)')
    return
end

%% nuisance covs

% regressors = filenames(fullfile(mdl_dir,strcat(mask_string, '*', PID, '*.mat')));
regressors = filenames(fullfile(mdl_dir,strcat('*', PID, '*.mat')));
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
    %job = 'MID_PPI_consumption_template.m';
    job = 'MID_SPM_consumption_template.m';
    %%
    spm('defaults', 'FMRI')
    spm_jobman('serial',job,'',in{:});

    cd(cwd);
end

end