% run first level models for sponpain task for a subject
% the URSI and whether to overwrite existing first levels.  Overwrite = 0
% (default) or 1. before running this file, must run
% read_timings_make_onsets.m

function run_subject_firstlevel_MID_quest(PID)


%% var set up
%PID = 10004;  
overwrite = 1;
ses = 2;
run = 2;
contrast1 = 'ant';
contrast2 = 'con'; %'ant' or 'con'
    
job1 = 'SPM_MID_anticipation_quest_template.m';    
job2 = 'SPM_MID_consumption_quest_template.m'; %'SPM_MID_anticipation_quest_template.m' or 'SPM_MID_consumption_quest_template'   
fl_dir1 = '/projects/b1108/projects/BrainMAPD_func_analysis/first_levels/first_level_output/run-2/anticipation';
fl_dir2 = '/projects/b1108/projects/BrainMAPD_func_analysis/first_levels/first_level_output/run-2/consumption';
preproc_dir = '/projects/b1108/projects/BrainMAPD_func_analysis/fmriprep';
raw_dir = '/projects/b1108/data/BrainMAPD';
timing_dir = '/projects/b1108/projects/BrainMAPD_func_analysis/timing_files/final/';
save_dir = '/projects/b1108/projects/BrainMAPD_func_analysis/first_levels/additional_files';

if nargin==1
    overwrite = 1;
end  

PID = strcat('sub-',num2str(PID));

%fprintf(['Preparing 1st level model for MID task for ' PID ' / ' ses], ['Overwrite = ' num2str(overwrite)]);


ndummies = 2;
TR = 2.05;

%% Model for MID task

% FL directory for saving 1st level results: beta images, SPM.mat, etc.
a_in{1} = {fullfile(fl_dir1, PID, strcat('ses-',num2str(ses)), strcat('run-', num2str(run)), 'MID')};
c_in{1} = {fullfile(fl_dir2, PID, strcat('ses-',num2str(ses)), strcat('run-', num2str(run)), 'MID')};
% preproc images
rundir = fullfile(preproc_dir, PID, strcat('ses-', num2str(ses)), 'func');
a_in{2} = cellstr(spm_select('ExtFPList', rundir, strcat('ssub.*task-MID_run-0',num2str(run),'_space-MNI152NLin2009cAsym_desc-preproc_bold.nii'), ndummies+1:9999));
c_in{2} = cellstr(spm_select('ExtFPList', rundir, strcat('ssub.*task-MID_run-0',num2str(run),'_space-MNI152NLin2009cAsym_desc-preproc_bold.nii'), ndummies+1:9999));

if isempty(a_in{2}{1})
    warning('No preprocd functional found')
    return
end

% onset files
sub_str = PID(5:9);
a_in{3} = filenames(fullfile(timing_dir, strcat('run', num2str(run)), strcat(sub_str,'*',contrast1,'*.mat')));
c_in{3} = filenames(fullfile(timing_dir, strcat('run', num2str(run)), strcat(sub_str,'*',contrast2,'*.mat')));
if isempty(a_in{3})
    warning('No modeling found (behav data might be missing)')
    return
end

%% nuisance covs

% fmriprep output
confound_fname = filenames(fullfile(preproc_dir, PID, strcat('ses-',num2str(ses)), 'func', '*task-MID*confounds*.tsv'));

% find the raw image file, for spike detection
rawrun = filenames(fullfile(raw_dir, PID, strcat('ses-',num2str(ses)), 'func', strcat('*task-MID_run-0',num2str(run),'_bold.nii*')));
cd(fullfile(preproc_dir, PID));
% get nuis covs
[Rfull, Rselected, n_spike_regs, FD, gsr_final] = make_nuisance_from_fmriprep_output_MID(confound_fname{run}, rawrun, TR);%, 4);
save(fullfile(save_dir, strcat('FD_',sub_str, '_ses', num2str(ses), '_run', num2str(run), '.mat')), 'FD')
save(fullfile(save_dir, strcat('GSR_',sub_str, '_ses', num2str(ses), '_run', num2str(run), '.mat')), 'gsr_final')
% choose which matrix to use
R = Rselected;

% its now possible that some of the spike regs are all zero, b/c the spikes
% were discarded in the step above. find all-zero regs and drop
R(:, ~any(table2array(R))) = [];
R = R(ndummies+1:end, :); %discard dummy vols

% put in SPM format: matrix called 'R', and 'names'
names = R.Properties.VariableNames;
R = table2array(R);

confoundFile = fullfile(fileparts(a_in{3}{1}),'MID_confounds.mat');
save(confoundFile,'R', 'names');
a_in{4} = {confoundFile};
c_in{4} = {confoundFile};
% checks
if any(cellfun( @(x) isempty(x{1}), a_in))
    a_in
    error('Some input to the model is missing')
end


% check for SPM.mat and overwrite if needed
skip = 0;
if exist(fullfile(a_in{1}{1},'SPM.mat'),'file')
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
    if ~isdir(a_in{1}{1}), mkdir(a_in{1}{1}); end
    if ~isdir(c_in{1}{1}), mkdir(c_in{1}{1}); end

    % run spm FL estimation
    cwd = pwd;
    %job = 'SPM_MID_anticipation_quest_template.m';
    %%
    spm('defaults', 'FMRI')
    spm_jobman('serial',job1,'',a_in{:});
    spm_jobman('serial',job2,'',c_in{:});

    cd(cwd);
end

end