% Need to run PPI. This involves rerunning first level models but with an
% interaction variable that I create. So the idea is to load data, extract
% a time course for an ROI using some of Tor's scripts and then take the
% product of that time course and the original design matrix. Right now I'm
% going to leave the tasks separate but I need to control for all the
% different kinds of trials.


% set up directories
datadir = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD/func_conn/ICA/MID_data';
spmdir = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD/func_conn/first_levels/first_level_output/first_level_output/consumption';
maskdir = '/Users/zaz3744/Documents/current_projects/ACNlab/masks/ROI_BrainMAPD_functional/consumption';
timecoursedir = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD/func_conn/PPI/time_courses/bi_VS';
ppimdldir = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD/func_conn/PPI/mdl_dir/consumption';
tempderivativefldir = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD/MID_contrasts_w_temp_der/nuisance_regressors';

mask = 'OFC*Oldham'; 
%% First, extract mean timecourses from a variable of interest
fnames = filenames(fullfile(datadir, '*.nii'));
maskobj = fmri_data(filenames(fullfile(maskdir,strcat(mask,'*'))));

for sub = 1:length(fnames)
    curr_id = fnames{sub}(88:92);
    dataobj = fmri_data(filenames(fullfile(datadir,strcat('*',curr_id,'*'))));
    [roiobj] = extract_roi_averages(dataobj,maskobj);
    mean_signal = roiobj.dat;
    curr_save_file = strcat(mask,'_',curr_id,'_time_series.mat');
    save(fullfile(timecoursedir,curr_save_file),'mean_signal');   
end
    

%% Now we'll load up SPM files. 
% Then we'll load in the time course data from above to create the
% interaction regressors. I'll concatenate all these regressors into a
% single variable R, which I should be able to use to rerun first levels.
fnames = filenames(fullfile(datadir, '*.nii'));

clear sub curr_id names
for sub = 1:length(fnames)
    curr_id = fnames{sub}(88:92);
    spm_fname = filenames(fullfile(spmdir,strcat('*',curr_id,'/ses-2/run-1/MID/SPM*')));
    if isempty(spm_fname) == 0
        load(spm_fname{1});
        timecourse_fname = filenames(fullfile(timecoursedir,strcat('*',curr_id,'*')));
        load(timecourse_fname{1})
        mean_signal = zscore(mean_signal);
        ppi_regressor1 = SPM.xX.X(:,1) .* mean_signal(3:281,1);
        ppi_regressor2 = SPM.xX.X(:,2) .* mean_signal(3:281,1);
        R = [ppi_regressor1, ppi_regressor2, mean_signal(3:281,:), SPM.xX.X(:,6:size(SPM.xX.X,2))];
        names{1} = 'loss_ppi'; names{2} = 'gain_ppi'; names{3} = 'time_course'; names = [names,SPM.xX.name(1,6:length(SPM.xX.name))];
        curr_save_file = strcat(mask,'_',curr_id,'_ppi_regressors.mat');
        save(fullfile(ppimdldir,curr_save_file),'R', 'names')
        clear names
    end
end

% At this point, go to: run_subject_firstlevel_PPI_MID_consumption.m and
% MID_PPI_consumption_template.m
%% I will write the section of batch script here
% This will run all subjects through the above 2 scripts

% directories
% first is where your stats files will be output to
directories{1} = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD/func_conn/PPI/ppi_fldir/consumption/OFC';
% next is where the preprocessed data is
directories{2} = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD/func_conn/ICA/MID_data';
% the timing files for modelling (onsets, durations, names)
directories{3} = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD/func_conn/final_timing_files/run-1/consumption_timing';
% where your extra covariates are including PPI regressors
directories{4} = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD/func_conn/PPI/mdl_dir/consumption';

% What run of your task are you looking at?
run = 1;
% What session appears in your raw filenames when in BIDS format?
ses = 2;
% Do you want to overwrite previously estimated first levels or just add to
% what you have? 1 overwrites, 0 adds
overwrite = 0;
% specify which mask you're looking at. It should just be the first
% characters of the file
mask_string = 'OFC'; % OFC, VS, HO_VMPFC

fnames = filenames(fullfile(directories{2}, '*.nii'));

for sub = 1:length(fnames)
    PID = fnames{sub}(88:92);
    spm_fname = filenames(fullfile(directories{4},strcat('*',PID,'*')));
    if isempty(spm_fname) == 0
        
        run_subject_firstlevel_PPI_MID_consumption(PID, ses, run, mask_string, overwrite)
    end
end

%% Group analysis
curr_analysis_dir = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD/func_conn/PPI/ppi_fldir/consumption/OFC';
loss_ppi_fnames = filenames(fullfile(curr_analysis_dir,'*/*/*/*/con_0001.nii')); %con_0003.nii'));
gain_ppi_fnames = filenames(fullfile(curr_analysis_dir,'*/*/*/*/con_0002.nii')); %con_0004.nii'));

loss_data = fmri_data(loss_ppi_fnames);
gain_data = fmri_data(gain_ppi_fnames);

% Load in clinical data for group analysis
clinicaldir = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD/clinical_data';

load(fullfile(clinicaldir,'BrainMAPD_clinical_diagnoses_final.mat'))
for sub = 1:length(loss_ppi_fnames)
    PID(sub,1) = str2num(loss_ppi_fnames{sub}(100:104));
    if isempty(find(clinical_info.PID(:) == PID(sub,1))) == 0
        curr = find(clinical_info.PID(:) == PID(sub,1));
        life_dep(sub,1) = clinical_info.dep_life_any(curr);
        life_anx(sub,1) = clinical_info.anx_life_any(curr);
        life_com(sub,1) = clinical_info.comorbid_life_dep_anx(curr);
    else
        disp(strcat(loss_ppi_fnames{sub}(100:104), ' missing clinical info'))
        life_dep(sub,1) = 0;%NaN;
        life_anx(sub,1) = 0;%NaN;
        life_com(sub,1) = 0;%NaN;
    end
end

R = [life_anx,life_dep,life_com];
life_healthy = zeros(length(R),1);
life_healthy(find(sum(R,2)==0)) = 1;
R = [R,ones(length(life_dep),1)];
% second_level_regressors = array2table(second_level_regressors); second_level_regressors.Properties.VariableNames = {'anx','dep','com','healthy'};
% save(fullfile(savedir,'temp_second_level_regressors.mat'),'R')
%% Seed to seed?
% 

target = mean(ha_output,2);
group = ones(size(loss_ppi_fnames,1),1);
group(R(:,1)==1) = 2;
group(R(:,2)==1) = 3;
group(R(:,3)==1) = 4;

group_strings = cell(size(group));
group_strings(group(:,1)==1) = {'Healthy'};
group_strings(group(:,1)==2) = {'Anxiety'};
group_strings(group(:,1)==3) = {'Depression'};
group_strings(group(:,1)==4) = {'Comorbidity'};

[p,tbl,stats] = anova1(target,group_strings);

% Extract Trilevel information
% Gotta see if we replicate the Anhedonia findings Rick and Robin are
% publishing

%%
clinicaldir = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD/clinical_data';

clear sub
 
load(fullfile(clinicaldir,'trilevel_factors.mat'));
trilevel_array = [trilevel.ID,trilevel.GenDis,trilevel.Anhedon,trilevel.Fears];
for sub = 1:length(loss_ppi_fnames)
    PID(sub,1) = str2num(loss_ppi_fnames{sub}(100:104));
    if isempty(find(trilevel.ID(:) == PID(sub,1))) == 0
        curr = find(trilevel.ID(:) == PID(sub,1));
        trilevel_regressors(sub,:) = trilevel_array(curr,:);
    else
        % for OFC loss_ppi_fnames{sub}(100:104); for HO_VMPFC loss_ppi_fnames{sub}(105:109)
        disp(strcat(loss_ppi_fnames{sub}(100:104), ' missing clinical info')) 
        trilevel_regressors(sub,:) = NaN;
        trilevel_regressors(sub,1) = PID(sub,1);
        loss_data.dat(:,sub) = [];
        gain_data.dat(:,sub) = [];
    end
end

GenDis = trilevel_regressors(:,2);
Anhedonia = trilevel_regressors(:,3);
Fears = trilevel_regressors(:,4);

GenDis(any(isnan(GenDis), 2), :) = [];
Anhedonia(any(isnan(Anhedonia), 2), :) = [];
Fears(any(isnan(Fears), 2), :) = [];

trilevel_analysis_table = [mean(loss_data.dat,1)', mean(gain_data.dat,1)',GenDis,Anhedonia,Fears];
trilevel_analysis_table = array2table(trilevel_analysis_table); trilevel_analysis_table.Properties.VariableNames = {'loss_activation','gain_activation','GenDis','Anhedonia','Fears'};
% loss_data.X = R;
% regress(loss_data, .01, 'unc')
% gain_data.X = R;
% regress(gain_data, .01, 'unc')
% fitlm(trilevel_analysis_table, 'loss_activation~GenDis')

%% try something a little crazy. Going to hyperalign conn matrices
% and then test for diagnosis related differences
clear aligned transforms ha_input ha_output corr_mat_unaligned corr_mat_aligned

maskdir = '/Users/zaz3744/Documents/current_projects/ACNlab/masks/ICA_generated_masks';
mask = fmri_data(filenames(fullfile(maskdir,'MID_mask1.nii')));


mask_logical = logical(mask.dat);
dat = gain_data.dat(mask_logical,:);
dat = dat';

clear sub
for sub = 1:size(dat,1)
    ha_input{sub} = dat(sub,:);
end

[aligned, transforms] = hyperalign(ha_input{:});
clear sub
for sub = 1:length(aligned)
    ha_output(sub,:) = aligned{sub};
end

corr_mat_unaligned = corrcoef(dat(:,:)');
corr_mat_aligned = corrcoef(ha_output');

heatmap(corr_mat_unaligned)
figure();heatmap(corr_mat_aligned)




