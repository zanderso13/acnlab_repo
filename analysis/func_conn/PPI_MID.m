% Need to run PPI. This involves rerunning first level models but with an
% interaction variable that I create. So the idea is to load data, extract
% a time course for an ROI using some of Tor's scripts and then take the
% product of that time course and the original design matrix. Right now I'm
% going to leave the tasks separate but I need to control for all the
% different kinds of trials.


% set up directories
datadir = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD/func_conn/ICA/MID_data';
spmdir = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD/func_conn/first_levels/first_level_output/anticipation';
maskdir = '/Users/zaz3744/Documents/current_projects/ACNlab/masks/ROI_BrainMAPD_functional/consumption';
%timecoursedir = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD/func_conn/PPI/time_courses/bi_VS';
%ppimdldir = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD/func_conn/PPI/mdl_dir/anticipation';
tempderivativefldir = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD/func_conn/first_levels/nuisance_regressors';

mask = 'LOFC'; % OFC, VS, VMPFC (gotta change line 89)

%% Quick motion check
% Hopefully this will be a one and done
datadir = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD/func_conn/ICA/MID_data';
motiondir = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD/func_conn/first_levels/additional_files';
fnames = filenames(fullfile(datadir, '*.nii'));

clear sub curr_id names
for sub = 1:length(fnames)
    curr_id = fnames{sub}(88:92);
    motion_fname = filenames(fullfile(motiondir,strcat('FD*',curr_id,'*')));
    if isempty(motion_fname) == 0
        load(motion_fname{1});
        if isempty(FD)
            disp(strcat('missing data:', num2str(curr_id)))
            continue
        else
            FD_final(:,sub) = table2array(FD(3:281,1));
            if mean(FD_final(:,sub)) > 0.3
                disp(strcat('motion problems:',num2str(curr_id)))
                mean(FD_final(:,sub))
            end
        end
    end
end


%% This section was written to redo first levels based on new timing files
% It's adapted from the original stuff I wrote to perform PPI but I needed
% to change it in order to create new regressors for rerunning first levels
% round 1 Going to copy and paste it in a new section below for PPI
fnames = filenames(fullfile(datadir, '*.nii'));

clear sub curr_id names
for sub = 1:length(fnames)
    curr_id = fnames{sub}(88:92);
    spm_fname = filenames(fullfile(spmdir,strcat('*',curr_id,'/ses-2/run-1/MID/SPM*')));
    if isempty(spm_fname) == 0
        load(spm_fname{1});
        %timecourse_fname = filenames(fullfile(timecoursedir,strcat('*',curr_id,'*')));
        %load(timecourse_fname{1})
        %mean_signal = zscore(mean_signal);
        %ppi_regressor1 = SPM.xX.X(:,1) .* mean_signal(3:281,1);
        %ppi_regressor2 = SPM.xX.X(:,2) .* mean_signal(3:281,1);
        %R = [ppi_regressor1, ppi_regressor2, mean_signal(3:281,:), SPM.xX.X(:,6:size(SPM.xX.X,2))];
        R = [SPM.xX.X(:,6:size(SPM.xX.X,2))];
        %names{1} = 'loss_ppi'; names{2} = 'gain_ppi'; names{3} = 'time_course'; 
        name = [SPM.xX.name(1,6:length(SPM.xX.name))];
        curr_save_file = strcat(curr_id,'_nuisance.mat'); % include mask name in this too, ppi_regressors.mat');
        
        save(fullfile(ppimdldir,curr_save_file),'R', 'names')
        clear names
    end
end

% At this point, go to: run_subject_firstlevel_PPI_MID_consumption.m and
% MID_PPI_consumption_template.m
%% I will write the section of batch script here
% Rerun first level models, NOT PPI
% You don't need to run the above sections, this stands alone. It makes
% sense to keep it though. We'll make the assumption that timing files will
% need to change again/there are some final subjects that Atrik is
% preprocessing

% first is where your stats files will be output to
directories{1} = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD/func_conn/first_levels/first_level_output/consumption';
% next is where the preprocessed data is
directories{2} = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD/func_conn/ICA/MID_data';
% the timing files for modelling (onsets, durations, names)
directories{3} = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD/func_conn/final_timing_files/run-1/consumption/spm_all_vs_0_timing';
% where your extra covariates are including PPI regressors
directories{4} = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD/func_conn/first_levels/nuisance_regressors';

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
        
        run_subject_firstlevel_MID_Mac(PID, ses, run, mask_string, overwrite)
    end
end

%% First, extract mean timecourses from a variable of interest
datadir = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD/func_conn/ICA/MID_data';
maskdir = '/Users/zaz3744/Documents/current_projects/ACNlab/masks/ROI_BrainMAPD_functional/anticipation/right';
timecoursedir = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD/func_conn/PPI/time_courses/R_VS';
fnames = filenames(fullfile(datadir, '*.nii'));
maskobj = fmri_data(filenames(fullfile(maskdir,strcat('R_VS_8mmsphere_Oldham_Loss.nii'))));

for sub = 1:length(fnames)
    curr_id = fnames{sub}(88:92);
    check_fname = filenames(fullfile(timecoursedir, strcat('R_VS_LossAnticipation_',curr_id,'_time_series.mat')));
    if isempty(check_fname) == 1
        dataobj = fmri_data(filenames(fullfile(datadir,strcat('*',curr_id,'*'))));
        [roiobj] = extract_roi_averages(dataobj,maskobj);
        mean_signal = roiobj.dat;
        curr_save_file = strcat('R_VS_LossAnticipation_',curr_id,'_time_series.mat');
        save(fullfile(timecoursedir,curr_save_file),'mean_signal');
    else
        continue
    end
end
    

%% Need to create regressor files for PPI: ROI - ROI
% This will be from a single ROI to whole brain
% Involves loading up SPM.mat files and then adding columns for 
fnames = filenames(fullfile(datadir, '*.nii'));
ppimdldir = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD/func_conn/PPI/mdl_dir/anticipation/L_VS_AntRew_wholebrain';
spmdir = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD/func_conn/first_levels/first_level_output/consumption';

timecoursedir1 = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD/func_conn/PPI/time_courses/L_VS';
%timecoursedir2 = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD/func_conn/PPI/time_courses/bi_VS';

clear sub curr_id names
for sub = 1:length(fnames)
    curr_id = fnames{sub}(88:92);
    spm_fname = filenames(fullfile(spmdir,strcat('*',curr_id,'/ses-2/run-1/MID/SPM*')));
    if isempty(spm_fname) == 0
        load(spm_fname{1});
        timecourse_fname1 = filenames(fullfile(timecoursedir1,strcat('*Rew*',curr_id,'*')));
        %timecourse_fname2 = filenames(fullfile(timecoursedir2,strcat('*',curr_id,'*')));
        roi1 = load(timecourse_fname1{1});
        %roi2 = load(timecourse_fname2{1});
        mean_signal1 = zscore(roi1.mean_signal);
        %mean_signal2 = zscore(roi2.mean_signal);
        ppi_regressor1 = SPM.xX.X(:,1) .* mean_signal1(3:281,1); % gain any, both for ant and con
        ppi_regressor2 = SPM.xX.X(:,2) .* mean_signal1(3:281,1); % loss, both for ant and con
        ppi_regressor3 = SPM.xX.X(:,3) .* mean_signal1(3:281,1); % gain 0, both for ant and con
        ppi_regressor4 = SPM.xX.X(:,4) .* mean_signal1(3:281,1); % loss 0, both for ant and con
        R = [ppi_regressor1, ppi_regressor2, ppi_regressor3, ppi_regressor4, mean_signal1(3:281,:), SPM.xX.X(:,6:size(SPM.xX.X,2))];
        %R = [SPM.xX.X(:,6:size(SPM.xX.X,2))];
        names{1} = 'gain_ppi_any'; names{2} = 'loss_ppi_any'; names{3} = 'gain_ppi_0'; names{4} = 'loss_ppi_0'; names{5} = 'timecourse_roi1'; 
        names = [names,SPM.xX.name(1,6:length(SPM.xX.name))];
        curr_save_file = strcat(curr_id,'_nuisance.mat'); % include mask name in this too, ppi_regressors.mat');
        
        save(fullfile(ppimdldir,curr_save_file),'R', 'names')
        clear names
    else
        disp(curr_id)
        continue
    end
end

% At this point, go to: run_subject_firstlevel_PPI_MID_consumption.m and
% MID_PPI_consumption_template.m


%% So now let's get back to PPI
% PPI first involves the estimation of task regressors which is done in the
% section above. Now I need to rerun first levels with two additional
% regressors 1. The time avg time course from a target ROI 2. An
% interaction variable avg time course * original task regressors.
% Ultimately this means this script is the exact same as above, the only
% differences is what regressors I'm calling

% first is where your stats files will be output to
directories{1} = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD/func_conn/PPI/ppi_fldir/anticipation/L_VS_AntRew_to_wholebrain';
% next is where the preprocessed data is
directories{2} = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD/func_conn/ICA/MID_data';
% the timing files for modelling (onsets, durations, names)
directories{3} = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD/func_conn/final_timing_files/run-1/anticipation/spm_all_vs_0_timing';
% where your extra covariates are including PPI regressors
directories{4} = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD/func_conn/PPI/mdl_dir/anticipation/L_VS_AntRew_wholebrain';

% What run of your task are you looking at?
run = 1;
% What session appears in your raw filenames when in BIDS format?
ses = 2;
% Do you want to overwrite previously estimated first levels or just add to
% what you have? 1 overwrites, 0 adds
overwrite = 0;
% specify which mask you're looking at. It should just be the first
% characters of the file
mask_string = 'L_VS_AntRew'; % OFC, VS, HO_VMPFC

fnames = filenames(fullfile(directories{2}, '*.nii'));

for sub = 1:length(fnames)
    PID = fnames{sub}(88:92);
    spm_fname = filenames(fullfile(directories{4},strcat('*',PID,'*')));
    if isempty(spm_fname) == 0
        
        run_subject_firstlevel_PPI_MID(PID, ses, run, mask_string, directories, overwrite)
    end
end

%% Fast Group analysis (see next section for analysis pulling clusters from deep storage)
curr_analysis_dir = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD/func_conn/PPI/ppi_fldir/anticipation';
gain_fnames = filenames(fullfile(curr_analysis_dir,'LHO_Accumbens_to_wholebrain/*/*/*/*/con_0001.nii')); %con_0003.nii'));
loss_fnames = filenames(fullfile(curr_analysis_dir,'LHO_Accumbens_to_wholebrain/*/*/*/*/con_0002.nii')); %con_0004.nii'));

loss_data = fmri_data(loss_fnames);
gain_data = fmri_data(gain_fnames);

% Load in clinical data for group analysis
clinicaldir = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD/clinical_data';

load(fullfile(clinicaldir,'BrainMAPD_clinical_diagnoses_final.mat'))
for sub = 1:length(loss_fnames)
    PID(sub,1) = str2num(loss_fnames{sub}(124:128));
    if isempty(find(clinical_info.PID(:) == PID(sub,1))) == 0
        curr = find(clinical_info.PID(:) == PID(sub,1));
        life_dep(sub,1) = clinical_info.dep_life_any(curr);
        life_anx(sub,1) = clinical_info.anx_life_any(curr);
        life_com(sub,1) = clinical_info.comorbid_life_dep_anx(curr);
    else
        disp(strcat(loss_fnames{sub}(124:128), ' missing clinical info'))
        life_dep(sub,1) = 0;%NaN;
        life_anx(sub,1) = 0;%NaN;
        life_com(sub,1) = 0;%NaN;
    end
end

R = [life_anx,life_dep,life_com];
life_healthy = zeros(length(R),1);
life_healthy(find(sum(R,2)==0)) = 1;
%R = [R,ones(length(life_dep),1)];
gain_data.X = R;
loss_data.X = R;
% regress(gain_data,'robust')
% regress(loss_data,'robust')
% second_level_regressors = array2table(second_level_regressors); second_level_regressors.Properties.VariableNames = {'anx','dep','com','healthy'};
% save(fullfile(savedir,'temp_second_level_regressors.mat'),'R')
%% Pull all results from deep storage DSM ANALYSIS
% Hope my external hard drive is ok with me pull all the data I'm about to.

% Robust_results1.mat: LHO_Accumbens, LOFC2, LOFC, RHO_Accumbens ROFC
% Robust_results2.mat: R/L VS Loss/Rew Anticipation

storedir = '/Volumes/ZachExternal/ACNlab/BrainMAPD/func_conn/PPI/ppi_fldir';
contrastdir = 'anticipation';
curr_analysis_dir = {'L_VS_AntLoss_to_wholebrain','L_VS_AntRew_to_wholebrain','R_VS_AntLoss_to_wholebrain','R_VS_AntRew_to_wholebrain'};%{'ROFC_to_wholebrain','RHO_Accumbens_to_wholebrain','LOFC_to_wholebrain','LOFC2_to_wholebrain','LHO_Accumbens_to_wholebrain'};

for roi = 1:length(curr_analysis_dir)
    curr_fnames_gain = filenames(fullfile(storedir,contrastdir,curr_analysis_dir{roi},'*/*/*/*/con_0001.nii'));
    curr_fnames_loss = filenames(fullfile(storedir,contrastdir,curr_analysis_dir{roi},'*/*/*/*/con_0002.nii'));
    
    curr_dat_gain = fmri_data(curr_fnames_gain);
    curr_dat_loss = fmri_data(curr_fnames_loss);
    
    % Load in clinical data for group analysis
    clinicaldir = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD/clinical_data';

    load(fullfile(clinicaldir,'BrainMAPD_clinical_diagnoses_final.mat'))
    if roi == 1
        for sub = 1:length(curr_fnames_gain)
            PID(sub,1) = str2num(curr_fnames_gain{sub}(104:108));
            if isempty(find(clinical_info.PID(:) == PID(sub,1))) == 0
                curr = find(clinical_info.PID(:) == PID(sub,1));
                life_dep(sub,1) = clinical_info.dep_life_any(curr);
                life_anx(sub,1) = clinical_info.anx_life_any(curr);
                life_com(sub,1) = clinical_info.comorbid_life_dep_anx(curr);
            else
                disp(strcat(curr_fnames_gain{sub}(104:108), ' missing clinical info'))
                life_dep(sub,1) = 0;%NaN;
                life_anx(sub,1) = 0;%NaN;
                life_com(sub,1) = 0;%NaN;
            end
        end        
    end
    
    R = [life_anx,life_dep,life_com];
    
    curr_dat_gain.X = R;
    curr_dat_loss.X = R;
    
    temp_results_gain = regress(curr_dat_gain,'robust');
    temp_results_loss = regress(curr_dat_loss,'robust');
    
    results_struct.gain.(curr_analysis_dir{roi}) = threshold(temp_results_gain.t,.001,'unc','k',50);
    results_struct.loss.(curr_analysis_dir{roi}) = threshold(temp_results_loss.t,.001,'unc','k',50);
    
    
end


%% quick visualization loop: DSM
analysis_dir = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD/func_conn/PPI/results/DSM_analysis/anticipation';
analysis_type = 'Robust'; % Also have OLS

curr_analysis_dir = {'L_VS_AntLoss_to_wholebrain','L_VS_AntRew_to_wholebrain','R_VS_AntLoss_to_wholebrain','R_VS_AntRew_to_wholebrain'};%{'ROFC_to_wholebrain','RHO_Accumbens_to_wholebrain','LOFC_to_wholebrain','LOFC2_to_wholebrain','LHO_Accumbens_to_wholebrain'};

load_fname = filenames(fullfile(analysis_dir,strcat(analysis_type,'*2.mat')));
load(load_fname{1})

for roi = 1:length(curr_analysis_dir)
    curr_analysis_dir{roi}
    cluster_out_gain = orthviews(results_struct.gain.(curr_analysis_dir{roi}))
    
    cluster_out_loss = orthviews(results_struct.loss.(curr_analysis_dir{roi}))
    fname = strcat(curr_analysis_dir{roi},'.mat');
    % save(fullfile(analysis_dir,fname),'cluster_out_gain','cluster_out_loss')
    data(1,:) = [cluster_out_gain{1}(1).mm_center,cluster_out_gain{1}(1).numVox,mean(cluster_out_gain{1}(1).Z(:))]
    data(2,:) = [cluster_out_gain{2}(1).mm_center,cluster_out_gain{2}(1).numVox,mean(cluster_out_gain{2}(1).Z(:))]
    data(3,:) = [cluster_out_gain{3}(1).mm_center,cluster_out_gain{3}(1).numVox,mean(cluster_out_gain{3}(1).Z(:))]
    data(4,:) = [cluster_out_loss{1}(1).mm_center,cluster_out_loss{1}(1).numVox,mean(cluster_out_loss{1}(1).Z(:))]
    data(5,:) = [cluster_out_loss{2}(1).mm_center,cluster_out_loss{2}(1).numVox,mean(cluster_out_loss{2}(1).Z(:))]
    data(6,:) = [cluster_out_loss{3}(1).mm_center,cluster_out_loss{3}(1).numVox,mean(cluster_out_loss{3}(1).Z(:))]
    keyboard
end

%% Pull all results from deep storage TRILEVEL ANALYSIS
% Hope my external hard drive is ok with me pull all the data I'm about to.

storedir = '/Volumes/ZachExternal/ACNlab/BrainMAPD/func_conn/PPI/ppi_fldir';
contrastdir = 'anticipation';
curr_analysis_dir = {'L_VS_AntLoss_to_wholebrain','L_VS_AntRew_to_wholebrain','R_VS_AntLoss_to_wholebrain','R_VS_AntRew_to_wholebrain'};%{'ROFC_to_wholebrain','RHO_Accumbens_to_wholebrain','LOFC_to_wholebrain','LOFC2_to_wholebrain','LHO_Accumbens_to_wholebrain'};

for roi = 1:length(curr_analysis_dir)
    clear curr_dat_gain curr_dat_loss trilevel_regressors GenDis Anhedonia Fears 
    curr_fnames_gain = filenames(fullfile(storedir,contrastdir,curr_analysis_dir{roi},'*/*/*/*/con_0001.nii'));
    curr_fnames_loss = filenames(fullfile(storedir,contrastdir,curr_analysis_dir{roi},'*/*/*/*/con_0002.nii'));
    
    curr_dat_gain = fmri_data(curr_fnames_gain);
    curr_dat_loss = fmri_data(curr_fnames_loss);
    
    % Load in clinical data for group analysis
    clinicaldir = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD/clinical_data';

    load(fullfile(clinicaldir,'trilevel_factors.mat'));
    trilevel_array = [trilevel.ID,trilevel.GenDis,trilevel.Anhedon,trilevel.Fears];
    if roi == 1
        for sub = 1:length(curr_fnames_gain)
            PID(sub,1) = str2num(curr_fnames_gain{sub}(95:99));
        end        
    end
    
    for sub = 1:length(curr_fnames_loss)
        
        if isempty(find(trilevel.ID(:) == PID(sub,1))) == 0
            curr = find(trilevel.ID(:) == PID(sub,1));
            trilevel_regressors(sub,:) = trilevel_array(curr,:);
        else
            % for OFC loss_ppi_fnames{sub}(100:104); for HO_VMPFC loss_ppi_fnames{sub}(105:109)
            disp(strcat(num2str(PID(sub,1)), ' missing clinical info')) 
            trilevel_regressors(sub,:) = NaN;
            trilevel_regressors(sub,1) = PID(sub,1);
            curr_dat_loss.dat(:,sub) = [];
            curr_dat_gain.dat(:,sub) = [];
        end
    end

    GenDis = trilevel_regressors(:,2);
    Anhedonia = trilevel_regressors(:,3);
    Fears = trilevel_regressors(:,4);

    GenDis(any(isnan(GenDis), 2), :) = [];
    Anhedonia(any(isnan(Anhedonia), 2), :) = [];
    Fears(any(isnan(Fears), 2), :) = [];
  
    curr_dat_gain.X = Anhedonia;
    curr_dat_loss.X = Anhedonia;
    
    temp_results_gain = regress(curr_dat_gain,'robust');
    temp_results_loss = regress(curr_dat_loss,'robust');
    
    results_struct.gain.(curr_analysis_dir{roi}) = threshold(temp_results_gain.t,.001,'unc','k',30);
    results_struct.loss.(curr_analysis_dir{roi}) = threshold(temp_results_loss.t,.001,'unc','k',30);
    
    
end

%% quick visualization loop: TRILEVEL
analysis_dir = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD/func_conn/PPI/results/Trilevel_analysis/consumption';
analysis_type = 'Robust'; % Also have OLS
outcome = 'Anhedonia';

curr_analysis_dir = {'ROFC_to_wholebrain','RHO_Accumbens_to_wholebrain','LOFC_to_wholebrain','LOFC2_to_wholebrain','LHO_Accumbens_to_wholebrain'};

load_fname = filenames(fullfile(analysis_dir,strcat(analysis_type,'*',outcome,'*')));
load(load_fname{1})

for roi = 1:length(curr_analysis_dir)
    curr_analysis_dir{roi}
    cluster_out_gain = orthviews(results_struct.gain.(curr_analysis_dir{roi}))
    
    cluster_out_loss = orthviews(results_struct.loss.(curr_analysis_dir{roi}))
    fname = strcat(curr_analysis_dir{roi},'_Anhedonia.mat');
    save(fullfile(analysis_dir,fname),'cluster_out_gain','cluster_out_loss')
    data(1,:) = [cluster_out_gain{1}(1).mm_center,cluster_out_gain{1}(1).numVox,mean(cluster_out_gain{1}(1).Z(:))]
    data(2,:) = [cluster_out_loss{1}(1).mm_center,cluster_out_loss{1}(1).numVox,mean(cluster_out_loss{1}(1).Z(:))]
    keyboard
end

%% OLD
clinicaldir = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD/clinical_data';

clear sub
 
load(fullfile(clinicaldir,'trilevel_factors.mat'));
trilevel_array = [trilevel.ID,trilevel.GenDis,trilevel.Anhedon,trilevel.Fears];
for sub = 1:length(loss_fnames)
    PID(sub,1) = str2num(loss_fnames{sub}(104:108));
    if isempty(find(trilevel.ID(:) == PID(sub,1))) == 0
        curr = find(trilevel.ID(:) == PID(sub,1));
        trilevel_regressors(sub,:) = trilevel_array(curr,:);
    else
        % for OFC loss_ppi_fnames{sub}(100:104); for HO_VMPFC loss_ppi_fnames{sub}(105:109)
        disp(strcat(loss_fnames{sub}(104:108), ' missing clinical info')) 
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

%% TEMPORARY
% Missing timing files from Brian's csv so I need to pull Ann's old files
% and make the appropriate changes. For anticipation onsets, I'm
% subtracting 1.9 seconds from each of Ann's onsets and changing durations
% to be 4 seconds. 

% load in txt that has missing subjects
txtdir = '/home/zach/Documents/current_projects/ACNlab/BrainMAPD';
ant_timing_dir = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD/Oldham_ROI_by_diagnosis/TimingFiles_082218/FSL_anticipation_072418/';
savedir = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD/func_conn/final_timing_files/run-1/anticipation/spm_all_vs_0_timing';
% load(fullfile(txtdir, 'subs_left_out.txt'));
% subs_left_out
clear sub
for sub = 1:length(subs_left_out)
    curr_sub = num2str(subs_left_out(sub));    
    trial_type_strings = {'antgain_all','antloss_all','antgain_0','antloss_0','motor_all'};
    
    if isfile(filenames(fullfile(ant_timing_dir,strcat(curr_sub,'*'),strcat(curr_sub,'_Anticipation_Win_Run1.txt')))) == 1
        for event = 1:length(trial_type_strings)
            if strcmp(trial_type_strings{event}, 'antgain_all') == 1
                fname = filenames(fullfile(ant_timing_dir,strcat(curr_sub,'*'),strcat(curr_sub,'_Anticipation_Win_Run1.txt')));
                temp = load(fname{1});
                onsets{event} = temp(:,1) - 1.9;
             
                name{event} = trial_type_strings{event};
                durations{event} = ones(length(onsets{event}),1) .* 4;
            elseif strcmp(trial_type_strings{event}, 'antloss_all') == 1
                fname = filenames(fullfile(ant_timing_dir,strcat(curr_sub,'*'),strcat(curr_sub,'_Anticipation_Loss_Run1.txt')));
                temp = load(fname{1});
                onsets{event} = temp(:,1) - 1.9;
                name{event} = trial_type_strings{event};
                durations{event} = ones(length(onsets{event}),1) .* 4;
            elseif strcmp(trial_type_strings{event}, 'antgain_0') == 1
                fname = filenames(fullfile(ant_timing_dir,strcat(curr_sub,'*'),strcat(curr_sub,'_Anticipation_Win0_Run1.txt')));
                temp = load(fname{1});
                onsets{event} = temp(:,1) - 1.9;
                name{event} = trial_type_strings{event};
                durations{event} = ones(length(onsets{event}),1) .* 4;
            elseif strcmp(trial_type_strings{event}, 'antloss_0') == 1
                fname = filenames(fullfile(ant_timing_dir,strcat(curr_sub,'*'),strcat(curr_sub,'_Anticipation_Loss0_Run1.txt')));
                temp = load(fname{1});
                onsets{event} = temp(:,1) - 1.9;
                name{event} = trial_type_strings{event};
                durations{event} = ones(length(onsets{event}),1) .* 4;
            elseif strcmp(trial_type_strings{event}, 'motor_all') == 1
                fname = filenames(fullfile(ant_timing_dir,strcat(curr_sub,'*'),strcat(curr_sub,'_Motor_Run1.txt')));
                temp = load(fname{1});
                onsets{event} = temp(:,1);
                name{event} = trial_type_strings{event};
                durations{event} = ones(length(onsets{event}),1) .* 4;
            end
        end
    
    temp_file_name = strcat(curr_sub,'_anticipation_timing.mat');
    save(fullfile(savedir, temp_file_name),'onsets','names','durations')
    else
        disp(curr_sub)
    end
    
end

%% temporary 2

% Missing timing files from Brian's csv so I need to pull Ann's old files
% and make the appropriate changes. For anticipation onsets, I'm
% subtracting 1.9 seconds from each of Ann's onsets and changing durations
% to be 4 seconds. 

% load in txt that has missing subjects
txtdir = '/home/zach/Documents/current_projects/ACNlab/BrainMAPD';
ant_timing_dir = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD/Oldham_ROI_by_diagnosis/TimingFiles_082218/FSL_consumption_110818';
savedir = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD/func_conn/final_timing_files/run-1/consumption/spm_all_vs_0_timing';
%load(fullfile(txtdir, 'subs_left_out.txt'));

clear sub
for sub = 1:length(subs_left_out)
    curr_sub = num2str(subs_left_out(sub));    
    trial_type_strings = {'congain_all','conloss_all','congain_0','conloss_0','motor_all'};
    
    if isfile(filenames(fullfile(ant_timing_dir,strcat(curr_sub,'*'),strcat(curr_sub,'_Hit_Win_Feedback_Run1.txt')))) == 1
        for event = 1:length(trial_type_strings)
            if strcmp(trial_type_strings{event}, 'congain_all') == 1
                fname = filenames(fullfile(ant_timing_dir,strcat(curr_sub,'*'),strcat(curr_sub,'_Hit_Win_Feedback_Run1.txt')));
                temp = load(fname{1});
                onsets{event} = temp(:,1);
             
                name{event} = trial_type_strings{event};
                durations{event} = ones(length(onsets{event}),1) .* 4;
            elseif strcmp(trial_type_strings{event}, 'conloss_all') == 1
                fname = filenames(fullfile(ant_timing_dir,strcat(curr_sub,'*'),strcat(curr_sub,'_Miss_Loss_Feedback_Run1.txt')));
                temp = load(fname{1});
                onsets{event} = temp(:,1);
                name{event} = trial_type_strings{event};
                durations{event} = ones(length(onsets{event}),1) .* 4;
            elseif strcmp(trial_type_strings{event}, 'congain_0') == 1
                fname = filenames(fullfile(ant_timing_dir,strcat(curr_sub,'*'),strcat(curr_sub,'_Hit_Win_0_Feedback_Run1.txt')));
                temp = load(fname{1});
                onsets{event} = temp(:,1);
                name{event} = trial_type_strings{event};
                durations{event} = ones(length(onsets{event}),1) .* 4;
            elseif strcmp(trial_type_strings{event}, 'conloss_0') == 1
                fname = filenames(fullfile(ant_timing_dir,strcat(curr_sub,'*'),strcat(curr_sub,'_Miss_Loss_0_Feedback_Run1.txt')));
                temp = load(fname{1});
                onsets{event} = temp(:,1);
                name{event} = trial_type_strings{event};
                durations{event} = ones(length(onsets{event}),1) .* 4;
            elseif strcmp(trial_type_strings{event}, 'motor_all') == 1
                fname = filenames(fullfile(ant_timing_dir,strcat(curr_sub,'*'),strcat(curr_sub,'_Motor_Run1.txt')));
                temp = load(fname{1});
                onsets{event} = temp(:,1);
                name{event} = trial_type_strings{event};
                durations{event} = ones(length(onsets{event}),1) .* 4;
            end
        end
    
    temp_file_name = strcat(curr_sub,'_consumption_timing.mat');
    save(fullfile(savedir, temp_file_name),'onsets','names','durations')
    else
        disp(curr_sub)
    end
    
end

