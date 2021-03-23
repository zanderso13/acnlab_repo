% 4 groups paper code
% You've been sitting on this and you know it bud. Pull your head out of
% your ass and let's get at er
maskdir = '/Users/zaz3744/Documents/current_projects/ACNlab/masks/ROI_BrainMAPD_functional';
datadir = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD/func_conn/first_levels/first_level_output/consumption';
figdir = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD/Oldham_ROI_by_diagnosis/second_level_output';
motiondir = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD/Oldham_ROI_by_diagnosis/motion';
clinicaldir = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD/clinical_data';
savedir = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD/Oldham_ROI_by_diagnosis/second_level_output';

% Define what mask you want to use. Here are the current options:
% 1 - bi_ofc_consumption_mask.dat
% 2 - left_ofc_consumption_mask.dat
% 3 - right_ofc_consumption_mask.dat
% 4 - bi_vs_consumption_mask.dat
% 5 - left_vs_consumption_mask.dat
% 6 - right_vs_consumption_mask.dat
% 7 - bi_vs_rew_anticipation_mask.dat
% 8 - left_vs_rew_anticipation_mask.dat
% 9 - right_vs_rew_anticipation_mask.dat
% 10 - bi_vs_loss_anticipation_mask.dat
% 11 - left_vs_loss_anticipation_mask.dat
% 12 - right_vs_loss_anticipation_mask.dat

% I can add more as we keep going. Should cause these analyses to go more
% rapidly than before. Or maybe they went quickly before. That's unclear to
% me. 

mask = 4;

% What motion cutoff do you want to use?

motion_cutoff = 0.3;
%% Load in data of all kinds
% gotta make sure my indices all line up

motion_fnames = filenames(fullfile(motiondir, '*run1.mat'));
con1_fnames = filenames(fullfile(datadir,'*/ses-2/run-1/MID/con*1.nii'));

for i = 1 : length(con1_fnames)
    sub_ids{i} = con1_fnames{i}(118:122);
end

for sub = 1:length(con1_fnames)
    curr_id = sub_ids{sub};
    curr_file = contains(motion_fnames,curr_id);
    load(motion_fnames{curr_file});
    sub_id{sub,1} = motion_fnames{sub}(94:98);
    subj_motion(sub,1) = mean(FD.framewise_displacement);
end

figure(); hist(subj_motion); title('Framewise displacement')


%% So now I need regressors for the group level contrasts
% Each row will be a contrast corresponding to one of the four groups. If I
% code this then I'll be able to add/remove subjects later.
% I'm saving the final version of this so only run it when you get more
% subjects. Manually going in and finding diagnoses for the following subs

% 10176 - healthy control
% 20461 - Not in most recent version
% 20674 - Not in most recent version

clear sub
 

load(fullfile(clinicaldir,'BrainMAPD_clinical_diagnoses_final.mat'))
for sub = 1:length(sub_id)
    PID(sub,1) = str2num(sub_id{sub});
    if isempty(find(clinical_info.PID(:) == str2num(sub_id{sub}))) == 0
        curr = find(clinical_info.PID(:) == str2num(sub_id{sub}));
        life_dep(sub,1) = clinical_info.dep_life_any(curr);
        life_anx(sub,1) = clinical_info.anx_life_any(curr);
        life_com(sub,1) = clinical_info.comorbid_life_dep_anx(curr);
    else
        disp(strcat(sub_id{sub}, ' missing clinical info'))
        life_dep(sub,1) = 0;%NaN;
        life_anx(sub,1) = 0;%NaN;
        life_com(sub,1) = 0;%NaN;
    end
end

R = [life_anx,life_dep,life_com];
life_healthy = zeros(length(R),1);
life_healthy(find(sum(R,2)==0)) = 1;
R = [R,ones(length(life_dep),1)];
%second_level_regressors = array2table(second_level_regressors); second_level_regressors.Properties.VariableNames = {'anx','dep','com','healthy'};
save(fullfile(savedir,'temp_second_level_regressors.mat'),'R')

%% Extract Trilevel information
% Gotta see if we replicate the Anhedonia findings Rick and Robin are
% publishing

clear sub
 
load(fullfile(clinicaldir,'trilevel_factors.mat'));
trilevel_array = [trilevel.ID,trilevel.GenDis,trilevel.Anhedon,trilevel.Fears];

for sub = 1:length(con1_fnames)
    curr_id = sub_ids{sub};
    PID(sub,1) = str2num(curr_id);
    if isempty(find(trilevel.ID(:) == str2num(curr_id))) == 0
        curr = find(trilevel.ID(:) == str2num(curr_id));
        curr_analysis_table(sub,:) = trilevel_array(curr,:);
    else
        disp(strcat(curr_id, ' missing clinical info')) 
        curr_analysis_table(sub,:) = NaN;
        curr_analysis_table(sub,1) = str2double(curr_id);
    end
end

GenDis = trilevel_array(:,2);
Anhedonia = trilevel_array(:,3);
Fears = trilevel_array(:,4);

%% Remove motion outliers based on FD

fnames_gain_consumption_temp = filenames(fullfile(datadir, '*/ses-2/run-1/MID/con_0002.nii'));
fnames_loss_consumption_temp = filenames(fullfile(datadir, '*/ses-2/run-1/MID/con_0001.nii'));


fnames_gain_vs_no_gain = fnames_gain_consumption_temp(subj_motion<motion_cutoff);
fnames_loss_vs_no_loss = fnames_loss_consumption_temp(subj_motion<motion_cutoff);

R = R(subj_motion<motion_cutoff,:);
trilevel_array = trilevel_array(subj_motion<motion_cutoff,:);

GenDis = GenDis(subj_motion<motion_cutoff,:);
Anhedonia = Anhedonia(subj_motion<motion_cutoff,:);
Fears = Fears(subj_motion<motion_cutoff,:);


%% BEGIN WHOLE BRAIN ANALYSIS SECTION FOR DSM DIAGNOSIS
% Current Gain Contrast
data_gain_con = fmri_data(fnames_gain_vs_no_gain)
data_gain_con.X = R;
out_gain_consump = regress(data_gain_con, .05, 'fdr')
out_gain_consump_coord = image2coordinates(out_gain_consump.t);

% out_gain_consump_coord_anx = tal2mni(out_gain_consump_coord{1});
% out_gain_consump_coord_dep = tal2mni(out_gain_consump_coord{2});
% out_gain_consump_coord_com = tal2mni(out_gain_consump_coord{3});
% out_gain_consump_coord_avg = tal2mni(out_gain_consump_coord{4});

%% Current Loss Contrast
data_loss_con = fmri_data(fnames_loss_vs_no_loss)
data_loss_con.X = R;
out_loss_consump = regress(data_loss_con, .05, 'fdr')
out_loss_consump_coord = image2coordinates(out_loss_consump.t);

% out_loss_consump_coord_anx = tal2mni(out_loss_consump_coord{1});
% out_loss_consump_coord_dep = tal2mni(out_loss_consump_coord{2});
% out_loss_consump_coord_com = tal2mni(out_loss_consump_coord{3});
% out_loss_consump_coord_avg = tal2mni(out_loss_consump_coord{4});

%% ROI ANALYSIS
% load in masks of interest
bi_ofc_consumption_mask = fmri_data(fullfile(maskdir,'consumption','OFC_8mmsphere_Oldham.nii'));
left_ofc_consumption_mask = fmri_data(fullfile(maskdir,'consumption/left','L_OFC_8mmsphere_Oldham.nii'));
left_ofc_consumption_mask2 = fmri_data(fullfile(maskdir,'consumption/left','L_OFC_8mmsphere_Oldham2.nii'));
right_ofc_consumption_mask = fmri_data(fullfile(maskdir,'consumption/right','R_OFC_8mmsphere_Oldham.nii'));

bi_vs_consumption_mask = fmri_data(fullfile(maskdir,'consumption','VS_8mmsphere_Oldham_Con.nii'));
left_vs_consumption_mask = fmri_data(fullfile(maskdir,'consumption/left','L_VS_8mmsphere_Oldham_Con.nii'));
right_vs_consumption_mask = fmri_data(fullfile(maskdir,'consumption/right','R_VS_8mmsphere_Oldham_Con.nii'));

bi_vs_rew_anticipation_mask = fmri_data(fullfile(maskdir,'anticipation','VS_8mmsphere_Oldham_Rew.nii'));
left_vs_rew_anticipation_mask = fmri_data(fullfile(maskdir,'anticipation/left','L_VS_8mmsphere_Oldham_Rew.nii'));
right_vs_rew_anticipation_mask = fmri_data(fullfile(maskdir,'anticipation/right','R_VS_8mmsphere_Oldham_Rew.nii'));

bi_vs_loss_anticipation_mask = fmri_data(fullfile(maskdir,'anticipation','VS_8mmsphere_Oldham_Loss.nii'));
left_vs_loss_anticipation_mask = fmri_data(fullfile(maskdir,'anticipation/left','L_VS_8mmsphere_Oldham_Loss.nii'));
right_vs_loss_anticipation_mask = fmri_data(fullfile(maskdir,'anticipation/right','R_VS_8mmsphere_Oldham_Loss.nii'));

%% Extract activation from contrasts OFC

mask_stash = {bi_ofc_consumption_mask.dat,left_ofc_consumption_mask.dat,right_ofc_consumption_mask.dat,...
    bi_vs_consumption_mask.dat,left_vs_consumption_mask.dat,right_vs_consumption_mask.dat,...
    bi_vs_rew_anticipation_mask.dat,left_vs_rew_anticipation_mask.dat,right_vs_rew_anticipation_mask.dat,...
    bi_vs_loss_anticipation_mask.dat,left_vs_loss_anticipation_mask.dat,right_vs_loss_anticipation_mask.dat};

curr_mask = mask_stash{mask};
% Gain consumption
masked_gain_con = data_gain_con.dat .* curr_mask;
masked_gain_con(masked_gain_con==0)=NaN;
masked_gain_con_bold_avg = nanmean(masked_gain_con);
% Loss consumption
masked_loss_con = data_loss_con.dat .* curr_mask;
masked_loss_con(masked_loss_con==0)=NaN;
masked_loss_con_bold_avg = nanmean(masked_loss_con);


%% Regression analysis
% Get regressors

anova_diagnoses = sum([R(:,1), (R(:,2)*2), (R(:,3)*3)],2);  

anova_regressors_strings = cell(size(anova_diagnoses));
anova_regressors_strings(anova_diagnoses(:,1)==0) = {'Healthy'};
anova_regressors_strings(anova_diagnoses(:,1)==1) = {'Depression'};
anova_regressors_strings(anova_diagnoses(:,1)==2) = {'Anxiety'};
anova_regressors_strings(anova_diagnoses(:,1)==3) = {'Comorbidity'};


[gain_p,gain_tbl,gain_stats] = anova1(masked_gain_con_bold_avg(:)',anova_regressors_strings);
[loss_p,loss_tbl,loss_stats] = anova1(masked_loss_con_bold_avg(:)',anova_regressors_strings);


%% BEGIN WHOLE BRAIN ANALYSIS SECTION FOR TRILEVEL 
% Current Gain Contrast and GenDis
data_gain_con = fmri_data(fnames_gain_vs_no_gain)
data_gain_con.X = GenDis;
out_gain_consump = regress(data_gain_con, .01, 'unc')
out_gain_consump_coord = image2coordinates(out_gain_consump.t);

% out_gain_consump_coord_anx = tal2mni(out_gain_consump_coord{1});
% out_gain_consump_coord_dep = tal2mni(out_gain_consump_coord{2});
% out_gain_consump_coord_com = tal2mni(out_gain_consump_coord{3});
% out_gain_consump_coord_avg = tal2mni(out_gain_consump_coord{4});

%% Current Gain Contrast and Anhedonia
data_gain_con = fmri_data(fnames_gain_vs_no_gain)
data_gain_con.X = Anhedonia;
out_gain_consump = regress(data_gain_con, .01, 'unc')
out_gain_consump_coord = image2coordinates(out_gain_consump.t);

%% Current Gain Contrast and Fears
data_gain_con = fmri_data(fnames_gain_vs_no_gain)
data_gain_con.X = Fears;
out_gain_consump = regress(data_gain_con, .01, 'unc')
out_gain_consump_coord = image2coordinates(out_gain_consump.t);


%% Current Loss Contrast and GenDis
data_loss_con = fmri_data(fnames_loss_vs_no_loss)
data_loss_con.X = GenDis;
out_loss_consump = regress(data_loss_con, .01, 'unc')
out_loss_consump_coord = image2coordinates(out_loss_consump.t);

% out_loss_consump_coord_anx = tal2mni(out_loss_consump_coord{1});
% out_loss_consump_coord_dep = tal2mni(out_loss_consump_coord{2});
% out_loss_consump_coord_com = tal2mni(out_loss_consump_coord{3});
% out_loss_consump_coord_avg = tal2mni(out_loss_consump_coord{4});

%% Current Loss Contrast and Anhedonia
data_loss_con = fmri_data(fnames_loss_vs_no_loss)
data_loss_con.X = Anhedonia;
out_loss_consump = regress(data_loss_con, .01, 'unc')
out_loss_consump_coord = image2coordinates(out_loss_consump.t);

%% Current Loss Contrast and Fears
data_loss_con = fmri_data(fnames_loss_vs_no_loss)
data_loss_con.X = Fears;
out_loss_consump = regress(data_loss_con, .01, 'unc')
out_loss_consump_coord = image2coordinates(out_loss_consump.t);

%% ROI analysis trilevel: GenDis

curr_analysis_table_gain = [masked_gain_con_bold_avg',trilevel_array(:,2:size(trilevel_array,2))];
curr_analysis_table_loss = [masked_loss_con_bold_avg',trilevel_array(:,2:size(trilevel_array,2))];

curr_analysis_table_gain = array2table(curr_analysis_table_gain); curr_analysis_table_gain.Properties.VariableNames = {'Brain','GenDis','Anhedonia','Fears'};
curr_analysis_table_loss = array2table(curr_analysis_table_loss); curr_analysis_table_loss.Properties.VariableNames = {'Brain','GenDis','Anhedonia','Fears'};

mdl_gain = fitlm(curr_analysis_table_gain,'Brain ~ GenDis')
anova(mdl_gain,'summary')

mdl_gain = fitlm(curr_analysis_table_loss,'Brain ~ GenDis')
anova(mdl_gain,'summary')

%% ROI analysis trilevel: Anhedonia

curr_analysis_table_gain = [masked_gain_con_bold_avg',trilevel_array(:,2:size(trilevel_array,2))];
curr_analysis_table_loss = [masked_loss_con_bold_avg',trilevel_array(:,2:size(trilevel_array,2))];

curr_analysis_table_gain = array2table(curr_analysis_table_gain); curr_analysis_table_gain.Properties.VariableNames = {'Brain','GenDis','Anhedonia','Fears'};
curr_analysis_table_loss = array2table(curr_analysis_table_loss); curr_analysis_table_loss.Properties.VariableNames = {'Brain','GenDis','Anhedonia','Fears'};

mdl_gain = fitlm(curr_analysis_table_gain,'Brain ~ Anhedonia')
anova(mdl_gain,'summary')

mdl_gain = fitlm(curr_analysis_table_loss,'Brain ~ Anhedonia')
anova(mdl_gain,'summary')

%% ROI analysis trilevel: Fears

curr_analysis_table_gain = [masked_gain_con_bold_avg',trilevel_array(:,2:size(trilevel_array,2))];
curr_analysis_table_loss = [masked_loss_con_bold_avg',trilevel_array(:,2:size(trilevel_array,2))];

curr_analysis_table_gain = array2table(curr_analysis_table_gain); curr_analysis_table_gain.Properties.VariableNames = {'Brain','GenDis','Anhedonia','Fears'};
curr_analysis_table_loss = array2table(curr_analysis_table_loss); curr_analysis_table_loss.Properties.VariableNames = {'Brain','GenDis','Anhedonia','Fears'};

mdl_gain = fitlm(curr_analysis_table_gain,'Brain ~ Fears')
anova(mdl_gain,'summary')

mdl_gain = fitlm(curr_analysis_table_loss,'Brain ~ Fears')
anova(mdl_gain,'summary')

%% Begin HA for Psychometrics assignment
% extract voxels
ica_dir = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD/func_conn/ICA';
ica_spatial_mask1 = fmri_data(fullfile(ica_dir,'gica_cmd_mean_component_ica_s_all_.nii'));
ica_spatial_mask1.dat = ica_spatial_mask1.dat(:,1);
ica_spatial_mask1.dat(ica_spatial_mask1.dat>-2 & ica_spatial_mask1.dat<2) = 0;
ica_spatial_mask1.dat(ica_spatial_mask1.dat~=0) = 1;

ica_spatial_mask4 = fmri_data(fullfile(ica_dir,'gica_cmd_mean_component_ica_s_all_.nii'));
ica_spatial_mask4.dat = ica_spatial_mask4.dat(:,21);
ica_spatial_mask4.dat(ica_spatial_mask4.dat>-2 & ica_spatial_mask4.dat<2) = 0;
ica_spatial_mask4.dat(ica_spatial_mask4.dat~=0) = 1;

comp_ica_mask = ica_spatial_mask1;
comp_ica_mask.dat = comp_ica_mask.dat + ica_spatial_mask4.dat; 
comp_ica_mask.dat(comp_ica_mask.dat>0) = 1;


ha_input{1} = data_gain_con.dat .* comp_ica_mask.dat;
ha_input{2} = data_gain_con.dat .* comp_ica_mask.dat;
ha_input{1} = ha_input{1}';
ha_input{2} = ha_input{2}';

[aligned,transforms] = hyperalign(ha_input{1},ha_input{2});

%% Extract most significant voxels
% That means snagging the max t scores, positive or negative which the loop
% does. Then I'm save the output I guess

for sub = 1:size(aligned{1},1)
    [corr_input_gain_aligned(sub,:),corr_input_gain_ind_aligned(sub,:)] = maxk(abs(aligned{1}(sub,:)),1000);
    [corr_input_loss_aligned(sub,:),corr_input_loss_ind_aligned(sub,:)] = maxk(abs(aligned{2}(sub,:)),1000);
    [corr_input_gain_unaligned(sub,:),corr_input_gain_ind_unaligned(sub,:)] = maxk(abs(ha_input{1}(sub,:)),1000);
    [corr_input_loss_unaligned(sub,:),corr_input_loss_ind_unaligned(sub,:)] = maxk(abs(ha_input{2}(sub,:)),1000);
    corr_matrix_gain_aligned(sub,:) = corr_input_gain_aligned(sub,:);
    corr_matrix_loss_aligned(sub,:) = corr_input_loss_aligned(sub,:);
    corr_matrix_gain_unaligned(sub,:) = corr_input_gain_unaligned(sub,:);
    corr_matrix_loss_unaligned(sub,:) = corr_input_loss_unaligned(sub,:);
end

corr_matrix_gain_aligned = corrcoef(corr_matrix_gain_aligned');
corr_matrix_loss_aligned = corrcoef(corr_matrix_loss_aligned');

corr_matrix_gain_unaligned = corrcoef(corr_matrix_gain_unaligned');
corr_matrix_loss_unaligned = corrcoef(corr_matrix_loss_unaligned');

corr_matrix_gain_aligned_to_save = corr_matrix_gain_aligned';
corr_matrix_loss_aligned_to_save = corr_matrix_loss_aligned';
corr_matrix_gain_unaligned_to_save = corr_matrix_gain_unaligned';
corr_matrix_loss_unaligned_to_save = corr_matrix_loss_unaligned';

