% 4 groups paper code
% You've been sitting on this and you know it bud. Pull your head out of
% your ass and let's get at er
maskdir = '/Users/zaz3744/Documents/current_projects/ACNlab/masks/ROI_BrainMAPD';
datadir = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD/func_conn/first_levels/first_level_output/consumption';
figdir = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD/Oldham_ROI_by_diagnosis/second_level_output';
motiondir = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD/Oldham_ROI_by_diagnosis/motion';
clinicaldir = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD/clinical_data';
savedir = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD/Oldham_ROI_by_diagnosis/second_level_output';

%% Load in data of all kinds
% gotta make sure my indices all line up

motion_fnames = filenames(fullfile(motiondir, '*run1.mat'));
con1_fnames = filenames(fullfile(datadir,'*/ses-2/run-1/MID/con*1.nii'));
for sub = 1:length(con1_fnames)
    curr_id = con1_fnames{sub}(118:122);
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

%% Remove motion outliers based on FD

fnames_gain_consumption_temp = filenames(fullfile(datadir, '*/ses-2/run-1/MID/con_0002.nii'));
fnames_loss_consumption_temp = filenames(fullfile(datadir, '*/ses-2/run-1/MID/con_0001.nii'));


fnames_gain_vs_no_gain = fnames_gain_consumption_temp(subj_motion<.3);
fnames_loss_vs_no_loss = fnames_loss_consumption_temp(subj_motion<.3);

R = R(subj_motion<.3,:);

%% BEGIN WHOLE BRAIN ANALYSIS SECTION
% Current Gain Contrast
data_gain_con = fmri_data(fnames_gain_vs_no_gain)
data_gain_con.X = R;
out_gain_consump = regress(data_gain_con, .001, 'unc')
out_gain_consump_coord = image2coordinates(out_gain_consump.t);

out_gain_consump_coord_anx = tal2mni(out_gain_consump_coord{1});
out_gain_consump_coord_dep = tal2mni(out_gain_consump_coord{2});
out_gain_consump_coord_com = tal2mni(out_gain_consump_coord{3});
out_gain_consump_coord_avg = tal2mni(out_gain_consump_coord{4});

%% Current Loss Contrast
data_loss_con = fmri_data(fnames_loss_vs_no_loss)
data_loss_con.X = R;
out_loss_consump = regress(data_loss_con, .001, 'unc')
out_loss_consump_coord = image2coordinates(out_loss_consump.t);

out_loss_consump_coord_anx = tal2mni(out_loss_consump_coord{1});
out_loss_consump_coord_dep = tal2mni(out_loss_consump_coord{2});
out_loss_consump_coord_com = tal2mni(out_loss_consump_coord{3});
out_loss_consump_coord_avg = tal2mni(out_loss_consump_coord{4});

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

% Gain consumption
bi_ofc_gain_con = data_gain_con.dat .* right_ofc_consumption_mask.dat;
bi_ofc_gain_con(bi_ofc_gain_con==0)=NaN;
bi_ofc_gain_con_bold_avg = nanmean(bi_ofc_gain_con);
% Loss consumption
bi_ofc_loss_con = data_loss_con.dat .* right_ofc_consumption_mask.dat;
bi_ofc_loss_con(bi_ofc_loss_con==0)=NaN;
bi_ofc_loss_con_bold_avg = nanmean(bi_ofc_loss_con);


%% Regression analysis
% Get regressors
anova_diagnoses = sum([R(:,1), (R(:,2)*2), (R(:,3)*3)],2);  

anova_regressors_strings = cell(size(anova_diagnoses));
anova_regressors_strings(anova_diagnoses(:,1)==0) = {'Healthy'};
anova_regressors_strings(anova_diagnoses(:,1)==1) = {'Depression'};
anova_regressors_strings(anova_diagnoses(:,1)==2) = {'Anxiety'};
anova_regressors_strings(anova_diagnoses(:,1)==3) = {'Comorbidity'};

[biOFC_gaincon_p,biOFC_gaincon_tbl,biOFC_gaincon_stats] = anova1(bi_ofc_gain_con_bold_avg(:),anova_regressors_strings);
[biOFC_losscon_p,biOFC_losscon_tbl,biOFC_losscon_stats] = anova1(bi_ofc_loss_con_bold_avg(:),anova_regressors_strings);

%% OLD
% %% Replace 0's with scaled contrast numbers so SPM doesn't freak out.
% % I think what's happening is that spm is not doing a traditional
% % regression through the gui. So inputing a matrix of 1's and 0's doesn't
% % cut it. Total sum of vector must = 0, so I need to replace 
% 
% conval_anx = sum(life_anx)/sum(life_anx==0);
% conval_dep = sum(life_dep)/sum(life_dep==0);
% conval_com = sum(life_com)/sum(life_com==0);
% 
% 
% life_anx(find(life_anx==0)) = conval_anx * -1;
% life_dep(find(life_dep==0)) = conval_dep * -1;
% life_com(find(life_com==0)) = conval_com * -1;
% 
% final = [life_anx';life_dep';life_com'];
% 
% save(fullfile(savedir,'temp_dx_regressors.mat'),'final')

