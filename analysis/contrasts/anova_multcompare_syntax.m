% addpaths
contrastdir = '~/Documents/repo/acnlab_repo/data';
addpath(genpath('~/Documents/repo'))
%% load the table in thar
% curr_analysis_table = readtable(fullfile(contrastdir,'data.txt'));
load('mid_significant_contrasts.mat')

%% General set up
lVS_loss = curr_analysis_table.lVS_Oldham_Loss_Loss;
lVS_gain = curr_analysis_table.lVS_Oldham_Loss_Gain;

dsm_diagnoses_regressors = [curr_analysis_table.Dep(:),curr_analysis_table.Anx(:),curr_analysis_table.Comorbid(:)];

anova_regressors = ones(height(curr_analysis_table),1);
anova_regressors(dsm_diagnoses_regressors(:,1)==1) = 2;
anova_regressors(dsm_diagnoses_regressors(:,2)==1) = 3;
anova_regressors(dsm_diagnoses_regressors(:,3)==1) = 4;

anova_regressors_strings = cell(size(anova_regressors));
anova_regressors_strings(anova_regressors(:,1)==1) = {'Healthy'};
anova_regressors_strings(anova_regressors(:,1)==2) = {'Depression'};
anova_regressors_strings(anova_regressors(:,1)==3) = {'Anxiety'};
anova_regressors_strings(anova_regressors(:,1)==4) = {'Comorbidity'};

sex = curr_analysis_table.sex;
psych = curr_analysis_table.PsychAny;
dop = curr_analysis_table.DopAny;
anhed = trilevel.Anhedon;

%% lVS model

[p,tbl,stats,terms] = anovan(lVS_loss,{anova_regressors_strings,anhed,lVS_gain,sex,psych},'varnames',{'Diagnosis','lVS_gain_anticipation','Sex','Psychotropic_meds'},'continuous',[3])

[c,m,h]=multcompare(stats)

%% for Oldham rOFC now
rOFC_loss = curr_analysis_table.Oldham_rOFC_loss;
rOFC_gain = curr_analysis_table.Oldham_rOFC_gain;

[p,tbl,stats,terms] = anovan(rOFC_gain,{anova_regressors_strings,rOFC_loss,sex,psych},'varnames',{'Diagnosis','rOFC loss receipt','Sex','Psychotropic_meds'},'continuous',[2])

[c,m,h]=multcompare(stats)
%% for Ng OFC bilateral
 
Ng_OFC_loss = curr_analysis_table.OFC_Ng_ConLoss_v_ConNoLoss_avg;
Ng_OFC_gain = curr_analysis_table.OFC_Ng_ConGain_v_ConNoGain_avg;
Ng_OFC_AntGain = curr_analysis_table.OFC_Ng_AntGain_v_AntNoGain_avg;
Ng_OFC_AntLoss = curr_analysis_table.OFC_Ng_AntLoss_v_AntNoLoss_avg;

% Gain controlling for Loss NOT SIG FOR BILATERAL
[p,tbl,stats,terms] = anovan(Ng_OFC_gain,{anova_regressors_strings,Ng_OFC_loss,sex,psych},'varnames',{'Diagnosis','OFC loss outcome','Sex','Psychotropic_meds'},'continuous',[2])

% AntGain controlling for AntLoss NOT SIG FOR BILATERAL
[p,tbl,stats,terms] = anovan(Ng_OFC_AntGain,{anova_regressors_strings,Ng_OFC_AntLoss,sex,psych},'varnames',{'Diagnosis','OFC loss antic','Sex','Psychotropic_meds'},'continuous',[2])

%% for Ng OFC right
 
R_Ng_OFC_loss = curr_analysis_table.R_OFC_Ng_ConLoss_v_ConNoLoss_avg;
R_Ng_OFC_gain = curr_analysis_table.R_OFC_Ng_ConGain_v_ConNoGain_avg;
R_Ng_OFC_AntGain = curr_analysis_table.R_OFC_Ng_AntGain_v_AntNoGain_avg;
R_Ng_OFC_AntLoss = curr_analysis_table.R_OFC_Ng_AntLoss_v_AntNoLoss_avg;

% Gain controlling for Loss NOT SIG FOR BILATERAL
[p,tbl,stats,terms] = anovan(R_Ng_OFC_gain,{anova_regressors_strings,R_Ng_OFC_loss,sex,psych},'varnames',{'Diagnosis','OFC loss outcome','Sex','Psychotropic_meds'},'continuous',[2])

% AntGain controlling for AntLoss NOT SIG FOR BILATERAL
[p,tbl,stats,terms] = anovan(R_Ng_OFC_AntGain,{anova_regressors_strings,R_Ng_OFC_AntLoss,sex,psych},'varnames',{'Diagnosis','OFC loss antic','Sex','Psychotropic_meds'},'continuous',[2])

%% for Ng OFC left

L_Ng_OFC_loss = curr_analysis_table.L_OFC_Ng_ConLoss_v_ConNoLoss_avg;
L_Ng_OFC_gain = curr_analysis_table.L_OFC_Ng_ConGain_v_ConNoGain_avg;
L_Ng_OFC_AntGain = curr_analysis_table.L_OFC_Ng_AntGain_v_AntNoGain_avg;
L_Ng_OFC_AntLoss = curr_analysis_table.L_OFC_Ng_AntLoss_v_AntNoLoss_avg;

% Gain controlling for Loss NOT SIG FOR BILATERAL
[p,tbl,stats,terms] = anovan(L_Ng_OFC_gain,{anova_regressors_strings,L_Ng_OFC_loss,sex,psych},'varnames',{'Diagnosis','OFC loss outcome','Sex','Psychotropic_meds'},'continuous',[2])

% AntGain controlling for AntLoss NOT SIG FOR BILATERAL
[p,tbl,stats,terms] = anovan(L_Ng_OFC_AntGain,{anova_regressors_strings,L_Ng_OFC_AntLoss,sex,psych},'varnames',{'Diagnosis','OFC loss antic','Sex','Psychotropic_meds'},'continuous',[2])

%% for Sphere VS bilateral
 
VS_loss = curr_analysis_table.VS_Sphere_ConLoss_v_ConNoLoss_avg;
VS_gain = curr_analysis_table.VS_Sphere_ConGain_v_ConNoGain_avg;
VS_AntGain = curr_analysis_table.VS_Sphere_AntGain_v_AntNoGain_avg;
VS_AntLoss = curr_analysis_table.VS_Sphere_AntLoss_v_AntNoLoss_avg;

% Gain controlling for Loss NOT SIG FOR BILATERAL
[p,tbl,stats,terms] = anovan(VS_gain,{anova_regressors_strings,VS_loss,sex,psych},'varnames',{'Diagnosis','VS loss outcome','Sex','Psychotropic_meds'},'continuous',[2])

% AntGain controlling for AntLoss NOT SIG FOR BILATERAL
[p,tbl,stats,terms] = anovan(VS_AntGain,{anova_regressors_strings,VS_AntLoss,sex,psych},'varnames',{'Diagnosis','VS loss antic','Sex','Psychotropic_meds'},'continuous',[2])

%% for Sphere VS right
 
R_VS_loss = curr_analysis_table.R_VS_Sphere_ConLoss_v_ConNoLoss_avg;
R_VS_gain = curr_analysis_table.R_VS_Sphere_ConGain_v_ConNoGain_avg;
R_VS_AntGain = curr_analysis_table.R_VS_Sphere_AntGain_v_AntNoGain_avg;
R_VS_AntLoss = curr_analysis_table.R_VS_Sphere_AntLoss_v_AntNoLoss_avg;

% Gain controlling for Loss 
[p,tbl,stats,terms] = anovan(R_VS_gain,{anova_regressors_strings,R_VS_loss,sex,psych},'varnames',{'Diagnosis','R VS loss outcome','Sex','Psychotropic_meds'},'continuous',[2])

% AntGain controlling for AntLoss 
[p,tbl,stats,terms] = anovan(R_VS_AntGain,{anova_regressors_strings,R_VS_AntLoss,sex,psych},'varnames',{'Diagnosis','R VS loss antic','Sex','Psychotropic_meds'},'continuous',[2])

%% for Sphere VS left

L_VS_loss = curr_analysis_table.L_VS_Sphere_ConLoss_v_ConNoLoss_avg;
L_VS_gain = curr_analysis_table.L_VS_Sphere_ConGain_v_ConNoGain_avg;
L_VS_AntGain = curr_analysis_table.L_VS_Sphere_AntGain_v_AntNoGain_avg;
L_VS_AntLoss = curr_analysis_table.L_VS_Sphere_AntLoss_v_AntNoLoss_avg;

% Gain controlling for Loss 
[p,tbl,stats,terms] = anovan(L_VS_gain,{anova_regressors_strings,L_VS_loss,sex,psych},'varnames',{'Diagnosis','L VS loss outcome','Sex','Psychotropic_meds'},'continuous',[2])

% AntGain controlling for AntLoss
[p,tbl,stats,terms] = anovan(L_VS_AntGain,{anova_regressors_strings,L_VS_AntLoss,sex,psych},'varnames',{'Diagnosis','L VS loss antic','Sex','Psychotropic_meds'},'continuous',[2])

