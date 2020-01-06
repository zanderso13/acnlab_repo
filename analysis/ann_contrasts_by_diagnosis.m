% For SOBP poster I need to break down initial results and then describe in
% what ways I will extend those findings before the conference. Luckily Ann
% is amazing and has done a lot of this for me. So I just need to select
% the correct contrasts and perform some t-tests to show blunted reward
% response in depression and elevated response in anxiety. So I'm really
% just going to load in Ann's work and then explore the relationships in
% the Win vs neutral condition with depression, anxiety, comorbidity, and
% then probably the tri-level features as well.

% First to set up my directories of course

repodir = '/home/zaz3744/repo';
spmdir = '/home/zaz3744/spm12';
contrastdir = '/projects/p30954/ann_contrasts';
clinicaldir = '/projects/p30954/clinical_data';
medicationdir = '/projects/p30954/demographic_data';

% I've extracted significant results. If you're going back to dig through
% more of Ann's contrasts set the following variable to 0. Otherwise
% setting it to 1 will run the anova's that were already significant so you
% can pull figures, stats, etc. Finally, setting it to 2 will let you look at VS specifically
% Have fun, Party on Wayne

significant = 1;

% I've cherry picked the significant results. 
load(fullfile(contrastdir, 'mid_significant_contrasts.mat'))

% Looking at Gain vs No Gain
% gain = load(fullfile(contrastdir, 'mid_contrasts_gain_vs_no_gain.mat'))

% Looking at Loss vs No Loss
% loss = load(fullfile(contrastdir, 'mid_contrasts_loss_vs_no_loss.mat'))

% Looking at Ant Gain vs Ant No Gain
% load(fullfile(contrastdir, 'mid_contrasts_ant_gain_vs_ant_no_gain.mat'))

% There are so many VS contrasts I've split them into a different file.
% This is just for Anticipation as well. Here are the variables I put in
% each 
% GainVsNoGain = [MIDROI091119Zach.L_VS_Oldham_Rew_AntGain_v_AntNoGain_avg, MIDROI091119Zach.R_VS_Oldham_Rew_AntGain_v_AntNoGain_avg, MIDROI091119Zach.VS_Oldham_Rew_AntGain_v_AntNoGain_avg, MIDROI091119Zach.L_VS_Oldham_Loss_AntGain_v_AntNoGain_avg, MIDROI091119Zach.R_VS_Oldham_Loss_AntGain_v_AntNoGain_avg, MIDROI091119Zach.VS_Oldham_Loss_AntGain_v_AntNoGain_avg];
% LossVsNoLoss = [MIDROI091119Zach.L_VS_Oldham_Rew_AntLoss_v_AntNoLoss_avg, MIDROI091119Zach.R_VS_Oldham_Rew_AntLoss_v_AntNoLoss_avg, MIDROI091119Zach.VS_Oldham_Rew_AntLoss_v_AntNoLoss_avg, MIDROI091119Zach.L_VS_Oldham_Loss_AntLoss_v_AntNoLoss_avg, MIDROI091119Zach.R_VS_Oldham_Loss_AntLoss_v_AntNoLoss_avg, MIDROI091119Zach.VS_Oldham_Loss_AntLoss_v_AntNoLoss_avg];
% load(fullfile(contrastdir, 'mid_contrasts_VS_Ant.mat'));


% Looking at Ant Loss vs Ant No Gain. For VS, I had to select the outcome
% of the trial as well. So I picked loss. This was arbitrary and I need to
% do the win outcome as well. Actually I really need to compare the two and
% then average them. If they're 
% load(fullfile(contrastdir, 'mid_contrasts_ant_loss_vs_ant_no_loss.mat'))

% drop exclusions briefly
% curr_analysis_table(curr_analysis_table.Exclusions==1,:) = [];
% curr_analysis_table.Properties.VariableNames = {'PID','Dep','Anx','Comorbid','GenDis','Anhed','Fears','rOFC_Oldham_GainVsNone','lOFC_Oldham_GainVsNone','rVS_Oldham_GainVsNone','lVS_Oldham_GainVsNone'};
% analyses are going to be a GLM for DSM diagnoses and trilevel model stuff
% separately

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
% R OFC is showing some significance here, let's plot the relationship of R
% OFC activity in each of the three diagnostic bins

% Based on the receipt of reward. These ROIs are all significant, show
% hyperactivity which is not expected. At least in terms of that meta
% analysis though the tasks were pretty different there so who knows.
%% Copy of analyses to run when only looking at the significant contrasts
% Just trying to make my life easier here
if significant == 1
    [rOFC_p,rOFC_tbl,rOFC_stats] = anova1(curr_analysis_table.Oldham_rOFC_gain(:),anova_regressors_strings);
    [lPutamen_p,lPutamen_tbl,lPutamen_stats] = anova1(curr_analysis_table.HO_lPutamen_gain,anova_regressors_strings);
    [biVS_Loss_Loss_p, biVS_Loss_Loss_tbl, biVS_Loss_Loss_stats] = anova1(curr_analysis_table.Oldham_biVS_Loss_Loss, anova_regressors_strings);
    [lVS_Loss_Loss_p, lVS_Loss_Loss_tbl, lVS_Loss_Loss_stats] = anova1(curr_analysis_table.Oldham_lVS_Loss_Loss, anova_regressors_strings);
    [lVS_Rew_Loss_p, lVS_Rew_Loss_tbl, lVS_Rew_Loss_stats] = anova1(curr_analysis_table.Oldham_lVS_Rew_Loss, anova_regressors_strings);
    % I want to make my own plots for these. 
elseif significant == 2 
    [biVS_Loss_Gain_p, biVS_Loss_Gain_tbl, biVS_Loss_Gain_stats] = anova1(curr_analysis_table.Oldham_biVS_Loss_Gain, anova_regressors);
    [biVS_Loss_Loss_p, biVS_Loss_Loss_tbl, biVS_Loss_Loss_stats] = anova1(curr_analysis_table.Oldham_biVS_Loss_Loss, anova_regressors);
    [biVS_Rew_Gain_p, biVS_Rew_Gain_tbl, biVS_Rew_Gain_stats] = anova1(curr_analysis_table.Oldham_biVS_Rew_Gain, anova_regressors);
    [biVS_Rew_Loss_p, biVS_Rew_Loss_tbl, biVS_Rew_Loss_stats] = anova1(curr_analysis_table.Oldham_biVS_Rew_Loss, anova_regressors);
    [lVS_Loss_Gain_p, lVS_Loss_Gain_tbl, lVS_Loss_Gain_stats] = anova1(curr_analysis_table.Oldham_lVS_Loss_Gain, anova_regressors);
    [lVS_Loss_Loss_p, lVS_Loss_Loss_tbl, lVS_Loss_Loss_stats] = anova1(curr_analysis_table.Oldham_lVS_Loss_Loss, anova_regressors);
    [lVS_Rew_Gain_p, lVS_Rew_Gain_tbl, lVS_Rew_Gain_stats] = anova1(curr_analysis_table.Oldham_lVS_Rew_Gain, anova_regressors);
    [lVS_Rew_Loss_p, lVS_Rew_Loss_tbl, lVS_Rew_Loss_stats] = anova1(curr_analysis_table.Oldham_lVS_Rew_Loss, anova_regressors);
    [rVS_Loss_Gain_p, rVS_Loss_Gain_tbl, rVS_Loss_Gain_stats] = anova1(curr_analysis_table.Oldham_rVS_Loss_Gain, anova_regressors);
    [rVS_Loss_Loss_p, rVS_Loss_Loss_tbl, rVS_Loss_Loss_stats] = anova1(curr_analysis_table.Oldham_rVS_Loss_Loss, anova_regressors);
    [rVS_Rew_Gain_p, rVS_Rew_Gain_tbl, rVS_Rew_Gain_stats] = anova1(curr_analysis_table.Oldham_rVS_Rew_Gain, anova_regressors);
    [rVS_Rew_Loss_p, rVS_Rew_Loss_tbl, rVS_Rew_Loss_stats] = anova1(curr_analysis_table.Oldham_rVS_Rew_Loss, anova_regressors);
   
    [biVS_AntLossNoLoss_p,biVS_AntLossNoLoss_tbl,biVS_AntLossNoLoss_stats] = anova1(curr_analysis_table.biVS_mean_LossNoLoss, anova_regressors_strings);
    [biVS_AntGainNoGain_p,biVS_AntGainNoGain_tbl,biVS_AntGainNoGain_stats] = anova1(curr_analysis_table.biVS_mean_GainNoGain, anova_regressors);
    
else

    % bilateral ROIs, no sig for loss
    % p = 0.0537 for gains
    [biPutamen_p, biPutamen_tbl, biPutamen_stats] = anova1(curr_analysis_table.HO_biPutamen, anova_regressors); 
    % p = 0.0845 for gains
    [bivmPFC_p, bivmPFC_tbl, bivmPFC_stats] = anova1(curr_analysis_table.HO_bivmPFC, anova_regressors); 
    [biAmyg_p, biAmyg_tbl, biAmyg_stats] = anova1(curr_analysis_table.HO_biAmyg, anova_regressors);
    [biNAcc_p, biNAcc_tbl, biNAcc_stats] = anova1(curr_analysis_table.HO_biNAcc, anova_regressors);
    % significant for anticipation before a loss during a loss trial
    [biOFC_p, biOFC_tbl, biOFC_stats] = anova1(curr_analysis_table.Oldham_biOFC, anova_regressors);
    [biVS_p, biVS_tbl, biOFC_stats] = anova1(curr_analysis_table.Oldham_biVS, anova_regressors);
    
    % Right unilateral ROIs, no sig for loss
    % significant for gains
    [rOFC_p,rOFC_tbl,rOFC_stats] = anova1(curr_analysis_table.Oldham_rOFC(:),anova_regressors);
    [rNAcc_p,rNAcc_tbl,rNAcc_stats] = anova1(curr_analysis_table.HO_rNAcc,anova_regressors);
    [rAmyg_p,rAmyg_tbl,rAmyg_stats] = anova1(curr_analysis_table.HO_rAmyg,anova_regressors);
    [rPutamen_p,rPutamen_tbl,rPutamen_stats] = anova1(curr_analysis_table.HO_rPutamen,anova_regressors);
    [rVS_p, rVS_tbl, rVS_stats] = anova1(curr_analysis_table.Oldham_rVS,anova_regressors);

    % Left unilateral ROIs, no sig for loss
    % significant for gains, p = 0.0648 for loss
    [lPutamen_p,lPutamen_tbl,lPutamen_stats] = anova1(curr_analysis_table.HO_lPutamen,anova_regressors);
    [lOFC_p, lOFC_tbl, lOFC_stats] = anova1(curr_analysis_table.Oldham_lOFC,anova_regressors);
    % significant for anticipation preceding a loss during a loss trial
    [lVS_p, lVS_tbl, lVS_gains_stats] = anova1(curr_analysis_table.Oldham_lVS,anova_regressors);
    [lAmygp,lAmyg_tbl,lAmyg_stats] = anova1(curr_analysis_table.HO_lAmyg,anova_regressors);
    [lNAcc_p,lNAcc_tbl,lNAcc_stats] = anova1(curr_analysis_table.HO_lNAcc,anova_regressors);
end
%% Archived code for generating new contrast mat files
% [p,tbl,stats] = anova1(curr_analysis_table.lNAcc,anova_regressors);
% multcompare(rOFC_gains_stats)

% None of the following are significant effects. 


% everything below this point is commented out because I finished
% outputting the contrasts and clinical data I needed. This is probably
% still helpful code though so I'm keeping it here as a sort of archive for
% this specific set of analyses.
% %% ok great, now remember to add your repo to the path if you haven't already
% % But now load in the contrast document and clinical data
% 
% load(fullfile(contrastdir, 'mid_contrasts_ann.mat'));
% trilevelfname = filenames(fullfile(clinicaldir, 'Multi*FS.mat'));
% load(trilevelfname{1});
% dsmfname = filenames(fullfile(clinicaldir, 'BrainMAPD_subject_diagnosis_table.mat'));
% load(dsmfname{1})
% %% Alright time to wrangle some doggies
% % Need to match up dsm, trilevel, and contrast data
% 
% curr_analysis_table(:,1) = MIDROI091119Zach.PID;
% curr_analysis_table(:,2:11) = NaN;
% 
% % loop through the PID's in curr_analysis_table and snag the contrasts I
% % want and clinical data I want
% 
% temp_analysis_table = zeros(224,4);
% for sub = 1:length(curr_analysis_table.PID(:))
%     % dsm
%     if isempty(find(medbin.PID(:) == curr_analysis_table.PID(sub))) == 1
%         disp(curr_analysis_table(sub,1))
%         
%     else
%         curr_med_index = find(medbin.PID(:) == curr_analysis_table.PID(sub));
%         temp_analysis_table(sub, 1:4) = [medbin.PID(curr_med_index),medbin.T1M1sqmedcurrdopyn1(curr_med_index),medbin.T1M1sqmedcurrdopyn2(curr_med_index),medbin.T1M1sqmedcurrdopyn3(curr_med_index)];%, T1MedicationsBF11S1.Anx_lifetime(curr_dsm_index), BrainMAPD_subject_diagnosis_table.Comorbid_lifetime(curr_dsm_index)];
%     end
% end    
%     
% %     % trilevel
%     if isempty(find(trilevelmultimethodFS.id(:) == curr_analysis_table(sub,1))) == 1
%         disp(curr_analysis_table(sub,1))
%         
%     else
%         curr_trilevel_index = find(trilevelmultimethodFS.id(:) == curr_analysis_table(sub,1));
%         curr_analysis_table(sub,5:7) = [trilevelmultimethodFS.GenDis(curr_trilevel_index), trilevelmultimethodFS.Anhedon(curr_trilevel_index), trilevelmultimethodFS.Fears(curr_trilevel_index)];
%     end
%     % Ann contrasts
%     curr_analysis_table(sub, 8) = MIDROI091119Zach.R_OFC_Oldham_ConGain_v_ConNoGain_avg(sub);
%     curr_analysis_table(sub, 9) = MIDROI091119Zach.L_OFC_Oldham2_ConGain_v_ConNoGain_avg(sub);
%     curr_analysis_table(sub, 10) = MIDROI091119Zach.R_VS_Oldham_Con_ConGain_v_ConNoGain_avg(sub);
%     curr_analysis_table(sub, 11) = MIDROI091119Zach.L_VS_Oldham_Con_ConGain_v_ConNoGain_avg(sub);
% end
% 
% 
% 
