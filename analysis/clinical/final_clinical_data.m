%% Man I hope this is the last time I worry about this
% Load the things, rename them for consistency across scripts
% pitter patter let's get at 'er

clinicaldir = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD/clinical_data';
clinical_file = 'SCID_T4.mat'; % SCID_T4_complete.mat

load(fullfile(clinicaldir,clinical_file));

SCID_data = BrainMAPD_T4(:,:);%BrainMAPDT4SCIDConsensusDiagnoses71019(:,:);

%% Now just pull all the obnoxiously long variable names out and name them something more reasonable.
PID = SCID_data.PID;
% Lifetime Dx
depression_life = [SCID_data.T4ANYMDElife, SCID_data.T4ANYPDDlife(:), SCID_data.T4ANYUPDlife(:), SCID_data.T4ANYCYClife(:)];
anxiety_life = [SCID_data.T4ANYPDLife(:), SCID_data.T4ANYAGOLife(:), SCID_data.T4ANYSADLife(:), SCID_data.T4ANYSPLife(:), SCID_data.T4ANYGADLife(:), SCID_data.T4ANYADlife(:), SCID_data.T4ANYOCDlife(:), SCID_data.T4ANYSEPlife(:), SCID_data.T4ANYHOAlife(:), SCID_data.T4ANYIADlife(:), SCID_data.T4ANYPTSDlife(:)];
substance_gambling_life = [SCID_data.T4ANYGAMlife(:), SCID_data.T4ANYAUDlife(:), SCID_data.T4ANYSUDlife(:)]; 
trauma_life = [SCID_data.T4ANYASDlife(:), SCID_data.T4ANYPTSDlife(:), SCID_data.T4ANYADJlife(:), SCID_data.T4ANYTRAlife(:)];
psychosis_mania_life = [SCID_data.T4ANYPSYlife(:), SCID_data.T4ANYBDIlife(:), SCID_data.T4ANYBDIIlife(:), SCID_data.T4ANYBSDlife(:)];
eating_life = [SCID_data.T4ANYANlife(:), SCID_data.T4ANYBNlife(:)];

% Curr Dx
depression_curr = [SCID_data.T4MDEcurr(:), SCID_data.T4PDDCurr(:), SCID_data.T4UPDcurr(:), SCID_data.T4CYCcurr(:)];
anxiety_curr = [SCID_data.T4PDcurr(:), SCID_data.T4AGOOScurr(:), SCID_data.T4SADcurr(:), SCID_data.T4SPcurr(:), SCID_data.T4GADcurr(:), SCID_data.T4ADcurr(:), SCID_data.T4OCDcurr(:), SCID_data.T4SEPcurr(:), SCID_data.T4HOAcurr(:), SCID_data.T4SSDcurr(:), SCID_data.T4IADcurr(:)]; 
substance_gambling_curr = [SCID_data.T4GAMcurr(:), SCID_data.T4AUDcurr(:), SCID_data.T4SUDcurr(:)]; 
trauma_curr = [SCID_data.T4ASDcurr(:), SCID_data.T4PTSDcurr(:), SCID_data.T4ADJcurr(:), SCID_data.T4TRAcurr(:)];
psychosis_mania_curr = [SCID_data.T4PSYcurr(:), SCID_data.T4BDIIcurr(:), SCID_data.T4MANcurr(:), SCID_data.T4HMANcurr(:), SCID_data.T4BSDcurr(:)];
eating_curr = [SCID_data.T4ANcurr(:), SCID_data.T4BNcurr(:)];

%%
diagnosis_table = [PID,depression_life,anxiety_life,substance_gambling_life,trauma_life,psychosis_mania_life,eating_life,depression_curr,anxiety_curr,substance_gambling_curr,trauma_curr,psychosis_mania_curr,eating_curr];
diagnosis_table = array2table(diagnosis_table);
diagnosis_table.Properties.VariableNames = {'PID','dep_mde_life','dep_pdd_life','dep_upd_life','dep_cyc_life',...
    'anx_pd_life','anx_ago_life','anx_sad_life','anx_sp_life','anx_gad_life','anx_ad_life','anx_ocd_life','anx_sep_life','anx_hoad_life','anx_ssd_life','anx_ias_life',...
    'gambling_life','alcohol_life','substance_life',...
    'trauma_asd_life','trauma_ptsd_life','trauma_adj_life','trauma_tra_life',...
    'psych_life','psych_BDI_life', 'psych_BDII_life','psych_bsd_life',...
    'eat_an_life','eat_bn_life',...
    'dep_mde_curr','dep_pdd_curr','dep_upd_curr','dep_cyc_curr',...
    'anx_pd_curr','anx_ago_curr','anx_sad_curr','anx_sp_curr','anx_gad_curr','anx_ad_curr','anx_ocd_curr','anx_sep_curr','anx_hoad_curr','anx_ssd_curr','anx_ias_curr',...
    'gambling_curr','alcohol_curr','substance_curr',...
    'trauma_asd_curr','trauma_ptsd_curr','trauma_adj_curr','trauma_tra_curr',...
    'psych_curr','psych_BDII_curr', 'psych_mania_curr', 'psych_hypo_curr','psych_bsd_curr',...
    'eat_an_curr','eat_bn_curr'};

dep_life_any = max(depression_life,[],2);
anx_life_any = max(anxiety_life,[],2);
substance_life_any = max(substance_gambling_life,[],2);
trauma_life_any = max(trauma_life,[],2);
psych_life_any = max(psychosis_mania_life,[],2);
eating_life_any = max(eating_life,[],2);

dep_curr_any = max(depression_curr,[],2);
anx_curr_any = max(anxiety_curr,[],2);
substance_curr_any = max(substance_gambling_curr,[],2);
trauma_curr_any = max(trauma_curr,[],2);
psych_curr_any = max(psychosis_mania_curr,[],2);
eating_curr_any = max(eating_curr,[],2);

compiled_table = [PID,dep_life_any,anx_life_any,substance_life_any,trauma_life_any,psych_life_any,eating_life_any,dep_curr_any,anx_curr_any,substance_curr_any,trauma_curr_any,psych_curr_any,eating_curr_any];
compiled_table = array2table(compiled_table);
compiled_table.Properties.VariableNames = {'PID','dep_life_any','anx_life_any','substance_life_any','trauma_life_any','psych_mania_life_any','eating_life_any','dep_curr_any','anx_curr_any','substance_curr_any','trauma_curr_any','psych_mania_curr_any','eating_curr_any'};

%% Create lifetime comorbid table

comorbid_life_anx_dep = compiled_table.dep_life_any + compiled_table.anx_life_any;
comorbid_life_anx_dep(comorbid_life_anx_dep<2)=0;comorbid_life_anx_dep(comorbid_life_anx_dep==2)=1;
comorbid_life_anx_dep = array2table(comorbid_life_anx_dep); comorbid_life_anx_dep.Properties.VariableNames = {'comorbid_life_anx_dep'};
clinical_info = [compiled_table(:,1:3),comorbid_life_anx_dep(:,1),compiled_table(:,4:13)];

save 'BrainMAPD_clinical_diagnoses_final_T4.mat' clinical_info



