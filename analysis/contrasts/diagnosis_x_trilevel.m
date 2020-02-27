% I need to resolve the issue of depression being related to elevated OFC
% activation while reduced OFC activation is apparently related to the tri
% level measure of anhedonia. Very strange. But now I get to see how
% symptoms relate to diagnoses so yay!!

% First! Set up repos, directories, the usual

repodir = '~/Documents/repo';
datadir = '/Users/zaz3744/Documents/repo/acnlab_repo/data';
clinicaldir = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD/clinical_data';
addpath(genpath(repodir))

%% Load in data files
load(fullfile(datadir,'tri_level_final.mat'));
load(fullfile(datadir,'mid_significant_contrasts.mat'));
load(fullfile(clinicaldir,'BrainMAPD_clinical_diagnoses_final.mat'))
%%
% Select cases 
for sub = 1:length(curr_analysis_table.PID)
    if isempty(find(trilevel.ID(:) == curr_analysis_table.PID(sub))) == 0  
        curr = find(trilevel.ID(:) == curr_analysis_table.PID(sub));
        GenDis(sub,1) = trilevel.GenDis(curr);
        GenDis(sub,2) = curr_analysis_table.PID(sub);
        Fears(sub,1) = trilevel.Fears(curr);
        Fears(sub,2) = curr_analysis_table.PID(sub);
        Anhed(sub,1) = trilevel.Anhedon(curr);
        Anhed(sub,2) = curr_analysis_table.PID(sub);            
    else
        disp('Missing Trilevel for:')
        disp(curr_analysis_table.PID(sub))
        GenDis(sub,1)=NaN;
        GenDis(sub,2) = curr_analysis_table.PID(sub);
        Fears(sub,1)=NaN;
        Fears(sub,2) = curr_analysis_table.PID(sub);
        Anhed(sub,1)=NaN;
        Anhed(sub,2) = curr_analysis_table.PID(sub);
    end
end

clear sub

for sub = 1:length(curr_analysis_table.PID)
    if isempty(find(clinical_info.PID(:) == curr_analysis_table.PID(sub))) == 0
        curr2 = find(clinical_info.PID(:) == curr_analysis_table.PID(sub));
        Dep(sub,1) = clinical_info.dep_life_any(curr2);
        Dep(sub,2) = curr_analysis_table.PID(sub);
        Anx(sub,1) = clinical_info.anx_life_any(curr2);
        Anx(sub,2) = curr_analysis_table.PID(sub);
        Com(sub,1) = clinical_info.comorbid_life_dep_anx(curr2);
        Com(sub,2) = curr_analysis_table.PID(sub);
        curr_dep(sub,1) = clinical_info.dep_curr_any(curr2);
        curr_dep(sub,2) = curr_analysis_table.PID(sub);
        curr_anx(sub,1) = clinical_info.anx_curr_any(curr2);
        curr_anx(sub,2) = curr_analysis_table.PID(sub);
        curr_com(sub,1) = clinical_info.comorbid_curr_dep_anx(curr2);
        curr_com(sub,2) = curr_analysis_table.PID(sub);
        dep_csr(sub,1) = clinical_info.max_dep_csr(curr2);
        dep_csr(sub,2) = curr_analysis_table.PID(sub);
        anx_csr(sub,1) = clinical_info.max_anx_csr(curr2);
        anx_csr(sub,2) = curr_analysis_table.PID(sub);
    else
        disp('Missing SCID for:')
        disp(curr_analysis_table.PID(sub))
        Dep(sub,1) = NaN;
        Dep(sub,2) = curr_analysis_table.PID(sub);
        Anx(sub,1) = NaN;
        Anx(sub,2) = curr_analysis_table.PID(sub);
        curr_dep(sub,1) = NaN;
        curr_dep(sub,2) = curr_analysis_table.PID(sub);
        curr_anx(sub,1) = NaN;
        curr_anx(sub,2) = curr_analysis_table.PID(sub);
        curr_com(sub,1) = NaN;
        curr_com(sub,2) = curr_analysis_table.PID(sub);
    end
end



%% Create anova input

lifetime_regressors = ones(length(Dep),1);
lifetime_regressors(Dep(:,1)==1) = 2;
lifetime_regressors(Anx(:,1)==1) = 3;
lifetime_regressors(Com(:,1)==1) = 4;

lifetime_diagnosis = cell(length(lifetime_regressors),1);
lifetime_diagnosis(lifetime_regressors(:,1)==1) = {'Healthy'};
lifetime_diagnosis(lifetime_regressors(:,1)==2) = {'Depression'};
lifetime_diagnosis(lifetime_regressors(:,1)==3) = {'Anxiety'};
lifetime_diagnosis(lifetime_regressors(:,1)==4) = {'Comorbidity'};

curr_regressors = ones(length(Dep),1);
curr_regressors(curr_dep(:,1)==1,:) = 2;
curr_regressors(curr_anx(:,1)==1,:) = 3;
curr_regressors(curr_com(:,1)==1,:) = 4;

curr_diagnosis = cell(length(curr_regressors),1);
curr_diagnosis(curr_regressors(:,1)==1) = {'Healthy'};
curr_diagnosis(curr_regressors(:,1)==2) = {'Current Depression'};
curr_diagnosis(curr_regressors(:,1)==3) = {'Current Anxiety'};
curr_diagnosis(curr_regressors(:,1)==4) = {'Current Comorbid'};
%% Stats and graphs baby
% Anhedonia
[p,tbl,stats] = anova1(Anhed(:,1),lifetime_diagnosis)

figure; [c,m,h]=multcompare(stats)

%% General Distress
[p,tbl,stats] = anova1(GenDis(:,1),lifetime_diagnosis)

figure; [c,m,h]=multcompare(stats)

%% Fears

[p,tbl,stats] = anova1(Fears(:,1),lifetime_diagnosis)

figure; [c,m,h]=multcompare(stats)

%% let's get some CSR's going
% load the file that, as it turns out, wasn't hard to create

figure();
subplot(1,2,1)
scatter(dep_csr(dep_csr(:,1) > 0,1),Anhed(dep_csr(:,1) > 0,1))
title("Dep CSR vs Anhedonia")
subplot(1,2,2) 
scatter(anx_csr(anx_csr(:,1) > 0,1),Anhed(anx_csr(:,1) > 0,1)) 
title("Anx CSR vs Anhedonia")

figure(); 
subplot(1,2,1)
scatter(dep_csr(dep_csr(:,1) > 0,1),GenDis(dep_csr(:,1) > 0,1))
title("Dep CSR vs General Distress")
subplot(1,2,2) 
scatter(anx_csr(anx_csr(:,1) > 0,1),GenDis(anx_csr(:,1) > 0,1)) 
title("Anx CSR vs General Distress")

figure();
subplot(1,2,1)
scatter(dep_csr(dep_csr(:,1) > 0,1),Fears(dep_csr(:,1) > 0,1))
title("Dep CSR vs Fears")
subplot(1,2,2) 
scatter(anx_csr(anx_csr(:,1) > 0,1),Fears(anx_csr(:,1) > 0,1)) 
title("Anx CSR vs Fears")




