% I need to resolve the issue of depression being related to elevated OFC
% activation while reduced OFC activation is apparently related to the tri
% level measure of anhedonia. Very strange. But now I get to see how
% symptoms relate to diagnoses so yay!!

% First! Set up repos, directories, the usual

repodir = '~/Documents/repo';
datadir = '/Users/zaz3744/Documents/repo/acnlab_repo/data';
clinicaldir = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD/clinical_data';
figdir = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD/clinical_data/figures';
%addpath(genpath(repodir))

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
% current 
figure();
subplot(1,2,1)
scatter(curr_dep(:,1),Anhed(:,1))
h1 = lsline();
h1.LineWidth = 5;
h1.Color = 'r';
r1 = corrcoef(curr_dep(:,1),Anhed(:,1),'rows','complete');
disp(r1(1,2));
str = [' r = ',num2str(r1(1,2))]
T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str); 
set(T, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');
title("Dep CSR Current vs Positive Affect")
subplot(1,2,2) 
scatter(curr_anx(:,1),Anhed(:,1)) 
h2 = lsline();
h2.LineWidth = 5;
h2.Color = 'r';
r2 = corrcoef(curr_anx(:,1),Anhed(:,1),'rows','complete');
disp(r2(1,2));
str = [' r = ',num2str(r2(1,2))]
T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str); 
set(T, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');
title("Anx CSR Current vs Positive Affect")
temp_file_name = strcat('CurrentCSR_PositiveAffect.jpg');
saveas(gcf,fullfile(figdir,temp_file_name))

figure(); 
subplot(1,2,1)
scatter(curr_dep(:,1),GenDis(:,1))
h3 = lsline();
h3.LineWidth = 5;
h3.Color = 'r';
r3 = corrcoef(curr_dep(:,1),GenDis(:,1),'rows','complete');
disp(r3(1,2));
str = [' r = ',num2str(r3(1,2))]
T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str); 
set(T, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');
title("Dep CSR Current vs General Distress")
subplot(1,2,2) 
scatter(curr_anx(:,1),GenDis(:,1))
h4 = lsline();
h4.LineWidth = 5;
h4.Color = 'r';
r4 = corrcoef(curr_anx(:,1),GenDis(:,1),'rows','complete');
disp(r4(1,2));
str = [' r = ',num2str(r4(1,2))]
T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str); 
set(T, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');
title("Anx CSR Current vs General Distress")
temp_file_name = strcat('CurrentCSR_GeneralDistress.jpg');
saveas(gcf,fullfile(figdir,temp_file_name))
    
figure();
subplot(1,2,1)
scatter(curr_dep(:,1),Fears(:,1))
h5 = lsline();
h5.LineWidth = 5;
h5.Color = 'r';
r5 = corrcoef(curr_dep(:,1),Fears(:,1),'rows','complete');
disp(r5(1,2));
str = [' r = ',num2str(r5(1,2))]
T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str); 
set(T, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');
title("Dep CSR Current vs Fears")
subplot(1,2,2) 
scatter(curr_anx(:,1),Fears(:,1)) 
h6 = lsline();
h6.LineWidth = 5;
h6.Color = 'r';
r6 = corrcoef(curr_anx(:,1),Fears(:,1),'rows','complete');
disp(r6(1,2));
str = [' r = ',num2str(r6(1,2))]
T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str); 
set(T, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');
title("Anx CSR Current vs Fears")
temp_file_name = strcat('CurrentCSR_Fears.jpg');
saveas(gcf,fullfile(figdir,temp_file_name))

% lifetime
figure();
subplot(1,2,1)
scatter(dep_csr(:,1),Anhed(:,1))
h1 = lsline();
h1.LineWidth = 5;
h1.Color = 'r';
r1 = corrcoef(dep_csr(:,1),Anhed(:,1),'rows','complete');
disp(r1(1,2));
str = [' r = ',num2str(r1(1,2))]
T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str); 
set(T, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');
title("Dep CSR Lifetime vs Positive Affect")
subplot(1,2,2) 
scatter(anx_csr(:,1),Anhed(:,1)) 
h2 = lsline();
h2.LineWidth = 5;
h2.Color = 'r';
r2 = corrcoef(anx_csr(:,1),Anhed(:,1),'rows','complete');
disp(r2(1,2));
str = [' r = ',num2str(r2(1,2))]
T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str); 
set(T, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');
title("Anx CSR Lifetime vs Positive Affect")
temp_file_name = strcat('LifetimeCSR_PositiveAffect.jpg');
saveas(gcf,fullfile(figdir,temp_file_name))

figure(); 
subplot(1,2,1)
scatter(dep_csr(:,1),GenDis(:,1))
h3 = lsline();
h3.LineWidth = 5;
h3.Color = 'r';
r3 = corrcoef(dep_csr(:,1),GenDis(:,1),'rows','complete');
disp(r3(1,2));
str = [' r = ',num2str(r3(1,2))]
T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str); 
set(T, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');
title("Dep CSR Lifetime vs General Distress")
subplot(1,2,2) 
scatter(anx_csr(:,1),GenDis(:,1))
h4 = lsline();
h4.LineWidth = 5;
h4.Color = 'r';
r4 = corrcoef(anx_csr(:,1),GenDis(:,1),'rows','complete');
disp(r4(1,2));
str = [' r = ',num2str(r4(1,2))]
T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str); 
set(T, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');
title("Anx CSR Lifetime vs General Distress")
temp_file_name = strcat('LifetimeCSR_GeneralDistress.jpg');
saveas(gcf,fullfile(figdir,temp_file_name))
    
figure();
subplot(1,2,1)
scatter(dep_csr(:,1),Fears(:,1))
h5 = lsline();
h5.LineWidth = 5;
h5.Color = 'r';
r5 = corrcoef(dep_csr(:,1),Fears(:,1),'rows','complete');
disp(r5(1,2));
str = [' r = ',num2str(r5(1,2))]
T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str); 
set(T, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');
title("Dep CSR Lifetime vs Fears")
subplot(1,2,2) 
scatter(anx_csr(:,1),Fears(:,1)) 
h6 = lsline();
h6.LineWidth = 5;
h6.Color = 'r';
r6 = corrcoef(anx_csr(:,1),Fears(:,1),'rows','complete');
disp(r6(1,2));
str = [' r = ',num2str(r6(1,2))]
T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str); 
set(T, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');
title("Anx CSR Lifetime vs Fears")
temp_file_name = strcat('LifetimeCSR_Fears.jpg');
saveas(gcf,fullfile(figdir,temp_file_name))
