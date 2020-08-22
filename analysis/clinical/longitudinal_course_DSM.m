% The file is made. See sections of code below for that.

clinicaldir = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD/clinical_data';
load(fullfile(clinicaldir,'longitudinal_diagnosis_table.mat'));


%% I have all four timepoints and I'd like to see who of our healthies
% develops psychopathology PREP

clinicaldir = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD/clinical_data';

t1 = load(fullfile(clinicaldir,'BrainMAPD_clinical_diagnoses_final_T1.mat'));
t2 = load(fullfile(clinicaldir,'BrainMAPD_clinical_diagnoses_final_T2.mat'));
t3 = load(fullfile(clinicaldir,'BrainMAPD_clinical_diagnoses_final_T3.mat'));
t4 = load(fullfile(clinicaldir,'BrainMAPD_clinical_diagnoses_final_T4.mat'));


%% PREP
for sub = 1:length(t4.clinical_info.PID)
    curr_sub = t4.clinical_info.PID(sub);
    PID(sub) = curr_sub;
    if isempty(t1.clinical_info.dep_life_any(t1.clinical_info.PID==t4.clinical_info.PID(sub)))==0
        dep1(sub) = [t1.clinical_info.dep_life_any(t1.clinical_info.PID==t4.clinical_info.PID(sub))];
        anx1(sub) = [t1.clinical_info.anx_life_any(t1.clinical_info.PID==t4.clinical_info.PID(sub))];
        com1(sub) = [t1.clinical_info.comorbid_life_anx_dep(t1.clinical_info.PID==t4.clinical_info.PID(sub))];
    else
        disp(curr_sub)
        dep1(sub) = [t2.clinical_info.dep_life_any(t2.clinical_info.PID==t4.clinical_info.PID(sub))];
        anx1(sub) = [t2.clinical_info.anx_life_any(t2.clinical_info.PID==t4.clinical_info.PID(sub))];
        com1(sub) = [t2.clinical_info.comorbid_life_anx_dep(t2.clinical_info.PID==t4.clinical_info.PID(sub))];
    end
    
    dep2(sub) = t4.clinical_info.dep_life_any(sub);
    anx2(sub) = t4.clinical_info.anx_life_any(sub);
    com2(sub) = t4.clinical_info.comorbid_life_anx_dep(sub);
end


%% PREP
clinical_info = [PID',anx1',anx2',dep1',dep2',com1',com2'];
clinical_info = array2table(clinical_info); clinical_info.Properties.VariableNames = {'PID','Anxiety_T1','Anxiety_T2','Depression_T1','Depression_T2','Comorbid_T1','Comorbid_T2'};

%% Taking a closer look

% The big question is how many people seem to lose lifetime diagnoses from
% T1 to T4 and what those diagnoses are.

% I think one way I could do this would be to find the difference between
% T4 and T1. 

% (-) will indicate problems, i.e. they had a lifetime diagnosis... now
% it's gone
% 0 is kinda what we're going for I guess. Suggests no change
% (+) this one will be important, indicates someone was diagnosed later and
% will be crucial for analysis

dep_diff = clinical_info.Depression_T2 - clinical_info.Depression_T1;
anx_diff = clinical_info.Anxiety_T2 - clinical_info.Anxiety_T1;
com_diff = clinical_info.Comorbid_T2 - clinical_info.Comorbid_T1;

% make binary regressors for subs that will develop psychopathology

dep_reg = dep_diff; dep_reg(dep_reg<1)=0;
anx_reg = anx_diff; anx_reg(anx_reg<1)=0;
com_reg = com_diff; com_reg(com_reg<1)=0;

R = [PID',anx_reg,dep_reg,com_reg]; R = array2table(R); R.Properties.VariableNames = {'PID','anx_dev','dep_dev','com_dev'};

save 'T4_diagnosis_onset.mat' R







