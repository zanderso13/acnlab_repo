% I have all four timepoints and I'd like to see who of our healthies
% develops psychopathology PREP

clinicaldir = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD/clinical_data';
load(fullfile('/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD/PID.mat'))
t1 = load(fullfile(clinicaldir,'BrainMAPD_clinical_diagnoses_final_T1.mat'));
t2 = load(fullfile(clinicaldir,'BrainMAPD_clinical_diagnoses_final_T2.mat'));
t3 = load(fullfile(clinicaldir,'BrainMAPD_clinical_diagnoses_final_T3.mat'));
t4 = load(fullfile(clinicaldir,'BrainMAPD_clinical_diagnoses_final_T4.mat'));


for sub = 1:length(PID)
    curr_sub = PID(sub);
    PID(sub) = curr_sub;
    if isempty(t2.clinical_info.dep_life_any(curr_sub==t2.clinical_info.PID(:)))==0
        dep2(sub) = [t2.clinical_info.dep_curr_any(curr_sub==t2.clinical_info.PID(:))];
        anx2(sub) = [t2.clinical_info.anx_curr_any(curr_sub==t2.clinical_info.PID(:))];
        com2(sub) = [t2.clinical_info.comorbid_life_anx_dep(curr_sub==t2.clinical_info.PID(:))];
    else
        disp(curr_sub)
        
    end
    
end

%% same as above for t1
for sub = 1:length(PID)
    curr_sub = PID(sub);
    PID(sub) = curr_sub;
    if isempty(t1.clinical_info.dep_life_any(curr_sub==t1.clinical_info.PID(:)))==0
        dep1(sub) = [t1.clinical_info.dep_curr_any(curr_sub==t1.clinical_info.PID(:))];
        anx1(sub) = [t1.clinical_info.anx_curr_any(curr_sub==t1.clinical_info.PID(:))];
        com1(sub) = [t1.clinical_info.comorbid_life_anx_dep(curr_sub==t1.clinical_info.PID(:))];
    else
        disp(curr_sub)
        
    end
    
end 

%% Just for shits, let's do all the func conn analyses but on diagnostic group
for sub = 1:length(PID)
    curr_sub = PID(sub);
    PID(sub) = curr_sub;
    if isempty(t2.clinical_info.dep_life_any(curr_sub==t2.clinical_info.PID(:)))==0
        deplife2(sub) = [t2.clinical_info.dep_life_any(curr_sub==t2.clinical_info.PID(:))];
        anxlife2(sub) = [t2.clinical_info.anx_life_any(curr_sub==t2.clinical_info.PID(:))];
        comlife2(sub) = [t2.clinical_info.comorbid_life_anx_dep(curr_sub==t2.clinical_info.PID(:))];
    elseif isempty(t1.clinical_info.dep_life_any(curr_sub==t1.clinical_info.PID(:)))==0
        deplife2(sub) = [t1.clinical_info.dep_curr_any(curr_sub==t1.clinical_info.PID(:))];
        anxlife2(sub) = [t1.clinical_info.anx_curr_any(curr_sub==t1.clinical_info.PID(:))];
        comlife2(sub) = [t1.clinical_info.comorbid_life_anx_dep(curr_sub==t1.clinical_info.PID(:))];
    else
        disp(curr_sub)
        
        
    end
    
end
