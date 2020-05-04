% 4 groups paper code
% You've been sitting on this and you know it bud. Pull your head out of
% your ass and let's get at er

datadir = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD/func_conn/first_levels/first_level_output/run-1';
figdir = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD/Oldham_ROI_by_diagnosis/second_level_output';
motiondir = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD/Oldham_ROI_by_diagnosis/motion';
clinicaldir = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD/clinical_data';
savedir = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD/Oldham_ROI_by_diagnosis/second_level_output';

%% Load in data of all kinds
% gotta make sure my indices all line up

motion_fnames = filenames(fullfile(motiondir, '*run1.mat'));
con1_fnames = filenames(fullfile(datadir,'*/ses-2/run-1/MID/con*1.nii'));

for sub = 1:length(con1_fnames)
    load(motion_fnames{sub});
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
R = [R,life_healthy];
%second_level_regressors = array2table(second_level_regressors); second_level_regressors.Properties.VariableNames = {'anx','dep','com','healthy'};
save(fullfile(savedir,'temp_second_level_regressors.mat'),'R')
%% Replace 0's with scaled contrast numbers so SPM doesn't freak out.
% I think what's happening is that spm is not doing a traditional
% regression through the gui. So inputing a matrix of 1's and 0's doesn't
% cut it. Total sum of vector must = 0, so I need to replace 

conval_anx = sum(life_anx)/sum(life_anx==0);
conval_dep = sum(life_dep)/sum(life_dep==0);
conval_com = sum(life_com)/sum(life_com==0);


life_anx(find(life_anx==0)) = conval_anx * -1;
life_dep(find(life_dep==0)) = conval_dep * -1;
life_com(find(life_com==0)) = conval_com * -1;

final = [life_anx';life_dep';life_com'];

save(fullfile(savedir,'temp_dx_regressors.mat'),'final')

