% group ica
gica_cmd --data gica_sub_list.txt --o '/Volumes/ADAPTlab/gica_output'
% gig-ica
gica_cmd --data gica_sub_list.txt --templates '/Volumes/ADAPTlab/gica_output/gica_cmd_agg__component_ica_.nii' --a gig-ica --o '/Volumes/ADAPTlab/gigica_output'
%load clinical data

clinicaldir = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD/clinical_data';
motiondir = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD/Oldham_ROI_by_diagnosis/motion';

% Bad subjects Schizconnect
%/Volumes/ADAPTlab/fmriprep/sub-A00018129/ses-20100101/func/sub-A00018129_ses-20100101_task-rest_space-MNI152NLin2009cAsym_desc-preproc_bold.nii
%/Volumes/ADAPTlab/fmriprep/sub-A00018129/ses-20110101/func/sub-A00018129_ses-20110101_task-rest_space-MNI152NLin2009cAsym_desc-preproc_bold.nii
%/Volumes/ADAPTlab/fmriprep/sub-A00017294/ses-20110101/func/sub-A00017294_ses-20110101_task-rest_space-MNI152NLin2009cAsym_desc-preproc_bold.nii

%import sub_list.csv
load(fullfile(clinicaldir,'BrainMAPD_clinical_diagnoses_final.mat'))
subnums = sublist.VarName5;
for sub = 1:length(subnums)
    PID(sub,1) = subnums(sub);
    if isempty(find(clinical_info.PID(:) == subnums(sub))) == 0
        curr = find(clinical_info.PID(:) == subnums(sub));
        life_dep(sub,1) = clinical_info.dep_life_any(curr);
        life_anx(sub,1) = clinical_info.anx_life_any(curr);
        life_com(sub,1) = clinical_info.comorbid_life_dep_anx(curr);
    else
        disp(strcat(string(subnums(sub)), ' missing clinical info'))
        life_dep(sub,1) = 0;%NaN;
        life_anx(sub,1) = 0;%NaN;
        life_com(sub,1) = 0;%NaN;
    end
end

clear sub
subnums = sublist.VarName5;
load(fullfile(clinicaldir,'trilevel_factors.mat'));
for sub = 1:length(subnums)
    PID(sub,1) = subnums(sub);
    if isempty(find(trilevel.ID(:) == subnums(sub))) == 0
        curr = find(trilevel.ID(:) == subnums(sub));
        curr_analysis_table(sub,:) = trilevel(curr,:);
    else
        disp(strcat(num2str(subnums(sub)), ' missing clinical info'))       
    end
end

% gotta make sure my indices all line up

motion_fnames = filenames(fullfile(motiondir, '*run1.mat'));

for sub = 1:length(subnums)
    curr_id = string(subnums(sub));
    curr_file = contains(motion_fnames,curr_id);
    load(motion_fnames{curr_file});
    sub_id{sub,1} = motion_fnames{sub}(94:98);
    subj_motion(sub,1) = mean(FD.framewise_displacement);
end

figure(); hist(subj_motion); title('Framewise displacement')
