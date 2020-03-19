function [curr_analysis_table] = BrainMAPD_dti_extract(datadir,clinicaldir,immunedir,motiondir,demodir,datamatdir)

% Helloooooooo campers. We gotta get some analyses going on these DTI data
% pronto. First, what are we looking at. You'll soon find out
% because we're starting with descriptives. After that we'll need to judge
% how good or bad everything looks. More info to come following those
% inevitable meetings and email chains. Party on Wayne


pull_from_txt = 1; % set to 1 if you want to load from txt files and not from mat files


%% Get framewise displacement to add as regressor
fnames_motion = filenames(fullfile(motiondir,'*txt'));
for sub_motion = 1:length(fnames_motion)
    motion = load(fnames_motion{sub_motion});
    if length(motion) == 129
        FD_sub_id(:,sub_motion) = str2num(fnames_motion{sub_motion}(80:84));
        FD_temp(:,sub_motion) = mean(motion,2);
    else
        FD_sub_id(:,sub_motion) = 0;
        FD_temp(:,sub_motion) = NaN;
    end  
    FD = FD_temp';   
    FD = nanmean(FD,2);
    
end
% FD_sub_id = FD_sub_id';
% FD(find(FD_sub_id==0),:)=[];
%% Snag those files
fnames = filenames(fullfile(datadir,'*.stat.txt'));
if pull_from_txt == 1
    % Load in stats txt files. Let's find out what we're looking at IF
    % NEEDED
    

    for sub = 1:length(fnames)   
        fid = fopen(fnames{sub});
        data = textscan(fid,'%s%s%s');
        fclose(fid);

        num_tracts(sub,1)=str2num(data{1,1}{2,1});
        fa_mean(sub,1)=str2num(data{1,3}{15,1});
        fa_sd(sub,1)=str2num(data{1,3}{16,1});
        qa_mean(sub,1)=str2num(data{1,3}{9,1});
        qa_sd(sub,1)=str2num(data{1,3}{10,1});
        rd_mean(sub,1)=str2num(data{1,3}{21,1});
        rd_sd(sub,1)=str2num(data{1,3}{22,1});
        ad_mean(sub,1)=str2num(data{1,3}{19,1});
        ad_sd(sub,1)=str2num(data{1,3}{20,1});
    end
else
    load(fullfile(datamatdir,track_file))    
end
% It worked! Wow was so scared of that when I saw Frank had spaces in all
% his variable names but matlab saves the day again!

%% Let's get some descriptive stats going on here
all_stats = [num_tracts,fa_mean,fa_sd,qa_mean,qa_sd,rd_mean,rd_sd,ad_mean,ad_sd];
mean(all_stats)

% Alright there are some tiny values here lol
all_stats = array2table(all_stats);
all_stats.Properties.VariableNames = {'num_tracts','fa_mean','fa_sd','qa_mean','qa_sd','rd_mean','rd_sd','ad_mean','ad_sd'};

figure();
histogram(all_stats.num_tracts)
title('Number of tracts')
figure();
histogram(all_stats.fa_mean)
title('Average FA')
figure();
histogram(all_stats.qa_mean)
title('Average QA')
figure();
histogram(all_stats.rd_mean)
title('Average RD')
figure();
histogram(all_stats.ad_mean)
title('Average AD')

%% Create list of subject IDs I have data for
clear sub
for sub = 1:length(fnames)
    sub_id(sub) = str2num(fnames{sub}(80:84));
end
sub_id = sub_id';
FD_final = [sub_id,FD];
% Grab clinical data for those people and check for people who are missing
% data

clear sub
load(fullfile(clinicaldir,'BrainMAPD_clinical_diagnoses_final.mat'))
for sub = 1:length(sub_id)
    if isempty(find(clinical_info.PID(:) == sub_id(sub))) == 0
        curr = find(clinical_info.PID(:) == sub_id(sub));
        curr_dep(sub) = clinical_info.dep_life_any(curr);
        curr_anx(sub) = clinical_info.anx_life_any(curr);
        curr_com(sub) = clinical_info.comorbid_life_dep_anx(curr);
    else
        disp(sub_id(sub))
        curr_dep(sub) = NaN;
        curr_anx(sub) = NaN;
        curr_com(sub) = NaN;
    end
end

    
%% create curr_analysis_table
% load immune data
load(fullfile(immunedir,'BrainMAPD_immune.mat'));

% Pull immune data corresponding to present DTI data

clear sub
for sub = 1:length(sub_id)
    if isempty(find(immune_data.PID(:) == sub_id(sub))) == 0
        curr = find(immune_data.PID(:) == sub_id(sub));
        immune(sub) = immune_data.T1BDicsavg(curr);
    else
        disp(sub_id(sub))
        immune(sub) = NaN;
    end
end

curr_analysis_table = [sub_id,curr_dep',curr_anx',curr_com',immune',fa_mean,qa_mean];    

curr_analysis_table = array2table(curr_analysis_table);
curr_analysis_table.Properties.VariableNames = {'PID','Dep','Anx','Com','composite_immune','FA','QA'};    

figure()
scatter(curr_analysis_table.composite_immune, curr_analysis_table.FA)
% title(strcat(datadir(64:69), ' ', datadir(68:72),' FA'))


corr_FA = corr(curr_analysis_table.composite_immune, curr_analysis_table.FA,'Rows','complete')
corr_QA = corr(curr_analysis_table.composite_immune, curr_analysis_table.QA,'Rows','complete')

%% Gonna get some analyses going here
% going to focus on FA, as per Robin's request. First thing is to get my
% covariates going. So pulling out gender from our original document. This
% is slow but my computer is great so I'm just gonna go for it
load(fullfile(datamatdir,'whole_brain_fa.mat'))
load(fullfile(demodir,'demographics.mat'));
clear sub
for sub = 1:length(curr_analysis_table.PID)
    if isempty(find(BrainMAPDT1S1Demo.PID(:) == curr_analysis_table.PID(sub))) == 0
        curr = find(BrainMAPDT1S1Demo.PID(:) == curr_analysis_table.PID(sub));
        gender(sub,1) = BrainMAPDT1S1Demo.sex(curr);
    else
        disp(sub_id(sub))
        gender(sub,1) = NaN;
    end
end

gender = array2table(gender);
gender.Properties.VariableNames = {'gender'};
curr_analysis_table = [curr_analysis_table,gender,whole_brain_fa];


FD = array2table(FD);
FD.Properties.VariableNames = {'FD'};
curr_analysis_table = [curr_analysis_table,FD];
%% Removing outliers first



curr_analysis_table(find(isoutlier(curr_analysis_table.FA,'mean')),:)=[];

dsm_diagnoses_regressors = [curr_analysis_table.Dep(:),curr_analysis_table.Anx(:),curr_analysis_table.Com(:)];

anova_regressors = ones(height(curr_analysis_table),1);
anova_regressors(dsm_diagnoses_regressors(:,1)==1) = 2;
anova_regressors(dsm_diagnoses_regressors(:,2)==1) = 3;
anova_regressors(dsm_diagnoses_regressors(:,3)==1) = 4;

anova_regressors_strings = cell(size(anova_regressors));
anova_regressors_strings(anova_regressors(:,1)==1) = {'Healthy'};
anova_regressors_strings(anova_regressors(:,1)==2) = {'Depression'};
anova_regressors_strings(anova_regressors(:,1)==3) = {'Anxiety'};
anova_regressors_strings(anova_regressors(:,1)==4) = {'Comorbidity'};

anova_regressors_strings = cell2table(anova_regressors_strings);

curr_analysis_table = [curr_analysis_table, anova_regressors_strings];





    

