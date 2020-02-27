% Helloooooooo campers. We gotta get some analyses going on these DTI data
% pronto. First, what are we looking at. You'll soon find out
% because we're starting with descriptives. After that we'll need to judge
% how good or bad everything looks. More info to come following those
% inevitable meetings and email chains. Party on Wayne

% Get your toolbox first and set up some paths
repodir = '~/Documents/repo';
datadir = '/Users/zaz3744/Documents/current_projects/ACNlab/dti_BrainMAPD/L-Amyg_NAcc';
datamatdir = '/Users/zaz3744/Documents/current_projects/ACNlab/dti_BrainMAPD/mat_files';
clinicaldir = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD/clinical_data';
immunedir = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD/immune_data';
track_file = 'RAmyg-RmOFC.mat';
pull_from_txt = 0;

addpath(genpath(repodir))

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

scatter(curr_analysis_table.composite_immune, curr_analysis_table.QA)
corr(curr_analysis_table.composite_immune, curr_analysis_table.QA,'Rows','complete')

%% Run the stats    

anova_regressors = ones(height(curr_analysis_table),1);
anova_regressors(curr_analysis_table.Dep==1) = 2;
anova_regressors(curr_analysis_table.Anx==1) = 3;
anova_regressors(curr_analysis_table.Com==1) = 4;

Diagnosis = cell(size(anova_regressors));
Diagnosis(anova_regressors(:,1)==1) = {'Healthy'};
Diagnosis(anova_regressors(:,1)==2) = {'Depression'};
Diagnosis(anova_regressors(:,1)==3) = {'Anxiety'};
Diagnosis(anova_regressors(:,1)==4) = {'Comorbidity'};

%% Gonna get some immune analyses going here











    

