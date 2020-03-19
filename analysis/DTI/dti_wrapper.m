% I wrote a long extraction script "BrainMAPD_dti_extract.m" that does all
% kinds of good stuff. But it's not really tailored toward what I want to
% do now, comoparison of all tracks in a single analysis, a repeated
% measures model. So that's what I'm going to do now. And I'm also going to
% revamp the extract script a bit so that no analyses are being done within
% it. It will exclusively extract useful information for useful analyses
% completed here. K, party on Wayne

% So what are our tracks, this really will reference folder names. 

track_list = {'L_Amyg_mOFC','L_Amyg_NAcc','L_NAcc_mOFC','R_Amyg_mOFC','R_Amyg_NAcc','R_NAcc_mOFC'};
suppress_figures = 1; % This is a little haphazard but this is easier than going back into the extraction script and removing the plots... so...

% Get your toolbox first and set up some paths
repodir = '~/Documents/repo';
addpath(genpath(repodir))
basedir = '/Users/zaz3744/Documents/current_projects/ACNlab/dti_BrainMAPD/'; % basedir, add from above track list to get datadir which is the var you want
datamatdir = '/Users/zaz3744/Documents/current_projects/ACNlab/dti_BrainMAPD/mat_files';

for track = 1:length(track_list)
    datadir = strcat(basedir,track_list{track});
    disp(strcat('Extracting stats for:', track_list{track}))
    clinicaldir = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD/clinical_data';
    immunedir = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD/immune_data';
    motiondir = '/Users/zaz3744/Documents/current_projects/ACNlab/dti_BrainMAPD/nii/data';
    demodir = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD';
    [curr_analysis_table] = BrainMAPD_dti_extract(datadir,clinicaldir,immunedir,motiondir,demodir,datamatdir);
    cumulative_analysis_struct.(track_list{track}) = curr_analysis_table;
    % add those paths though
    if suppress_figures == 1
        close all
    end
end

%% Run the model: This will recreate all the analyses you ran before
% models with DSM disorders non significant though all disorders are
% showing decreases. I'm going to write a wrapper that will pull variables
% from this for the repeated measures model.
clear track
for track = 1:length(track_list)

    [p,tbl,stats_fa,terms] = anovan(cumulative_analysis_struct.(track_list{track}).FA,{cumulative_analysis_struct.(track_list{track}).anova_regressors_strings,cumulative_analysis_struct.(track_list{track}).FD,cumulative_analysis_struct.(track_list{track}).gender,cumulative_analysis_struct.(track_list{track}).whole_brain_fa},'varnames',{'Diagnosis','Motion','gender','whole_brain_fa'},'continuous',[2,4])
    
    figure();
    [c,m,h]=multcompare(stats_fa)
end
    
%% But those didn't show anything! So let's try a repeated measures model, which honestly makes more sense anyway


    
    
    
    
    
    
    
    
    
