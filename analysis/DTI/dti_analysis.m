% I wrote a long extraction script "BrainMAPD_dti_extract.m" that does all
% kinds of good stuff. But it's not really tailored toward what I want to
% do now, comoparison of all tracks in a single analysis, a repeated
% measures model. So that's what I'm going to do now. And I'm also going to
% revamp the extract script a bit so that no analyses are being done within
% it. It will exclusively extract useful information for useful analyses
% completed here. K, party on Wayne

% So what are our tracks, this really will reference folder names. 

track_list = {'L_Amyg_mOFC','L_NAcc_Amyg','L_NAcc_mOFC','R_Amyg_mOFC','R_NAcc_Amyg','R_NAcc_mOFC'};
suppress_figures = 1; % This is a little haphazard but this is easier than going back into the extraction script and removing the plots... so...
% I should say the above will suppress graphs of the distributions of each
% of the white matter outcomes we extract for these tracks.

% Get your toolbox first and set up some paths
repodir = '~/Documents/repo';
%addpath(genpath(repodir))
basedir = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD/dti/'; % basedir, add from above track list to get datadir which is the var you want
datamatdir = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD/dti/mat_files';
figdir = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD/dti/figures';

for track = 1:length(track_list)
    datadir = strcat(basedir,track_list{track});
    disp(strcat('Extracting stats for:', track_list{track}))
    clinicaldir = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD/clinical_data';
    immunedir = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD/immune_data';
    motiondir = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD/dti/nii/data';
    demodir = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD';
    [curr_analysis_table] = BrainMAPD_dti_extract(datadir,clinicaldir,immunedir,motiondir,demodir,datamatdir);
    cumulative_analysis_struct.(track_list{track}) = curr_analysis_table;
    % add those paths though
    if suppress_figures == 1
        close all
    end
end
clear curr_analysis_table
%% Run the model for DSM: This will recreate all the analyses you ran before
% models with DSM disorders non significant though all disorders are
% showing decreases. I'm going to write a wrapper that will pull variables
% from this for the repeated measures model.
clear track
plot_titles = {'L Amyg -> mOFC','L NAcc -> Amyg','L NAcc -> mOFC','R Amyg -> mOFC','R NAcc -> Amyg','R NAcc -> mOFC'};
for track = 1:length(track_list)

    [p,tbl,stats_fa,terms] = anovan(cumulative_analysis_struct.(track_list{track}).FA,{cumulative_analysis_struct.(track_list{track}).anova_regressors_strings,cumulative_analysis_struct.(track_list{track}).FD,cumulative_analysis_struct.(track_list{track}).gender,cumulative_analysis_struct.(track_list{track}).whole_brain_fa},'varnames',{'Diagnosis','Motion','gender','whole_brain_fa'},'continuous',[2,4])
    
    figure();
    [c,m,h]=multcompare(stats_fa)
    title(plot_titles{track},'FontSize', 20)
    temp_file_name = strcat('Diagnosis_',track_list{track},'.jpg');
    saveas(gcf,fullfile(figdir,temp_file_name))
end
%% Run the model for Trilevel factor scores: 
% Ya none of this is significant either. Motion estimates really screwed us
% here.

clear track
plot_titles = {'L Amyg -> mOFC','L NAcc -> Amyg','L NAcc -> mOFC','R Amyg -> mOFC','R NAcc -> Amyg','R NAcc -> mOFC'};
for track = 1:length(track_list)
    trilevel_mdl.(track_list{track}) = fitlm(cumulative_analysis_struct.(track_list{track}),'anhedon ~ FA + FD + gender + whole_brain_fa')
    trilevel_mdl.(track_list{track}).Coefficients
    anova(trilevel_mdl.(track_list{track}),'summary')
    
    disp('Show correlations between each individual')
    figure();
    subplot(1,3,1)
    scatter(cumulative_analysis_struct.(track_list{track}).gendis(:),cumulative_analysis_struct.(track_list{track}).FA(:))
    xlabel('General Distress', 'FontSize', 15);
    ylabel('Mean FA', 'FontSize', 15);
    h5 = lsline();
    h5.LineWidth = 5;
    h5.Color = 'r';
    r5 = corrcoef(cumulative_analysis_struct.(track_list{track}).gendis(:),cumulative_analysis_struct.(track_list{track}).FA(:),'rows','complete');
    disp(r5(1,2));
    str = [' r = ',num2str(r5(1,2))]
    T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str); 
    set(T, 'fontsize', 15, 'verticalalignment', 'top', 'horizontalalignment', 'left');
    %title(plot_titles{track},'FontSize', 30)
    subplot(1,3,2)
    scatter(cumulative_analysis_struct.(track_list{track}).anhedon(:),cumulative_analysis_struct.(track_list{track}).FA(:))
    xlabel('Anhedonia', 'FontSize', 15);
    ylabel('Mean FA', 'FontSize', 15);
    h5 = lsline();
    h5.LineWidth = 5;
    h5.Color = 'r';
    r5 = corrcoef(cumulative_analysis_struct.(track_list{track}).anhedon(:),cumulative_analysis_struct.(track_list{track}).FA(:),'rows','complete');
    disp(r5(1,2));
    str = [' r = ',num2str(r5(1,2))]
    T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str); 
    set(T, 'fontsize', 15, 'verticalalignment', 'top', 'horizontalalignment', 'left');
    title(plot_titles{track},'FontSize', 20)
    subplot(1,3,3)
    scatter(cumulative_analysis_struct.(track_list{track}).fears(:),cumulative_analysis_struct.(track_list{track}).FA(:))
    xlabel('Fears', 'FontSize', 15);
    ylabel('Mean FA', 'FontSize', 15);
    h5 = lsline();
    h5.LineWidth = 5;
    h5.Color = 'r';
    r5 = corrcoef(cumulative_analysis_struct.(track_list{track}).fears(:),cumulative_analysis_struct.(track_list{track}).FA(:),'rows','complete');
    disp(r5(1,2));
    str = [' r = ',num2str(r5(1,2))]
    T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str); 
    set(T, 'fontsize', 15, 'verticalalignment', 'top', 'horizontalalignment', 'left');
    %title(plot_titles{track},'FontSize', 30)
    temp_file_name = strcat('Trilevel_',track_list{track},'.jpg');
    saveas(gcf,fullfile(figdir,temp_file_name))
end

%% Run the model for immune measures:

clear track
for track = 1:length(track_list)
    immune_mdl.(track_list{track}) = fitlm(cumulative_analysis_struct.(track_list{track}),'FA ~ composite_immune + FD + gender + whole_brain_fa')
    immune_mdl.(track_list{track}).Coefficients
    anova(immune_mdl.(track_list{track}),'summary')
    
    figure();
    scatter(cumulative_analysis_struct.(track_list{track}).composite_immune(:),cumulative_analysis_struct.(track_list{track}).FA(:))
    xlabel('Composite Immune Score', 'FontSize', 15);
    ylabel('Mean FA', 'FontSize', 15);
    h5 = lsline();
    h5.LineWidth = 5;
    h5.Color = 'r';
    r5 = corrcoef(cumulative_analysis_struct.(track_list{track}).composite_immune(:),cumulative_analysis_struct.(track_list{track}).FA(:),'rows','complete');
    disp(r5(1,2));
    str = [' r = ',num2str(r5(1,2))]
    T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str); 
    set(T, 'fontsize', 15, 'verticalalignment', 'top', 'horizontalalignment', 'left');
    title(plot_titles{track},'FontSize', 20)
    temp_file_name = strcat('Immune_',track_list{track},'.jpg');
    saveas(gcf,fullfile(figdir,temp_file_name))
end

%% Run the model for neuroticism and BAS

clear track
for track = 1:length(track_list)
    immune_mdl.(track_list{track}) = fitlm(cumulative_analysis_struct.(track_list{track}),'FA ~ BAS + FD + gender + whole_brain_fa')
    immune_mdl.(track_list{track}).Coefficients
    anova(immune_mdl.(track_list{track}),'summary')
    
    figure();
    scatter(cumulative_analysis_struct.(track_list{track}).BAS(:),cumulative_analysis_struct.(track_list{track}).FA(:))
    xlabel('Behavioral Activation Scale', 'FontSize', 15);
    ylabel('Mean FA', 'FontSize', 15);
    h5 = lsline();
    h5.LineWidth = 5;
    h5.Color = 'r';
    r5 = corrcoef(cumulative_analysis_struct.(track_list{track}).BAS(:),cumulative_analysis_struct.(track_list{track}).FA(:),'rows','complete');
    disp(r5(1,2));
    str = [' r = ',num2str(r5(1,2))]
    T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str); 
    set(T, 'fontsize', 15, 'verticalalignment', 'top', 'horizontalalignment', 'left');
    title(plot_titles{track},'FontSize', 20)
    temp_file_name = strcat('BAS_',track_list{track},'.jpg');
    saveas(gcf,fullfile(figdir,temp_file_name))
end

clear track
for track = 1:length(track_list)
    immune_mdl.(track_list{track}) = fitlm(cumulative_analysis_struct.(track_list{track}),'FA ~ neuroticism + FD + gender + whole_brain_fa')
    immune_mdl.(track_list{track}).Coefficients
    anova(immune_mdl.(track_list{track}),'summary')
    
    figure();
    scatter(cumulative_analysis_struct.(track_list{track}).neuroticism(:),cumulative_analysis_struct.(track_list{track}).FA(:))
    xlabel('Behavioral Activation Scale', 'FontSize', 15);
    ylabel('Mean FA', 'FontSize', 15);
    h5 = lsline();
    h5.LineWidth = 5;
    h5.Color = 'r';
    r5 = corrcoef(cumulative_analysis_struct.(track_list{track}).neuroticism(:),cumulative_analysis_struct.(track_list{track}).FA(:),'rows','complete');
    disp(r5(1,2));
    str = [' r = ',num2str(r5(1,2))]
    T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str); 
    set(T, 'fontsize', 15, 'verticalalignment', 'top', 'horizontalalignment', 'left');
    title(plot_titles{track},'FontSize', 20)
    temp_file_name = strcat('Neuroticism_',track_list{track},'.jpg');
    saveas(gcf,fullfile(figdir,temp_file_name))
end

%% But those didn't show anything! So let's try a repeated measures model, which honestly makes more sense anyway
% Gotta get everything in order first. All my white matter tables are
% different lengths. Super obnoxious.

% Got rid of all my find commands in this version too! The world is
% changing and my code is the better for it.
for sub = 1:length(cumulative_analysis_struct.R_NAcc_mOFC.PID)
    
    % curr_outcome_table(sub,1) = cumulative_analysis_struct.R_NAcc_mOFC.PID(sub);
    % Left side tracks
    if sum(cumulative_analysis_struct.L_Amyg_mOFC.PID(:)==cumulative_analysis_struct.R_NAcc_mOFC.PID(sub)) > 0
        curr_outcome_table(sub,1) = cumulative_analysis_struct.L_Amyg_mOFC.FA(cumulative_analysis_struct.L_Amyg_mOFC.PID(:)==cumulative_analysis_struct.R_NAcc_mOFC.PID(sub));
    else
        curr_outcome_table(sub,1) = NaN;
    end

    if sum(cumulative_analysis_struct.L_NAcc_Amyg.PID(:)==cumulative_analysis_struct.R_NAcc_mOFC.PID(sub)) > 0
        curr_outcome_table(sub,2) = cumulative_analysis_struct.L_NAcc_Amyg.FA(cumulative_analysis_struct.L_NAcc_Amyg.PID(:)==cumulative_analysis_struct.R_NAcc_mOFC.PID(sub));
    else
        curr_outcome_table(sub,2) = NaN;
    end
    
    if sum(cumulative_analysis_struct.L_NAcc_mOFC.PID(:)==cumulative_analysis_struct.R_NAcc_mOFC.PID(sub)) > 0
        curr_outcome_table(sub,3) = cumulative_analysis_struct.L_NAcc_mOFC.FA(cumulative_analysis_struct.L_NAcc_mOFC.PID(:)==cumulative_analysis_struct.R_NAcc_mOFC.PID(sub));
    else
        curr_outcome_table(sub,3) = NaN;
    end
    
    % Right side tracks
    if sum(cumulative_analysis_struct.R_Amyg_mOFC.PID(:)==cumulative_analysis_struct.R_NAcc_mOFC.PID(sub)) > 0
        curr_outcome_table(sub,4) = cumulative_analysis_struct.R_Amyg_mOFC.FA(cumulative_analysis_struct.R_Amyg_mOFC.PID(:)==cumulative_analysis_struct.R_NAcc_mOFC.PID(sub));
    else
        curr_outcome_table(sub,4) = NaN;
    end

    if sum(cumulative_analysis_struct.R_NAcc_Amyg.PID(:)==cumulative_analysis_struct.R_NAcc_mOFC.PID(sub)) > 0
        curr_outcome_table(sub,5) = cumulative_analysis_struct.R_NAcc_Amyg.FA(cumulative_analysis_struct.R_NAcc_Amyg.PID(:)==cumulative_analysis_struct.R_NAcc_mOFC.PID(sub));
    else
        curr_outcome_table(sub,5) = NaN;
    end
    
    if sum(cumulative_analysis_struct.R_NAcc_mOFC.PID(:)==cumulative_analysis_struct.R_NAcc_mOFC.PID(sub)) > 0
        curr_outcome_table(sub,6) = cumulative_analysis_struct.R_NAcc_mOFC.FA(cumulative_analysis_struct.R_NAcc_mOFC.PID(:)==cumulative_analysis_struct.R_NAcc_mOFC.PID(sub));
    else
        curr_outcome_table(sub,6) = NaN;
    end

end
%curr_outcome_table = [cumulative_analysis_struct.L_Amyg_mOFC.FA,cumulative_analysis_struct.L_NAcc_Amyg.FA,cumulative_analysis_struct.L_NAcc_mOFC.FA,cumulative_analysis_struct.R_Amyg_mOFC.FA,cumulative_analysis_struct.R_NAcc_Amyg.FA,cumulative_analysis_struct.R_NAcc_mOFC.FA];
curr_outcome_table = array2table(curr_outcome_table); curr_outcome_table.Properties.VariableNames = {'track1','track2','track3','track4','track5','track6'};

diagnosis = cumulative_analysis_struct.R_NAcc_mOFC.anova_regressors_strings;diagnosis = array2table(diagnosis); diagnosis.Properties.VariableNames = {'diagnosis'};

gender = cumulative_analysis_struct.R_NAcc_mOFC.gender; gender = array2table(gender); gender.Properties.VariableNames = {'gender'};

motion = cumulative_analysis_struct.R_NAcc_mOFC.FD; motion = array2table(motion); motion.Properties.VariableNames = {'motion'};

whole_brain_fa = cumulative_analysis_struct.R_NAcc_mOFC.whole_brain_fa; whole_brain_fa = array2table(whole_brain_fa); whole_brain_fa.Properties.VariableNames = {'whole_brain_fa'};

curr_regressors = [diagnosis,gender,motion,whole_brain_fa];

curr_outcome_table = [curr_outcome_table,curr_regressors];

%% This section specifies a repeated measures model taking into account all the tracks at once
% key to tracks:
% 1: L_Amyg_mOFC
% 2: L_NAcc_Amyg
% 3: L_NAcc_mOFC
% 4: R_Amyg_mOFC
% 5: R_NAcc_Amyg
% 6: R_NAcc_mOFC
Tracks = table([1 2 3 4 5 6]','VariableNames',{'WM_Tracks'});
rm = fitrm(curr_outcome_table,'track1-track6~diagnosis+gender+motion+whole_brain_fa','WithinDesign',Tracks);

ranovatbl = ranova(rm,'WithinModel','separatemeans')  
plot(rm)    
    
    
