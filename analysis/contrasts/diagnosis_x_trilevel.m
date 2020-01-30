% I need to resolve the issue of depression being related to elevated OFC
% activation while reduced OFC activation is apparently related to the tri
% level measure of anhedonia. Very strange. But now I get to see how
% symptoms relate to diagnoses so yay!!

% First! Set up repos, directories, the usual

repodir = '~/Documents/repo';
datadir = '/Users/zaz3744/Documents/repo/acnlab_repo/data';

addpath(genpath(repodir))

%% Load in data files
load(fullfile(datadir,'tri_level_final.mat'));
load(fullfile(datadir,'mid_significant_contrasts.mat'));
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
        disp(curr_analysis_table.PID(sub))
        GenDis(sub,1)=NaN;
        GenDis(sub,2) = curr_analysis_table.PID(sub);
        Fears(sub,1)=NaN;
        Fears(sub,2) = curr_analysis_table.PID(sub);
        Anhed(sub,1)=NaN;
        Anhed(sub,2) = curr_analysis_table.PID(sub);
    end
end

Diagnosis = curr_analysis_table.Diagnosis;



%% Stats and graphs baby
% Anhedonia
[p,tbl,stats] = anova1(Anhed(:,1),Diagnosis)

figure; [c,m,h]=multcompare(stats)

%% General Distress
[p,tbl,stats] = anova1(GenDis(:,1),Diagnosis)

figure; [c,m,h]=multcompare(stats)

%% Fears

[p,tbl,stats] = anova1(Fears(:,1),Diagnosis)

figure; [c,m,h]=multcompare(stats)

