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
    if exist(find(trilevel.ID(:) == curr_analysis_table.PID(sub)))
        curr = find(trilevel.ID(:) == curr_analysis_table.PID(sub));
        GenDis(sub,1) = trilevel.GenDis(curr);
        Fears(sub,1) = trilevel.Fears(curr);
        Anhed(sub,1) = trilevel.Anhedon(curr);
    else
        disp(curr_analysis_table.PID(sub))
        GenDis(sub,1)=NaN;
        Fears(sub,1)=NaN;
        Anhed(sub,1)=NaN;
    end
end

Diagnosis = curr_analysis_table.Diagnosis;



%% Stats and graphs baby
[p,tbl,stats] = anova1(Anhed,Diagnosis)

[c,m,h]=multcompare(stats)

[p,tbl,stats] = anova1(AnhedMult,Diagnosis)

[c,m,h]=multcompare(stats)

