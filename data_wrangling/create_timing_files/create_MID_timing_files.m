%% User spec
% specify run 
run_name = 'Run2'; %or Run2
% specify how many subs you should have
% run1 has 266, run2 has 
nsubs = 266;
% where are the timing files
basedir = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD/TimingFiles_082218';
% what contrasts are you hoping to perform? Gainany or specific kinds of
% gain?
gainany_vs_0 = 1; % set to zero if you want to split out 150 and 5
% where do you want SPM timing files to end up?
savedir = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD/func_conn/final_timing_files';
%%


% split up the timing files into separate variables
% anticipation filenames
afnames_win0 = filenames(fullfile(basedir,'*anticipation*',strcat('*',run_name,'*'),'*Win0*txt'));
afnames_win150 = filenames(fullfile(basedir,'*anticipation*',strcat('*',run_name,'*'),'*Win150*txt'));
afnames_win5 = filenames(fullfile(basedir,'*anticipation*',strcat('*',run_name,'*'),'*Win5*txt'));
afnames_loss0 = filenames(fullfile(basedir,'*anticipation*',strcat('*',run_name,'*'),'*Loss0*txt'));
afnames_loss150 = filenames(fullfile(basedir,'*anticipation*',strcat('*',run_name,'*'),'*Loss150*txt'));
afnames_loss5 = filenames(fullfile(basedir,'*anticipation*',strcat('*',run_name,'*'),'*Loss5*txt'));
% consumption filenames for getting a gain or losing a loss
cfnames_hitwin0 = filenames(fullfile(basedir,'*consumption*',strcat('*',run_name,'*'),'*Hit*Win_0*txt'));
cfnames_hitwin = filenames(fullfile(basedir,'*consumption*',strcat('*',run_name,'*'),'*Hit*Win_Feedback*txt'));
cfnames_missloss0 = filenames(fullfile(basedir,'*consumption*',strcat('*',run_name,'*'),'*Miss*Loss_0*txt'));
cfnames_missloss = filenames(fullfile(basedir,'*consumption*',strcat('*',run_name,'*'),'*Miss*Loss_Feedback*txt'));
% consumption filenames for losing a gain or avoiding a loss
cfnames_misswin0 = filenames(fullfile(basedir,'*consumption*',strcat('*',run_name,'*'),'*Miss*Win_0*txt'));
cfnames_misswin = filenames(fullfile(basedir,'*consumption*',strcat('*',run_name,'*'),'*Miss*Win_Feedback*txt'));
cfnames_hitloss0 = filenames(fullfile(basedir,'*consumption*',strcat('*',run_name,'*'),'*Hit*Loss_0*txt'));
cfnames_hitloss = filenames(fullfile(basedir,'*consumption*',strcat('*',run_name,'*'),'*Hit*Loss_Feedback*txt'));
% motor filenames
mfnames = filenames(fullfile(basedir,'*anticipation*',strcat('*',run_name,'*'),'*Motor*txt'));

%% load files
% load in files for anticipation
for sub = 1:nsubs
    adat_win0 = load(afnames_win0{sub});
    adat_win150 = load(afnames_win150{sub});
    adat_win5 = load(afnames_win5{sub});
    adat_loss0 = load(afnames_loss0{sub});
    adat_loss150 = load(afnames_loss150{sub});
    adat_loss5 = load(afnames_loss5{sub});
    motor = load(mfnames{sub});
    
    if gainany_vs_0 == 1  
        % for anticipation, we want everything to start with the cue. Old
        % files started with the fixation so all onsets start 1.9 seconds
        % into the trial. We subtract this 1.9 to model the beginning of
        % the onset
        onsets = {sort([adat_win5(:,1);adat_win150(:,1)]-1.9,1),adat_win0(:,1)-1.9,sort([adat_loss5(:,1);adat_loss150(:,1)]-1.9,1),adat_loss0(:,1)-1.9,motor(:,1)};
        names = {'gain','gain0','loss','loss0','motor'};
        durations = {ones(length(onsets{1}),1) .* 4,ones(length(onsets{2}),1) .* 4,ones(length(onsets{3}),1) .* 4,ones(length(onsets{4}),1) .* 4,ones(length(onsets{5}),1) .* 4};
        tmp_filename = fullfile(savedir,run_name,strcat(afnames_loss0{sub}(103:107),'_ant_timing.mat'));
        save(tmp_filename, 'onsets','durations','names')
    else
        onsets = {adat_win5(:,1)-1.9,adat_win150(:,1)-1.9,adat_win0(:,1)-1.9,adat_loss5(:,1)-1.9,adat_loss150(:,1)-1.9,adat_loss0(:,1)-1.9,motor(:,1)};
        names = {'gain5','gain150','gain0','loss5','loss150','loss0','motor'};
        durations = {ones(length(adat_win5),1) .* 4,ones(length(adat_win150),1) .* 4,ones(length(adat_win0),1) .* 4,ones(length(adat_loss5),1) .* 4,ones(length(adat_loss150),1) .* 4,ones(length(adat_loss0),1) .* 4};
        tmp_filename = fullfile(savedir,run_name,strcat(afnames_loss0{sub}(103:107),'_ant_timing.mat'));
        save(tmp_filename, 'onsets','durations','names')
    end
end    
    
%% load in files for consumption    
for sub = 1:nsubs
    cdat_hitwin = load(cfnames_hitwin{sub});
    cdat_hitwin0 = load(cfnames_hitwin0{sub});
    cdat_missloss = load(cfnames_missloss{sub});
    cdat_missloss0 = load(cfnames_missloss0{sub});
    
    cdat_hitloss = load(cfnames_hitloss{sub});
    cdat_hitloss0 = load(cfnames_hitloss0{sub});
    cdat_misswin = load(cfnames_misswin{sub});
    cdat_misswin0 = load(cfnames_misswin0{sub});
    
    motor = load(mfnames{sub});
    
    if gainany_vs_0 == 1  
        onsets1 = cdat_hitwin;%cdat_hitwin0];
        onsets2 = cdat_misswin;%cdat_misswin0];
        onsets3 = cdat_missloss;%cdat_hitloss0];
        onsets4 = cdat_hitloss;%cdat_missloss0];
        if isempty(onsets1) == 0 && isempty(onsets2) == 0 && isempty(onsets3) == 0 && isempty(onsets4) == 0
            onsets = {sort(onsets1(:,1),1),sort(onsets2(:,1),1),sort(onsets3(:,1),1),sort(onsets4(:,1),1),motor(:,1)};
            names = {'gain','nogain','loss','noloss','motor'};
            durations = {ones(length(onsets{1}),1) .* 4,ones(length(onsets{2}),1) .* 4,ones(length(onsets{3}),1) .* 4,ones(length(onsets{4}),1) .* 4,ones(length(onsets{5}),1) .* 4};
            tmp_filename = fullfile(savedir,run_name,strcat(afnames_loss0{sub}(103:107),'_con_timing.mat'));
            save(tmp_filename, 'onsets','durations','names')
        else
            disp(afnames_loss0{sub}(103:107))
            continue
        end
    else
        % you need to edit this before you can run it
        onsets = {cdat_hitwin(:,1),cdat_hitwin0(:,1),cdat_hitloss(:,1),cdat_hitloss0(:,1),cdat_misswin(:,1),cdat_misswin0(:,1),cdat_missloss(:,1),cdat_missloss0(:,1),motor(:,1)};
        names = {'hitgain','hitgain0','hitloss','hitloss0','missgain','missgain0','missloss','missloss0','motor'};
        durations = {ones(length(onsets{1}),1) .* 4,ones(length(onsets{2}),1) .* 4,ones(length(onsets{3}),1) .* 4,ones(length(onsets{4}),1) .* 4,ones(length(onsets{5}),1) .* 4,ones(length(onsets{5}),1) .* 4,ones(length(onsets{5}),1) .* 4,ones(length(onsets{5}),1) .* 4,ones(length(onsets{5}),1) .* 4};
        tmp_filename = fullfile(savedir,run_name,strcat(afnames_loss0{sub}(103:107),'_con_timing.mat'));
        save(tmp_filename, 'onsets','durations','names')
    end
end        
    