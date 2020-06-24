% gotta change timing files

eprime_dir = '/Users/zaz3744/Desktop/MID_Timing_Files';
cd(eprime_dir)
% import BrainMAPD_MID_T1_EPrime.csv using the GUI

%% create a working table for run 1

curr_table = [BrainMAPDMIDT1EPrime(:,2),BrainMAPDMIDT1EPrime(:,85),BrainMAPDMIDT1EPrime(:,120),BrainMAPDMIDT1EPrime(:,124),BrainMAPDMIDT1EPrime(:,103)];
curr_table.Properties.VariableNames = {'PID','ant','rwd','type','fbk'};

curr_table(isnan(curr_table.ant),:) = [];
%% run 2 timing
trial_ind1 = 1;
trial_ind2 = 48;

for sub = 1:242
    curr_table.fbk(trial_ind1:trial_ind2) = curr_table.fbk(trial_ind1:trial_ind2) - curr_table.ant(trial_ind1);
    curr_table.ant(trial_ind1:trial_ind2) = curr_table.ant(trial_ind1:trial_ind2) - curr_table.ant(trial_ind1);
    if range(curr_table.PID(trial_ind1:trial_ind2)) ~= 0
        curr_table.PID(trial_ind1:trial_ind2)
    end
    trial_ind1 = trial_ind1 + 48;
    trial_ind2 = trial_ind2 + 48;
end

curr_table.ant = curr_table.ant ./ 1000;
curr_table.fbk = curr_table.fbk ./ 1000;

% convert to TRs which make soooo much more sense
curr_table.ant = floor(curr_table.ant ./ 2.05);
curr_table.fbk = floor(curr_table.fbk ./ 2.05);

%% Save new timing files for anticipation
savedir = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD/func_conn/final_timing_files/run-2/anticipation/all_trial_types';

trial_ind1 = 1;
trial_ind2 = 48;
% this part is gonna be a little tricky with different indices
for sub = 1:242
    temp_table = curr_table(trial_ind1:trial_ind2,:); 
    if range(curr_table.PID(trial_ind1:trial_ind2)) ~= 0
        curr_table.PID(trial_ind1:trial_ind2)
    end
    % start with anticipation because it's a little easier
    trial_type_strings = {'antgain_1.5','antgain_5','antgain_0','antloss_1.5','antloss_5','antloss_0', 'antgain_all','antloss_all'};
    for files = 1:length(trial_type_strings)
        if strcmp(trial_type_strings{files}, 'antgain_1.5') == 1
            onsets{files} = temp_table.ant(temp_table.type(:) == 'Win' & temp_table.rwd(:) == 1.5);
            names{files} = trial_type_strings{files};
            durations{files} = ones(length(onsets{files}),1) .* 2;
        elseif strcmp(trial_type_strings{files}, 'antgain_5') == 1
            onsets{files} = temp_table.ant(temp_table.type(:) == 'Win' & temp_table.rwd(:) == 5);
            names{files} = trial_type_strings{files};
            durations{files} = ones(length(onsets{files}),1) .* 2;
        elseif strcmp(trial_type_strings{files}, 'antgain_0') == 1
            onsets{files} = temp_table.ant(temp_table.type(:) == 'Win' & temp_table.rwd(:) == 0);
            names{files} = trial_type_strings{files};
            durations{files} = ones(length(onsets{files}),1) .* 2;
        elseif strcmp(trial_type_strings{files}, 'antloss_1.5') == 1
            onsets{files} = temp_table.ant(temp_table.type(:) == 'Lose' & temp_table.rwd(:) == 1.5);
            names{files} = trial_type_strings{files};
            durations{files} = ones(length(onsets{files}),1) .* 2;
        elseif strcmp(trial_type_strings{files}, 'antloss_5') == 1
            onsets{files} = temp_table.ant(temp_table.type(:) == 'Lose' & temp_table.rwd(:) == 5);
            names{files} = trial_type_strings{files};
            durations{files} = ones(length(onsets{files}),1) .* 2;
        elseif strcmp(trial_type_strings{files}, 'antloss_0') == 1
            onsets{files} = temp_table.ant(temp_table.type(:) == 'Lose' & temp_table.rwd(:) == 0);
            names{files} = trial_type_strings{files};
            durations{files} = ones(length(onsets{files}),1) .* 2;
        elseif strcmp(trial_type_strings{files}, 'antgain_all') == 1
            onsets{files} = temp_table.ant(temp_table.type(:) == 'Win' & temp_table.rwd(:) > 0);
            names{files} = trial_type_strings{files};
            durations{files} = ones(length(onsets{files}),1) .* 2;
        elseif strcmp(trial_type_strings{files}, 'antloss_all') == 1
            onsets{files} = temp_table.ant(temp_table.type(:) == 'Lose' & temp_table.rwd(:) > 0);
            names{files} = trial_type_strings{files};
            durations{files} = ones(length(onsets{files}),1) .* 2;
        end
    end
    temp_file_name = strcat(num2str(curr_table.PID(trial_ind1)),'_anticipation_timing.mat');
    save(fullfile(savedir, temp_file_name), 'onsets', 'names', 'durations');
    trial_ind1 = trial_ind1 + 48;
    trial_ind2 = trial_ind2 + 48;
end

%% Save new timing files for anticipation, specifically for the SPM contrasts I'm hoping to run

savedir = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD/func_conn/final_timing_files/run-2/anticipation/spm_all_vs_0_timing';
clear onsets names durations

trial_ind1 = 1;
trial_ind2 = 48;
% this part is gonna be a little tricky with different indices
for sub = 1:242
    temp_table = curr_table(trial_ind1:trial_ind2,:); 
    if range(curr_table.PID(trial_ind1:trial_ind2)) ~= 0
        curr_table.PID(trial_ind1:trial_ind2)
    end
    % start with anticipation because it's a little easier
    trial_type_strings = {'antgain_all','antloss_all','antgain_0','antloss_0'};
    for files = 1:length(trial_type_strings)   
        if strcmp(trial_type_strings{files}, 'antgain_all') == 1
            onsets{files} = temp_table.ant(temp_table.type(:) == 'Win' & temp_table.rwd(:) > 0);
            names{files} = trial_type_strings{files};
            durations{files} = ones(length(onsets{files}),1) .* 2;
        elseif strcmp(trial_type_strings{files}, 'antloss_all') == 1
            onsets{files} = temp_table.ant(temp_table.type(:) == 'Lose' & temp_table.rwd(:) > 0);
            names{files} = trial_type_strings{files};
            durations{files} = ones(length(onsets{files}),1) .* 2;
        elseif strcmp(trial_type_strings{files}, 'antgain_0') == 1
            onsets{files} = temp_table.ant(temp_table.type(:) == 'Win' & temp_table.rwd(:) == 0);
            names{files} = trial_type_strings{files};
            durations{files} = ones(length(onsets{files}),1) .* 2;
        elseif strcmp(trial_type_strings{files}, 'antloss_0') == 1
            onsets{files} = temp_table.ant(temp_table.type(:) == 'Lose' & temp_table.rwd(:) == 0);
            names{files} = trial_type_strings{files};
            durations{files} = ones(length(onsets{files}),1) .* 2;      
        end
    end
    temp_file_name = strcat(num2str(curr_table.PID(trial_ind1)),'_anticipation_timing.mat');
    save(fullfile(savedir, temp_file_name), 'onsets', 'names', 'durations');
    trial_ind1 = trial_ind1 + 48;
    trial_ind2 = trial_ind2 + 48;
end

%% Save new timing files for consumption

savedir = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD/func_conn/final_timing_files/run-2/consumption/all_trial_types';
clear onsets names durations

trial_ind1 = 1;
trial_ind2 = 48;
% this part is gonna be a little tricky with different indices
for sub = 1:242
    temp_table = curr_table(trial_ind1:trial_ind2,:); 
    if range(curr_table.PID(trial_ind1:trial_ind2)) ~= 0
        curr_table.PID(trial_ind1:trial_ind2)
    end
    % start with anticipation because it's a little easier
    trial_type_strings = {'congain_1.5','congain_5','congain_0','conloss_1.5','conloss_5','conloss_0', 'congain_all','conloss_all'};
    for files = 1:length(trial_type_strings)
        if strcmp(trial_type_strings{files}, 'congain_1.5') == 1
            onsets{files} = temp_table.fbk(temp_table.type(:) == 'Win' & temp_table.rwd(:) == 1.5);
            names{files} = trial_type_strings{files};
            durations{files} = ones(length(onsets{files}),1) .* 2;
        elseif strcmp(trial_type_strings{files}, 'congain_5') == 1
            onsets{files} = temp_table.fbk(temp_table.type(:) == 'Win' & temp_table.rwd(:) == 5);
            names{files} = trial_type_strings{files};
            durations{files} = ones(length(onsets{files}),1) .* 2;
        elseif strcmp(trial_type_strings{files}, 'congain_0') == 1
            onsets{files} = temp_table.fbk(temp_table.type(:) == 'Win' & temp_table.rwd(:) == 0);
            names{files} = trial_type_strings{files};
            durations{files} = ones(length(onsets{files}),1) .* 2;
        elseif strcmp(trial_type_strings{files}, 'conloss_1.5') == 1
            onsets{files} = temp_table.fbk(temp_table.type(:) == 'Lose' & temp_table.rwd(:) == 1.5);
            names{files} = trial_type_strings{files};
            durations{files} = ones(length(onsets{files}),1) .* 2;
        elseif strcmp(trial_type_strings{files}, 'conloss_5') == 1
            onsets{files} = temp_table.fbk(temp_table.type(:) == 'Lose' & temp_table.rwd(:) == 5);
            names{files} = trial_type_strings{files};
            durations{files} = ones(length(onsets{files}),1) .* 2;
        elseif strcmp(trial_type_strings{files}, 'conloss_0') == 1
            onsets{files} = temp_table.fbk(temp_table.type(:) == 'Lose' & temp_table.rwd(:) == 0);
            names{files} = trial_type_strings{files};
            durations{files} = ones(length(onsets{files}),1) .* 2;
        elseif strcmp(trial_type_strings{files}, 'congain_all') == 1
            onsets{files} = temp_table.fbk(temp_table.type(:) == 'Win' & temp_table.rwd(:) > 0);
            names{files} = trial_type_strings{files};
            durations{files} = ones(length(onsets{files}),1) .* 2;
        elseif strcmp(trial_type_strings{files}, 'conloss_all') == 1
            onsets{files} = temp_table.fbk(temp_table.type(:) == 'Lose' & temp_table.rwd(:) > 0);
            names{files} = trial_type_strings{files};
            durations{files} = ones(length(onsets{files}),1) .* 2;
        end
    end
    temp_file_name = strcat(num2str(curr_table.PID(trial_ind1)),'_consumption_timing.mat');
    save(fullfile(savedir, temp_file_name), 'onsets', 'names', 'durations');
    trial_ind1 = trial_ind1 + 48;
    trial_ind2 = trial_ind2 + 48;
end

%% Save new timing files for consumption, specifically for the SPM contrasts I'm hoping to run

savedir = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD/func_conn/final_timing_files/run-2/consumption/spm_all_vs_0_timing';
clear onsets names durations

trial_ind1 = 1;
trial_ind2 = 48;
% this part is gonna be a little tricky with different indices
for sub = 1:242
    temp_table = curr_table(trial_ind1:trial_ind2,:); 
    if range(curr_table.PID(trial_ind1:trial_ind2)) ~= 0
        curr_table.PID(trial_ind1:trial_ind2)
    end
    % start with anticipation because it's a little easier
    trial_type_strings = {'congain_all','conloss_all','congain_0','conloss_0'};
    for files = 1:length(trial_type_strings)   
        if strcmp(trial_type_strings{files}, 'congain_all') == 1
            onsets{files} = temp_table.fbk(temp_table.type(:) == 'Win' & temp_table.rwd(:) > 0);
            names{files} = trial_type_strings{files};
            durations{files} = ones(length(onsets{files}),1) .* 2;
        elseif strcmp(trial_type_strings{files}, 'conloss_all') == 1
            onsets{files} = temp_table.fbk(temp_table.type(:) == 'Lose' & temp_table.rwd(:) > 0);
            names{files} = trial_type_strings{files};
            durations{files} = ones(length(onsets{files}),1) .* 2;
        elseif strcmp(trial_type_strings{files}, 'congain_0') == 1
            onsets{files} = temp_table.fbk(temp_table.type(:) == 'Win' & temp_table.rwd(:) == 0);
            names{files} = trial_type_strings{files};
            durations{files} = ones(length(onsets{files}),1) .* 2;
        elseif strcmp(trial_type_strings{files}, 'conloss_0') == 1
            onsets{files} = temp_table.fbk(temp_table.type(:) == 'Lose' & temp_table.rwd(:) == 0);
            names{files} = trial_type_strings{files};
            durations{files} = ones(length(onsets{files}),1) .* 2;      
        end
    end
    temp_file_name = strcat(num2str(curr_table.PID(trial_ind1)),'_consumption_timing.mat');
    save(fullfile(savedir, temp_file_name), 'onsets', 'names', 'durations');
    trial_ind1 = trial_ind1 + 48;
    trial_ind2 = trial_ind2 + 48;
end
