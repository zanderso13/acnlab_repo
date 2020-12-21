% gotta change timing files

eprime_dir = '/home/zach/Downloads/MID_Timing_Files';
cd(eprime_dir)
% import BrainMAPD_MID_T1_EPrime.csv using the GUI

%% create a working table for run 1

curr_table = [allsubs(:,2),allsubs(:,58),allsubs(:,124),allsubs(:,128),allsubs(:,76),allsubs(:,82),allsubs(:,78)];
curr_table.Properties.VariableNames = {'PID','ant','rwd','type','fbk','motor','acc'};

curr_table(isnan(curr_table.ant),:) = [];
%% run 1 timing
trial_ind1 = 1;
trial_ind2 = 48;
num_subs = round(length(curr_table.PID)/48);
for sub = 1:num_subs
    curr_table.fbk(trial_ind1:trial_ind2) = curr_table.fbk(trial_ind1:trial_ind2) - curr_table.ant(trial_ind1);
    curr_table.motor(trial_ind1:trial_ind2) = curr_table.motor(trial_ind1:trial_ind2) - curr_table.ant(trial_ind1);
    curr_table.ant(trial_ind1:trial_ind2) = curr_table.ant(trial_ind1:trial_ind2) - curr_table.ant(trial_ind1);
    if range(curr_table.PID(trial_ind1:trial_ind2)) ~= 0
        curr_table.PID(trial_ind1:trial_ind2)
    end
    trial_ind1 = trial_ind1 + 48;
    trial_ind2 = trial_ind2 + 48;
end

curr_table.ant = round(curr_table.ant ./ 1000,2);
curr_table.fbk = round(curr_table.fbk ./ 1000,2);
curr_table.motor = round(curr_table.motor ./ 1000,2);

% convert to TRs which make soooo much more sense
% curr_table.ant = floor(curr_table.ant ./ 2.05);
% curr_table.fbk = floor(curr_table.fbk ./ 2.05);
% curr_table.motor = ceil(curr_table.motor ./ 2.05);


%% Save new timing files for anticipation
savedir = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD/func_conn/final_timing_files/run-1/anticipation/separate_trial_types';

trial_ind1 = 1;
trial_ind2 = 48;
% this part is gonna be a little tricky with different indices
for sub = 1:round(length(curr_table.PID)/48);
    temp_table = curr_table(trial_ind1:trial_ind2,:); 
    if range(curr_table.PID(trial_ind1:trial_ind2)) ~= 0
        curr_table.PID(trial_ind1:trial_ind2)
    end
    % start with anticipation because it's a little easier
    trial_type_strings = {'antgain_1.5','antgain_5','antgain_0','antloss_1.5','antloss_5','antloss_0', 'antgain_all','antloss_all','motor_all'};
    for files = 1:length(trial_type_strings)
        if strcmp(trial_type_strings{files}, 'antgain_1.5') == 1
            onsets{files} = temp_table.ant(temp_table.type(:) == 'Win' & temp_table.rwd(:) == 1.5);
            names{files} = trial_type_strings{files};
            durations{files} = ones(length(onsets{files}),1) .* 4;
        elseif strcmp(trial_type_strings{files}, 'antgain_5') == 1
            onsets{files} = temp_table.ant(temp_table.type(:) == 'Win' & temp_table.rwd(:) == 5);
            names{files} = trial_type_strings{files};
            durations{files} = ones(length(onsets{files}),1) .* 4;
        elseif strcmp(trial_type_strings{files}, 'antgain_0') == 1
            onsets{files} = temp_table.ant(temp_table.type(:) == 'Win' & temp_table.rwd(:) == 0);
            names{files} = trial_type_strings{files};
            durations{files} = ones(length(onsets{files}),1) .* 4;
        elseif strcmp(trial_type_strings{files}, 'antloss_1.5') == 1
            onsets{files} = temp_table.ant(temp_table.type(:) == 'Lose' & temp_table.rwd(:) == 1.5);
            names{files} = trial_type_strings{files};
            durations{files} = ones(length(onsets{files}),1) .* 4;
        elseif strcmp(trial_type_strings{files}, 'antloss_5') == 1
            onsets{files} = temp_table.ant(temp_table.type(:) == 'Lose' & temp_table.rwd(:) == 5);
            names{files} = trial_type_strings{files};
            durations{files} = ones(length(onsets{files}),1) .* 4;
        elseif strcmp(trial_type_strings{files}, 'antloss_0') == 1
            onsets{files} = temp_table.ant(temp_table.type(:) == 'Lose' & temp_table.rwd(:) == 0);
            names{files} = trial_type_strings{files};
            durations{files} = ones(length(onsets{files}),1) .* 4;
        elseif strcmp(trial_type_strings{files}, 'antgain_all') == 1
            onsets{files} = temp_table.ant(temp_table.type(:) == 'Win' & temp_table.rwd(:) > 0);
            names{files} = trial_type_strings{files};
            durations{files} = ones(length(onsets{files}),1) .* 4;
        elseif strcmp(trial_type_strings{files}, 'antloss_all') == 1
            onsets{files} = temp_table.ant(temp_table.type(:) == 'Lose' & temp_table.rwd(:) > 0);
            names{files} = trial_type_strings{files};
            durations{files} = ones(length(onsets{files}),1) .* 4;
        elseif strcmp(trial_type_strings{files}, 'motor_all') == 1
            onsets{files} = temp_table.motor;
            names{files} = trial_type_strings{files};
            durations{files} = ones(length(onsets{files}),1) .* 4;
        end
    end
    temp_file_name = strcat(num2str(curr_table.PID(trial_ind1)),'_anticipation_timing.mat');
    save(fullfile(savedir, temp_file_name), 'onsets', 'names', 'durations');
    trial_ind1 = trial_ind1 + 48;
    trial_ind2 = trial_ind2 + 48;
end

%% Save new timing files for anticipation, specifically for the SPM contrasts I'm hoping to run

savedir = '/home/zach/Documents/current_projects/ACNlab/BrainMAPD/timing/anticipation/spm_all_vs_0_timing';
clear onsets names durations

trial_ind1 = 1;
trial_ind2 = 48;
% this part is gonna be a little tricky with different indices
for sub = 1:round(length(curr_table.PID)/48);
    temp_table = curr_table(trial_ind1:trial_ind2,:); 
    if range(curr_table.PID(trial_ind1:trial_ind2)) ~= 0
        curr_table.PID(trial_ind1:trial_ind2)
    end
    % start with anticipation because it's a little easier
    trial_type_strings = {'antgain_all','antloss_all','antgain_0','antloss_0','motor_all'};
    for files = 1:length(trial_type_strings)   
        if strcmp(trial_type_strings{files}, 'antgain_all') == 1
            onsets{files} = temp_table.ant(temp_table.type(:) == 'Win' & temp_table.rwd(:) > 0);
            names{files} = trial_type_strings{files};
            durations{files} = ones(length(onsets{files}),1) .* 4;
        elseif strcmp(trial_type_strings{files}, 'antloss_all') == 1
            onsets{files} = temp_table.ant(temp_table.type(:) == 'Lose' & temp_table.rwd(:) > 0);
            names{files} = trial_type_strings{files};
            durations{files} = ones(length(onsets{files}),1) .* 4;
        elseif strcmp(trial_type_strings{files}, 'antgain_0') == 1
            onsets{files} = temp_table.ant(temp_table.type(:) == 'Win' & temp_table.rwd(:) == 0);
            names{files} = trial_type_strings{files};
            durations{files} = ones(length(onsets{files}),1) .* 4;
        elseif strcmp(trial_type_strings{files}, 'antloss_0') == 1
            onsets{files} = temp_table.ant(temp_table.type(:) == 'Lose' & temp_table.rwd(:) == 0);
            names{files} = trial_type_strings{files};
            durations{files} = ones(length(onsets{files}),1) .* 4;
        elseif strcmp(trial_type_strings{files}, 'motor_all') == 1
            onsets{files} = temp_table.motor;
            names{files} = trial_type_strings{files};
            durations{files} = ones(length(onsets{files}),1) .* 4
        end
    end
    temp_file_name = strcat(num2str(curr_table.PID(trial_ind1)),'_anticipation_timing.mat');
    save(fullfile(savedir, temp_file_name), 'onsets', 'names', 'durations');
    trial_ind1 = trial_ind1 + 48;
    trial_ind2 = trial_ind2 + 48;
end

%% Save new timing files for consumption

savedir = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD/func_conn/final_timing_files/run-1/consumption/separate_trial_types';
clear onsets names durations

trial_ind1 = 1;
trial_ind2 = 48;
% this part is gonna be a little tricky with different indices
for sub = 1:round(length(curr_table.PID)/48);
    temp_table = curr_table(trial_ind1:trial_ind2,:); 
    if range(curr_table.PID(trial_ind1:trial_ind2)) ~= 0
        curr_table.PID(trial_ind1:trial_ind2)
    end
    % start with anticipation because it's a little easier
    trial_type_strings = {'congain_5','congain_1.5','congain_0','conloss_5','conloss_1.5','conloss_0', 'motor_all'}; %'congain_all','conloss_all',
    for files = 1:length(trial_type_strings)
        if strcmp(trial_type_strings{files}, 'congain_5') == 1
            onsets{files} = temp_table.fbk(temp_table.type(:) == 'Win' & temp_table.rwd(:) == 1.5);
            names{files} = trial_type_strings{files};
            durations{files} = ones(length(onsets{files}),1) .* 4;
        elseif strcmp(trial_type_strings{files}, 'congain_1.5') == 1
            onsets{files} = temp_table.fbk(temp_table.type(:) == 'Win' & temp_table.rwd(:) == 5);
            names{files} = trial_type_strings{files};
            durations{files} = ones(length(onsets{files}),1) .* 4;
        elseif strcmp(trial_type_strings{files}, 'congain_0') == 1
            onsets{files} = temp_table.fbk(temp_table.type(:) == 'Win' & temp_table.rwd(:) == 0);
            names{files} = trial_type_strings{files};
            durations{files} = ones(length(onsets{files}),1) .* 4;
        elseif strcmp(trial_type_strings{files}, 'conloss_5') == 1
            onsets{files} = temp_table.fbk(temp_table.type(:) == 'Lose' & temp_table.rwd(:) == 1.5);
            names{files} = trial_type_strings{files};
            durations{files} = ones(length(onsets{files}),1) .* 4;
        elseif strcmp(trial_type_strings{files}, 'conloss_1.5') == 1
            onsets{files} = temp_table.fbk(temp_table.type(:) == 'Lose' & temp_table.rwd(:) == 5);
            names{files} = trial_type_strings{files};
            durations{files} = ones(length(onsets{files}),1) .* 4;
        elseif strcmp(trial_type_strings{files}, 'conloss_0') == 1
            onsets{files} = temp_table.fbk(temp_table.type(:) == 'Lose' & temp_table.rwd(:) == 0);
            names{files} = trial_type_strings{files};
            durations{files} = ones(length(onsets{files}),1) .* 4;
%         elseif strcmp(trial_type_strings{files}, 'congain_all') == 1
%             onsets{files} = temp_table.fbk(temp_table.type(:) == 'Win' & temp_table.rwd(:) > 0);
%             names{files} = trial_type_strings{files};
%             durations{files} = ones(length(onsets{files}),1) .* 4;
%         elseif strcmp(trial_type_strings{files}, 'conloss_all') == 1
%             onsets{files} = temp_table.fbk(temp_table.type(:) == 'Lose' & temp_table.rwd(:) > 0);
%             names{files} = trial_type_strings{files};
%             durations{files} = ones(length(onsets{files}),1) .* 4;
        elseif strcmp(trial_type_strings{files}, 'motor_all') == 1
            onsets{files} = temp_table.motor;
            names{files} = trial_type_strings{files};
            durations{files} = ones(length(onsets{files}),1) .* 4;
        end
    end
    temp_file_name = strcat(num2str(curr_table.PID(trial_ind1)),'_consumption_timing.mat');
    save(fullfile(savedir, temp_file_name), 'onsets', 'names', 'durations');
    trial_ind1 = trial_ind1 + 48;
    trial_ind2 = trial_ind2 + 48;
end

%% Save new timing files for consumption, specifically for the SPM contrasts I'm hoping to run

savedir = '/home/zach/Documents/current_projects/ACNlab/BrainMAPD/timing/consumption/spm_all_vs_0_timing';
clear onsets names durations

trial_ind1 = 1;
trial_ind2 = 48;
% this part is gonna be a little tricky with different indices
for sub = 1:round(length(curr_table.PID)/48);
    temp_table = curr_table(trial_ind1:trial_ind2,:); 
    if range(curr_table.PID(trial_ind1:trial_ind2)) ~= 0
        curr_table.PID(trial_ind1:trial_ind2)
    end
    % start with anticipation because it's a little easier
    trial_type_strings = {'congain_all','conloss_all','congain_0','conloss_0','motor_all'};
    for files = 1:length(trial_type_strings)   
        if strcmp(trial_type_strings{files}, 'congain_all') == 1
            onsets{files} = temp_table.fbk(temp_table.type(:) == 'Win' & temp_table.rwd(:) > 0);
            names{files} = trial_type_strings{files};
            durations{files} = ones(length(onsets{files}),1) .* 4;
        elseif strcmp(trial_type_strings{files}, 'conloss_all') == 1
            onsets{files} = temp_table.fbk(temp_table.type(:) == 'Lose' & temp_table.rwd(:) > 0);
            names{files} = trial_type_strings{files};
            durations{files} = ones(length(onsets{files}),1) .* 4;
        elseif strcmp(trial_type_strings{files}, 'congain_0') == 1
            onsets{files} = temp_table.fbk(temp_table.type(:) == 'Win' & temp_table.rwd(:) == 0);
            names{files} = trial_type_strings{files};
            durations{files} = ones(length(onsets{files}),1) .* 4;
        elseif strcmp(trial_type_strings{files}, 'conloss_0') == 1
            onsets{files} = temp_table.fbk(temp_table.type(:) == 'Lose' & temp_table.rwd(:) == 0);
            names{files} = trial_type_strings{files};
            durations{files} = ones(length(onsets{files}),1) .* 4;  
        elseif strcmp(trial_type_strings{files}, 'motor_all') == 1
            onsets{files} = temp_table.motor;
            names{files} = trial_type_strings{files};
            durations{files} = ones(length(onsets{files}),1) .* 4;
        end
    end
    temp_file_name = strcat(num2str(curr_table.PID(trial_ind1)),'_consumption_timing.mat');
    save(fullfile(savedir, temp_file_name), 'onsets', 'names', 'durations');
    trial_ind1 = trial_ind1 + 48;
    trial_ind2 = trial_ind2 + 48;
end

%% Save new timing files for James project
% This focuses on the consumption phase of the project and is intended to
% break down common versus rare stimuli. This requires that I, not only

savedir = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD/MID_all_trial_types/timing_files';
clear onsets names durations

trial_ind1 = 1;
trial_ind2 = 48;
% this part is gonna be a little tricky with different indices
for sub = 1:round(length(curr_table.PID)/48)
    temp_table = curr_table(trial_ind1:trial_ind2,:); 
    if range(curr_table.PID(trial_ind1:trial_ind2)) ~= 0
        curr_table.PID(trial_ind1:trial_ind2)
    end
    % start with anticipation because it's a little easier
    trial_type_strings = {'congain_5W','congain_5L','congain_1.5W', 'congain_1.5L','congain_0W','congain_0L','conloss_5W','conloss_5L','conloss_1.5W','conloss_1.5L','conloss_0W','conloss_0L', 'motor_all'}; %'congain_all','conloss_all',
    for files = 1:length(trial_type_strings)
        if strcmp(trial_type_strings{files}, 'congain_5W') == 1
            onsets{files} = temp_table.fbk(temp_table.type(:) == 'Win' & temp_table.rwd(:) == 5 & temp_table.acc == 1);
            names{files} = trial_type_strings{files};
            durations{files} = ones(length(onsets{files}),1) .* 4;
        elseif strcmp(trial_type_strings{files}, 'congain_5L') == 1
            onsets{files} = temp_table.fbk(temp_table.type(:) == 'Win' & temp_table.rwd(:) == 5 & temp_table.acc == 0);
            names{files} = trial_type_strings{files};
            durations{files} = ones(length(onsets{files}),1) .* 4;
        elseif strcmp(trial_type_strings{files}, 'congain_1.5W') == 1
            onsets{files} = temp_table.fbk(temp_table.type(:) == 'Win' & temp_table.rwd(:) == 1.5 & temp_table.acc == 1);
            names{files} = trial_type_strings{files};
            durations{files} = ones(length(onsets{files}),1) .* 4;
        elseif strcmp(trial_type_strings{files}, 'congain_1.5L') == 1
            onsets{files} = temp_table.fbk(temp_table.type(:) == 'Win' & temp_table.rwd(:) == 1.5 & temp_table.acc == 0);
            names{files} = trial_type_strings{files};
            durations{files} = ones(length(onsets{files}),1) .* 4;
        elseif strcmp(trial_type_strings{files}, 'congain_0W') == 1
            onsets{files} = temp_table.fbk(temp_table.type(:) == 'Win' & temp_table.rwd(:) == 0 & temp_table.acc == 1);
            names{files} = trial_type_strings{files};
            durations{files} = ones(length(onsets{files}),1) .* 4;
        elseif strcmp(trial_type_strings{files}, 'congain_0L') == 1
            onsets{files} = temp_table.fbk(temp_table.type(:) == 'Win' & temp_table.rwd(:) == 0 & temp_table.acc == 0);
            names{files} = trial_type_strings{files};
            durations{files} = ones(length(onsets{files}),1) .* 4;
        elseif strcmp(trial_type_strings{files}, 'conloss_5W') == 1
            onsets{files} = temp_table.fbk(temp_table.type(:) == 'Lose' & temp_table.rwd(:) == 5 & temp_table.acc == 1);
            names{files} = trial_type_strings{files};
            durations{files} = ones(length(onsets{files}),1) .* 4;
        elseif strcmp(trial_type_strings{files}, 'conloss_5L') == 1
            onsets{files} = temp_table.fbk(temp_table.type(:) == 'Lose' & temp_table.rwd(:) == 5 & temp_table.acc == 0);
            names{files} = trial_type_strings{files};
            durations{files} = ones(length(onsets{files}),1) .* 4;
        elseif strcmp(trial_type_strings{files}, 'conloss_1.5W') == 1
            onsets{files} = temp_table.fbk(temp_table.type(:) == 'Lose' & temp_table.rwd(:) == 1.5 & temp_table.acc == 1);
            names{files} = trial_type_strings{files};
            durations{files} = ones(length(onsets{files}),1) .* 4;
        elseif strcmp(trial_type_strings{files}, 'conloss_1.5L') == 1
            onsets{files} = temp_table.fbk(temp_table.type(:) == 'Lose' & temp_table.rwd(:) == 1.5 & temp_table.acc == 0);
            names{files} = trial_type_strings{files};
            durations{files} = ones(length(onsets{files}),1) .* 4;
        elseif strcmp(trial_type_strings{files}, 'conloss_0W') == 1
            onsets{files} = temp_table.fbk(temp_table.type(:) == 'Lose' & temp_table.rwd(:) == 0 & temp_table.acc == 1);
            names{files} = trial_type_strings{files};
            durations{files} = ones(length(onsets{files}),1) .* 4;
        elseif strcmp(trial_type_strings{files}, 'conloss_0L') == 1
            onsets{files} = temp_table.fbk(temp_table.type(:) == 'Lose' & temp_table.rwd(:) == 0 & temp_table.acc == 0);
            names{files} = trial_type_strings{files};
            durations{files} = ones(length(onsets{files}),1) .* 4;
%         elseif strcmp(trial_type_strings{files}, 'congain_all') == 1
%             onsets{files} = temp_table.fbk(temp_table.type(:) == 'Win' & temp_table.rwd(:) > 0);
%             names{files} = trial_type_strings{files};
%             durations{files} = ones(length(onsets{files}),1) .* 4;
%         elseif strcmp(trial_type_strings{files}, 'conloss_all') == 1
%             onsets{files} = temp_table.fbk(temp_table.type(:) == 'Lose' & temp_table.rwd(:) > 0);
%             names{files} = trial_type_strings{files};
%             durations{files} = ones(length(onsets{files}),1) .* 4;
        elseif strcmp(trial_type_strings{files}, 'motor_all') == 1
            onsets{files} = temp_table.motor;
            names{files} = trial_type_strings{files};
            durations{files} = ones(length(onsets{files}),1) .* 4;
        end
    end
    temp_file_name = strcat(num2str(curr_table.PID(trial_ind1)),'_consumption_timing.mat');
    save(fullfile(savedir, temp_file_name), 'onsets', 'names', 'durations');
    trial_ind1 = trial_ind1 + 48;
    trial_ind2 = trial_ind2 + 48;
end
