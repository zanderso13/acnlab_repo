% batch script to submit many jobs. Initial version of this will be for my
% local test of adding variables to the run_subject_firstlevel_MID.m
% script. Then I will add components that will allow for integration with
% the cluster.

% Dependencies: spm12, CanlabCore, acnlab_repo

% add paths
% repodir = '~/Documents/repo';
% addpath(genpath(repodir))

% Define input to single sub first levels script
% I think the best way to do this is going to be to find which files exist
% and then defining variables relative to those paths. It's going to be
% clunky and IT WILL CHANGE EVERYTIME YOU MOVE FILES AROUND. But since
% we're in BIDS format... Maybe I can actually make this more dynamic. 

%%%%%%% USER DEFINED %%%%%%%%%%
% Define some paths
% This is going to generate a first level script to be submitted to the
% cluster with each run. Where do you want all these .sh scripts saved?
scriptdir = '/projects/b1108/projects/BrainMAPD_func_analysis/PPI/quest_submission';

% Where are all your scripts saved for first levels? i.e. where is the
% acnlab_repo folder? Also where is spm12... you need spm

repodir = '~/repo';

% directories
% first is where your stats files will be output to
fldir = '/projects/b1108/projects/BrainMAPD_func_analysis/PPI/ppi_fldir/anticipation/HOAmyg';
% next is where the preprocessed data is
datadir = '/projects/b1108/projects/BrainMAPD_func_analysis/fmriprep';
% where the raw data lives (raw meaning before preprocessing)
%directories{3} = '/projects/b1108/data/BrainMAPD';
% the timing files for modelling (onsets, durations, names)
%directories{4} = '/projects/b1108/projects/BrainMAPD_func_analysis/timing_files/final/';
% where framewise displacement files will be saved
%directories{5} = '/projects/b1108/projects/BrainMAPD_func_analysis/first_levels/additional_files';

% What run of your task are you looking at?
run = 1;
% What session appears in your raw filenames when in BIDS format?
ses = 2;
% Do you want to overwrite previously estimated first levels or just add to
% what you have? 1 overwrites, 0 adds
overwrite = 0;
% Last thing is janky but bear with me. How long are your participant ID's?
% i.e. 10234 would correspond with a 5 for this variable
ID_length = 5;
% timing files are named such that you need the string below to specify
% what you're doing
contrast = 'ant'; %'ant' or 'con'
%%%%%%% END USER DEFINED %%%%%%%%%%

% This will give you PID, don't change this.
% foldernames = char(filenames(fullfile(directories{2},'sub*/')));
% sublist = foldernames(:,size(foldernames,2)-(ID_length):size(foldernames,2)-1);
% 
% sublist = string(sublist);
% I like this better than the overwrite option. This way, the list of
% subjects to be run will change each time I run the script. The overwrite
% option will still be helpful here to turn things on or off.
file_list = filenames(fullfile(datadir,strcat('sub-*/ses-',num2str(ses),'/func/ssub*MID*run-0',num2str(run),'*preproc_bold.nii')));
for i = 1:length(file_list)
    sublist{i} = file_list{i}(63:67);
end
sublist = string(sublist);


if overwrite == 0
    fl_list = filenames(fullfile(fldir,strcat('*/ses-',num2str(ses)),strcat('run-',num2str(run)),'/MID/con_0001.nii'));
    counter = 1;
    for sub = 1:length(sublist)
        curr_sub = num2str(sublist(sub));
        if isempty(find(contains(fl_list,curr_sub)))
            new_list(counter) = sublist(sub);
            counter = counter + 1;
        else
            continue
        end
    end
end

% Run/submit first level script

cd(scriptdir)
keyboard
for sub = 1:length(new_list)
    PID = new_list(sub);
%     ses = 2;
%     run = 2;
%     overwrite = 0;
%     run_subject_firstlevel_MID_quest(PID,ses,run,directories,overwrite)
%    run_subject_firstlevel_MID_quest(PID, ses, run,overwrite)
        s = ['#!/bin/bash\n\n'...
     '#SBATCH -A p30954\n'...
     '#SBATCH -p short\n'...
     '#SBATCH -t 00:20:00\n'...  
     '#SBATCH --mem=30G\n\n'...
     'matlab -nodisplay -nosplash -nodesktop -r "addpath(genpath(''' repodir ''')); run_subject_firstlevel_PPI_MID_quest(' num2str(PID) ', ' num2str(ses) ',' num2str(run) ',' num2str(overwrite) '); quit"\n\n'];
%   
     scriptfile = fullfile(scriptdir, 'first_level_script.sh');
     fout = fopen(scriptfile, 'w');
     fprintf(fout, s);
%     
%     
     !chmod 777 first_level_script.sh
     !sbatch first_level_script.sh
%     
end
% I probably want to add flags or warnings that will be easy to reference
% later with respect to motion problems, missing data, corrupted files,
% etc.



% for i=1:length(URSIs)
%     
%     % to make jobs, to be submitted w/ sbatch
%     overwrite = 1;
%     ses = 2;
%     
%     if ~overwrite % if not overwriting, can skip right here, before job submission
%         [subID, ~, ~, sub_str, ses_str] = get_subjID_from_URSI(URSIs{i});    
%         if exist(fullfile(fl_dir, sub_str, ses_str(ses,:), 'sponpain','SPM.mat'),'file') % skip if SPM.mat exists
%             fprintf('skipping %d\n', i)
%             continue
%         end
%     end
%     
%     
%     % to run directly (interactively)
%     %run_subject_firstlevel_acute(URSIs{i}, ses, 1)    
%     %continue
% 
% 
%        s = ['#!/bin/bash\n\n'...
%     '#SBATCH -p blanca-ics\n'...
%     '#SBATCH -n 1\n'...
%     '#SBATCH --time 2:00:00\n'...  
%     '#SBATCH --mem=30G\n\n'...
%     'matlab -nodisplay -nosplash -nodesktop -r "addpath(genpath(''' repodir ''')); run_subject_firstlevel_sponpain(''' URSIs{i} ''', ' num2str(ses) ',' num2str(overwrite) '); quit"\n\n'];
%   
%     scriptfile = fullfile(scriptdir, 'first_level_script.sh');
%     fout = fopen(scriptfile, 'w');
%     fprintf(fout, s);
%     
%     !sbatch first_level_script.sh
% end
