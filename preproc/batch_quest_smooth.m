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
scriptdir = '/projects/b1108/projects/BrainMAPD_func_conn/first_levels/quest_submission';

% Where are all your scripts saved for first levels? i.e. where is the
% acnlab_repo folder? Also where is spm12... you need spm

repodir = '~/repo';

% directories
% first is where your stats files will be output to
directories{1} = '/projects/b1108/projects/BrainMAPD_func_conn/first_levels/first_level_output';
% next is where the preprocessed data is
directories{2} = '/projects/b1108/projects/BrainMAPD_func_conn/fmriprep';
% where the raw data lives (raw meaning before preprocessing)
directories{3} = '/projects/b1108/data/BrainMAPD';
% the timing files for modelling (onsets, durations, names)
directories{4} = '/projects/b1108/projects/BrainMAPD_func_conn/timing_files';
% where framewise displacement files will be saved
directories{5} = '/projects/b1108/projects/BrainMAPD_func_conn/additional_files';

% What run of your task are you looking at?
run = 2;
% What session appears in your raw filenames when in BIDS format?
ses = 2;
% Do you want to overwrite previously estimated first levels or just add to
% what you have? 1 overwrites, 0 adds
overwrite = 0;
% Last thing is janky but bear with me. How long are your participant ID's?
% i.e. 10234 would correspond with a 5 for this variable
ID_length = 5;

%%%%%%% END USER DEFINED %%%%%%%%%%

% This will give you PID, don't change this.
% foldernames = char(filenames(fullfile(directories{2},'sub*/')));
% sublist = foldernames(:,size(foldernames,2)-(ID_length):size(foldernames,2)-1);
% 
% sublist = string(sublist);
% I like this better than the overwrite option. This way, the list of
% subjects to be run will change each time I run the script. The overwrite
% option will still be helpful here to turn things on or off.
file_list = filenames(fullfile(directories{2},strcat('*/ses-',num2str(ses),'/func/sub*MID*run',num2str(run),'*preproc_bold.nii')));
for i = 1:length(file_list)
    sublist{i} = file_list{i}(59:63);
end
sublist = string(sublist);
if overwrite == 0
    smooth_list = filenames(fullfile(directories{2},strcat('*/ses-',num2str(ses),'/func/ssub*MID*run*',num2str(run),'*')));
    counter = 1;
    for sub = 1:length(sublist)
        curr_sub = num2str(sublist(sub));
        if isempty(find(contains(smooth_list,curr_sub)))
            new_list(counter) = sublist(sub);
            counter = counter + 1;
        else
            continue
        end
    end
end

% Run/submit first level script

cd(scriptdir)
for sub = 1:length(new_list)
    PID = sublist(sub);
%     ses = 2;
%     run = 2;
%    overwrite = 0;
%     single_sub_smooth(PID,ses,run,overwrite)

        s = ['#!/bin/bash\n\n'...
     '#SBATCH -A p30954\n'...
     '#SBATCH -p short\n'...
     '#SBATCH -t 00:10:00\n'...  
     '#SBATCH --mem=30G\n\n'...
     'matlab -nodisplay -nosplash -nodesktop -r "addpath(genpath(''' repodir ''')); single_sub_smooth(' num2str(PID) ', ' num2str(ses) ',' num2str(run) ',' num2str(overwrite) '); quit"\n\n'];
   
     scriptfile = fullfile(scriptdir, 'smoothing_script.sh');
     fout = fopen(scriptfile, 'w');
     fprintf(fout, s);
     
     
     !chmod 777 smoothing_script.sh
     !sbatch smoothing_script.sh
     
end