% batch script to submit many jobs. Initial version of this will be for my
% local test of adding variables to the run_subject_firstlevel_MID.m
% script. Then I will add components that will allow for integration with
% the cluster.

% Dependencies: spm12, CanlabCore, acnlab_repo

% add paths
% repodir = '/projects/p30985/repos';
% addpath(genpath(repodir))

% Define input to single sub first levels script
% I think the best way to do this is going to be to find which files exist
% and then defining variables relative to those paths. It's going to be
% clunky and IT WILL CHANGE EVERYTIME YOU MOVE FILES AROUND. But since
% we're in BIDS format... Maybe I can actually make this more dynamic. 

%%%%%%% USER DEFINED %%%%%%%%%%
% Define some paths
basedir = '/projects/b1108/projects/MWMH_project';

% This is going to generate a first level script to be submitted to the
% cluster with each run. Where do you want all these .sh scripts saved?
scriptdir = fullfile(basedir,'/first_levels/quest_submission');

% Where are all your scripts saved for first levels? i.e. where is the
% acnlab_repo folder? Also where is spm12... you need spm

repodir = '~/repo';

% directories
% first is where your stats files will be output to
directories{1} = fullfile(basedir,'first_levels');
% next is where the preprocessed data is
directories{2} = fullfile(basedir,'fmriprep');
% where the raw data lives (raw meaning before preprocessing)
% directories{3} = '/projects/b1108/data/BrainMAPD';
directories{3} = '/home/zaz3744/ACNlab/data/MWMH';
% where framewise displacement files will be saved
directories{4} = fullfile(basedir,'first_levels/FD');

% What run of your task are you looking at?
run = 1;
% What session appears in your raw filenames when in BIDS format?
ses = 1;
% Do you want to overwrite previously estimated first levels or just add to
% what you have? 1 overwrites, 0 adds
overwrite = 1;
% Last thing is janky but bear with me. How long are your participant ID's?
% i.e. 10234 would correspond with a 5 for this variable
ID_length = 3;

%%%%%%% END USER DEFINED %%%%%%%%%%

% This will give you PID, don't change this.
foldernames = char(filenames(fullfile(directories{2},'sub*/')));
sublist = foldernames(:,size(foldernames,2)-(ID_length):size(foldernames,2)-1);

sublist = string(sublist);
% I like this better than the overwrite option. This way, the list of
% subjects to be run will change each time I run the script. The overwrite
% option will still be helpful here to turn things on or off.

if overwrite == 0
    fl_list = filenames(fullfile(directories{1},'*',strcat('ses-',num2str(ses)),strcat('run-',num2str(run)),'/rest/SPM.mat'));
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
for sub = 1:length(sublist)
    PID = sublist(sub);
%    ses = 2;
%    run = 1;
%    overwrite = 0;
%     run_subject_firstlevel_MWMH_rest(PID,ses,run,overwrite,directories)

        s = ['#!/bin/bash\n\n'...
     '#SBATCH -A p30985\n'...
     '#SBATCH -p short\n'...
     '#SBATCH -t 00:10:00\n'...  
     '#SBATCH --mem=30G\n\n'...
     'matlab -nodisplay -nosplash -nodesktop -r "addpath(genpath(''' repodir ''')); run_subject_firstlevel_MWMH_rest(' num2str(PID) ', ' num2str(ses) ',' num2str(run) ',' num2str(overwrite) '); quit"\n\n'];
   
     scriptfile = fullfile(scriptdir, 'first_level_script.sh');
     fout = fopen(scriptfile, 'w');
     fprintf(fout, s);
     
     
     !chmod 777 first_level_script.sh
     !sbatch first_level_script.sh
     
end
