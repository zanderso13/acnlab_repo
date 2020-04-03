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
directories{5} = '/projects/b1108/projects/BrainMAPD_func_conn/framewise_displacement';

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

%%%%%%% END USER DEFINED %%%%%%%%%%


% This will give you PID, don't change this.
foldernames = char(filenames(fullfile(directories{2},'sub*/')));
sublist = foldernames(:,size(foldernames,2)-(ID_length):size(foldernames,2)-1);

sublist = string(sublist);
% Run/submit first level script
fl_list = filenames(fullfile(directories{1},'*/ses-2/run-1/MID/con_0001.nii'));
counter = 1;
for sub = 1:length(sublist)
    curr_sub = num2str(sublist(sub));
    if isempty(find(contains(fl_list,curr_sub)));
        new_list(counter) = sublist(sub);
        counter = counter + 1;
    else
        continue
    end
end



    
    