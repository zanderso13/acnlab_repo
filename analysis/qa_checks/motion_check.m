%% This script is created with the structure of my current first level output in mind
% cd into the first level output directory where sub* folders are. Print
% the result of the find command of one of the highest residual files and
% call it sublist.txt

% You also need to download all the mat files from the FD output you create
% but that should be obvs. 

motiondir = '/Users/zaz3744/Documents/current_projects/ACNlab/MWMH/motion';
fnames = filenames(fullfile(motiondir,'*.mat'));
sublist = importdata(fullfile(motiondir,'sublist_T1.txt'));



for sub = 1:length(sublist)
    curr_id = sublist{sub}(11:13);
    if sum(contains(fnames,curr_id)) == 1
        ind = contains(fnames,curr_id);
        load(fnames{ind});
        PID{sub} = curr_id;
        FD_temp = table2array(FD);
        if length(FD_temp) == 1110
            FD_all(:,sub) = FD_temp;
        else
            fprintf(strcat('Not enough volumes: ',curr_id))
        end
    else
        disp(curr_id)
        fprintf(strcat('No motion data found: ', curr_id))
    end
end



%% Visualize 
% For MWMH rest, if we make the motion cutoff .3, likely to lose 6-15 subs.
% Not bad in the grand scheme and probably justifies the number of motion
% regressors that you used.

hist(mean(FD_all))


