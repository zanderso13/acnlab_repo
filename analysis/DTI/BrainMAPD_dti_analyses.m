% Helloooooooo campers. We gotta get some analyses going on these DTI data
% pronto. First, what are we looking at. You'll soon find out
% because we're starting with descriptives. After that we'll need to judge
% how good or bad everything looks. More info to come following those
% inevitable meetings and email chains. Party on Wayne

% Get your toolbox first and set up some paths
repodir = '~/Documents/repo';
datadir = '/Users/zaz3744/Documents/current_projects/dti_BrainMAPD/Right-Amyg_mOFC_len100';

addpath(genpath(repodir))
%% Snag those files
% Load in stats txt files. Let's find out what we're looking at

fnames = filenames(fullfile(datadir,'*.stat.txt'));

for sub = 1:length(fnames)   
    fid = fopen(fnames{sub});
    data = textscan(fid,'%s%s%s');
    fclose(fid);
    
    num_tracts(sub,1)=str2num(data{1,1}{2,1});
    fa_mean(sub,1)=str2num(data{1,3}{15,1});
    fa_sd(sub,1)=str2num(data{1,3}{16,1});
    qa_mean(sub,1)=str2num(data{1,3}{9,1});
    qa_sd(sub,1)=str2num(data{1,3}{10,1});
    rd_mean(sub,1)=str2num(data{1,3}{21,1});
    rd_sd(sub,1)=str2num(data{1,3}{22,1});
    ad_mean(sub,1)=str2num(data{1,3}{19,1});
    ad_sd(sub,1)=str2num(data{1,3}{20,1});
    
end

% It worked! Wow was so scared of that when I saw Frank had spaces in all
% his variable names but matlab saves the day again!

%% Let's get some descriptive stats going on here
all_stats = [num_tracts,fa_mean,fa_sd,qa_mean,qa_sd,rd_mean,rd_sd,ad_mean,ad_sd];
mean(all_stats)

% Alright there are some tiny values here lol
all_stats = array2table(all_stats);
all_stats.Properties.VariableNames = {'num_tracts','fa_mean','fa_sd','qa_mean','qa_sd','rd_mean','rd_sd','ad_mean','ad_sd'};

figure();
histogram(all_stats.num_tracts)
title('Number of tracts')
figure();
histogram(all_stats.fa_mean)
title('Average FA')
figure();
histogram(all_stats.qa_mean)
title('Average QA')
figure();
histogram(all_stats.rd_mean)
title('Average RD')
figure();
histogram(all_stats.ad_mean)
title('Average AD')



