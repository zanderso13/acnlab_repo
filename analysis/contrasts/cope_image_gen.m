% This script will deal with the first level files generated in FSL by
% someone before me. This is more of a temporary progress report script so
% I can show Robin brain figures. I need images for each of the following:

% Run1 GainNoGain (cope1) and NoGainGain (cope2)
% Run1 LossNoLoss (cope1) and NoLossLoss (cope2)
% Run2 GainNoGain (cope1) and NoGainGain (cope2)
% Run2 LossNoLoss (cope1) and NoLossLoss (cope2)

% define directories
datadir = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD/MID_FSL_contrasts';
figdir = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD/MID_FSL_contrasts/figures';

%% define list of files
% Run 1 filenames
fnames_cope1_run1_gain = filenames(fullfile(datadir,'MID_FSL_Run1/*/*/*/GainvNoGain/cope1_standard.nii'));
fnames_cope2_run1_gain = filenames(fullfile(datadir,'MID_FSL_Run1/*/*/*/GainvNoGain/cope2_standard.nii'));
fnames_cope1_run1_loss = filenames(fullfile(datadir,'MID_FSL_Run1/*/*/*/LossvNoLoss/cope1_standard.nii'));
fnames_cope2_run1_loss = filenames(fullfile(datadir,'MID_FSL_Run1/*/*/*/LossvNoLoss/cope2_standard.nii'));

% Run 2 filenames
fnames_cope1_run2_gain = filenames(fullfile(datadir,'MID_FSL_Run2/*/*/*/GainvNoGain_Run2/cope1_standard.nii'));
fnames_cope2_run2_gain = filenames(fullfile(datadir,'MID_FSL_Run2/*/*/*/GainvNoGain_Run2/cope2_standard.nii'));
fnames_cope1_run2_loss = filenames(fullfile(datadir,'MID_FSL_Run2/*/*/*/LossvNoLoss_Run2/cope1_standard.nii'));
fnames_cope2_run2_loss = filenames(fullfile(datadir,'MID_FSL_Run2/*/*/*/LossvNoLoss_Run2/cope2_standard.nii'));


%% of course there are uneven numbers of files here.
% This is outrageous lol oh well. Gain - No Gain first

% first load the data in 
data = fmri_data(fnames_cope1_run1_gain);
% threshold the data
thresh_data = threshold(ttest(data),.001,'unc');
% then plot it 
figure();    
montage(thresh_data)
title('Gain - No Gain Run 1','FontSize', 20)
temp_file_name = 'Gain-NoGain_Run1.jpg';
saveas(gcf,fullfile(figdir,temp_file_name))
%% No gain - gain. Sanity check, should just be the opposite pattern of activation
% first load the data in 
data = fmri_data(fnames_cope2_run1_gain);
% threshold the data
thresh_data = threshold(ttest(data),.001,'unc');
% then plot it 
figure();    
montage(thresh_data)
title('No Gain - Gain Run 1','FontSize', 20)
temp_file_name = 'NoGain-Gain_Run1.jpg';
saveas(gcf,fullfile(figdir,temp_file_name))

% sanity check passed. It is, in fact No Gain - Gain contrast
%% Loss - No Loss
% first load the data in 
data = fmri_data(fnames_cope1_run1_loss);
% threshold the data
thresh_data = threshold(ttest(data),.001,'unc');
% then plot it 
figure();    
montage(thresh_data)
title('Loss - No Loss Run 1','FontSize', 20)
temp_file_name = 'Loss-NoLoss_Run1.jpg';
saveas(gcf,fullfile(figdir,temp_file_name))

%% Repeat for session 2
% Gain - No Gain 
% first load the data in 
data = fmri_data(fnames_cope1_run2_gain);
% threshold the data
thresh_data = threshold(ttest(data),.001,'unc');
% then plot it 
figure();    
montage(thresh_data)
title('Gain - No Gain Run 2','FontSize', 20)
temp_file_name = 'Gain-NoGain_Run2.jpg';
saveas(gcf,fullfile(figdir,temp_file_name))

%% Loss - No Loss
% first load the data in 
data = fmri_data(fnames_cope1_run2_loss);
% threshold the data
thresh_data = threshold(ttest(data),.001,'unc');
% then plot it 
figure();    
montage(thresh_data)
title('Loss - No Loss Run 2','FontSize', 20)
temp_file_name = 'Loss-NoLoss_Run2.jpg';
saveas(gcf,fullfile(figdir,temp_file_name))
