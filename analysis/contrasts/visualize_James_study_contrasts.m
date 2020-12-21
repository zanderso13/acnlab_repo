datadir = '/Volumes/ZachExternal/ACNlab/BrainMAPD/James_study/flout/consumption';
savedir = '/Users/zaz3744/Documents/repo/acnlab_repo/analysis/contrasts/figures';

fnames1 = filenames(fullfile(datadir,'*/*/*/MID/con_0001.nii'));
dat_con1 = fmri_data(fnames1);
thresh_dat1 = threshold(ttest(dat_con1),.05,'fdr');

fnames2 = filenames(fullfile(datadir,'*/*/*/MID/con_0002.nii'));
dat_con2 = fmri_data(fnames2);
thresh_dat2 = threshold(ttest(dat_con2),.05,'fdr');

fnames3 = filenames(fullfile(datadir,'*/*/*/MID/con_0003.nii'));
dat_con3 = fmri_data(fnames3);
thresh_dat3 = threshold(ttest(dat_con3),.05,'fdr');

fnames4 = filenames(fullfile(datadir,'*/*/*/MID/con_0004.nii'));
dat_con4 = fmri_data(fnames4);
thresh_dat4 = threshold(ttest(dat_con4),.05,'fdr');

fnames5 = filenames(fullfile(datadir,'*/*/*/MID/con_0005.nii'));
dat_con5 = fmri_data(fnames5);
thresh_dat5 = threshold(ttest(dat_con5),.05,'fdr');

% t = tiledlayout(5,1)
% brain1 = nexttile;
montage(thresh_dat1)
keyboard
% saveas(gcf, fullfile(savedir, 'RareVsCommon.pdf'))
% close all

% brain2 = nexttile;
montage(thresh_dat2)
keyboard
% saveas(gcf, fullfile(savedir, 'IncentiveRareVsCommon.pdf'))
% close all

% brain3 = nexttile;
montage(thresh_dat3)
keyboard
% saveas(gcf, fullfile(savedir, 'RewardRareVsCommon.pdf'))
% close all

% brain4 = nexttile;
montage(thresh_dat4)
keyboard
% saveas(gcf, fullfile(savedir, 'PunishmentRareVsCommon.pdf'))
% close all

% brain5 = nexttile;
montage(thresh_dat5)

% saveas(gcf, fullfile(savedir, 'NeutralRareVsCommon.pdf'))
% close all
