function [sub_dir] = AFQ_single_sub(visualize)

if vargin=0
    visualize=0;
end

repodir = '/home/zaz3744/repo';
savedir = '/home/zaz3744/projects';

spm_dir=char(strcat(repodir, filesep, filesep, 'spm8')); 
addpath(spm_dir); 
afq_dir=char(strcat(repodir, filesep, filesep, 'AFQ')); 
p=genpath(afq_dir); 
addpath(p); 
vista_dir=char(strcat(repodir, filesep, filesep, 'vistasoft')); 
v=genpath(vista_dir); 
addpath(v); 
AFQbase=afq_dir; 
AFQdata=strcat('/home/zaz3744/projects/AFQ_MWMH113CV1_20150520', filesep, 'tpl'); 
AFQfunc=strcat(AFQbase, filesep, 'functions'); 
AFQutil=strcat(AFQbase, filesep, 'utilities'); 
AFQdoc=strcat(AFQbase, filesep, 'doc'); 
AFQgui=strcat(AFQbase, filesep, 'gui'); 
cd(AFQdata); 

roi_dir = '/projects/p30954/AFQ_MWMH113CV1_20150520/tpl/ROIs';
%%
% The directory path to the first example subject within the AFQdata folder
% will be:
sub_dir = fullfile(AFQdata);

% Load the subject's dt6 file (generated from dtiInit).
dt = dtiLoadDt6(fullfile(sub_dir,'dt6.mat'));

%%
% Track every fiber from a mask of white matter voxels. Use 'test' mode to
% track fewer fibers and make the example run quicker.
wholebrainFG = AFQ_WholebrainTractography(dt,'test');

%% Trying out some boutique tracts
roi1 = fullfile(roi_dir,'HO_Amygdala_50prob.nii');
roi2 = fullfile(roi_dir,'OFC_8mmsphere_Oldham.nii');

[afq] = AFQ_run(AFQdata,[0]) ; % just creates a struct with all default params
afq = AFQ_AddNewFiberGroup(afq,'test2',roi1,roi2,true,true,true)



%%
% Visualize the wholebrain fiber group.  Because there are a few hundred
% thousand fibers we will use the 'numfibers' input to AFQ_RenderFibers to
% randomly select 1,000 fibers to render. The 'color' input is used to set
% the rgb values that specify the desired color of the fibers.


if visualize == 1
    AFQ_RenderFibers(wholebrainFG, 'numfibers',1000, 'color', [1 .6 .2]);

    
    % Add a sagittal slice from the subject's b0 image to the plot. First load
    % the b0 image.
    b0 = readFileNifti(fullfile(sub_dir,'t1.nii.gz'));

    % Then add the slice X = -2 to the 3d rendering.
    AFQ_AddImageTo3dPlot(b0,[-2, 0, 0]);
end

save(fullfile(savedir,'AFQ_test.mat'))