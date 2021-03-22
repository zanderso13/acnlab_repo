% First get a list of subs
% Then need to identify all the subs: spm_select 
% Then I need to run generate the template for the estimates: spm_realign
% Finally, reslice actually generates the motion parameters: spm_reslice


fnames = filenames('/projects/b1108/data/MWMH/*/ses-1/dwi/*.nii.gz');

for sub = 1:length(fnames)
    
    sub_folder = fnames{sub}(31:37);
    datadir = strcat('',sub_folder,'ses-2/dwi');
    sub_num = strcat('sub-',num2str(fnames{sub}(77:81)));
    P = spm_select('FPList', pwd, strcat('^',sub_num,'.*\.nii$'));
    spm_realign(P)
    spm_reslice(P)
end