% First get a list of subs
% Then need to identify all the subs: spm_select 
% Then I need to run generate the template for the estimates: spm_realign
% Finally, reslice actually generates the motion parameters: spm_reslice


fnames = filenames('/Users/zaz3744/Documents/current_projects/ACNlab/dti_BrainMAPD/nii/data/*');

for sub = 1:length(fnames)
    sub_folder = fnames{sub}(77:81);
    datadir = strcat('/Users/zaz3744/Documents/current_projects/ACNlab/dti_BrainMAPD/nii/home/zaz3744/ACNlab/data/BrainMAPD/',sub_folder,'ses-2/dwi');
    sub_num = strcat('sub-',num2str(fnames{sub}(77:81)));
    P = spm_select('FPList', pwd, strcat('^',sub_num,'.*\.nii$'));
    spm_realign(P)
    spm_reslice(P)
end