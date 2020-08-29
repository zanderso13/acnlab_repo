% Outline
% 1. Load in .mat files that contain component information
% 2. load in clinical information
% 3. Visualize data

% 1.
datadir = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD/func_conn/ICA/gigica_output';
fnames = filenames(fullfile(datadir, '*sub*component*nii'));
for sub = 1:length(fnames)
    dat = fmri_data(fnames{sub});
    dat1 = dat;
    dat1.dat(:,sub) = dat.dat(:,1);
    dat2 = dat;
    dat2.dat(:,sub) = dat.dat(:,2);
    dat5 = dat;
    dat5.dat(:,sub) = dat.dat(:,5);
    dat17 = dat;
    dat17.dat(:,sub) = dat.dat(:,17);
    dat21 = dat;
    dat21.dat(:,sub) = dat.dat(:,21);
    
end

% template_fname = filenames(fullfile(datadir,'gica_cmd_mean_component_ica_s1_.nii'));
% 
% dat = fmri_data(template_fname{1});
% dat.dat = net1;




