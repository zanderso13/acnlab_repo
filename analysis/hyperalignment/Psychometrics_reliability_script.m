% Need to run ICA on subjects smoothed data. I need to do it in a way that
% allows me to corrlate the temporal components that are output with the
% design matrices.

% This line does group ICA for a group of subjects based on an inputted txt
% with subject paths. It's relative, gonna need to play around with it.
% Kind of weird how it trys to take into account what folder you're running
% the script from vs what paths are included in the txt

% gica_cmd --data sublist.txt

% It is hilariously simple but this runs ICA with all the default settings.
% There are a lot of things I need to play with. For instance the -a
% gig-ica option is likely going to be needed for the class discrimination
% we want to do.

%% Now, on to the reliability (kind of for 405) 
% But really, this reliability metric is what you want to input into
% hyperalignment. We'll select voxels based on their correlations with the
% original design matrix of a subject.

% load the group average component. This will be your mask and acts as an
% initial reduction of data. The selection of these IC's appears to be
% manual and involves looking at the avg time course in addition to the
% usual warning signs of noise components. 
ica_dir = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD/func_conn/ICA';
con_dir = '/Users/zaz3744/Documents/current_projects/ACNlab/hyperalignment/first_level_temp';
spm_dir = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD/func_conn/first_levels/first_level_output/first_level_output/';
datadir = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD/func_conn/fmriprep';
timing_dir = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD/func_conn/final_timing_files/run-1';
%%
ica_spatial_mask1 = fmri_data(fullfile(ica_dir,'gica_cmd_mean_component_ica_s_all_.nii'));
ica_spatial_mask1.dat = ica_spatial_mask1.dat(:,5);
ica_spatial_mask1.dat(ica_spatial_mask1.dat>-2 & ica_spatial_mask1.dat<2) = 0;
ica_spatial_mask1.dat(ica_spatial_mask1.dat~=0) = 1;

%%
ica_spatial_mask2 = fmri_data(fullfile(ica_dir,'gica_cmd_mean_component_ica_s_all_.nii'));
ica_spatial_mask2.dat = ica_spatial_mask2.dat(:,7);
ica_spatial_mask2.dat(ica_spatial_mask2.dat>-2 & ica_spatial_mask2.dat<2) = 0;
ica_spatial_mask2.dat(ica_spatial_mask2.dat~=0) = 1;

%% Not very many significant voxels here. Consider dropping
ica_spatial_mask3 = fmri_data(fullfile(ica_dir,'gica_cmd_mean_component_ica_s_all_.nii'));
ica_spatial_mask3.dat = ica_spatial_mask3.dat(:,11);
ica_spatial_mask3.dat(ica_spatial_mask3.dat>-2 & ica_spatial_mask3.dat<2) = 0;
ica_spatial_mask3.dat(ica_spatial_mask3.dat~=0) = 1;

%% Some weird wrap around on the frontal lobe in this image
ica_spatial_mask4 = fmri_data(fullfile(ica_dir,'gica_cmd_mean_component_ica_s_all_.nii'));
ica_spatial_mask4.dat = ica_spatial_mask4.dat(:,21);
ica_spatial_mask4.dat(ica_spatial_mask4.dat>-2 & ica_spatial_mask4.dat<2) = 0;
ica_spatial_mask4.dat(ica_spatial_mask4.dat~=0) = 1;

% %% Not crazy about the ventricles here but maybe hits
% ica_spatial_mask5 = fmri_data(fullfile(ica_dir,'gica_cmd_mean_component_ica_s_all_.nii'));
% ica_spatial_mask5.dat = ica_spatial_mask5.dat(:,23);
% ica_spatial_mask5.dat(ica_spatial_mask5.dat>-2 & ica_spatial_mask5.dat<2) = 0;
% ica_spatial_mask5.dat(ica_spatial_mask5.dat~=0) = 1;

%% create a composite mask with most stable components from prev 4
comp_ica_mask = ica_spatial_mask1;
comp_ica_mask.dat = comp_ica_mask.dat + ica_spatial_mask4.dat; % + ica_spatial_mask3.dat + ica_spatial_mask4.dat;

comp_ica_mask.dat(comp_ica_mask.dat>0) = 1;

%% calculate reliability for each contrast
full_data_fnames = filenames(fullfile(datadir,'*/ses-2/func/ssub*MID*run-01*preproc_bold.nii'));

% The goal is to create a 3D matrix voxel x event x subject that will
% contain the correlation of each voxel timecourse with the predicted time
% course generated from first level analysis for every type of event (loss
% ant, gain ant, gain consump, etc) for every subject. I also need to pull
% the 1000 best voxels for each event as we go. I should do this by their
% index and not by the value itself so that I end with a 1 x event x
% subject set of arrays. Each array will refer to the index of the voxel in
% brain space that is providing the most reliable information in each task

% Then, I need to collapse across all those events. The goal is having a
% single array of 1000 correlations for each subject in addition to a
% second 1000 element array that will link them back to their position in
% brain space.

for sub = 1:length(full_data_fnames)
    subid = full_data_fnames{sub}(105:109);
    % load in a subject brain data
    sub_data = fmri_data(full_data_fnames{sub});
    % apply the mask
    sub_data.dat = sub_data.dat .* comp_ica_mask.dat;
    % load their SPM file for consumption
    spm_file_con = filenames(fullfile(spm_dir,'consumption',strcat('sub-',subid),'ses-2/run-1/MID/SPM.mat'));
    load(spm_file_con{1});
    design_dat_con = [SPM.xX.nKX(:,1),SPM.xX.nKX(:,2)];
    clear SPM
    % load their SPM file for anticipation
    spm_file_ant = filenames(fullfile(spm_dir,'anticipation',strcat('sub-',subid),'ses-2/run-1/MID/SPM.mat'));
    load(spm_file_ant{1});
    design_dat_ant = [SPM.xX.nKX(:,1),SPM.xX.nKX(:,3)];
    % ok now let's correlate time series with BOLD signal from masked data
    corr_data_index = find(sub_data.dat(:,sum(sub_data.dat,1)>0));
    corr_data_index(corr_data_index>242953)=[];
    for contrast = 1:2
        for vox = 1:length(corr_data_index)
            corr_pair = [zscore(sub_data.dat(corr_data_index(vox),3:281))',design_dat_ant(:,contrast)];
            R_temp = corrcoef(corr_pair);
            R_ant_temp(vox,contrast) = R_temp(2,1);
        end   
    end
    % Select the best 1000 voxels for each contrast and create a new
    % variable that references where that voxel is in brain space. 
    [max_corr_ant{sub},max_corr_ant_ind{sub}] = maxk(R_ant_temp,300);
    
    % add the correlation matrix you have for this subject containing all
    % time course correlations with all contrasts to the larger all subject
    % matrix
    R_ant(:,:,sub) = R_ant_temp;
    % Now for consumption
    clear vox
    clear contrast
    for contrast = 1:2
        for vox = 1:length(corr_data_index)       
            corr_pair = [sub_data.dat(corr_data_index(vox),3:281)',design_dat_con(:,contrast)];
            R_temp = corrcoef(corr_pair);
            R_con_temp(vox,contrast) = R_temp(2,1);
        end
    end
    [max_corr_con{sub},max_corr_con_ind{sub}] = maxk(R_con_temp,300);  
    
    % add the correlation matrix you have for this subject containing all
    % time course correlations with all contrasts to the larger all subject
    % matrix
    R_con(:,:,sub) = R_con_temp;
    
    % Combine the time courses to actually pull out 1000 unique voxels
    design_dat_total = [design_dat_con, design_dat_ant];
    design_dat_total = sum(design_dat_total,2);
    clear vox
    clear contrast
    
    for vox = 1:length(corr_data_index)       
        corr_pair = [sub_data.dat(corr_data_index(vox),3:281)',design_dat_total(:,1)];
        R_temp = corrcoef(corr_pair);
        R_total_temp(vox) = R_temp(2,1);
    end

    [max_corr_total{sub},max_corr_total_ind{sub}] = maxk(R_total_temp,1000);  
    
    % add the correlation matrix you have for this subject containing all
    % time course correlations with all contrasts to the larger all subject
    % matrix
    R_total(:,:,sub) = R_total_temp;
    
    %% New section! Same loop 
    % I need to snag the average image for each instance of each event
    % while I have it all here. So I'm thinking that I load those big
    % timing files that I made that have all contrasts. No I just need to
    % use the onsets and take the average of each image. I don't want this
    % to be modeled data. Ok so let's load the broader onsets file.
    
    temp_timing_file = filenames(fullfile(timing_dir,strcat(subid,'.mat')));
    load(temp_timing_file{1})
    
end    
 
%% Variance/covariance matrix
% Bill had a great though in that I should really be combining my trial
% types. There are going to be voxels that reliably activate to both,
% one or neither type of trial. So let's do that. But instead of just
% adding the time courses (which is a mess... multicolinearity
% everywhere). So instead, I extracted the top 1000 voxels with
% respect to each contrast of interest (8000 total). If voxels are
% represented in more than one contrast I'll pull those out first. Pull
% these from the max_corr and max_corr_ind cell matrices

% I'm going to take 500 voxels from all anticipation events and 500 from all consumption events
clear sub
for sub = 1:length(max_corr_ant)
    corr_total(:,sub) = [max_corr_ant{sub}(:,1);max_corr_ant{sub}(:,2);max_corr_con{sub}(:,1);max_corr_con{sub}(:,2)];    
    corr_total_index(:,sub) = [max_corr_ant_ind{sub}(:,1);max_corr_ant_ind{sub}(:,2);max_corr_con_ind{sub}(:,1);max_corr_con_ind{sub}(:,2)];
end

%% create corr matrices

clear sub
for sub = 1:size(corr_total_index,2)
    sub_data = fmri_data(full_data_fnames{sub});
    corr_input = sub_data.dat(unique(corr_total_index(:,sub)),3:281); 
    final_brain_coordinates_by_contrast{sub} = corr_data_index(unique(corr_total_index(:,sub)));
    correlation_matrices_by_contrast{sub} = corrcoef(corr_input');
end

%% create corr matrices for the R values where I combined the time courses
clear corr_input
clear sub
for sub = 1:size(max_corr_total_ind,2)
    sub_data = fmri_data(full_data_fnames{sub});
    corr_input(:,:,sub) = sub_data.dat(corr_data_index(max_corr_total_ind{sub}),3:281); 
    final_brain_coordinates_total{sub} = corr_data_index(max_corr_total_ind{sub});
    correlation_matrices_total{sub} = corrcoef(corr_input(:,:,sub)');
end

%% Visualize brain regions projected back onto brain
indiv_overlap_mask = comp_ica_mask;
control_overlap_mask = comp_ica_mask;
dsm_overlap_mask = comp_ica_mask;
patient_vs_control_overlap_mask = comp_ica_mask;

control_subs = [0,1,1,1,0,0,0,1,1,0];

for sub = 1:length(final_brain_coordinates_total)
    indiv_overlap_mask.dat(final_brain_coordinates_total{sub}) = indiv_overlap_mask.dat(final_brain_coordinates_total{sub})+10;
    if control_subs == 1
        control_overlap_mask.dat(final_brain_coordinates_total{sub}) = control_overlap_mask.dat(final_brain_coordinates_total{sub})+1;
    else
        dsm_overlap_mask.dat(final_brain_coordinates_total{sub}) = dsm_overlap_mask.dat(final_brain_coordinates_total{sub})-1;
    end
end

figure(); montage(indiv_overlap_mask)

patient_vs_control_overlap_mask.dat = patient_vs_control_overlap_mask.dat + control_overlap_mask.dat;
patient_vs_control_overlap_mask.dat = patient_vs_control_overlap_mask.dat + dsm_overlap_mask.dat;
patient_vs_control_overlap_mask.dat(patient_vs_control_overlap_mask.dat==1) = 0;
figure(); montage(patient_vs_control_overlap_mask)
%% I need to align the time points and then hyperalign.
% I'll be able to compare the hyperaligned matrix and the original data in
% two ways. I can explore within sub differences between their matrices and
% will also be able to see how well each matrix corresponds across
% individuals.

% Average voxel time series data from the repeated runs of the same kind to
% generate one dataset for each condition. I think this happens after first
% levels. So I'm actually hyperaligning voxels that I extract from contrast
% images, not time series data. So now that I've found reliable voxels as a
% function of their relatedness to the task, I'll use the extracted
% coordinates to pick out the voxels I modeled in first levels. For future
% folks, I know this is circular. I'll extract networks from resting state
% data in the future. 

gain_consump_data_fnames = filenames(fullfile(con_dir,'consumption/*','ses-2/run-1/MID/con_0002.nii'));
loss_consump_data_fnames = filenames(fullfile(con_dir,'consumption/*','ses-2/run-1/MID/con_0001.nii'));
gain_anticip_data_fnames = filenames(fullfile(con_dir,'anticipation/*','ses-2/run-1/MID/con_0002.nii'));
loss_anticip_data_fnames = filenames(fullfile(con_dir,'anticipation/*','ses-2/run-1/MID/con_0001.nii'));

clear sub_data
for sub = 1:size(corr_input,3)
    sub_data = fmri_data(loss_consump_data_fnames{sub});
    ha_input(sub,:,1) = sub_data.dat(final_brain_coordinates_total{sub})';
end
clear sub_data
for sub = 1:size(corr_input,3)
    sub_data = fmri_data(gain_consump_data_fnames{sub});
    ha_input(sub,:,2) = sub_data.dat(final_brain_coordinates_total{sub})';
end
clear sub_data
for sub = 1:size(corr_input,3)
    sub_data = fmri_data(loss_anticip_data_fnames{sub});
    ha_input(sub,:,3) = sub_data.dat(final_brain_coordinates_total{sub})';
end
clear sub_data
for sub = 1:size(corr_input,3)
    sub_data = fmri_data(gain_anticip_data_fnames{sub});
    ha_input(sub,:,4) = sub_data.dat(final_brain_coordinates_total{sub})';
end

%% hyperalign subjects
[aligned, transforms] = hyperalign(ha_input(:,:,1),ha_input(:,:,2),ha_input(:,:,3),ha_input(:,:,4));

%% Visualize corr matrices across subs 
% This is my eyeballing estimate of increases in reliability
% side by side for loss consumption
figure();
subplot(1,2,1)
heatmap(corrcoef(ha_input(:,:,1)'))
title('Unaligned subjects loss consumption')
subplot(1,2,2)
heatmap(corrcoef(aligned{1}'))
title('Aligned subjects loss consumption')

%% side by side gain consumption
figure();
subplot(1,2,1)
heatmap(corrcoef(ha_input(:,:,2)'))
title('Unaligned subjects gain consumption')
subplot(1,2,2)
heatmap(corrcoef(aligned{2}'))
title('Aligned subjects gain consumption')

%% side by side loss anticipation
figure();
subplot(1,2,1)
heatmap(corrcoef(ha_input(:,:,3)'))
title('Unaligned subjects loss anticipation')
subplot(1,2,2)
heatmap(corrcoef(aligned{3}'))
title('Aligned subjects loss anticipation')

%% side by side gain consumption
figure();
subplot(1,2,1)
heatmap(corrcoef(ha_input(:,:,4)'))
title('Unaligned subjects gain anticipation')
subplot(1,2,2)
heatmap(corrcoef(aligned{4}'))
title('Aligned subjects gain anticipation')

%% To end off my project for psychometrics...
% I really just need to try and predict Trilevel from brain and be done
% with it. It's not publishable content whatever the result (or maybe it
% is?) but I need to knock it out so I can get to writing

clinicaldir = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD/clinical_data';

clear sub
 
load(fullfile(clinicaldir,'trilevel_factors.mat'));
for sub = 1:length(full_data_fnames)
    subid = full_data_fnames{sub}(105:109);
    PID(sub,1) = str2num(subid);
    if isempty(find(trilevel.ID(:) == str2num(subid))) == 0
        curr = find(trilevel.ID(:) == str2num(subid));
        curr_analysis_table(sub,:) = trilevel(curr,:);
    else
        disp(strcat(subid{sub}, ' missing clinical info'))       
    end
end

% Now we'll predict trilevel score from brain for each contrast separately
clear contrast 
for contrast = 1:length(aligned)
    regressors = aligned{contrast};
    % let's loop through all the trilevel symptoms to see if my little data
    % set predicts any of them
    for symp = 1:size(curr_analysis_table,2)
        curr_symp_name = curr_analysis_table.Properties.VariableNames{symp};
        regression_results.(curr_symp_name) = fitlm(curr_analysis_table.(curr_symp_name),regressors);
    end
end
                
        
%% Machine learning data prep 
% So it's visually apparent that there were some major gains in between sub
% correspondance between subjects. But now, I want to relate that
% improvement to something a little more tangible/sensitive. That means
% some kind of MVPA and writing this code now (in the context of a class so
% Robin can't be mad) is the move. So first, I need to go back to the time
% course data. I need separate matrices of each event and not just the
% average contrast that I used to show general improvements above

% Create array with the first field being a key that tells you what
% event the corresponding image in the second field of the array is
% referring to

% This first run will be a control condition where we try to classify each
% study event based on non hyperaligned voxel time courses. Then we'll need
% to hyperalign... I'll sort that out in the next section.
clear sub
clear event
clear onset
for sub = 1:length(full_data_fnames) 
    subid = full_data_fnames{sub}(105:109);
    % load in a subject brain data
    sub_data = fmri_data(full_data_fnames{sub}); 
    sub_data.dat = zscore(sub_data.dat);
    % loop through consump and anticip events. Each sub will have a
    % different number of these
    for event = 1:length(names)
        rsa_labels(sub,event) = names(event);
        % Since there's a different number of events, there's a different
        % number of onsets for each sub. So we have to store in a cell
        % array
        for onset = 1:length(onsets{event})
            TR = round(onsets{event}(onset) / 2.05);
            % Brian brought up a good point about trying to model that last
            % trial. I've done the best I could there. It makes a strong
            % case for concatenating runs actually... I should think about
            % that. But I think the scanner stops in between runs... So
            % actually I think we're shit outta luck there
            if TR < 279
                rsa_data{sub,event}(:,onset) = mean(sub_data.dat(final_brain_coordinates_total{sub},TR:TR+3),2);
            else
                rsa_data{sub,event}(:,onset) = mean(sub_data.dat(final_brain_coordinates_total{sub},TR:TR+2),2);
            end
        end
    end
end

% 

