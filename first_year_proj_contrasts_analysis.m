% Description: 
% 1. find filenames associated with MID anticipation period for gain and
% loss
% 2. load filenames into a data obj
% 3. load clinical data (trilevel data)
% 4. set dat.X = clinical symptoms
% 5. regress the data object to test for linear effects of trilevel
% symptoms on reward/loss anticipation data
% 6. Visualize results
% 7. Create region objs
% 8. Region of interest analysis with all ROIs

% 1.

analysisdir = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD/func_conn/first_levels/first_level_output';
condir = 'consumption';%/bi_OFC_to_wholebrain';

cd(fullfile(analysisdir,condir))

D = dir; D(1:2,:) = [];
for sub = 1:length(D)
    curr_folder = D(sub).name;
    curr_fnames_gain{sub} = char(filenames(fullfile(curr_folder,'ses-2/run-1/MID/con_0001.nii')));
    curr_fnames_loss{sub} = char(filenames(fullfile(curr_folder,'ses-2/run-1/MID/con_0002.nii')));
end


%% 2.
curr_dat_gain = fmri_data(curr_fnames_gain');
curr_dat_loss = fmri_data(curr_fnames_loss');


%% 3. and 4. 
% Load in clinical data for group analysis
clinicaldir = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD/clinical_data';
med_dir = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD';
demo_dir = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD';
load(fullfile(clinicaldir,'first_year_project_trilevel_T1.mat'));
load(fullfile(med_dir,'Medication_T2.mat'));
load(fullfile(demo_dir,'demographics.mat'));

trilevel_array = [trilevel_T1.id,trilevel_T1.GenDis,trilevel_T1.Anhedon,trilevel_T1.Fears];
med_array = [T2MedicationInventory3(:,1),T2MedicationInventory3(:,88)];
med_array = table2array(med_array);
dem_array = [BrainMAPDT1S1Demo.PID,BrainMAPDT1S1Demo.sex];
%dem_array = table2array(dem_array);

for sub = 1:length(curr_fnames_gain)
    PID(sub,1) = str2num(curr_fnames_gain{sub}(1:5));
end        


for sub = 1:length(curr_fnames_gain)

    if isempty(find(trilevel_T1.id(:) == PID(sub,1))) == 0
        curr = find(trilevel_T1.id(:) == PID(sub,1));
        trilevel_regressors(sub,:) = trilevel_array(curr,:);
    else
        % for OFC loss_ppi_fnames{sub}(100:104); for HO_VMPFC loss_ppi_fnames{sub}(105:109)
        disp(strcat(num2str(PID(sub,1)), ' missing clinical info')) 
        trilevel_regressors(sub,:) = NaN;
        trilevel_regressors(sub,1) = PID(sub,1);
    end
    if isempty(find(med_array(:,1) == PID(sub,1))) == 0
        curr2 = find(med_array(:,1) == PID(sub,1));
        med_regressors(sub,:) = med_array(curr2,2);
    else
        disp(strcat(num2str(PID(sub,1)), ' missing medication info'))
        med_regressors(sub,:) = [0];
    end
    if isempty(find(dem_array(:,1) == PID(sub,1))) == 0
        curr2 = find(dem_array(:,1) == PID(sub,1));
        demographic_regressors(sub,:) = dem_array(curr2,2);
    else
        disp(strcat(num2str(PID(sub,1)), ' missing demographic info'))
        demographic_regressors(sub,:) = 0;
    end

end

site_regressors = (PID < 20000);

GenDis = trilevel_regressors(:,2);
Anhedonia = trilevel_regressors(:,3);
Fears = trilevel_regressors(:,4);

% CHANGE THIS
R = [GenDis,Anhedonia,Fears,med_regressors,site_regressors]; %demographic_regressors = sex
intercept = ones(length(R),1);

%R = [R,intercept];

R_gain_temp = [R,curr_dat_gain.dat'];
R_loss_temp = [R,curr_dat_loss.dat'];

R_gain_temp(any(isnan(R_gain_temp),2),:)=[];
R_loss_temp(any(isnan(R_loss_temp),2),:)=[];
PID(any(isnan(R_loss_temp),2),:)=[];

R_final_gain = R_gain_temp(:,1:size(R,2));
R_final_loss = R_loss_temp(:,1:size(R,2));

dat_final_gain = R_gain_temp(:,size(R,2)+1:size(R_gain_temp,2));
dat_final_loss = R_loss_temp(:,size(R,2)+1:size(R_loss_temp,2));

curr_dat_gain.X = R_final_gain;
curr_dat_loss.X = R_final_loss;

curr_dat_gain.dat = dat_final_gain';
curr_dat_loss.dat = dat_final_loss';

%% 5. 
temp_results_gain= regress(curr_dat_gain); %,'robust');
temp_results_loss= regress(curr_dat_loss); %,'robust');
results_struct.gain= threshold(temp_results_gain.t,.001,'unc','k',10);
results_struct.loss= threshold(temp_results_loss.t,.001,'unc','k',10);
%% 6. gain

% orthviews(results_struct.gain)
%% 6. loss
% orthviews(results_struct.loss)

%% 7. 
% CHANGE THIS
symptom_names = {'GenDis','Anhedonia','Fears'};

for symptom = 1:length(symptom_names)
    r_gain.(symptom_names{symptom}) = region(select_one_image(results_struct.gain,symptom));
    r_loss.(symptom_names{symptom}) = region(select_one_image(results_struct.loss,symptom));
end

% GenDis Anhedonia Fears r_gain r_loss


%% 8. 
roidir = '/Users/zaz3744/Documents/current_projects/ACNlab/masks/ROI_BrainMAPD_functional';
condir = 'anticipation';
bilateral_roi_fnames = filenames(fullfile(roidir,condir,'*VS*Oldham*.nii'));
right_roi_fnames = filenames(fullfile(roidir,condir,'right','*VS*Oldham*.nii'));
left_roi_fnames = filenames(fullfile(roidir,condir,'left','*VS*Oldham*.nii'));
region_name_list = {'OFC*Oldham*','Amyg*Ng*','BA9BA46*','bi_vs_sphere*','Caudate*Ng*','HO_Accumbens*','HO_Amyg*','HO_Caudate*','HO_Pallidum*','HO_Putamen*','HO_VMPFC*','*MPFC*Knutson*','OFC*Knutson*','OFC*Ng*','VS*Oldham_Loss*','VS*Oldham_Rew*','VS_Sphere*'};
region_name_list_for_struct = {'bi_OFC_Oldham','Amyg_Ng','BA9BA46','bi_vs_sphere','Caudate_Ng','HO_Accumbens','HO_Amyg','HO_Caudate','HO_Pallidum','HO_Putamen','HO_VMPFC','mPFC_Knutson','OFC_Knutson','OFC_Ng','VS_Oldham_Loss','VS_Oldham_Rew','VS_sphere'};


for i = 1:length(region_name_list_for_struct)
    disp(region_name_list_for_struct{i})
    roi_fname = filenames(fullfile(roidir,condir,region_name_list{i}));
    roi = fmri_data(roi_fname{1});
    roi_avg_gain.(region_name_list_for_struct{i}) = extract_roi_averages(curr_dat_gain,roi);
    roi_avg_gain.(region_name_list_for_struct{i}).title = region_name_list_for_struct{i}; 
    roi_gain_table(:,i) = roi_avg_gain.(region_name_list_for_struct{i}).dat;
    % Note that the structure that is created below will have two
    % layers. The first refers to the seed region for this analysis.
    % The second will refer to the end region in this seed to seed func
    % conn analysis
    %gain_mdl.(region_name_list_for_struct{i}) = fitlm(curr_dat_gain.X(:,1:5),roi_avg_gain(i,:).dat);
    %fitlm(curr_dat_gain.X(:,1:5),roi_avg_gain(i,:).dat)
    gain_mdl.(region_name_list_for_struct{i}) = fitlm(R_final_gain,roi_avg_gain.(region_name_list_for_struct{i}).dat);
    fitlm(R_final_gain,roi_avg_gain.(region_name_list_for_struct{i}).dat)
    roi_avg_loss.(region_name_list_for_struct{i}) = extract_roi_averages(curr_dat_loss,roi);
    roi_avg_loss.(region_name_list_for_struct{i}).title = region_name_list_for_struct{i};
    roi_loss_table(:,i) = roi_avg_loss.(region_name_list_for_struct{i}).dat;
    %loss_mdl.(region_name_list_for_struct{i}) = fitlm(curr_dat_loss.X(:,1:5),roi_avg_loss(i,:).dat);
    %fitlm(curr_dat_gain.X(:,1:5),roi_avg_loss(i,:).dat)
    loss_mdl.(region_name_list_for_struct{i}) = fitlm(R_final_loss,roi_avg_loss.(region_name_list_for_struct{i}).dat);
    fitlm(R_final_loss,roi_avg_loss.(region_name_list_for_struct{i}).dat)
    disp(strcat('ROI: ',region_name_list_for_struct{i}))
    
    %keyboard
end


R_region_name_list = {'R_OFC*Oldham*','R_VS_Sphere*','*Caudate*Ng*','HO*_Accumbens*','HO*_Amyg*','HO*_Caudate*','HO*_Pallidum*','HO*_Putamen*','*OFC*Ng*','*VS*Oldham_Loss*','*VS*Oldham_Rew*'};
R_region_name_list_for_struct = {'R_OFC_Oldham','R_VS_sphere','R_Caudate_Ng','R_HO_Accumbens','R_HO_Amyg','R_HO_Caudate','R_HO_Pallidum','R_HO_Putamen','R_OFC_Ng','R_VS_Oldham_Loss','R_VS_Oldham_Rew'};

% Right lateralized
for i = 1:length(R_region_name_list_for_struct)
    disp(R_region_name_list_for_struct{i})
    roi_fname = filenames(fullfile(roidir,condir,'right',R_region_name_list{i}));
    roi = fmri_data(roi_fname{1});
    roi_avg_gain.(R_region_name_list_for_struct{i}) = extract_roi_averages(curr_dat_gain,roi);
    roi_avg_gain.(R_region_name_list_for_struct{i}).title = R_region_name_list_for_struct{i}; 
    roi_gain_table(:,i) = roi_avg_gain.(R_region_name_list_for_struct{i}).dat;
    % Note that the structure that is created below will have two
    % layers. The first refers to the seed region for this analysis.
    % The second will refer to the end region in this seed to seed func
    % conn analysis
    %gain_mdl.(region_name_list_for_struct{i}) = fitlm(curr_dat_gain.X(:,1:5),roi_avg_gain(i,:).dat);
    %fitlm(curr_dat_gain.X(:,1:5),roi_avg_gain(i,:).dat)
    gain_mdl.(R_region_name_list_for_struct{i}) = fitlm(R_final_gain,roi_avg_gain.(R_region_name_list_for_struct{i}).dat);
    fitlm(R_final_gain,roi_avg_gain.(R_region_name_list_for_struct{i}).dat)
    roi_avg_loss.(R_region_name_list_for_struct{i}) = extract_roi_averages(curr_dat_loss,roi);
    roi_avg_loss.(R_region_name_list_for_struct{i}).title = R_region_name_list_for_struct{i};
    roi_loss_table(:,i) = roi_avg_loss.(R_region_name_list_for_struct{i}).dat;
    %loss_mdl.(region_name_list_for_struct{i}) = fitlm(curr_dat_loss.X(:,1:5),roi_avg_loss(i,:).dat);
    %fitlm(curr_dat_gain.X(:,1:5),roi_avg_loss(i,:).dat)
    loss_mdl.(R_region_name_list_for_struct{i}) = fitlm(R_final_loss,roi_avg_loss.(R_region_name_list_for_struct{i}).dat);
    fitlm(R_final_loss,roi_avg_loss.(R_region_name_list_for_struct{i}).dat)
    disp(strcat('ROI: ',R_region_name_list_for_struct{i}))
    
    %keyboard
end

L_region_name_list = {'L_OFC_*Oldham.nii','L_OFC_*Oldham2.nii','L*Amyg*Ng.nii','L*Amyg*Ng_2*','L_VS_Sphere*','*Caudate*Ng*','HO*_Accumbens*','HO*_Amyg*','HO*_Caudate*','HO*_Pallidum*','HO*_Putamen*','*OFC*Ng*','*VS*Oldham_Loss*','*VS*Oldham_Rew*'};
L_region_name_list_for_struct = {'L_OFC1_Oldham','L_OFC2_Oldham','L_Amyg_Ng','L_Amyg_Ng_2','L_VS_Sphere','L_Caudate_Ng','L_HO_Accumbens','L_HO_Amyg','L_HO_Caudate','L_HO_Pallidum','L_HO_Putamen','L_OFC_Ng','L_VS_Oldham_Loss','L_VS_Oldham_Rew'};

% left lateralized

for i = 1:length(L_region_name_list_for_struct)
    disp(L_region_name_list_for_struct{i})
    roi_fname = filenames(fullfile(roidir,condir,'left',L_region_name_list{i}));
    roi = fmri_data(roi_fname{1});
    roi_avg_gain.(L_region_name_list_for_struct{i}) = extract_roi_averages(curr_dat_gain,roi);
    roi_avg_gain.(L_region_name_list_for_struct{i}).title = L_region_name_list_for_struct{i}; 
    roi_gain_table(:,i) = roi_avg_gain.(L_region_name_list_for_struct{i}).dat;
    % Note that the structure that is created below will have two
    % layers. The first refers to the seed region for this analysis.
    % The second will refer to the end region in this seed to seed func
    % conn analysis
    %gain_mdl.(region_name_list_for_struct{i}) = fitlm(curr_dat_gain.X(:,1:5),roi_avg_gain(i,:).dat);
    %fitlm(curr_dat_gain.X(:,1:5),roi_avg_gain(i,:).dat)
    gain_mdl.(L_region_name_list_for_struct{i}) = fitlm(R_final_gain,roi_avg_gain.(L_region_name_list_for_struct{i}).dat);
    fitlm(R_final_gain,roi_avg_gain.(L_region_name_list_for_struct{i}).dat)
    roi_avg_loss.(L_region_name_list_for_struct{i}) = extract_roi_averages(curr_dat_loss,roi);
    roi_avg_loss.(L_region_name_list_for_struct{i}).title = L_region_name_list_for_struct{i};
    roi_loss_table(:,i) = roi_avg_loss.(L_region_name_list_for_struct{i}).dat;
    %loss_mdl.(region_name_list_for_struct{i}) = fitlm(curr_dat_loss.X(:,1:5),roi_avg_loss(i,:).dat);
    %fitlm(curr_dat_gain.X(:,1:5),roi_avg_loss(i,:).dat)
    loss_mdl.(L_region_name_list_for_struct{i}) = fitlm(R_final_loss,roi_avg_loss.(L_region_name_list_for_struct{i}).dat);
    fitlm(R_final_loss,roi_avg_loss.(L_region_name_list_for_struct{i}).dat)
    disp(strcat('ROI: ',L_region_name_list_for_struct{i}))
    
    %keyboard
end


%% Same analyses as above but controlling for loss/gain condition
% bilateral
for i = 1:length(region_name_list_for_struct)
    disp(region_name_list_for_struct{i})
    roi_fname = filenames(fullfile(roidir,condir,region_name_list{i}));
    roi = fmri_data(roi_fname{1});

    % Note that the structure that is created below will have two
    % layers. The first refers to the seed region for this analysis.
    % The second will refer to the end region in this seed to seed func
    % conn analysis

    % extract roi data for gain and loss
    roi_avg_gain2.(region_name_list_for_struct{i}) = extract_roi_averages(curr_dat_gain,roi);
    roi_avg_gain2.(region_name_list_for_struct{i}).title = region_name_list_for_struct{i}; 
    roi_avg_loss2.(region_name_list_for_struct{i}) = extract_roi_averages(curr_dat_loss,roi);
    roi_avg_loss2.(region_name_list_for_struct{i}).title = region_name_list_for_struct{i};

    % add the loss condition to the gain model to control for
    % activation during loss anticipation
    R_final_gain2 = [R_final_gain,roi_avg_loss2.(region_name_list_for_struct{i}).dat];
    gain_mdl2.(region_name_list_for_struct{i}) = fitlm(R_final_gain2,roi_avg_gain2.(region_name_list_for_struct{i}).dat);
    fitlm(R_final_gain2,roi_avg_gain2.(region_name_list_for_struct{i}).dat)

    % add the gain condition to the loss model to control for
    % activation during gain anticipation
    R_final_loss2 = [R_final_loss,roi_avg_gain.(region_name_list_for_struct{i}).dat];
    loss_mdl2.(region_name_list_for_struct{i}) = fitlm(R_final_loss2,roi_avg_loss2.(region_name_list_for_struct{i}).dat);
    fitlm(R_final_loss2,roi_avg_loss2.(region_name_list_for_struct{i}).dat)
    disp(strcat('ROI: ',region_name_list_for_struct{i}))
    %keyboard
end

% right lateralized
for i = 1:length(R_region_name_list_for_struct)
    disp(R_region_name_list_for_struct{i})
    roi_fname = filenames(fullfile(roidir,condir,'right',R_region_name_list{i}));
    roi = fmri_data(roi_fname{1});

    % Note that the structure that is created below will have two
    % layers. The first refers to the seed region for this analysis.
    % The second will refer to the end region in this seed to seed func
    % conn analysis

    % extract roi data for gain and loss
    roi_avg_gain2.(R_region_name_list_for_struct{i}) = extract_roi_averages(curr_dat_gain,roi);
    roi_avg_gain2.(R_region_name_list_for_struct{i}).title = R_region_name_list_for_struct{i}; 
    roi_avg_loss2.(R_region_name_list_for_struct{i}) = extract_roi_averages(curr_dat_loss,roi);
    roi_avg_loss2.(R_region_name_list_for_struct{i}).title = R_region_name_list_for_struct{i};

    % add the loss condition to the gain model to control for
    % activation during loss anticipation
    R_final_gain2 = [R_final_gain,roi_avg_loss2.(R_region_name_list_for_struct{i}).dat];
    gain_mdl2.(R_region_name_list_for_struct{i}) = fitlm(R_final_gain2,roi_avg_gain2.(R_region_name_list_for_struct{i}).dat);
    fitlm(R_final_gain2,roi_avg_gain2.(R_region_name_list_for_struct{i}).dat)

    % add the gain condition to the loss model to control for
    % activation during gain anticipation
    R_final_loss2 = [R_final_loss,roi_avg_gain.(R_region_name_list_for_struct{i}).dat];
    loss_mdl2.(R_region_name_list_for_struct{i}) = fitlm(R_final_loss2,roi_avg_loss2.(R_region_name_list_for_struct{i}).dat);
    fitlm(R_final_loss2,roi_avg_loss2.(R_region_name_list_for_struct{i}).dat)
    disp(strcat('ROI: ',R_region_name_list_for_struct{i}))
    %keyboard
end

% left lateralized
for i = 1:length(L_region_name_list_for_struct)
    disp(L_region_name_list_for_struct{i})
    roi_fname = filenames(fullfile(roidir,condir,'left',L_region_name_list{i}));
    roi = fmri_data(roi_fname{1});

    % Note that the structure that is created below will have two
    % layers. The first refers to the seed region for this analysis.
    % The second will refer to the end region in this seed to seed func
    % conn analysis

    % extract roi data for gain and loss
    roi_avg_gain2.(L_region_name_list_for_struct{i}) = extract_roi_averages(curr_dat_gain,roi);
    roi_avg_gain2.(L_region_name_list_for_struct{i}).title = L_region_name_list_for_struct{i}; 
    roi_avg_loss2.(L_region_name_list_for_struct{i}) = extract_roi_averages(curr_dat_loss,roi);
    roi_avg_loss2.(L_region_name_list_for_struct{i}).title = L_region_name_list_for_struct{i};

    % add the loss condition to the gain model to control for
    % activation during loss anticipation
    R_final_gain2 = [R_final_gain,roi_avg_loss2.(L_region_name_list_for_struct{i}).dat];
    gain_mdl2.(L_region_name_list_for_struct{i}) = fitlm(R_final_gain2,roi_avg_gain2.(L_region_name_list_for_struct{i}).dat);
    fitlm(R_final_gain2,roi_avg_gain2.(L_region_name_list_for_struct{i}).dat)

    % add the gain condition to the loss model to control for
    % activation during gain anticipation
    R_final_loss2 = [R_final_loss,roi_avg_gain.(L_region_name_list_for_struct{i}).dat];
    loss_mdl2.(L_region_name_list_for_struct{i}) = fitlm(R_final_loss2,roi_avg_loss2.(L_region_name_list_for_struct{i}).dat);
    fitlm(R_final_loss2,roi_avg_loss2.(L_region_name_list_for_struct{i}).dat)
    disp(strcat('ROI: ',L_region_name_list_for_struct{i}))
    %keyboard
end


