%% Started as an Atrik script, now it's going to form the backbone of my project probs
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

analysisdir = '/Volumes/ZachExternal/ACNlab/BrainMAPD/func_conn/PPI/ppi_fldir';
condir = 'anticipation';%/bi_OFC_to_wholebrain';
concode = 1; % 1 - anticipation, other - consumption

% if concode == 1
%     region_list = {'bi_Amyg_to_wholebrain','biOFC_to_wholebrain','L_VS_AntRew_to_wholebrain','L_VS_AntLoss_to_wholebrain','LHO_Accumbens_to_wholebrain','LOFC2_to_wholebrain','LOFC_to_wholebrain','R_VS_AntLoss_to_wholebrain','R_VS_AntRew_to_wholebrain','RHO_Accumbens_to_wholebrain','ROFC_to_wholebrain'};
% else
%     region_list = {'bi_Amyg_to_wholebrain','biOFC_to_wholebrain','L_VS_to_wholebrain','LHO_Accumbens_to_wholebrain','LOFC2_to_wholebrain','LOFC_to_wholebrain','R_VS_to_wholebrain','RHO_Accumbens_to_wholebrain','ROFC_to_wholebrain'};
% end
region_list = {'biOFC_to_wholebrain','ROFC_to_wholebrain','LOFC_to_wholebrain','LOFC2_to_wholebrain'};
for seed_region = 1:length(region_list)

    cd(fullfile(analysisdir,condir,region_list{seed_region}))

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

    load(fullfile(clinicaldir,'trilevel_longitudinal.mat'));
    load(fullfile(med_dir,'Medication_T2.mat'));

    trilevel_array = [trilevel_T1.ID,trilevel_T1.GenDis,trilevel_T1.Anhedonia,trilevel_T1.Fears,trilevel_T1.Narrow];
    med_array = [T2MedicationInventory3(:,1),T2MedicationInventory3(:,88),T2MedicationInventory3(:,89)];
    med_array = table2array(med_array);

    for sub = 1:length(curr_fnames_gain)
        PID(sub,1) = str2num(curr_fnames_gain{sub}(1:5));
    end        


    for sub = 1:length(curr_fnames_gain)

        if isempty(find(trilevel_T1.ID(:) == PID(sub,1))) == 0
            curr = find(trilevel_T1.ID(:) == PID(sub,1));
            trilevel_regressors(sub,:) = trilevel_array(curr,:);
        else
            % for OFC loss_ppi_fnames{sub}(100:104); for HO_VMPFC loss_ppi_fnames{sub}(105:109)
            disp(strcat(num2str(PID(sub,1)), ' missing clinical info')) 
            trilevel_regressors(sub,:) = 0;%NaN;
            trilevel_regressors(sub,1) = PID(sub,1);
        end
        if isempty(find(med_array(:,1) == PID(sub,1))) == 0
            curr2 = find(med_array(:,1) == PID(sub,1));
            med_regressors(sub,:) = [med_array(curr2,2),med_array(curr2,3)];
        else
            disp(strcat(num2str(PID(sub,1)), ' missing medication info'))
            med_regressors(sub,:) = [0,0];
        end

    end

    site_regressors = (PID < 20000);

    GenDis = trilevel_regressors(:,2);
    Anhedonia = trilevel_regressors(:,3);
    Fears = trilevel_regressors(:,4);
    Narrow = trilevel_regressors(:,5);


    % CHANGE THIS
    R = [GenDis,Anhedonia,Fears,med_regressors,site_regressors];
    intercept = ones(length(R),1);

    R = [R,intercept];

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
    temp_results_gain = regress(curr_dat_gain);%,'robust');
    temp_results_loss = regress(curr_dat_loss);%,'robust');

    results_struct.gain.(region_list{seed_region}) = threshold(temp_results_gain.t,.001,'unc','k',30);
    results_struct.loss.(region_list{seed_region}) = threshold(temp_results_loss.t,.001,'unc','k',30);

    %% 6. gain

    orthviews(results_struct.gain.(region_list{seed_region}))
    %% 6. loss
    orthviews(results_struct.loss.(region_list{seed_region}))

    %% 7. 
    % CHANGE THIS
    symptom_names = {'GenDis','Anhedonia','Fears'};

    for symptom = 1:length(symptom_names)
        r_gain.(region_list{seed_region}).(symptom_names{symptom}) = region(select_one_image(results_struct.gain.(region_list{seed_region}),symptom));
        r_loss.(region_list{seed_region}).(symptom_names{symptom}) = region(select_one_image(results_struct.loss.(region_list{seed_region}),symptom));
    end

    %% GenDis Anhedonia Fears r_gain r_loss
     region_list{seed_region}
     orthviews(r_gain.(region_list{seed_region}).Fears) %, 'colormap', 'regioncenters');
 
     table(r_gain.(region_list{seed_region}).Fears)

    %% 8. 
    roidir = '/Users/zaz3744/Documents/current_projects/ACNlab/masks/ROI_BrainMAPD_functional';
    condir = 'anticipation';
    bilateral_roi_fnames = filenames(fullfile(roidir,condir,'*VS*Oldham*.nii'));
    right_roi_fnames = filenames(fullfile(roidir,condir,'right','*VS*Oldham*.nii'));
    left_roi_fnames = filenames(fullfile(roidir,condir,'left','*VS*Oldham*.nii'));
    %region_name_list = {'Ng_Amyg','BA9BA46','bi_vs_sphere','Ng_Caudate','HO_Accumbens','HO_Amyg','HO_Caudate','HO_Pallidum','HO_Putamen','HO_vmPFC','Knutson_mPFC','Knutson_OFC','Ng_OFC','Oldham_loss_VS','Oldham_gain_VS','VS_sphere'};
    region_name_list = {'VS_loss','VS_reward'};
    for i = 1:length(left_roi_fnames)
        disp(region_name_list{i})
        roi_fname = filenames(fullfile(left_roi_fnames{i}));
        roi = fmri_data(roi_fname);
        roi_avg_gain(i,:) = extract_roi_averages(curr_dat_gain,roi);
        roi_avg_gain(i,:).title = region_name_list{i}; 
        roi_gain_table(:,i) = roi_avg_gain(i,:).dat;
        % Note that the structure that is created below will have two
        % layers. The first refers to the seed region for this analysis.
        % The second will refer to the end region in this seed to seed func
        % conn analysis
        gain_mdl.(region_list{seed_region}).(region_name_list{i}) = fitlm(curr_dat_gain.X(:,1:6),roi_avg_gain(i,:).dat);
        fitlm(curr_dat_gain.X(:,1:6),roi_avg_gain(i,:).dat)
        roi_avg_loss(i,:) = extract_roi_averages(curr_dat_loss,roi);
        roi_avg_loss(i,:).title = region_name_list{i};
        roi_loss_table(:,i) = roi_avg_loss(i,:).dat;
        loss_mdl.(region_list{seed_region}).(region_name_list{i}) = fitlm(curr_dat_loss.X(:,1:6),roi_avg_loss(i,:).dat);
        fitlm(curr_dat_gain.X(:,1:6),roi_avg_loss(i,:).dat)
        disp(strcat('Seed region is: ',region_list{seed_region}))
        disp(strcat('End region is: ', region_name_list{i}))
        %keyboard
    end
    param_struct.(region_list{seed_region}).gain = roi_avg_gain;
    param_struct.(region_list{seed_region}).loss = roi_avg_loss;
%     roi_gain_table = [curr_dat_gain.X,roi_gain_table]; roi_gain_table = array2table(roi_gain_table); 
%     roi_gain_table.Properties.VariableNames = {'GenDis','Anhedonia','Fears','Ng_Amyg','BA9BA46','bi_vs_sphere','Ng_Caudate','HO_Accumbens','HO_Amyg','HO_Caudate','HO_Pallidum','HO_Putamen','HO_vmPFC','Knutson_mPFC','Knutson_OFC','Ng_OFC','Oldham_loss_VS','Oldham_gain_VS','VS_sphere'};
%     roi_loss_table = [curr_dat_gain.X,roi_loss_table]; roi_loss_table = array2table(roi_loss_table); 
%     roi_loss_table.Properties.VariableNames = {'GenDis','Anhedonia','Fears','Ng_Amyg','BA9BA46','bi_vs_sphere','Ng_Caudate','HO_Accumbens','HO_Amyg','HO_Caudate','HO_Pallidum','HO_Putamen','HO_vmPFC','Knutson_mPFC','Knutson_OFC','Ng_OFC','Oldham_loss_VS','Oldham_gain_VS','VS_sphere'};

    %% SVM test

%     svm_in_X = [curr_dat_gain.dat',R_final_gain(:,2:7)];
%     svm_in_Y = R_final_gain(:,1);
% 
%     svm_mdl.(region_list{seed_region}) = fitrsvm(svm_in_X,svm_in_Y);
% 
%     % Cross validate the model generated above
% 
%     svm_mdl_cv.(region_list{seed_region}) = crossval(svm_mdl.(region_list{seed_region}),'KFold',10);
% 
%     % Once you have a cross validated model, predict responses based on that model
% 
%     yfit.(region_list{seed_region}) = kfoldPredict(svm_mdl_cv.(region_list{seed_region}));
% 
%     % RMSE
% 
%     RMSE.(region_list{seed_region}) = sqrt(mean((svm_in_Y - yfit.(region_list{seed_region})).^2));  
    

end


% helpful for stuff

X = [param_struct.biOFC_to_wholebrain.gain(2).dat,param_struct.biOFC_to_wholebrain.gain(1).dat,param_struct.ROFC_to_wholebrain.gain(2).dat,param_struct.ROFC_to_wholebrain.gain(1).dat,param_struct.LOFC_to_wholebrain.gain(2).dat,param_struct.LOFC_to_wholebrain.gain(1).dat,param_struct.LOFC2_to_wholebrain.gain(2).dat,param_struct.LOFC2_to_wholebrain.gain(1).dat];
X = [param_struct.biOFC_to_wholebrain.loss(2).dat,param_struct.biOFC_to_wholebrain.loss(1).dat,param_struct.ROFC_to_wholebrain.loss(2).dat,param_struct.ROFC_to_wholebrain.loss(1).dat,param_struct.LOFC_to_wholebrain.loss(2).dat,param_struct.LOFC_to_wholebrain.loss(1).dat,param_struct.LOFC2_to_wholebrain.loss(2).dat,param_struct.LOFC2_to_wholebrain.loss(1).dat];
