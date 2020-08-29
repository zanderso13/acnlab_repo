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
condir = 'consumption';%/bi_OFC_to_wholebrain';
region_list = {'bi_OFC*','L_VS_AntRew*','L_VS_AntLoss*'


for seed_region = 1:length(region_list)

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
            trilevel_regressors(sub,:) = NaN;
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

    R_final_gain = R_gain_temp(:,1:size(R,2));
    R_final_loss = R_loss_temp(:,1:size(R,2));

    dat_final_gain = R_gain_temp(:,size(R,2)+1:size(R_gain_temp,2));
    dat_final_loss = R_loss_temp(:,size(R,2)+1:size(R_loss_temp,2));

    curr_dat_gain.X = R_final_gain;
    curr_dat_loss.X = R_final_loss;

    curr_dat_gain.dat = dat_final_gain';
    curr_dat_loss.dat = dat_final_loss';


    %% 5. 
    temp_results_gain = regress(curr_dat_gain,'robust');
    temp_results_loss = regress(curr_dat_loss,'robust');

    results_struct.gain = threshold(temp_results_gain.t,.001,'unc','k',30);
    results_struct.loss = threshold(temp_results_loss.t,.001,'unc','k',30);

    %% 6. gain

    orthviews(results_struct.gain)
    %% 6. loss
    orthviews(results_struct.loss)

    %% 7. 
    % CHANGE THIS
    symptom_names = {'GenDis','Anhedonia','Fears'};

    for symptom = 1:length(symptom_names)
        r_gain.(symptom_names{symptom}) = region(select_one_image(results_struct.gain,symptom));
        r_loss.(symptom_names{symptom}) = region(select_one_image(results_struct.loss,symptom));
    end

    %% GenDis Anhedonia Fears r_gain r_loss
    montage(r_loss.Fears) %, 'colormap', 'regioncenters');

    table(r_loss.Fears)

    %% 8. 
    roidir = '/Users/zaz3744/Documents/current_projects/ACNlab/masks/ROI_BrainMAPD_functional';
    condir = 'anticipation'
    bilateral_roi_fnames = filenames(fullfile(roidir,condir,'*.nii'));
    right_roi_fnames = filenames(fullfile(roidir,condir,'*.nii'));
    left_roi_fnames = filenames(fullfile(roidir,condir,'*.nii'));
    region_name_list = {'Ng_Amyg','BA9BA46','bi_vs_sphere','Ng_Caudate','HO_Accumbens','HO_Amyg','HO_Caudate','HO_Pallidum','HO_Putamen','HO_vmPFC','Knutson_mPFC','Knutson_OFC','Ng_OFC','Oldham_loss_VS','Oldham_gain_VS','VS_sphere'};
    for region = 1:length(bilateral_roi_fnames)
        disp(region_name_list{region})
        roi_fname = filenames(fullfile(bilateral_roi_fnames{region}));
        roi = fmri_data(roi_fname);
        roi_avg_gain(region,:) = extract_roi_averages(curr_dat_gain,roi);
        roi_avg_gain(region,:).title = region_name_list{region}; 
        roi_gain_table(:,region) = roi_avg_gain(region,:).dat;
        gain_mdl.(region_name_list{region}) = fitlm(curr_dat_gain.X,roi_avg_gain(region,:).dat)
        roi_avg_loss(region,:) = extract_roi_averages(curr_dat_loss,roi);
        roi_avg_loss(region,:).title = region_name_list{region};
        roi_loss_table(:,region) = roi_avg_loss(region,:).dat;
        loss_mdl.(region_name_list{region}) = fitlm(curr_dat_gain.X,roi_avg_gain(region,:).dat)
    end

    % roi_gain_table = [curr_dat_gain.X,roi_gain_table]; roi_gain_table = array2table(roi_gain_table); 
    % roi_gain_table.Properties.VariableNames = {'GenDis','Anhedonia','Fears','Ng_Amyg','BA9BA46','bi_vs_sphere','Ng_Caudate','HO_Accumbens','HO_Amyg','HO_Caudate','HO_Pallidum','HO_Putamen','HO_vmPFC','Knutson_mPFC','Knutson_OFC','Ng_OFC','Oldham_loss_VS','Oldham_gain_VS','VS_sphere'};
    % roi_loss_table = [curr_dat_gain.X,roi_loss_table]; roi_loss_table = array2table(roi_loss_table); 
    % roi_loss_table.Properties.VariableNames = {'GenDis','Anhedonia','Fears','Ng_Amyg','BA9BA46','bi_vs_sphere','Ng_Caudate','HO_Accumbens','HO_Amyg','HO_Caudate','HO_Pallidum','HO_Putamen','HO_vmPFC','Knutson_mPFC','Knutson_OFC','Ng_OFC','Oldham_loss_VS','Oldham_gain_VS','VS_sphere'};
end





