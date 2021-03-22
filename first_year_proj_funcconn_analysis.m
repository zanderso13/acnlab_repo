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
% What kind of analysis do you want to run? Are you interested in trilevel
% symptom dimensions or dsm diagnostic categories? Only ONE of the below
% two options can be set to 1
trilevel = 1;
dsm = 0;
life_stress = 0;

analysisdir = '/Volumes/ZachExternal/ACNlab/BrainMAPD/func_conn/PPI/ppi_fldir';
condir = 'consumption';
% Generate a model to plot with the scatter
savedir = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD/func_conn/PPI/results/Trilevel_analysis/figures';
constring = 'Consumption';

% OPTIONS FOR MODELS TO RUN
% symptom_names = {'GenDis','Anhedonia','Fears'};
% symptom_names = {'Anhedonia'};
% symptom_names = {'GenDis'};
% symptom_names = {'Fears'};
% symptom_names = {'StressxAnhedonia','Anhedonia','Stress','GenDis','Fears'};
% symptom_names = {'StressxAnhedonia','Anhedonia','Stress'};
% symptom_names = {'FearsxAnhedonia','Anhedonia','GenDis','Fears'};
% symptom_names = {'GenDisxFearsxAnhedonia','FearsxAnhedonia','Anhedonia','GenDisxFears','GenDisxAnhedonia','Fears','GenDis'};

if trilevel == 1
    symptom_names = {'GenDis','Anhedonia','Fears'};
elseif dsm == 1
    symptom_names = {'Comorbid','Depression','Anxiety'};
end

% concode = 1; % 1 - anticipation, other - consumption
% 
% if concode == 1
%     region_list = {'bi_VS_AntRew','bi_VS_AntLoss','bi_vmPFC_to_wholebrain','bi_Amyg_to_wholebrain','biOFC_to_wholebrain','L_VS_AntRew_to_wholebrain','L_VS_AntLoss_to_wholebrain','LHO_Accumbens_to_wholebrain','LOFC2_to_wholebrain','LOFC_to_wholebrain','R_VS_AntLoss_to_wholebrain','R_VS_AntRew_to_wholebrain','RHO_Accumbens_to_wholebrain','ROFC_to_wholebrain'};
% else
%     region_list = {'bi_vmPFC_to_wholebrain','bi_Amyg_to_wholebrain','biOFC_to_wholebrain','L_VS_to_wholebrain','LHO_Accumbens_to_wholebrain','LOFC2_to_wholebrain','LOFC_to_wholebrain','R_VS_to_wholebrain','RHO_Accumbens_to_wholebrain','ROFC_to_wholebrain'};
% end
if strcmp(condir, 'anticipation')
    % anticipation list for first year project region_list = {'biOFC_to_wholebrain','bi_VS_AntRew','bi_VS_AntLoss','ROFC_to_wholebrain','LOFC2_to_wholebrain','LOFC_to_wholebrain','L_VS_AntRew_to_wholebrain','L_VS_AntLoss_to_wholebrain','R_VS_AntLoss_to_wholebrain','R_VS_AntRew_to_wholebrain'};
    region_list = {'biOFC_to_wholebrain','bi_VS_AntLoss','bi_VS_AntRew','bi_Amyg_to_wholebrain','ROFC_to_wholebrain','LOFC_to_wholebrain','LOFC2_to_wholebrain'};%,'LOFC2_to_wholebrain','LOFC_to_wholebrain','R_VS_AntRew_to_wholebrain','L_VS_AntRew_to_wholebrain'};
    % region_list = {'biOFC_to_wholebrain','bi_VS_AntLoss','bi_VS_AntRew','HOAmygdala','HOAccumbens','HOvmPFC'};
else
    % consumption
    region_list = {'biOFC_to_wholebrain','bi_VS_Oldham','bi_Amyg_to_wholebrain','ROFC_to_wholebrain','LOFC_to_wholebrain','LOFC2_to_wholebrain'};%{'biOFC_to_wholebrain','bi_VS_Oldham','ROFC_to_wholebrain','LOFC2_to_wholebrain','LOFC_to_wholebrain'};%,'R_VS_to_wholebrain','L_VS_to_wholebrain'};
    %region_list = {'biOFC_to_wholebrain','bi_VS_Oldham','HOAmygdala','HOAccumbens','HOvmPFC'};
    
end


for seed_region = 1:length(region_list)

    cd(fullfile(analysisdir,condir,region_list{seed_region}))

    D = dir; D(1:2,:) = [];
    for sub = 1:length(D)
        curr_folder = D(sub).name;
        curr_fnames_gain{sub} = char(filenames(fullfile(curr_folder,'ses-2/run-1/MID/con_0001.nii')));
        curr_fnames_loss{sub} = char(filenames(fullfile(curr_folder,'ses-2/run-1/MID/con_0002.nii')));
    end


    %% 2.
    curr_dat_gain.(region_list{seed_region}) = fmri_data(curr_fnames_gain');
    curr_dat_loss.(region_list{seed_region}) = fmri_data(curr_fnames_loss');



    %% 3. and 4. 
    % Load in clinical data for group analysis
    clinicaldir = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD/clinical_data';
    med_dir = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD';
    demo_dir = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD';
    load(fullfile(clinicaldir,'first_year_project_trilevel_T1.mat'));
    load(fullfile(med_dir,'Medication_T2.mat'));
    load(fullfile(demo_dir,'demographics.mat'));
    
    %% load in clinical categorical diagnosis
    % I have all four timepoints and I'd like to see who of our healthies
    % develops psychopathology. In the next two sections, I'm checking for
    % how many participants have current clinical diagnoses

    clinicaldir = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD/clinical_data';
    load(fullfile('/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD/PID.mat'))
    load(fullfile(clinicaldir,'first_year_project_BrainMAPD_clinical_diagnoses_final.mat'))
    
    stressdir = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD/life_stress';
    load(fullfile(stressdir,'LSI_T1.mat'))
    

    %% calculate number currently in episode
    for sub = 1:length(PID)
        curr_sub = PID(sub);
        PID(sub) = curr_sub;
        if isempty(clinical_info.dep_life_any(curr_sub==clinical_info.PID(:)))==0
            dep1(sub) = [clinical_info.dep_curr_any(curr_sub==clinical_info.PID(:))];
            anx1(sub) = [clinical_info.anx_curr_any(curr_sub==clinical_info.PID(:))];
            com1(sub) = [clinical_info.comorbid_life_dep_anx(curr_sub==clinical_info.PID(:))];
        else
            disp(curr_sub)

        end

    end 

    %%  code below takes different kinds of regressors into account
    
    for sub = 1:length(PID)
        curr_sub = PID(sub);
        PID(sub) = curr_sub;
        if isempty(clinical_info.dep_life_any(curr_sub==clinical_info.PID(:)))==0
            deplife2(sub) = [clinical_info.dep_life_any(curr_sub==clinical_info.PID(:))];
            anxlife2(sub) = [clinical_info.anx_life_any(curr_sub==clinical_info.PID(:))];
            comlife2(sub) = [clinical_info.comorbid_life_dep_anx(curr_sub==clinical_info.PID(:))];
        
        else
            disp(curr_sub)


        end

    end

    
    
    trilevel_array = [trilevel_T1.id,trilevel_T1.GenDis,trilevel_T1.Anhedon,trilevel_T1.Fears,trilevel_T1.narrow];
    %trilevel_array = trilevel_array(:,2:size(trilevel_array,2));
    med_array = [T2MedicationInventory3(:,1),T2MedicationInventory3(:,88)];
    med_array = table2array(med_array);
    dem_array = [BrainMAPDT1S1Demo.PID,BrainMAPDT1S1Demo.sex,BrainMAPDT1S1Demo.ethnicity,BrainMAPDT1S1Demo.race,BrainMAPDT1S1Demo.race2];
    lsi_array = [LSI_T1(:,1),LSI_T1(:,6),LSI_T1(:,8:9),LSI_T1(:,11:13),LSI_T1(:,15),LSI_T1(:,17:19)];
    lsi_array = table2array(lsi_array);
    lsi_array(any(isnan(lsi_array),2),:)=[];
    %dem_array = table2array(dem_array);
    clear sub PID
    for sub = 1:length(curr_fnames_gain)
        PID(sub,1) = str2num(curr_fnames_gain{sub}(1:5));
    end        


    for sub = 1:length(curr_fnames_gain)
        % trilevel 
        if isempty(find(trilevel_T1.id(:) == PID(sub,1))) == 0
            curr = find(trilevel_T1.id(:) == PID(sub,1));
            trilevel_regressors(sub,:) = trilevel_array(curr,:);
        else
            % for OFC loss_ppi_fnames{sub}(100:104); for HO_VMPFC loss_ppi_fnames{sub}(105:109)
            disp(strcat(num2str(PID(sub,1)), ' missing clinical info')) 
            trilevel_regressors(sub,:) = NaN;
            trilevel_regressors(sub,1) = PID(sub,1);
        end
        % medication
        if isempty(find(med_array(:,1) == PID(sub,1))) == 0
            curr2 = find(med_array(:,1) == PID(sub,1));
            med_regressors(sub,:) = med_array(curr2,2);
        else
            disp(strcat(num2str(PID(sub,1)), ' missing medication info'))
            med_regressors(sub,:) = [0];
        end
        % demographics
        if isempty(find(dem_array(:,1) == PID(sub,1))) == 0
            curr2 = find(dem_array(:,1) == PID(sub,1));
            demographic_regressors(sub,:) = dem_array(curr2,2);%:5);
        else
            disp(strcat(num2str(PID(sub,1)), ' missing demographic info'))
            demographic_regressors(sub,:) = 0;
        end
        % life stress
        if isempty(find(lsi_array(:,1) == PID(sub,1))) == 0
            curr2 = find(lsi_array(:,1) == PID(sub,1));
            lsi_regressors(sub,:) = sum(lsi_array(curr2,2:11));
        else
            disp(strcat(num2str(PID(sub,1)), ' missing life stress info'))
            lsi_regressors(sub,1) = NaN;
        end
    end

    
    site_regressors = (PID < 20000);

    GenDis = trilevel_regressors(:,2);
    Anhedonia = trilevel_regressors(:,3);
    Fears = trilevel_regressors(:,4);
    Narrow = trilevel_regressors(:,5);
    

    % CHANGE THIS
    if trilevel == 1
        % Will also run models separately for each diagnosis. 
        if length(symptom_names) == 5
            R = [Anhedonia.*lsi_regressors,Anhedonia,lsi_regressors,GenDis,Fears,Narrow,med_regressors,site_regressors,demographic_regressors]; %demographic_regressors = sex
        elseif strcmp(symptom_names{1},'GenDis')
            R = [GenDis,Anhedonia,Fears,med_regressors, site_regressors,demographic_regressors]; %Anhedonia,Fears,,med_regressors
        elseif strcmp(symptom_names{1},'Anhedonia')
            R = [Anhedonia,med_regressors,site_regressors,demographic_regressors];
        elseif strcmp(symptom_names{1},'StressxAnhedonia')
            R = [Anhedonia.*lsi_regressors,Anhedonia,lsi_regressors,med_regressors,site_regressors,demographic_regressors];
        elseif strcmp(symptom_names{1},'FearsxAnhedonia')
            R = [Fears.*Anhedonia, Anhedonia, GenDis, Fears, med_regressors,site_regressors,demographic_regressors];
        elseif strcmp(symptom_names{1},'GenDisxFearsxAnhedonia')
            R = [GenDis.*Fears.*Anhedonia,Fears.*Anhedonia,Anhedonia,GenDis.*Fears,GenDis.*Anhedonia,Fears,GenDis,med_regressors,site_regressors,demographic_regressors];
        elseif strcmp(symptom_names{1},'Fears')
            R = [Fears,med_regressors,site_regressors,demographic_regressors];
        end
        intercept = ones(length(R),1);

        %R = [R,intercept];

        R_gain_temp = [R,curr_dat_gain.(region_list{seed_region}).dat'];
        R_loss_temp = [R,curr_dat_loss.(region_list{seed_region}).dat'];

        R_gain_temp(any(isnan(R_gain_temp),2),:)=[];
        R_loss_temp(any(isnan(R_loss_temp),2),:)=[];
        PID(any(isnan(R_loss_temp),2),:)=[];
        R_gain_temp(:,3) = zscore(R_gain_temp(:,3));
        R_loss_temp(:,3) = zscore(R_loss_temp(:,3));

        R_final_gain = R_gain_temp(:,1:size(R,2));
        R_final_loss = R_loss_temp(:,1:size(R,2));

        dat_final_gain = R_gain_temp(:,size(R,2)+1:size(R_gain_temp,2));
        dat_final_loss = R_loss_temp(:,size(R,2)+1:size(R_loss_temp,2));

        curr_dat_gain.(region_list{seed_region}).X = R_final_gain;
        curr_dat_loss.(region_list{seed_region}).X = R_final_loss;

        curr_dat_gain.(region_list{seed_region}).dat = dat_final_gain';
        curr_dat_loss.(region_list{seed_region}).dat = dat_final_loss';
    elseif dsm == 1
        R = [comlife2',deplife2',anxlife2',med_regressors,site_regressors]; %demographic_regressors = sex
        intercept = ones(length(R),1);
        R(R(:,1)==1,2) = 0;
        R(R(:,1)==1,3) = 0;
        %R = [R,intercept];

        R_gain_temp = [R,curr_dat_gain.(region_list{seed_region}).dat'];
        R_loss_temp = [R,curr_dat_loss.(region_list{seed_region}).dat'];

        R_gain_temp(any(isnan(R_gain_temp),2),:)=[];
        R_loss_temp(any(isnan(R_loss_temp),2),:)=[];
        PID(any(isnan(R_loss_temp),2),:)=[];

        R_final_gain = R_gain_temp(:,1:size(R,2));
        R_final_loss = R_loss_temp(:,1:size(R,2));

        dat_final_gain = R_gain_temp(:,size(R,2)+1:size(R_gain_temp,2));
        dat_final_loss = R_loss_temp(:,size(R,2)+1:size(R_loss_temp,2));

        curr_dat_gain.(region_list{seed_region}).X = R_final_gain;
        curr_dat_loss.(region_list{seed_region}).X = R_final_loss;

        curr_dat_gain.(region_list{seed_region}).dat = dat_final_gain';
        curr_dat_loss.(region_list{seed_region}).dat = dat_final_loss';
    elseif life_stress == 1
        R = [lsi_regressors(:,1),med_regressors,site_regressors,demographic_regressors]; %demographic_regressors = sex
        intercept = ones(length(R),1);
        
        R_gain_temp = [R,curr_dat_gain.(region_list{seed_region}).dat'];
        R_loss_temp = [R,curr_dat_loss.(region_list{seed_region}).dat'];

        R_gain_temp(any(isnan(R_gain_temp),2),:)=[];
        R_loss_temp(any(isnan(R_loss_temp),2),:)=[];
        PID(any(isnan(R_loss_temp),2),:)=[];

        R_final_gain = R_gain_temp(:,1:size(R,2));
        R_final_loss = R_loss_temp(:,1:size(R,2));

        dat_final_gain = R_gain_temp(:,size(R,2)+1:size(R_gain_temp,2));
        dat_final_loss = R_loss_temp(:,size(R,2)+1:size(R_loss_temp,2));

        curr_dat_gain.(region_list{seed_region}).X = R_final_gain;
        curr_dat_loss.(region_list{seed_region}).X = R_final_loss;

        curr_dat_gain.(region_list{seed_region}).dat = dat_final_gain';
        curr_dat_loss.(region_list{seed_region}).dat = dat_final_loss';
    end
    %% 5. 
    temp_results_gain.(region_list{seed_region}) = regress(curr_dat_gain.(region_list{seed_region}));%,'robust');
    temp_results_loss.(region_list{seed_region}) = regress(curr_dat_loss.(region_list{seed_region}));%,'robust');
    results_struct.gain.(region_list{seed_region}) = threshold(temp_results_gain.(region_list{seed_region}).t,.05,'fdr','k',1);
    results_struct.loss.(region_list{seed_region}) = threshold(temp_results_loss.(region_list{seed_region}).t,.05,'fdr','k',1);
    %% 6. gain

    % orthviews(results_struct.gain.(region_list{seed_region}))
    %% 6. loss
    % orthviews(results_struct.loss.(region_list{seed_region}))

    %% 7. 
    % CHANGE THIS
    if trilevel == 1
        for symptom = 1:length(symptom_names)
            r_gain.(region_list{seed_region}).(symptom_names{symptom}) = region(select_one_image(results_struct.gain.(region_list{seed_region}),symptom));
            r_loss.(region_list{seed_region}).(symptom_names{symptom}) = region(select_one_image(results_struct.loss.(region_list{seed_region}),symptom));
        end
    end
    % GenDis Anhedonia Fears r_gain r_loss
     

    %% 8. 
    if trilevel == 1
        roidir = '/Users/zaz3744/Documents/current_projects/ACNlab/masks/ROI_BrainMAPD_functional';

        bilateral_roi_fnames = filenames(fullfile(roidir,condir,'*VS*Oldham*.nii'));
        right_roi_fnames = filenames(fullfile(roidir,condir,'right','*VS*Oldham*.nii'));
        left_roi_fnames = filenames(fullfile(roidir,condir,'left','*VS*Oldham*.nii'));
        if strcmp(condir, 'consumption')
            region_name_list = {'Amyg*Ng*','BA9BA46*','bi_vs_sphere*','Caudate*Ng*','HO_Accumbens*','HO_Amyg*','HO_Caudate*','HO_Pallidum*','HO_Putamen*','HO_VMPFC*','OFC*Ng*','VS*Oldham*'}; %'VS*Oldham_Loss*','VS*Oldham_Rew*','VS_Sphere*'};
            %region_name_list = {'HO_Accumbens*','HO_Amyg*','HO_Caudate*','HO_Pallidum*','HO_Putamen*','HO_VMPFC*'}; %'VS*Oldham_Loss*','VS*Oldham_Rew*','VS_Sphere*'};
            
            region_name_list_for_struct = {'Amyg_Ng','BA9BA46','bi_vs_sphere','Caudate_Ng','HO_Accumbens','HO_Amyg','HO_Caudate','HO_Pallidum','HO_Putamen','HO_VMPFC','OFC_Ng','bi_VS_Oldham'};%'VS_Oldham_Loss','VS_Oldham_Rew','VS_sphere'};
            %region_name_list_for_struct = {'HO_Amyg','HO_Caudate','HO_Pallidum','HO_Putamen','HO_VMPFC'};%'VS_Oldham_Loss','VS_Oldham_Rew','VS_sphere'};
        
        else
            region_name_list = {'Amyg*Ng*','BA9BA46*','bi_vs_sphere*','Caudate*Ng*','HO_Accumbens*','HO_Amyg*','HO_Caudate*','HO_Pallidum*','HO_Putamen*','HO_VMPFC*','OFC*Ng*','VS*Oldham_Loss*','VS*Oldham_Rew*','VS_Sphere*'};
            %region_name_list = {'HO_Accumbens*','HO_Amyg*','HO_Caudate*','HO_Pallidum*','HO_Putamen*','HO_VMPFC*'};
            
            region_name_list_for_struct = {'Amyg_Ng','BA9BA46','bi_vs_sphere','Caudate_Ng','HO_Accumbens','HO_Amyg','HO_Caudate','HO_Pallidum','HO_Putamen','HO_VMPFC','OFC_Ng','VS_Oldham_Loss','VS_Oldham_Rew','VS_sphere'};
            %region_name_list_for_struct = {'HO_Amyg','HO_Caudate','HO_Pallidum','HO_Putamen','HO_VMPFC'};%'VS_Oldham_Loss','VS_Oldham_Rew','VS_sphere'};
        end
        % region_name_list = {'VS_loss','VS_reward'};
        for i = 1:length(region_name_list)
            disp(region_name_list_for_struct{i})
            roi_fname = filenames(fullfile(roidir,condir,region_name_list{i}));
            roi = fmri_data(roi_fname{1});
            roi_avg_gain.(region_list{seed_region}).(region_name_list_for_struct{i}) = extract_roi_averages(curr_dat_gain.(region_list{seed_region}),roi);
            roi_avg_gain.(region_list{seed_region}).(region_name_list_for_struct{i}).title = region_name_list_for_struct{i}; 
            roi_gain_table(:,i) = roi_avg_gain.(region_list{seed_region}).(region_name_list_for_struct{i}).dat;
            % Note that the structure that is created below will have two
            % layers. The first refers to the seed region for this analysis.
            % The second will refer to the end region in this seed to seed func
            % conn analysis
            %gain_mdl.(region_list{seed_region}).(region_name_list_for_struct{i}) = fitlm(curr_dat_gain.(region_list{seed_region}).X(:,1:5),roi_avg_gain(i,:).dat);
            %fitlm(curr_dat_gain.(region_list{seed_region}).X(:,1:5),roi_avg_gain(i,:).dat)
            gain_mdl.(region_list{seed_region}).(region_name_list_for_struct{i}) = fitlm(R_final_gain,roi_avg_gain.(region_list{seed_region}).(region_name_list_for_struct{i}).dat);
            fitlm(R_final_gain,roi_avg_gain.(region_list{seed_region}).(region_name_list_for_struct{i}).dat)
            roi_avg_loss.(region_list{seed_region}).(region_name_list_for_struct{i}) = extract_roi_averages(curr_dat_loss.(region_list{seed_region}),roi);
            roi_avg_loss.(region_list{seed_region}).(region_name_list_for_struct{i}).title = region_name_list_for_struct{i};
            roi_loss_table(:,i) = roi_avg_loss.(region_list{seed_region}).(region_name_list_for_struct{i}).dat;
            %loss_mdl.(region_list{seed_region}).(region_name_list_for_struct{i}) = fitlm(curr_dat_loss.(region_list{seed_region}).X(:,1:5),roi_avg_loss(i,:).dat);
            %fitlm(curr_dat_gain.(region_list{seed_region}).X(:,1:5),roi_avg_loss(i,:).dat)
            loss_mdl.(region_list{seed_region}).(region_name_list_for_struct{i}) = fitlm(R_final_loss,roi_avg_loss.(region_list{seed_region}).(region_name_list_for_struct{i}).dat);
            fitlm(R_final_loss,roi_avg_loss.(region_list{seed_region}).(region_name_list_for_struct{i}).dat)
            disp(strcat('Seed region is: ',region_list{seed_region}))
            disp(strcat('End region is: ', region_name_list_for_struct{i}))
            %keyboard
        end
        param_struct.(region_name_list_for_struct{seed_region}).gain = roi_avg_gain.(region_list{seed_region}).(region_name_list_for_struct{i});
        param_struct.(region_name_list_for_struct{seed_region}).loss = roi_avg_loss.(region_list{seed_region}).(region_name_list_for_struct{i});
    %     roi_gain_table = [curr_dat_gain.X,roi_gain_table]; roi_gain_table = array2table(roi_gain_table); 
    %     roi_gain_table.Properties.VariableNames = {'GenDis','Anhedonia','Fears','Ng_Amyg','BA9BA46','bi_vs_sphere','Ng_Caudate','HO_Accumbens','HO_Amyg','HO_Caudate','HO_Pallidum','HO_Putamen','HO_vmPFC','Knutson_mPFC','Knutson_OFC','Ng_OFC','Oldham_loss_VS','Oldham_gain_VS','VS_sphere'};
    %     roi_loss_table = [curr_dat_gain.X,roi_loss_table]; roi_loss_table = array2table(roi_loss_table); 
    %     roi_loss_table.Properties.VariableNames = {'GenDis','Anhedonia','Fears','Ng_Amyg','BA9BA46','bi_vs_sphere','Ng_Caudate','HO_Accumbens','HO_Amyg','HO_Caudate','HO_Pallidum','HO_Putamen','HO_vmPFC','Knutson_mPFC','Knutson_OFC','Ng_OFC','Oldham_loss_VS','Oldham_gain_VS','VS_sphere'};


        %% Same analyses as above but controlling for loss/gain condition

        for i = 1:length(region_name_list)
            disp(region_name_list_for_struct{i})
            roi_fname = filenames(fullfile(roidir,condir,region_name_list{i}));
            roi = fmri_data(roi_fname{1});

            % Note that the structure that is created below will have two
            % layers. The first refers to the seed region for this analysis.
            % The second will refer to the end region in this seed to seed func
            % conn analysis

            % extract roi data for gain and loss
            roi_avg_gain2.(region_list{seed_region}).(region_name_list_for_struct{i}) = extract_roi_averages(curr_dat_gain.(region_list{seed_region}),roi);
            roi_avg_gain2.(region_list{seed_region}).(region_name_list_for_struct{i}).title = region_name_list_for_struct{i}; 
            roi_avg_loss2.(region_list{seed_region}).(region_name_list_for_struct{i}) = extract_roi_averages(curr_dat_loss.(region_list{seed_region}),roi);
            roi_avg_loss2.(region_list{seed_region}).(region_name_list_for_struct{i}).title = region_name_list_for_struct{i};

            % add the loss condition to the gain model to control for
            % activation during loss anticipation
            R_final_gain2 = [R_final_gain,roi_avg_loss2.(region_list{seed_region}).(region_name_list_for_struct{i}).dat];
            gain_mdl2.(region_list{seed_region}).(region_name_list_for_struct{i}) = fitlm(R_final_gain2,roi_avg_gain2.(region_list{seed_region}).(region_name_list_for_struct{i}).dat);
            fitlm(R_final_gain2,roi_avg_gain2.(region_list{seed_region}).(region_name_list_for_struct{i}).dat)

            % add the gain condition to the loss model to control for
            % activation during gain anticipation
            R_final_loss2 = [R_final_loss,roi_avg_gain.(region_list{seed_region}).(region_name_list_for_struct{i}).dat];
            loss_mdl2.(region_list{seed_region}).(region_name_list_for_struct{i}) = fitlm(R_final_loss2,roi_avg_loss2.(region_list{seed_region}).(region_name_list_for_struct{i}).dat);
            fitlm(R_final_loss2,roi_avg_loss2.(region_list{seed_region}).(region_name_list_for_struct{i}).dat)
            disp(strcat('Seed region is: ',region_list{seed_region}))
            disp(strcat('End region is: ', region_name_list_for_struct{i}))
            %keyboard
            % plot residual scatter and model for presentations
    %        symptoms = [GenDis,Anhedonia,Fears];
    %         for symp = 1:length(symptom_names)
    %             % for gain
    %             resid_mdl_R = [R_final_gain2(:,1),R_final_gain2(:,3:6)];
    %             resid_mdl_gain = fitlm(resid_mdl_R,roi_avg_gain2.(region_list{seed_region}).(region_name_list_for_struct{i}).dat);
    %             plot_mdl = fitlm(symptoms(:,symp),resid_mdl_gain.Residuals.Raw)
    %             figure();scatter(symptoms(:,symp),resid_mdl_gain.Residuals.Raw)
    %             h1 = lsline();
    %             h1.LineWidth = 2;
    %             h1.Color = 'r';
    %             r1 = corrcoef(symptoms(:,symp),resid_mdl_gain.Residuals.Raw,'rows','complete');
    %             disp(r1(1,2));
    %             str = [' r = ',num2str(r1(1,2))];
    %             T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str);
    %             set(T, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');
    %             hold on; plot(plot_mdl); title(strcat("Gain ",condir,":",region_list{seed_region},"-",region_name_list_for_struct{i}));xlabel(symptom_names{symp});ylabel("Functional Connectivity");legend off
    %             drawnow, snapnow
    %             filename = strcat(symptom_names{symp},"_Gain",constring,":",region_list{seed_region},"-",region_name_list_for_struct{i},'.jpg');
    %             saveas(gcf,fullfile(savedir,filename))
    %             % for loss
    %             resid_mdl_R = [R_final_loss2(:,1),R_final_loss2(:,3:6)];
    %             resid_mdl_loss = fitlm(resid_mdl_R,roi_avg_loss2.(region_list{seed_region}).(region_name_list_for_struct{i}).dat);
    %             plot_mdl = fitlm(symptoms(:,symp),resid_mdl_loss.Residuals.Raw)
    %             figure();scatter(symptoms(:,symp),resid_mdl_loss.Residuals.Raw)
    %             h1 = lsline();
    %             h1.LineWidth = 2;
    %             h1.Color = 'r';
    %             r1 = corrcoef(symptoms(:,symp),resid_mdl_loss.Residuals.Raw,'rows','complete');
    %             disp(r1(1,2));
    %             str = [' r = ',num2str(r1(1,2))];
    %             T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str);
    %             set(T, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');
    %             hold on; plot(plot_mdl); title(strcat("Loss ",condir,":",region_list{seed_region},"-",region_name_list_for_struct{i}));xlabel(symptom_names{symp});ylabel("Functional Connectivity");legend off
    %             drawnow, snapnow
    %             filename = strcat(symptom_names{symp},"_Loss",constring,":",region_list{seed_region},"-",region_name_list_for_struct{i},'.jpg');
    %             saveas(gcf,fullfile(savedir,filename))
    %        end
            close all
        end
    end
    %% Same analyses but for dsm diagnoses and controlling for loss/gain condition
    if dsm == 1
        roidir = '/Users/zaz3744/Documents/current_projects/ACNlab/masks/ROI_BrainMAPD_functional';

        bilateral_roi_fnames = filenames(fullfile(roidir,condir,'*VS*Oldham*.nii'));
        right_roi_fnames = filenames(fullfile(roidir,condir,'right','*VS*Oldham*.nii'));
        left_roi_fnames = filenames(fullfile(roidir,condir,'left','*VS*Oldham*.nii'));
        if strcmp(condir, 'consumption')
            region_name_list = {'Amyg*Ng*','BA9BA46*','bi_vs_sphere*','Caudate*Ng*','HO_Accumbens*','HO_Amyg*','HO_Caudate*','HO_Pallidum*','HO_Putamen*','HO_VMPFC*','OFC*Ng*','VS*Oldham*'}; %'VS*Oldham_Loss*','VS*Oldham_Rew*','VS_Sphere*'};
            region_name_list_for_struct = {'Amyg_Ng','BA9BA46','bi_vs_sphere','Caudate_Ng','HO_Accumbens','HO_Amyg','HO_Caudate','HO_Pallidum','HO_Putamen','HO_VMPFC','OFC_Ng','bi_VS_Oldham'};%'VS_Oldham_Loss','VS_Oldham_Rew','VS_sphere'};
        else
            region_name_list = {'Amyg*Ng*','BA9BA46*','bi_vs_sphere*','Caudate*Ng*','HO_Accumbens*','HO_Amyg*','HO_Caudate*','HO_Pallidum*','HO_Putamen*','HO_VMPFC*','OFC*Ng*','VS*Oldham_Loss*','VS*Oldham_Rew*','VS_Sphere*'};
            region_name_list_for_struct = {'Amyg_Ng','BA9BA46','bi_vs_sphere','Caudate_Ng','HO_Accumbens','HO_Amyg','HO_Caudate','HO_Pallidum','HO_Putamen','HO_VMPFC','OFC_Ng','VS_Oldham_Loss','VS_Oldham_Rew','VS_sphere'};
        end
        % region_name_list = {'VS_loss','VS_reward'};
        for i = 1:length(region_name_list)
            disp(region_name_list_for_struct{i})
            roi_fname = filenames(fullfile(roidir,condir,region_name_list{i}));
            roi = fmri_data(roi_fname{1});
            roi_avg_gain.(region_list{seed_region}).(region_name_list_for_struct{i}) = extract_roi_averages(curr_dat_gain.(region_list{seed_region}),roi);
            roi_avg_gain.(region_list{seed_region}).(region_name_list_for_struct{i}).title = region_name_list_for_struct{i}; 
            roi_gain_table(:,i) = roi_avg_gain.(region_list{seed_region}).(region_name_list_for_struct{i}).dat;
            % Note that the structure that is created below will have two
            % layers. The first refers to the seed region for this analysis.
            % The second will refer to the end region in this seed to seed func
            % conn analysis
            %gain_mdl.(region_list{seed_region}).(region_name_list_for_struct{i}) = fitlm(curr_dat_gain.(region_list{seed_region}).X(:,1:5),roi_avg_gain(i,:).dat);
            %fitlm(curr_dat_gain.(region_list{seed_region}).X(:,1:5),roi_avg_gain(i,:).dat)
            gain_mdl.(region_list{seed_region}).(region_name_list_for_struct{i}) = fitlm(R_final_gain,roi_avg_gain.(region_list{seed_region}).(region_name_list_for_struct{i}).dat);
            fitlm(R_final_gain,roi_avg_gain.(region_list{seed_region}).(region_name_list_for_struct{i}).dat)
            roi_avg_loss.(region_list{seed_region}).(region_name_list_for_struct{i}) = extract_roi_averages(curr_dat_loss.(region_list{seed_region}),roi);
            roi_avg_loss.(region_list{seed_region}).(region_name_list_for_struct{i}).title = region_name_list_for_struct{i};
            roi_loss_table(:,i) = roi_avg_loss.(region_list{seed_region}).(region_name_list_for_struct{i}).dat;
            %loss_mdl.(region_list{seed_region}).(region_name_list_for_struct{i}) = fitlm(curr_dat_loss.(region_list{seed_region}).X(:,1:5),roi_avg_loss(i,:).dat);
            %fitlm(curr_dat_gain.(region_list{seed_region}).X(:,1:5),roi_avg_loss(i,:).dat)
            loss_mdl.(region_list{seed_region}).(region_name_list_for_struct{i}) = fitlm(R_final_loss,roi_avg_loss.(region_list{seed_region}).(region_name_list_for_struct{i}).dat);
            fitlm(R_final_loss,roi_avg_loss.(region_list{seed_region}).(region_name_list_for_struct{i}).dat)
            disp(strcat('Seed region is: ',region_list{seed_region}))
            disp(strcat('End region is: ', region_name_list_for_struct{i}))
            %keyboard
        end
        
        % same analysis controlling for opposite condition
        for i = 1:length(region_name_list)
            disp(region_name_list_for_struct{i})
            roi_fname = filenames(fullfile(roidir,condir,region_name_list{i}));
            roi = fmri_data(roi_fname{1});

            % Note that the structure that is created below will have two
            % layers. The first refers to the seed region for this analysis.
            % The second will refer to the end region in this seed to seed func
            % conn analysis

            % extract roi data for gain and loss
            roi_avg_gain2.(region_list{seed_region}).(region_name_list_for_struct{i}) = extract_roi_averages(curr_dat_gain.(region_list{seed_region}),roi);
            roi_avg_gain2.(region_list{seed_region}).(region_name_list_for_struct{i}).title = region_name_list_for_struct{i}; 
            roi_avg_loss2.(region_list{seed_region}).(region_name_list_for_struct{i}) = extract_roi_averages(curr_dat_loss.(region_list{seed_region}),roi);
            roi_avg_loss2.(region_list{seed_region}).(region_name_list_for_struct{i}).title = region_name_list_for_struct{i};

            % add the loss condition to the gain model to control for
            % activation during loss anticipation
            R_final_gain2 = [R_final_gain,roi_avg_loss2.(region_list{seed_region}).(region_name_list_for_struct{i}).dat];
            gain_mdl2.(region_list{seed_region}).(region_name_list_for_struct{i}) = fitlm(R_final_gain2,roi_avg_gain2.(region_list{seed_region}).(region_name_list_for_struct{i}).dat);
            fitlm(R_final_gain2,roi_avg_gain2.(region_list{seed_region}).(region_name_list_for_struct{i}).dat)

            % add the gain condition to the loss model to control for
            % activation during gain anticipation
            R_final_loss2 = [R_final_loss,roi_avg_gain.(region_list{seed_region}).(region_name_list_for_struct{i}).dat];
            loss_mdl2.(region_list{seed_region}).(region_name_list_for_struct{i}) = fitlm(R_final_loss2,roi_avg_loss2.(region_list{seed_region}).(region_name_list_for_struct{i}).dat);
            fitlm(R_final_loss2,roi_avg_loss2.(region_list{seed_region}).(region_name_list_for_struct{i}).dat)
            disp(strcat('Seed region is: ',region_list{seed_region}))
            disp(strcat('End region is: ', region_name_list_for_struct{i}))
            %keyboard
            % plot residual scatter and model for presentations
    %        symptoms = [GenDis,Anhedonia,Fears];
    %         for symp = 1:length(symptom_names)
    %             % for gain
    %             resid_mdl_R = [R_final_gain2(:,1),R_final_gain2(:,3:6)];
    %             resid_mdl_gain = fitlm(resid_mdl_R,roi_avg_gain2.(region_list{seed_region}).(region_name_list_for_struct{i}).dat);
    %             plot_mdl = fitlm(symptoms(:,symp),resid_mdl_gain.Residuals.Raw)
    %             figure();scatter(symptoms(:,symp),resid_mdl_gain.Residuals.Raw)
    %             h1 = lsline();
    %             h1.LineWidth = 2;
    %             h1.Color = 'r';
    %             r1 = corrcoef(symptoms(:,symp),resid_mdl_gain.Residuals.Raw,'rows','complete');
    %             disp(r1(1,2));
    %             str = [' r = ',num2str(r1(1,2))];
    %             T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str);
    %             set(T, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');
    %             hold on; plot(plot_mdl); title(strcat("Gain ",condir,":",region_list{seed_region},"-",region_name_list_for_struct{i}));xlabel(symptom_names{symp});ylabel("Functional Connectivity");legend off
    %             drawnow, snapnow
    %             filename = strcat(symptom_names{symp},"_Gain",constring,":",region_list{seed_region},"-",region_name_list_for_struct{i},'.jpg');
    %             saveas(gcf,fullfile(savedir,filename))
    %             % for loss
    %             resid_mdl_R = [R_final_loss2(:,1),R_final_loss2(:,3:6)];
    %             resid_mdl_loss = fitlm(resid_mdl_R,roi_avg_loss2.(region_list{seed_region}).(region_name_list_for_struct{i}).dat);
    %             plot_mdl = fitlm(symptoms(:,symp),resid_mdl_loss.Residuals.Raw)
    %             figure();scatter(symptoms(:,symp),resid_mdl_loss.Residuals.Raw)
    %             h1 = lsline();
    %             h1.LineWidth = 2;
    %             h1.Color = 'r';
    %             r1 = corrcoef(symptoms(:,symp),resid_mdl_loss.Residuals.Raw,'rows','complete');
    %             disp(r1(1,2));
    %             str = [' r = ',num2str(r1(1,2))];
    %             T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str);
    %             set(T, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');
    %             hold on; plot(plot_mdl); title(strcat("Loss ",condir,":",region_list{seed_region},"-",region_name_list_for_struct{i}));xlabel(symptom_names{symp});ylabel("Functional Connectivity");legend off
    %             drawnow, snapnow
    %             filename = strcat(symptom_names{symp},"_Loss",constring,":",region_list{seed_region},"-",region_name_list_for_struct{i},'.jpg');
    %             saveas(gcf,fullfile(savedir,filename))
    %        end
            close all
        end
    end
    
    %% I guess I'm doing one more round for life stress
    if life_stress == 1
        roidir = '/Users/zaz3744/Documents/current_projects/ACNlab/masks/ROI_BrainMAPD_functional';

        bilateral_roi_fnames = filenames(fullfile(roidir,condir,'*VS*Oldham*.nii'));
        right_roi_fnames = filenames(fullfile(roidir,condir,'right','*VS*Oldham*.nii'));
        left_roi_fnames = filenames(fullfile(roidir,condir,'left','*VS*Oldham*.nii'));
        if strcmp(condir, 'consumption')
            region_name_list = {'Amyg*Ng*','BA9BA46*','bi_vs_sphere*','Caudate*Ng*','HO_Accumbens*','HO_Amyg*','HO_Caudate*','HO_Pallidum*','HO_Putamen*','HO_VMPFC*','OFC*Ng*','VS*Oldham*'}; %'VS*Oldham_Loss*','VS*Oldham_Rew*','VS_Sphere*'};
            region_name_list_for_struct = {'Amyg_Ng','BA9BA46','bi_vs_sphere','Caudate_Ng','HO_Accumbens','HO_Amyg','HO_Caudate','HO_Pallidum','HO_Putamen','HO_VMPFC','OFC_Ng','bi_VS_Oldham'};%'VS_Oldham_Loss','VS_Oldham_Rew','VS_sphere'};
        else
            region_name_list = {'Amyg*Ng*','BA9BA46*','bi_vs_sphere*','Caudate*Ng*','HO_Accumbens*','HO_Amyg*','HO_Caudate*','HO_Pallidum*','HO_Putamen*','HO_VMPFC*','OFC*Ng*','VS*Oldham_Loss*','VS*Oldham_Rew*','VS_Sphere*'};
            region_name_list_for_struct = {'Amyg_Ng','BA9BA46','bi_vs_sphere','Caudate_Ng','HO_Accumbens','HO_Amyg','HO_Caudate','HO_Pallidum','HO_Putamen','HO_VMPFC','OFC_Ng','VS_Oldham_Loss','VS_Oldham_Rew','VS_sphere'};
        end
        % region_name_list = {'VS_loss','VS_reward'};
        for i = 1:length(region_name_list)
            disp(region_name_list_for_struct{i})
            roi_fname = filenames(fullfile(roidir,condir,region_name_list{i}));
            roi = fmri_data(roi_fname{1});
            roi_avg_gain.(region_list{seed_region}).(region_name_list_for_struct{i}) = extract_roi_averages(curr_dat_gain.(region_list{seed_region}),roi);
            roi_avg_gain.(region_list{seed_region}).(region_name_list_for_struct{i}).title = region_name_list_for_struct{i}; 
            roi_gain_table(:,i) = roi_avg_gain.(region_list{seed_region}).(region_name_list_for_struct{i}).dat;
            % Note that the structure that is created below will have two
            % layers. The first refers to the seed region for this analysis.
            % The second will refer to the end region in this seed to seed func
            % conn analysis
            %gain_mdl.(region_list{seed_region}).(region_name_list_for_struct{i}) = fitlm(curr_dat_gain.(region_list{seed_region}).X(:,1:5),roi_avg_gain(i,:).dat);
            %fitlm(curr_dat_gain.(region_list{seed_region}).X(:,1:5),roi_avg_gain(i,:).dat)
            gain_mdl.(region_list{seed_region}).(region_name_list_for_struct{i}) = fitlm(R_final_gain,roi_avg_gain.(region_list{seed_region}).(region_name_list_for_struct{i}).dat);
            fitlm(R_final_gain,roi_avg_gain.(region_list{seed_region}).(region_name_list_for_struct{i}).dat)
            roi_avg_loss.(region_list{seed_region}).(region_name_list_for_struct{i}) = extract_roi_averages(curr_dat_loss.(region_list{seed_region}),roi);
            roi_avg_loss.(region_list{seed_region}).(region_name_list_for_struct{i}).title = region_name_list_for_struct{i};
            roi_loss_table(:,i) = roi_avg_loss.(region_list{seed_region}).(region_name_list_for_struct{i}).dat;
            %loss_mdl.(region_list{seed_region}).(region_name_list_for_struct{i}) = fitlm(curr_dat_loss.(region_list{seed_region}).X(:,1:5),roi_avg_loss(i,:).dat);
            %fitlm(curr_dat_gain.(region_list{seed_region}).X(:,1:5),roi_avg_loss(i,:).dat)
            loss_mdl.(region_list{seed_region}).(region_name_list_for_struct{i}) = fitlm(R_final_loss,roi_avg_loss.(region_list{seed_region}).(region_name_list_for_struct{i}).dat);
            fitlm(R_final_loss,roi_avg_loss.(region_list{seed_region}).(region_name_list_for_struct{i}).dat)
            disp(strcat('Seed region is: ',region_list{seed_region}))
            disp(strcat('End region is: ', region_name_list_for_struct{i}))
            %keyboard
        end
        
        % same analysis controlling for opposite condition
        for i = 1:length(region_name_list)
            disp(region_name_list_for_struct{i})
            roi_fname = filenames(fullfile(roidir,condir,region_name_list{i}));
            roi = fmri_data(roi_fname{1});

            % Note that the structure that is created below will have two
            % layers. The first refers to the seed region for this analysis.
            % The second will refer to the end region in this seed to seed func
            % conn analysis

            % extract roi data for gain and loss
            roi_avg_gain2.(region_list{seed_region}).(region_name_list_for_struct{i}) = extract_roi_averages(curr_dat_gain.(region_list{seed_region}),roi);
            roi_avg_gain2.(region_list{seed_region}).(region_name_list_for_struct{i}).title = region_name_list_for_struct{i}; 
            roi_avg_loss2.(region_list{seed_region}).(region_name_list_for_struct{i}) = extract_roi_averages(curr_dat_loss.(region_list{seed_region}),roi);
            roi_avg_loss2.(region_list{seed_region}).(region_name_list_for_struct{i}).title = region_name_list_for_struct{i};

            % add the loss condition to the gain model to control for
            % activation during loss anticipation
            R_final_gain2 = [R_final_gain,roi_avg_loss2.(region_list{seed_region}).(region_name_list_for_struct{i}).dat];
            gain_mdl2.(region_list{seed_region}).(region_name_list_for_struct{i}) = fitlm(R_final_gain2,roi_avg_gain2.(region_list{seed_region}).(region_name_list_for_struct{i}).dat);
            fitlm(R_final_gain2,roi_avg_gain2.(region_list{seed_region}).(region_name_list_for_struct{i}).dat)

            % add the gain condition to the loss model to control for
            % activation during gain anticipation
            R_final_loss2 = [R_final_loss,roi_avg_gain.(region_list{seed_region}).(region_name_list_for_struct{i}).dat];
            loss_mdl2.(region_list{seed_region}).(region_name_list_for_struct{i}) = fitlm(R_final_loss2,roi_avg_loss2.(region_list{seed_region}).(region_name_list_for_struct{i}).dat);
            fitlm(R_final_loss2,roi_avg_loss2.(region_list{seed_region}).(region_name_list_for_struct{i}).dat)
            disp(strcat('Seed region is: ',region_list{seed_region}))
            disp(strcat('End region is: ', region_name_list_for_struct{i}))
            %keyboard
            % plot residual scatter and model for presentations
    %        symptoms = [GenDis,Anhedonia,Fears];
    %         for symp = 1:length(symptom_names)
    %             % for gain
    %             resid_mdl_R = [R_final_gain2(:,1),R_final_gain2(:,3:6)];
    %             resid_mdl_gain = fitlm(resid_mdl_R,roi_avg_gain2.(region_list{seed_region}).(region_name_list_for_struct{i}).dat);
    %             plot_mdl = fitlm(symptoms(:,symp),resid_mdl_gain.Residuals.Raw)
    %             figure();scatter(symptoms(:,symp),resid_mdl_gain.Residuals.Raw)
    %             h1 = lsline();
    %             h1.LineWidth = 2;
    %             h1.Color = 'r';
    %             r1 = corrcoef(symptoms(:,symp),resid_mdl_gain.Residuals.Raw,'rows','complete');
    %             disp(r1(1,2));
    %             str = [' r = ',num2str(r1(1,2))];
    %             T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str);
    %             set(T, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');
    %             hold on; plot(plot_mdl); title(strcat("Gain ",condir,":",region_list{seed_region},"-",region_name_list_for_struct{i}));xlabel(symptom_names{symp});ylabel("Functional Connectivity");legend off
    %             drawnow, snapnow
    %             filename = strcat(symptom_names{symp},"_Gain",constring,":",region_list{seed_region},"-",region_name_list_for_struct{i},'.jpg');
    %             saveas(gcf,fullfile(savedir,filename))
    %             % for loss
    %             resid_mdl_R = [R_final_loss2(:,1),R_final_loss2(:,3:6)];
    %             resid_mdl_loss = fitlm(resid_mdl_R,roi_avg_loss2.(region_list{seed_region}).(region_name_list_for_struct{i}).dat);
    %             plot_mdl = fitlm(symptoms(:,symp),resid_mdl_loss.Residuals.Raw)
    %             figure();scatter(symptoms(:,symp),resid_mdl_loss.Residuals.Raw)
    %             h1 = lsline();
    %             h1.LineWidth = 2;
    %             h1.Color = 'r';
    %             r1 = corrcoef(symptoms(:,symp),resid_mdl_loss.Residuals.Raw,'rows','complete');
    %             disp(r1(1,2));
    %             str = [' r = ',num2str(r1(1,2))];
    %             T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str);
    %             set(T, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');
    %             hold on; plot(plot_mdl); title(strcat("Loss ",condir,":",region_list{seed_region},"-",region_name_list_for_struct{i}));xlabel(symptom_names{symp});ylabel("Functional Connectivity");legend off
    %             drawnow, snapnow
    %             filename = strcat(symptom_names{symp},"_Loss",constring,":",region_list{seed_region},"-",region_name_list_for_struct{i},'.jpg');
    %             saveas(gcf,fullfile(savedir,filename))
    %        end
            close all
        end
    end
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

%% Seed to seed while controlling for the additional condition
% region_list = {'bi_VS_AntRew','bi_VS_AntLoss','bi_vmPFC_to_wholebrain','bi_Amyg_to_wholebrain','biOFC_to_wholebrain','L_VS_AntRew_to_wholebrain','L_VS_AntLoss_to_wholebrain','LHO_Accumbens_to_wholebrain','LOFC2_to_wholebrain','LOFC_to_wholebrain','R_VS_AntLoss_to_wholebrain','R_VS_AntRew_to_wholebrain','RHO_Accumbens_to_wholebrain','ROFC_to_wholebrain'};
% 
% new_region_list = {
% 
% for 
% 
% temp_results_gain.(region_list{seed_region})
% temp_results_loss.(region_list{seed_region})
% 
% task_regressor
% original_regressor
% current_outcome
% 
% data_tbl = [current_outcome, task_regressor, original_regressor];
% 
% loss_mdl_final = fitlm(data_tbl,'current_outcome ~ task_regressor + original_regressor');



% helpful for stuff

% X = [param_struct.biOFC_to_wholebrain.gain(2).dat,param_struct.biOFC_to_wholebrain.gain(1).dat,param_struct.ROFC_to_wholebrain.gain(2).dat,param_struct.ROFC_to_wholebrain.gain(1).dat,param_struct.LOFC_to_wholebrain.gain(2).dat,param_struct.LOFC_to_wholebrain.gain(1).dat,param_struct.LOFC2_to_wholebrain.gain(2).dat,param_struct.LOFC2_to_wholebrain.gain(1).dat];
% X = [param_struct.biOFC_to_wholebrain.loss(2).dat,param_struct.biOFC_to_wholebrain.loss(1).dat,param_struct.ROFC_to_wholebrain.loss(2).dat,param_struct.ROFC_to_wholebrain.loss(1).dat,param_struct.LOFC_to_wholebrain.loss(2).dat,param_struct.LOFC_to_wholebrain.loss(1).dat,param_struct.LOFC2_to_wholebrain.loss(2).dat,param_struct.LOFC2_to_wholebrain.loss(1).dat];

%% Final Results Report

% region_list = {'bi_VS_AntRew','bi_VS_AntLoss','bi_vmPFC_to_wholebrain','bi_Amyg_to_wholebrain','biOFC_to_wholebrain','L_VS_AntRew_to_wholebrain',...
% 'L_VS_AntLoss_to_wholebrain','LHO_Accumbens_to_wholebrain','LOFC2_to_wholebrain','LOFC_to_wholebrain','R_VS_AntLoss_to_wholebrain','R_VS_AntRew_to_wholebrain',...
% 'RHO_Accumbens_to_wholebrain','ROFC_to_wholebrain'};

% specifiy the regions of interest you actually want to do for wholebrain
% analysis

%% Final Results Report

% The code below is out of date and all final results reporting occurs with
% canlab_help_publish('publish_first_year_project_func_conn.m')
% - ZA (10/26/20)



% % region_list = {'bi_VS_AntRew','bi_VS_AntLoss','bi_vmPFC_to_wholebrain','bi_Amyg_to_wholebrain','biOFC_to_wholebrain','L_VS_AntRew_to_wholebrain',...
% % 'L_VS_AntLoss_to_wholebrain','LHO_Accumbens_to_wholebrain','LOFC2_to_wholebrain','LOFC_to_wholebrain','R_VS_AntLoss_to_wholebrain','R_VS_AntRew_to_wholebrain',...
% % 'RHO_Accumbens_to_wholebrain','ROFC_to_wholebrain'};
% 
% % specifiy the regions of interest you actually want to do for wholebrain
% % analysis
% 
% specific_roi_list = {'bi_VS_AntRew','bi_VS_AntLoss','biOFC_to_wholebrain','LOFC2_to_wholebrain','LOFC_to_wholebrain','ROFC_to_wholebrain','L_VS_AntRew_to_wholebrain','L_VS_AntLoss_to_wholebrain','R_VS_AntLoss_to_wholebrain','R_VS_AntRew_to_wholebrain'};
% analyze = 1;
% if analyze == 1
%     clear r_gain r_loss
%     % report whole brain clusters
%     for i = 1:length(specific_roi_list)
%         results_struct.gain.(specific_roi_list{i}) = threshold(temp_results_gain.(specific_roi_list{i}).t,.0005,'unc','k',10);
%         results_struct.loss.(specific_roi_list{i}) = threshold(temp_results_loss.(specific_roi_list{i}).t,.0005,'unc','k',10);
%         
%         for s = 1:length(symptom_names)
%             
%             % Current Seed Region
%             disp('GAIN ANTICIPTATION')
%             disp(specific_roi_list{i})
%             % Current Symptom
%             disp(symptom_names{s})
%             
%             disp('_________________________________________________________________________')
%             
%             r_gain.(specific_roi_list{i}).(symptom_names{s}) = region(select_one_image(results_struct.gain.(specific_roi_list{i}),s));
%             info_gain = table(r_gain.(specific_roi_list{i}).(symptom_names{s}));
%             
%             if isempty(info_gain) == 1
%                 continue
%             else
%                 create_figure(strcat('Gain anticipation montage Seed:',specific_roi_list{i},' Symptom:',symptom_names{s})); axis off;
%                 canlab_results_fmridisplay(r_gain.(specific_roi_list{i}).(symptom_names{s}), 'regioncenters')
%                 drawnow, snapnow
%                 filename = strcat('Gain_anticipation_montage_Seed_',specific_roi_list{i},'_Symptom_',symptom_names{s},'.pdf');
%                 disp('_________________________________________________________________________')
%                 %keyboard
%                 %saveas(gcf,fullfile(figdir,filename))
%             end
%             
%             
%             % Current Seed Region
%             disp('LOSS ANTICIPATION')
%             disp(specific_roi_list{i})
%             % Current Symptom
%             disp(symptom_names{s})
%             
%             disp('_________________________________________________________________________')
%             
%             r_loss.(specific_roi_list{i}).(symptom_names{s}) = region(select_one_image(results_struct.loss.(specific_roi_list{i}),s));
%             info_loss = table(r_loss.(specific_roi_list{i}).(symptom_names{s}));
%             
%             if isempty(info_loss) == 1
%                 continue
%             else
%                 
%                 create_figure(strcat('Loss anticipation montage Seed:',specific_roi_list{i},' Symptom:',symptom_names{s})); axis off;
%                 canlab_results_fmridisplay(r_loss.(specific_roi_list{i}).(symptom_names{s}), 'regioncenters');
%                 drawnow, snapnow
%                 filename = strcat('Loss_anticipation_montage_Seed_',specific_roi_list{i},'_Symptom_',symptom_names{s},'.pdf');
%                 disp('_________________________________________________________________________')
%                 %keyboard
%                 %saveas(gcf,fullfile(figdir,filename))
%             end
%         end
%     end
%     
%     % report region of interest analysis
%     for i = 1:length(specific_roi_list)
%         for r = 1:length(region_name_list)
%             % So basically we're looping through every seed region, and
%             % then kind of every seed region again... But one refers to the
%             % seed and the next refers to the region of interest extracted
%             % for results. Then we'll check the p value associated with the
%             % relationship and will pause to explore things further if
%             % there's a significant relationship. I hope this works.
%             test_gain = gain_mdl.(specific_roi_list{i}).(region_name_list_for_struct{r});
%             if test_gain.Coefficients.pValue(2) < 0.05 || test_gain.Coefficients.pValue(3) < 0.05 || test_gain.Coefficients.pValue(4) < 0.05
%                 % This is going to create a results structure that will be
%                 % helpful for me to try and parse all the results I have.
%                 % Controlling for sex really opened the box
%                 if test_gain.Coefficients.pValue(2) < 0.05
%                     final_seed_to_seed_results.gain.GenDis.(specific_roi_list{i}).(region_name_list_for_struct{r}) = test_gain;
%                     disp('During Gain Anticipation')
%                     disp('Symptom: General Distress')
%                     disp(strcat('Seed region is: ',region_list{i}))
%                     disp(strcat('End region is: ', region_name_list_for_struct{r}))
%                     disp(strcat('t = ',num2str(test_gain.Coefficients.tStat(2)),'p = ',num2str(test_gain.Coefficients.pValue(2))))
%                     disp('________________________________________________________________')
%                 elseif test_gain.Coefficients.pValue(3) < 0.05
%                     final_seed_to_seed_results.gain.Anhedonia.(specific_roi_list{i}).(region_name_list_for_struct{r}) = test_gain;
%                     disp('During Gain Anticipation')
%                     disp('Symptom: Anhedonia')
%                     disp(strcat('Seed region is: ',region_list{i}))
%                     disp(strcat('End region is: ', region_name_list_for_struct{r}))
%                     disp(strcat('t = ',num2str(test_gain.Coefficients.tStat(3)),'p = ',num2str(test_gain.Coefficients.pValue(3))))
%                     disp('________________________________________________________________')
% 
%                 elseif test_gain.Coefficients.pValue(4) < 0.05
%                     final_seed_to_seed_results.gain.Fears.(specific_roi_list{i}).(region_name_list_for_struct{r}) = test_gain;
%                     disp('During Gain Anticipation') 
%                     disp('Symptom: Fears')
%                     disp(strcat('Seed region is: ',region_list{i}))
%                     disp(strcat('End region is: ', region_name_list_for_struct{r}))
%                     disp(strcat('t = ',num2str(test_gain.Coefficients.tStat(4)),'p = ',num2str(test_gain.Coefficients.pValue(4))))
%                     disp('________________________________________________________________')
% 
%                 end
%                     %keyboard
%             else
% %                     disp(strcat('Seed region is: ',specific_roi_list{i}))
% %                     disp(strcat('End region is: ', region_name_list_for_struct{r}))
% %                     disp('NO SIG GAIN RESULTS')
%             end
%             test_loss = loss_mdl.(region_list{i}).(region_name_list_for_struct{r});
%             if test_loss.Coefficients.pValue(2) < 0.05 || test_loss.Coefficients.pValue(3) < 0.05 || test_loss.Coefficients.pValue(4) < 0.05
%                 if test_loss.Coefficients.pValue(2) < 0.05
%                     final_seed_to_seed_results.loss.GenDis.(region_list{i}).(region_name_list_for_struct{r}) = test_loss;
%                     disp('During Loss Anticipation')
%                     disp('Symptom: General Distress')
%                     disp(strcat('Seed region is: ',region_list{i}))
%                     disp(strcat('End region is: ', region_name_list_for_struct{r}))
%                     disp(strcat('t = ',num2str(test_loss.Coefficients.tStat(2)),'p = ',num2str(test_loss.Coefficients.pValue(2))))
%                     disp('________________________________________________________________')
%                     
%                 elseif test_loss.Coefficients.pValue(3) < 0.05
%                     final_seed_to_seed_results.loss.Anhedonia.(region_list{i}).(region_name_list_for_struct{r}) = test_loss;
%                     disp('During Loss Anticipation')
%                     disp('Symptom: Anhedonia')
%                     disp(strcat('Seed region is: ',region_list{i}))
%                     disp(strcat('End region is: ', region_name_list_for_struct{r}))
%                     disp(strcat('t = ',num2str(test_loss.Coefficients.tStat(3)),'p = ',num2str(test_loss.Coefficients.pValue(3))))
%                     disp('________________________________________________________________')
%                     
%                 elseif test_loss.Coefficients.pValue(4) < 0.05
%                     final_seed_to_seed_results.loss.Fears.(region_list{i}).(region_name_list_for_struct{r}) = test_loss;
%                     disp('During Loss Anticipation')
%                     disp('Symptom: Fears')
%                     disp(strcat('Seed region is: ',region_list{i}))
%                     disp(strcat('End region is: ', region_name_list_for_struct{r}))
%                     disp(strcat('t = ',num2str(test_loss.Coefficients.tStat(4)),'p = ',num2str(test_loss.Coefficients.pValue(4))))
%                     disp('________________________________________________________________')
% 
%                 end
%             else
% %                     disp(strcat('Seed region is: ',region_list{i}))
% %                     disp(strcat('End region is: ', region_name_list_for_struct{r}))
% %                     disp('NO SIG LOSS RESULTS')
%             end
%         end
%         
%     end
% end
