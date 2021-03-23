%% I want to use RSA
% And I want to do it to more precisely test if patterns of dissimilarity
% in specitic regrions of interest across people correspond with trilevel
% factors. 

% Round 1 of this will look at average activation in a set of ROI's. I'm
% thinking bilateral masks first. 

% PPI: mOFC - Ng Amyg, Oldham VS loss/rew - Ng Amyg, mOFC - Oldham VS loss/rew

% 1. Run func conn script from first year project
% 2. Concatenate subject data from regions. This should take the form of a
% cell array. This will involve creating picking one of the above PPI
% associations and lining up the four contrasts for a single subject next
% to one another. We will calculate dissimilarity across these conditions
% for each subject with respect to these specific ROIs
% 3. Calculate dissimilarity 1 - corr matrix in every person for each of
% the associations above. 
% 4. Pull out overall dissimilrity estimates and put them into a single
% matrix for all subjects
% 5. We'll see what we have at that point, struggling to visualize, but the
% idea will then be to predict symptom on the basis of this dissimilarity
whole_brain_analysis = 0;
roi_analysis = 1;

contrastdir = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD/func_conn/PPI/ppi_fldir';
clinicaldir = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD/clinical_data';
maskdir = '/Users/zaz3744/Documents/current_projects/ACNlab/masks/ROI_BrainMAPD_functional';
med_dir = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD';
demo_dir = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD';
stressdir = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD/life_stress';


eventdir = 'anticipation';
brain_region_dir = 'biOFC';
seed_to_seed_end = 'OFC*Oldham*'; % other ones: VS*Oldham_Loss*, VS*Oldham_Rew*, HO*Amygdala*, Amyg*Ng*
iterations = 1000;
alpha = 0.001;

%% Calculate whole brain disimilarity for Gain all - Gain 0, Loss all - Loss 0
cd(fullfile(contrastdir, eventdir))
fnames_wholebrain_gain_all = filenames(fullfile('*/*/*/*/con_0001.nii'));
fnames_wholebrain_gain_0 = filenames(fullfile('*/*/*/*/con_0003.nii'));
fnames_wholebrain_loss_all = filenames(fullfile('*/*/*/*/con_0002.nii'));
fnames_wholebrain_loss_0 = filenames(fullfile('*/*/*/*/con_0004.nii'));

region_list = {'Amyg*Ng*','VS*Oldham_Loss*','VS*Oldham_Rew*'};
region_names = {'Amygdala','VS_Loss','VS_Rew'};

for sub = 1:length(fnames_wholebrain_gain_all)
    
    % load each of four file types. Each was generated as the unique effect
    % of a single trial type on the brain. ie [1 0 0 0];
    ppi_data_gain_all = fmri_data(fnames_wholebrain_gain_all{sub});
    ppi_data_gain_0 = fmri_data(fnames_wholebrain_gain_0{sub});
    ppi_data_loss_all = fmri_data(fnames_wholebrain_loss_all{sub});
    ppi_data_loss_0 = fmri_data(fnames_wholebrain_loss_0{sub});        

    % run correlations to get whole brain similarity generating RSA "contrasts"
    ppi_corr_mat_gain = corrcoef(ppi_data_gain_all.dat,ppi_data_gain_0.dat);
    ppi_corr_mat_loss = corrcoef(ppi_data_loss_all.dat,ppi_data_loss_0.dat);

    % generate whole brain disimilarity
    ppi_disim_gain = 1 - ppi_corr_mat_gain;
    ppi_disim_loss = 1 - ppi_corr_mat_loss;
    ppi_disim_dat_gain(:,sub) = ppi_data_gain_all.dat;
    ppi_disim_dat_gain_neutral(:,sub) = ppi_data_gain_0.dat;
    ppi_disim_dat_loss(:,sub) = ppi_data_loss_all.dat;
    ppi_disim_dat_loss_neutral(:,sub) = ppi_data_loss_0.dat;

    % store each whole brain disimilarity value for subs in a single vector
    PID(sub,1) = str2num(fnames_wholebrain_gain_all{sub}(1:5));
    results_vector_gain(sub,1) = ppi_disim_gain(2,1);
    results_vector_loss(sub,1) = ppi_disim_loss(2,1);
    
    if roi_analysis == 1
        for roi = 1:length(region_list)
            % let's get seed to seed PPI RSA going
            % load in a specific seed
            disp(region_names{roi})
            roi_fname = filenames(fullfile(maskdir,eventdir,region_list{roi}));
            roi_struct.(region_names{roi}) = fmri_data(roi_fname{1});
            
            % grab the data for that seed from each data structure
            curr_region_gain.(region_names{roi}) = extract_roi_averages(ppi_data_gain_all,roi_struct.(region_names{roi}));
            curr_region_gain0.(region_names{roi}) = extract_roi_averages(ppi_data_gain_all,roi_struct.(region_names{roi}));
            curr_region_loss.(region_names{roi}) = extract_roi_averages(ppi_data_gain_all,roi_struct.(region_names{roi}));
            curr_region_loss0.(region_names{roi}) = extract_roi_averages(ppi_data_gain_all,roi_struct.(region_names{roi}));
            
            % get correlation matrices for these seed to seed "contrasts"
            roi_ppi_corr_mat_gain = corrcoef(curr_region_gain.(region_names{roi}).all_data,curr_region_gain0.(region_names{roi}).all_data);
            roi_ppi_corr_mat_loss = corrcoef(curr_region_loss.(region_names{roi}).all_data,curr_region_loss0.(region_names{roi}).all_data);
            
            % generate roi disimilarity
            roi_ppi_disim_gain = 1 - roi_ppi_corr_mat_gain;
            roi_ppi_disim_loss = 1 - roi_ppi_corr_mat_loss;
            roi_ppi_disim_dat_gain{:,sub} = curr_region_gain.(region_names{roi}).all_data';
            roi_ppi_disim_dat_gain_neutral{:,sub} = curr_region_gain0.(region_names{roi}).all_data';
            roi_ppi_disim_dat_loss{:,sub} = curr_region_loss.(region_names{roi}).all_data';
            roi_ppi_disim_dat_loss_neutral{:,sub} = curr_region_loss0.(region_names{roi}).all_data';
            
            % store each roi disimilarity value for subs in a single vector
            PID(sub,1) = str2num(fnames_wholebrain_gain_all{sub}(1:5));
            roi_results_vector_gain(sub,1) = roi_ppi_disim_gain(2,1);
            roi_results_vector_loss(sub,1) = roi_ppi_disim_loss(2,1);
        end
    end         
end

% compare_actual_gain = mean(results_vector_gain);
% compare_actual_loss = mean(results_vector_loss);

% %% Generate null case to compare mean disimilarity to
% 
% % Only need to create a single bootstrapped distribution I think and will
% % plot average disimilarity on it. This is really just a sanity check to
% % make sure we're getting any kind of difference at the whole brain level.
% 
% for i = 1:iterations
%     index = randperm(size(ppi_disim_dat_gain,2));
%     bootstrapped_brain = ppi_disim_dat_gain(:,index);
%     bootstrapped_ppi_corr = corr2(bootstrapped_brain,ppi_disim_dat_gain_neutral);
%     bootstrapped_disim(i) = 1 - bootstrapped_ppi_corr;
%     disp(strcat('Number of bootstraps complete:',num2str(i)))
% end
% 
% figure(); histogram(bootstrapped_disim); title("Whole brain: Gain")          
% str = [' disim = ',num2str(compare_actual_gain)];
% T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str);
% set(T, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');
% 

%% load in clinical data and other relevant measures
load(fullfile(clinicaldir,'first_year_project_trilevel_T1.mat'));
load(fullfile(med_dir,'Medication_T2.mat'));
load(fullfile(demo_dir,'demographics.mat'));
load(fullfile(stressdir,'LSI_T1.mat'))

trilevel_array = [trilevel_T1.id,trilevel_T1.GenDis,trilevel_T1.Anhedon,trilevel_T1.Fears,trilevel_T1.narrow];
%trilevel_array = trilevel_array(:,2:size(trilevel_array,2));
med_array = [T2MedicationInventory3(:,1),T2MedicationInventory3(:,88)];
med_array = table2array(med_array);
dem_array = [BrainMAPDT1S1Demo.PID,BrainMAPDT1S1Demo.sex,BrainMAPDT1S1Demo.ethnicity,BrainMAPDT1S1Demo.race,BrainMAPDT1S1Demo.race2];
lsi_array = [LSI_T1(:,1),LSI_T1(:,6),LSI_T1(:,8:9),LSI_T1(:,11:13),LSI_T1(:,15),LSI_T1(:,17:19)];
lsi_array = table2array(lsi_array);

for sub = 1:length(fnames_wholebrain_gain_all)
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
        lsi_regressors(sub,1) = mean(lsi_array(curr2,2:11));
    else
        disp(strcat(num2str(PID(sub,1)), ' missing life stress info'))
        lsi_regressors(sub,1) = NaN;
    end
end
site_regressor = ones(length(PID),1); site_regressor(PID>=20000,1) = 0;
%% Analysis portion
if whole_brain_analysis == 1
    curr_analysis_table_symp = [trilevel_regressors,med_regressors,demographic_regressors,site_regressor,results_vector_gain,results_vector_loss];
    curr_analysis_table_stress = [lsi_regressors,med_regressors,demographic_regressors,results_vector_gain,site_regressor,results_vector_loss];

    curr_analysis_table_symp = array2table(curr_analysis_table_symp);curr_analysis_table_symp.Properties.VariableNames = {'PID','GenDis','Anhedonia','Fears','Narrow','Meds','Sex','Site','Disim_gain','Disim_loss'};
    curr_analysis_table_stress = array2table(curr_analysis_table_stress);curr_analysis_table_stress.Properties.VariableNames = {'Stress','Meds','Sex','Site','Disim_gain','Disim_loss'};

    fitlm(curr_analysis_table_symp,'Disim_gain ~ GenDis + Anhedonia + Fears + Narrow + Meds + Sex + Site')
    fitlm(curr_analysis_table_symp,'Disim_loss ~ GenDis + Anhedonia + Fears + Narrow + Meds + Sex + Site')

    fitlm(curr_analysis_table_stress,'Disim_gain ~ Stress + Meds + Sex + Site')
    fitlm(curr_analysis_table_stress,'Disim_loss ~ Stress + Meds + Sex + Site')
elseif roi_analysis == 1
    curr_analysis_table_symp = [trilevel_regressors,med_regressors,demographic_regressors,site_regressor,roi_results_vector_gain,roi_results_vector_loss];
    curr_analysis_table_stress = [lsi_regressors,med_regressors,demographic_regressors,roi_results_vector_gain,site_regressor,roi_results_vector_loss];

    curr_analysis_table_symp = array2table(curr_analysis_table_symp);curr_analysis_table_symp.Properties.VariableNames = {'PID','GenDis','Anhedonia','Fears','Narrow','Meds','Sex','Site','Disim_gain','Disim_loss'};
    curr_analysis_table_stress = array2table(curr_analysis_table_stress);curr_analysis_table_stress.Properties.VariableNames = {'Stress','Meds','Sex','Site','Disim_gain','Disim_loss'};

    fitlm(curr_analysis_table_symp,'Disim_gain ~ GenDis + Anhedonia + Fears + Narrow + Meds + Sex + Site')
    fitlm(curr_analysis_table_symp,'Disim_loss ~ GenDis + Anhedonia + Fears + Narrow + Meds + Sex + Site')

    fitlm(curr_analysis_table_stress,'Disim_gain ~ Stress + Meds + Sex + Site')
    fitlm(curr_analysis_table_stress,'Disim_loss ~ Stress + Meds + Sex + Site')
end
%% OLD CODE, MAY BE HANDY LATER
% %% Extended Bilateral VS - bilateral Amygdala vs Big three symptoms
% 
% ppi_data2 = [roi_avg_gain2.bi_VS_AntRew.HO_Amyg.all_data,roi_avg_loss2.bi_VS_AntLoss.HO_Amyg.all_data];
% ppi_corr_mat2 = corrcoef(ppi_data2');
% ppi_disim2 = 1 - ppi_corr_mat2;
% 
% symp_data = [GenDis, Fears, Anhedonia];
% symp_corr_mat1 = corrcoef(symp_data');
% symp_disim1 = 1 - symp_corr_mat1;
% 
% heatmap_input2 = tril(symp_disim1) + triu(ppi_disim2);
% 
% figure(); heatmap(heatmap_input2)
% saveas(gcf,fullfile(contrastdir,'heatmap2.jpg'))
% 
% lower_half_ppi1 = tril(ppi_disim2);
% lower_half_symp = tril(symp_disim1);
% 
% compare_actual = 1 - corr2(lower_half_ppi1,lower_half_symp);
% % Need to bootstrap to get a comparison condition
% 
% for i = 1:iterations
%     index = randperm(length(ppi_disim2));
%     bootstrapped_brain = ppi_data2(:,index);
%     bootstrapped_ppi_corr = corrcoef(bootstrapped_brain');
%     bootstrapped_disim = 1 - bootstrapped_ppi_corr;
%     lower_half_ppi_dist = tril(bootstrapped_disim);
%     bootstrapped_symp = symp_disim1(:,:);
%     lower_half_symp_dist = tril(bootstrapped_symp);
%     compare_dist(i) = 1 - corr2(lower_half_ppi_dist,lower_half_symp_dist);
% end
% 
% figure(); histogram(compare_dist); title("General Distress, Fears, Anhedonia: bilateral VS - bilateral Amygdala")          
% str = [' r = ',num2str(compare_actual)];
% T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str);
% set(T, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');
% 
% saveas(gcf,fullfile(contrastdir,'RSA2.jpg'))
%     
% %% mOFC - Ng Amyg vs all symptoms
% 
% ppi_data3 = [roi_avg_gain2.biOFC_to_wholebrain.Amyg_Ng.all_data,roi_avg_loss2.biOFC_to_wholebrain.Amyg_Ng.all_data];
% ppi_corr_mat3 = corrcoef(ppi_data3');
% ppi_disim3 = 1 - ppi_corr_mat3;
% 
% symp_data2 = trilevel_regressors;
% symp_corr_mat2 = corrcoef(symp_data2');
% symp_disim2 = 1 - symp_corr_mat2;
% 
% heatmap_input3 = tril(symp_disim2) + triu(ppi_disim3);
% 
% figure(); heatmap(heatmap_input3)
% saveas(gcf,fullfile(contrastdir,'heatmap3.jpg'))
% 
% lower_half_ppi3 = tril(ppi_disim3);
% lower_half_symp2 = tril(symp_disim2);
% 
% compare_actual = 1 - corr2(lower_half_ppi3,lower_half_symp2);
% % Need to bootstrap to get a comparison condition
% 
% for i = 1:iterations
%     index = randperm(length(ppi_disim3));
%     bootstrapped_brain = ppi_data3(:,index);
%     bootstrapped_ppi_corr = corrcoef(bootstrapped_brain');
%     bootstrapped_disim = 1 - bootstrapped_ppi_corr;
%     lower_half_ppi_dist = tril(bootstrapped_disim);
%     bootstrapped_symp = symp_disim2(:,:);
%     lower_half_symp_dist = tril(bootstrapped_symp);
%     compare_dist(i) = 1 - corr2(lower_half_ppi_dist,lower_half_symp_dist);
% end
% 
% figure(); histogram(compare_dist); title("All Symptoms: mOFC - Extended Amygdala")          
% str = [' r = ',num2str(compare_actual)];
% T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str);
% set(T, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');
% 
% saveas(gcf,fullfile(contrastdir,'RSA3.jpg'))
%     
% %% bilateral VS - bilateral Amyg vs all symptoms
% 
% ppi_data4 = [roi_avg_gain2.bi_VS_AntRew.HO_Amyg.all_data,roi_avg_loss2.bi_VS_AntLoss.HO_Amyg.all_data];
% ppi_corr_mat4 = corrcoef(ppi_data4');
% ppi_disim4 = 1 - ppi_corr_mat4;
% 
% symp_data2 = trilevel_regressors;
% symp_corr_mat2 = corrcoef(symp_data2');
% symp_disim2 = 1 - symp_corr_mat2;
% 
% heatmap_input4 = tril(symp_disim2) + triu(ppi_disim4);
% 
% figure(); heatmap(heatmap_input4)
% saveas(gcf,fullfile(contrastdir,'heatmap4.jpg'))
% 
% lower_half_ppi4 = tril(ppi_disim4);
% lower_half_symp2 = tril(symp_disim2);
% 
% compare_actual = 1 - corr2(lower_half_ppi4,lower_half_symp2);
% % Need to bootstrap to get a comparison condition
% 
% for i = 1:iterations
%     index = randperm(length(ppi_disim4));
%     bootstrapped_brain = ppi_data4(:,index);
%     bootstrapped_ppi_corr = corrcoef(bootstrapped_brain');
%     bootstrapped_disim = 1 - bootstrapped_ppi_corr;
%     lower_half_ppi_dist = tril(bootstrapped_disim);
%     bootstrapped_symp = symp_disim2(:,:);
%     lower_half_symp_dist = tril(bootstrapped_symp);
%     compare_dist(i) = 1 - corr2(lower_half_ppi_dist,lower_half_symp_dist);
% end
% 
% figure(); histogram(compare_dist); title("All symptoms: bilateral VS - bilateral Amygdala")          
% str = [' r = ',num2str(compare_actual)];
% T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str);
% set(T, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');
% 
% saveas(gcf,fullfile(contrastdir,'RSA4.jpg'))
%     
% 
% 
% 
% 
% 
