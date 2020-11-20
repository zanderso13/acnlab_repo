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

datadir = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD/func_conn/RSA';
load(fullfile(datadir,'PPI_roi_data_for_RSA.mat'))
load(fullfile(datadir,'symptoms_complete.mat'))
iterations = 10000;
alpha = 0.001;

%% mOFC - Ng Amyg vs Big three symptoms

ppi_data1 = [roi_avg_gain2.biOFC_to_wholebrain.Amyg_Ng.all_data,roi_avg_loss2.biOFC_to_wholebrain.Amyg_Ng.all_data];
ppi_corr_mat1 = corrcoef(ppi_data1');
ppi_disim1 = 1 - ppi_corr_mat1;

symp_data = [GenDis, Fears, Anhedonia];
symp_corr_mat1 = corrcoef(symp_data');
symp_disim1 = 1 - symp_corr_mat1;

heatmap_input1 = tril(symp_disim1) + triu(ppi_disim1);

figure(); heatmap(heatmap_input1)
saveas(gcf,fullfile(datadir,'heatmap1.jpg'))

lower_half_ppi1 = tril(ppi_disim1);
lower_half_symp = tril(symp_disim1);

compare_actual = 1 - corr2(lower_half_ppi1,lower_half_symp);
% Need to bootstrap to get a comparison condition

for i = 1:iterations
    index = randperm(length(ppi_disim1));
    bootstrapped_brain = ppi_data1(:,index);
    bootstrapped_ppi_corr = corrcoef(bootstrapped_brain');
    bootstrapped_disim = 1 - bootstrapped_ppi_corr;
    lower_half_ppi_dist = tril(bootstrapped_disim);
    bootstrapped_symp = symp_disim1(:,:);
    lower_half_symp_dist = tril(bootstrapped_symp);
    compare_dist(i) = 1 - corr2(lower_half_ppi_dist,lower_half_symp_dist);
end

figure(); histogram(compare_dist); title("General Distress, Fears, Anhedonia: mOFC - Extended Amygdala")          
str = [' r = ',num2str(compare_actual)];
T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str);
set(T, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');

saveas(gcf,fullfile(datadir,'RSA1.jpg'))
%% mOFC - bilateral VS vs Big three symptoms

ppi_data2 = [roi_avg_gain2.biOFC_to_wholebrain.bi_vs_sphere.all_data,roi_avg_loss2.biOFC_to_wholebrain.bi_vs_sphere.all_data];
ppi_corr_mat2 = corrcoef(ppi_data2');
ppi_disim2 = 1 - ppi_corr_mat2;

symp_data = [GenDis, Fears, Anhedonia];
symp_corr_mat1 = corrcoef(symp_data');
symp_disim1 = 1 - symp_corr_mat1;

heatmap_input2 = tril(symp_disim1) + triu(ppi_disim2);

figure(); heatmap(heatmap_input2)
saveas(gcf,fullfile(datadir,'heatmap2.jpg'))

lower_half_ppi1 = tril(ppi_disim2);
lower_half_symp = tril(symp_disim1);

compare_actual = 1 - corr2(lower_half_ppi1,lower_half_symp);
% Need to bootstrap to get a comparison condition

for i = 1:iterations
    index = randperm(length(ppi_disim2));
    bootstrapped_brain = ppi_data2(:,index);
    bootstrapped_ppi_corr = corrcoef(bootstrapped_brain');
    bootstrapped_disim = 1 - bootstrapped_ppi_corr;
    lower_half_ppi_dist = tril(bootstrapped_disim);
    bootstrapped_symp = symp_disim1(:,:);
    lower_half_symp_dist = tril(bootstrapped_symp);
    compare_dist(i) = 1 - corr2(lower_half_ppi_dist,lower_half_symp_dist);
end

figure(); histogram(compare_dist); title("General Distress, Fears, Anhedonia: mOFC - bilateral VS")          
str = [' r = ',num2str(compare_actual)];
T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str);
set(T, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');

saveas(gcf,fullfile(datadir,'RSA2.jpg'))
    
%% mOFC - Ng Amyg vs all symptoms

ppi_data3 = [roi_avg_gain2.biOFC_to_wholebrain.Amyg_Ng.all_data,roi_avg_loss2.biOFC_to_wholebrain.Amyg_Ng.all_data];
ppi_corr_mat3 = corrcoef(ppi_data3');
ppi_disim3 = 1 - ppi_corr_mat1;

symp_data2 = trilevel_regressors;
symp_corr_mat2 = corrcoef(symp_data2');
symp_disim2 = 1 - symp_corr_mat2;

heatmap_input3 = tril(symp_disim2) + triu(ppi_disim3);

figure(); heatmap(heatmap_input3)
saveas(gcf,fullfile(datadir,'heatmap3.jpg'))

lower_half_ppi3 = tril(ppi_disim3);
lower_half_symp2 = tril(symp_disim2);

compare_actual = 1 - corr2(lower_half_ppi1,lower_half_symp2);
% Need to bootstrap to get a comparison condition

for i = 1:iterations
    index = randperm(length(ppi_disim3));
    bootstrapped_brain = ppi_data1(:,index);
    bootstrapped_ppi_corr = corrcoef(bootstrapped_brain');
    bootstrapped_disim = 1 - bootstrapped_ppi_corr;
    lower_half_ppi_dist = tril(bootstrapped_disim);
    bootstrapped_symp = symp_disim2(:,:);
    lower_half_symp_dist = tril(bootstrapped_symp);
    compare_dist(i) = 1 - corr2(lower_half_ppi_dist,lower_half_symp_dist);
end

figure(); histogram(compare_dist); title("All Symptoms: mOFC - Extended Amygdala")          
str = [' r = ',num2str(compare_actual)];
T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str);
set(T, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');

saveas(gcf,fullfile(datadir,'RSA3.jpg'))
    
%% bilateral VS - bilateral Amyg vs all symptoms

ppi_data4 = [roi_avg_gain2.bi_VS_AntRew.HO_Amyg.all_data,roi_avg_loss2.bi_VS_AntLoss.HO_Amyg.all_data];
ppi_corr_mat4 = corrcoef(ppi_data1');
ppi_disim4 = 1 - ppi_corr_mat1;

symp_data2 = trilevel_regressors;
symp_corr_mat2 = corrcoef(symp_data2');
symp_disim2 = 1 - symp_corr_mat2;

heatmap_input4 = tril(symp_disim2) + triu(ppi_disim4);

figure(); heatmap(heatmap_input4)
saveas(gcf,fullfile(datadir,'heatmap4.jpg'))

lower_half_ppi4 = tril(ppi_disim4);
lower_half_symp2 = tril(symp_disim2);

compare_actual = 1 - corr2(lower_half_ppi4,lower_half_symp2);
% Need to bootstrap to get a comparison condition

for i = 1:iterations
    index = randperm(length(ppi_disim4));
    bootstrapped_brain = ppi_data1(:,index);
    bootstrapped_ppi_corr = corrcoef(bootstrapped_brain');
    bootstrapped_disim = 1 - bootstrapped_ppi_corr;
    lower_half_ppi_dist = tril(bootstrapped_disim);
    bootstrapped_symp = symp_disim2(:,:);
    lower_half_symp_dist = tril(bootstrapped_symp);
    compare_dist(i) = 1 - corr2(lower_half_ppi_dist,lower_half_symp_dist);
end

figure(); histogram(compare_dist); title("General Distress, Fears, Anhedonia: mOFC - bilateral VS")          
str = [' r = ',num2str(compare_actual)];
T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str);
set(T, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');

saveas(gcf,fullfile(datadir,'RSA4.jpg'))
    





