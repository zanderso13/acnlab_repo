sub_spm_dir = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD/func_conn/first_levels/sub-20309/ses-2/run-1/MID';
roi_dir = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD/func_conn/roi';
preproc_dir = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD/func_conn/fmriprep/sub-20309/ses-2/func';
save_dir = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD/func_conn/figures/PPI_task_x_timecourse_regressors';

run1 = fmri_data(filenames(fullfile(preproc_dir,'*MID*run-01*nii')));
roi = fmri_data(filenames(fullfile(roi_dir,'*OFC*nii')));
load(fullfile(sub_spm_dir,'SPM.mat'));

time_series = run1.dat;
roi_dat = roi.dat;
roi_time_series = time_series .* roi_dat;

mean_roi_bold = sum(roi_time_series,1) ./ 281;

std_mean_roi_bold = zscore(mean_roi_bold');

%%

roi_time_regressor = std_mean_roi_bold(5:283);
regressors = SPM.xX.X;


% convolve task*HRF with Bold time series

number_of_contrasts = 9;
con_names = {'Loss Anticipation','Loss Anticipation Zero','Gain Anticipation','Gain Anticipation Zero','Loss Neutral','Loss Consumption','Gain Consumption','Gain Neutral','Motor'};
con_names_for_save = {'Loss_Anticipation','Loss_Anticipation_Zero','Gain_Anticipation','Gain_Anticipation_Zero','Loss_Neutral','Loss_Consumption','Gain_Consumption','Gain_Neutral','Motor'};
for con = 1:number_of_contrasts
    con_temp = regressors(:,con) .* roi_time_regressor;
    
    figure();
    plot(con_temp)
    title(con_names{con})
    temp_file_name = strcat('Contrast_',con_names_for_save{con},'.jpg');
    saveas(gcf,fullfile(save_dir,temp_file_name))
end
close all



