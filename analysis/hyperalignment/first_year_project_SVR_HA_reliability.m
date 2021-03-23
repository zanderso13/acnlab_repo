% With the replication movement in mind. I'm going to use SVR as suggested
% by (Marek et al, 2020). I'll do this on the specific seed to seed output
% from my PPI analysis. Currently, this seed to seed approach is showing a
% weak association with clinical symptoms. This is going to be my cross
% validated stress test.

clear dat ha_input ha_output corr_mat_unaligned corr_mat_aligned
dat = [roi_avg_gain2.biOFC_to_wholebrain.HO_Amyg.all_data,roi_avg_gain2.bi_VS_Oldham.Amyg_Ng.all_data];
% cols are subjects
dat = dat';

clear sub
for sub = 1:size(dat,2)
    ha_input{sub} = dat(:,sub);
end

[aligned,transforms] = hyperalign(ha_input{:});
clear sub
for sub = 1:length(aligned)
    ha_output(sub,:) = aligned{sub};
end

corr_mat_unaligned = corrcoef(dat(:,:));
corr_mat_aligned = corrcoef(ha_output(:,:));

heatmap(corr_mat_unaligned)
figure();heatmap(corr_mat_aligned)


%% SVR test unaligned data

svm_in_X = [dat'];
svm_in_Y = Anhedonia;

svm_mdl = fitrlinear(svm_in_X,svm_in_Y,'KFold',10,'Regularization','lasso');

% Once you have a cross validated model, predict responses based on that model

%yfit = predict(svm_mdl);

% RMSE

mse_unaligned = sqrt(kfoldLoss(svm_mdl));


%% SVR test aligned data
clear svm_mdl
svm_in_X = [ha_output];
svm_in_Y = Anhedonia;

svm_mdl = fitrlinear(svm_in_X,svm_in_Y,'KFold',10,'Regularization','lasso');

% Once you have a cross validated model, predict responses based on that model

%yfit = predict(svm_mdl);

% RMSE

mse_aligned = sqrt(kfoldLoss(svm_mdl));

%% Generate bootstrapped data
% and generate mse with each bootstrap to demonstrate significance of above
% mse finding






