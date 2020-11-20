% Generate a model to plot with the scatter
savedir = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD/func_conn/PPI/results/Trilevel_analysis/figures';
%% Gain Anticipation
% biOFC - biHOAmyg
resid_mdl_R = [R_final_gain2(:,1),R_final_gain2(:,3:size(R_final_gain2,2))];
resid_mdl = fitlm(resid_mdl_R,roi_avg_gain.biOFC_to_wholebrain.HO_Amyg.dat);
plot_mdl = fitlm(Anhedonia,resid_mdl.Residuals.Raw)
figure();scatter(Anhedonia,resid_mdl.Residuals.Raw)
h1 = lsline();
h1.LineWidth = 2;
h1.Color = 'r';
r1 = corrcoef(Anhedonia(:,1),resid_mdl.Residuals.Raw,'rows','complete');
disp(r1(1,2));
str = [' r = ',num2str(r1(1,2))];
T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str);
set(T, 'fontsize', 20, 'verticalalignment', 'top', 'horizontalalignment', 'left');
hold on; plot(plot_mdl); title("biOFC - HO Amyg");xlabel("");ylabel("");legend off
drawnow, snapnow
saveas(gcf,fullfile(savedir,'Residual_Anhedonia_Gain_Anticipation1.jpg'))

% biOFC - Ng Amyg
resid_mdl_R = [R_final_gain2(:,1),R_final_gain2(:,3:size(R_final_gain2,2))];
resid_mdl = fitlm(resid_mdl_R,roi_avg_gain.biOFC_to_wholebrain.Amyg_Ng.dat);
plot_mdl = fitlm(Anhedonia,resid_mdl.Residuals.Raw)
figure();scatter(Anhedonia,resid_mdl.Residuals.Raw)
h1 = lsline();
h1.LineWidth = 2;
h1.Color = 'r';
r1 = corrcoef(Anhedonia(:,1),resid_mdl.Residuals.Raw,'rows','complete');
disp(r1(1,2));
str = [' r = ',num2str(r1(1,2))];
T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str);
set(T, 'fontsize', 20, 'verticalalignment', 'top', 'horizontalalignment', 'left');
hold on; plot(plot_mdl); title("");xlabel("");ylabel("");legend off
drawnow, snapnow
saveas(gcf,fullfile(savedir,'Residual_Anhedonia_Gain_Anticipation1.jpg'))


% biVS - biHOAmyg

resid_mdl_R = [R_final_gain2(:,1),R_final_gain2(:,3:size(R_final_gain2,2))];
resid_mdl = fitlm(resid_mdl_R,roi_avg_gain.bi_VS_AntLoss.HO_Amyg.dat);
plot_mdl = fitlm(Anhedonia,resid_mdl.Residuals.Raw)
figure();scatter(Anhedonia,resid_mdl.Residuals.Raw)
h1 = lsline();
h1.LineWidth = 2;
h1.Color = 'r';
r1 = corrcoef(Anhedonia(:,1),resid_mdl.Residuals.Raw,'rows','complete');
disp(r1(1,2));
str = [' r = ',num2str(r1(1,2))];
T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str);
set(T, 'fontsize', 20, 'verticalalignment', 'top', 'horizontalalignment', 'left');
hold on; plot(plot_mdl); title("");xlabel("");ylabel("");legend off
drawnow, snapnow
saveas(gcf,fullfile(savedir,'Residual_Anhedonia_Gain_Anticipation2.jpg'))

% biVS - biOFC

resid_mdl_R = [R_final_gain2(:,1),R_final_gain2(:,3:size(R_final_gain2,2))];
resid_mdl = fitlm(resid_mdl_R,roi_avg_gain.bi_VS_AntRew.HO_VMPFC.dat);
plot_mdl = fitlm(Anhedonia,resid_mdl.Residuals.Raw)
figure();scatter(Anhedonia,resid_mdl.Residuals.Raw)
h1 = lsline();
h1.LineWidth = 2;
h1.Color = 'r';
r1 = corrcoef(Anhedonia(:,1),resid_mdl.Residuals.Raw,'rows','complete');
disp(r1(1,2));
str = [' r = ',num2str(r1(1,2))];
T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str);
set(T, 'fontsize', 20, 'verticalalignment', 'top', 'horizontalalignment', 'left');
hold on; plot(plot_mdl); title("bi OFC - bi VS");xlabel("");ylabel("");legend off
drawnow, snapnow
saveas(gcf,fullfile(savedir,'Residual_Anhedonia_Gain_Anticipation3.jpg'))

%% Gain Consumption
% Anhedonia
% biOFC - biHOAmyg
resid_mdl_R = [R_final_gain2(:,1),R_final_gain2(:,3:size(R_final_gain2,2))];
resid_mdl = fitlm(resid_mdl_R,roi_avg_gain.biOFC_to_wholebrain.HO_Amyg.dat);
plot_mdl = fitlm(Anhedonia,resid_mdl.Residuals.Raw)
figure();scatter(Anhedonia,resid_mdl.Residuals.Raw)
h1 = lsline();
h1.LineWidth = 2;
h1.Color = 'r';
r1 = corrcoef(Anhedonia(:,1),resid_mdl.Residuals.Raw,'rows','complete');
disp(r1(1,2));
str = [' r = ',num2str(r1(1,2))];
T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str);
set(T, 'fontsize', 20, 'verticalalignment', 'top', 'horizontalalignment', 'left');
hold on; plot(plot_mdl); title("");xlabel("");ylabel("");legend off
drawnow, snapnow
saveas(gcf,fullfile(savedir,'Residual_Anhedonia_Gain_Consumption1.jpg'))

% biVS - Ng Amyg

resid_mdl_R = [R_final_gain2(:,1),R_final_gain2(:,3:size(R_final_gain2,2))];
resid_mdl = fitlm(resid_mdl_R,roi_avg_gain.bi_VS_Oldham.Amyg_Ng.dat);
plot_mdl = fitlm(Anhedonia,resid_mdl.Residuals.Raw)
figure();scatter(Anhedonia,resid_mdl.Residuals.Raw)
h1 = lsline();
h1.LineWidth = 2;
h1.Color = 'r';
r1 = corrcoef(Anhedonia(:,1),resid_mdl.Residuals.Raw,'rows','complete');
disp(r1(1,2));
str = [' r = ',num2str(r1(1,2))];
T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str);
set(T, 'fontsize', 20, 'verticalalignment', 'top', 'horizontalalignment', 'left');
hold on; plot(plot_mdl); title("");xlabel("");ylabel("");legend off
drawnow, snapnow
saveas(gcf,fullfile(savedir,'Residual_Anhedonia_Gain_Consumption2.jpg'))

% VS - HO Amyg

resid_mdl_R = [R_final_gain2(:,1),R_final_gain2(:,3:size(R_final_gain2,2))];
resid_mdl = fitlm(resid_mdl_R,roi_avg_gain.bi_VS_Oldham.HO_Amyg.dat);
plot_mdl = fitlm(Anhedonia,resid_mdl.Residuals.Raw)
figure();scatter(Anhedonia,resid_mdl.Residuals.Raw)
h1 = lsline();
h1.LineWidth = 2;
h1.Color = 'r';
r1 = corrcoef(Anhedonia(:,1),resid_mdl.Residuals.Raw,'rows','complete');
disp(r1(1,2));
str = [' r = ',num2str(r1(1,2))];
T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str);
set(T, 'fontsize', 20, 'verticalalignment', 'top', 'horizontalalignment', 'left');
hold on; plot(plot_mdl); title("");xlabel("");ylabel("");legend off
drawnow, snapnow
saveas(gcf,fullfile(savedir,'Residual_Anhedonia_Gain_Consumption3.jpg'))

% bi OFC - VS

resid_mdl_R = [R_final_gain2(:,1),R_final_gain2(:,3:6)];
resid_mdl = fitlm(resid_mdl_R,roi_avg_gain.biOFC_to_wholebrain.bi_VS_Oldham.dat);
plot_mdl = fitlm(Anhedonia,resid_mdl.Residuals.Raw)
figure();scatter(Anhedonia,resid_mdl.Residuals.Raw)
h1 = lsline();
h1.LineWidth = 2;
h1.Color = 'r';
r1 = corrcoef(Anhedonia(:,1),resid_mdl.Residuals.Raw,'rows','complete');
disp(r1(1,2));
str = [' r = ',num2str(r1(1,2))];
T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str);
set(T, 'fontsize', 20, 'verticalalignment', 'top', 'horizontalalignment', 'left');
hold on; plot(plot_mdl); title("");xlabel("");ylabel("");legend off
drawnow, snapnow
saveas(gcf,fullfile(savedir,'Residual_Anhedonia_Gain_Consumption4.jpg'))


%% General Distress

% biVS - biHOvmPFC

resid_mdl_R = [R_final_gain2(:,2:size(R_final_gain2,2))];
resid_mdl = fitlm(resid_mdl_R,roi_avg_gain.bi_VS_Oldham.HO_VMPFC.dat);
plot_mdl = fitlm(GenDis,resid_mdl.Residuals.Raw)
figure();scatter(GenDis,resid_mdl.Residuals.Raw)
h1 = lsline();
h1.LineWidth = 2;
h1.Color = 'r';
r1 = corrcoef(GenDis(:,1),resid_mdl.Residuals.Raw,'rows','complete');
disp(r1(1,2));
str = [' r = ',num2str(r1(1,2))];
T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str);
set(T, 'fontsize', 20, 'verticalalignment', 'top', 'horizontalalignment', 'left');
hold on; plot(plot_mdl); title("bi VS - bi vmPFC");xlabel("");ylabel("");legend off
drawnow, snapnow
saveas(gcf,fullfile(savedir,'Residual_GenDis_Gain_Consumption1.jpg'))

% L OFC - biVS
resid_mdl_R = [R_final_gain2(:,2:size(R_final_gain2,2))];
resid_mdl = fitlm(resid_mdl_R,roi_avg_gain.LOFC_to_wholebrain.bi_VS_Oldham.dat);
plot_mdl = fitlm(GenDis,resid_mdl.Residuals.Raw)
figure();scatter(GenDis,resid_mdl.Residuals.Raw)
h1 = lsline();
h1.LineWidth = 2;
h1.Color = 'r';
r1 = corrcoef(GenDis(:,1),resid_mdl.Residuals.Raw,'rows','complete');
disp(r1(1,2));
str = [' r = ',num2str(r1(1,2))];
T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str);
set(T, 'fontsize', 20, 'verticalalignment', 'top', 'horizontalalignment', 'left');
hold on; plot(plot_mdl); title("");xlabel("");ylabel("");legend off
drawnow, snapnow
saveas(gcf,fullfile(savedir,'Residual_GenDis_Gain_Consumption2.jpg'))

%% Loss Consumption
% Fears
% biVS - Ng Amyg

resid_mdl_R = [R_final_loss2(:,1:2),R_final_loss2(:,4:size(R_final_gain2,2))];
resid_mdl = fitlm(resid_mdl_R,roi_avg_loss.biOFC_to_wholebrain.bi_VS_Oldham.dat);
plot_mdl = fitlm(Fears,resid_mdl.Residuals.Raw)
figure();scatter(Fears,resid_mdl.Residuals.Raw)
h1 = lsline();
h1.LineWidth = 2;
h1.Color = 'r';
r1 = corrcoef(Fears(:,1),resid_mdl.Residuals.Raw,'rows','complete');
disp(r1(1,2));
str = [' r = ',num2str(r1(1,2))];
T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str);
set(T, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');
hold on; plot(plot_mdl); title("");xlabel("");ylabel("");legend off
drawnow, snapnow
saveas(gcf,fullfile(savedir,'Residual_Fears_Loss_Consumption1.jpg'))

% L OFC2 - biVS

resid_mdl_R = [R_final_loss2(:,1:2),R_final_loss2(:,4:size(R_final_gain2,2))];
resid_mdl = fitlm(resid_mdl_R,roi_avg_loss.LOFC2_to_wholebrain.bi_VS_Oldham.dat);
plot_mdl = fitlm(Fears,resid_mdl.Residuals.Raw)
figure();scatter(Fears,resid_mdl.Residuals.Raw)
h1 = lsline();
h1.LineWidth = 2;
h1.Color = 'r';
r1 = corrcoef(Fears(:,1),resid_mdl.Residuals.Raw,'rows','complete');
disp(r1(1,2));
str = [' r = ',num2str(r1(1,2))];
T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str);
set(T, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');
hold on; plot(plot_mdl); title("");xlabel("");ylabel("");legend off
drawnow, snapnow
saveas(gcf,fullfile(savedir,'Residual_Fears_Loss_Consumption2.jpg'))