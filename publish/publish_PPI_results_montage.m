function publish_PPI_results_montage(curr_ROI,curr_title,analysisdir,condir)

%% PPI results today('datetime')
% Current results are generated using a 0.001 primary threshold and a
% cluster threshold of 10 voxels. This is a very liberal cluster extent
% thresholding procedure and we might think of making this more stringent
% if we choose to publish

% Results are generated using a robust regression model. Downweighting
% outliers in this case will reduce issues related to motion related spikes

% First level models includes 12 motion regressors (transverse, rotation,
% and all their derivatives) in addition to a spike detection scheme
% generated using CanlabCore tools. Images are similarly generated using
% CanlabCore tools.
currdate = today('datetime')
% figuredir = 'figures';
% savedir = fullfile(analysisdir,condir,figuredir);

curr_analysis_fname = filenames(fullfile(analysisdir,condir,strcat('longitudinal*',curr_ROI,'*')));
load(curr_analysis_fname{1})
for gain_or_loss = 1:2
    if gain_or_loss == 1
        % First timepoint, gain
        figure();
        montage(obj_struct.gain.life_dep1_dat_gain)
        title(strcat(curr_title,' T1 depression: gain contrast'))
        %saveas(gcf,fullfile(savedir,strcat(curr_analysis{roi},'T1_dep_gain.jpg')));
        figure();
        montage(obj_struct.gain.life_anx1_dat_gain)
        title(strcat(curr_title,' T1 anxiety: gain contrast')) 
        %saveas(gcf,fullfile(savedir,strcat(curr_analysis{roi},'T1_anx_gain.jpg')));
        figure();
        montage(obj_struct.gain.life_dep1_dat_gain)
        title(strcat(curr_title,' T1 comorbid: gain contrast'))
        %saveas(gcf,fullfile(savedir,strcat(curr_analysis{roi},'T1_com_gain.jpg')));
        % Second timepoint, gain
        figure();
        montage(obj_struct.gain.life_dep2_dat_gain)
        title(strcat(curr_title,' T2 depression: gain contrast'))
        %saveas(gcf,fullfile(savedir,strcat(curr_analysis{roi},'T2_dep_gain.jpg')));
        figure();
        montage(obj_struct.gain.life_anx2_dat_gain)
        title(strcat(curr_title,' T2 anxiety: gain contrast'))
        %saveas(gcf,fullfile(savedir,strcat(curr_analysis{roi},'T2_anx_gain.jpg')));
        figure();
        montage(obj_struct.gain.life_dep2_dat_gain)
        title(strcat(curr_title,' T2 comorbid: gain contrast'))
        %saveas(gcf,fullfile(savedir,strcat(curr_analysis{roi},'T2_com_gain.jpg')));
    else
        % First timepoint, loss 
        figure();
        montage(obj_struct.loss.life_dep1_dat_loss)
        title(strcat(curr_title,' T1 depression: loss contrast'))
        %saveas(gcf,fullfile(savedir,strcat(curr_analysis{roi},'T1_dep_loss.jpg')));
        figure();
        montage(obj_struct.loss.life_anx1_dat_loss)
        title(strcat(curr_title,' T1 anxiety: loss contrast'))
        %saveas(gcf,fullfile(savedir,strcat(curr_analysis{roi},'T1_anx_loss.jpg')));
        montage(obj_struct.loss.life_dep1_dat_loss)
        title(strcat(curr_title,' T1 comorbid: loss contrast'))
        %saveas(gcf,fullfile(savedir,strcat(curr_analysis{roi},'T1_com_loss.jpg')));
        % Second timepoint, loss
        figure();
        montage(obj_struct.loss.life_dep2_dat_loss)
        title(strcat(curr_title,' T2 depression: loss contrast'))
        %saveas(gcf,fullfile(savedir,strcat(curr_analysis{roi},'T2_dep_loss.jpg')));
        figure();
        montage(obj_struct.loss.life_anx2_dat_loss)
        title(strcat(curr_title,' T2 anxiety: loss contrast'))
        %saveas(gcf,fullfile(savedir,strcat(curr_analysis{roi},'T2_anx_loss.jpg')));
        figure();
        montage(obj_struct.loss.life_dep2_dat_loss)
        title(strcat(curr_title,' T2 comorbid: loss contrast'))
        %saveas(gcf,fullfile(savedir,strcat(curr_analysis{roi},'T2_com_loss.jpg')));
    end
end

end