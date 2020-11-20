% specs for whole brain analysis

primary_threshold = 0.001;
cluster_size = 25;
analysis_string = ' consumption';
do_all = 1; % Check out line 15 for this. I don't have sig results but I'd like to get all the plots that will at least tell a story into my presentation
publish_region_name_list_for_struct = [region_name_list_for_struct, R_region_name_list_for_struct, L_region_name_list_for_struct];
for r = 1:length(publish_region_name_list_for_struct)
% So basically we're looping through every seed region, and
% then kind of every seed region again... But one refers to the
% seed and the next refers to the region of interest extracted
% for results. Then we'll check the p value associated with the
% relationship and will pause to explore things further if
% there's a significant relationship. I hope this works.
    test_gain = gain_mdl.(publish_region_name_list_for_struct{r});
    if test_gain.Coefficients.pValue(2) < 0.05 || test_gain.Coefficients.pValue(3) < 0.05 || test_gain.Coefficients.pValue(4) < 0.05
    % This is going to create a results structure that will be
    % helpful for me to try and parse all the results I have.

    % Robin would like me to include graphs that show the
    % specific relationships between connectivity output. That
    % will probably work best if I pull out average estimates
    % for the regions here. Then I'll create specific plots in
    % the following if statement


        if test_gain.Coefficients.pValue(2) < 0.05 
            final_seed_to_seed_results.gain.GenDis.(publish_region_name_list_for_struct{r}) = test_gain;
            disp(strcat('During Gain ',analysis_string))
            disp('Symptom: General Distress')
            disp(strcat('ROI: ', publish_region_name_list_for_struct{r}))
            disp(strcat('t = ',num2str(test_gain.Coefficients.tStat(2)),'p = ',num2str(test_gain.Coefficients.pValue(2)))) 
            mdl = fitlm(GenDis,roi_avg_gain.(publish_region_name_list_for_struct{r}).dat);
            figure(); scatter(GenDis,roi_avg_gain.(publish_region_name_list_for_struct{r}).dat)
            title('Gain Anticipation')
            h1 = lsline();
            h1.LineWidth = 2;
            h1.Color = 'r';
            r1 = corrcoef(GenDis(:,1),roi_avg_gain.(publish_region_name_list_for_struct{r}).dat(:,1),'rows','complete');
            disp(r1(1,2));
            str = [' r = ',num2str(r1(1,2))];
            T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str);
            set(T, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');
            hold on; plot(mdl); legend off
            drawnow, snapnow
            clear mdl
            disp('________________________________________________________________')
        elseif test_gain.Coefficients.pValue(3) < 0.05
            final_seed_to_seed_results.gain.Anhedonia.(publish_region_name_list_for_struct{r}) = test_gain;
            disp(strcat('During Gain ',analysis_string))
            disp('Symptom: Anhedonia')
            disp(strcat('ROI: ', publish_region_name_list_for_struct{r}))
            disp(strcat('t = ',num2str(test_gain.Coefficients.tStat(3)),'p = ',num2str(test_gain.Coefficients.pValue(3))))
            mdl = fitlm(Anhedonia,roi_avg_gain.(publish_region_name_list_for_struct{r}).dat);
            figure(); scatter(Anhedonia,roi_avg_gain.(publish_region_name_list_for_struct{r}).dat);
            title('Gain Anticipation')
            h1 = lsline();
            h1.LineWidth = 2;
            h1.Color = 'r';
            r1 = corrcoef(Anhedonia(:,1),roi_avg_gain.(publish_region_name_list_for_struct{r}).dat(:,1),'rows','complete');
            disp(r1(1,2));
            str = [' r = ',num2str(r1(1,2))]
            T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str);
            set(T, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');
            hold on; plot(mdl); legend off
            drawnow, snapnow
            clear mdl
            disp('________________________________________________________________')
        elseif test_gain.Coefficients.pValue(4) < 0.05 
            final_seed_to_seed_results.gain.Fears.(publish_region_name_list_for_struct{r}) = test_gain;
            disp(strcat('During Gain ',analysis_string)) 
            disp('Symptom: Fears')
            disp(strcat('ROI: ', publish_region_name_list_for_struct{r}))
            disp(strcat('t = ',num2str(test_gain.Coefficients.tStat(4)),'p = ',num2str(test_gain.Coefficients.pValue(4))))
            mdl = fitlm(Fears,roi_avg_gain.(publish_region_name_list_for_struct{r}).dat);
            figure(); scatter(Fears,roi_avg_gain.(publish_region_name_list_for_struct{r}).dat);
            title('Gain Anticipation')
            h1 = lsline();
            h1.LineWidth = 2;
            h1.Color = 'r';
            r1 = corrcoef(Fears(:,1),roi_avg_gain.(publish_region_name_list_for_struct{r}).dat(:,1),'rows','complete');
            disp(r1(1,2));
            str = [' r = ',num2str(r1(1,2))]
            T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str);
            set(T, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');
            hold on; plot(mdl); legend off
            drawnow, snapnow
            clear mdl
            disp('________________________________________________________________')

        end
    end

    test_loss = loss_mdl.(publish_region_name_list_for_struct{r});
    if test_loss.Coefficients.pValue(2) < 0.05 || test_loss.Coefficients.pValue(3) < 0.05 || test_loss.Coefficients.pValue(4) < 0.05 

        if test_loss.Coefficients.pValue(2) < 0.05 
            final_seed_to_seed_results.loss.GenDis.(publish_region_name_list_for_struct{r}) = test_loss;
            disp(strcat('During Loss ',analysis_string))
            disp('Symptom: General Distress')
            disp(strcat('ROI: ', publish_region_name_list_for_struct{r}))
            disp(strcat('t = ',num2str(test_loss.Coefficients.tStat(2)),'p = ',num2str(test_loss.Coefficients.pValue(2))))
            mdl = fitlm(GenDis,roi_avg_loss.(publish_region_name_list_for_struct{r}).dat);
            figure(); scatter(GenDis,roi_avg_loss.(publish_region_name_list_for_struct{r}).dat);
            title('Loss Anticipation')
            h1 = lsline();
            h1.LineWidth = 2;
            h1.Color = 'r';
            r1 = corrcoef(GenDis(:,1),roi_avg_loss.(publish_region_name_list_for_struct{r}).dat(:,1),'rows','complete');
            disp(r1(1,2));
            str = [' r = ',num2str(r1(1,2))]
            T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str);
            set(T, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');
            hold on; plot(mdl); legend off
            drawnow, snapnow
            clear mdl
            disp('________________________________________________________________')
            clear mdl

        elseif test_loss.Coefficients.pValue(3) < 0.05 
            final_seed_to_seed_results.loss.Anhedonia.(publish_region_name_list_for_struct{r}) = test_loss;
            disp(strcat('During Loss ',analysis_string))
            disp('Symptom: Anhedonia')
            disp(strcat('ROI: ', publish_region_name_list_for_struct{r}))
            disp(strcat('t = ',num2str(test_loss.Coefficients.tStat(3)),'p = ',num2str(test_loss.Coefficients.pValue(3))))
            mdl = fitlm(Anhedonia,roi_avg_loss.(publish_region_name_list_for_struct{r}).dat);
            figure(); scatter(Anhedonia,roi_avg_loss.(publish_region_name_list_for_struct{r}).dat);
            title('Loss Anticipation')
            h1 = lsline();
            h1.LineWidth = 2;
            h1.Color = 'r';
            r1 = corrcoef(Anhedonia(:,1),roi_avg_loss.(publish_region_name_list_for_struct{r}).dat(:,1),'rows','complete');
            disp(r1(1,2));
            str = [' r = ',num2str(r1(1,2))]
            T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str);
            set(T, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');
            hold on; plot(mdl); legend off
            drawnow, snapnow
            clear mdl
            disp('________________________________________________________________')

        elseif test_loss.Coefficients.pValue(4) < 0.05
            final_seed_to_seed_results.loss.Fears.(publish_region_name_list_for_struct{r}) = test_loss;
            disp(strcat('During Loss ',analysis_string))
            disp('Symptom: Fears')
            disp(strcat('ROI: ', publish_region_name_list_for_struct{r}))
            disp(strcat('t = ',num2str(test_loss.Coefficients.tStat(4)),'p = ',num2str(test_loss.Coefficients.pValue(4))))
            mdl = fitlm(Fears,roi_avg_loss.(publish_region_name_list_for_struct{r}).dat);
            figure(); scatter(Fears,roi_avg_loss.(publish_region_name_list_for_struct{r}).dat);
            title('Loss Anticipation')
            h1 = lsline();
            h1.LineWidth = 2;
            h1.Color = 'r';
            r1 = corrcoef(Fears(:,1),roi_avg_loss.(publish_region_name_list_for_struct{r}).dat(:,1),'rows','complete');
            disp(r1(1,2));
            str = [' r = ',num2str(r1(1,2))]
            T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str);
            set(T, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');
            hold on; plot(mdl); legend off
            drawnow, snapnow
            clear mdl
            disp('________________________________________________________________')

        end
    end
end


clear test_gain test_loss
disp('STARTING SEED TO SEED ANALYSIS: CONTROL FOR OPPOSITE CONDITION')
disp('_________________________________________________________________________')
for r = 1:length(publish_region_name_list_for_struct)
    % report region of interest analysis

    % So basically we're looping through every seed region, and
    % then kind of every seed region again... But one refers to the
    % seed and the next refers to the region of interest extracted
    % for results. Then we'll check the p value associated with the
    % relationship and will pause to explore things further if
    % there's a significant relationship. I hope this works.
    test_gain = gain_mdl2.(publish_region_name_list_for_struct{r});
    if test_gain.Coefficients.pValue(2) < 0.05 || test_gain.Coefficients.pValue(3) < 0.05 || test_gain.Coefficients.pValue(4) < 0.05
        % This is going to create a results structure that will be
        % helpful for me to try and parse all the results I have.
        % Controlling for sex really opened the box
        if test_gain.Coefficients.pValue(2) < 0.05
            final_seed_to_seed_results.gain.GenDis.(publish_region_name_list_for_struct{r}) = test_gain;
            disp(strcat('During Gain ',analysis_string))
            disp('Symptom: General Distress')
            disp(strcat('ROI: ', publish_region_name_list_for_struct{r}))
            disp(strcat('t = ',num2str(test_gain.Coefficients.tStat(2)),'p = ',num2str(test_gain.Coefficients.pValue(2))))
            mdl = fitlm(GenDis,roi_avg_gain.(publish_region_name_list_for_struct{r}).dat);
            figure(); scatter(GenDis,roi_avg_gain.(publish_region_name_list_for_struct{r}).dat)
            title('Gain Anticipation')
            h1 = lsline();
            h1.LineWidth = 2;
            h1.Color = 'r';
            r1 = corrcoef(GenDis(:,1),roi_avg_gain.(publish_region_name_list_for_struct{r}).dat(:,1),'rows','complete');
            disp(r1(1,2));
            str = [' r = ',num2str(r1(1,2))]
            T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str);
            set(T, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');
            hold on; plot(mdl); legend off
            drawnow, snapnow
            disp('________________________________________________________________')
        elseif test_gain.Coefficients.pValue(3) < 0.05
            final_seed_to_seed_results.gain.Anhedonia.(publish_region_name_list_for_struct{r}) = test_gain;
            disp(strcat('During Gain ',analysis_string))
            disp('Symptom: Anhedonia')
            disp(strcat('ROI: ', publish_region_name_list_for_struct{r}))
            disp(strcat('t = ',num2str(test_gain.Coefficients.tStat(3)),'p = ',num2str(test_gain.Coefficients.pValue(3))))
            mdl = fitlm(Anhedonia,roi_avg_gain.(publish_region_name_list_for_struct{r}).dat);
            figure(); scatter(Anhedonia,roi_avg_gain.(publish_region_name_list_for_struct{r}).dat)
            title('Gain Anticipation')
            h1 = lsline();
            h1.LineWidth = 2;
            h1.Color = 'r';
            r1 = corrcoef(Anhedonia(:,1),roi_avg_gain.(publish_region_name_list_for_struct{r}).dat(:,1),'rows','complete');
            disp(r1(1,2));
            str = [' r = ',num2str(r1(1,2))]
            T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str);
            set(T, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');
            hold on; plot(mdl); legend off
            drawnow, snapnow
            clear mdl
            disp('________________________________________________________________')

        elseif test_gain.Coefficients.pValue(4) < 0.05
            final_seed_to_seed_results.gain.Fears.(publish_region_name_list_for_struct{r}) = test_gain;
            disp(strcat('During Gain ',analysis_string)) 
            disp('Symptom: Fears')
            disp(strcat('ROI: ', publish_region_name_list_for_struct{r}))
            disp(strcat('t = ',num2str(test_gain.Coefficients.tStat(4)),'p = ',num2str(test_gain.Coefficients.pValue(4))))
            mdl = fitlm(Fears,roi_avg_gain.(publish_region_name_list_for_struct{r}).dat);
            figure(); scatter(Fears,roi_avg_gain.(publish_region_name_list_for_struct{r}).dat)
            title('Gain Anticipation')
            h1 = lsline();
            h1.LineWidth = 2;
            h1.Color = 'r';
            r1 = corrcoef(Fears(:,1),roi_avg_gain.(publish_region_name_list_for_struct{r}).dat(:,1),'rows','complete');
            disp(r1(1,2));
            str = [' r = ',num2str(r1(1,2))]
            T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str);
            set(T, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');
            hold on; plot(mdl); legend off
            drawnow, snapnow
            clear mdl
            disp('________________________________________________________________')

        end
    end
    test_loss = loss_mdl2.(publish_region_name_list_for_struct{r});
    if test_loss.Coefficients.pValue(2) < 0.05 || test_loss.Coefficients.pValue(3) < 0.05 || test_loss.Coefficients.pValue(4) < 0.05
        if test_loss.Coefficients.pValue(2) < 0.05
            final_seed_to_seed_results.loss.GenDis.(publish_region_name_list_for_struct{r}) = test_loss;
            disp(strcat('During Loss ',analysis_string))
            disp('Symptom: General Distress')
            disp(strcat('ROI: ', publish_region_name_list_for_struct{r}))
            disp(strcat('t = ',num2str(test_loss.Coefficients.tStat(2)),'p = ',num2str(test_loss.Coefficients.pValue(2))))
            mdl = fitlm(GenDis,roi_avg_loss.(publish_region_name_list_for_struct{r}).dat);
            figure(); scatter(GenDis,roi_avg_loss.(publish_region_name_list_for_struct{r}).dat)
            title('Loss Anticipation')
            h1 = lsline();
            h1.LineWidth = 2;
            h1.Color = 'r';
            r1 = corrcoef(GenDis(:,1),roi_avg_loss.(publish_region_name_list_for_struct{r}).dat(:,1),'rows','complete');
            disp(r1(1,2));
            str = [' r = ',num2str(r1(1,2))]
            T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str);
            set(T, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');
            hold on; plot(mdl); legend off
            drawnow, snapnow
            clear mdl
            disp('________________________________________________________________')

        elseif test_loss.Coefficients.pValue(3) < 0.05
            final_seed_to_seed_results.loss.Anhedonia.(publish_region_name_list_for_struct{r}) = test_loss;
            disp(strcat('During Loss ',analysis_string))
            disp('Symptom: Anhedonia')
            disp(strcat('ROI: ', publish_region_name_list_for_struct{r}))
            disp(strcat('t = ',num2str(test_loss.Coefficients.tStat(3)),'p = ',num2str(test_loss.Coefficients.pValue(3))))
            mdl = fitlm(Anhedonia,roi_avg_loss.(publish_region_name_list_for_struct{r}).dat);
            figure(); scatter(Anhedonia,roi_avg_loss.(publish_region_name_list_for_struct{r}).dat)
            title('Loss Anticipation')
            h1 = lsline();
            h1.LineWidth = 2;
            h1.Color = 'r';
            r1 = corrcoef(Anhedonia(:,1),roi_avg_loss.(publish_region_name_list_for_struct{r}).dat(:,1),'rows','complete');
            disp(r1(1,2));
            str = [' r = ',num2str(r1(1,2))]
            T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str);
            set(T, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');
            hold on; plot(mdl); legend off
            drawnow, snapnow
            clear mdl
            disp('________________________________________________________________')

        elseif test_loss.Coefficients.pValue(4) < 0.05
            final_seed_to_seed_results.loss.Fears.(publish_region_name_list_for_struct{r}) = test_loss;
            disp(strcat('During Loss ',analysis_string))
            disp('Symptom: Fears')
            disp(strcat('ROI: ', publish_region_name_list_for_struct{r}))
            disp(strcat('t = ',num2str(test_loss.Coefficients.tStat(4)),'p = ',num2str(test_loss.Coefficients.pValue(4))))
            mdl = fitlm(Fears,roi_avg_loss.(publish_region_name_list_for_struct{r}).dat);
            figure(); scatter(Fears,roi_avg_loss.(publish_region_name_list_for_struct{r}).dat)
            title('Loss Anticipation')
            h1 = lsline();
            h1.LineWidth = 2;
            h1.Color = 'r';
            r1 = corrcoef(Fears(:,1),roi_avg_loss.(publish_region_name_list_for_struct{r}).dat(:,1),'rows','complete');
            disp(r1(1,2));
            str = [' r = ',num2str(r1(1,2))]
            T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str);
            set(T, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');
            hold on; plot(mdl); legend off
            drawnow, snapnow
            clear mdl
            disp('________________________________________________________________')

        end
    end
end
clear r_gain r_loss
disp('STARTING WHOLE BRAIN ANALYSIS AND VISUALIZATION')
disp('_________________________________________________________________________')

results_struct.gain = threshold(temp_results_gain.t, primary_threshold,'unc','k',cluster_size);
results_struct.loss = threshold(temp_results_loss.t, primary_threshold,'unc','k',cluster_size);

for s = 1:length(symptom_names)

    % Current Seed Region
    disp(strcat('Gain ',analysis_string))
    % Current Symptom
    disp(symptom_names{s})

    disp('_________________________________________________________________________')

    r_gain.(symptom_names{s}) = region(select_one_image(results_struct.gain,s));
    info_gain = table(r_gain.(symptom_names{s}));

    if isempty(info_gain) == 1
        disp('I think the continue broke things before')
    else
        create_figure(strcat(strcat('Gain ',analysis_string), 'Symptom:',symptom_names{s})); axis off;
        canlab_results_fmridisplay(r_gain.(symptom_names{s}), 'regioncenters');
        drawnow, snapnow
        disp('_________________________________________________________________________')
        %keyboard
        %saveas(gcf,fullfile(figdir,filename))
    end


    % Current Seed Region
    disp(strcat('LOSS ',analysis_string))
    % Current Symptom
    disp(symptom_names{s})

    disp('_________________________________________________________________________')

    r_loss.(symptom_names{s}) = region(select_one_image(results_struct.loss,s));
    info_loss = table(r_loss.(symptom_names{s}));

    if isempty(info_loss) == 1
        disp('I think the continue broke things before')
    else

        create_figure(strcat('Loss , ', 'Symptom:',symptom_names{s})); axis off;
        canlab_results_fmridisplay(r_loss.(symptom_names{s}), 'regioncenters');
        drawnow, snapnow
        disp('_________________________________________________________________________')
        %keyboard
        %saveas(gcf,fullfile(figdir,filename))
    end
end
