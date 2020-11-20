%% Final Results Report

% region_list = {'bi_mPFC_Knutson','bi_VS_AntRew','bi_VS_AntLoss','bi_vmPFC_to_wholebrain','bi_Amyg_to_wholebrain','biOFC_to_wholebrain','L_VS_AntRew_to_wholebrain',...
% 'L_VS_AntLoss_to_wholebrain','LHO_Accumbens_to_wholebrain','LOFC2_to_wholebrain','LOFC_to_wholebrain','R_VS_AntLoss_to_wholebrain','R_VS_AntRew_to_wholebrain',...
% 'RHO_Accumbens_to_wholebrain','ROFC_to_wholebrain'};

% specifiy the regions of interest you actually want to do for wholebrain
% analysis
% oldham_list = {'biOFC_to_wholebrain','bi_VS_AntRew','bi_VS_AntLoss','ROFC_to_wholebrain','LOFC2_to_wholebrain','LOFC_to_wholebrain','L_VS_AntRew_to_wholebrain','L_VS_AntLoss_to_wholebrain','R_VS_AntLoss_to_wholebrain','R_VS_AntRew_to_wholebrain'};

%specific_roi_list = {'bi_Amyg_to_wholebrain','bi_VS_AntRew','bi_VS_AntLoss','L_VS_AntRew_to_wholebrain','L_VS_AntLoss_to_wholebrain','R_VS_AntLoss_to_wholebrain','R_VS_AntRew_to_wholebrain'};%{'bi_Amyg_to_wholebrain','bi_VS_AntRew','bi_VS_AntLoss','biOFC_to_wholebrain','LOFC2_to_wholebrain','LOFC_to_wholebrain','ROFC_to_wholebrain','L_VS_AntRew_to_wholebrain','L_VS_AntLoss_to_wholebrain','R_VS_AntLoss_to_wholebrain','R_VS_AntRew_to_wholebrain'};
primary_threshold = 0.05;
cluster_size = 1;
correction = 'fdr';

if strcmp(condir, 'anticipation')
    % anticipation list for first year project region_list = {'biOFC_to_wholebrain','bi_VS_AntRew','bi_VS_AntLoss','ROFC_to_wholebrain','LOFC2_to_wholebrain','LOFC_to_wholebrain','L_VS_AntRew_to_wholebrain','L_VS_AntLoss_to_wholebrain','R_VS_AntLoss_to_wholebrain','R_VS_AntRew_to_wholebrain'};
    specific_roi_list = {'biOFC_to_wholebrain','bi_VS_AntLoss','bi_VS_AntRew','ROFC_to_wholebrain'};%,'LOFC2_to_wholebrain','LOFC_to_wholebrain'};%,'R_VS_AntRew_to_wholebrain','L_VS_AntRew_to_wholebrain'};
else
    % consumption
    specific_roi_list = {'biOFC_to_wholebrain','bi_VS_Oldham','ROFC_to_wholebrain','LOFC2_to_wholebrain','LOFC_to_wholebrain'};%,'R_VS_to_wholebrain','L_VS_to_wholebrain'};
end
analysis = condir;
analyze = 1;
if analyze == 1
    % report region of interest analysis
    disp('STARTING SEED TO SEED ANALYSIS')
    disp('_________________________________________________________________________')
    for i = 1:length(specific_roi_list)
        for r = 1:length(region_name_list)
            % So basically we're looping through every seed region, and
            % then kind of every seed region again... But one refers to the
            % seed and the next refers to the region of interest extracted
            % for results. Then we'll check the p value associated with the
            % relationship and will pause to explore things further if
            % there's a significant relationship. I hope this works.
            test_gain = gain_mdl.(specific_roi_list{i}).(region_name_list_for_struct{r});
            if test_gain.Coefficients.pValue(2) < 0.05 || test_gain.Coefficients.pValue(3) < 0.05 || test_gain.Coefficients.pValue(4) < 0.05
                % This is going to create a results structure that will be
                % helpful for me to try and parse all the results I have.
                
                % Robin would like me to include graphs that show the
                % specific relationships between connectivity output. That
                % will probably work best if I pull out average estimates
                % for the regions here. Then I'll create specific plots in
                % the following if statement
                                          
                
                if test_gain.Coefficients.pValue(2) < 0.05
                    final_seed_to_seed_results.gain.(symptom_names{1}).(specific_roi_list{i}).(region_name_list_for_struct{r}) = test_gain;
                    disp(strcat('During Gain_',analysis))
                    disp(strcat('Symptom: ',symptom_names{1}))
                    disp(strcat('Seed region is: ',specific_roi_list{i}))
                    disp(strcat('End region is: ', region_name_list_for_struct{r}))
                    disp(strcat('t = ',num2str(test_gain.Coefficients.tStat(2)),'p = ',num2str(test_gain.Coefficients.pValue(2)))) 
                    
                    disp('________________________________________________________________')
                    
                elseif test_gain.Coefficients.pValue(3) < 0.05 & length(symptom_names) > 1
                    final_seed_to_seed_results.gain.Anhedonia.(specific_roi_list{i}).(region_name_list_for_struct{r}) = test_gain;
                    disp(strcat('During Gain_',analysis))
                    disp('Symptom: Anhedonia')
                    disp(strcat('Seed region is: ',specific_roi_list{i}))
                    disp(strcat('End region is: ', region_name_list_for_struct{r}))
                    disp(strcat('t = ',num2str(test_gain.Coefficients.tStat(3)),'p = ',num2str(test_gain.Coefficients.pValue(3))))
                    figure(); scatter(Anhedonia,roi_avg_gain.(specific_roi_list{i}).(region_name_list_for_struct{r}).dat)
                    drawnow, snapnow
                    disp('________________________________________________________________')
                elseif test_gain.Coefficients.pValue(4) < 0.05 & length(symptom_names) > 1
                    final_seed_to_seed_results.gain.Fears.(specific_roi_list{i}).(region_name_list_for_struct{r}) = test_gain;
                    disp(strcat('During Gain_',analysis))
                    disp('Symptom: Fears')
                    disp(strcat('Seed region is: ',specific_roi_list{i}))
                    disp(strcat('End region is: ', region_name_list_for_struct{r}))
                    disp(strcat('t = ',num2str(test_gain.Coefficients.tStat(4)),'p = ',num2str(test_gain.Coefficients.pValue(4))))
                    figure(); scatter(Fears,roi_avg_gain.(specific_roi_list{i}).(region_name_list_for_struct{r}).dat)
                    drawnow, snapnow
                    disp('________________________________________________________________')

                end
            end
            test_loss = loss_mdl.(specific_roi_list{i}).(region_name_list_for_struct{r});
            if test_loss.Coefficients.pValue(2) < 0.05 || test_loss.Coefficients.pValue(3) < 0.05 || test_loss.Coefficients.pValue(4) < 0.05
                
                if test_loss.Coefficients.pValue(2) < 0.05
                    final_seed_to_seed_results.loss.(symptom_names{1}).(specific_roi_list{i}).(region_name_list_for_struct{r}) = test_loss;
                    disp(strcat('During Loss_',analysis))
                    disp(strcat('Symptom: ',symptom_names{1}))
                    disp(strcat('Seed region is: ',specific_roi_list{i}))
                    disp(strcat('End region is: ', region_name_list_for_struct{r}))
                    disp(strcat('t = ',num2str(test_loss.Coefficients.tStat(2)),'p = ',num2str(test_loss.Coefficients.pValue(2))))
                    figure(); scatter(GenDis,roi_avg_loss.(specific_roi_list{i}).(region_name_list_for_struct{r}).dat)
                    drawnow, snapnow
                    disp('________________________________________________________________')
                  
                elseif test_loss.Coefficients.pValue(3) < 0.05 & length(symptom_names) > 1
                    final_seed_to_seed_results.loss.Anhedonia.(specific_roi_list{i}).(region_name_list_for_struct{r}) = test_loss;
                    disp(strcat('During Loss_',analysis))
                    disp('Symptom: Anhedonia')
                    disp(strcat('Seed region is: ',specific_roi_list{i}))
                    disp(strcat('End region is: ', region_name_list_for_struct{r}))
                    disp(strcat('t = ',num2str(test_loss.Coefficients.tStat(3)),'p = ',num2str(test_loss.Coefficients.pValue(3))))
                    figure(); scatter(Anhedonia,roi_avg_loss.(specific_roi_list{i}).(region_name_list_for_struct{r}).dat)
                    drawnow, snapnow
                    disp('________________________________________________________________')
                    
                elseif test_loss.Coefficients.pValue(4) < 0.05 & length(symptom_names) > 1
                    final_seed_to_seed_results.loss.Fears.(specific_roi_list{i}).(region_name_list_for_struct{r}) = test_loss;
                    disp(strcat('During Loss_',analysis))
                    disp('Symptom: Fears')
                    disp(strcat('Seed region is: ',specific_roi_list{i}))
                    disp(strcat('End region is: ', region_name_list_for_struct{r}))
                    disp(strcat('t = ',num2str(test_loss.Coefficients.tStat(4)),'p = ',num2str(test_loss.Coefficients.pValue(4))))
                    figure(); scatter(Fears,roi_avg_loss.(specific_roi_list{i}).(region_name_list_for_struct{r}).dat)
                    drawnow, snapnow
                    disp('________________________________________________________________')

                end
            end
        end        
    end
    % report region of interest analysis
    clear test_gain test_loss
    disp('STARTING SEED TO SEED ANALYSIS: CONTROL FOR OPPOSITE CONDITION')
    disp('_________________________________________________________________________')
    for i = 1:length(specific_roi_list)
        for r = 1:length(region_name_list)
            % So basically we're looping through every seed region, and
            % then kind of every seed region again... But one refers to the
            % seed and the next refers to the region of interest extracted
            % for results. Then we'll check the p value associated with the
            % relationship and will pause to explore things further if
            % there's a significant relationship. I hope this works.
            test_gain = gain_mdl2.(specific_roi_list{i}).(region_name_list_for_struct{r});
            if test_gain.Coefficients.pValue(2) < 0.05 || test_gain.Coefficients.pValue(3) < 0.05 || test_gain.Coefficients.pValue(4) < 0.05
                % This is going to create a results structure that will be
                % helpful for me to try and parse all the results I have.
                % Controlling for sex really opened the box
                if test_gain.Coefficients.pValue(2) < 0.05
                    final_seed_to_seed_results.gain.(symptom_names{1}).(specific_roi_list{i}).(region_name_list_for_struct{r}) = test_gain;
                    disp(strcat('During Gain_',analysis))
                    disp(strcat('Symptom: ',symptom_names{1}))
                    disp(strcat('Seed region is: ',specific_roi_list{i}))
                    disp(strcat('End region is: ', region_name_list_for_struct{r}))
                    disp(strcat('t = ',num2str(test_gain.Coefficients.tStat(2)),'p = ',num2str(test_gain.Coefficients.pValue(2))))
                    figure(); scatter(GenDis,roi_avg_gain2.(specific_roi_list{i}).(region_name_list_for_struct{r}).dat)
                    drawnow, snapnow
                    disp('________________________________________________________________')
                elseif test_gain.Coefficients.pValue(3) < 0.05 & length(symptom_names) > 1
                    final_seed_to_seed_results.gain.Anhedonia.(specific_roi_list{i}).(region_name_list_for_struct{r}) = test_gain;
                    disp(strcat('During Gain_',analysis))
                    disp('Symptom: Anhedonia')
                    disp(strcat('Seed region is: ',specific_roi_list{i}))
                    disp(strcat('End region is: ', region_name_list_for_struct{r}))
                    disp(strcat('t = ',num2str(test_gain.Coefficients.tStat(3)),'p = ',num2str(test_gain.Coefficients.pValue(3))))
                    figure(); scatter(Anhedonia,roi_avg_gain2.(specific_roi_list{i}).(region_name_list_for_struct{r}).dat)
                    drawnow, snapnow
                    disp('________________________________________________________________')

                elseif test_gain.Coefficients.pValue(4) < 0.05 & length(symptom_names) > 1
                    final_seed_to_seed_results.gain.Fears.(specific_roi_list{i}).(region_name_list_for_struct{r}) = test_gain;
                    disp(strcat('During Gain_',analysis))
                    disp('Symptom: Fears')
                    disp(strcat('Seed region is: ',specific_roi_list{i}))
                    disp(strcat('End region is: ', region_name_list_for_struct{r}))
                    disp(strcat('t = ',num2str(test_gain.Coefficients.tStat(4)),'p = ',num2str(test_gain.Coefficients.pValue(4))))
                    figure(); scatter(Fears,roi_avg_gain2.(specific_roi_list{i}).(region_name_list_for_struct{r}).dat)
                    drawnow, snapnow
                    disp('________________________________________________________________')

                end
            end
            test_loss = loss_mdl2.(specific_roi_list{i}).(region_name_list_for_struct{r});
            if test_loss.Coefficients.pValue(2) < 0.05 || test_loss.Coefficients.pValue(3) < 0.05 || test_loss.Coefficients.pValue(4) < 0.05
                if test_loss.Coefficients.pValue(2) < 0.05
                    final_seed_to_seed_results.loss.(symptom_names{1}).(specific_roi_list{i}).(region_name_list_for_struct{r}) = test_loss;
                    disp(strcat('During Loss_',analysis))
                    disp(strcat('Symptom: ',symptom_names{1}))
                    disp(strcat('Seed region is: ',specific_roi_list{i}))
                    disp(strcat('End region is: ', region_name_list_for_struct{r}))
                    disp(strcat('t = ',num2str(test_loss.Coefficients.tStat(2)),'p = ',num2str(test_loss.Coefficients.pValue(2))))
                    figure(); scatter(GenDis,roi_avg_loss2.(specific_roi_list{i}).(region_name_list_for_struct{r}).dat)
                    drawnow, snapnow
                    disp('________________________________________________________________')
                    
                elseif test_loss.Coefficients.pValue(3) < 0.05 & length(symptom_names) > 1
                    final_seed_to_seed_results.loss.Anhedonia.(specific_roi_list{i}).(region_name_list_for_struct{r}) = test_loss;
                    disp(strcat('During Loss_',analysis))
                    disp('Symptom: Anhedonia')
                    disp(strcat('Seed region is: ',specific_roi_list{i}))
                    disp(strcat('End region is: ', region_name_list_for_struct{r}))
                    disp(strcat('t = ',num2str(test_loss.Coefficients.tStat(3)),'p = ',num2str(test_loss.Coefficients.pValue(3))))
                    figure(); scatter(Anhedonia,roi_avg_loss2.(specific_roi_list{i}).(region_name_list_for_struct{r}).dat)
                    drawnow, snapnow
                    disp('________________________________________________________________')
                    
                elseif test_loss.Coefficients.pValue(4) < 0.05 & length(symptom_names) > 1
                    final_seed_to_seed_results.loss.Fears.(specific_roi_list{i}).(region_name_list_for_struct{r}) = test_loss;
                    disp(strcat('During Loss_',analysis))
                    disp('Symptom: Fears')
                    disp(strcat('Seed region is: ',specific_roi_list{i}))
                    disp(strcat('End region is: ', region_name_list_for_struct{r}))
                    disp(strcat('t = ',num2str(test_loss.Coefficients.tStat(4)),'p = ',num2str(test_loss.Coefficients.pValue(4))))
                    figure(); scatter(Fears,roi_avg_loss2.(specific_roi_list{i}).(region_name_list_for_struct{r}).dat)
                    drawnow, snapnow
                    disp('________________________________________________________________')

                end
            end
        end        
    end
    
    clear r_gain r_loss
    disp('STARTING WHOLE BRAIN ANALYSIS AND VISUALIZATION')
    disp('_________________________________________________________________________')
    for i = 1:length(specific_roi_list)
        results_struct.gain.(specific_roi_list{i}) = threshold(temp_results_gain.(specific_roi_list{i}).t,primary_threshold,correction,'k',cluster_size);
        results_struct.loss.(specific_roi_list{i}) = threshold(temp_results_loss.(specific_roi_list{i}).t,primary_threshold,correction,'k',cluster_size);
        
        for s = 1:length(symptom_names)
            
            % Current Seed Region
            disp('GAIN ANTICIPTATION')
            disp(specific_roi_list{i})
            % Current Symptom
            disp(symptom_names{s})
            
            disp('_________________________________________________________________________')
            
            r_gain.(specific_roi_list{i}).(symptom_names{s}) = region(select_one_image(results_struct.gain.(specific_roi_list{i}),s));
            info_gain = table(r_gain.(specific_roi_list{i}).(symptom_names{s}));
            
            if isempty(info_gain) == 1
                disp('I think the continue broke things before')
            else
                create_figure(strcat('Gain anticipation montage Seed:',specific_roi_list{i},' Symptom:',symptom_names{s})); axis off;
                canlab_results_fmridisplay(r_gain.(specific_roi_list{i}).(symptom_names{s}), 'regioncenters');
                drawnow, snapnow
                filename = strcat('Gain_anticipation_montage_Seed_',specific_roi_list{i},'_Symptom_',symptom_names{s},'.pdf');
                disp('_________________________________________________________________________')
                %keyboard
                %saveas(gcf,fullfile(figdir,filename))
            end
            
            
            % Current Seed Region
            disp('LOSS ANTICIPATION')
            disp(specific_roi_list{i})
            % Current Symptom
            disp(symptom_names{s})
            
            disp('_________________________________________________________________________')
            
            r_loss.(specific_roi_list{i}).(symptom_names{s}) = region(select_one_image(results_struct.loss.(specific_roi_list{i}),s));
            info_loss = table(r_loss.(specific_roi_list{i}).(symptom_names{s}));
            
            if isempty(info_loss) == 1
                disp('I think the continue broke things before')
            else
                
                create_figure(strcat('Loss anticipation montage Seed:',specific_roi_list{i},' Symptom:',symptom_names{s})); axis off;
                canlab_results_fmridisplay(r_loss.(specific_roi_list{i}).(symptom_names{s}), 'regioncenters');
                drawnow, snapnow
                filename = strcat('Loss_anticipation_montage_Seed_',specific_roi_list{i},'_Symptom_',symptom_names{s},'.pdf');
                disp('_________________________________________________________________________')
                %keyboard
                %saveas(gcf,fullfile(figdir,filename))
            end
        end
    end
end
