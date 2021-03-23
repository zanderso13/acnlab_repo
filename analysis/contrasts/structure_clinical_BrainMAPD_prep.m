% This used to be a script to compare insula structure and clinical
% measures... Eventually it turned into a workbench that I repurposed to
% the point of insanity. SO! I've deleted a bunch of stuff. (I never really
% wanted to study GMV anyway...

% set up you subjects directory
basedir = '/projects/p30954/struct_data/session1';
clinicaldir = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD/clinical_data';
repodir = '/home/zaz3744/repo';
demodir = '/projects/p30954/demographic_data';

%% load demographics
demo_file = filenames(fullfile(demodir,'BrainMAPDT1S1Demo.mat'));
load(demo_file{1});
    
%% Load in some Tri level factor scores and plot them. That's the goal here
clinicalfname = filenames(fullfile(clinicaldir, 'Multi*FS.mat'));
load(clinicalfname{1});

% Plot distribution of z-scored factor scores
% They're normally distributed... so ya
figure(); 
subplot(1,3,1)
violinplot(trilevelmultimethodFS.GenDis); legend off
title('Subplot 1:General Distress')
subplot(1,3,2)
violinplot(trilevelmultimethodFS.Anhedon); legend off
title('Subplot 2: Anhedonia')
subplot(1,3,3)
violinplot(trilevelmultimethodFS.Fears); legend off
title('Subplot 3: Fears')

%% Load in the SCID and output some frequency charts looking at hom many of each disorder we have. 
%scidfname = filenames(fullfile(clinicaldir, 'SCID_data.mat'));
%load(scidfname{1});

% Want to get counts for relevant DSM disorders
% Past and curr(:)ent depression
% Past and curr(:)ent anxiety
% Life(:)time comorbidity of depression and anxiety
% Other category (for now, likely want to flesh this out more but may just
% exclude anything that falls in here)

% Depression variables list: SCID_data.T1MDElife(:), SCID_data.T1MDEOSCurr(:), SCID_data.T1MDEOSlife(:),
% SCID_data.T1PDDlife(:), SCID_data.T1PDDOSlife(:), SCID_data.T1UPDDxlife(:), SCID_data.T1UPDOSlife(:), SCID_data.T1CYCDxlife(:),
% SCID_data.T1CYCOSlife(:)
% Mania/BD/psychosis variables list: SCID_data.T1PSYDxlife(:), SCID_data.T1PSYlife(:), SCID_data.T1BDIIDxlife(:), SCID_data.T1BDIIOSlife(:),SCID_data.T1MANDxlife(:), SCID_data.T1MANOSlife(:),
% SCID_data.T1HMANDxlife(:), SCID_data.T1HMANOSlife(:), SCID_data.T1BSDDxlife(:), SCID_data.T1BSDOSlife(:)
% Anxiety related variables list: SCID_data.T1PDDxLife(:), SCID_data.T1PDOSLife(:), SCID_data.T1AGODxLife(:),
% SCID_data.T1AGOOSLife(:), SCID_data.T1SADDxLife(:), SCID_data.T1SADOSLife(:), SCID_data.T1SPDxLife(:), SCID_data.T1SPOSLife(:),
% SCID_data.T1GADDxLife(:), SCID_data.T1GADOSLife(:), SCID_data.T1ADDxlife(:), SCID_data.T1ADOSlife(:),
% SCID_data.T1OCDDxlife(:), SCID_data.T1OCDOSlife(:), SCID_data.T1SEPDxlife(:), SCID_data.T1SEPOSlife(:), SCID_data.T1HOADxlife(:),
% SCID_data.T1HOAOSlife(:), SCID_data.T1SSDDxlife(:), SCID_data.T1SSDOSlife(:), SCID_data.T1IADDxlife(:), SCID_data.T1IADOSlife(:),
% SCID_data.T1GAMDxlife(:), SCID_data.T1GAMOSlife(:), SCID_data.T1AUDDxlife(:), SCID_data.T1AUDOSlife(:), SCID_data.T1SUDDxlife(:),
% SCID_data.T1SUDOSlife(:) 
% Eating disorders and psychosis: SCID_data.T1ANDxlife(:), SCID_data.T1ANOSlife(:), SCID_data.T1BNDxlife(:),
% SCID_data.T1BNOSlife(:)
% Trauma based disorders: SCID_data.T1ASDDxlife(:), SCID_data.T1ASDOSlife(:), SCID_data.T1PTSDDxlife(:),
% SCID_data.T1PTSDOSlife(:), SCID_data.T1ADJDxlife(:), SCID_data.T1ADJOSlife(:), SCID_data.T1TRADxlife(:), SCID_data.T1TRAOSlife(:)

% Unclear what to do with this sum(SCID_data.T1MDEOSCurr(:)),  and other
% variables combining OS and Dx. I'm thinking nothing for now.

depression_all = [sum(SCID_data.T1MDElife(:)), sum(SCID_data.T1MDEOSlife(:)), sum(SCID_data.T1PDDlife(:)), sum(SCID_data.T1PDDOSlife(:)), sum(SCID_data.T1UPDDxlife(:)), sum(SCID_data.T1UPDOSlife(:)), sum(SCID_data.T1CYCDxlife(:)), sum(SCID_data.T1CYCOSlife(:))];
anxiety_all = [sum(SCID_data.T1PDDxLife(:)), sum(SCID_data.T1PDOSLife(:)), sum(SCID_data.T1AGODxLife(:)), sum(SCID_data.T1AGOOSLife(:)), sum(SCID_data.T1SADDxLife(:)), sum(SCID_data.T1SADOSLife(:)), sum(SCID_data.T1SPDxLife(:)), sum(SCID_data.T1SPOSLife(:)), sum(SCID_data.T1GADDxLife(:)), sum(SCID_data.T1GADOSLife(:)), sum(SCID_data.T1ADDxlife(:)), sum(SCID_data.T1ADOSlife(:)), sum(SCID_data.T1OCDDxlife(:)), sum(SCID_data.T1OCDOSlife(:)), sum(SCID_data.T1SEPDxlife(:)), sum(SCID_data.T1SEPOSlife(:)), sum(SCID_data.T1HOADxlife(:)), sum(SCID_data.T1HOAOSlife(:)), sum(SCID_data.T1SSDDxlife(:)), sum(SCID_data.T1SSDOSlife(:)), sum(SCID_data.T1IADDxlife(:)), sum(SCID_data.T1IADOSlife(:))]; 
substance_gambling_all = [sum(SCID_data.T1GAMDxlife(:)), sum(SCID_data.T1GAMOSlife(:)), sum(SCID_data.T1AUDDxlife(:)), sum(SCID_data.T1AUDOSlife(:)), sum(SCID_data.T1SUDDxlife(:)), sum(SCID_data.T1SUDOSlife(:))]; 
trauma_all = [sum(SCID_data.T1ASDDxlife(:)), sum(SCID_data.T1ASDOSlife(:)), sum(SCID_data.T1PTSDDxlife(:)), sum(SCID_data.T1PTSDOSlife(:)), sum(SCID_data.T1ADJDxlife(:)), sum(SCID_data.T1ADJOSlife(:)), sum(SCID_data.T1TRADxlife(:)), sum(SCID_data.T1TRAOSlife(:))];
psychosis_mania_all = [sum(SCID_data.T1PSYDxlife(:)), sum(SCID_data.T1PSYlife(:)), sum(SCID_data.T1BDIIDxlife(:)), sum(SCID_data.T1BDIIOSlife(:)), sum(SCID_data.T1MANDxlife(:)), sum(SCID_data.T1MANOSlife(:)), sum(SCID_data.T1HMANDxlife(:)), sum(SCID_data.T1HMANOSlife(:)), sum(SCID_data.T1BSDDxlife(:)), sum(SCID_data.T1BSDOSlife(:))];
eating_all = [sum(SCID_data.T1ANDxlife(:)), sum(SCID_data.T1ANOSlife(:)), sum(SCID_data.T1BNDxlife(:)), sum(SCID_data.T1BNOSlife(:))];

depression_dx = depression_all(1, 1:2:end);
anxiety_dx = anxiety_all(1, 1:2:end);
substance_gambling_dx = substance_gambling_all(1, 1:2:end);
trauma_dx = trauma_all(1, 1:2:end);
psychosis_mania_dx = psychosis_mania_all(1, 1:2:end);
eating_dx = eating_all(1, 1:2:end);

depression_os = depression_all(1, 2:2:end);
anxiety_os = anxiety_all(1, 2:2:end);
substance_gambling_os = substance_gambling_all(1, 2:2:end);
trauma_os = trauma_all(1, 2:2:end);
psychosis_mania_os = psychosis_mania_all(1, 2:2:end);
eating_os = eating_all(1, 2:2:end);

depression_label = categorical({'MDE', 'PDD', 'UPD', 'CYC'});
anxiety_label = categorical({'Panic', 'AGO', 'SAD', 'Phobia', 'GAD', 'Any Anxiety', 'OCD', 'Separation', 'Hoarding', 'Somatic', 'Illness'});
gamb_label = categorical({'Gambling', 'Alcohol', 'Substance'});
trauma_label = categorical({'Acute Stress', 'PTSD', 'Adjustment', 'Any lifetime trauma'});
psychosis_label = categorical({'Psychosis', 'BDII', 'Mania', 'Hypomania', 'Borderline'});
eating_label = categorical({'Anorexia', 'Bulimia'});

% dx plots for DSM disorders
figure();
subplot(1,6,1)
title('Depression Dx')
bar(depression_label, depression_dx)
ylim([0 100])
    
subplot(1,6,2)
title('Anxiety Dx')
bar(anxiety_label, anxiety_dx)
ylim([0 100])

subplot(1,6,3)
title('Gambling and Substance Use Dx')
bar(gamb_label, substance_gambling_dx)
ylim([0 100])

subplot(1,6,4)
title('Trauma Dx')
bar(trauma_label, trauma_dx)
ylim([0 100])

subplot(1,6,5)
title('Eating Disorder Dx')
bar(eating_label, eating_dx)
ylim([0 100])

subplot(1,6,6)
title('Psychosis and Mania Dx')
bar(psychosis_label, psychosis_mania_dx)
ylim([0 100])

% os plots for DSM disorders
figure();
subplot(1,6,1)
title('Depression OS')
bar(depression_label, depression_os)
ylim([0 100])
    
subplot(1,6,2)
title('Anxiety OS')
bar(anxiety_label, anxiety_os)
ylim([0 100])

subplot(1,6,3)
title('Gambling and Substance Use OS')
bar(gamb_label, substance_gambling_os)
ylim([0 100])

subplot(1,6,4)
title('Trauma OS')
bar(trauma_label, trauma_os)
ylim([0 100])

subplot(1,6,5)
title('Eating Disorder OS')
bar(eating_label, eating_os)
ylim([0 100])

subplot(1,6,6)
title('Psychosis and Mania OS')
bar(psychosis_label, psychosis_mania_os)
ylim([0 100])

%% find comorbidity between anxiety and depression. Also create freq diagrams for any depression, any anxiety, all comorbidity
% first concatenate all depression and anxiety related disorders into
% single binary table.
dep_all = [SCID_data.T1ANYPDDlife(:), SCID_data.T1ANYUPDlife(:), SCID_data.T1ANYCYClife(:)];
anx_all = [SCID_data.T1ANYPDLife(:), SCID_data.T1ANYAGOLife(:), SCID_data.T1ANYSADLife(:), SCID_data.T1ANYSPLife(:), SCID_data.T1ANYGADLife(:), SCID_data.T1ANYADlife(:), SCID_data.T1ANYOCDlife(:), SCID_data.T1ANYSEPlife(:), SCID_data.T1ANYHOAlife(:), SCID_data.T1ANYIADlife(:), SCID_data.T1ANYPTSDlife(:)];

% Now reduce that table to single columns. Each participant only needs to
% have one of the above diagnoses to be a 1 in this list
dep_general(:,1) = SCID_data.PID;
dep_general(:,2) = max(dep_all, [], 2);
anx_general(:,1) = SCID_data.PID;
anx_general(:,2) = max(anx_all, [], 2);

% Now that I know all the subs with some kind of depression or some kind of
% anxiety. I need to find the subs that have both
comorbidity(:,1) = SCID_data.PID;
comorbidity(:,2) = dep_general(:,2) == 1 & anx_general(:,2) == 1;

%% Gonna force some counts of current depression vs anxiety in here
mde_curr(:,1) = SCID_data.PID;
mde_curr(:,2) = [SCID_data.T1MDEcurr(:)];
pdd_curr(:,1) = SCID_data.PID;
pdd_curr(:,2) = [SCID_data.T1PDDCurr(:)];
% anx_curr = [SCID_data.T1PDcurr(:), SCID_data.T1AGOcurr(:), SCID_data.T1SADcurr(:), SCID_data.T1SPcurr(:), SCID_data.T1GADcurr(:), SCID_data.T1ADcurr(:), SCID_data.T1OCDcurr(:), SCID_data.T1SEPcurr(:), SCID_data.T1HOAcurr(:), SCID_data.T1IADcurr(:), SCID_data.T1PTSDcurr(:)];
anx_curr(:,1) = SCID_data.PID;
anx_curr(:,2) = [SCID_data.T1ADcurr(:)];

%%
% I have too many people and the likely cause of this is my sample is 330
% instead of 265... so let's fix that by looking for a corresponding
% structural file and then deleting the row if a struct file doesn't exist
for sub = 1:length(SCID_data.PID)

    if exist(fullfile(basedir, strcat(num2str(SCID_data.PID(sub)) ,'_MR1'))) == 7
        continue
    else
        
        missing_sub = find(dep_general(:,1) == SCID_data.PID(sub));
        dep_general(missing_sub,:) = [];
        anx_general(missing_sub,:) = [];
        comorbidity(missing_sub,:) = [];
        mde_curr(missing_sub,:) = [];
        pdd_curr(missing_sub,:) = [];
        anx_curr(missing_sub,:) = [];
    end
end

disorder_counts = [sum(dep_general(:,2)), sum(anx_general(:,2)), sum(comorbidity(:,2))];
counts_label = categorical({'Depression', 'Anxiety', 'Comorbid'});

figure()
title('Overall Disorder Count')
bar(counts_label, disorder_counts)
ylim([0 110])

curr_freq = [sum(mde_curr(:,2)), sum(pdd_curr(:,2)), sum(anx_curr(:,2))];
curr_counts = categorical({'Current MDE', 'Current PDD', 'Current Anxiety'});

figure();
title('Freq of current anxiety and depression')
bar(curr_counts, curr_freq)


%% time to get counts for the screener bins

load(fullfile(clinicaldir,'BrainMAPD_screener_bins.mat'));
% I have all these saved data files but I probably need to start from
% scratch to build the comprehensive mat file. Yoni would probably be able
% to do this no problem but the table syntax is giving me a headache so I'm
% just going to make a double and then convert it to a table before saving
% the comprehensive mat file. I only hope, one day, the yonestar will
% forgive me for my indiscretions

%final_table_bins(:,1) = clinical_bins.PID(:);
%final_table_bins(:,2) = clinical_bins.Bin(:);

% this is janky but the lh_tabl and rh_tabl are pissing me off. So now
% they're doubles. I'll need to copy the properties into the final table I
% create. But thankfully I'll be done with all this after I finish this
% last portion. Then I can start on a fresh script and maybe get to the
% paper June just emailed me about.... Sigh... Update on Mar 25! You
% finally got to the paper! And resubmitted! Yay!
clear sub
for sub = 1:length(lh_tabl.subject_number)
    temp_lh_struct_data(sub,1) = lh_tabl.subject_number{sub};
    temp_lh_struct_data(sub,2:width(lh_tabl)) = lh_tabl{sub,2:width(lh_tabl)};
end
clear sub
for sub = 1:length(rh_tabl.subject_number)
    temp_rh_struct_data(sub,1) = rh_tabl.subject_number{sub};
    temp_rh_struct_data(sub,2:width(rh_tabl)) = rh_tabl{sub,2:width(rh_tabl)};
end

final_table_bins(:,1) = temp_lh_struct_data(:,1);

clear sub
for sub = 1:length(final_table_bins)
    % the clinical bins table is the largest which is why it is the
    % beginning of this exercise. Now I need to add elements. This will
    % include the dep_general, anx_general, and comorbidity variables so I
    % can create DSM counts for each clinical bin. Then I'm putting the
    % struct data in there as well. Once I have this I'll convert it back
    % to a table. 
    % insert structural data
    if isempty(find(final_table_bins(sub,1) == temp_lh_struct_data(:,1))) == 1
        final_table_bins(sub,14:width(lh_tabl)+13) = NaN;
        final_table_bins(sub,2:8) = 0;
        continue
    else
        curr_index_struct_data = find(final_table_bins(sub,1) == temp_lh_struct_data(:,1));
        final_table_bins(sub,2:8) = 0;
        final_table_bins(sub,14:width(lh_tabl)+12) = table2array(lh_tabl(curr_index_struct_data,2:width(lh_tabl)));
    end
    
    % insert dsm data
    if isempty(find(final_table_bins(sub,1) == dep_general(:,1))) 
        final_table_bins(sub,3:5) = NaN;
        continue
    else
        curr_index_dsm_data = find(final_table_bins(sub,1) == dep_general(:,1));
        final_table_bins(sub,3:5) = [dep_general(curr_index_dsm_data,2), anx_general(curr_index_dsm_data,2), comorbidity(curr_index_dsm_data,2)];
    end
    
    % insert gender data
    if isempty(find(final_table_bins(sub,1) == BrainMAPDT1S1Demo.PID(:)))
        final_table_bins(sub,9) = NaN;
        continue
    else
        curr_index_demo_data = find(final_table_bins(sub,1) == BrainMAPDT1S1Demo.PID(:));
        final_table_bins(sub,9) = BrainMAPDT1S1Demo.sex(curr_index_demo_data);
    end
    
    % insert current disorder count
    if isempty(find(final_table_bins(sub,1) == anx_curr(:,1)))
        final_table_bins(sub,6:8) = NaN;
        continue
    else
        curr_index_curr_disorder = find(final_table_bins(sub,1) == anx_curr(:,1));
        final_table_bins(sub,6:8) = [mde_curr(curr_index_curr_disorder,2), pdd_curr(curr_index_curr_disorder,2), anx_curr(curr_index_curr_disorder,2)];
    end
    
    % insert bin data
    if isempty(find(final_table_bins(sub,1) == clinical_bins.PID(:)))
        final_table_bins(sub,2) = NaN;
        continue
    else
        curr_index_bin_data = find(final_table_bins(sub,1) == clinical_bins.PID(:));
        final_table_bins(sub,2) = clinical_bins.Bin(curr_index_bin_data);
    end
    
    % insert trilevel data
    if isempty(find(final_table_bins(sub,1) == trilevelmultimethodFS.id(:))) 
        final_table_bins(sub,10:12) = NaN;
        continue
    else
        curr_index_tri_data = find(final_table_bins(sub,1) == trilevelmultimethodFS.id(:));
        final_table_bins(sub,10:13) = [trilevelmultimethodFS.GenDis(curr_index_tri_data), trilevelmultimethodFS.Anhedon(curr_index_tri_data), trilevelmultimethodFS.Fears(curr_index_tri_data), trilevelmultimethodFS.Arousal(curr_index_tri_data)];
    end
    
    
    
end
    
% Doing another janky edit. The above code (although I don't like how I did it) gets the job done.
% I now have a single table that incorporates structural, dsm, and screener
% data. I'll need to check how much data I lost. I don't know how well this
% all lines up. Sigh, I hate this part of science. But ok, so now let's
% change the comorbid variable so that I'm not double counting people

clear sub
for sub = 1:length(final_table_bins)
    if final_table_bins(sub,5) == 1
        final_table_bins(sub,4) = 0;
        final_table_bins(sub,3) = 0;
    else
        continue
    end
end

%% Final set of plots
% Frequency diagrams for each of the 9 bins. Identify the bins first, then
% get a count of depression, anxiety, comorbid in each one. 

high_high_index = find(final_table_bins(:,2) == 1);
high_low_index = find(final_table_bins(:,2) == 2);
med_med_index = find(final_table_bins(:,2) == 3);
low_high_index = find(final_table_bins(:,2) == 4);
low_low_index = find(final_table_bins(:,2) == 5);
high_med_index = find(final_table_bins(:,2) == 6);
low_med_index = find(final_table_bins(:,2) == 7);
med_high_index = find(final_table_bins(:,2) == 8);
med_low_index = find(final_table_bins(:,2) == 9);

% depression frequencies for each bin
high_high_dep = nansum(final_table_bins(high_high_index, 3));
high_low_dep = nansum(final_table_bins(high_low_index, 3));
med_med_dep = nansum(final_table_bins(med_med_index, 3));
low_high_dep = nansum(final_table_bins(low_high_index, 3));
low_low_dep = nansum(final_table_bins(low_low_index, 3));
high_med_dep = nansum(final_table_bins(high_med_index, 3));
low_med_dep = nansum(final_table_bins(low_med_index, 3));
med_high_dep = nansum(final_table_bins(med_high_index, 3));
med_low_dep = nansum(final_table_bins(med_low_index, 3));
dep_hist_input = [high_high_dep, high_low_dep, med_med_dep, low_high_dep, low_low_dep, high_med_dep, low_med_dep, med_high_dep, med_low_dep]; 

% anxiety frequencies for each bin
high_high_anx = nansum(final_table_bins(high_high_index, 4));
high_low_anx = nansum(final_table_bins(high_low_index, 4));
med_med_anx = nansum(final_table_bins(med_med_index, 4));
low_high_anx = nansum(final_table_bins(low_high_index, 4));
low_low_anx = nansum(final_table_bins(low_low_index, 4));
high_med_anx = nansum(final_table_bins(high_med_index, 4));
low_med_anx = nansum(final_table_bins(low_med_index, 4));
med_high_anx = nansum(final_table_bins(med_high_index, 4));
med_low_anx = nansum(final_table_bins(med_low_index, 4));
anx_hist_input = [high_high_anx, high_low_anx, med_med_anx, low_high_anx, low_low_anx, high_med_anx, low_med_anx, med_high_anx, med_low_anx]; 


% comorbid frequencies for each bin
high_high_com = nansum(final_table_bins(high_high_index, 5));
high_low_com = nansum(final_table_bins(high_low_index, 5));
med_med_com = nansum(final_table_bins(med_med_index, 5));
low_high_com = nansum(final_table_bins(low_high_index, 5));
low_low_com = nansum(final_table_bins(low_low_index, 5));
high_med_com = nansum(final_table_bins(high_med_index, 5));
low_med_com = nansum(final_table_bins(low_med_index, 5));
med_high_com = nansum(final_table_bins(med_high_index, 5));
med_low_com = nansum(final_table_bins(med_low_index, 5));
com_hist_input = [high_high_com, high_low_com, med_med_com, low_high_com, low_low_com, high_med_com, low_med_com, med_high_com, med_low_com]; 
    
%% Time to make the real final table
% Need to put final_table variable into a format that everyone can
% access/add to if they find new subjects I missed, all that good stuff.

BrainMAPD_subject_struct_table = array2table(final_table_bins);
BrainMAPD_subject_struct_table = [BrainMAPD_subject_struct_table,array2table(temp_rh_struct_data(:,1:75))];

BrainMAPD_subject_struct_table.Properties.VariableNames = {'PID', 'Clinical_bins', 'Dep_lifetime', 'Anx_lifetime', 'Comorbid_lifetime', 'MDE_current', 'PDD_current', 'Anx_current', 'Sex', 'GenDis', 'Anhedon', 'Fears','Arousal', lh_tabl.Properties.VariableNames{2:75}, 'BrainSegVolNotVent', 'eTIV' 'PID2', rh_tabl.Properties.VariableNames{2:75}};
BrainMAPD_subject_struct_table(find(BrainMAPD_subject_struct_table.Sex == 2), :) = [];
BrainMAPD_subject_diagnosis_table = BrainMAPD_subject_struct_table(:,1:8);

%% Set up regression variables
rh_ant_insula = sum([BrainMAPD_subject_struct_table.rh_G_insular_short_volume, BrainMAPD_subject_struct_table.rh_S_circular_insula_ant_volume],2);
lh_ant_insula = sum([BrainMAPD_subject_struct_table.lh_G_insular_short_volume, BrainMAPD_subject_struct_table.lh_S_circular_insula_ant_volume],2);

rh_post_insula = sum([BrainMAPD_subject_struct_table.rh_G_Ins_lg_S_cent_ins_volume, BrainMAPD_subject_struct_table.rh_S_circular_insula_inf_volume, BrainMAPD_subject_struct_table.rh_S_circular_insula_sup_volume],2);
lh_post_insula = sum([BrainMAPD_subject_struct_table.lh_G_Ins_lg_S_cent_ins_volume, BrainMAPD_subject_struct_table.lh_S_circular_insula_inf_volume, BrainMAPD_subject_struct_table.lh_S_circular_insula_sup_volume],2);

regressors = [BrainMAPD_subject_struct_table.Dep_lifetime, BrainMAPD_subject_struct_table.Anx_lifetime, BrainMAPD_subject_struct_table.Comorbid_lifetime, BrainMAPD_subject_struct_table.Sex, BrainMAPD_subject_struct_table.eTIV];
% ones(length(rh_ant_insula),1),

%% Set up anovan variables
% diagnostic criteria
for sub = 1:size(BrainMAPD_subject_struct_table,1)
    if BrainMAPD_subject_struct_table.Dep_lifetime(sub) == 1
        anova_input_dsm(sub) = 1;
    elseif BrainMAPD_subject_struct_table.Anx_lifetime(sub) == 1
        anova_input_dsm(sub) = 2;
    elseif BrainMAPD_subject_struct_table.Comorbid_lifetime(sub) == 1
        anova_input_dsm(sub) = 3;
    else
        anova_input_dsm(sub) = 4;
    end
end
anova_input_dsm = anova_input_dsm';
anova_input_sex = BrainMAPD_subject_struct_table.Sex; 
anova_input_headsize = BrainMAPD_subject_struct_table.eTIV;

%% Model with GLM
mdl1_lh_post = fitlm(regressors, lh_post_insula, 'Categorical', logical([1 1 1 1 0]))
mdl2_rh_post = fitlm(regressors, rh_post_insula, 'Categorical', logical([1 1 1 1 0]))
mdl3_lh_ant = fitlm(regressors, lh_ant_insula, 'Categorical', logical([1 1 1 1 0]))
mdl4_rh_ant = fitlm(regressors, rh_ant_insula, 'Categorical', logical([1 1 1 1 0]))

%% Hierarchical model 
% Need to model my nuisance regressors first
% Then I'll include diagnoses in the model
% I'll finish up by manually subtracting R^2 values

level_one_regressors = [BrainMAPD_subject_struct_table.Sex, BrainMAPD_subject_struct_table.eTIV];

level_one_model_rh_ant = fitlm(level_one_regressors, rh_ant_insula, 'Categorical', logical([1,0]))
level_one_model_lh_ant = fitlm(level_one_regressors, lh_ant_insula, 'Categorical', logical([1,0]))
level_one_model_rh_post = fitlm(level_one_regressors, rh_post_insula, 'Categorical', logical([1,0]))
level_one_model_lh_post = fitlm(level_one_regressors, lh_post_insula, 'Categorical', logical([1,0]))

% Subtract R^2 values to get the variance accounted for by internalizing
% disorders generally
% whole lot of nuthin

Rsquared_internalizing_RAnt = level_one_model_rh_ant.Rsquared.Ordinary - mdl4_rh_ant.Rsquared.Ordinary;
Rsquared_internalizing_LAnt = level_one_model_lh_ant.Rsquared.Ordinary - mdl3_lh_ant.Rsquared.Ordinary;
Rsquared_internalizing_RPost = level_one_model_rh_post.Rsquared.Ordinary - mdl2_rh_post.Rsquared.Ordinary;
Rsquared_internalizing_LPost = level_one_model_lh_post.Rsquared.Ordinary - mdl1_lh_post.Rsquared.Ordinary;

%% New set of models for Trilevel factors
% Still going to do this hierarchically. If I find an effect for one or all
% of the tri-level factor scores then I'm going to need to put them into
% the same model with DSM to make sure their accounting for variance in the
% data above and beyond DSM diagnoses

trilevel_regressors = [BrainMAPD_subject_struct_table.GenDis, BrainMAPD_subject_struct_table.Anhedon, BrainMAPD_subject_struct_table.Fears, BrainMAPD_subject_struct_table.Sex, BrainMAPD_subject_struct_table.eTIV];

level_two_trilevel_rh_ant = fitlm(trilevel_regressors, rh_ant_insula, 'Categorical', logical([0,0,0,1,0]))
level_two_trilevel_lh_ant = fitlm(trilevel_regressors, lh_ant_insula, 'Categorical', logical([0,0,0,1,0]))
level_two_trilevel_rh_post = fitlm(trilevel_regressors, rh_post_insula, 'Categorical', logical([0,0,0,1,0]))
level_two_trilevel_lh_post = fitlm(trilevel_regressors, lh_post_insula, 'Categorical', logical([0,0,0,1,0]))

Rsquared_trilevel_symptom_RAnt = level_one_model_rh_ant.Rsquared.Ordinary - level_two_trilevel_rh_ant.Rsquared.Ordinary;
Rsquared_trilevel_symptom_LAnt = level_one_model_lh_ant.Rsquared.Ordinary - level_two_trilevel_lh_ant.Rsquared.Ordinary;
Rsquared_trilevel_symptom_RPost = level_one_model_rh_post.Rsquared.Ordinary - level_two_trilevel_rh_post.Rsquared.Ordinary;
Rsquared_trilevel_symptom_LPost = level_one_model_lh_post.Rsquared.Ordinary - level_two_trilevel_lh_post.Rsquared.Ordinary;

%% Now to model everything together to show trilevel symptoms predict above and beyond DSM

whole_model_regressors = [BrainMAPD_subject_struct_table.Dep_lifetime, BrainMAPD_subject_struct_table.Anx_lifetime, BrainMAPD_subject_struct_table.Comorbid_lifetime, BrainMAPD_subject_struct_table.GenDis, BrainMAPD_subject_struct_table.Anhedon, BrainMAPD_subject_struct_table.Fears, BrainMAPD_subject_struct_table.Sex, BrainMAPD_subject_struct_table.eTIV];

overall_mdl_rh_ant = fitlm(whole_model_regressors, rh_ant_insula, 'Categorical', logical([1,1,1,0,0,0,1,0,]))

%% in level_two_trilevel_rh_post we see a significant effect of BrainMAPD_subject_struct_table.Anhedon on R posterior insula volume (t = -2.0595, p = .040447)
% This could make sense. Need to dig into what anhedonia really is. But now
% we're going to relate this specific ROI to Anxious-Misery to see if there
% is a relationship with a symptom specific profile. 
anx_misery_regressors = [BrainMAPD_subject_struct_table.Arousal, BrainMAPD_subject_struct_table.Sex, BrainMAPD_subject_struct_table.eTIV];

anx_misery_mdl_rh_post = fitlm(anx_misery_regressors, rh_post_insula, 'Categorical', logical([0,1,0]))
anx_misery_mdl_lh_post = fitlm(anx_misery_regressors, lh_post_insula, 'Categorical', logical([0,1,0]))

