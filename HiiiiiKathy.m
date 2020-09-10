% Outline
% 1. GICA
% 2. Manual selection of components
% 3. Run gig-ica
% 4. Load in .mat files that contain component information
% 5. load in clinical information
% 6. remove brain data for subs that don't have clinical data
% 7. Visualize what we have
% 8. regress function with clinical data and brain data, may also do this
% in the gift software
% 9. Once we've visualized where junk is, then support vector regression,
% to predict level of trilevel symptom

% 1.
datadir = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD/func_conn/ICA/gigica_output';
fnames = filenames(fullfile(datadir, 'gica_cmd_ica_br*mat'));
for sub = 1:length(fnames)
    load(fnames{sub});
    dat1(:,sub) = compSet.ic(1,:)';
    dat2(:,sub) = compSet.ic(2,:)';
    dat4(:,sub) = compSet.ic(4,:)';
    dat5(:,sub) = compSet.ic(5,:)';
    dat6(:,sub) = compSet.ic(6,:)';
    dat12(:,sub) = compSet.ic(12,:)';
    dat17(:,sub) = compSet.ic(17,:)';
    dat21(:,sub) = compSet.ic(21,:)';   
end

%% 5.
% load in sublist.txt to reference subject ID's
txtdir = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD/func_conn/ICA';
txt_fname = fullfile(txtdir, 'sublist.mat');
load(txt_fname);
sub_list = sublist.datassub10242 * -1;

% Load in clinical data for group analysis
clinicaldir = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD/clinical_data';
med_dir = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD';

load(fullfile(clinicaldir,'trilevel_longitudinal.mat'));
load(fullfile(med_dir,'Medication_T2.mat'));

trilevel_array = [trilevel_T1.ID,trilevel_T1.GenDis,trilevel_T1.Anhedonia,trilevel_T1.Fears,trilevel_T1.Narrow];
med_array = [T2MedicationInventory3(:,1),T2MedicationInventory3(:,88),T2MedicationInventory3(:,89)];
med_array = table2array(med_array);

for sub = 1:length(sub_list)
    PID(sub,1) = sub_list(sub);
end        


for sub = 1:length(sub_list)

    if isempty(find(trilevel_T1.PID(:) == sub_list(sub,1))) == 0
        curr = find(trilevel_T1.PID(:) == sub_list(sub,1));
        trilevel_regressors(sub,:) = trilevel_T1(curr,:);
    else
        % for OFC loss_ppi_fnames{sub}(100:104); for HO_VMPFC loss_ppi_fnames{sub}(105:109)
        disp(strcat(num2str(PID(sub,1)), ' missing clinical info')) 
        trilevel_regressors(sub,:) = NaN;
        trilevel_regressors(sub,1) = PID(sub,1);
    end
    if isempty(find(med_array(:,1) == PID(sub,1))) == 0
        curr2 = find(med_array(:,1) == PID(sub,1));
        med_regressors(sub,:) = [med_array(curr2,2),med_array(curr2,3)];
    else
        disp(strcat(num2str(PID(sub,1)), ' missing medication info'))
        med_regressors(sub,:) = [0,0];
    end
    
end

site_regressors = (PID < 20000);

GenDis = trilevel_regressors(:,2);
Anhedonia = trilevel_regressors(:,3);
Fears = trilevel_regressors(:,4);
Narrow = trilevel_regressors(:,5);


%% SVM model

GenDis_mdl = fitrsvm(dat1',GenDis);

%% Cross validate the model generated above

CV_GenDis_mdl = crossval(GenDis_mdl,'KFold',10);

%% Once you have a cross validated model, predict responses based on that model

yfit = kfoldPredict(CV_GenDis_mdl);

%%
for sub = 1:length(clinical_info.PID)
    curr = find(clinical_info.PID



