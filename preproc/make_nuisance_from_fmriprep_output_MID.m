% this function takes a confounds file genervated by fmriprep, selects some 
% desired regressors, and computes and adds some other desired regressors.
% It returns a matrix of all possible nuisance regressors, as well as a 
% matrix of a reasonable subset of regressors decided up on in lab mtg -- see below.
%
% This function also saves the output back to the fmriprep confounds file,
% so they can stay with the dataset. If the regressors it is looking for
% are already in the fmriprep confounds file, this script assumes this function was already
% run and so it does not regenerate them, and just returns what was in the
% file.
%
% in a Jan 2019 Wager lab mtg, we decided that its sensible to include the
% following in a 1st level model for task data: 24 motion regressors, CSF
% (esp. a degraded/conservative CSF mask), canlab spike detection on the
% raw data, spikes for initial volumes (5 sec), spikes as determined by DVARS / RMSQ (Zscore > 2.5), framewise displacement.  We decided _not_ to
% include WM signal, as this often can contain BOLD signal of neuronal
% origin.  -- Yoni Ashar
%
% Adding inputs and code for extending the number of images affected by the
% scn_spike_detection script. Based on Power (2019), we're assuming that
% TR's immediately following a transient signal spike will be affected and
% should therefore be filtered in a similar manner to the spike itself. If
% the 'spike_additional_vols' variable is left blank, script will proceed
% as it did before, generating nuisance regressors identifying only the
% original spikes.
% --Zach Anderson (4/19/19)

function [Rfull, Rselected, n_spike_regs, framewise_displacement_final, gsr_final] = make_nuisance_covs_from_fmriprep_output_MID(fmriprep_confounds_fname, raw_img_fname, TR, spike_additional_vols)

R = readtable(fmriprep_confounds_fname, 'TreatAsEmpty', 'n/a', 'filetype', 'text', 'ReadVariableNames', true);

% replace NaNs in first row with Os
wh_replace = ismissing(R(1,:));
if any(wh_replace)
    R{1, wh_replace} = zeros(1, sum(wh_replace)); % make array of zeros of the right size
end


% compute 24 motion regs
mot_names = {'trans_x','trans_y','trans_z','rot_x','rot_y','rot_z'};

motion = R{:,mot_names};
diffs = [zeros(1,6); diff(motion)];
mo_sq = motion .^ 2;
mo_sq_diff = [zeros(1,6); diff(mo_sq)];
motion18 = [diffs]; %mo_sq mo_sq_diff];

% make a table of the 18 additional ones (ZA: 5/28/20 I've commented out 12
% regressors worth of motion. I have a theory. I'm getting weak signal on
% BrainMAPD. We're working with late adolescence though and I'm wondering
% if this isn't similar to pain research in that higher pain is correlating
% with task in ways that are killing the signal. Ann used 12 regressors for
% motion in her pipeline. Given my weak results, I'm thinking I overstepped
% in redesigning a whole pipeline and not trusting the work that had been
% done previously.

mot_names18 = [cellfun( @(x) [x '_diff'], mot_names, 'UniformOutput',false)]; % cellfun( @(x) [x '_sq'], mot_names, 'UniformOutput',false) cellfun( @(x) [x '_sq_diff'], mot_names, 'UniformOutput',false) ];
motion18 = array2table(motion18, 'VariableNames', mot_names18);

% remove previously saved motion cols 1) in case there was an error, and 2)
% so i can re-add them without conflict
R(:,contains(R.Properties.VariableNames, 'diff')) = [];
R(:,contains(R.Properties.VariableNames, 'sq')) = [];

% add them in (table join)
R = [R motion18];


% add spikes for initial volumes. drop any existing one first, i.e., if
% user changes the number of initial vols to dr op
R(:,contains(R.Properties.VariableNames, 'initial_vols')) = [];
nvols = round(5/TR);  % add spikes for the first 5 seconds
initial_vols = zeros(height(R), nvols);
for i=1:nvols
    initial_vols(i,i) = 1;
end
R = [R array2table(initial_vols)];


% before dealing w/ nuisance covs, drop any "additional spike regs". this
% should be done anyway, and should be done at this point in the script to
% address a bug in an earlier version that can still impact code if run on
% older confounds files
additional_spike_cols = contains(R.Properties.VariableNames,'nuisance_covs_additional_spikes'); 
R(:,additional_spike_cols)=[];

% has scn_session_spike_id already been run on this data? if so, results
% will be saved back into the confounds file, with variables named
% nuisance_covsXX
spike_cols = contains(R.Properties.VariableNames,'nuisance_covs'); 

if sum(spike_cols) == 0 % have not yet computed and added these
    
    % add in canlab spike detection (Mahalanobis distance)
    [g, spikes, gtrim, nuisance_covs, snr] = scn_session_spike_id(raw_img_fname, 'doplot', 0);
    % add in canlab spike detection (Mahalanobis distance)
    % TRs
    nuisance_covs = nuisance_covs{1};
    nuisance_covs(:,1) = []; %drop gtrim 
    R = [R array2table(nuisance_covs)];
    % find updated spike cols
    spike_cols = contains(R.Properties.VariableNames,'nuisance_covs'); 
else
    
    % recreate nuisance_covs and spikes variables for later use
    nuisance_covs = R{:,spike_cols};
    spikes{1} = find(sum(nuisance_covs,2)); % TODO: this will only work for single session data!
end

% make spike regs from dvars. we dont expect a reliable signal in the brain
% that tracks dvars, so less sensible to include as a parametric regressor.
% better to use to identify outliers
dvarsZ = [ 0; zscore(R.dvars(2:end))]; % first element of dvars always = 0, drop it from zscoring and set it to Z=0
dvars_spikes = find(dvarsZ > 3); % arbitrary cutoff -- Z > 2.5
    
% make regs from spike indices
dvars_spikes_regs = zeros(height(R),length(dvars_spikes));
for i=1:length(dvars_spikes)
    dvars_spikes_regs(dvars_spikes(i),i) = 1;
end

% find dvars_spike_regs that are non-redundant with mahal spikes, and
% include them. compare TRs at which spike happens
same = ismember(dvars_spikes, find(sum(R{:,spike_cols},2)));
dvars_spikes_regs(:,same) = []; % drop the redundant ones

% remove any previous dvars_spikes_regs, and add the ones i just made
dvars_cols = contains(R.Properties.VariableNames,'dvars_spikes'); 
R(:,dvars_cols) = [];
dvars_spikes_regs = array2table(dvars_spikes_regs);
R = [R dvars_spikes_regs];


% Motion can create artifacts lasting longer than the single image we
% usually account for using spike id scripts. we're also going to flag the
% following TRs, the number of which is defined by the user. If
% 'spike_additional_vols' remains unspecified, everything will proceed as
% it did before, meaning spikes will be identified and flagged in the
% creation of nuisance regressors without considering the following TRs
% Add them if user requested, for both nuisance_covs and dvars_spikes_regs

% if exist('spike_additional_vols')
%         
%     % concatenate generated spike nuisance covs and also dvars regs. We
%     % would like to flag subsequent TR's with respect to both of these
%     % measures.
%     nuisance_covs_with_timing_adjustment = [nuisance_covs, table2array(dvars_spikes_regs)];
%     spikes = [spikes{1};dvars_spikes];
%     nuisance_covs_additional_spikes = zeros(length(nuisance_covs),length(spikes)*spike_additional_vols);
%     
%     % This loop will create a separate column with ones in each row (TR) 
%     % we would like to consider a nuisance regressor
%     % Performs this function for dvars and spikes. From now on we'll
%     % consider the two as a single set of regressors
%     for i = 1:length(spikes) 
%         nuisance_covs_additional_spikes(spikes(i)+1 : spikes(i)+spike_additional_vols,(i*spike_additional_vols-(spike_additional_vols-1)):(i*spike_additional_vols)) = eye(spike_additional_vols);
%     end
% 
%     % if any spikes went beyond the end, trim it down
%     nuisance_covs_additional_spikes = nuisance_covs_additional_spikes(1:height(R),:);
% 
%     % Add the additional spikes to the larger covariance matrix
%     % if any already exist in R, drop them first so can (re)add them
%     % without issue
%     additional_spike_cols = contains(R.Properties.VariableNames,'additional_spikes'); 
%     R(:,additional_spike_cols) = []; 
%     R = [R array2table(nuisance_covs_additional_spikes)];
% end


% Now, remove redundant spike regressors
regs = R.Properties.VariableNames;
dvars_cols = contains(regs,'dvars_spikes'); 
spike_cols = contains(regs,'nuisance_covs'); 
% additional_spike_cols = contains(regs,'additional_spikes'); 
initial_vols = contains(regs,'initial_vols'); 

[duplicate_rows, ~] = find(sum(R{:, spike_cols | dvars_cols | initial_vols}, 2)>1); %  additional_spike_cols |
for i = 1:length(duplicate_rows) %This loop sets duplicate values to zero; drops them later (to keep indices the same during the loop)
    [~,curr_cols] = find(R{duplicate_rows(i),:}==1);
    R{duplicate_rows(i), curr_cols(2:end)} = 0;
end
R = R(1:length(nuisance_covs), any(table2array(R)));

% Select a subset of regressors to return for use in GLM to return to user
% regs = R.Properties.VariableNames;
% dvars_cols = contains(regs,'dvars_spikes'); 
% spike_cols = contains(regs,'nuisance_covs'); 
% motion_cols = contains(regs,'rot') | contains(regs,'trans') | contains(regs,'diff'); 
% % additional_spike_cols = contains(regs,'nuisance_covs_additional_spikes'); 
% 
if size(R,2)==length(regs) == 0
    % Now, remove redundant spike regressors
    regs = R.Properties.VariableNames;
    dvars_cols = contains(regs,'dvars_spikes'); 
    spike_cols = contains(regs,'nuisance_covs'); 
    additional_spike_cols = contains(regs,'additional_spikes'); 
    initial_vols = contains(regs,'initial_vols'); 

    [duplicate_rows, ~] = find(sum(R{:, spike_cols | dvars_cols | additional_spike_cols | initial_vols}, 2)>1);
    for i = 1:length(duplicate_rows) %This loop sets duplicate values to zero; drops them later (to keep indices the same during the loop)
        [~,curr_cols] = find(R{duplicate_rows(i),:}==1);
        R{duplicate_rows(i), curr_cols(2:end)} = 0;
    end
    R = R(1:length(nuisance_covs), any(table2array(R)));

    % Select a subset of regressors to return for use in GLM to return to user
    % regs = R.Properties.VariableNames;
    dvars_cols = contains(regs,'dvars_spikes'); 
    spike_cols = contains(regs,'nuisance_covs'); 
end
    
    
% Rselected = R(:,spike_cols | dvars_cols); %| additional_spike_cols);
motion_final = [R.trans_x,R.trans_y,R.trans_z,R.rot_x,R.rot_y,R.rot_z,R.trans_x_derivative1,R.trans_y_derivative1,R.trans_z_derivative1,R.rot_x_derivative1,R.rot_y_derivative1,R.rot_z_derivative1]; %,R.trans_x_power2,R.trans_y_power2,R.trans_z_power2,R.rot_x_power2,R.rot_y_power2,R.rot_z_power2,R.trans_x_derivative1_power2,R.trans_y_derivative1_power2,R.trans_z_derivative1_power2,R.rot_x_derivative1_power2,R.rot_y_derivative1_power2,R.rot_z_derivative1_power2];
motion_final = array2table(motion_final);motion_final.Properties.VariableNames = {'trans_x','trans_y','trans_z','rot_x','rot_y','rot_z','trans_x_derivative1','trans_y_derivative1','trans_z_derivative1','rot_x_derivative1','rot_y_derivative1','rot_z_derivative1'}; %,'trans_x_power2','trans_y_power2','trans_z_power2','rot_x_power2','rot_y_power2','rot_z_power2','trans_x_derivative1_power2','trans_y_derivative1_power2','trans_z_derivative1_power2','rot_x_derivative1_power2','rot_y_derivative1_power2','rot_z_derivative1_power2'};
framewise_displacement_final = R.framewise_displacement;
framewise_displacement_final = array2table(framewise_displacement_final);framewise_displacement_final.Properties.VariableNames = {'framewise_displacement'};
csf_final = R.csf;
csf_final = array2table(csf_final);csf_final.Properties.VariableNames = {'csf'};
gsr_final = R.global_signal;
gsr_final = array2table(gsr_final);gsr_final.Properties.VariableNames = {'global_signal'};

Rselected = [motion_final,csf_final,R(:,spike_cols)]; %,Rselected, framewise_displacement_final];
% compute and output how many spikes total
n_spike_regs = sum(spike_cols) % dvars_cols | additional_spike_cols)
n_spike_regs_percent = n_spike_regs / height(R)


% write back to file
writetable(R, fmriprep_confounds_fname, 'filetype', 'text');

Rfull = R;

end
