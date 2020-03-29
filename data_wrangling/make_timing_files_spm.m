% I have lots of timing files from Ann (who's amazing as we all know). But
% they're split into consumption and anticipation and I don't want to go
% through what I am right now in FSL. So I need to concatenate all subject
% timing files so that I can define all contrasts at once in spm


datadir = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD/Oldham_ROI_by_diagnosis/TimingFiles_082218';
antdir = '/FSL_anticipation_072418';
condir = '/FSL_consumption_110818';
savedir = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD/func_conn/final_timing_files/run-1';

fnames_antloss_any = filenames(fullfile(datadir,antdir,'*/*Loss_Run1.txt'));
fnames_antloss_zero = filenames(fullfile(datadir,antdir,'*/*Loss0_Run1.txt'));
fnames_antgain_any = filenames(fullfile(datadir,antdir,'*/*Win_Run1.txt'));
fnames_antgain_zero = filenames(fullfile(datadir,antdir,'*/*Win0_Run1.txt'));
fnames_hitconloss = filenames(fullfile(datadir,condir,'*/*Hit*Loss_Feedback_Run1.txt'));
fnames_hitcongain = filenames(fullfile(datadir,condir,'*/*Hit*Win_Feedback_Run1.txt'));
fnames_missconloss = filenames(fullfile(datadir,condir,'*/*Miss*Loss_Feedback_Run1.txt'));
fnames_misscongain = filenames(fullfile(datadir,condir,'*/*Miss*Win_Feedback_Run1.txt'));
fnames_motor = filenames(fullfile(datadir,antdir,'*/*Motor_Run1.txt'));

for sub = 1:length(fnames_antloss_any)
    
    % define all the different kinds of contrasts
    id = fnames_antgain_any{sub}(155:159);
    dat_antloss_any = load(fnames_antloss_any{sub});
    dat_antloss_zero = load(fnames_antloss_zero{sub});
    dat_antgain_any = load(fnames_antgain_any{sub});
    dat_antgain_zero = load(fnames_antgain_zero{sub});
    dat_hitconloss = load(fnames_hitconloss{sub});
    dat_missconloss = load(fnames_missconloss{sub});
    dat_hitcongain = load(fnames_hitcongain{sub});
    dat_misscongain = load(fnames_misscongain{sub});
    dat_motor = load(fnames_motor{sub});    
    
    if isempty(dat_antloss_any) | isempty(dat_antloss_zero) | isempty(dat_antgain_any) | isempty(dat_antgain_zero) | isempty(dat_hitconloss) | isempty(dat_missconloss) | isempty(dat_hitcongain) | isempty(dat_misscongain) | isempty(dat_motor)
        disp(id)
    else        
        onsets_temp = [dat_antloss_any(:,1);dat_antloss_zero(:,1);dat_antgain_any(:,1);dat_antgain_zero(:,1);dat_hitconloss(:,1);dat_missconloss(:,1);dat_hitcongain(:,1);dat_misscongain(:,1);dat_motor(:,1)];
        names_temp = cellstr([string(repmat('antloss',size(dat_antloss_any,1),1));string(repmat('antlosszero',size(dat_antloss_zero,1),1));string(repmat('antgain',size(dat_antgain_any,1),1));string(repmat('antgainzero',size(dat_antgain_zero,1),1));string(repmat('hitconloss',size(dat_hitconloss,1),1));string(repmat('hitcongain',size(dat_hitcongain,1),1));string(repmat('missconloss',size(dat_missconloss,1),1));string(repmat('misscongain',size(dat_misscongain,1),1));string(repmat('motor',size(dat_motor,1),1))]);
        durations_temp = ones(length(onsets_temp),1)*4;
        names = {'antloss','antlosszero','antgain','antgainzero','hitconloss','missconloss','hitcongain','misscongain','motor'};
        onsets = {onsets_temp(string(names_temp)=='antloss',:)',onsets_temp(string(names_temp)=='antlosszero',:)',onsets_temp(string(names_temp)=='antgain',:)',onsets_temp(string(names_temp)=='antgainzero',:)',onsets_temp(string(names_temp)=='hitconloss',:)',onsets_temp(string(names_temp)=='missconloss',:)',onsets_temp(string(names_temp)=='hitcongain',:)',onsets_temp(string(names_temp)=='misscongain',:)',onsets_temp(string(names_temp)=='motor',:)'};
        durations = {durations_temp(string(names_temp)=='antloss',:)',durations_temp(string(names_temp)=='antlosszero',:)',durations_temp(string(names_temp)=='antgain',:)',durations_temp(string(names_temp)=='antgainzero',:)',durations_temp(string(names_temp)=='hitconloss',:)',durations_temp(string(names_temp)=='missconloss',:)',durations_temp(string(names_temp)=='hitcongain',:)',durations_temp(string(names_temp)=='misscongain',:)',durations_temp(string(names_temp)=='motor',:)'};
        filename = fullfile(savedir,strcat(id,'.mat'));
        save(filename,'onsets','names','durations')
    end
    clear id dat_antgain dat_antgain_zero dat_antloss dat_antloss_zero dat_hitcongain dat_hitconloss dat_misscongain dat_missconloss dat_motor  filename T type duration onsets
end


% Will DEFINITELY be helpful later: T.onset(string(T.trial)=='antloss',:)

% These might be useful later. 
% fnames_antgain{1}(155:159);
% fnames_hitcongain{1}(153:157);
