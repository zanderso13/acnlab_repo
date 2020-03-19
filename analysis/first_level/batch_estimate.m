% this script writes out job files for 1st levels, for submission w/ sbatch
% to toggle between bladder and acute tasks, simply change name of the
% function between run_subject_firstlevel_acute or _bladder on line 32 and
% the folder name on line 28

fl_dir = '/work/ics/data/projects/wagerlab/labdata/projects/OLP4CBP/first_level';
repodir = '/projects/zaan8774/repo';
addpath(genpath('/projects/zaan8774/repo'))
scriptdir = fullfile(fileparts(fl_dir), 'script_tmp', datestr(now));
mkdir(scriptdir)

URSIs = get_URSIs();


%% for each S, write a job submission file, then submit it

cd(scriptdir); %slurm out files will be saved here

for i=1:length(URSIs)
    
    % to make jobs, to be submitted w/ sbatch
    overwrite = 1;
    ses = 2;
    
    if ~overwrite % if not overwriting, can skip right here, before job submission
        [subID, ~, ~, sub_str, ses_str] = get_subjID_from_URSI(URSIs{i});    
        if exist(fullfile(fl_dir, sub_str, ses_str(ses,:), 'sponpain','SPM.mat'),'file') % skip if SPM.mat exists
            fprintf('skipping %d\n', i)
            continue
        end
    end
    
    
    % to run directly (interactively)
    %run_subject_firstlevel_acute(URSIs{i}, ses, 1)    
    %continue


       s = ['#!/bin/bash\n\n'...
    '#SBATCH -p blanca-ics\n'...
    '#SBATCH -n 1\n'...
    '#SBATCH --time 2:00:00\n'...  
    '#SBATCH --mem=30G\n\n'...
    'matlab -nodisplay -nosplash -nodesktop -r "addpath(genpath(''' repodir ''')); run_subject_firstlevel_sponpain(''' URSIs{i} ''', ' num2str(ses) ',' num2str(overwrite) '); quit"\n\n'];
  
    scriptfile = fullfile(scriptdir, 'first_level_script.sh');
    fout = fopen(scriptfile, 'w');
    fprintf(fout, s);
    
    !sbatch first_level_script.sh
end