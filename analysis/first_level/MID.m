% List of open inputs
nrun = X; % enter the number of runs here
jobfile = {'/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD/func_conn/MID_job.m'};
jobs = repmat(jobfile, 1, nrun);
inputs = cell(0, nrun);
for crun = 1:nrun
end
spm('defaults', 'FMRI');
spm_jobman('run', jobs, inputs{:});
