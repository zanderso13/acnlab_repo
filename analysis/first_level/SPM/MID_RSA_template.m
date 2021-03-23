% spm template code

matlabbatch{1}.spm.stats.fmri_spec.dir = '<UNDEFINED>';
matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
matlabbatch{1}.spm.stats.fmri_spec.timing.RT = 2.05;
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 8;
matlabbatch{1}.spm.stats.fmri_spec.sess.scans = '<UNDEFINED>';
matlabbatch{1}.spm.stats.fmri_spec.sess.cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {});
matlabbatch{1}.spm.stats.fmri_spec.sess.multi = '<UNDEFINED>';
matlabbatch{1}.spm.stats.fmri_spec.sess.regress = struct('name', {}, 'val', {});
matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = '<UNDEFINED>';
matlabbatch{1}.spm.stats.fmri_spec.sess.hpf = 128;
matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
matlabbatch{1}.spm.stats.fmri_spec.mthresh = -Inf;
matlabbatch{1}.spm.stats.fmri_spec.mask = {''};
matlabbatch{1}.spm.stats.fmri_spec.cvi = 'None';

matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('fMRI model specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;

matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'GainAll';
matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = [1];
matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.consess{2}.tcon.name = 'LossAll';
matlabbatch{3}.spm.stats.con.consess{2}.tcon.weights = [0 1];
matlabbatch{3}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.consess{3}.tcon.name = 'Gain0';
matlabbatch{3}.spm.stats.con.consess{3}.tcon.weights = [0 0 1];
matlabbatch{3}.spm.stats.con.consess{3}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.consess{4}.tcon.name = 'Loss0';
matlabbatch{3}.spm.stats.con.consess{4}.tcon.weights = [0 0 0 1];
matlabbatch{3}.spm.stats.con.consess{4}.tcon.sessrep = 'none';

% When I do all trial types
% matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'Reward5Success';
% matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = [1];
% matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
% matlabbatch{3}.spm.stats.con.consess{2}.tcon.name = 'Reward5Failure';
% matlabbatch{3}.spm.stats.con.consess{2}.tcon.weights = [0 1];
% matlabbatch{3}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
% matlabbatch{3}.spm.stats.con.consess{3}.tcon.name = 'Reward150Success';
% matlabbatch{3}.spm.stats.con.consess{3}.tcon.weights = [0 0 1];
% matlabbatch{3}.spm.stats.con.consess{3}.tcon.sessrep = 'none';
% matlabbatch{3}.spm.stats.con.consess{4}.tcon.name = 'Reward150Loss';
% matlabbatch{3}.spm.stats.con.consess{4}.tcon.weights = [0 0 0 1];
% matlabbatch{3}.spm.stats.con.consess{4}.tcon.sessrep = 'none';
% matlabbatch{3}.spm.stats.con.consess{5}.tcon.name = 'Reward0Success';
% matlabbatch{3}.spm.stats.con.consess{5}.tcon.weights = [0 0 0 0 1];
% matlabbatch{3}.spm.stats.con.consess{5}.tcon.sessrep = 'none';
% matlabbatch{3}.spm.stats.con.consess{6}.tcon.name = 'Reward0Failure';
% matlabbatch{3}.spm.stats.con.consess{6}.tcon.weights = [0 0 0 0 0 1];
% matlabbatch{3}.spm.stats.con.consess{6}.tcon.sessrep = 'none';
% 
% matlabbatch{3}.spm.stats.con.consess{7}.tcon.name = 'Loss5Success';
% matlabbatch{3}.spm.stats.con.consess{7}.tcon.weights = [0 0 0 0 0 0 1];
% matlabbatch{3}.spm.stats.con.consess{7}.tcon.sessrep = 'none';
% matlabbatch{3}.spm.stats.con.consess{8}.tcon.name = 'Loss5Failure';
% matlabbatch{3}.spm.stats.con.consess{8}.tcon.weights = [0 0 0 0 0 0 0 1];
% matlabbatch{3}.spm.stats.con.consess{8}.tcon.sessrep = 'none';
% matlabbatch{3}.spm.stats.con.consess{9}.tcon.name = 'Loss150Success';
% matlabbatch{3}.spm.stats.con.consess{9}.tcon.weights = [0 0 0 0 0 0 0 0 1];
% matlabbatch{3}.spm.stats.con.consess{9}.tcon.sessrep = 'none';
% matlabbatch{3}.spm.stats.con.consess{10}.tcon.name = 'Loss150Loss';
% matlabbatch{3}.spm.stats.con.consess{10}.tcon.weights = [0 0 0 0 0 0 0 0 0 1];
% matlabbatch{3}.spm.stats.con.consess{10}.tcon.sessrep = 'none';
% matlabbatch{3}.spm.stats.con.consess{11}.tcon.name = 'Loss0Success';
% matlabbatch{3}.spm.stats.con.consess{11}.tcon.weights = [0 0 0 0 0 0 0 0 0 0 1];
% matlabbatch{3}.spm.stats.con.consess{11}.tcon.sessrep = 'none';
% matlabbatch{3}.spm.stats.con.consess{12}.tcon.name = 'Loss0Failure';
% matlabbatch{3}.spm.stats.con.consess{12}.tcon.weights = [0 0 0 0 0 0 0 0 0 0 0 1];
% matlabbatch{3}.spm.stats.con.consess{12}.tcon.sessrep = 'none';

matlabbatch{3}.spm.stats.con.delete = 1;