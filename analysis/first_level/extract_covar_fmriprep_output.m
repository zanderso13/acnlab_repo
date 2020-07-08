basedir = '/home/zach/Documents/ACNlab/current_projects/func_conn';
datadir = 'data/';

fnames = filenames(fullfile(basedir,datadir,'*tsv'));
for sub = 1:length(fnames)
    fid = fopen(fnames{sub});
    data = textscan(fid,'%s%s%s');
    fclose(fid);
    keyboard
end