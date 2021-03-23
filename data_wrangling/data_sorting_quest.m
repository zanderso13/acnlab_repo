% I'm dumb and deleted a bunch of files. Thanksfully I had them backed up
% but now I need to put them back where they belong
sortdir = '/projects/b1108/data/BrainMAPD/NEED_TO_SORT/T4_rest_to_sort';
datadir = '/projects/b1108/data/BrainMAPD';
cd(sortdir)

files_to_move = filenames('*/ses-4/func/*');
for file = 1:length(files_to_move)
    curr = files_to_move{file};
    subid = curr(5:9);
    sess = curr(19);
    movetofolder = filenames(fullfile(datadir,strcat('*',subid,'/*',sess,'/func')));
    movefile(curr, movetofolder{1});
end