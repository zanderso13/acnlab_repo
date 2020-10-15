% I'm dumb and deleted a bunch of files. Thanksfully I had them backed up
% but now I need to put them back where they belong
basedir = '/home/zaz3744/projects';
% Starting folder
datadir = '/projects/p30954/Schizconnect/bids';
% Ending folder
sortdir = '/projects/p30954/Schizconnect_final/COBRE_BIDS';
cd(basedir)

files_to_move = filenames(fullfile(datadir, '*/*/*.nii.gz'));
anat_list = filenames(fullfile(sortdir,'*/*/anat/*echo-01_T1w.nii.gz'));
func_list = filenames(fullfile(sortdir,'sub*/ses*/*.nii.gz'));
for file = 1:length(func_list)
    subid = anat_list{file}(80:92);
    sesid = anat_list{file}(94:105);
    %movefromfolder1 = filenames(fullfile(datadir,strcat(subid,'*json*'),'func'));
    %movefromfolder2 = filenames(fullfile(datadir,strcat(subid,'*nii.gz'),'func/*'));
    %movetofolder1 = fullfile(sortdir,subid,sesid);
    filetomove = func_list{file};
    movetofolder2 = fullfile(sortdir,subid,sesid,'func');
    keyboard
    if exist(movetofolder2) == 0
        if isempty(movefromfolder1)==0
            %keyboard
            %movefile(movefromfolder2{1}, movetofolder1)
            %keyboard
            %movefile(movefromfolder1{1}, movetofolder2) 
            %disp(subid)
        else
            %disp(subid)
        end
    
    end
end


