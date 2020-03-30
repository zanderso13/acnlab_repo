

ftpserver = ftp('ftp.box.com', 'briankraus2024@u.northwestern.edu', 'btklm555!@');
cd '/projects/b1081/Test_Transfer'
cd(ftpserver,'DATA/HCP/BOLD_DATA/100206_rfMRI_REST1_LR');
ftpdir = dir(ftpserver);

%hs = struct(ftpserver);
%hs.jobject.setConnectionTimeout(50);

for x = 3:length(ftpdir)    %% Can only do 1 at a time in MATLAB

    mget(ftpserver,ftpdir(x).name)
    
end
    
close(ftpserver)



system(['ftp -inv ftp.box.com <<EOF;'...
       'lcd /projects/b1081/Test_Transfer;'...
       'cd DATA/HCP/BOLD_DATA/100206_rfMRI_REST1_LR;'...
       'mget *;'...
       'quit;'...
       'EOF'])
   
   
   
% Run BASH script
   
pathToScript = fullfile('/home/btk2142/MATLABFTPTransfer.sh');  % assumes script is in curent directory
inputfilepath = '/DATA/HCP/BOLD_DATA/100206_rfMRI_REST1_LR';
%inputfolder = '100206_rfMRI_REST1_LR';
outputfilepath = '/projects/b1081/Test_Transfer';
%cmdStr = [pathToScript ' ' inputfilepath ' ' outputfilepath ' ' inputfolder];
cmdStr = [pathToScript ' ' inputfilepath ' ' outputfilepath];
system(cmdStr);
   

