%% Before the dcm2bids,need to delete unnecessary files in all image folders;
workingDir = '~/2020-SCZ';
ddd = [workingDir];

fileFolder = fullfile(ddd, '/control');
dirOutput = dir(fullfile(fileFolder));
dirOutput(1:2) = [];
subNames = {dirOutput.name};
subNum = length(subNames);
% new folder to contain the only .dcm image files of patients;
mkdir([ddd, filesep], 'controlRaw');

for i = 1: subNum
    mkdir([ddd, '/controlRaw'], subNames{i});
    subiDir = [ddd, '/control', filesep, subNames{i}];
    subiFolder = fullfile(subiDir);
    dirOutput1 = dir(fullfile(subiFolder));
    dirOutput1(1:2) = [];
    fileNames = {dirOutput1.name};
    fileNum = length(fileNames);
    for j = 1: fileNum
        mkdir([ddd, '/controlRaw', filesep, subNames{i}], fileNames{j});
        AimDcmFolder = [ddd, '/controlRaw', filesep, subNames{i},...
            filesep, fileNames{j}];
        dcmfolder = [ddd, '/control', filesep, subNames{i},...
            filesep, fileNames{j}, filesep];
        fname = dir([dcmfolder, '*.dcm']);
        if(length(fname)==1)
           copyfile([dcmfolder, fname(1).name], [AimDcmFolder,...
                '/0001.dcm']);
        end
    end
end

%% dcm2bids preparing;
filetext = fileread('~/DATA210725/SZ2010data/2020-SCZ/BIDS_configPatient.json');
matches = regexp(filetext,'SidecarFilename','end');
restID = matches(2)+5;
taskIDrun(1) = matches(3)+5;
taskIDrun(2) = matches(4)+5;
taskIDrun(3) = matches(5)+5;

DCMfolder = '~/DATA210725/SZ2010data/2020-SCZ/patientsRaw';
Folder1 = '~/DATA210725/SZ2010data/2020-SCZ';
dd = [Folder1];
folders = dir(DCMfolder);
folders(1:2) = [];
clear subName
clear filetextnew
for i = 1: % length(folders)
    clear EPIfolders
    clear s
    subName{i} = folders(i).name;
    EPIfolders = dir([DCMfolder,filesep,folders(i).name, filesep, '*FE_EPI*']);
    for j = 1 : length(EPIfolders)
        clear x
        x = split(EPIfolders(j).name, '_');
        s(j) = str2num(x{1});      
    end
    snew = sort(s);
    % rest    
    filetextnew = [filetext(1:restID-1), num2str(snew(1)), filetext(restID+3:end)];   
    for k = 2 : length(snew) %patiens: 3; contols: 4
        indexChange = length(num2str(snew(k-1))) - 3;
        taskIDrun = taskIDrun + indexChange;
        filetextnew = [filetextnew(1:taskIDrun(k-1)-1), ...
            num2str(snew(k)), ...
            filetextnew(taskIDrun(k-1)+3:end)];  
    end 
    % define the configration file for each subject
    configurationfile = [DCMfolder, filesep, folders(i).name, filesep, 'BIDS_config', subName{i}, '.json'];
    fid = fopen(configurationfile, 'w');
    fprintf(fid, '%s', filetextnew);
    fclose(fid);
    
    %%% dcm2bids -d patientsRaw/104/ -p 104 -c BIDS_config.json -o data/patients/104
    command = ['dcm2bids -d ', dd, '/patientsRaw/', subName{i}, ...
              ' -p ', subName{i}, ' -c ', configurationfile, ...
               ' -o ', dd, '/data/patients/', subName{i}];
    unix(command);
end

%% %%%% run preprocesssing with fmriprep
wfolders = ' ~/patients';
FolderP = '~/patients';
% license_file = '/share/inspurStorge/apps/freesurfer/license.txt';

subfolders = dir(FolderP);
subfolders(1:2) = [];

subNum = length(subfolders)-2;
for i = 1: subNum
    clear bids_dir;
    clear output_dir;
    clear subsubname;
    subname{i} = subfolders(i).name;
    bids_dir = [FolderP, filesep, subname{i}];
    subsubname = ['sub-', subname{i}];
    output_dir = [' ', FolderP, filesep, subname{i}, filesep, subsubname];
    
    % remove the first 6 volumes for the task-fmri and 5 volumes for restfMRI;
    niigzFileP = [FolderP, filesep, subname{i}, filesep, subsubname, '/func'];
    niigzFile = dir([niigzFileP, filesep, '*.nii.gz']);
    
    for a = 1
        niigzName{a} = niigzFile(a).name;
        command = ['fslroi ',  niigzFileP, filesep, niigzName{a}, ' ',...
            niigzFileP, filesep,  niigzName{a}, ' 5 201'];
        unix(command);
    end
    for a = 2: length(niigzFile)
        niigzName{a} = niigzFile(a).name;
        command = ['fslroi ', niigzFileP, filesep, niigzName{a},' ',...
            niigzFileP, filesep, niigzName{a}, ' 6 220'];
        unix(command);
    end
       command = ['fmriprep ', bids_dir, output_dir,...
              ' participant',...
              ' -w', wfolders,...
              ' --skip_bids_validation',...
              ' --output-spaces MNI152NLin2009cAsym:res-2 --fs-no-reconall'];
       unix(command);
       disp(command);
end
