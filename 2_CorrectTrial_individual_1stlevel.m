%% perform 1st level contrast analysis only selecting the correct trial
%  Run3 for the following analysis;
% Step1: resample, smoothing and highpass filtering (>=0.01Hz)using WM task fMRI after fmriprep prepro, resample 3, smooth 6mm; 
%% resample the 2*2*2 files to 3*3*3
% 3dresample -dxyz 3 3 3 -prefix output.nii.gz -input input.nii.gz
WKD = '~/DATAfMRIprepAGAIN2112';
filesP = [WKD, '/bidsfMRIprepSub/fmriprepHCrun3'];
files = dir([filesP, '/gunzipHCrun3/*.nii']);
for i = 1:length(files)
    newname = ['R3mm', files(i).name];
    command = ['3dresample -dxyz 3 3 3 -prefix ', filesP, '/R3mmgunzipHCrun3/', newname,...
        ' -input ', filesP, '/gunzipHCrun3/', files(i).name];
    unix(command);
end
% sz
WKD = '~/DATAfMRIprepAGAIN2112';
filesP = [WKD, '/bidsfMRIprepSub/fmriprepSZrun3'];
files = dir([filesP, '/gunzipSZrun3/*.nii']);
for i = 1:length(files)
    newname = ['R3mm', files(i).name];
    command = ['3dresample -dxyz 3 3 3 -prefix ', filesP, '/R3mmgunzipSZrun3/', newname,...
        ' -input ', filesP, '/gunzipSZrun3/', files(i).name];
    unix(command);
end
%% smooth the files from step4 after resample
WKD = '~/DATAfMRIprepAGAIN2112';
filesP = [WKD, '/bidsfMRIprepSub/fmriprepSZrun3'];
mkdir([filesP, filesep], 'S6R3mmgunzipSZrun3') % 6mm smooth folder
files = dir([filesP, '/R3mmgunzipSZrun3/*.nii']);
for i = 1:length(files)
    disp(['smooth 6mm processing the ', num2str(i), '-th subject...']);
    newname = ['S6', files(i).name];
    command = ['3dmerge -1blur_fwhm 6 -doall -prefix ', filesP, '/S6R3mmgunzipSZrun3/', newname,...
        ' ', filesP, '/R3mmgunzipSZrun3/', files(i).name];
    unix(command);
end
% hc
WKD = '/home1/Fengnana/DATA210725/SZ2010data/2020-SCZ/DATAfMRIprepAGAIN2112';
filesP = [WKD, '/bidsfMRIprepSub/fmriprepHCrun3'];
mkdir([filesP, filesep], 'S6R3mmgunzipHCrun3') % 6mm smooth folder
files = dir([filesP, '/R3mmgunzipHCrun3/*.nii']);
for i = 1:length(files)
    disp(['smooth 6mm processing the ', num2str(i), '-th subject...']);
    newname = ['S6', files(i).name];
    command = ['3dmerge -1blur_fwhm 6 -doall -prefix ', filesP, '/S6R3mmgunzipHCrun3/', newname,...
        ' ', filesP, '/R3mmgunzipHCrun3/', files(i).name];
    unix(command);
end
%% filtering the images Again;
% fslmaths sXXX.nii.gz -Tmean sXXX_mean
% fslmaths sXXX.nii.gz -bptf 25 -1 -add XXX_mean fsXXX.nii.gz
% HC
filePath = [WKD, '/S6R3mmgunzipHCrun3']; 
meanPath = [WKD, '/F100smeanHCrun3'];
aimPath = [WKD, '/F100sS6R3mmgunzipHCrun3A'];
subfiles = dir([filePath, '/*.nii']);
% filtering 100s
for i = 1: length(subfiles)
    display(['Computing the ', num2str(i), '-th subjects!']);
    filename = subfiles(i).name;
    meanname = ['mean_', filename];
    newname = ['F_', filename];
    command = ['fslmaths ', filePath, '/', filename, ' -Tmean ',...
        meanPath, '/', meanname];
    unix(command);
    % hp_sigma = 100 / TR / 2;
    command = ['fslmaths ', filePath, '/', filename, ' -bptf 25 -1 -add ',...
        meanPath, '/', meanname, ' ', aimPath, '/', newname];
    unix(command);
end
% unzip file
filePath = [WKD, '/F100sS6R3mmgunzipHCrun3A']; 
files = dir([filePath, '/*.nii.gz']);
aimPath = [WKD, '/F100sS6R3mmgunzipHCrun3Aun'];
for i = 1: length(files)
    splitCell = split(files(i).name, '.');
    gunzipN = [splitCell{1}, '.', splitCell{2}];
    gunzip([filePath, '/', files(i).name]);
    movefile([filePath, '/', gunzipN], ...
            [aimPath, '/', gunzipN]);
end
% SZ
filePath = [WKD, '/S6R3mmgunzipSZrun3']; 
meanPath = [WKD, '/F100smeanSZrun3'];
aimPath = [WKD, '/F100sS6R3mmgunzipSZrun3A'];
subfiles = dir([filePath, '/*.nii']);
for i = 1: length(subfiles)
    display(['Computing the ', num2str(i), '-th subjects!']);
    filename = subfiles(i).name;
    meanname = ['mean_', filename];
    newname = ['F_', filename];% *.nii
    command = ['fslmaths ', filePath, '/', filename, ' -Tmean ',...
        meanPath, '/', meanname]; % result mean_*.nii.gz
    unix(command);
    % hp_sigma = 100 / TR / 2;
    command = ['fslmaths ', filePath, '/', filename, ' -bptf 25 -1 -add ',...
        meanPath, '/', meanname, ' ', aimPath, '/', newname];
    unix(command); % result F_*.nii.gz
    
    gunzip([aimPath, '/', [newname, '.gz']]);
    movefile([aimPath, '/', newname], ...
            [WKD, '/F100sS6R3mmgunzipSZrun3Aun/', newname]);
end

%% step2.prepare for rp.txt---extract the motion tranX, tranY, tranZ, rotX, rotY, rotZ 
% WM, CSF, as the .txt file for each sub;
% HC
WKD = '~/bidsfMRIprepSub';
filesP = [WKD, '/DATA_FS6R3load1n10RightTrailRun3/fmriprepHCrun3/confounds_timseries'];

seriesfile = dir([filesP, '/*.tsv']);
for subi =1: length(seriesfile)
    clear tsvname;
    clear newname;
    clear tabletran;
    clear motiondataN;
    display(['Computing the ', num2str(subi), '-th subjects!']);
    tsvname = seriesfile(subi).name;
    tsvnamesplit = split(tsvname, '_');
    newname = ['rp', tsvnamesplit{1}, '_', tsvnamesplit{3},...
        '_timeseries.txt'];  % newname for txt
    confound = importdata([filesP, filesep, tsvname]);%struct
    tabletran = array2table(confound.textdata);%table
    [m, n] = size(tabletran);
    cellcolu = {'trans_x', 'trans_y', 'trans_z', 'rot_x', 'rot_y', 'rot_z',...
        'csf', 'white_matter'};
    for co =1: length(cellcolu)
        for cxyz = 1: n
            if strcmp(cellcolu{co}, tabletran{1, cxyz})
                number(subi, co) = cxyz; % find the column number of every motion data
            end
        end
    end
    motiondata = [tabletran{2:221, number(subi, 1)}, tabletran{2:221, number(subi, 2)}, tabletran{2:221,number(subi, 3)},...
        tabletran{2:221, number(subi, 4)}, tabletran{2:221, number(subi, 5)}, tabletran{2:221, number(subi, 6)},...
        tabletran{2:221, number(subi, 7)}, tabletran{2:221, number(subi, 8)}]; % cell
    % cell to num
    for i=1: 220
        for j=1: 8 % 6 motion variables+csf+wm;
            motiondataN(i, j)  = str2num(motiondata{i, j});
        end
    end
    save([filesP,'/confoundtimestxt/',  newname], 'motiondataN', '-ascii');
end
%% step3.extract the onset duration for each subject;
%% set event times table files
% 1.select the trial data for the included subject here;
% WKD = working space;
rawPath = [WKD, '/control'];
goalPath = [WKD, '/HC26'];
txtfile = dir([rawPath, '/*.txt']);
for i=31: 56
    % read the sub ID enrolled in this study;
    ID = SZ30HC26matchRightTrial(i, 1);
    for j = 1: length(txtfile)
        txtname = txtfile(j).name;
        txtnum = txtname(8:10);
        if ID == str2num(txtnum)
            copyfile([rawPath, filesep, txtname],...
                [goalPath, filesep, [txtnum, '.txt']]);
        end
    end
end
% SZ
rawPath = [WKD, '/patients'];
goalPath = [WKD, '/SZ30'];
txtfile = dir([rawPath, '/*.txt']);
for i=1: 30
    display(['Computing the ', num2str(i), '-th subjects!']);
    ID = SZ30HC26matchRightTrial(i, 1);
    for j = 1: length(txtfile)
        txtname = txtfile(j).name;
        txtnum = txtname(8:10);
        if ID == str2num(txtnum)
            copyfile([rawPath, filesep, txtname],...
                [goalPath, filesep, [txtnum, '.txt']]);
        end
    end
end
% 2. read the trail txt and construct the individual event files;
WKD = ['~/bidsfMRIprepSub/DATA_FS6R3load1n10RightTrailRun3'];
load([WKD, '/Run3eventTable.mat']);
load([WKD, '/eventNames.mat']);
rawPath = [WKD, '/digitWM/SZ30'];
file  = dir([rawPath, '/*.txt']);
for i = 1: length(file)
    display(['Computing the ', num2str(i), '-th subjects!']);
    eventFile = importdata([rawPath, filesep, file(i).name]);%struct
    IDcondition(1, i) = str2num(file(i).name(1:3)); % save ID
    subicomb = [onlyRun3eventTimes{:, 1:9}, eventFile.data(81:120, 1:3)];
    subicombta = array2table(subicomb,'VariableNames',...
        {'Run', 'condition','interTime','trialStart','trialfinal',...
        'delayOnset','probeOnset','ITIOnset','trialorder',...
        'probeACC', 'probeRT', 'shumu'});
    subicomA = subicombta(find(subicombta{:, 10}>0), :);% select the correct trail;
    subicomAsort = sortrows(subicomA, 2); % ascending sort condition;
    num = unique(subicomAsort{:, 2}, 'rows'); % [0;1;2;3;4]
    IDcondition(2: length(num)+1, i) = num; % save condition ;
    % condition 0/1/2/3/4
    for condi =0: length(num)-1
        subicom0 = subicomAsort(find(subicomAsort{:, 2}==num(condi+1)), :);
        [trailsum, n] = size(subicom0);% trailsum: the number of trail of special conndit
        colu = [4,6,7]; % 
        for stage = 1:3
            onsetset = subicom0{:, colu(stage)}';
            onsets{1, condi*3+stage} = onsetset;
        end
        % duration
        dur = [2, 4, 2]; % s
        for sta = 1:3
            durat = repmat(dur(sta), 1, trailsum);
            durations{1, condi*3+sta} = durat;
        end
    end
    IDname = [file(i).name(1:3), '_eventtimes.mat'];
    save([WKD, '/fmriprepSZrun3/eventtimes/', IDname],...
        'onsets', 'durations', 'names');
end
% HC
WKD = ['~/bidsfMRIprepSub/DATA_FS6R3load1n10RightTrailRun3'];
load([WKD, '/Run3eventTable.mat']);
load([WKD, '/eventNames.mat']);
rawPath = [WKD, '/digitWM/HC27'];
file  = dir([rawPath, '/*.txt']);
for i = 1: length(file)
    display(['Computing the ', num2str(i), '-th subjects!']);
    eventFile = importdata([rawPath, filesep, file(i).name]);%struct
    IDcondition(1, i) = str2num(file(i).name(1:3)); % save ID
    subicomb = [onlyRun3eventTimes{:, 1:9}, eventFile.data(81:120, 1:3)];
    subicombta = array2table(subicomb,'VariableNames',...
        {'Run', 'condition','interTime','trialStart','trialfinal',...
        'delayOnset','probeOnset','ITIOnset','trialorder',...
        'probeACC', 'probeRT', 'shumu'});
    subicomA = subicombta(find(subicombta{:, 10}>0), :);% select the correct trail;
    subicomAsort = sortrows(subicomA, 2); % ascending sort condition;
    num = unique(subicomAsort{:, 2}, 'rows'); % [0;1;2;3;4]
    IDcondition(2: length(num)+1, i) = num; % save condition ;
    % condition 0/1/2/3/4
    for condi =0: length(num)-1
        subicom0 = subicomAsort(find(subicomAsort{:, 2}==num(condi+1)), :);
        [trailsum, n] = size(subicom0);% trailsum: the number of trail of special conndit
        colu = [4,6,7]; % 
        for stage = 1:3
            onsetset = subicom0{:, colu(stage)}';
            onsets{1, condi*3+stage} = onsetset;
        end
        % duration
        dur = [2, 4, 2]; % s
        for sta = 1:3
            durat = repmat(dur(sta), 1, trailsum);
            durations{1, condi*3+sta} = durat;
        end
    end
    IDname = [file(i).name(1:3), '_eventtimes.mat'];
    save([WKD, '/fmriprepHCrun3/eventtimes/', IDname],...
        'onsets', 'durations', 'names'); % save as struct;
end

%% 1st level Contrast copmutation -- contrast Analysis-2022-4-16
%%%%%%%%%%%%%%%%%%%%% sz %%%%%%%%%%%%%%%%%%%%%
WKD = ['~/bidsfMRIprepSub/DATA_FS6R3load1n10RightTrailRun3'];
rpdirFiles = [WKD, '/fmriprepSZrun3/confounds_timseries/confoundtimestxt']; % motion para
rpfile = dir([rpdirFiles, '/rpsub-*.txt']);  % rp files for each sub;
datadir = [WKD, '/fmriprepSZrun3/F100sS6R3mmgunzipSZrun3Aun'];  %preprocessing results    
subfile = dir([datadir, '/F_S6R3mmsub-*.nii']);
eventDir = [WKD, '/fmriprepSZrun3/eventtimes']; % event files for each sub
eventFiles = dir([eventDir, '/*_eventtimes.mat']);

spm_get_defaults              
% global defaults
spm_jobman('initcfg');

nsub=length(subfile);
% nses=5; 
ncon=23; % ncon represent how many conditions you have in model(as in .mat),including names,onsets and durations

mkdir([datadir, filesep], '1levelResults'); 
cwd = [datadir, filesep, '1levelResults']; % path to Results(1st-level);

for sub=1:nsub
%     if sub==2||sub==3;   % skip invaild subjects
%         continue;
%     end
    display(['Computing the ', num2str(sub), '-th subjects!']);
    subname = subfile(sub).name;
    subnamecell = split(subname, '_');
    mkdir([cwd, filesep], [subnamecell{2}, '_', subnamecell{4}]); 
    %%%%%%%%%%%%%model specification and estimation%%%%%%%%%%%%%%%
    ses = 1;
    subfolder = [cwd, filesep, [subnamecell{2}, '_', subnamecell{4}]];
    jobs{1}.stats{1}.fmri_spec.dir=cellstr(subfolder);   %results for individual model
    jobs{1}.stats{1}.fmri_spec.timing.units='secs';    %Unit for design:secs/scans
    jobs{1}.stats{1}.fmri_spec.timing.RT=2;   % RT:repetition time (TR)
    jobs{1}.stats{1}.fmri_spec.timing.fmri_t = 36;  % if performed slicetiming, the number of slice; otherwise, by default.
    jobs{1}.stats{1}.fmri_spec.timing.fmri_t0 = 18; %setting t0=t/2 is a good compromise
    
    % muldirFiles=[datadir subExpID(sub).name, '/', runnames{ses}, '/'];     %folders after preprocessing
    muldirFiles=[datadir, '/'];
    mulfiles=spm_get('Files', muldirFiles, subfile(sub).name);        %the last step files of preprossing
    mulmotionfiles=spm_get('Files', rpdirFiles, rpfile(sub).name) ;   %head movement files for each run
    % develop individual onsets file
    % mulcondition_file_name=['sub',num2str(sub),'-',num2str(ses),'.mat'];  %behavior name for earch run of each subject
    
    mulcondpath=[WKD, '/fmriprepSZrun3/eventtimes'];    %behavior path
    mulcondition_file_name = eventFiles(sub).name;  
    mulconditions=load(fullfile(mulcondpath,mulcondition_file_name));  %load behavior results(.mat) for earch run of each subject
    jobs{1}.stats{1}.fmri_spec.sess(ses).scans=cellstr(mulfiles);    
    
    for j=1:ncon-8 
            jobs{1}.stats{1}.fmri_spec.sess(ses).cond(j).name = mulconditions.names{j};
            jobs{1}.stats{1}.fmri_spec.sess(ses).cond(j).onset = mulconditions.onsets{j};
            jobs{1}.stats{1}.fmri_spec.sess(ses).cond(j).duration = mulconditions.durations{j};
    end
    
    jobs{1}.stats{1}.fmri_spec.bases.hrf.derivs = [0 0];  % add time derivative = [1 0] no time derivative  = [0 0]; if changed, change the contrast file above and below
    jobs{1}.stats{1}.fmri_spec.sess(ses).multi_reg=cellstr(mulmotionfiles); % add head movements in model
     
    resultpath=[cwd, '/', [subnamecell{2}, '_', subnamecell{4}], '/'];
    jobs{1}.stats{2}.fmri_est.spmmat=cellstr(fullfile(resultpath,'SPM.mat')); % Save model settings and estimate results in SPM.mat
    save(fullfile(resultpath,'modelspecification.mat'),'jobs'); % save jobs
    
    spm_jobman('run', jobs);
    clear jobs;
    
    %%%%%%%%%%%%%%%%%%%seting contrast conditions%%%%%%%%%%%%%%%%%%%%%%%%%
    resultpath=[cwd, '/', [subnamecell{2}, '_', subnamecell{4}], '/'];
    jobs{1}.stats{1}.con.spmmat=cellstr(fullfile(resultpath,'SPM.mat'));  % write contrast information into SPM.mat
    % set contrast matrix for contrast conditions for individual subjexts
    % develop a pipline if you have too many contrasts
    jobs{1}.stats{1}.con.consess{1}.tcon=struct('name', 'load1Conload0','convec',...
        [0 -1 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0], 'sessrep', 'none');
    jobs{1}.stats{1}.con.consess{2}.tcon=struct('name', 'load2Conload0','convec',...
        [0 -1 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0], 'sessrep', 'none');
    jobs{1}.stats{1}.con.consess{3}.tcon=struct('name', 'load3Conload0','convec',...
        [0 -1 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0], 'sessrep', 'none');
    jobs{1}.stats{1}.con.consess{4}.tcon=struct('name', 'load4Conload0','convec',...
        [0 -1 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0], 'sessrep', 'none');
    
    spm_jobman('run',jobs);            % activation map is 'ps' format, which can be opened with GhostView and can also be converted to PDF format
    clear jobs;
end    % end of individu
%% HC 1st analysis
%%%%%%%%%%%%%%%%%  hc  %%%%%%%%%%%%%%%%%%%
WKD = ['~/bidsfMRIprepSub/DATA_FS6R3load1n10RightTrailRun3'];
rpdirFiles = [WKD, '/fmriprepHCrun3/confounds_timseries/confoundtimestxt']; % motion para
rpfile = dir([rpdirFiles, '/rpsub-*.txt']);  % rp files for each sub;
datadir = [WKD, '/fmriprepHCrun3/F100sS6R3mmgunzipHCrun3Aun'];  %preprocessing results    
subfile = dir([datadir, '/F_S6R3mmsub-*.nii']);
eventDir = [WKD, '/fmriprepHCrun3/eventtimes']; % event files for each sub
eventFiles = dir([eventDir, '/*_eventtimes.mat']);

spm_get_defaults              
% global defaults
spm_jobman('initcfg');

nsub=length(subfile);
% nses=5; 
ncon=23; % ncon represent how many conditions you have in model(as in .mat),including names,onsets and durations
% runnames={'RUN1' 'RUN2' 'RUN3' 'RUN4' 'RUN5'};    
mkdir([datadir, filesep], '1levelResults'); 
cwd = [datadir, filesep, '1levelResults']; % path to Results(1st-level);

for sub=1:nsub
%     if sub==2||sub==3;   % skip invaild subjects
%         continue;
%     end
    display(['Computing the ', num2str(sub), '-th subjects!']);
    subname = subfile(sub).name;
    subnamecell = split(subname, '_');
    mkdir([cwd, filesep], [subnamecell{2}, '_', subnamecell{4}]); 
    %%%%%%%%%%%%%model specification and estimation%%%%%%%%%%%%%%%
    ses = 1;
    subfolder = [cwd, filesep, [subnamecell{2}, '_', subnamecell{4}]];
    jobs{1}.stats{1}.fmri_spec.dir=cellstr(subfolder);   %results for individual model
    jobs{1}.stats{1}.fmri_spec.timing.units='secs';    %Unit for design:secs/scans
    jobs{1}.stats{1}.fmri_spec.timing.RT=2;   % RT:repetition time (TR)
    jobs{1}.stats{1}.fmri_spec.timing.fmri_t = 36;  % if performed slicetiming, the number of slice; otherwise, by default.
    jobs{1}.stats{1}.fmri_spec.timing.fmri_t0 = 18; %setting t0=t/2 is a good compromise
    
    % muldirFiles=[datadir subExpID(sub).name, '/', runnames{ses}, '/'];     %folders after preprocessing
    muldirFiles=[datadir, '/'];
    mulfiles=spm_get('Files', muldirFiles, subfile(sub).name);        %the last step files of preprossing
    mulmotionfiles=spm_get('Files', rpdirFiles, rpfile(sub).name) ;   %head movement files for each run
    % develop individual onsets file
    % mulcondition_file_name=['sub',num2str(sub),'-',num2str(ses),'.mat'];  %behavior name for earch run of each subject
    
    mulcondpath=[WKD, '/fmriprepHCrun3/eventtimes'];    %behavior path
    mulcondition_file_name = eventFiles(sub).name;  
    mulconditions=load(fullfile(mulcondpath,mulcondition_file_name));  %load behavior results(.mat) for earch run of each subject
    jobs{1}.stats{1}.fmri_spec.sess(ses).scans=cellstr(mulfiles);    
    
    for j=1:ncon-8  
            jobs{1}.stats{1}.fmri_spec.sess(ses).cond(j).name = mulconditions.names{j};
            jobs{1}.stats{1}.fmri_spec.sess(ses).cond(j).onset = mulconditions.onsets{j};
            jobs{1}.stats{1}.fmri_spec.sess(ses).cond(j).duration = mulconditions.durations{j};
    end
    
    jobs{1}.stats{1}.fmri_spec.bases.hrf.derivs = [0 0];  % add time derivative = [1 0] no time derivative  = [0 0]; if changed, change the contrast file above and below
    jobs{1}.stats{1}.fmri_spec.sess(ses).multi_reg=cellstr(mulmotionfiles); % add head movements in model
     
    resultpath=[cwd, '/', [subnamecell{2}, '_', subnamecell{4}], '/'];
    jobs{1}.stats{2}.fmri_est.spmmat=cellstr(fullfile(resultpath,'SPM.mat')); % Save model settings and estimate results in SPM.mat
    save(fullfile(resultpath,'modelspecification.mat'),'jobs'); % save jobs
    
    spm_jobman('run', jobs);
    clear jobs;
    
    %%%%%%%%%%%%%%%%%%%seting contrast conditions%%%%%%%%%%%%%%%%%%%%%%%%%
    resultpath=[cwd, '/', [subnamecell{2}, '_', subnamecell{4}], '/'];
    jobs{1}.stats{1}.con.spmmat=cellstr(fullfile(resultpath,'SPM.mat'));  % write contrast information into SPM.mat
    % set contrast matrix for contrast conditions for individual subjexts
    % develop a pipline if you have too many contrasts
    jobs{1}.stats{1}.con.consess{1}.tcon=struct('name', 'load1Conload0','convec',...
        [0 -1 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0], 'sessrep', 'none');
    jobs{1}.stats{1}.con.consess{2}.tcon=struct('name', 'load2Conload0','convec',...
        [0 -1 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0], 'sessrep', 'none');
    jobs{1}.stats{1}.con.consess{3}.tcon=struct('name', 'load3Conload0','convec',...
        [0 -1 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0], 'sessrep', 'none');
    jobs{1}.stats{1}.con.consess{4}.tcon=struct('name', 'load4Conload0','convec',...
        [0 -1 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0], 'sessrep', 'none');
 
    spm_jobman('run',jobs);            % activation map is 'ps' format, which can be opened with GhostView and can also be converted to PDF format
    clear jobs;
end    % end of individual processing

%% select the contrast files and revised the file name; 
%%%%%%%%%%%%%% sz %%%%%%%%%%%%%%%%%
WKD = ['~/bidsfMRIprepSub/DATA_FS6R3load1n10RightTrailRun3'];
rawConP = [WKD, '/fmriprepSZrun3/F100sS6R3mmgunzipSZrun3Aun/1levelResults'];% revise
SUBFOLDER = dir(rawConP);
SUBFOLDER(1:2) = [];
% make new folder
foldername = {'2ndlevelAnalysisFS6R3load1n10', 'Run3CONTRAST'};
% 2ndlevelAnalysisS6R3load1n10;
for nu = 1
    mkdir([WKD, filesep], foldername{nu}); 
end
% /2ndlevelAnalysisS6R3load1n10/Run3CONTRAST;
for nu = 2: length(foldername) 
    mkdir([WKD, filesep, foldername{1}, filesep], foldername{nu});
end
% /2ndlevelAnalysisS6R3load1n10/Run3CONTRAST/load1Conload0;
ConFoldname = {'load1Conload0', 'load2Conload0', 'load3Conload0',...
    'load4Conload0'}; % use for new name
for nu1 = 1: length(ConFoldname)
    mkdir([WKD, filesep, foldername{1}, filesep, foldername{2}, filesep],...
        ConFoldname{nu1});
    mkdir([WKD, filesep, foldername{1}, filesep, foldername{2}, filesep,...
        ConFoldname{nu1}, filesep], 'SZ');
    mkdir([WKD, filesep, foldername{1}, filesep, foldername{2}, filesep,...
        ConFoldname{nu1}, filesep], 'HC');
end
rawname = {'con_0001.nii', 'con_0002.nii', 'con_0003.nii',...
        'con_0004.nii'};
for i =1: length(SUBFOLDER)
    display(['Computing the ', num2str(i), '-th subjects!']);
    clear rawConPP;
    subi = SUBFOLDER(i).name;% name
    rawConPP = [rawConP, filesep, SUBFOLDER(i).name];
    
    for j = 1: length(rawname)
        aimP = [WKD, filesep, foldername{1}, filesep, foldername{2}, filesep ,...
            ConFoldname{j}];
        clear newname;
        newname = ['1_', subi, '_', ConFoldname{j}, '.nii'];
        copyfile([rawConPP, filesep, rawname{j}], ...
            [aimP, '/SZ/', newname]);
    end
end
%%%%%%%%%%%%%% HC %%%%%%%%%%%%%%%%%
WKD = ['~/bidsfMRIprepSub/DATA_FS6R3load1n10RightTrailRun3'];
rawConP = [WKD, '/fmriprepHCrun3/F100sS6R3mmgunzipHCrun3Aun/1levelResults'];% revise
SUBFOLDER = dir(rawConP);
SUBFOLDER(1:2) = [];
% make new folder
foldername = {'2ndlevelAnalysisFS6R3load1n10', 'Run3CONTRAST'};
% /2ndlevelAnalysisS6R3load1n10/Run3CONTRAST/load1Conload0;
ConFoldname = {'load1Conload0', 'load2Conload0', 'load3Conload0',...
    'load4Conload0'}; % use for new name

rawname = {'con_0001.nii', 'con_0002.nii', 'con_0003.nii',...
        'con_0004.nii'};

for i =1: length(SUBFOLDER)
    display(['Computing the ', num2str(i), '-th subjects!']);
    clear rawConPP;
    subi = SUBFOLDER(i).name;% name
    rawConPP = [rawConP, filesep, SUBFOLDER(i).name];
    
    for j = 1: length(rawname)
        aimP = [WKD, filesep, foldername{1}, filesep, foldername{2}, filesep ,...
            ConFoldname{j}];
        clear newname;
        newname = ['2_', subi, '_', ConFoldname{j}, '.nii'];
        copyfile([rawConPP, filesep, rawname{j}], ...
            [aimP, '/HC/', newname]);
    end
end

