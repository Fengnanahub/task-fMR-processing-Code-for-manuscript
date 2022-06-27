%% compute the inverted U fitting of extracted contrast data ACC cluster;
load('sz30M4fitlmefdrsignACC_signal.mat');
load('M4fitlmefdrsignACC_signal.mat');
% SZ30
dataROICov = szSIGNALACC(:, 1:4);% load1-4;
[m, n] = size(dataROICov);
dataforModel =  [[zeros(m, 1); ones(m, 1); ones(m, 1)*2; ones(m, 1)*3],...
    repmat((1: m)', 4, 1), [dataROICov(:, 1);  dataROICov(:, 2);  dataROICov(:, 3);  dataROICov(:, 4)]];
dataforModelTa = array2table(dataforModel, 'VariableNames', {'loading',	'subject',...
    'activation'});
lme = fitlme(dataforModelTa,...
    'activation~loading^2+(1+loading^2|subject)'); % activation ~ 1 + loading + loading^2 + (1 | subject)
beta=  fixedEffects(lme);
[~,~,STAT] = randomEffects(lme);
 
% Organize the STAT(random effects) for each sub;
for i = 1: m
    randomCBA(1, i) = STAT{i*3-2, 4}; % intercept
    randomCBA(2, i) = STAT{i*3-1, 4}; % loading
    randomCBA(3, i) = STAT{i*3: i*3, 4};% loading^2
    fixRandCBA(1:3, i) = beta +randomCBA(1:3, i);
    % compute U measure
    applyM = @(t)(fixRandCBA(3, i)*(t.*t) +fixRandCBA(2, i)*t+fixRandCBA(1, i));
    axisX = - fixRandCBA(2, i)./(2*fixRandCBA(3, i));
    fixRandCBA(4, i)  = axisX; % axis
    fixRandCBA(5, i)  = applyM(axisX); %peakY
end
szfixRand = fixRandCBA';
szfixRandABCpeak = array2table(szfixRand, 'VariableNames', {'mixUCoeC',...
    'mixUCoeB', 'mixUCoeA', 'mixUCoeAxis', 'mixUCoepeakY'});
% %%%%%%%%%compute correlation U measure--ACC/RT
load('SZ30HC26hehavMatchImagTable.mat'); % 2ndlevelanalysisSZ30hc26/
load('ACCdataMatchImageSZ28hc21Ta.mat'); % 2ndlevelanalysisSZ30hc26/
indexsz = [1:10, 11:20, 22, 23, 25:30];  
for  i =4: 5
    for j = 11: 12
        % covariates age, sex, eduyear;
        ageCov = SZ30HC26hehavMatchImagTable{indexsz, 4};
        sexCov = SZ30HC26hehavMatchImagTable{indexsz, 5};
        eduCov = SZ30HC26hehavMatchImagTable{indexsz, 6};
        
        colu = (j-2)*2;
        [R, P] = partialcorr(szfixRandABCpeak{indexsz, i},...
            ACCdataMatchImageSZ28hc21Ta{1:28, j},...
            [ageCov, sexCov, eduCov], 'Type', 'Spearman');
        RP(i-3, colu-1) = R;
        RP(i-3, colu) = P;
    end
end
%% bootstrap %%% 
% compare the correlations between (1) sz28 axis with load4_ACC;
indexsz = [1:10, 11:20, 22, 23, 25:30]; 
axis4ACCagesexedu = [szfixRandABCpeak{indexsz, 4}, ACCdataMatchImageSZ28hc21Ta{1:28, 11},...
    SZ30HC26hehavMatchImagTable{indexsz, 4:6}];
rank  = [1: 1: 28]';
[bootstat, bootsam] = bootstrp(3000, @(x) [mean(x)], rank);
for i=1: 3000
    newrank = bootsam(:, i);
    for j = 1: 28
        newdata(j, 1) = axis4ACCagesexedu(newrank(j), 1); % axis
        newdata(j, 2) = axis4ACCagesexedu(newrank(j), 2); % load4_acc
        
        newdata(j, 3) = axis4ACCagesexedu(newrank(j), 3); % age
        newdata(j, 4) = axis4ACCagesexedu(newrank(j), 4); % sex1
        newdata(j, 5) = axis4ACCagesexedu(newrank(j), 5); % eduyear
    end
    bootstrapRPsz(i, :) = particalCorrela(newdata); % 3000*4
    disp(['computing the  ', num2str(i) , '-th bootstrap correlation !']);
end
bootstrapRPsz33Table = array2table(bootstrapRPsz, 'VariableNames',...
    {'axisload4ACCr', 'axisload4ACCp'});
Y = prctile(bootstrapRPsz(:, 1), [2.5 97.5], 'all'); % ==>>([-0.0908, -0.6633])
% %%%%%%%%%% HC26
dataROICov = hcSIGNALACC(:, 1:4);% load1-4;
[m, ~] = size(dataROICov);
dataforModel =  [[zeros(m, 1); ones(m, 1); ones(m, 1)*2; ones(m, 1)*3],...
    repmat((1: m)', 4, 1), [dataROICov(:, 1);  dataROICov(:, 2);  dataROICov(:, 3);  dataROICov(:, 4)]];
dataforModelTa = array2table(dataforModel, 'VariableNames', {'loading',	'subject',...
    'activation'});
lme = fitlme(dataforModelTa,...
    'activation~loading^2+(1+loading^2|subject)'); % activation ~ 1 + loading + loading^2 + (1 | subject)
beta=  fixedEffects(lme);
[~,~,STAT] = randomEffects(lme);
% Organize the STAT(random effects) for each sub;
for i = 1: m
    randomCBA(1, i) = STAT{i*3-2, 4}; % intercept
    randomCBA(2, i) = STAT{i*3-1, 4}; % loading
    randomCBA(3, i) = STAT{i*3: i*3, 4};% loading^2
    fixRandCBA(1:3, i) = beta +randomCBA(1:3, i);
    % compute U measure
    applyM = @(t)(fixRandCBA(3, i)*(t.*t) +fixRandCBA(2, i)*t+fixRandCBA(1, i));
    axisX = - fixRandCBA(2, i)./(2*fixRandCBA(3, i));
    fixRandCBA(4, i)  = axisX; % axis
    fixRandCBA(5, i)  = applyM(axisX); %peakY
end
hcfixRand = fixRandCBA';
hcfixRandABCpeak = array2table(hcfixRand, 'VariableNames', {'mixUCoeC',...
    'mixUCoeB', 'mixUCoeA', 'mixUCoeAxis', 'mixUCoepeakY'});

% %%%%%%%%%compute correlation U measure--ACC/RT hc26
load('SZ30HC26hehavMatchImagTable.mat'); % 2ndlevelanalysisSZ30hc26/
load('ACCdataMatchImageSZ28hc21Ta.mat'); % 2ndlevelanalysisSZ30hc26/
fixrandBAC56 = [szfixRandABCpeak; hcfixRandABCpeak];
indexhc = [31:35, 37:39, 42:46, 48:49, 51:56];  % remove? sz sub149(row10);
for  i =4: 5
    for j = 11: 12
        % covariates age, sex, eduyear;
        ageCov = SZ30HC26hehavMatchImagTable{indexhc, 4};
        sexCov = SZ30HC26hehavMatchImagTable{indexhc, 5};
        eduCov = SZ30HC26hehavMatchImagTable{indexhc, 6};
        % load0ACC = ACCdataMatchImageSZ28hc21Ta{index, 3};
        colu = (j-2)*2;
        [R, P] = partialcorr(fixrandBAC56{indexhc, i},...
            ACCdataMatchImageSZ28hc21Ta{29:49, j},...
            [ageCov, sexCov, eduCov], 'Type', 'Spearman');
        RP(i-3, colu-1) = R;
        RP(i-3, colu) = P;
    end
end

%% between-group anovan analysis
anova1([hcfixRandABCpeak{:, 3}; szfixRandABCpeak{:, 3}],...
    [zeros(1,26), ones(1,30)])% A
anova1([hcfixRandABCpeak{:, 2}; szfixRandABCpeak{:, 2}],...
    [zeros(1,26), ones(1,30)])% B
anova1([hcfixRandABCpeak{:, 1}; szfixRandABCpeak{:, 1}],...
    [zeros(1,26), ones(1,30)])% C
anova1([hcfixRandABCpeak{:, 4}; szfixRandABCpeak{:, 4}],...
    [zeros(1,26), ones(1,30)])% axis
anova1([hcfixRandABCpeak{:, 5}; szfixRandABCpeak{:, 5}],...
    [zeros(1,26), ones(1,30)])% peak Y
%%anova WITH covariates--
for i = 1: 5
    data2comp = [szfixRandABCpeak{:, i};hcfixRandABCpeak{:, i}];
    grouplabel = [ones(30, 1);zeros(26, 1)];
    % age,sex,edu
    groups = [grouplabel,SZ30HC26hehavMatchImagTable{:, 3},...
        SZ30HC26hehavMatchImagTable{:,4},SZ30HC26hehavMatchImagTable{:, 5}];
    p= anovan(data2comp, groups, 'continuous', [2, 4]);
    pf(i, 1) = p(1);
end
%% %%%%%%%%% axis/peakY - TMTA-RT/TMTB-RT
% SZ30 -20220621add
ageCov = SZ30HC26hehavMatchImagTable{1:30, 3};
sexCov = SZ30HC26hehavMatchImagTable{1:30, 4};
eduCov = SZ30HC26hehavMatchImagTable{1:30, 5};
col = [7, 8];% col index of TMTA-RT, TMTB-RT;
for  i =4: 5 % Cluster peak Signal
    for j = 1:length(col)
        colu = j*2;
        [R, P] = partialcorr(szfixRandABCpeak{1:30, i},...
            SZ30HC26hehavMatchImagTable{1:30, col(j)},...
            [ageCov, sexCov, eduCov], 'Type', 'Spearman', 'Rows', 'Complete');
        RPSZ(i, colu-1) = R;
        RPSZ(i, colu) = P;
    end
end
% HC21
ageCov = SZ30HC26hehavMatchImagTable{31:56, 3};
sexCov = SZ30HC26hehavMatchImagTable{31:56, 4};
eduCov = SZ30HC26hehavMatchImagTable{31:56, 5};
col = [7, 8];% col index of TMTA-RT, TMTB-RT;
for  i =4: 5 % axis/peak
    for j = 1:length(col)
        colu = j*2;
        [R, P] = partialcorr(hcfixRandABCpeak{:, i},...
            SZ30HC26hehavMatchImagTable{31:56, col(j)},...
            [ageCov, sexCov, eduCov], 'Type', 'Spearman', 'Rows', 'Complete');
        RPhc(i, colu-1) = R;
        RPhc(i, colu) = P;
    end
end
%% %%%%%%%%%%%%%%%%%%%%%%%OTHER correlation
%% extract the laod4-load0 one and whole Ttest fdr cluster signals
% compute the mean signal from each cluster area  each sub;
%%%%%%%%%%%%%%% sz30hc26 load4-load0 SZ OneT
WKD = ['~DATA_FS6R3load1n10RightTrailRun3/2ndlevelAnalysisFS6R3load1n10',...
    '/2ndlevelanalysisSZ30hc26'];
cellname = {'SZ', 'HC'};
cellclumask =  {'cluster1_mask.nii', 'cluster2_mask.nii',...
    'cluster3_mask.nii', 'cluster4_mask.nii', 'cluster5_mask.nii',...
    'cluster6_mask.nii', 'cluster7_mask.nii'};
% clusrer 1 -7
for j= 1: length(cellname)
    clear subfiles;
    clear signal;
    subfileW = [WKD, '/load4Conload0/', cellname{j}];
     subfiles = dir([subfileW, '/*.nii']);
    for i = 1: length(cellclumask)
        maskvol = spm_vol([WKD, '/CorrelationAnaClustBehav/load4load0SZoneT/',...
            cellclumask{i}]);
        maskAAL = spm_read_vols(maskvol);
        
        for subi = 1: length(subfiles)
            V = spm_vol([subfileW, filesep, subfiles(subi).name]);
            [con1, XYZ] = spm_read_vols(V);  %  spm read image
            % con1--67*77*65double; XYZ--3*325325double;
            con1_data = con1(maskAAL>0); %  this is a vector
            signal(subi, i) = nanmean(con1_data); %-->>30sub * 7mask
        end
    end
    szhcSignal{j, 1} = signal; 
end

% %% compute the correlations (cluster1-7 signals with ACC;) %%
szhcSignal56 = [szhcSignal{1, 1}; szhcSignal{2, 1}];
indexsz = [1:10, 11:20, 22, 23, 25:30];
indexhc = [31:35, 37:39, 42:46, 48:49, 51:56];
index = [indexsz, indexhc];
load('SZ30HC26hehavMatchImagTable.mat');
load('ACCdataMatchImageSZ28hc21Ta.mat');
% SZ
ageCov = SZ30HC26hehavMatchImagTable{indexsz, 4};
sexCov = SZ30HC26hehavMatchImagTable{indexsz, 5};
eduCov = SZ30HC26hehavMatchImagTable{indexsz, 6};
for  i =1: 7 % clusterSignal
    for j = 11:12 % load4-ACC/RT
        colu = (j-2)*2;
        [R, P] = partialcorr(szhcSignal56(indexsz, i),...
            ACCdataMatchImageSZ28hc21Ta{1:28, j},...
            [ageCov, sexCov, eduCov], 'Type', 'Spearman');
        RPSZ(i, colu-1) = R;
        RPSZ(i, colu) = P;
    end
end
% HC
ageCov = SZ30HC26hehavMatchImagTable{indexhc, 4};
sexCov = SZ30HC26hehavMatchImagTable{indexhc, 5};
eduCov = SZ30HC26hehavMatchImagTable{indexhc, 6};
for  i =1:7 % clusterSignal
    for j = 11:12 % load4-ACC/RT
        colu = (j-2)*2;
        [R, P] = partialcorr(szhcSignal56(indexhc, i),...
            ACCdataMatchImageSZ28hc21Ta{29:49, j},...
            [ageCov, sexCov, eduCov], 'Type', 'Spearman');
        RPSZ(i, colu-1) = R;
        RPSZ(i, colu) = P;
    end
end

% %%%%%%%%% PANSS-P, N, G %%%%%%%
ageCov = SZ30HC26hehavMatchImagTable{1:30, 4};
sexCov = SZ30HC26hehavMatchImagTable{1:30, 5};
eduCov = SZ30HC26hehavMatchImagTable{1:30, 6};
for  i =1: 7 % cluster
    for j = 68: 70 % panss-P,N,G
        colu = (j-66)*2;
        [R, P] = partialcorr(szhcSignal56(1:30, i),...
            SZ30HC26hehavMatchImagTable{1:30, j},...
            [ageCov, sexCov, eduCov], 'Type', 'Spearman', 'Rows', 'complete');
        RPszpanss(i, colu-1) = R;
        RPszpanss(i, colu) = P;
    end
end

% %%%%%%%%%%%%%%%%%%partial correlation --TMT-RT;
% SZ30 -20220620add
ageCov = SZ30HC26hehavMatchImagTable{1:30, 4};
sexCov = SZ30HC26hehavMatchImagTable{1:30, 5};
eduCov = SZ30HC26hehavMatchImagTable{1:30, 6};
col = [21, 23];% col index of TMTA-RT, TMTB-RT;
for  i =1: 7 % clusterSignal
    for j = 1:length(col)
        colu = j*2;
        [R, P] = partialcorr(szhcSignal56(1:30, i),...
            SZ30HC26hehavMatchImagTable{1:30, col(j)},...
            [ageCov, sexCov, eduCov], 'Type', 'Spearman', 'Rows', 'Complete');
        RPSZ(i, colu-1) = R;
        RPSZ(i, colu) = P;
    end
end
% HC26 -20220620
ageCov = SZ30HC26hehavMatchImagTable{31:56, 4};
sexCov = SZ30HC26hehavMatchImagTable{31:56, 5};
eduCov = SZ30HC26hehavMatchImagTable{31:56, 6};
col = [21, 23];% col index of TMTA-RT, TMTB-RT;
for  i =1: 7 % clusterSignal
    for j = 1:length(col)
        colu = j*2;
        [R, P] = partialcorr(szhcSignal56(31:56, i),...
            SZ30HC26hehavMatchImagTable{31:56, col(j)},...
            [ageCov, sexCov, eduCov], 'Type', 'Spearman', 'Rows', 'Complete');
        RPhc(i, colu-1) = R;
        RPhc(i, colu) = P;
    end
end
%% extract the laod4-load0 whole Ttest fdr cluster signals
% compute the mean signal from each cluster area  each sub;
%%%%%%%%%%%%%%% sz30hc26 load4-load0 SZHC whole group T
WKD = ['~/DATA_FS6R3load1n10RightTrailRun3/2ndlevelAnalysisFS6R3load1n10',...
    '/2ndlevelanalysisSZ30hc26'];
cellname = {'SZ', 'HC'};
for cl = 1: 19
    cellclumask{1, cl} = ['cluster', num2str(cl), '_mask.nii'];
end
% clusrer 1 -19
for j= 1: length(cellname)
    clear subfiles;
    clear signal;
    subfileW = [WKD, '/load4Conload0/', cellname{j}];
     subfiles = dir([subfileW, '/*.nii']);
    for i = 1: length(cellclumask)
        maskvol = spm_vol([WKD, '/CorrelationAnaClustBehav/load4load0SZHCwholeT/',...
            cellclumask{i}]);
        maskAAL = spm_read_vols(maskvol);
        
        for subi = 1: length(subfiles)
            V = spm_vol([subfileW, filesep, subfiles(subi).name]);
            [con1, XYZ] = spm_read_vols(V);  %  spm read image
            % con1--67*77*65double; XYZ--3*325325double;
            con1_data = con1(maskAAL>0); %  this is a vector
            signal(subi, i) = nanmean(con1_data); %-->>sub * 19mask
        end
    end
    szhcSignal{j, 1} = signal; 
end

% %%%%%%%%%%%%%% partial correlation --TMT-RT %%%%%%%%%%%
% SZ30 -20220620add
szhcSignal56 = [szhcSignal{1, 1}; szhcSignal{2, 1}];
ageCov = SZ30HC26hehavMatchImagTable{1:30, 4};
sexCov = SZ30HC26hehavMatchImagTable{1:30, 5};
eduCov = SZ30HC26hehavMatchImagTable{1:30, 6};
col = [21, 23];% col index of TMTA-RT, TMTB-RT;
for  i =1: 19 % clusterSignal
    for j = 1:length(col)
        colu = j*2;
        [R, P] = partialcorr(szhcSignal56(1:30, i),...
            SZ30HC26hehavMatchImagTable{1:30, col(j)},...
            [ageCov, sexCov, eduCov], 'Type', 'Spearman', 'Rows', 'Complete');
        RPSZ(i, colu-1) = R;
        RPSZ(i, colu) = P;
    end
end
% HC26 -20220620add
ageCov = SZ30HC26hehavMatchImagTable{31:56, 4};
sexCov = SZ30HC26hehavMatchImagTable{31:56, 5};
eduCov = SZ30HC26hehavMatchImagTable{31:56, 6};
col = [21, 23];% col index of TMTA-RT, TMTB-RT;
for  i =1: 19 % clusterSignal
    for j = 1:length(col)
        colu = j*2;
        [R, P] = partialcorr(szhcSignal56(31:56, i),...
            SZ30HC26hehavMatchImagTable{31:56, col(j)},...
            [ageCov, sexCov, eduCov], 'Type', 'Spearman', 'Rows', 'Complete');
        RPhc(i, colu-1) = R;
        RPhc(i, colu) = P;
    end
end

%%%%%%%%%%% correlation  signal-ACC/RT %%%%%%%%%%%%%%%%%
szhcSignal56 = [szhcSignal{1, 1}; szhcSignal{2, 1}];
indexsz = [1:10, 11:20, 22, 23, 25:30];  
indexhc = [31:35, 37:39, 42:46, 48:49, 51:56];
index = [indexsz, indexhc];
load('SZ30HC26hehavMatchImagTable.mat');
load('ACCdataMatchImageSZ28hc21Ta.mat');
% SZ
ageCov = SZ30HC26hehavMatchImagTable{indexsz, 4};
sexCov = SZ30HC26hehavMatchImagTable{indexsz, 5};
eduCov = SZ30HC26hehavMatchImagTable{indexsz, 6};
for  i =1: 19 % clusterSignal
    for j = 11: 12 % ACC
        colu = (j-2)*2;
        [R, P] = partialcorr(szhcSignal56(indexsz, i),...
            ACCdataMatchImageSZ28hc21Ta{1:28, j},...
            [ageCov, sexCov, eduCov], 'Type', 'Spearman');
        RPSZ(i, colu-1) = R;
        RPSZ(i, colu) = P;
    end
end
% HC
ageCov = SZ30HC26hehavMatchImagTable{indexhc, 4};
sexCov = SZ30HC26hehavMatchImagTable{indexhc, 5};
eduCov = SZ30HC26hehavMatchImagTable{indexhc, 6};
for  i =1:19 % clusterSignal
    for j = 11:12 % Load4-ACC/RT
        colu = (j-2)*2;
        [R, P] = partialcorr(szhcSignal56(indexhc, i),...
            ACCdataMatchImageSZ28hc21Ta{29:49, j},...
            [ageCov, sexCov, eduCov], 'Type', 'Spearman');
        RPSZ(i, colu-1) = R;
        RPSZ(i, colu) = P;
    end
end
% sz28+HC21
ageCov = SZ30HC26hehavMatchImagTable{index, 4};
sexCov = SZ30HC26hehavMatchImagTable{index, 5};
eduCov = SZ30HC26hehavMatchImagTable{index, 6};
group = [ones(28,1); zeros(21, 1)];
for  i =1:19 % clusterSignal
    for j = 3:13 % ACC
        colu = (j-2)*2;
        [R, P] = partialcorr(szhcSignal56(index, i),...
            ACCdataMatchImageSZ28hc21Ta{:, j},...
            [ageCov, sexCov, eduCov, group], 'Type', 'Spearman');
        RPSZ(i, colu-1) = R;
        RPSZ(i, colu) = P;
    end
end

%%%%%%%%%%% correlation  signal-PANSS %%%%%%%%%%%%%%%%%
ageCov = SZ30HC26hehavMatchImagTable{1:30, 3};
sexCov = SZ30HC26hehavMatchImagTable{1:30, 4};
eduCov = SZ30HC26hehavMatchImagTable{1:30, 5};
for  i =1: 19 % cluster
    for j = 10: 12 % panss
        colu = (j-9)*2;
        [R, P] = partialcorr(szhcSignal56(1:30, i),...
            SZ30HC26hehavMatchImagTable{1:30, j},...
            [ageCov, sexCov, eduCov], 'Type', 'Spearman', 'Rows', 'complete');
        RPszpanss(i, colu-1) = R;
        RPszpanss(i, colu) = P;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% M2 SZ loading FDR cluster--correlation;%%%%%%%%%%%%%
%% fitting the fitlme model for each location for the randomeffects Cofi of loading;
[m, n] = size(volumes); %30*4
for k = 1: length(volumes{1,1}) % for each location--47634
    display(['Computing the ', num2str(k), '-th location point!']);
    clear subidata;
    for colu = 1: n
        for row = 1: m
            subidata(row, colu) = volumes{row, colu}(k);
        end
    end
    LOCdata{k, 1} = subidata;
    % fitlme model 1
    subidata(isnan(subidata) == 1)=[];
    Coeffic(k, 1) = 0;
    fixedCoefficSave{k, 1} = [];
    randomCoefficSave{k, 1} = [];
    if rank(subidata)==4
        clear lme;
        dataROICov = subidata(:, 1:4);% load1-4
        dataforModel =  [[zeros(30, 1); ones(30, 1); ones(30, 1)*2; ones(30, 1)*3],...
            repmat((1: 30)', 4, 1),...
            [dataROICov(:, 1);  dataROICov(:, 2);  dataROICov(:, 3);  dataROICov(:, 4)]];
        dataforModelTa = array2table(dataforModel, 'VariableNames', {'loading',	'subject',...
            'activation'});
        lme = fitlme(dataforModelTa,...
            'activation~loading+(1+loading|subject)');
        Coeffic(k, 1) = lme.Coefficients.Estimate(2); % loading Coeffic;
         % save the random effects Coeffic titleddataset of each location
        fixedCoefficSave{k, 1}=  fixedEffects(lme); % fixed coeffic;
        [~, ~,stats] = randomEffects(lme);
        randomCoefficSave{k, 1} = stats; % random coeffic;
    end
end

FixRandomLoading = zeros(m, length(volumes{1,1})); % 2*30sub * 47634location;
for loc = 1: length(volumes{1,1})
    if ~isempty(randomCoefficSave{loc, 1})
        allrandom1(:, loc) = randomCoefficSave{loc, 1}.Estimate;
        for ro = 1: m
            % random +fixed for each location each sub;30*47634
            FixRandomLoading(ro, loc) = allrandom1(2*ro, loc)+Coeffic(loc, 1);
        end
    end
end
% %%write back the fixed+random loading Coffic  value to image
clear Fcon;
for i = 1: m
    V.fname  = ['model2fitlmeSZ_', num2str(i), '_FixRandomCoffic.nii'];
    
    Fcon = con1;
    Fcon(maskAAL<1) = 0;
    Fcon(maskAAL>0) = FixRandomLoading(i, :) ;
    spm_write_vol(V, Fcon);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% extract the signal of cluster1-16 Coffic value(fix+random) %%
WKD = ['/share/inspurStorage/home1/Fengnana/DATA210725/SZ2010data',...
    '/2020-SCZ/DATAfMRIprepAGAIN2112/bidsfMRIprepSub',...
    '/DATA_FS6R3load1n10RightTrailRun3/2ndlevelAnalysisFS6R3load1n10',...
    '/2ndlevelanalysisSZ30hc26'];
for cl = 1: 16
    cellclumask{1, cl} = ['cluster', num2str(cl), '_mask.nii'];
end
% clusrer 1 -16
subfileW = [WKD, '/SZ30load1234Ftest/M2loadingCofficfixPrandomImage'];
subfiles = dir([subfileW, '/*.nii']);
for i = 1: length(cellclumask)
    maskvol = spm_vol([WKD, '/CorrelationAnaClustBehav/M2SZloadingFDRcluster/',...
        cellclumask{i}]);
    maskAAL = spm_read_vols(maskvol);
    
    for subi = 1: length(subfiles)
        V = spm_vol([subfileW, filesep, subfiles(subi).name]);
        [con1, XYZ] = spm_read_vols(V);  %  spm read image
        % con1--67*77*65double; XYZ--3*325325double;
        con1_data = con1(maskAAL>0); %  this is a vector
        signal(subi, i) = nanmean(con1_data); %-->>30sub * 16mask
    end
end
M2szCoefSignal = signal;
%%%%%%%%%%%%%% parial correlation -TMT %%%%%%%%%%%
% SZ30 -20220620add
ageCov = SZ30HC26hehavMatchImagTable{1:30, 4};
sexCov = SZ30HC26hehavMatchImagTable{1:30, 5};
eduCov = SZ30HC26hehavMatchImagTable{1:30, 6};
col = [21, 23];% col index of TMTA-RT, TMTB-RT;
for  i =1: 16 % clusterSignal
    for j = 1:length(col)
        colu = j*2;
        [R, P] = partialcorr(M2szCoefSignal(1:30, i),...
            SZ30HC26hehavMatchImagTable{1:30, col(j)},...
            [ageCov, sexCov, eduCov], 'Type', 'Spearman', 'Rows', 'Complete');
        RPSZ(i, colu-1) = R;
        RPSZ(i, colu) = P;
    end
end

% %%%%compute the correlations (cluster1-16 Coffic value with ACC/RT)
indexsz = [1:10, 11:20, 22, 23, 25:30];  
indexhc = [31:35, 37:39, 42:46, 48:49, 51:56];
index = [indexsz, indexhc];
load('SZ30HC26hehavMatchImagTable.mat');
load('ACCdataMatchImageSZ28hc21Ta.mat');
% SZ
ageCov = SZ30HC26hehavMatchImagTable{indexsz, 4};
sexCov = SZ30HC26hehavMatchImagTable{indexsz, 5};
eduCov = SZ30HC26hehavMatchImagTable{indexsz, 6};
for  i =1: 16 % clusterSignal;
    for j = 11: 12 % load4 ACC;
        colu = (j-2)*2;
        [R, P] = partialcorr(M2szCoefSignal(indexsz, i),...
            ACCdataMatchImageSZ28hc21Ta{1:28, j},...
            [ageCov, sexCov, eduCov], 'Type', 'Spearman');
        RPSZ(i, colu-1) = R;
        RPSZ(i, colu) = P;
    end
end

% %%%%%%%%% correlation -PANSS %%%%%%%
ageCov = SZ30HC26hehavMatchImagTable{1:30, 4};
sexCov = SZ30HC26hehavMatchImagTable{1:30, 5};
eduCov = SZ30HC26hehavMatchImagTable{1:30, 6};
for  i =1: 16 % cluster
    for j = 68: 70 % panss-P,N,G
        colu = (j-66)*2;
        [R, P] = partialcorr(M2szCoefSignal(1:30, i),...
            SZ30HC26hehavMatchImagTable{1:30, j},...
            [ageCov, sexCov, eduCov], 'Type', 'Spearman', 'Rows', 'complete');
        RPszpanss(i, colu-1) = R;
        RPszpanss(i, colu) = P;
    end
end

%% attention-related analysis (TMTA-WM/TMTB)
load('SZ30HC26hehavMatchImagTable.mat'); % 2ndlevelanalysisSZ30hc26/
load('ACCdataMatchImageSZ28hc21Ta.mat'); % 2ndlevelanalysisSZ30hc26/
indexsz = [1:10, 11:20, 22, 23, 25:30]; 
indexhc = [31:35, 37:39, 42:46, 48:49, 51:56];

% compute mean and SD value;
colu = [21,22,23,24];
for i = 1:length(colu)
    SZmeansd(i,1) = nanmean(SZ30HC26hehavMatchImagTable{1:30, colu(i)});
    SZmeansd(i,2) = nanstd(SZ30HC26hehavMatchImagTable{1:30, colu(i)});
    HCmeansd(i,1) = nanmean(SZ30HC26hehavMatchImagTable{31:56, colu(i)});
    HCmeansd(i,2) = nanstd(SZ30HC26hehavMatchImagTable{31:56, colu(i)});
end

anova1([SZ30HC26hehavMatchImagTable{31:56, 23};...
    SZ30HC26hehavMatchImagTable{1:30, 23}],...
    [zeros(1,26), ones(1,30)])% TMTA/B-RT
anova1([SZ30HC26hehavMatchImagTable{31:56, 24};...
    SZ30HC26hehavMatchImagTable{1:30, 24}],...
    [zeros(1,26), ones(1,30)])% TMTA/B-WN
%%anova WITH covariates--
for i = 23:24
    data2comp = [SZ30HC26hehavMatchImagTable{1:30, i};...
        SZ30HC26hehavMatchImagTable{31:56, i}];
    grouplabel = [ones(30, 1);zeros(26, 1)];
    groups = [grouplabel,SZ30HC26hehavMatchImagTable{:, 4},...
        SZ30HC26hehavMatchImagTable{:, 5},SZ30HC26hehavMatchImagTable{:, 6}];
    p= anovan(data2comp, groups, 'continuous', [2, 4]);
    pf(i, 1) = p(1);
end

% correlation-SZ
for  i =23:24
    for j = 11: 12
        % covariates age, sex, eduyear;
        ageCov = SZ30HC26hehavMatchImagTable{indexsz, 4};
        sexCov = SZ30HC26hehavMatchImagTable{indexsz, 5};
        eduCov = SZ30HC26hehavMatchImagTable{indexsz, 6};
        
        colu = (j-2)*2;
        [R, P] = partialcorr(SZ30HC26hehavMatchImagTable{indexsz, i},...
            ACCdataMatchImageSZ28hc21Ta{1:28, j},...
            [ageCov, sexCov, eduCov],'Type', 'Spearman','Rows','Complete');
        RP(i-20, colu-1) = R;
        RP(i-20, colu) = P;
    end
end

% HC
indexhc = [31:35, 37:39, 42:46, 48:49, 51:56];
for  i =21
    for j = 3: 13
        % covariates age, sex, eduyear;
        ageCov = SZ30HC26hehavMatchImagTable{indexhc, 4};
        sexCov = SZ30HC26hehavMatchImagTable{indexhc, 5};
        eduCov = SZ30HC26hehavMatchImagTable{indexhc, 6};
        
        colu = (j-2)*2;
        [R, P] = partialcorr(SZ30HC26hehavMatchImagTable{indexhc, i},...
            ACCdataMatchImageSZ28hc21Ta{29:49, j},...
            [ageCov, sexCov, eduCov], 'Type', 'Spearman','Rows','Complete');
        RP(i-20, colu-1) = R;
        RP(i-20, colu) = P;
    end
end
%%scatter
indexhc = [31:35, 37:39, 42:46, 48:49, 51:56];
figure
X = [ones(21, 1) SZ30HC26hehavMatchImagTable{indexhc, 3:5}];
[b,bint,r] = regress(SZ30HC26hehavMatchImagTable{indexhc, 7}, X); % TMTA-RT
[b1,bint1,r1] = regress(ACCdataMatchImageSZ28hc21Ta{29:49, 4}, X);% Load4-ACC

theta = glmfit(r, r1);
yCalc1 = theta(2)*r +theta(1);
plot(r(:, 1), r1(:, 1), 'Color', sixteen2ten('#B6825F')/255,...
    'Marker','o', 'Markersize',10, 'LineStyle', 'none', 'LineWidth', 1);
hold on;

plot(r(:, 1), yCalc1, 'Color', sixteen2ten('#FACB66')/255, 'LineWidth', 2);
ylim([-8, 5]);
yticks([-8 0 5]);
xlim([-15, 30]);
xticks([-15 0 30]);
% xticklabels({'-4', '0', '4'});
ylabel('Residual accuracy of load4', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('Residual TMT-RT', 'FontSize', 12, 'FontWeight', 'bold');
% title(['Correlation of slope with load', num2str(i), ' ACC']);
grid on