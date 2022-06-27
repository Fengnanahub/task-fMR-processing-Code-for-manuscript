%% mixed effect model fitlme fitting ,sz30hc26 matching AGE
%% 1-HC: x ~ loadings + age + sex + edu + meanFD + (1 | subID )
% HC26
WKD = ['~/DATA_FS6R3load1n10RightTrailRun3',...
    '/2ndlevelAnalysisFS6R3load1n10/2ndlevelanalysisSZ30hc26/HC26load1234Ftest'];

cellload = {'load1Conload0HC', 'load2Conload0HC', 'load3Conload0HC', 'load4Conload0HC'};
% read the image data Cell{sub26*load1234}
maskvol = spm_vol([WKD, '/aal3_3mmBinnocerellu.nii']);
maskAAL = spm_read_vols(maskvol); % AAL Mask
for i = 1:length(cellload)
    clear niifiles;
    loadpath = [WKD, filesep, cellload{i}];
    niifiles = dir([loadpath, '/*.nii']);
    for j = 1:length(niifiles)
        clear V;
        clear volume;
        V = spm_vol([loadpath, filesep, niifiles(j).name]);
        [con1, XYZ] = spm_read_vols(V);  %  spm read image
        % con1--67*77*65double; XYZ--3*325325double;
        
        con1_data = con1(maskAAL>0); % get the voxles within the AAL mask only, this is a vector
        volumes{j, i} = con1_data; % Cell{sub26*load1234}--47634*1double
    end
end
[m, n] = size(volumes); %26*4
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
    Tvalue(k, 1) = 0;
    pvalueD(k, 1) = 1;
    if rank(subidata)==4
        clear lme;
        dataROICov = [subidata(:, 1:4), SZ30HC26hehavMatchImagTable{31:56, 4},...
            SZ30HC26hehavMatchImagTable{31:56, 5}, SZ30HC26hehavMatchImagTable{31:56, 6},...
            SZ30HC26hehavMatchImagTable{31:56, 9}];% load1-4,age, sex,edu, fd;
        dataforModel =  [[zeros(26, 1); ones(26, 1); ones(26, 1)*2; ones(26, 1)*3],...
            repmat((1: 26)', 4, 1), [dataROICov(:, 1);  dataROICov(:, 2);  dataROICov(:, 3);  dataROICov(:, 4)],...
            repmat(dataROICov(:, 5), 4, 1), repmat(dataROICov(:, 6), 4, 1), ...
            repmat(dataROICov(:, 7), 4, 1), repmat(dataROICov(:, 8), 4, 1)];
        dataforModelTa = array2table(dataforModel, 'VariableNames', {'loading',	'subject',...
            'activation', 'ageCov', 'sexCov', 'eduCov', 'fdCov'});
        lme = fitlme(dataforModelTa,...
            'activation~loading+ageCov+sexCov+eduCov+fdCov+(1|subject) ');
         Tvalue(k, 1) = lme.Coefficients.tStat(2); % 47634*1cell
         pvalueD(k, 1) = lme.Coefficients.pValue(2); % 47634*1double
    end
end
fdr1 = mafdr(pvalueD, 'BHFDR', true);
sign1 = sum(fdr1<0.05);% no sign

%% 2-SCZ:   x ~ loadings + age + sex + edu + meanFD + (1 | subjID )
% sz30 20:08-23:10
WKD = ['~/DATA_FS6R3load1n10RightTrailRun3',...
    '/2ndlevelAnalysisFS6R3load1n10/2ndlevelanalysisSZ30hc26/SZ30load1234Ftest'];
cellload = {'load1Conload0SZ', 'load2Conload0SZ', 'load3Conload0SZ', 'load4Conload0SZ'};
% read the image data Cell{sub30*load1234}
maskvol = spm_vol([WKD, '/aal3_3mmBinnocerellu.nii']);
maskAAL = spm_read_vols(maskvol); % AAL Mask
for i = 1:length(cellload)
    
    clear niifiles;
    loadpath = [WKD, filesep, cellload{i}];
    niifiles = dir([loadpath, '/*.nii']);
    for j = 1:length(niifiles)
        clear V;
        clear volume;
        V = spm_vol([loadpath, filesep, niifiles(j).name]);
        [con1, XYZ] = spm_read_vols(V);  %  spm read image
        % con1--67*77*65double; XYZ--3*325325double;
        
        con1_data = con1(maskAAL>0); % get the voxles within the AAL mask only, this is a vector
        volumes{j, i} = con1_data; % Cell{sub30*load1234}--47634*1double
    end
end
% fitting the fitlme model for each location;
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
    Tvalue(k, 1) = 0;
    pvalueD(k, 1) = 1;
    if rank(subidata)==4
        clear lme;
        dataROICov = [subidata(:, 1:4), SZ30HC26hehavMatchImagTable{1:30, 4},...
            SZ30HC26hehavMatchImagTable{1:30, 5}, SZ30HC26hehavMatchImagTable{1:30, 6},...
            SZ30HC26hehavMatchImagTable{1:30, 9}];% load1-4,age, sex,edu, fd;
        dataforModel =  [[zeros(30, 1); ones(30, 1); ones(30, 1)*2; ones(30, 1)*3],...
            repmat((1: 30)', 4, 1), [dataROICov(:, 1);  dataROICov(:, 2);  dataROICov(:, 3);  dataROICov(:, 4)],...
            repmat(dataROICov(:, 5), 4, 1), repmat(dataROICov(:, 6), 4, 1), ...
            repmat(dataROICov(:, 7), 4, 1), repmat(dataROICov(:, 8), 4, 1)];
        dataforModelTa = array2table(dataforModel, 'VariableNames', {'loading',	'subject',...
            'activation', 'ageCov', 'sexCov', 'eduCov', 'fdCov'});
        lme = fitlme(dataforModelTa,...
            'activation~loading+ageCov+sexCov+eduCov+fdCov+(1|subject) ');
         Tvalue(k, 1) = lme.Coefficients.tStat(2); % loadings
         pvalueD(k, 1) = lme.Coefficients.pValue(2); % 47634*1double
    end
end
fdr1 = mafdr(pvalueD, 'BHFDR', true);
sign1 = sum(fdr1<0.05); %1841
figure;
hist(Tvalue(fdr1<0.05));

%% write the data back to the brain (1-pvalue or 1-fdr pvalue)
signT = Tvalue(fdr1<0.05);
va = min(signT(signT>0));
vaN = max(signT(signT<0));% -3.1760-3.1742
% Tvalue-
V.fname = 'model2fitlmeTvalueSign.nii';
Fcon = con1;
siTvalue = Tvalue;
siTvalue(fdr1>=0.05) = 0;
Fcon(maskAAL<1)= 0;
Fcon(maskAAL>0)= siTvalue;
spm_write_vol(V, Fcon);
% extract sign (clustersize 10)mask;
% p-value
V.fname = 'model2fitlmePvalueSign.nii';
Fcon = con1;
siPvalue = pvalueD;
siPvalue(fdr1>=0.05) = 0;
Fcon(maskAAL<1)= 0;
Fcon(maskAAL>0)= siPvalue;
spm_write_vol(V, Fcon);

%% 3-HC vs SCZ: x ~ group : loadings:group + loadings + age + sex + edu + meanFD + (1 | subjID )
%11:17--18:10
WKD = ['~/DATA_FS6R3load1n10RightTrailRun3',...
    '/2ndlevelAnalysisFS6R3load1n10/2ndlevelanalysisSZ30hc26'];
WKD1 = [WKD, '/SZ30load1234Ftest'];
cellload1 = {'load1Conload0SZ', 'load2Conload0SZ', 'load3Conload0SZ', 'load4Conload0SZ'};
% read the image data Cell{sub30*load1234}
maskvol = spm_vol([WKD, '/aal3_3mmBinnocerellu.nii']);
maskAAL = spm_read_vols(maskvol); % AAL Mask
for i = 1:length(cellload1)
    clear niifiles;
    loadpath = [WKD1, filesep, cellload1{i}];
    niifiles = dir([loadpath, '/*.nii']);
    for j = 1:length(niifiles)
        clear V;
        clear volume;
        V = spm_vol([loadpath, filesep, niifiles(j).name]);
        [con1, XYZ] = spm_read_vols(V);  %  spm read image
        % con1--67*77*65double; XYZ--3*325325double;
        
        con1_data = con1(maskAAL>0); % get the voxles within the AAL mask only, this is a vector
        volumessz{j, i} = con1_data; % Cell{sub30*load1234}--47634*1double
    end
end
WKD2 = [WKD, '/HC26load1234Ftest'];
cellload2 = {'load1Conload0HC', 'load2Conload0HC', 'load3Conload0HC', 'load4Conload0HC'};
% read the image data Cell{sub26*load1234}
maskvol = spm_vol([WKD, '/aal3_3mmBinnocerellu.nii']);
maskAAL = spm_read_vols(maskvol); % AAL Mask
for i = 1:length(cellload2)
    clear niifiles;
    loadpath = [WKD2, filesep, cellload2{i}];
    niifiles = dir([loadpath, '/*.nii']);
    for j = 1:length(niifiles)
        clear V;
        clear volume;
        V = spm_vol([loadpath, filesep, niifiles(j).name]);
        [con2, XYZ] = spm_read_vols(V);  %  spm read image
        % con2--67*77*65double; XYZ--3*325325double;
        
        con2_data = con2(maskAAL>0); % get the voxles within the AAL mask only, this is a vector
        volumeshc{j, i} = con2_data; % Cell{sub26*load1234}--47634*1double
    end
end
volumes = [volumessz; volumeshc]; %(sz30+hc26)*4cell--64*4cell
% fitting the fitlme model for each location;
[m, n] = size(volumes); %56*4
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
    Tvalue(k, 1) = 0;
    pvalueD(k, 1) = 1;
    if rank(subidata)==4
        clear lme;
        dataROICov = [subidata(:, 1:4), SZ30HC26hehavMatchImagTable{:, 4},...
            SZ30HC26hehavMatchImagTable{:, 5}, SZ30HC26hehavMatchImagTable{:, 6},...
            SZ30HC26hehavMatchImagTable{:, 9}];% load1-4,age, sex,edu, fd;
        dataforModel =  [[zeros(56, 1); ones(56, 1); ones(56, 1)*2; ones(56, 1)*3],...
            repmat([ones(30, 1); zeros(26, 1)], 4, 1), repmat((1: 56)', 4, 1),... 
            [dataROICov(:, 1);  dataROICov(:, 2);  dataROICov(:, 3);  dataROICov(:, 4)],...
            repmat(dataROICov(:, 5), 4, 1), repmat(dataROICov(:, 6), 4, 1), ...
            repmat(dataROICov(:, 7), 4, 1), repmat(dataROICov(:, 8), 4, 1)];
        dataforModelTa = array2table(dataforModel, 'VariableNames', {'loadings', 'group',...
            'subject', 'activation', 'ageCov', 'sexCov', 'eduCov', 'fdCov'});
        lme = fitlme(dataforModelTa,...
            'activation~group:loadings+group+loadings+ageCov+sexCov+eduCov+fdCov+(1|subject) ');
         
         Tvalue(k, 1) = lme.Coefficients.tStat(8); % group:loadings
         pvalueD(k, 1) = lme.Coefficients.pValue(8); % 47634*1double
    end
end
fdr1 = mafdr(pvalueD, 'BHFDR', true);
sign1 = sum(fdr1<0.05); % no sign;
%% 4-HC Ushape: x ~ loandings^2 + loadings + age + sex + edu + meanFD + (1|subjID)
% 2022-4-20-16:16-19:21
WKD = ['~/DATA_FS6R3load1n10RightTrailRun3',...
    '/2ndlevelAnalysisFS6R3load1n10/2ndlevelanalysisSZ30hc26/HC26load1234Ftest'];

cellload = {'load1Conload0HC', 'load2Conload0HC', 'load3Conload0HC', 'load4Conload0HC'};
% read the image data Cell{sub26*load1234}
maskvol = spm_vol([WKD, '/aal3_3mmBinnocerellu.nii']);
maskAAL = spm_read_vols(maskvol); % AAL Mask
for i = 1:length(cellload)
    clear niifiles;
    loadpath = [WKD, filesep, cellload{i}];
    niifiles = dir([loadpath, '/*.nii']);
    for j = 1:length(niifiles)
        clear V;
        clear volume;
        V = spm_vol([loadpath, filesep, niifiles(j).name]);
        [con1, XYZ] = spm_read_vols(V);  %  spm read image
        % con1--67*77*65double; XYZ--3*325325double;
        
        con1_data = con1(maskAAL>0); % get the voxles within the AAL mask only, this is a vector
        volumes{j, i} = con1_data; % Cell{sub26*load1234}--47634*1double
    end
end
[m, n] = size(volumes); %26*4
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
    
    Tvalue(k, 1) = 0;
    pvalueD(k, 1) = 1;
    if rank(subidata)==4
        clear lme;
        dataROICov = [subidata(:, 1:4), SZ30HC26hehavMatchImagTable{31:56, 4},...
            SZ30HC26hehavMatchImagTable{31:56, 5}, SZ30HC26hehavMatchImagTable{31:56, 6},...
            SZ30HC26hehavMatchImagTable{31:56, 9}];% load1-4,age, sex,edu, fd;
        dataforModel =  [[zeros(26, 1); ones(26, 1); ones(26, 1)*2; ones(26, 1)*3],...
            repmat((1: 26)', 4, 1), [dataROICov(:, 1);  dataROICov(:, 2);  dataROICov(:, 3);  dataROICov(:, 4)],...
            repmat(dataROICov(:, 5), 4, 1), repmat(dataROICov(:, 6), 4, 1), ...
            repmat(dataROICov(:, 7), 4, 1), repmat(dataROICov(:, 8), 4, 1)];
        dataforModelTa = array2table(dataforModel, 'VariableNames', {'loading',	'subject',...
            'activation', 'ageCov', 'sexCov', 'eduCov', 'fdCov'});
        lme = fitlme(dataforModelTa,...
            'activation~loading^2+ageCov+sexCov+eduCov+fdCov+(1|subject) ');% including loading term;
        
         Tvalue(k, 1) = lme.Coefficients.tStat(7);
         pvalueD(k, 1) = lme.Coefficients.pValue(7); % 47634*1double
    end
end
fdr1 = mafdr(pvalueD, 'BHFDR', true);
sign1 = sum(fdr1<0.05);
figure;
min(abs((Tvalue(fdr1<0.05))));

%% write the data back to the brain (p-value or T-value)
signT = Tvalue(fdr1<0.05);
va = min(signT(signT>0));
vaN = max(signT(signT<0));% -4.53
% Tvalue£»
V.fname = 'model4fitlmeTvalueSign.nii';
Fcon = con1;
siTvalue = Tvalue;
siTvalue(fdr1>=0.05) = 0;
Fcon(maskAAL<1)= 0;
Fcon(maskAAL>0)= siTvalue;
spm_write_vol(V, Fcon);
% p-value Sign
V.fname = 'M4fitlmePvalueSign.nii';
Fcon = con1;
siPvalue = pvalueD;
siPvalue(fdr1>=0.05) = 0;
Fcon(maskAAL<1)= 0;
Fcon(maskAAL>0)= siPvalue;
spm_write_vol(V, Fcon);
% p-value1 raw
V.fname = 'M4fitlmelPvalue1.nii';
Fcon = con1;
siPvalue = pvalueD;
Fcon(maskAAL<1)= 1;
Fcon(maskAAL>0)= siPvalue;
spm_write_vol(V, Fcon);
% T-value raw
V.fname = 'M4fitlmeTvalue.nii';
Fcon = con1;
siPvalue = Tvalue;
Fcon(maskAAL<1)= 0;
Fcon(maskAAL>0)= siPvalue;
spm_write_vol(V, Fcon);
%% using 2-sz fdr significant area(clustersize 10) as mask for correction of model4 pvalue;
% 2022-4-24;
WKD = ['~/DATA_FS6R3load1n10RightTrailRun3',...
    '/2ndlevelAnalysisFS6R3load1n10/2ndlevelanalysisSZ30hc26'];
% read mask and the resulted model3 pvalue;
% mask2vol = spm_vol([WKD, '/SZ30load1234Ftest/model2fitlmePvalueSignClustersize10Mask_mask.nii']);
% mask2AAL = spm_read_vols(mask2vol); % Mask
maskvol = spm_vol([WKD, '/aal3_3mmBinnocerellu.nii']);
maskAAL = spm_read_vols(maskvol); % AAL Mask
Tvol = spm_vol([WKD, '/SZ30load1234Ftest/model2fitlmeTvalueSign.nii']);
TAAL = spm_read_vols(Tvol);
nozero = TAAL*0;
nozero(TAAL~=0 & maskAAL>0) = 1; % 1841 1

clear con3_data;
clear con3;
V3 = spm_vol([WKD, '/HC26load1234Ftest/M4fitlmelPvalue1.nii']);
[con3, XYZ3] = spm_read_vols(V3);  
% con3--67*77*65double; XYZ--3*325325double;
con3_data = con3(nozero~=0); %1841
% fdr
fdr3 = mafdr(con3_data, 'BHFDR', true);
sign3 = sum(fdr3<0.05); % sign3=0

%% 5-SZ Ushape: x ~ loandings^2 + loadings + age + sex + edu + meanFD + (1|subjID)
% 2022-4-20-16:16-19:21
WKD = ['~/DATA_FS6R3load1n10RightTrailRun3',...
    '/2ndlevelAnalysisFS6R3load1n10/2ndlevelanalysisSZ30hc26/SZ30load1234Ftest'];
cellload = {'load1Conload0SZ', 'load2Conload0SZ', 'load3Conload0SZ', 'load4Conload0SZ'};
% read the image data Cell{sub30*load1234}
maskvol = spm_vol([WKD, '/aal3_3mmBinnocerellu.nii']);
maskAAL = spm_read_vols(maskvol); % AAL Mask
for i = 1:length(cellload)
    clear niifiles;
    loadpath = [WKD, filesep, cellload{i}];
    niifiles = dir([loadpath, '/*.nii']);
    for j = 1:length(niifiles)
        clear V;
        clear volume;
        V = spm_vol([loadpath, filesep, niifiles(j).name]);
        [con1, XYZ] = spm_read_vols(V);  %  spm read image
        % con1--67*77*65double; XYZ--3*325325double;
        
        con1_data = con1(maskAAL>0); % get the voxles within the AAL mask only, this is a vector
        volumes{j, i} = con1_data; % Cell{sub30*load1234}--47634*1double
    end
end
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
    
    Tvalue(k, 1) = 0;
    pvalueD(k, 1) = 1;
    if rank(subidata)==4
        clear lme;
        dataROICov = [subidata(:, 1:4), SZ30HC26hehavMatchImagTable{1:30, 4},...
            SZ30HC26hehavMatchImagTable{1:30, 5}, SZ30HC26hehavMatchImagTable{1:30, 6},...
            SZ30HC26hehavMatchImagTable{1:30, 9}];% load1-4,age, sex,edu, fd;
        dataforModel =  [[zeros(30, 1); ones(30, 1); ones(30, 1)*2; ones(30, 1)*3],...
            repmat((1: 30)', 4, 1), [dataROICov(:, 1);  dataROICov(:, 2);  dataROICov(:, 3);  dataROICov(:, 4)],...
            repmat(dataROICov(:, 5), 4, 1), repmat(dataROICov(:, 6), 4, 1), ...
            repmat(dataROICov(:, 7), 4, 1), repmat(dataROICov(:, 8), 4, 1)];
        dataforModelTa = array2table(dataforModel, 'VariableNames', {'loading',	'subject',...
            'activation', 'ageCov', 'sexCov', 'eduCov', 'fdCov'});
        lme = fitlme(dataforModelTa,...
            'activation~loading^2+ageCov+sexCov+eduCov+fdCov+(1|subject) ');
        
         Tvalue(k, 1) = lme.Coefficients.tStat(7);
         pvalueD(k, 1) = lme.Coefficients.pValue(7); % 47634*1double
    end
end
fdr1 = mafdr(pvalueD, 'BHFDR', true);
sign1 = sum(fdr1<0.05); % 4--angular_R
figure;
hist(Tvalue(fdr1<0.05));
% write back
V.fname = 'model5fitlmeTvalueSign.nii';
Fcon = con1;
siTvalue = Tvalue;
siTvalue(fdr1>=0.05) = 0;
Fcon(maskAAL<1)= 0;
Fcon(maskAAL>0)= siTvalue;
spm_write_vol(V, Fcon);

%% 6-HC vs SCZ: x ~ loadings^2:group + loadings + age + sex + edu + meanFD + (1 | subjID )
% 11:17--18:10
WKD = ['~/DATA_FS6R3load1n10RightTrailRun3',...
    '/2ndlevelAnalysisFS6R3load1n10/2ndlevelanalysisSZ30hc26'];
WKD1 = [WKD, '/SZ30load1234Ftest'];
cellload1 = {'load1Conload0SZ', 'load2Conload0SZ', 'load3Conload0SZ', 'load4Conload0SZ'};
% read the image data Cell{sub30*load1234}
maskvol = spm_vol([WKD, '/aal3_3mmBinnocerellu.nii']);
maskAAL = spm_read_vols(maskvol); % AAL Mask
for i = 1:length(cellload1)
    clear niifiles;
    loadpath = [WKD1, filesep, cellload1{i}];
    niifiles = dir([loadpath, '/*.nii']);
    for j = 1:length(niifiles)
        clear V;
        clear volume;
        V = spm_vol([loadpath, filesep, niifiles(j).name]);
        [con1, XYZ] = spm_read_vols(V);  %  spm read image
        % con1--67*77*65double; XYZ--3*325325double;
        
        con1_data = con1(maskAAL>0); % get the voxles within the AAL mask only, this is a vector
        volumessz{j, i} = con1_data; % Cell{sub30*load1234}--47634*1double
    end
end
WKD2 = [WKD, '/HC26load1234Ftest'];
cellload2 = {'load1Conload0HC', 'load2Conload0HC', 'load3Conload0HC', 'load4Conload0HC'};
% read the image data Cell{sub26*load1234}
maskvol = spm_vol([WKD, '/aal3_3mmBinnocerellu.nii']);
maskAAL = spm_read_vols(maskvol); % AAL Mask
for i = 1:length(cellload2)
    clear niifiles;
    loadpath = [WKD2, filesep, cellload2{i}];
    niifiles = dir([loadpath, '/*.nii']);
    for j = 1:length(niifiles)
        clear V;
        clear volume;
        V = spm_vol([loadpath, filesep, niifiles(j).name]);
        [con2, XYZ] = spm_read_vols(V);  %  spm read image
        % con2--67*77*65double; XYZ--3*325325double;
        
        con2_data = con2(maskAAL>0); % get the voxles within the AAL mask only, this is a vector
        volumeshc{j, i} = con2_data; % Cell{sub26*load1234}--47634*1double
    end
end
volumes = [volumessz; volumeshc]; %(sz30+hc26)*4cell--64*4cell
% fitting the fitlme model for each location;
[m, n] = size(volumes); %56*4
Tvalue = zeros(length(volumes{1,1}), 1);
pvalueD = ones(length(volumes{1,1}), 1);
subidata = zeros(m, n);
parfor k = 1: length(volumes{1,1}) % for each location--47634
    display(['Computing the ', num2str(k), '-th location point!']);
   
    for colu = 1: n
        for row = 1: m
            subidata(row, colu) = volumes{row, colu}(k);
        end
    end
    subidata(isnan(subidata) == 1)=[];
  
    if rank(subidata)==4  % k=686;

        dataROICov = [subidata(:, 1:4), SZ30HC26hehavMatchImagTable{:, 4},...
            SZ30HC26hehavMatchImagTable{:, 5}, SZ30HC26hehavMatchImagTable{:, 6},...
            SZ30HC26hehavMatchImagTable{:, 9}];% load1-4,age, sex,edu, fd;
        dataforModel =  [[zeros(56, 1); ones(56, 1); ones(56, 1)*2; ones(56, 1)*3],...
            repmat([ones(30, 1); zeros(26, 1)], 4, 1), repmat((1: 56)', 4, 1),... 
            [dataROICov(:, 1);  dataROICov(:, 2);  dataROICov(:, 3);  dataROICov(:, 4)],...
            repmat(dataROICov(:, 5), 4, 1), repmat(dataROICov(:, 6), 4, 1), ...
            repmat(dataROICov(:, 7), 4, 1), repmat(dataROICov(:, 8), 4, 1)];
        dataforModelTa = array2table(dataforModel, 'VariableNames', {'loadings', 'group',...
            'subject', 'activation', 'ageCov', 'sexCov', 'eduCov', 'fdCov'});
        lme = fitlme(dataforModelTa,...
            'activation~group:loadings^2-group:loadings+group+loadings^2+ageCov+sexCov+eduCov+fdCov+(1|subject) ');
         
         Tvalue(k, 1) = lme.Coefficients.tStat(9); % group:loadings^2
         pvalueD(k, 1) = lme.Coefficients.pValue(9); % 47634*1double
    end
end
fdr1 = mafdr(pvalueD, 'BHFDR', true);
sign1 = sum(fdr1<0.05); % no sign;

%% extract the mean signal from M4 fdr sign 1 cluster ACC(clustersize=10);
% compute the mean signal from each cluster area  each sub;
% 2022-4-24
%%%%%%%%%%%%%%% hc26 %%%%%%%%%%%%%%%%%%%
WKD = ['~/DATA_FS6R3load1n10RightTrailRun3',...
    '/2ndlevelAnalysisFS6R3load1n10/2ndlevelanalysisSZ30hc26/HC26load1234Ftest'];
cellload = {'load1Conload0HC', 'load2Conload0HC', 'load3Conload0HC', 'load4Conload0HC'};
maskvol1 = spm_vol([WKD, '/model4fitlmeTvalueSignACC_mask.nii']);
maskAAL1 = spm_read_vols(maskvol1); 

for i = 1:length(cellload)
    clear niifiles;
    loadpath = [WKD, filesep, cellload{i}];
    niifiles = dir([loadpath, '/*.nii']);
    for j = 1:length(niifiles)
        clear V;
        clear volume;
        V = spm_vol([loadpath, filesep, niifiles(j).name]);
        [con1, XYZ] = spm_read_vols(V);  %  spm read image
        % con1--67*77*65double; XYZ--3*325325double;
        % cluster
        con1_data = con1(maskAAL1>0); %  this is a vector
        hcSIGNALACC(j, i) = nanmean(con1_data);
        
    end
end
%%%%%%%%%%%%%%% sz30 %%%%%%%%%%%%%%%%%%%
WKD = ['~/DATA_FS6R3load1n10RightTrailRun3',...
    '/2ndlevelAnalysisFS6R3load1n10/2ndlevelanalysisSZ30hc26/SZ30load1234Ftest'];
cellload = {'load1Conload0SZ', 'load2Conload0SZ', 'load3Conload0SZ', 'load4Conload0SZ'};
for i = 1:length(cellload)
    clear niifiles;
    loadpath = [WKD, filesep, cellload{i}];
    niifiles = dir([loadpath, '/*.nii']);
    for j = 1:length(niifiles)
        clear V;
        clear volume;
        V = spm_vol([loadpath, filesep, niifiles(j).name]);
        [con1, XYZ] = spm_read_vols(V);  %  spm read image
        % con1--67*77*65double; XYZ--3*325325double;
        % cluster
        con1_data = con1(maskAAL1>0); %  this is a vector
        szSIGNALACC(j, i) = nanmean(con1_data);
        
    end
end
%% mixed model for the cluster signal;
% standard U model fitted by the HC26 ACC cluster signals;
load('M4fitlmefdrsignACC_signal.mat');

signal1 = hcSIGNALACC(:, 1);% column vector
signal2 = hcSIGNALACC(:, 2);
signal3 = hcSIGNALACC(:, 3);
signal4 = hcSIGNALACC(:, 4);
meansignal1 = mean(signal1);
meansignal2 = mean(signal2);
meansignal3 = mean(signal3);
meansignal4 = mean(signal4);

loadx = [0; 1; 2; 3]'; % column vector;
loady = [signal1; signal2; signal3; signal4]; 

figure % bar image
bar( [0, 1, 2, 3], [meansignal1, meansignal2,meansignal3,...
    meansignal4], 0.4, 'facecolor', 'w');
hold on
set(gca, 'xTickLabel', {'load1','load2','load3', 'load4'});
ylim([-2, 2])
xlim([-1, 4])
hold on
% fitting the mixed model accounting for  subject
loadysub = [signal1, signal2, signal3, signal4];  % 26*4
mixplot = plot(loadx, loadysub, 'o', 'LineWidth', 2);
legend([repmat('sub ', 26, 1), num2str((1: 26)')],...
       'Location', 'NW')
grid on
hold on
model = @(phi, t)(phi(1)*(t.*t) +phi(2)*t+phi(3));% model
LOADX = repmat(loadx, 26, 1); % 26row*(0 1 2 3)---26*4
NUMS = repmat((1: 26)', size(loadx)); % (1...26)*4column--26*4
beta0 = [1 1 1];

[beta1,PSI1,stats1, b] = nlmefit(LOADX(:), loadysub(:), NUMS(:),...
                              [], model, beta0);
disp(char('Nonlinear mixed model of hc group: ', ['y= ', poly2str(beta1, 'x')], ...
    ['a = ', num2str(beta1(1)),'   b=', num2str(beta1(2)),...
    '   c = ', num2str(beta1(3))]));
phi = repmat(beta1, 1, 26) + ...          % Fixed effects
        [b(1, :); b(2, :); b(3, :)];               % random effects        
xplot = -1:0.1: 4;  
% plot(xplot, model(phi,xplot), 'k', 'LineWidth', 2)
applymodel = @(t)(beta1(1)*(t.*t) +beta1(2)*t+beta1(3)); 
plot(xplot, applymodel(xplot), 'g')
ylabel('load dlay contrast value mean ACC');
title('HC26 U curver: mixed model')
hold off;

% model1 y1 = ax+b, linear fitting
model1 = @(phi1, x)(phi1(1)*x+phi1(2)); % model1
beta0 = [1 1];
[beta2,PSI2,stats2, b1] = nlmefit(LOADX(:), loadysub(:),NUMS(:),...
                              [], model1, beta0);
disp(char('linear mixed model1 of hc group: ', ...
    ['a = ', num2str(beta2(1)),'   b=', num2str(beta2(2))]));
phi1 = repmat(beta2, 1, 26) + ...          % Fixed effects
        [b1(1, :); b1(2, :)];               % random effects

uLogL = stats1.logl; % x^2model
rLogL = stats2.logl; % bx+c model, restrict the x^2 factor;
dof = 1; % 
% Likelihood ratio test of model specification
[h,pValue,stat] = lratiotest(uLogL, rLogL, dof);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%standard U model fitted by the sz30 ACC cluster signals;
load('sz30M4fitlmefdrsignACC_signal.mat');

signal1 = szSIGNALACC(:, 1);% column vector
signal2 = szSIGNALACC(:, 2);
signal3 = szSIGNALACC(:, 3);
signal4 = szSIGNALACC(:, 4);
meansignal1 = mean(signal1);
meansignal2 = mean(signal2);
meansignal3 = mean(signal3);
meansignal4 = mean(signal4);

loadx = [0; 1; 2; 3]'; % column vector;
loady = [signal1; signal2; signal3; signal4]; 

figure % bar image
bar( [0, 1, 2, 3], [meansignal1, meansignal2,meansignal3,...
    meansignal4], 0.4, 'facecolor', 'w');
hold on
set(gca, 'xTickLabel', {'load1','load2','load3', 'load4'});
ylim([-2, 2])
xlim([-1, 4])
hold on
% fitting the mixed model accounting for  subject
loadysub = [signal1, signal2, signal3, signal4];  % 26*4
mixplot = plot(loadx, loadysub, 'o', 'LineWidth', 2);
legend([repmat('sub ', 30, 1), num2str((1: 30)')],...
       'Location', 'NW')
grid on
hold on
model = @(phi, t)(phi(1)*(t.*t) +phi(2)*t+phi(3));% model
LOADX = repmat(loadx, 30, 1); % 26row*(0 1 2 3)---26*4
NUMS = repmat((1: 30)', size(loadx)); % (1...26)*4column--26*4
beta0 = [1 1 1];

[beta1,PSI1,stats1, b] = nlmefit(LOADX(:), loadysub(:), NUMS(:),...
                              [], model, beta0);
disp(char('Nonlinear mixed model of hc group: ', ['y= ', poly2str(beta1, 'x')], ...
    ['a = ', num2str(beta1(1)),'   b=', num2str(beta1(2)),...
    '   c = ', num2str(beta1(3))]));
phi = repmat(beta1, 1, 30) + ...          % Fixed effects
        [b(1, :); b(2, :); b(3, :)];               % random effects        
xplot = -1:0.1: 4;  
% plot(xplot, model(phi,xplot), 'k', 'LineWidth', 2)
applymodel = @(t)(beta1(1)*(t.*t) +beta1(2)*t+beta1(3)); 
plot(xplot, applymodel(xplot), 'g')
ylabel('load dlay contrast value mean ACC');
title('SZ30 U curver: mixed model')
hold off;

% model1 y1 = ax+b, linear fitting
model1 = @(phi1, x)(phi1(1)*x+phi1(2)); % model1
beta0 = [1 1];
[beta2,PSI2,stats2, b1] = nlmefit(LOADX(:), loadysub(:),NUMS(:),...
                              [], model1, beta0);
disp(char('linear mixed model1 of hc group: ', ...
    ['a = ', num2str(beta2(1)),'   b=', num2str(beta2(2))]));
phi1 = repmat(beta2, 1, 30) + ...          % Fixed effects
        [b1(1, :); b1(2, :)];               % random effects
    
uLogL = stats1.logl; % x^2model
rLogL = stats2.logl; % bx+c model, restrict the x^2 factor;
dof = 1; % 
% Likelihood ratio test of model specification
[h,pValue,stat] = lratiotest(uLogL, rLogL, dof);
