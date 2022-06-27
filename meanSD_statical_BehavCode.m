%% the resulted SZ30 HC26 behav 
% statistical analysis again examine£»
% 2022-6-21
% compute mean/SD of behavioral data
colu = [4,6,9,67,68,69,70];
for i = 1:length(colu)
    SZmeansd(i,1) = nanmean(SZ30HC26hehavMatchImagTable{1:30, colu(i)});
    SZmeansd(i,2) = nanstd(SZ30HC26hehavMatchImagTable{1:30, colu(i)});
    HCmeansd(i,1) = nanmean(SZ30HC26hehavMatchImagTable{31:56, colu(i)});
    HCmeansd(i,2) = nanstd(SZ30HC26hehavMatchImagTable{31:56, colu(i)});
end
% t-test
colut = [4,6,9];
for j= 1: length(colut)
    [h,p,ci,stats] = ttest2(SZ30HC26hehavMatchImagTable{1:30, colut(j)},...
        SZ30HC26hehavMatchImagTable{31:56, colut(j)});
    tt(j, 1) = stats.tstat; % ttest T value;
    tt(j, 2) = p; % ttest P value
end
% anova1
anova1([SZ30HC26hehavMatchImagTable{31:56, 9};...
    SZ30HC26hehavMatchImagTable{1:30, 9}],...
    [zeros(1,26), ones(1,30)])

%% organize ACC data matching imaging subID;
% 2022-2-17
[m, n] = size(TASKWMaccuracyHCSCZ);
[a, b] = size(SZ30HC27hehavMatchImagTable);
for i = 1: a
    for j = 1: m
    if SZ30HC27hehavMatchImagTable{i, 1} == TASKWMaccuracyHCSCZ{j,1};
        ACCdataMatchImageSZ30hc22(i,:) = TASKWMaccuracyHCSCZ{j,:};
    end
    end
end
ACCdataMatchImageSZ30hc22Ta = array2table(ACCdataMatchImageSZ30hc22, 'VariableNames',...
    {'subID',	'group', 'cond0ACC', 'cond0RT',  'cond1ACC', 'cond1RT', 'cond2ACC', 'cond2RT',...
     'cond3ACC', 'cond3RT', 'cond4ACC', 'cond4RT','meanACC'});
%% compute mean/SD of ACC behavior data
for i = 3: 13
    SZaccmeansd(i-2, 1) = nanmean(ACCdataMatchImageSZ28hc21Ta{1:28, i});
    SZaccmeansd(i-2, 2) = nanstd(ACCdataMatchImageSZ28hc21Ta{1:28, i});
    HCaccmeansd(i-2, 1) = nanmean(ACCdataMatchImageSZ28hc21Ta{29:49, i});
    HCaccmeansd(i-2, 2) = nanstd(ACCdataMatchImageSZ28hc21Ta{29:49, i});
end
% test without covariates
for ki= 3: 13
    [h,p,ci,stats] = ttest2(ACCdataMatchImageSZ28hc21Ta{1:28, ki},...
       ACCdataMatchImageSZ28hc21Ta{29:49, ki});
    tt(ki-2, 1) = stats.tstat; % ttest T value;
    tt(ki-2, 2) = p; % ttest P value
end
% anova with covariates(age, sex, eduyear)
index = [1:20,22, 23, 25:30, 31:35, 37:39, 42:46, 48:49, 51:56]; % index sz30hc26 to sz28hc21;
for kj =3: 13
    data2comp = [ACCdataMatchImageSZ28hc21Ta{1:28, kj};...
        ACCdataMatchImageSZ28hc21Ta{29:49, kj}];
    grouplabel = [zeros(28,1); ones(21,1)];
    
    groups = [grouplabel,SZ30HC26hehavMatchImagTable{index, 4},...
       SZ30HC26hehavMatchImagTable{index, 5},...
       SZ30HC26hehavMatchImagTable{index, 6}];...
    
    [p,tb1] = anovan(data2comp, groups, 'continuous', [2, 4]);
    anovanP(kj-2, :) = p'; % extract p value
    anovanf(kj-2, :) = tb1{2, 6}';
end 

% anova with covariates(age, sex, eduyear,load0ACC)
% 2022-2-19
for kj =4: 13
    data2comp = [ACCdataMatchImageSZ28hc21Ta{1:28, kj};...
        ACCdataMatchImageSZ28hc21Ta{29:49, kj}];
    grouplabel = [zeros(28,1); ones(21,1)];
    
    groups = [grouplabel,SZ30HC26hehavMatchImagTable{index, 4},...
       SZ30HC26hehavMatchImagTable{index, 5},...
       SZ30HC26hehavMatchImagTable{index, 6},...
        ACCdataMatchImageSZ28hc21Ta{:, 3}];
    [p, tb1,stats] = anovan(data2comp, groups, 'continuous', [2, 4, 5]);
    anovanPload0accCov(kj-3, :) = p'; % extract p value
    anovanfload0accCov(kj-3, :) = tb1{2, 6}';
end

%% compute correlation of ACC and behavior data(PANSS) 
% sz33 hc22
indexSZ = [1:20,22, 23, 25:30];
for i =3: 13 %¡¡£Á£Ã£Ã
    for j = 67:70 % panss P N G
        jco =(j-66)*2;
        % parCorrRP--11*8
        [parCorrRP(i-2, jco-1), parCorrRP(i-2, jco)]=partialcorr(ACCdataMatchImageSZ28hc21Ta{1:28, i},...
            SZ30HC26hehavMatchImagTable{indexSZ,j}, [SZ30HC26hehavMatchImagTable{indexSZ, 4},...
            SZ30HC26hehavMatchImagTable{indexSZ, 5}, SZ30HC26hehavMatchImagTable{indexSZ, 6}]); % Cov Age, sex, edu;
    end
end
%% attention-related analysis (TMTA-WM/TMTB)
indexsz = [1:10, 11:20, 22, 23, 25:30]; 
indexhc = [31:35, 37:39, 42:46, 48:49, 51:56];

% compute mean and SD value;
colu = [21,23];
for i = 1:length(colu)
    SZmeansd(i,1) = nanmean(SZ30HC26hehavMatchImagTable{1:30, colu(i)});
    SZmeansd(i,2) = nanstd(SZ30HC26hehavMatchImagTable{1:30, colu(i)});
    HCmeansd(i,1) = nanmean(SZ30HC26hehavMatchImagTable{31:56, colu(i)});
    HCmeansd(i,2) = nanstd(SZ30HC26hehavMatchImagTable{31:56, colu(i)});
end

anova1([SZ30HC26hehavMatchImagTable{31:56, 21};...
    SZ30HC26hehavMatchImagTable{1:30, 21}],...
    [zeros(1,26), ones(1,30)])% TMTA/B-RT
anova1([SZ30HC26hehavMatchImagTable{31:56, 23};...
    SZ30HC26hehavMatchImagTable{1:30, 23}],...
    [zeros(1,26), ones(1,30)])% TMTA/B-WN
%%anova WITH covariates--
for i = 21:24
    data2comp = [SZ30HC26hehavMatchImagTable{1:30, i};...
        SZ30HC26hehavMatchImagTable{31:56, i}];
    grouplabel = [ones(30, 1);zeros(26, 1)];
    groups = [grouplabel,SZ30HC26hehavMatchImagTable{:, 4},...
        SZ30HC26hehavMatchImagTable{:, 5},SZ30HC26hehavMatchImagTable{:, 6}];
    p= anovan(data2comp, groups, 'continuous', [2, 4]);
    pf(i, 1) = p(1);
end
% correlation-SZ
for  i =21:22
    for j = 3: 13
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
% correlation-HC
for  i =21:22
    for j = 3: 13
        % covariates age, sex, eduyear;
        ageCov = SZ30HC26hehavMatchImagTable{indexhc, 4};
        sexCov = SZ30HC26hehavMatchImagTable{indexhc, 5};
        eduCov = SZ30HC26hehavMatchImagTable{indexhc, 6};
       
        colu = (j-2)*2;
        [R, P] = partialcorr(SZ30HC26hehavMatchImagTable{indexhc, i},...
            ACCdataMatchImageSZ28hc21Ta{29:49, j},...
            [ageCov, sexCov, eduCov],'Type', 'Spearman','Rows','Complete');
        RP(i-20, colu-1) = R;
        RP(i-20, colu) = P;
    end
end
%%scatter:TMTA-RT---load4 ACC;
indexhc = [31:35, 37:39, 42:46, 48:49, 51:56];
figure
X = [ones(21, 1) SZ30HC26hehavMatchImagTable{indexhc, 4:6}];
[b,bint,r] = regress(SZ30HC26hehavMatchImagTable{indexhc, 21}, X); % TMTA-RT
[b1,bint1,r1] = regress(ACCdataMatchImageSZ28hc21Ta{29:49, 11}, X);% Load4-ACC

theta = glmfit(r, r1);
yCalc1 = theta(2)*r +theta(1);
plot(r(:, 1), r1(:, 1), 'Color', sixteen2ten('#B6825F')/255,...
    'Marker','o', 'Markersize',10, 'LineStyle', 'none', 'LineWidth', 1);
hold on;

plot(r(:, 1), yCalc1, 'Color', sixteen2ten('#FACB66')/255, 'LineWidth', 2);
ylim([-10, 5]);
yticks([-10 0 5]);
xlim([-15, 30]);
xticks([-15 0 30]);
% xticklabels({'-4', '0', '4'});
ylabel('Residual accuracy of load4', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('Residual TMT-RT', 'FontSize', 12, 'FontWeight', 'bold');
% title(['Correlation of slope with load', num2str(i), ' ACC']);
grid on

%% FDR correction;
fdr1 = mafdr(pvalue, 'BHFDR', true);
sign1 = sum(fdr1<0.05); 
%% compute particial correlation;
[parCorr(1, 1), parCorr(1, 2)]=partialcorr(ACCdataMatchImageSZ30hc22Ta{:, 3},...
            ACCdataMatchImageSZ30hc22Ta{:,12}, [SZ28HC22hehavMatchACC{:, 4},...
            SZ28HC22hehavMatchACC{:, 5}, SZ28HC22hehavMatchACC{:, 6}]); %