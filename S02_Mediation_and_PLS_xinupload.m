clear;
clc;

% cd('F:\work_dir\AHBAenrich/');
work_dir =['~/argenbrain/PLSenrichZscore'];
load([work_dir, '/AHBAenrichcao/AHBA_Mean_scaled_reanote.mat']); 
load ([work_dir, '/AHBAenrichcao/AHBA_data_reanote_slim.mat']);
load([work_dir, '/AHBAenrichcao/AHBA_ROI_index_80.mat']);
%% M4 extract T value;
Yraw = readtable([work_dir, '/M4pointSignalTvalue.xlsx']); 
Yraw.Properties.VariableNames
% % cortical
% X=cort_l_expMS;%1105samples*15408 genes
% Y=table2array(Yraw(Ncort_l_80(:,1),5));%1105 rows.
% subcortical
X=sub_l_expMS;% 645samples(6 subject total subcortex)*15408 genes
Y=table2array(Yraw(Nsub_l_80(:,1),5));%645 rows, select 5th column

nozeroID = find(Y~=0);
X = X(nozeroID, :);
Y = Y(nozeroID, :);

X = zscore(X);
Y = zscore(Y);
%% PLS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLS variance picture start
%perform full PLS and plot variance in Y explained by top 20 components
%(hopefully!)
dim=20; % weishu (more larger more better)
[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(X,Y,dim);
% XL:15408*20;YL:1:20; XS: 645*20; YS: 645*20;PCTVAR: 2*20;
%PCTVAR containing the percentage of variance explained by the model. 
%The first row of PCTVAR contains the percentage of variance explained in X by each PLS component, 
%and the second row contains the percentage of variance explained in Y
figure
plot(1:dim,cumsum(100*PCTVAR(2,1:dim)),'-o','LineWidth',1.5,'Color',[140/255,0,0]);
set(gca,'Fontsize',14)
xlabel('Number of PLS components','FontSize',14);
ylabel('Percent Variance Explained in Y','FontSize',14);
text(1:dim,cumsum(100*PCTVAR(2,1:dim))-3,num2str(cumsum(100*PCTVAR(2,1:dim)).','%.1f'))
grid on
PCTVAR

figure
plot(1:dim,cumsum(100*PCTVAR(1,1:dim)),'-o','LineWidth',1.5,'Color',[140/255,0,0]);
set(gca,'Fontsize',14)
xlabel('Number of PLS components','FontSize',14);
ylabel('Percent Variance Explained in X','FontSize',14);

grid on

%%% plot correlation of PLS component 1 with t-statistic:
figure
plot(XS(:,1),Y,'r.')% the predictor scores XS, that is, the PLS components that are linear combinations of the variables in X
[R,P] = corrcoef(XS(:,1),Y);% select the component explained the most y;
xlabel('XS scores for PLS component 1','FontSize',14);
ylabel('Y','FontSize',14);
grid on
l = lsline;
l.Color = 'k';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLS variance picture end
%% permutation testing to assess significance of PLS result as a function of
% the number of components (dim) included:
rep=5000; 
dim = 5;
[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(X,Y,dim);
PCTVAR
Rsquared_real_y = cumsum(100*PCTVAR(2,1:dim));
Rsq_y = [];
Pvalue_perm_y = ones(1,dim);
for i = 1:dim % component
    i
    parfor j = 1:rep % automatic parallel calculate
    order=randperm(size(Y,1)); % row
    Xp=X(order,:);
    [XL_p,YL_p,XS_p,YS_p,BETA_p,PCTVAR_p,MSE_p,stats_p]=plsregress(Xp,Y,dim);
    temp_y=cumsum(100*PCTVAR_p(2,1:dim));
    Rsq_y(i,j) = temp_y(i);
    end
    Pvalue_perm_y(i) = length(find(Rsq_y(i,:)>=Rsquared_real_y(i)))/rep;
end

fig = figure;
plot(1:dim, Pvalue_perm_y,'ok','MarkerSize',8,'MarkerFaceColor','r');
xlabel('Number of PLS components','FontSize',14);
ylabel('p-value','FontSize',14);
grid on
% img = frame2in(getframe(fig));
% imwrite(img,'test.jpg');

%% Bootstrap to get the gene list:
%in order to estimate the error in estimating each gene's PLS1 weight, then
%to caculate the z score
%number of bootstrap iterations:
bootnum=5000;

% Do PLS in 2 dimensions (with 2 components):
dim=2;
[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(X,Y,dim);

% temp_dim = 4;
%store regions' IDs and weights in descending order of weight for both components:
[R1,p1]=corr([XS(:,1),XS(:,2)],Y);%XS: predicted X score

%align PLS components with desired direction for interpretability 
if R1(1,1)<0  %this is specific to the data shape we were using - will need ammending
    stats.W(:,1)=-1*stats.W(:,1);
    XS(:,1)=-1*XS(:,1);
end
if R1(2,1)<0 %this is specific to the data shape we were using - will need ammending
    stats.W(:,2)=-1*stats.W(:,2);
    XS(:,2)=-1*XS(:,2);
end

% save PLS_XandY X Y
%real weight and sorted by weight
[PLS1w,x1] = sort(stats.W(:,1),'descend');%W: A p-by-ncomp matrix of PLS weights so that XS = X0*W.
PLS1ids=genes(x1);
[PLS2w,x2] = sort(stats.W(:,2),'descend');
PLS2ids=genes(x2);

%define variables for storing the (ordered) weights from all bootstrap runs
PLS1weights=[];
PLS2weights=[];

%start bootstrap
parfor i=1:bootnum
    myresample = randsample(size(X,1),size(X,1),1);%replacement
    res(i,:)=myresample; %store resampling out of interest
    Xr=X(myresample,:); % define X for resampled subjects
    Yr=Y(myresample,:); % define X for resampled subjects
    [XL_b,YL_b,XS_b,YS_b,BETA_b,PCTVAR_b,MSE_b,stats_b]=plsregress(Xr,Yr,dim); %perform PLS for resampled data
      
    temp1=stats_b.W(:,1);%extract PLS1 weights
    newW1=temp1(x1); %order the newly obtained weights the same way as initial PLS 
    if corr(PLS1w,newW1)<0 % the sign of PLS components is arbitrary - make sure this aligns between runs
        newW1=-1*newW1;
    end
    PLS1weights=[PLS1weights,newW1];%store (ordered) weights from this bootstrap run
    
    temp2=stats_b.W(:,2);%extract PLS2 weights
    newW2=temp2(x2); %order the newly obtained weights the same way as initial PLS 
    if corr(PLS2w,newW2)<0 % the sign of PLS components is arbitrary - make sure this aligns between runs
        newW2=-1*newW2;
    end
    PLS2weights=[PLS2weights,newW2];%store (ordered) weights from this bootstrap run

end

%get standard deviation of weights from bootstrap runs
PLS1sw=std(PLS1weights');
PLS2sw=std(PLS2weights');

%get bootstrap weights (Z)
PLS1Z=PLS1w./PLS1sw';
PLS2Z=PLS2w./PLS2sw';

% save AHBA_PLS_bootstrap PLS1ids PLS1weights PLS2weights PLS1sw PLS2sw PLS2ids PLS1Z PLS2Z -v7.3;

PLS_out_genes = table(genes, geneSymbol, PLS1Z);
PLS2_out_genes = table(genes, geneSymbol, PLS2Z);

cd (work_dir)
writetable(sortrows(PLS_out_genes,-3),'PLS_out_genes_M4_subcort_PLS_variance1.csv');
writetable(sortrows(PLS2_out_genes,-3),'PLS_out_genes_M4_subcort_PLS_variance2.csv');
cd '..'


%% leave one donor out;-20220704 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% leave one donor out;-202207017 XINNNNNNNN%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cd('F:\work_dir\AHBAenrich/');
work_dir =['/share/inspurStorage/home1/Fengnana/DATA210725/SZ2010data/2020-SCZ',...
    '/DATAfMRIprepAGAIN2112/bidsfMRIprepSub/DATA_FS6R3load1n10RightTrailRun3/argenbrain/PLSenrichZscore'];
load([work_dir, '/AHBAenrichcao/AHBA_Mean_scaled_reanote.mat']); % /AHBAenrichcao;
load ([work_dir, '/AHBAenrichcao/AHBA_data_reanote_slim.mat']);
load([work_dir, '/AHBAenrichcao/AHBA_ROI_index_80.mat']);

%% M4 extract T value;
Yraw = readtable([work_dir, '/M2pointSignalTvalue.xlsx']); 
Yraw.Properties.VariableNames
% cortical
Xall=cort_l_expMS;%1105samples*15408 genes
Yall=table2array(Yraw(Ncort_l_80(:,1),5));%1105 rows.
% donorN= [307, 247, 110, 148, 157, 136];
Donor = {[308:1105], [1:307,555:1105], [1:554, 665:1105], [1:664, 813:1105],...
    [1:812, 970:1105], [1:969]};

nozeroID = find(Yall~=0);%938/950
zeroID = find(Yall==0);

Xall1 = Xall(nozeroID, :);
Yall1 = Yall(nozeroID, :);

Xall1 = zscore(Xall1);%938/950
Yall1 = zscore(Yall1);% no zero for M4/M2;

dim=20; % weishu (more larger more better)
% PLS(all donor)
[XLa,YLa,XSa,YSa,BETAa,PCTVARa,MSEa,statsa]=plsregress(Xall1, Yall1, dim);
% back write 938/950 zscore -nozero X/Y to num 1105, default set the inital y=0 location as 0;
for i = 1:length(Yall1)
    Xallback(nozeroID(i), :) = Xall1(i, :);% 1105*15408
    Yallback(nozeroID(i), 1) = Yall1(i, 1); % 1105*1
end

for d = 1: 6
    clear X;
    clear Y;
    % leave one donor each step;
    X = Xallback(Donor{1, d}, :);
    Y = Yallback(Donor{1, d}, :);
    
    nozeroID = find(Y~=0);
    X = X(nozeroID, :);
    Y = Y(nozeroID, :);

    %% PLS
    clear stats;
    [XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(X,Y,dim);
    
    figure(1)
    subplot(2,3,d);
    plot(1:dim,cumsum(100*PCTVAR(2,1:dim)),'-o','LineWidth',1,'Color',[140/255,0,0]);
    set(gca,'Fontsize',14)
    xlabel(['Num of PLS compon (leave ', num2str(d), '-th)'],'FontSize',11);
    ylabel('Percent Variance Explained in Y','FontSize',11);
    % text(1:dim,cumsum(100*PCTVAR(2,1:dim))-3,num2str(cumsum(100*PCTVAR(2,1:dim)).','%.1f'))
    grid on
    PCTVAR
    
    %%% plot correlation of PLS(leave one donor out) coefficient PLS (all donor) coefficient:
    figure(2)
    subplot(2, 3, d)
    plot(stats.W(:,1),statsa.W(:,1),'r.')
    [R,P] = corrcoef(stats.W(:,1),statsa.W(:,1));% first component;
    RP(1,d) = R(1,2);
    RP(2,d) = P(1,2);
    xlabel('PLS(all donor)coef for component1','FontSize',11);
    ylabel(['PLS(leave ', num2str(d), '-th)coef for component1'],'FontSize',11);
    grid on
    l = lsline;
    l.Color = 'k';
end
% %%%%%%%% subcortical XIN%%%%%%%%%%
Yraw = readtable([work_dir, '/M2pointSignalTvalue.xlsx']); 
Yraw.Properties.VariableNames
Xall=sub_l_expMS;% 645samples(6 subject total subcortex)*15408 gene;
Yall=table2array(Yraw(Nsub_l_80(:,1),5));%645 rows, select 5th column;
% donorN= [195, 199, 49, 78, 64, 60];
Donor = {[196:645], [1:195,395:645], [1:394, 444:645], [1:443, 522:645],...
    [1:521, 586:645], [1:585]};

nozeroID = find(Yall~=0);
zeroID = find(Yall==0);

Xall1 = Xall(nozeroID, :);
Yall1 = Yall(nozeroID, :);

Xall1 = zscore(Xall1);
Yall1 = zscore(Yall1);% no zero for M4/M2;

dim=20; % weishu (more larger more better)
[XLa,YLa,XSa,YSa,BETAa,PCTVARa,MSEa,statsa]=plsregress(Xall1, Yall1, dim);

for i = 1:length(Yall1)
    Xallback(nozeroID(i), :) = Xall1(i, :);% 645*15408
    Yallback(nozeroID(i), 1) = Yall1(i, 1); % 645*1
end

for d = 1: 6
    clear X;
    clear Y;
    X = Xallback(Donor{1, d}, :);
    Y = Yallback(Donor{1, d}, :);
    
    nozeroID = find(Y~=0);
    X = X(nozeroID, :);
    Y = Y(nozeroID, :);
    %% PLS
    clear stats;
    [XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(X,Y,dim);
    
    figure(1)
    subplot(2,3,d);
    plot(1:dim,cumsum(100*PCTVAR(2,1:dim)),'-o','LineWidth',1,'Color',[140/255,0,0]);
    set(gca,'Fontsize',14)
    xlabel(['Num of PLS compon (leave ', num2str(d), '-th)'],'FontSize',11);
    ylabel('Percent Variance Explained in Y','FontSize',11);
    % text(1:dim,cumsum(100*PCTVAR(2,1:dim))-3,num2str(cumsum(100*PCTVAR(2,1:dim)).','%.1f'))
    grid on
    PCTVAR
    
    %%% plot correlation of PLS(leave one donor out) coefficient PLS (all donor) coefficient:
    figure(2)
    subplot(2, 3, d)
    plot(stats.W(:,1),statsa.W(:,1),'r.')
    [R,P] = corrcoef(stats.W(:,1),statsa.W(:,1));
    RP(1,d) = R(1,2);
    RP(2,d) = P(1,2);
    xlabel('PLS(all donor)coef for component1','FontSize',11);
    ylabel(['PLS(leave ', num2str(d), '-th)coef for component1'],'FontSize',11);
    grid on
    l = lsline;
    l.Color = 'k';
end