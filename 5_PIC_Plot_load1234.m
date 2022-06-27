%% 2022-5-19 code for picture of SZ U-shape__manuscript2
%% %%1_plot the load0 &load4 ACC-2022-3-20%%%%%%%
signal1 = ACCdataMatchImageSZ28hc21Ta{1:28, 3};% load0
signal2 = ACCdataMatchImageSZ28hc21Ta{29:49, 3};
signal3 = ACCdataMatchImageSZ28hc21Ta{1:28, 5};% load1
signal4 = ACCdataMatchImageSZ28hc21Ta{29:49, 5};
signal5 = ACCdataMatchImageSZ28hc21Ta{1:28, 7};% load2
signal6 = ACCdataMatchImageSZ28hc21Ta{29:49, 7};
signal7 = ACCdataMatchImageSZ28hc21Ta{1:28, 9};% load3
signal8 = ACCdataMatchImageSZ28hc21Ta{29:49, 9};
signal9 = ACCdataMatchImageSZ28hc21Ta{1:28, 11};% load4
signal10 = ACCdataMatchImageSZ28hc21Ta{29:49, 11};

figure(1)
%plot(2, meanACCSCZ, '-ro', 1, meanACCHC, '-ro');
b = bar( [2, 3, 4, 5, 6], [mean(signal1), mean(signal2); mean(signal3), mean(signal4);...
    mean(signal5), mean(signal6); mean(signal7), mean(signal8); mean(signal9), mean(signal10);...
    ], 0.9, 'FaceColor', 'flat');
set(b(1), 'FaceColor', sixteen2ten('#82ADD6')/255,...
    'EdgeColor',sixteen2ten('#82ADD6')/255);
set(b(2), 'FaceColor', sixteen2ten('#FACB66')/255,...
    'EdgeColor',sixteen2ten('#FACB66')/255);

hold on
plot(ones(28,1)*1.85 + rand(28,1)*0.2-0.1, signal1, ....
     'Color', sixteen2ten('#06397D')/255, 'Marker', 'o', 'LineStyle', 'none'); %ro,go
hold on
plot(ones(21,1)*2.15 + rand(21,1)*0.2-0.1, signal2,...
     'Color', sixteen2ten('#9A3B4E')/255, 'Marker', 'o', 'LineStyle', 'none');
plot(ones(28,1)*2.85 + rand(28,1)*0.2-0.1, signal3,...
     'Color', sixteen2ten('#06397D')/255, 'Marker', 'o', 'LineStyle', 'none'); %ro,go
hold on
plot(ones(21,1)*3.15 + rand(21,1)*0.2-0.1, signal4,...
    'Color', sixteen2ten('#9A3B4E')/255, 'Marker', 'o', 'LineStyle', 'none');
hold on;
plot(ones(28,1)*3.85 + rand(28,1)*0.2-0.1, signal5, ....
     'Color', sixteen2ten('#06397D')/255, 'Marker', 'o', 'LineStyle', 'none'); %ro,go
hold on
plot(ones(21,1)*4.15 + rand(21,1)*0.2-0.1, signal6,...
     'Color', sixteen2ten('#9A3B4E')/255, 'Marker', 'o', 'LineStyle', 'none');
plot(ones(28,1)*4.85 + rand(28,1)*0.2-0.1, signal7,...
     'Color', sixteen2ten('#06397D')/255, 'Marker', 'o', 'LineStyle', 'none'); %ro,go
hold on
plot(ones(21,1)*5.15 + rand(21,1)*0.2-0.1, signal8,...
    'Color', sixteen2ten('#9A3B4E')/255, 'Marker', 'o', 'LineStyle', 'none');
hold on;
plot(ones(28,1)*5.85 + rand(28,1)*0.2-0.1, signal9, ....
     'Color', sixteen2ten('#06397D')/255, 'Marker', 'o', 'LineStyle', 'none'); %ro,go
hold on
plot(ones(21,1)*6.15 + rand(21,1)*0.2-0.1, signal10,...
     'Color', sixteen2ten('#9A3B4E')/255, 'Marker', 'o', 'LineStyle', 'none');

hold on;
ylabel('Accuracy (%)', 'FontSize', 12);
set(gca, 'xTickLabel', {'Load0', 'Load1', 'Load2', 'Load3', 'Load4'}, 'FontSize', 12);
ylim([50, 105])
yticks([50 80 100]);
xlim([1, 7])
hold off
legend('SZ', 'HC');

%% %%%%U shape mixed effects models
load('M4fitlmefdrsignACC_signal.mat');
load('sz30M4fitlmefdrsignACC_signal.mat');
figure(1);

subplot(2, 2, 3);
hold on;
signal1 = hcSIGNALACC(:, 1);% column vector
signal2 = hcSIGNALACC(:, 2);
signal3 = hcSIGNALACC(:, 3);
signal4 = hcSIGNALACC(:, 4);

loadx = [0; 1; 2; 3]'; % column vector;
loady = [signal1; signal2; signal3; signal4]; 

 % bar image
bar( [0, 1, 2, 3], [mean(signal1), mean(signal2),mean(signal3),...
    mean(signal4)], 0.4, 'facecolor', 'w', 'EdgeColor',sixteen2ten('#47484C')/255);
hold on
set(gca, 'xTickLabel', {'Load1', 'Load2','Load3', 'Load4'},...
    'FontSize', 12, 'FontWeight', 'bold');
ylim([-2, 2]);
yticks([-2 0 2]);
xlim([-1, 4]);

hold on
% fitting the mixed model accounting for  subject
loadysub = [signal1, signal2, signal3, signal4];  % 30*4
loadxran = [ones(26,1)*0 + rand(26,1)*0.4-0.2, ones(26,1)*1 + rand(26,1)*0.4-0.2,...
    ones(26,1)*2 + rand(26,1)*0.4-0.2, ones(26,1)*3 + rand(26,1)*0.4-0.2];
mixplot = plot(loadxran, loadysub, 'o', 'Color', sixteen2ten('#5E616D')/255, 'LineWidth', 1); 

grid on
hold on
model = @(phi, t)(phi(1)*(t.*t) +phi(2)*t+phi(3));% model
LOADX = repmat(loadx, 26, 1); % 30row*(0 1 2 3)---26*4
NUMS = repmat((1: 26)', size(loadx)); % (1...29)*4column--26*4
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
plot(xplot, applymodel(xplot), 'Color', sixteen2ten('#FACB66')/255,...
    'LineWidth', 1.5);
ylabel('HC','FontSize', 12, 'FontWeight', 'bold');
hold on

subplot(2, 2, 4);
signal1 = szSIGNALACC(:, 1);% column vector
signal2 = szSIGNALACC(:, 2);
signal3 = szSIGNALACC(:, 3);
signal4 = szSIGNALACC(:, 4);

loadx = [0; 1; 2; 3]'; % column vector;
loady = [signal1; signal2; signal3; signal4]; 

bar( [0, 1, 2, 3], [mean(signal1), mean(signal2),mean(signal3),...
    mean(signal4)], 0.4, 'facecolor', 'w');
hold on
set(gca, 'xTickLabel', {'Load1', 'Load2','Load3', 'Load4'},...
    'FontSize', 12, 'FontWeight', 'bold');
ylim([-2, 2])
yticks([-2 0 2]);
xlim([-1, 4])
hold on
% fitting the mixed model accounting for  subject
loadysub = [signal1, signal2, signal3, signal4];  % 30*4
loadxran = [ones(30,1)*0 + rand(30,1)*0.4-0.2, ones(30,1)*1 + rand(30,1)*0.4-0.2,...
    ones(30,1)*2 + rand(30,1)*0.4-0.2, ones(30,1)*3 + rand(30,1)*0.4-0.2];
mixplot = plot(loadxran, loadysub, 'o', 'Color', sixteen2ten('#5E616D')/255, 'LineWidth', 1); 

grid on
hold on
model = @(phi, t)(phi(1)*(t.*t) +phi(2)*t+phi(3));% model
LOADX = repmat(loadx, 30, 1); % 30row*(0 1 2 3)---30*4
NUMS = repmat((1: 30)', size(loadx)); % (1...30)*4column--30*4
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
plot(xplot, applymodel(xplot), 'Color', sixteen2ten('#82ADD6')/255,...
    'LineWidth', 1.5)
ylabel('SZ', 'FontSize', 12, 'FontWeight', 'bold');
% title('HC U curver: mixed model-right angular')
hold off;

%%%%%%%%%%%%%%SZ28 axis of symmetry--load 4 ACC;
indexsz = [1:10, 11:20, 22, 23, 25:30]; 
ACCcol = [11];
figure
% set(0, 'defaultfigurecolor', 'w'); 
for i = 1: length(ACCcol) % axis-acc Residual
   
    X = [ones(28, 1) SZ30HC26hehavMatchImagTable{indexsz, 4:6}];
    [b,bint,r] = regress(szfixRandABCpeak{indexsz, 4}, X); % axis
    [b1,bint1,r1] = regress(ACCdataMatchImageSZ28hc21Ta{1:28, ACCcol(i)},X);% ACC
    
    theta = glmfit(r, r1);
    yCalc1 = theta(2)*r +theta(1);
    plot(r(1:28, 1), r1(1:28, 1), 'Color', sixteen2ten('#66A9C9')/255,...
        'Marker','o', 'Markersize',10, 'LineStyle', 'none', 'LineWidth', 1);
    hold on;
%     plot(r(29:49, 1), r1(29:49, 1), 'Color', sixteen2ten('#D1C2D3')/255,...
%         'Marker','o', 'LineStyle', 'none');
%     hold on
    plot(r(:, 1), yCalc1, 'Color', sixteen2ten('#82ADD6')/255, 'LineWidth', 2.5);
    ylim([-30, 25]);
    yticks([-30 0 25]);
    xlim([-3, 4]);
    xticks([-3 0 4]);
    xticklabels({'-3', '0', '4'});
    ylabel('Residual of load 4 accuracy', 'FontSize', 12, 'FontWeight', 'bold');
    xlabel('Residual axis of symmetry', 'FontSize', 12, 'FontWeight', 'bold');
    grid on
end
%%%%%%%%%%%%%%SZ30 panss--who T cluster(fusiform, angular);
indexrow = [1:23, 25:30]; 
clustercol = [11, 17];
% set(0,‘defaultfigurecolor’,‘w’);
figure
for i = 1: length(clustercol) % mean in cluster Residual
    subplot(1, 2, i);
    
    X = [ones(29, 1) SZ30HC26hehavMatchImagTable{indexrow, 3:5}];
    [b,bint,r] = regress(szhcSignal56(indexrow, clustercol(i)), X); % mean in cluster 
    [b1,bint1,r1] = regress(SZ30HC26hehavMatchImagTable{indexrow, 10},X);% panss-P
    
    theta = glmfit(r, r1);
    yCalc1 = theta(2)*r +theta(1);
 
    plot(r(1:29, 1), r1(1:29, 1), 'Color', sixteen2ten('#66A9C9')/255,...
        'Marker','o', 'Markersize',10, 'LineStyle', 'none', 'LineWidth', 1);
    hold on;

    plot(r(:, 1), yCalc1, 'Color', sixteen2ten('#82ADD6')/255, 'LineWidth', 2.5);
    ylim([-10, 10]);
    yticks([-10 0 10]);
    xlim([-3, 2]);
    xticks([-3 0 2]);
    xticklabels({'-3', '0', '2'});
    ylabel('Residual of PANSS-P', 'FontSize', 12, 'FontWeight', 'bold');
    xlabel('Residual axis of symmetry', 'FontSize', 12, 'FontWeight', 'bold');
    grid on
end