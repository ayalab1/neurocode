% two sample t test: 
% null: vectors x and y comes from independent random samples from normal 
%       distributions with equal means and equal but unknown variances, using the 
%       two-sample t-test. 
% The alternative hypothesis is that the data in x and y 
%       comes from populations with unequal means. The result h is 1 if the test rejects 
%       the null hypothesis at the 5% significance level, and 0 otherwise.
% [h,p,ci,stats] = ttest2(depth_total_E14_2M2,depth_total_E14_2M3)
%% plot normalized depths - separate
close all
figure
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.05, 0.6, 0.96]);
subplot(1,2,1)
title('CA1 Pyramidal cell depth (histogram)')
hold on
nBins=50;
edges=linspace(-3,3,nBins);
[fraction, edges] = histcounts(depth_total_E13_norm, edges, 'Normalization', 'probability');plot(fraction,edges(1:end-1)+diff( edges ) / 2, 'linewidth', 2)
[fraction, edges] = histcounts(depth_total_E14_2M3_norm, edges, 'Normalization', 'probability');plot(fraction,edges(1:end-1)+diff( edges ) / 2, 'linewidth', 2)
[fraction, edges] = histcounts(depth_total_E14_2M2_norm, edges, 'Normalization', 'probability');plot(fraction,edges(1:end-1)+diff( edges ) / 2, 'linewidth', 2)
[fraction, edges] = histcounts(depth_total_E15_11M1_norm, edges, 'Normalization', 'probability');plot(fraction,edges(1:end-1)+diff( edges ) / 2, 'linewidth', 2)
[fraction, edges] = histcounts(depth_total_E16_1F1_norm, edges, 'Normalization', 'probability');plot(fraction,edges(1:end-1)+diff( edges ) / 2, 'linewidth', 2)

% histogram(depth_total_E13_norm,'Orientation','horizontal','BinWidth',0.1, 'linestyle', 'none','FaceColor',[0.4 0.5940 0.6],'Normalization','probability');
% histogram(depth_total_E14_2M3_norm,'Orientation','horizontal','BinWidth',0.1, 'linestyle', 'none','FaceColor',[0.2 0.3940 0.8],'Normalization','probability');
% histogram(depth_total_E14_2M2_norm,'Orientation','horizontal','BinWidth',0.1, 'linestyle', 'none','FaceColor',[0.5 0.240 0.8],'Normalization','probability');
% histogram(depth_total_E15_11M1_norm,'Orientation','horizontal','BinWidth',0.1, 'linestyle', 'none','FaceColor',[0.3 0.540 0.2],'Normalization','probability');
% histogram(depth_total_E16_1F1_norm,'Orientation','horizontal','BinWidth',0.1, 'linestyle', 'none','FaceColor',[0.9290 0.6940 0.1250],'Normalization','probability');
set(gca, 'XAxisLocation', 'bottom')
xlabel('Percentage');
ylabel('Depth (normalized)');
ylim([-1,2])
yline(1)
yline(0)
text(0.21,1.2,'Oriens')
text(0.21,0.5,'Pyramidale')
text(0.21,-0.2,'Radiatum')
legend('E13','E14_{2M3}','E14_{2M2}','E15_{11M1}','E16_{1F1}')
% legend('E13','E14_{2M2}','E14_{2M3}','E16')
legend boxoff

%% box plot - separate
% loop through allllll var in the workspace
disp('please clean the workspace first! Only leave depth data')
data_all = whos;
length_max=0;
for ii = 1:numel(data_all)
    if length(eval(data_all(ii).name))>length_max
        length_max=length(eval(data_all(ii).name));
    end
end
data_Mat = NaN(length_max,5);
data_Mat(1:length(depth_total_E13_norm),1)=depth_total_E13_norm;
data_Mat(1:length(depth_total_E14_2M3_norm),2)=depth_total_E14_2M3_norm;
data_Mat(1:length(depth_total_E14_2M2_norm),3)=depth_total_E14_2M2_norm;
data_Mat(1:length(depth_total_E15_11M1_norm),4)=depth_total_E15_11M1_norm;
data_Mat(1:length(depth_total_E16_1F1_norm),5)=depth_total_E16_1F1_norm;
% figure
subplot(1,2,2)

% boxplot(data_Mat,'Notch','off','Labels',{'E13','E14_2M2','E14_2M3','E16'},'Whisker',1,'OutlierSize',3)
boxplot(data_Mat,'Notch','off','Labels',{'E13','E14_2M3','E14_2M2','E15_11M1','E16_1F1'},'Whisker',1,'OutlierSize',2)
box off
title('CA1 Pyramidal cell depth (box plot)')

%% plot normalized depths - combined
close all
figure
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.05, 0.6, 0.96]);
subplot(1,2,1)
title('CA1 Pyramidal cell depth (histogram)')
hold on
nBins=50;
edges=linspace(-3,3,nBins);


[fraction, edges] = histcounts([depth_total_E13_1M1_norm',depth_total_E13_1F1_norm',depth_total_E13_16F1_norm'], edges, 'Normalization', 'probability');plot(fraction,edges(1:end-1)+diff( edges ) / 2, 'linewidth', 2)
[fraction, edges] = histcounts([depth_total_E14_2M2_norm',depth_total_E14_2M3_norm',depth_total_E14_1B1_norm',depth_total_E14_2B1_norm'], edges, 'Normalization', 'probability');plot(fraction,edges(1:end-1)+diff( edges ) / 2, 'linewidth', 2)
[fraction, edges] = histcounts([depth_total_E15_11M1_norm',depth_total_E15_131_norm',depth_total_E15_132_norm'], edges, 'Normalization', 'probability');plot(fraction,edges(1:end-1)+diff( edges ) / 2, 'linewidth', 2)
[fraction, edges] = histcounts([depth_total_E16_1F1_norm',depth_total_E16_2M1_norm', depth_total_E16_3F1_norm',depth_total_E16_3F2_norm'], edges, 'Normalization', 'probability');plot(fraction,edges(1:end-1)+diff( edges ) / 2, 'linewidth',2)

% histogram(depth_total_E13_norm,'Orientation','horizontal','BinWidth',0.1, 'linestyle', 'none','FaceColor',[0.4 0.5940 0.6],'Normalization','probability');
% histogram([depth_total_E14_2M2_norm',depth_total_E14_2M3_norm'],'Orientation','horizontal','BinWidth',0.1, 'linestyle', 'none','FaceColor',[0.2 0.3940 0.8],'Normalization','probability');
% % histogram(depth_total_E14_2M2_norm,'Orientation','horizontal','BinWidth',0.1, 'linestyle', 'none','FaceColor',[0.5 0.240 0.8],'Normalization','probability');
% histogram(depth_total_E15_11M1_norm,'Orientation','horizontal','BinWidth',0.1, 'linestyle', 'none','FaceColor',[0.3 0.540 0.2],'Normalization','probability');
% histogram(depth_total_E16_1F1_norm,'Orientation','horizontal','BinWidth',0.1, 'linestyle', 'none','FaceColor',[0.9290 0.6940 0.1250],'Normalization','probability');
set(gca, 'XAxisLocation', 'bottom')
xlabel('Percentage');
ylabel('Depth (normalized)');
yline(1)
yline(0)
ylim([-1,2])
text(0.21,1.2,'Oriens')
text(0.21,0.5,'Pyramidale')
text(0.21,-0.2,'Radiatum')
legend('E13','E14','E15','E16')
% legend('E13','E14_{2M2}','E14_{2M3}','E16')
legend boxoff

%% box plot - combined
% loop through allllll var in the workspace
disp('please clean the workspace first! Only leave depth data')
data_all = whos;
length_max=0;
for ii = 1:numel(data_all)
    if length(eval(data_all(ii).name))>length_max
        length_max=length(eval(data_all(ii).name));
    end
end
data_Mat = NaN(length_max,4);
data_Mat(1:length([depth_total_E13_1M1_norm',depth_total_E13_1F1_norm',depth_total_E13_16F1_norm']),1)=[depth_total_E13_1M1_norm',depth_total_E13_1F1_norm',depth_total_E13_16F1_norm'];
data_Mat(1:length([depth_total_E14_2M2_norm',depth_total_E14_2M3_norm',depth_total_E14_1B1_norm',depth_total_E14_2B1_norm']),2)=[depth_total_E14_2M2_norm',depth_total_E14_2M3_norm',depth_total_E14_1B1_norm',depth_total_E14_2B1_norm'];
data_Mat(1:length([depth_total_E15_11M1_norm',depth_total_E15_131_norm',depth_total_E15_132_norm']),3)=[depth_total_E15_11M1_norm',depth_total_E15_131_norm',depth_total_E15_132_norm'];
data_Mat(1:length([depth_total_E16_1F1_norm',depth_total_E16_2M1_norm', depth_total_E16_3F1_norm',depth_total_E16_3F2_norm']),4)=[depth_total_E16_1F1_norm',depth_total_E16_2M1_norm', depth_total_E16_3F1_norm',depth_total_E16_3F2_norm'];
% figure
subplot(1,2,2)

% boxplot(data_Mat,'Notch','off','Labels',{'E13','E14_2M2','E14_2M3','E16'},'Whisker',1,'OutlierSize',3)
boxplot(data_Mat,'Notch','off','Labels',{'E13','E14','E15','E16'},'Whisker',1,'OutlierSize',2)
ylim([-1,2])
box off
title('CA1 Pyramidal cell depth (box plot)')

%% plot cdf - separate
figure
cdfplot(medial_lateral_position_total_E13)
hold on
cdfplot(medial_lateral_position_total_E14_2M3)
cdfplot(medial_lateral_position_total_E15_11M1)
cdfplot(medial_lateral_position_total_E16_1F1)
title("CDF of mediolateral distribution")
xlabel('Mediolateral index');
ylabel('Percentage');
legend('E13_{all}','E14_{2M3}','E15_{11M1}','E16_{1F1}')
xlim([-0.1,1.1])
box off

%% plot cdf - combined
figure
cdfplot(depth_total_E13_norm)
hold on
cdfplot([depth_total_E14_2M3_norm' depth_total_E14_2M2_norm'])
cdfplot(depth_total_E15_11M1_norm)
cdfplot(depth_total_E16_1F1_norm)
xlabel('Depth (normalized)');
ylabel('Percentage');
legend('E13_{1F1}','E14_{combined}','E15_{11M1}','E16_{1F1}')
xlim([-0.5,2])
box off
% %% plot absolute depths
% close all
% figure
% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.05, 0.6, 0.96]);
% subplot(1,2,1)
% title('CA1 Pyramidal cell depth (histogram)')
% hold on
% histogram(depth_total_13,'Orientation','horizontal','BinWidth',5, 'linestyle', 'none','FaceColor',[0.4 0.5940 0.6],'Normalization','probability');
% histogram(depth_total_E14_2M2,'Orientation','horizontal','BinWidth',5, 'linestyle', 'none','FaceColor',[0.2 0.3940 0.8],'Normalization','probability');
% histogram(depth_total_E14_2M3,'Orientation','horizontal','BinWidth',5, 'linestyle', 'none','FaceColor',[0.3 0.540 0.2],'Normalization','probability');
% histogram(depth_total_16,'Orientation','horizontal','BinWidth',5, 'linestyle', 'none','FaceColor',[0.9290 0.6940 0.1250],'Normalization','probability');
% set(gca, 'XAxisLocation', 'bottom')
% xlabel('Percentage');
% ylabel('Depth (um)');
% yline(50)
% yline(0)
% text(0.1,75,'Oriens')
% text(0.1,25,'Pyramidale')
% text(0.1,-25,'Radiatum')
% legend('E13','E14_{2M2}','E14_{2M3}','E16')
% legend boxoff
% 
% % end
% 
% %% box plot absolute depth
% % loop through allllll var in the workspace
% data_all = whos;
% length_max=0;
% for ii = 1:numel(data_all)
%     if length(eval(data_all(ii).name))>length_max
%         length_max=length(eval(data_all(ii).name));
%     end
% end
% data_Mat = NaN(length_max,4);
% data_Mat(1:length(depth_total_13),1)=depth_total_13;
% data_Mat(1:length(depth_total_E14_2M2),2)=depth_total_E14_2M2;
% data_Mat(1:length(depth_total_E14_2M3),3)=depth_total_E14_2M3;
% data_Mat(1:length(depth_total_16),4)=depth_total_16;
% % figure
% subplot(1,2,2)
% 
% boxplot(data_Mat,'Notch','off','Labels',{'E13','E14_2M2','E14_2M3','E16'},'Whisker',1,'OutlierSize',3)
% box off
% title('CA1 Pyramidal cell depth (box plot)')