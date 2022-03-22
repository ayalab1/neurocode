%% Single Cell Plots
% For plotting individual cell metrics

%% Settings
basepath = 'Y:\SMproject\AZ1\day13';
basename = basenameFromBasepath(basepath);
animName = animalFromBasepath(basepath);
% load([basepath '\' 'Barrage_Files' '\' basename '.HSE.mat']);
% load([basepath '\' 'Barrage_Files' '\' basename '.HSE.mat']);
load([basepath '\' basename '.spikes.cellinfo.mat']);
cellNum = 56;

if ~exist(strcat('Z:\home\Lindsay\Barrage\Cell\SingleCellPlots\', animName))
    mkdir(strcat('Z:\home\Lindsay\Barrage\Cell\SingleCellPlots\', animName));
end
saveSchem = strcat('Z:\home\Lindsay\Barrage\Cell\SingleCellPlots\',animName,'\',basename,'.',num2str(cellNum));

%%
%single cell ISIs distribution
BinSize_ISIs=1e-3; %in seconds
ntotal_ISIs=100; %counts
[~, n, ~] = ISIGrams(spikes.times, spikes.UID, BinSize_ISIs, ntotal_ISIs);
%plot for one random cell - eventually, we group them by region and average
figure();
plot(n(:,cellNum)./sum(n(:,cellNum)));%normalize per the total n of events
title(['ISI distribution for unit #' num2str(cellNum)]);
ylabel('Count normalized to total number of events');
xlabel('ISI (s)');
saveas(gcf,[saveSchem '.ISIdist.png']);