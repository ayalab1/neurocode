%% Barrage Analysis
% This is currently for use in singular sessions, but can easily be
% expanded using the outline from MassSessRun.m
%
%
% TO DO LIST:
% -- Create folders and files for saving the metrics run from this (likely
%    should be something in a central location, notated then by animal and
%    day)
% -- Create a folder with plots notated in the name by animal and day
% -- Integrate python to create nicer plots
% -- Actually write code for the metrics

%% Settings
basepath = pwd;
basename = basenameFromBasepath(basepath);
savePath = strcat(basepath, '\Barrage_Files');

% We should run this analysis with ALL spikes and pull using UIDs later if
% we want to look closer at regions or cell types
load([basename '.spikes.cellinfo.mat']);
load([basepath '\Barrage_Files\' basename '.HSEfutEVT.mat']);
load([basepath '\' basename '.cell_metrics.cellinfo.mat']);
fc = 0;

%% Initial runs
detection = 8;
% Get baseline metrics during and between events using pre-made code
[spkEventTimes] = getRipSpikes('basepath',pwd,'events',evtSave{detection,1},'spikes',spikes,'savePath',savePath,'padding',0);

[unitBar] = unitSWRmetrics(spkEventTimes,spikes);
save(strcat(savePath,'\',basename,'.unitBar.mat'),'unitBar');

%%% BURST PROPERTIES [SINGLE CELL, 10ms] %%%

%% # spikes in 10 ms windows
% Compute per unit, plot
% Histogram, then smooth - spkrthist?
%%%%% should check to see what they're envisioning bc I don't think that
%%%%% the histogram smoothed is very useful
binsz = 10/1000;
tSmooth = 150*binsz;
avgRate = NaN(1,length(spikes.times));
% fc = fc+1
% figure(fc);
% hold on
for i = 1:length(spikes.times)
    [spkhist spkmean spkstd ts] = spkRtHist(spikes.times{i}, tSmooth, binsz);
%     plot(ts, spkhist);
    avgRate(i) = mean(spkhist);
end
% hold off
fc = fc+1;
figure(fc);
plot(sort(avgRate, 'descend'));

%% ISI per unit
% Distribution
ISI = NaN(1,length(spikes.times));
for i = 1:length(spikes.times)
    num = 0;
    tempTime = spikes.times{i};
    for j = 1:length(tempTime)-1
        num = num+(tempTime(j+1)-tempTime(j));
    end
    ISI(i) = num/(length(tempTime)-1);
end
fc = fc+1;
figure(fc);
plot(1:length(ISI), ISI, '.');
xlabel('Unit #');
ylabel('ISI (s)');

[regID cellID] = getSubs(cell_metrics);
[ISIs, pI] = sort(ISI);
for i = 1:length(ISIs)
    regIDs(i) = regID(pI(i));
    cellIDs(i) = cellID(pI(i));
end 
fc=fc+1;
plotSubs(1:length(ISIs), ISIs, regIDs, cellIDs, fc);
xlabel('Unit #');
ylabel('ISI (s)');
hold off
%% FR per unit
% Distribution

%% Fraction of events that each unit fires in for barrages
% Plot (unit# by fraction of participation)
[regID cellID] = getSubs(cell_metrics);

[particip, pI] = sort(unitBar.particip, 'descend');
for i = 1:length(particip)
    regIDs(i) = regID(pI(i));
    cellIDs(i) = cellID(pI(i));
end 

fc = fc+1;
plotSubs(1:length(particip),particip, regIDs, cellIDs, fc);
ylim([0 1]);
title('Fraction of participation per unit');
ylabel('Fraction of parcipation (max 1)');
xlabel('Sorted Unit # (descending participation)');
hold off
%%% POPULATION BURST PROPERTIES [POPULATION, 100ms] %%%
%% # units per event

%% % units from different subregions active during these events

%% # spikes each unit fires in each event

%% FR per unit (per event)


%%% BARRAGE PROPERTIES [POPULATION AND SINGLE CELL, 1s] %%%
%% # units per event

%% % units from different subregions active during these events

%% # spikes each unit fires in each event

%% FR per unit (per event)



%%%%%%% FUNCTIONS %%%%%%%
function [regID cellID] = getSubs(cell_metrics)
regions = cat(1,cell_metrics.brainRegion{:});
regID = ones(length(regions),1);
check = ["CA1" "CA2" "CA3" "CTX" "DG"];
regYes = zeros(length(check),length(regions));

for i = 1:size(regYes,1)
    for j = 1:size(regYes,2)
        if contains(regions(j,:),convertStringsToChars(check(i)))
            regID(j) = (i+1);
        end
    end
end

for i = 1:length(cell_metrics.putativeCellType)
    cellType(i) = convertCharsToStrings(cell_metrics.putativeCellType{i});
end
cellID = ones(length(cellType),1);
check = ["Pyramidal" "Interneuron"];
cellYes = zeros(length(check),length(cellType));

for i = 1:size(cellYes,1)
    for j = 1:size(cellYes,2)
        if contains(cellType(j),convertStringsToChars(check(i)))
            cellID(j) = (i+1);
        end
    end
end
end

function plotSubs(x, y, regID, cellID, fc)
color = ['m';'b';'g';'r';'c';'k'];
type = ['*';'^';'o'];

figure(fc)
hold on
for i = 1:length(cellID)
    plot(x(i), y(i), cat(1,color(regID(i)),type(cellID(i))));
end
h = zeros(8, 1);
h(1) = plot(NaN,NaN,'b');
h(2) = plot(NaN,NaN,'g');
h(3) = plot(NaN,NaN,'r');
h(4) = plot(NaN,NaN,'c');
h(5) = plot(NaN,NaN,'k');
h(6) = plot(NaN,NaN,'k^');
h(7) = plot(NaN,NaN,'ko');
h(8) = plot(NaN,NaN,'*m');
legend(h, 'CA1','CA2','CA3','CTX','DG','Pyramidal','Interneuron','Unknown');
end