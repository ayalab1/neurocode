%% Barrage Analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% [STEP 0] Settings
close all
savePath = ('Z:\home\Lindsay\Barrage\');
load('Z:\home\Lindsay\Barrage\combinedPaths.mat');
ses = paths(16); %can change this to an iteratable variable for mass runs
cd(ses);
basepath = pwd;
basename = basenameFromBasepath(basepath);
animName = animalFromBasepath(basepath);
% Load in session information
load([basename '.spikes.cellinfo.mat']);
load([basename '.cell_metrics.cellinfo.mat']);
load([basename '.session.mat']);
showPlt = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% [STEP 1]
%%% BURST PROPERTIES [SINGLE CELL, 10ms] %%%
specPath = strcat(savePath,'Cell\');
saveSchem = strcat(animName,'.',basename);
if exist(strcat(specPath,'props.mat'))
    load([specPath 'props.mat']);
end
[regID cellID] = getSubs(cell_metrics);
cellProp.regID = regID;
cellProp.cellID = cellID;

%% [1.1] ISI per unit
% Distribution
[ISI,ISIc,t] = ISIGrams(spikes.times, spikes.UID, 1/1000, 100);
cellProp.ISI = ISI;
cellProp.ISIhistT = t;
cellProp.ISIhist = ISIc;
for i = 1:length(ISI)
    ISIavg(i) = mean(ISI{i});
end
cellProp.ISIavg = ISIavg;

% Plotting
figure('Position', get(0, 'Screensize'));
plot(1:length(ISIavg), ISIavg*1e3, '.');
xlabel('Unit #');
ylabel('Log of avg ISI (ms)');
title('Average ISI per cell');
set(gca, 'YScale', 'log')
saveas(gcf,[specPath 'Plots\' saveSchem '.ISIavg.png']);

% By type plotting
% [ISIs, pI] = sort(ISIavg);
% for i = 1:length(ISIs)
%     regIDs(i) = regID(pI(i));
%     cellIDs(i) = cellID(pI(i));
% end 
plotSubs(1:length(ISIavg), ISIavg*1e3, regID, cellID);
xlabel('Unit #');
ylabel('Log of avg ISI (ms)');
title('Average ISI per cell by type');
set(gca, 'YScale', 'log')
hold off
saveas(gcf,[specPath 'Plots\' saveSchem '.ISIavgSort.png']);

if ~showPlt
    close all
end

%% [1.2] # spikes in 10 ms windows
% Bin ISI's across session (maybe also make plots per neuron type and per
% region)
% Per cell
figure('Position', get(0, 'Screensize'));
hold on
for i = 1:size(ISIc,2)
    plot(t*1000, ISIc(:,i));
end
xlim([t(1)*1000 t(end)*1000]);
title('Histogram of ISIs per cell');
ylabel('Count');
xlabel('ISI (ms)');
hold off
saveas(gcf,[specPath 'Plots\' saveSchem '.ISIdist.png']);

% Per cell with type designations
plotSubs(t*1000, ISIc, regID, cellID, '-');
xlim([t(1)*1000 t(end)*1000]);
title('Histogram of ISIs by cell type and region');
ylabel('Count');
xlabel('ISI (ms)');
saveas(gcf,[specPath 'Plots\' saveSchem '.ISIdistSort.png']);

% Per cell with type designations, shortened
plotSubs(t*1000, ISIc, regID, cellID, '-');
xlim([t(1)*1000 20]);
title('Histogram of ISIs by cell type and region');
ylabel('Count');
xlabel('ISI (ms)');
saveas(gcf,[specPath 'Plots\' saveSchem '.ISIdistSortShort.png']);

% Summary box plot per region
check = ["CA1" "CA2" "CA3" "CTX" "DG"]; %this is embedded in getSubs()
regBox = cell(length(check),1);
boxC = cell(length(check),1);
for i = 2:(length(check)+1)
    if find(regID==i, 1)
        [keep] = includeType(regID,i);
        inds = find(keep==1);
        temp = [];
        for j = 1:length(inds)
            temp = [temp; ISI{inds(j)}];
        end
        regBox{i-1} = temp;
        boxC{i-1} = (i-1)*ones(length(temp),1);
    end
end
boxx = [];
boxg = [];
for i = 1:length(regBox)
    boxx = [boxx; regBox{i}];
    boxg = [boxg; boxC{i}];
end
figure('Position', get(0, 'Screensize'));  
% boxplot(boxx*1000, boxg,'symbol', ''); %use symbol '' to suppress outliers
boxplot(boxx*1000, boxg);
useName = unique(boxg);
nameInd = ismember(1:length(check),useName);
presentIND = find(nameInd==1);
presentName = convertStringsToChars(check(nameInd));
xticklabels(presentName);
title('Box plot per region, outliers cut off');
ylabel('ISI (ms)');
xlabel('Region');
ylim([-1000 7000]);
saveas(gcf,[specPath 'Plots\' saveSchem '.ISIboxReg.png']);



% Summary box plot per region and per type

if ~showPlt
    close all
end

%% [1.3] FR per unit (we could also pull this from cellmetrics if it takes too long)
% Distribution
FR = cell(length(spikes.times),1);
avgFR = NaN(length(spikes.times),1);
for i = 1:length(spikes.times)
    spkHist = spkRtHist(sort(spikes.times{i}), 'ifz', false);
    FR{i} = spkHist/0.001; %get rate from our count/bin size
    avgFR(i) = sum(FR{i})/length(FR{i});
end
cellProp.FR = FR;
cellProp.avgFR = avgFR;

figure('Position', get(0, 'Screensize'));
plot(1:length(avgFR), avgFR, '.');
title('Average Firing Rate');
ylabel('Firing Rate (Hz)');
xlabel('Unit Number');
saveas(gcf,[specPath 'Plots\' saveSchem '.avgFR.png']);

plotSubs(1:length(avgFR), avgFR, regID, cellID);
title('Average Firing Rate by Type');
ylabel('Firing Rate (Hz)');
xlabel('Unit Number');
saveas(gcf,[specPath 'Plots\' saveSchem '.avgFRSort.png']);

if ~showPlt
    close all
end

%% [1.4] Burst Index per cell
binssum = 6; %6ms
burstIndex = NaN(length(spikes.times),1);
for i = 1:length(spikes.times)
    burstIndex(i) = sum(ISIc(1:binssum,i))/sum(ISIc(:,i));
end
cellProp.burstIndex = burstIndex;

%% [1.5] Cell SAVE
save(strcat(specPath,'props.mat'),'cellProp', '-v7.3');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% [STEP 2]
%%% BARRAGE DETECTION PROPERTIES %%%
% Need to keep event duration in mind
specPath = strcat(savePath,'Population\');
load([basepath '\Barrage_Files\' basename '.HSEfutEVT.mat']);
load([basepath '\Barrage_Files\' basename '.HSEmetrics.mat']);
detection = 8; %%%%% need to pick the detection we want
saveSchem = strcat(animName,'.',basename,'.',num2str(detection));
if exist(strcat(specPath,'props.mat'))
    load([specPath 'props.mat']);
end

%% [2.0] Initial runs
[spkEventTimes] = getRipSpikes('basepath',pwd,'events',evtSave{detection,1},'spikes',spikes,'savePath',savePath,'saveNum',detection,'padding',0,'saveMat',false);
save(strcat(specPath,saveSchem,'.spkEventTimes.mat'),'spkEventTimes', '-v7.3');
[unitBar] = unitSWRmetrics(spkEventTimes,spikes);
save(strcat(specPath,saveSchem,'.unitBar.mat'),'unitBar', '-v7.3');

%% [2.1] # units per event
% This should be calculated within unitBar
barProp.untPerEvt = unitBar.nCellsEvent;
barProp.untPerEvtID = unitBar.nCellsEventUID;

[sortCell] = sort(unitBar.nCellsEvent,'descend');
figure('Position', get(0, 'Screensize'));
bar(1:length(sortCell),sortCell);
title('Number of units per event');
xlabel('Event # (sorted descending)');
ylabel('Number of Units');
saveas(gcf,[specPath 'Plots\' saveSchem '.untPerEvt.png']);

[~,edges,bins] = histcounts(spkEventTimes.EventDuration,25);
durCount = zeros(max(bins),1);
binct = zeros(max(bins),1);
for i = 1:length(bins)
    durCount(bins(i)) = durCount(bins(i))+unitBar.nCellsEvent(i);
end
for i = 1:length(binct)
    binct(i) = length(find(bins==i));
    if binct(i)==0
        binct(i) = 1;
    end
end
durCount = durCount./(binct);
for i = 1:length(edges)-1
    durX(i) = edges(i) + ((edges(i+1)-edges(i))/2);
end
figure('Position', get(0, 'Screensize'));
bar(durX,durCount);
title('Average number of units per event, sorted by event duration');
xlabel('Event duration (s)');
xticks([durX]);
ylabel('Average number of units');
saveas(gcf,[specPath 'Plots\' saveSchem '.untPerEvtDur.png']);

if ~showPlt
    close all
end

%% [2.2] % units from different subregions active during these events

regPerc = zeros(length(presentIND),length(unitBar.nCellsEventUID));
for i = 1:length(unitBar.nCellsEventUID) %iterate through each event
    temp = [];
    temp = unitBar.nCellsEventUID{i}; %pull unit IDs that fire during event i
    for j = 1:length(temp)
        useInd = find(presentIND==(regID(temp(j))-1)); %find
        regPerc(useInd,i) = regPerc(useInd,i)+1;
    end
    regSum = sum(regPerc(:,i));
    for k = 1:size(regPerc,1)
        regPerc(k,i) = 100*(regPerc(k,i)/regSum);
    end
end

barProp.regPerc = regPerc;

figure('Position', get(0, 'Screensize'));
bar(regPerc.','stacked');
title('Percent of units from each region per event');
xlabel('Event #');
ylabel('Percent');
legend(presentName, 'Location', 'eastoutside');
ylim([0 100]);
saveas(gcf,[specPath 'Plots\' saveSchem '.percEvtReg.png']);

% CA2 will be recorded everywhere theoretically, but just make sure
if ismember(2,presentIND)
    [regSortCA2, rI] = sort(regPerc((presentIND==2),:));
    for i = 1:length(regSortCA2)
        for j = 1:length(presentIND)
            regSort(j,i) = regPerc(j,(rI(i)));
            regSort(j,i) = regPerc(j,(rI(i)));
            regSort(j,i) = regPerc(j,(rI(i)));
        end
    end
    figure('Position', get(0, 'Screensize'));
    bar(regSort.','stacked');
    title(['Percent of units from each region per event, sorted by increasing CA2 involvment']);
    xlabel('Event #');
    ylabel('Percent');
    legend(presentName, 'Location', 'eastoutside');
    ylim([0 100]);
    saveas(gcf,[specPath 'Plots\' saveSchem '.percEvtRegCA2sort.png']);
end

% Sort by duration
[~,edges,bins] = histcounts(spkEventTimes.EventDuration);
durCount = zeros(max(bins),length(presentIND));
binct = zeros(max(bins),1);
for i = 1:length(bins)
    for j = 1:length(presentIND)
        durCount(bins(i),j) = durCount(bins(i),j)+regPerc(j,i);
    end
end
for i = 1:length(binct)
    binct(i) = length(find(bins==i));
    if binct(i)==0
        binct(i) = 1;
    end
end
for i = 1:size(durCount,1)
    durCount(i,:) = durCount(i,:)/(binct(i));
end
durX = [];
for i = 1:length(edges)-1
    durX(i) = edges(i) + ((edges(i+1)-edges(i))/2);
end
figure('Position', get(0, 'Screensize'));
bar(durX,durCount,'stacked');
title('Region participation percentage by event duration');
xlabel('Event duration (s)');
xticks([durX]);
ylabel('Percent participation in events');
ylim([0 100]);
legend(presentName, 'Location', 'eastoutside');
saveas(gcf,[specPath 'Plots\' saveSchem '.untPercentDur.png']);

% Summary box plot per region

boxx = [];
boxg = [];
for i = 1:size(regPerc,1)
    boxx = [boxx; regPerc(i,:)];
    boxg = [boxg; i*ones(size(regPerc,2),1)];
end
figure('Position', get(0, 'Screensize'));  
boxplot(boxx, boxg);
title('Percent of units which make up events for each region');
ylabel('Percent of units in events');
xticklabels(presentName);
xlabel('Region');
ylim([-5 105]);
saveas(gcf,[specPath 'Plots\' saveSchem '.PercEvtBox.png']);

if ~showPlt
    close all
end

%% [2.3] # spikes each unit fires in each event
barProp.nSpkEach = unitBar.nSpkEach;
% Plot average # spikes by each unit by event (sum across events)
avgSpkUn = NaN(size(unitBar.nSpkEach,1),1);
for i = 1:size(unitBar.nSpkEach,1)
    avgSpkUn(i) = sum(unitBar.nSpkEach(i,:));
    tempCnt = length(find(unitBar.nSpkEach(i,:)~=0));
    if tempCnt
        avgSpkUn(i) = avgSpkUn(i)/tempCnt;
    else
        avgSpkUn(i) = 0;
    end
end
figure('Position', get(0, 'Screensize')); 
plot(avgSpkUn);
xlabel('Unit #');
ylabel('Average # of Spikes');
title('Average # of spikes during events when the unit is active');
saveas(gcf,[specPath 'Plots\' saveSchem '.avgSpkUn.png']);
barProp.avgSpkUn = avgSpkUn;
% maybe add region and type identification here so we know that the
% outliers are***********************

% Plot average # of spikes per event by each event (sum across units)
avgSpkEvt = NaN(size(unitBar.nSpkEach,2),1);
for i = 1:size(unitBar.nSpkEach,2)
    avgSpkEvt(i) = sum(unitBar.nSpkEach(:,i));
    tempCnt = length(find(unitBar.nSpkEach(:,i)~=0));
    if tempCnt
        avgSpkEvt(i) = avgSpkEvt(i)/tempCnt;
    else
        avgSpkEvt(i) = 0;
    end
end
figure('Position', get(0, 'Screensize')); 
plot(avgSpkEvt);
xlabel('Event #');
ylabel('Average # of Spikes');
xlim([1 size(unitBar.nSpkEach,2)]);
title('Average # of spikes per event from active units');
saveas(gcf,[specPath 'Plots\' saveSchem '.avgSpkEvt.png']);
barProp.avgSpkEvt = avgSpkEvt;

% Grouped bar plot per region (want to do this, maybe with python?)
regBox = cell(length(presentIND),1);
boxC = cell(length(presentIND),1);
for i = 1:length(presentIND)
    [keep] = includeType(regID,presentIND(i)+1);
    inds = find(keep==1);
    temp = [];
    for j = 1:length(inds)
        temp = [temp; avgSpkUn(inds(j))];
    end
    regBox{i} = temp;
    boxC{i} = (i)*ones(length(temp),1);
end
boxx = [];
boxg = [];
for i = 1:length(regBox)
    boxx = [boxx; regBox{i}];
    boxg = [boxg; boxC{i}];
end
figure('Position', get(0, 'Screensize'));  
boxplot(boxx, boxg);
xticklabels(presentName);
title('Average number of spikes per event per region');
ylabel('# Spikes');
xlabel('Region');
saveas(gcf,[specPath 'Plots\' saveSchem '.avgSpkBox.png']);

if ~showPlt
    close all
end

%% [2.4] FR per unit (per event)
barProp.FReach = unitBar.FReach;
% Plot average FR by each unit by event (sum across events)
avgFRUn = NaN(size(unitBar.FReach,1),1);
for i = 1:size(unitBar.FReach,1)
    avgFRUn(i) = sum(unitBar.FReach(i,:));
    tempCnt = length(find(unitBar.FReach(i,:)~=0));
    if tempCnt
        avgFRUn(i) = avgFRUn(i)/tempCnt;
    else
        avgFRUn(i) = 0;
    end
end
figure('Position', get(0, 'Screensize')); 
plot(avgFRUn);
xlabel('Unit #');
ylabel('Average Firing Rate');
title('Average Firing Rate during events when the unit is active');
saveas(gcf,[specPath 'Plots\' saveSchem '.avgFRUn.png']);
barProp.avgFRUn = avgFRUn;

% Plot average # of spikes per event by each event (sum across units)
avgFREvt = NaN(size(unitBar.FReach,2),1);
for i = 1:size(unitBar.FReach,2)
    avgFREvt(i) = sum(unitBar.FReach(:,i));
    tempCnt = length(find(unitBar.FReach(:,i)~=0));
    if tempCnt
        avgFREvt(i) = avgFREvt(i)/tempCnt;
    else
        avgFREvt(i) = 0;
    end
end
figure('Position', get(0, 'Screensize')); 
plot(avgFREvt);
xlabel('Event #');
ylabel('Average Firing Rate');
xlim([1 size(unitBar.FReach,2)]);
title('Average Firing Rate per event from active units');
saveas(gcf,[specPath 'Plots\' saveSchem '.avgFREvt.png']);
barProp.avgFREvt = avgFREvt;

% Grouped bar plot per region (want to do this, maybe with python?)
regBox = cell(length(presentIND),1);
boxC = cell(length(presentIND),1);
for i = 1:length(presentIND)
    [keep] = includeType(regID,presentIND(i)+1);
    inds = find(keep==1);
    temp = [];
    for j = 1:length(inds)
        temp = [temp; avgFRUn(inds(j))];
    end
    regBox{i} = temp;
    boxC{i} = (i)*ones(length(temp),1);
end
boxx = [];
boxg = [];
for i = 1:length(regBox)
    boxx = [boxx; regBox{i}];
    boxg = [boxg; boxC{i}];
end
figure('Position', get(0, 'Screensize'));  
boxplot(boxx, boxg);
xticklabels(presentName);
title('Average Firing Rate per event per region');
ylabel('# Spikes');
xlabel('Region');
saveas(gcf,[specPath 'Plots\' saveSchem '.avgFRBox.png']);

if ~showPlt
    close all
end


%% [2.5] Population SAVE
save(strcat(specPath,'props.mat'),'barProp', '-v7.3');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% [3.5] Participation per cell (extra) - similar to 3.2
% Plot (unit# by fraction of participation)
% [particip, pI] = sort(unitMet.particip, 'descend');
% for i = 1:length(particip)
%     regIDs(i) = regID(pI(i));
%     cellIDs(i) = cellID(pI(i));
% end 
% 
% plotSubs(1:length(particip),particip, regIDs, cellIDs);
% ylim([0 1]);
% title('Fraction of participation per unit');
% ylabel('Fraction of parcipation (max 1)');
% xlabel('Sorted Unit # (descending participation)');
% hold off

%% [APPENDIX]
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [nx, ny] = excludeType(x, y, ID, ex)
nx = []; ny = []; nc = 1; flg = 0;
for i = 1:length(ID)
    for j=1:length(ex)
        if ID(i) == ex(j)
            flg = flg+1;
        end
    end
    if ~flg
        nx(nc) = x(i);
        ny(nc) = y(i);
        nc = nc+1;
    end
    flg = 0;
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [keep] = includeType(ID, inc)
keep = zeros(length(ID),1); flg = 0;
for i = 1:length(ID)
    for j=1:length(inc)
        if ID(i) == inc(j)
            flg = flg+1;
        end
    end
    keep(i) = logical(flg);
    flg = 0;
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotSubs(x, y, regID, cellID, line)
if nargin < 5
    line = '';
end
color = ['m';'b';'g';'r';'c';'k'];
type = ['*';'^';'o'];

figure('Position', get(0, 'Screensize'))
hold on
for i = 1:length(cellID)
    if isempty(line)
        plot(x(i), y(i), cat(1,color(regID(i)),type(cellID(i))));
    else
        plot(x, y(:,i), cat(1,color(regID(i)),type(cellID(i)),line));
    end
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
legend(h, 'CA1','CA2','CA3','CTX','DG','Pyramidal','Interneuron','Unknown', 'Location', 'eastoutside');
end