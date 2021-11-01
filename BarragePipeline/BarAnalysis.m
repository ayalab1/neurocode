%% Barrage Analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is currently for use in singular sessions, but can easily be
% expanded using the outline from MassSessRun.m
%
%
% TO DO LIST:
% -- Clean up and sort functions
%    -- BoxPlot
%    -- 
% -- Integrate python to create nicer plots
% -- Implement sleep states***
% -- Implement ripple comparisons

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function BarAnalysis(ses)
%% [STEP 0] Settings
close all
savePath = ('Z:\home\Lindsay\Barrage\');
load('Z:\home\Lindsay\Barrage\combinedPaths.mat');
cd(ses);
basepath = pwd;
basename = basenameFromBasepath(basepath);
animName = animalFromBasepath(basepath);
% Load in session information
load(['Barrage_Files\' basename '.allpyr.cellinfo.mat']);
load([basename '.cell_metrics.cellinfo.mat']);
load([basename '.session.mat']);
showPlt = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% [STEP 1]
%%% BURST PROPERTIES [SINGLE CELL, 10ms] %%%
specPath = strcat(savePath,'Cell\',animName,'\');
plotPath = strcat(savePath,'Cell\','Plots\',animName,'\');
if ~exist(specPath)
    mkdir(specPath);
end
if ~exist(plotPath)
    mkdir(plotPath);
end
saveSchem = basename;
if exist(strcat(specPath,'props.mat'))
    load([specPath 'props.mat']);
end
[regID, cellID] = getSubs(cell_metrics);
cellProp.regID = regID;
cellProp.cellID = cellID;

%% [1.1] ISI per unit
% Distribution
[ISI,ISIc,t] = ISIGrams(spikes.times, spikes.UID, 1/1000, 1000);
cellProp.ISI = ISI;
cellProp.ISIhistT = t;
cellProp.ISIhist = ISIc;
ISIavg = NaN(length(ISI),1);
for i = 1:length(ISI)
    ISIavg(i) = mean(ISI{i});
end
cellProp.ISIavg = ISIavg;

% Plotting
figure('Position', get(0, 'Screensize'));
plot(t*1e3,sum(ISIc,2));
xlabel('Log of ISI (ms)');
ylabel('Count');
xlim([0 1000]);
title('Distribution of ISIs in Log scale');
set(gca, 'XScale', 'log')
saveas(gcf,[plotPath saveSchem '.ISIdistLog.png']);

cumMet.ISIx = t*1e3;
cumMet.ISIy = sum(ISIc,2);

% figure('Position', get(0, 'Screensize'));
% plot(spikes.UID, ISIavg, '.');
% xlabel('Unit #');
% ylabel('Log of avg ISI (ms)');
% title('Average ISI per cell');
% set(gca, 'YScale', 'log')
% saveas(gcf,[plotPath saveSchem '.ISIavg.png']);

% By type plotting
plotSubs(spikes.UID, ISIavg*1e3, regID, cellID, spikes.UID);
xlabel('Unit #');
ylabel('Log of avg ISI (ms)');
title('Average ISI per cell by type');
set(gca, 'YScale', 'log')
hold off
saveas(gcf,[plotPath saveSchem '.ISIavgSort.png']);

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
saveas(gcf,[plotPath saveSchem '.ISIdist.png']);

% Per cell with type designations
plotSubs(t*1000, ISIc, regID, cellID, spikes.UID, '-');
xlim([t(1)*1000 t(end)*1000]);
title('Histogram of ISIs by cell type and region');
ylabel('Count');
xlabel('ISI (ms)');
saveas(gcf,[plotPath saveSchem '.ISIdistSort.png']);

% Per cell with type designations, shortened
plotSubs(t*1000, ISIc, regID, cellID, spikes.UID, '-');
xlim([t(1)*1000 20]);
title('Histogram of ISIs by cell type and region');
ylabel('Count');
xlabel('ISI (ms)');
saveas(gcf,[plotPath saveSchem '.ISIdistSortShort.png']);

% Summary box plot per region
check = ["CA1" "CA2" "CA3" "CTX" "DG"]; %this is embedded in getSubs()
regBox = cell(length(check),1);
boxC = cell(length(check),1);
for i = 2:(length(check)+1) %shift because regID=1 is unknown
    if find(regID==i, 1) %if region is present
        [keep] = includeType(regID,i); %pull from ALL cells, need to check they're actually in our spike subset
        % now we have a logical index of which regions fit our region from the FULL recording
        % We only want to keep those indices which are found within spikes.UID
        indsTemp = find(keep==1);
        [~,~,inds] = intersect(indsTemp,spikes.UID);
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
% boxplot(boxx, boxg,'symbol', ''); %use symbol '' to suppress outliers
boxplot(boxx, boxg);
useName = unique(boxg);
nameInd = ismember(1:length(check),useName);
presentIND = find(nameInd==1);
presentName = convertStringsToChars(check(nameInd));
xticklabels(presentName);
title('ISI per region, outliers cut off');
ylabel('ISI (s)');
xlabel('Region');
ylim([-1 7]);
saveas(gcf,[plotPath saveSchem '.ISIboxReg.png']);

cumMet.boxxISI = boxx;
cumMet.boxgISI = boxg;

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

% figure('Position', get(0, 'Screensize'));
% plot(spikes.UID, avgFR, '.');
% title('Average Firing Rate');
% ylabel('Firing Rate (Hz)');
% xlabel('Unit Number');
% saveas(gcf,[plotPath saveSchem '.avgFR.png']);

plotSubs(spikes.UID, avgFR, regID, cellID, spikes.UID);
title('Average Firing Rate by Type');
ylabel('Firing Rate (Hz)');
xlabel('Unit Number');
saveas(gcf,[plotPath saveSchem '.avgFRSort.png']);

if ~showPlt
    close all
end

%% [1.4] Burst Index per cell
binssum = 6; %6ms
burstIndex = sum(ISIc(1:binssum,:),1)./sum(ISIc,1);
cellProp.burstIndex = burstIndex;

plotSubs(spikes.UID, burstIndex, regID, cellID, spikes.UID);
title('Burst Index by Cell Type');
ylabel('Burst Index');
xlabel('Unit Number');
saveas(gcf,[plotPath saveSchem '.burstIndSort.png']);

% Box plot
regBox = cell(length(presentIND),1);
boxB = cell(length(presentIND),1);
for i = 1:length(presentIND)
    [keep] = includeType(regID,presentIND(i)+1);
    indsTemp = find(keep==1);
    [~,~,inds] = intersect(indsTemp,spikes.UID);
    temp = [];
    for j = 1:length(inds)
        temp = [temp; burstIndex(inds(j))];
    end
    regBox{i} = temp;
    boxB{i} = (i)*ones(length(temp),1);
end
boxx = [];
boxg = [];
for i = 1:length(regBox)
    boxx = [boxx; regBox{i}];
    boxg = [boxg; boxB{i}];
end
figure('Position', get(0, 'Screensize'));  
boxplot(boxx, boxg);
xticklabels(presentName);
title('Burstiness Index per region');
ylabel('Burst Index');
xlabel('Region');
saveas(gcf,[plotPath saveSchem '.BurstBox.png']);

cumMet.boxxBurst = boxx;
cumMet.boxgBurst = boxg;


if ~showPlt
    close all
end

%% [1.5] Average burst length by region
burstEvts = cell(size(ISI,1),1);
avgBurstSz = NaN(size(ISI,1),1);
for i = 1:size(ISI,1)
    evtstart = spikes.times{i};
    evtstart = evtstart(1:length(evtstart)-1);
    evtstop = spikes.times{i};
    evtstop = evtstop(2:end);
    evtdur = evtstop-evtstart;
    evtpeak = evtstart + (evtdur/2);
    evtamp = zeros(length(evtstart),1);
    flagConc = (evtdur <= (6/1000));
    [start, stop, dur,~] = CatCon(evtstart,evtstop,evtpeak,evtamp,flagConc);
    burstEvts{i} = [start; stop];
    avgBurstSz(i) = sum(dur)/length(dur);
end

% Plotting
% Box plot
regBox = cell(length(presentIND),1);
boxB = cell(length(presentIND),1);
for i = 1:length(presentIND)
    [keep] = includeType(regID,presentIND(i)+1);
    indsTemp = find(keep==1);
    [~,~,inds] = intersect(indsTemp,spikes.UID);
    temp = [];
    for j = 1:length(inds)
        temp = [temp avgBurstSz(inds(j))];
    end
    regBox{i} = temp';
    boxB{i} = (i)*ones(length(temp),1);
end
boxx = [];
boxg = [];
for i = 1:length(regBox)
    boxx = [boxx; regBox{i}];
    boxg = [boxg; boxB{i}];
end
figure('Position', get(0, 'Screensize'));  
boxplot(boxx, boxg);
xticklabels(presentName);
title('Avg number of spikes per burst per region');
ylabel('Burst Length (# spikes)');
xlabel('Region');
saveas(gcf,[plotPath saveSchem '.BurstSzBox.png']);

cellProp.burstSz = avgBurstSz;
cellProp.burstEvts = burstEvts;
cumMet.boxxBurstSz = boxx;
cumMet.boxgBurstSz = boxg;

if ~showPlt
    close all
end

%% [1.6] Cell SAVE
save(strcat(specPath,basename,'.props.mat'),'cellProp', '-v7.3');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% [STEP 2]
%%% BARRAGE DETECTION PROPERTIES %%%
% Need to keep event duration in mind
specPath = strcat(savePath,'Population\',animName,'\');
plotPath = strcat(savePath,'Population\','Plots\',animName,'\');
if ~exist(specPath)
    mkdir(specPath);
end
if ~exist(plotPath)
    mkdir(plotPath);
end
load([basepath '\Barrage_Files\' basename '.HSEfutEVT.mat']);
load([basepath '\Barrage_Files\' basename '.HSEmetrics.mat']);
detection = find(HSEmetrics.Notes == "All pyr",1,'last');
saveSchem = strcat(basename,'.',num2str(detection));
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

% [sortCell] = sort(unitBar.nCellsEvent,'descend');
% figure('Position', get(0, 'Screensize'));
% bar(1:length(sortCell),sortCell);
% title('Number of units per event');
% xlabel('Event # (sorted descending)');
% ylabel('Number of Units');
% saveas(gcf,[plotPath saveSchem '.untPerEvt.png']);

[~,edges,bins] = histcounts(spkEventTimes.EventDuration);
durCount = zeros(max(bins),1);
binct = zeros(max(bins),1);
durX = [];
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
saveas(gcf,[plotPath saveSchem '.untPerEvtDur.png']);

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
saveas(gcf,[plotPath saveSchem '.percEvtReg.png']);

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
    saveas(gcf,[plotPath saveSchem '.percEvtRegCA2sort.png']);
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
saveas(gcf,[plotPath saveSchem '.untPercentDur.png']);

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
saveas(gcf,[plotPath saveSchem '.PercEvtBox.png']);

cumMet.boxxPerc = boxx;
cumMet.boxgPerc = boxg;

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
plot(spikes.UID, avgSpkUn);
xlabel('Unit #');
ylabel('Average # of Spikes');
title('Average # of spikes during events when the unit is active');
saveas(gcf,[plotPath saveSchem '.avgSpkUn.png']);
barProp.avgSpkUn = avgSpkUn;

% By region
plotSubs(spikes.UID, avgSpkUn, regID, cellID,spikes.UID);
xlabel('Unit #');
ylabel('Average # of Spikes');
title('Average # of spikes during events when the unit is active');
saveas(gcf,[plotPath saveSchem '.avgSpkUnSort.png']);

% Plot average # of spikes per event by each event (sum across units)
avgSpkEvt = NaN(size(unitBar.nSpkEach,2),1);
for i = 1:size(unitBar.nSpkEach,2)
    avgSpkEvt(i) = sum(unitBar.nSpkEach(:,i));
    tempCnt(i) = length(find(unitBar.nSpkEach(:,i)~=0));
    if tempCnt(i)
        avgSpkEvt(i) = avgSpkEvt(i)/tempCnt(i);
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
saveas(gcf,[plotPath saveSchem '.avgSpkEvt.png']);
barProp.avgSpkEvt = avgSpkEvt;

% Sort by event duration
[~,edges,bins] = histcounts(spkEventTimes.EventDuration);
durCount = zeros(max(bins),1);
durSum = zeros(max(bins),1);
binct = zeros(max(bins),1);
for i = 1:length(bins)
    durCount(bins(i)) = durCount(bins(i))+sum(unitBar.nSpkEach(:,i));
    durSum(bins(i)) = durSum(bins(i)) + tempCnt(i);
end
for i = 1:length(durSum)
    if durSum(i) == 0
        durSum(i) = 1;
    end
end
durAvg = durCount./durSum;
for i = 1:length(edges)-1
    durX(i) = edges(i) + ((edges(i+1)-edges(i))/2);
end
figure('Position', get(0, 'Screensize')); 
bar(durX,durAvg);
xlabel('Event Duration (s)');
ylabel('Average # of Spikes');
title('Average # of spikes per event from active units by event duration');
saveas(gcf,[plotPath saveSchem '.avgSpkEvtDur.png']);

% Grouped bar plot per region (want to do this, maybe with python?)
regBox = cell(length(presentIND),1);
boxC = cell(length(presentIND),1);
for i = 1:length(presentIND)
    [keep] = includeType(regID,presentIND(i)+1);
    indsTemp = find(keep==1);
    [~,~,inds] = intersect(indsTemp,spikes.UID);
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
saveas(gcf,[plotPath saveSchem '.avgSpkBox.png']);

cumMet.boxxAvgSpk = boxx;
cumMet.boxgAvgSpk = boxg;

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
plot(spikes.UID,avgFRUn);
xlabel('Unit #');
ylabel('Average Firing Rate');
title('Average Firing Rate during events when the unit is active');
saveas(gcf,[plotPath saveSchem '.avgFRUn.png']);
barProp.avgFRUn = avgFRUn;

% Plot average # of spikes per event by each event (sum across units)
avgFREvt = NaN(size(unitBar.FReach,2),1);
for i = 1:size(unitBar.FReach,2)
    avgFREvt(i) = sum(unitBar.FReach(:,i));
    tempCnt(i) = length(find(unitBar.FReach(:,i)~=0));
    if tempCnt(i)
        avgFREvt(i) = avgFREvt(i)/tempCnt(i);
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
saveas(gcf,[plotPath saveSchem '.avgFREvt.png']);
barProp.avgFREvt = avgFREvt;

% Sort by event duration
[~,edges,bins] = histcounts(spkEventTimes.EventDuration);
durCount = zeros(max(bins),1);
durSum = zeros(max(bins),1);
binct = zeros(max(bins),1);
for i = 1:length(bins)
    durCount(bins(i)) = durCount(bins(i))+sum(unitBar.FReach(:,i));
    durSum(bins(i)) = durSum(bins(i)) + tempCnt(i);
end
for i = 1:length(durSum)
    if durSum(i) == 0
        durSum(i) = 1;
    end
end
durAvg = durCount./durSum;
for i = 1:length(edges)-1
    durX(i) = edges(i) + ((edges(i+1)-edges(i))/2);
end
figure('Position', get(0, 'Screensize')); 
bar(durX,durAvg);
xlabel('Event Duration (s)');
ylabel('Average FR of Spikes');
title('Average FR of spikes per event from active units by event duration');
saveas(gcf,[plotPath saveSchem '.avgFREvtDur.png']);

% Grouped bar plot per region (want to do this, maybe with python?)
regBox = cell(length(presentIND),1);
boxC = cell(length(presentIND),1);
for i = 1:length(presentIND)
    [keep] = includeType(regID,presentIND(i)+1);
    indsTemp = find(keep==1);
    [~,~,inds] = intersect(indsTemp,spikes.UID);
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
ylabel('Firing Rate (Hz)');
xlabel('Region');
saveas(gcf,[plotPath saveSchem '.avgFRBox.png']);

cumMet.boxxAvgFR = boxx;
cumMet.boxgAvgFR = boxg;

if ~showPlt
    close all
end

%% [2.5] Event Duration Characterization
bins = [0:50e-3:5]; %in seconds
evtDur = histc(spkEventTimes.EventDuration,bins);%distribute counts per time bins
% duration2plot=smooth(duration2plot,'sgolay',3);%smooth if you want
figure('Position', get(0, 'Screensize')); 
plot(bins, evtDur./sum(evtDur));%normalize per the total n of events
title('Event Duration normalized to number of events');
ylabel('Normalized Event Count');
xlabel('Event duration (s)');
saveas(gcf,[plotPath saveSchem '.evtDur.png']);

if ~showPlt
    close all
end

%% [2.6] Population SAVE
save(strcat(specPath,saveSchem,'.props.mat'),'barProp', '-v7.3');
save(strcat('Z:\home\Lindsay\Barrage\CumMet\',animName,'.',basename,'.cumMet.mat'),'cumMet', '-v7.3');
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

end