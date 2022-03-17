%% Barrage Analysis (remade)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Updated with unifying save structures to minimize errors
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function BarAnalysis(ses)
if nargin <1
    warning('No session path input, defaulting to current folder');
    ses = pwd;
end
%% [STEP 0] Settings
% warning('off','all');
close all
savePath = strcat(ses, '\Barrage_Files\'); savePath = convertCharsToStrings(savePath);
load('Z:\home\Lindsay\Barrage\combinedPaths.mat');
cd(ses);
basepath = pwd;
basename = basenameFromBasepath(basepath);
animName = animalFromBasepath(basepath);
% Load in session information
load(['Barrage_Files\' basename '.allpyr.cellinfo.mat']);
load(['Barrage_Files\' basename '.useSpk.UIDkeep.mat']);
load([basename '.cell_metrics.cellinfo.mat']);
load([basename '.session.mat']);
showPlt = 0;
reRun = 1; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% [STEP 1]
%%% SINGLE CELL PROPERTIES
specPath = strcat(savePath,'Cell\');
if ~exist(specPath)
    mkdir(specPath);
end
saveSchem = basename;
if exist(strcat(specPath,'props.mat'))
    load([specPath 'props.mat']);
end

[regID, modID, regKey, modKey] = getSubs(cell_metrics);
cellProp.regID = regID;
cellProp.modID = modID;
cellProp.regKey = regKey;
cellProp.modKey = modKey;
i = 1;
for r = 1:size(regKey,2)
    for m = 1:size(modKey,2)
        comKey(1,i) = strcat(regKey(1,r),'.',modKey(1,m));
        comKey(2,i) = i;
        i = i+1;
    end
end
cellProp.comKey = comKey;

cumMet.regID = regID;
cumMet.modID = modID;
cumMet.regKey = regKey;
cumMet.modKey = modKey;
cumMet.comKey = comKey;

%% [1.0] Sort cells by subtype
% Put list of UIDs into a structure size nReg x mMod
UIDfull = cell(size(regKey,2),size(modKey,2));
UIDsort = cell(size(regKey,2),size(modKey,2));
for r = 1:size(regKey,2)
    keepR = find(regID == r); %find unit #s in our region of interest
    for m = 1:size(modKey,2)
        keepM = find(modID == m); %find unit #s in our mod type of interest
        useM = ismember(keepM,keepR);
        useU = ismember(spikes.UID,keepM(useM));
        UIDfull{r,m} = spikes.UID(useU); %referenced to full set of spikes
        UIDsort{r,m} = find(useU == 1); %referenced to pulled spike set
    end
end

%% [1.1] ISI per unit
% Get distribution
[ISI,ISIc,t] = ISIGrams(spikes.times, spikes.UID, 1/1000, 1000);
cellProp.ISI = ISI;
cellProp.ISIt = t;
cellProp.ISIcnt = ISIc;
ISIavg = NaN(length(ISI),1);
for i = 1:length(ISI)
    ISIavg(i) = mean(ISI{i}); %average per cell
end
cellProp.ISIavg = ISIavg;

% Plotting
figure('Position', get(0, 'Screensize'));
plot(t*1e3,sum(ISIc,2));
xlabel('Log of ISI (ms)');
ylabel('Count');
xlim([0 500]);
title('Distribution of ISIs in Log scale');
set(gca, 'XScale', 'log')
saveas(gcf,strcat(specPath, saveSchem, '.ISIdistLog.png'));

cumMet.ISIx = t*1e3;
cumMet.ISIy = sum(ISIc,2);
cumMet.ISIc = ISIc;

% Plot by cells together by subtype
plotSubs(spikes.UID, ISIavg*1e3, regID, modID, spikes.UID, UIDkeep);
xlabel('Unit #');
ylabel('Log of avg ISI (ms)');
title('Average ISI per cell by type');
set(gca, 'YScale', 'log')
hold off
saveas(gcf,strcat(specPath, saveSchem, '.ISIavgType.png'));

if ~showPlt
    close all
end

%% [1.2] # spikes in 10 ms windows
% Per cell with type designations
plotSubs(t*1000, ISIc, regID, modID, spikes.UID, UIDkeep, '-');
xlim([0 500]);
title('Distribution of ISIs by cell type and region');
ylabel('Count');
xlabel('ISI (ms)');
saveas(gcf,strcat(specPath, saveSchem, '.ISIdistType.png'));

% Per cell with type designations, shortened
plotSubs(t*1000, ISIc, regID, modID, spikes.UID, UIDkeep, '-');
xlim([0 20]);
title('Shortened distribution of ISIs by cell type and region');
ylabel('Count');
xlabel('ISI (ms)');
saveas(gcf,strcat(specPath, saveSchem, '.ISIdistTypeShort.png'));

% Summary box plot per region
boxx = [];
boxg = [];
for r = 1:size(regKey,2)
    useR = cat(2,UIDsort{r,:});
    useU = ismember(spikes.UID,useR);
    Uind = find(useU==1);
    for g = 1:length(Uind)
        boxx = [boxx; ISI{spikes.UID(Uind(g))}];
        boxg = [boxg; r*ones(length(ISI{spikes.UID(Uind(g))}),1)];
    end
end
    
figure('Position', get(0, 'Screensize'));  
boxplot(boxx, boxg);
regInd = ismember(1:size(regKey,2),unique(boxg));
xticklabels(regKey(1,regInd));
title('ISI per region, outliers cut off');
ylabel('ISI (s)');
xlabel('Region');
ylim([-1 7]);
saveas(gcf,strcat(specPath, saveSchem, '.ISIboxReg.png'));    

cumMet.boxxISI = boxx;
cumMet.boxgISI = boxg;

% Summary box plot per region and per mod type
c = 1;
boxx = [];
boxg = [];
for r = 1:size(regKey,2)
    for m = 1:size(modKey,2)
        useR = UIDsort{r,m};
        useU = ismember(spikes.UID,useR);
        Uind = find(useU==1);
        for g = 1:length(Uind)
            boxx = [boxx; ISI{spikes.UID(Uind(g))}];
            boxg = [boxg; c*ones(length(ISI{spikes.UID(Uind(g))}),1)];
        end
        c = c+1;
    end
end
figure('Position', get(0, 'Screensize'));  
boxplot(boxx, boxg);
comInd = ismember(1:size(comKey,2),unique(boxg));
xticklabels(comKey(1,comInd));
title('ISI per region and mod type, outliers cut off');
ylabel('ISI (s)');
xlabel('Region');
ylim([-1 7]);
saveas(gcf,strcat(specPath, saveSchem, '.ISIboxAllType.png'));  

if ~showPlt
    close all
end

%% [1.3] FR per unit
FR = cell(length(spikes.times),1);
avgFR = NaN(length(spikes.times),1);
tSmooth = 0.02;
binSz = 0.005;
for i = 1:length(spikes.times)
    spkHist = spkRtHist(sort(spikes.times{i}),'tSmooth',tSmooth,'binSz',binSz,'ifz', false);
    FR{i} = spkHist/binSz; %get rate from our count/bin size
    avgFR(i) = sum(FR{i})/length(FR{i});
end
cellProp.FR = FR;
cellProp.avgFR = avgFR;

% FR by type
plotSubs(spikes.UID, avgFR, regID, modID, spikes.UID, UIDkeep);
title('Average Firing Rate by Type');
ylabel('Firing Rate (Hz)');
xlabel('Unit Number');
saveas(gcf,strcat(specPath, saveSchem, '.avgFRtype.png'));

if ~showPlt
    close all
end

%% [1.4] Burst Index per cell
binssum = 6; %6ms
burstIndex = sum(ISIc(1:binssum,:),1)./sum(ISIc,1);
cellProp.burstIndex = burstIndex;

% Plot 
plotSubs(spikes.UID, burstIndex, regID, modID, spikes.UID, UIDkeep);
title('Burst Index by Cell Type');
ylabel('Burst Index');
xlabel('Unit Number');
saveas(gcf,strcat(specPath, saveSchem, '.burstIndType.png'));

% Plot per region (box plot)
boxx = [];
boxg = [];
for r = 1:size(regKey,2)
    useR = cat(2,UIDsort{r,:});
    useU = ismember(spikes.UID,useR);
    Uind = find(useU==1);
    for g = 1:length(Uind)
        boxx = [boxx; burstIndex(spikes.UID(Uind(g)))];
        boxg = [boxg; r*ones(length(burstIndex(spikes.UID(Uind(g)))),1)];
    end
end

figure('Position', get(0, 'Screensize'));  
boxplot(boxx, boxg);
regInd = ismember(1:size(regKey,2),unique(boxg));
xticklabels(regKey(1,regInd));
title('Burstiness Index per region');
ylabel('Burst Index');
xlabel('Region');
saveas(gcf,strcat(specPath, saveSchem, '.burstBox.png'));

cumMet.boxxBurst = boxx;
cumMet.boxgBurst = boxg;


if ~showPlt
    close all
end

%% [1.5] Burst Size metrics
burstEvts = cell(size(ISI,1),1);
burstSz = cell(size(ISI,1),1);
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
    [start, stop, ~,~,numCat] = CatCon(evtstart,evtstop,evtpeak,evtamp,flagConc);
    keep_burst = find(numCat > 1);
    burstEvts{i} = [start(keep_burst); stop(keep_burst)];
    
    avgBurstSz(i) = sum(numCat(keep_burst))/length(numCat(keep_burst));
    burstSz{i} = numCat(keep_burst);
end

% Plot by region
boxx = [];
boxg = [];
for r = 1:size(regKey,2)
    useR = cat(2,UIDsort{r,:});
    useU = ismember(spikes.UID,useR);
    Uind = find(useU==1);
    for g = 1:length(Uind)
        boxx = [boxx; avgBurstSz(spikes.UID(Uind(g)))];
        boxg = [boxg; r*ones(length(burstIndex(spikes.UID(Uind(g)))),1)];
    end
end
figure('Position', get(0, 'Screensize'));  
boxplot(boxx, boxg);
regInd = ismember(1:size(regKey,2),unique(boxg));
xticklabels(regKey(1,regInd));
title('Average number of spikes per burst per region');
ylabel('Burst Length (# spikes)');
xlabel('Region');
saveas(gcf,strcat(specPath, saveSchem, '.burstSzBox.png'));

cellProp.avgBurstSz = avgBurstSz;
cellProp.burstEvts = burstEvts;
cellProp.burstSz = burstSz;
cumMet.boxxAvgBurstSz = boxx;
cumMet.boxgAvgBurstSz = boxg;

% Plot by region and type
c = 1;
boxx = [];
boxg = [];
resCell = cell(size(regKey,2),size(modKey,2));
for r = 1:size(regKey,2)
    for m = 1:size(modKey,2)
        useR = UIDsort{r,m};
        useU = ismember(spikes.UID,useR);
        Uind = find(useU==1);
        resTemp = [];
        for g = 1:length(Uind)
            resTemp = [resTemp; burstSz(spikes.UID(Uind(g)))];
            boxx = [boxx; avgBurstSz(spikes.UID(Uind(g)))];
            boxg = [boxg; c*ones(length(avgBurstSz(spikes.UID(Uind(g)))),1)];
        end
        resCell{r,m} = resTemp;
        c = c+1;
    end
end
figure('Position', get(0, 'Screensize'));  
boxplot(boxx, boxg);
comInd = ismember(1:size(comKey,2),unique(boxg));
xticklabels(comKey(1,comInd));
title('Average number of spikes per burst per region per type');
ylabel('Burst Length (# spikes)');
xlabel('Region');
saveas(gcf,strcat(specPath, saveSchem, '.burstSzBoxAllType.png'));

% Burst size histogram per region and per type
figure('Position', get(0, 'Screensize'));
hold on;
title('Histograms of burst length per region types');
c = 2;
i = 1;
for r = 1:(size(resCell,1))
    tempRes = [];
    if ~isempty(resCell{r,c})
        tempRes = cat(2,resCell{r,c}{:});
    end
    subplot(size(resCell,1),2,i); histogram(tempRes,[2:2:10]);hold on; title(strcat(regKey(1,r),' ',modKey(1,c)));
    if r == size(resCell,1)
        ylabel('Burst count');
        xlabel('Burst length');
    end
    i=i+1;
    tempRes = [];
    if ~isempty(resCell{r,c+1})
        tempRes = cat(2,resCell{r,c+1}{:});
    end
    subplot(size(resCell,1),2,i); histogram(tempRes,[2:2:10]);hold on; title(strcat(regKey(1,r),' ',modKey(1,c+1)));
    i=i+1;
end
hold off;
saveas(gcf,strcat(specPath, saveSchem, '.burstSzHist.png'));

cumMet.burstCnt = resCell;

% Number of bursts per region and per type with > 4 spikes
figure('Position', get(0, 'Screensize'));
hold on;
title('Number of bursts per region and type with more than 4 spikes');
c = 2;
i = 1;
naming = []; barPx = []; barPy = [];
for r = 1:size(resCell,1)
    tempRes = [];
    if ~isempty(resCell{r,c})
        tempRes = cat(2,resCell{r,c}{:});
    end
    barPy = [barPy (sum(tempRes>4)/length(tempRes))];
    barPx = [barPx i];
    naming = [naming; strcat(regKey(1,r),'.',modKey(1,c))];
    i = i+1;
    tempRes = [];
    if ~isempty(resCell{r,c+1})
        tempRes = cat(2,resCell{r,c+1}{:});
    end
    barPy = [barPy (sum(tempRes>4)/length(tempRes))];
    barPx = [barPx i];
    naming = [naming; strcat(regKey(1,r),'.',modKey(1,c+1))];
    i = i+1;
end
bar(barPx, barPy);
xticks(barPx);
xticklabels(naming);
ylabel('Normalized burst length count > 4 spikes');
xlabel('Region and Ripple Modulation Type');
ylim([0 0.02]);
hold off;
saveas(gcf,strcat(specPath, saveSchem, '.burstSzThresh.png'));

% Spread of burst lengths (imagesc)
figure('Position', get(0, 'Screensize'));
hold on;
title('Spread of burst sizes per cell type');
c = 2;
i = 1;
for r = 1:size(resCell,1)
    tempRes = [];
    if ~isempty(resCell{r,c})
        for j = 1:length(resCell{r,c})
            tempRes(j,:) = hist(resCell{r,c}{j},[2:2:10])/sum(hist(resCell{r,c}{j},[2:2:10]));
        end
    end
    subplot(size(resCell,1),2,i); imagesc(tempRes,[0 0.05]);hold on; title(strcat(regKey(1,r),' ',modKey(1,c)));
    xticks([1 2 3 4 5]);
    xticklabels({num2str(2) num2str(4) num2str(6) num2str(8) num2str(10)});
    if r == size(resCell,1)
        ylabel('Unit');
        xlabel('Burst length');
    end
    i = i+1;
    tempRes = [];
    if ~isempty(resCell{r,c+1})
        for j = 1:length(resCell{r,c+1})
            tempRes(j,:) = hist(resCell{r,c+1}{j},[2:2:10])/sum(hist(resCell{r,c+1}{j},[2:2:10]));
        end
    end
    subplot(size(resCell,1),2,i); imagesc(tempRes,[0 0.05]);hold on; title(strcat(regKey(1,r),' ',modKey(1,c+1)));
    xticks([1 2 3 4 5]);
    xticklabels({num2str(2) num2str(4) num2str(6) num2str(8) num2str(10)});
    i = i+1;
end
hold off;
saveas(gcf,strcat(specPath, saveSchem, '.burstSzSC.png'));

if ~showPlt
    close all
end

%% [1.5] Cell SAVE
saveStructs(cellProp,specPath,basename);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% [STEP 2]
%%% "BARRAGE" DETECTION PROPERTIES %%%
specPath = strcat(savePath,'Population\');
if ~exist(specPath)
    mkdir(specPath);
end
load([basepath '\Barrage_Files\' basename '.HSE.mat']);
load([basepath '\Barrage_Files\' basename '.HSEmetrics.mat']);
load([basepath '\' basename '.SleepState.states.mat']);
saveSchem = strcat(basename,'.useSpk');
if exist(strcat(specPath,'props.mat'))
    load([specPath 'props.mat']);
end

if isfield(HSE, 'keep')
    HSE.timestamps = HSE.timestamps(HSE.keep,:);
    HSE.peaks = HSE.peaks(HSE.keep);
    HSE.amplitudes = HSE.amplitudes(HSE.keep);
    HSE.duration = HSE.duration(HSE.keep);
end
HSEnrem = eventIntervals(HSE,SleepState.ints.NREMstate,1);
%% [2.0] Initial runs
if reRun
    [spkEventTimes] = getRipSpikes('basepath',pwd,'events',HSEnrem.timestamps,'spikes',spikes,'savePath',savePath,'padding',0,'saveMat',false);
    save(strcat(specPath,saveSchem,'.spkEventTimes.mat'),'spkEventTimes', '-v7.3');
    [unitBar] = unitSWRmetrics(spkEventTimes,spikes);
    save(strcat(specPath,saveSchem,'.unitBar.mat'),'unitBar', '-v7.3');
else
    load(strcat(specPath,saveSchem,'.spkEventTimes.mat'));
    load(strcat(specPath,saveSchem,'.unitBar.mat'));
end

%% [2.1] Event Duration Characteristics
bins = [0:50e-3:5]; %in seconds
evtDur = histc(spkEventTimes.EventDuration,bins);%distribute counts per time bins
% duration2plot=smooth(duration2plot,'sgolay',3);%smooth if you want
figure('Position', get(0, 'Screensize')); 
plot(bins, evtDur./sum(evtDur));%normalize per the total n of events
title('Event Duration normalized to number of events');
ylabel('Normalized Event Count');
xlabel('Event duration (s)');
saveas(gcf,strcat(specPath, saveSchem, '.evtDur.png'));

if ~showPlt
    close all
end

cumMet.evtDur = spkEventTimes.EventDuration;

%% [2.2] # units per event
barProp.untPerEvt = unitBar.nCellsEvent;
barProp.untPerEvtID = unitBar.nCellsEventUID;
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
saveas(gcf,strcat(specPath, saveSchem, '.untPerEvtDur.png'));

if ~showPlt
    close all
end

%% [2.3] % units from different subregions active during events
percEvt = zeros(length(unitBar.nCellsEventUID),size(regKey,2));
for e = 1:length(unitBar.nCellsEventUID)
    regSum = 0;
    for r = 1:size(regKey,2)
        useR = cat(2, UIDsort{r,:});
        useU = ismember(unitBar.nCellsEventUID{e},useR);
        Uind = find(useU==1);
        percEvt(e,r) = length(Uind);
        regSum = regSum + length(Uind);
    end
    if regSum == 0; regSum = 1; end
    percEvt(e,:) = 100*percEvt(e,:)/regSum;
end
barProp.percEvt = percEvt;

figure('Position', get(0, 'Screensize'));
bar(percEvt,'stacked');
title('Percent of units from each region per event');
xlabel('Event #');
ylabel('Percent');
legend(regKey(1,:), 'Location', 'eastoutside');
ylim([0 100]);
saveas(gcf,strcat(specPath, saveSchem, '.percEvtReg.png'));
    % we have unitBar.nCellsEventUID - just needs to be sorted by UIDsort
        
% Plot sorted by increasing CA2 involvment
if ismember(3, regID)
   [regSortCA2, rI] = sort(percEvt(:,(regKey(2,:)==num2str(3))));
   for i = 1:length(regSortCA2)
       for j = 1:size(regKey,2)
           regSort(i,j) = percEvt(rI(i),j);
       end
   end
   figure('Position', get(0, 'Screensize'));
   bar(regSort,'stacked');
   title(['Percent of units from each region per event, sorted by increasing CA2 involvment']);
   xlabel('Events Sorted');
   ylabel('Percent');
   legend(regKey(1,:), 'Location', 'eastoutside');
   ylim([0 100]);
   saveas(gcf,strcat(specPath, saveSchem, '.percEvtRegCA2sort.png')); 
end

% Sort by duration
[~,edges,bins] = histcounts(spkEventTimes.EventDuration);
durCount = zeros(max(bins),size(regKey,2));
binct = zeros(max(bins),1);
for i = 1:length(bins)
    for j = 1:size(regKey,2)
        durCount(bins(i),j) = durCount(bins(i),j)+percEvt(i,j);
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
legend(regKey(1,:), 'Location', 'eastoutside');
saveas(gcf,strcat(specPath, saveSchem, '.percDurReg.png'));

% Box plot
boxx = [];
boxg = [];
for r = 1:size(regKey,2)
    useP = cat(2,percEvt(:,r));
    boxx = [boxx; useP];
    boxg = [boxg; r*ones(length(useP),1)];
end
figure('Position', get(0, 'Screensize'));  
boxplot(boxx, boxg);
regInd = ismember(1:size(regKey,2),unique(boxg));
xticklabels(regKey(1,regInd));
title('Percent of units which make up events for each region');
ylabel('Percent of units in events');
ylim([0 100]);
xlabel('Region');
saveas(gcf,strcat(specPath, saveSchem, '.percEvtBox.png'));

cumMet.boxxPerc = boxx;
cumMet.boxgPerc = boxg;

if ~showPlt
    close all
end

%% [2.4] # spikes each unit fires in each event
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
saveas(gcf,strcat(specPath, saveSchem, '.avgSpkUn.png'));
barProp.avgSpkUn = avgSpkUn;

% By region
plotSubs(spikes.UID, avgSpkUn, regID, modID,spikes.UID, UIDkeep);
xlabel('Unit #');
ylabel('Average # of Spikes');
title('Average # of spikes during events when the unit is active');
saveas(gcf,strcat(specPath, saveSchem, '.avgSpkUnSort.png'));

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
saveas(gcf,strcat(specPath, saveSchem, '.avgSpkEvt.png'));
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
saveas(gcf,strcat(specPath, saveSchem, '.avgSpkEvtDur.png'));

% ImageSC of #spikes per event
cumSpk = cell(size(regKey,2),size(modKey,2));
figure('Position', get(0, 'Screensize'));
hold on;
title('Spread of number of spikes per event per cell type (when it fires)');
c = 2;
i = 1;
for r = 1:size(regKey,2)
    tempHist = [];
    tempHistNon = [];
    for j = 1:length(UIDsort{r,c})
        tempUse = unitBar.nSpkEach(UIDsort{r,c}(j)) > 0;
        tempHist(j,:) = hist(unitBar.nSpkEach(UIDsort{r,c}(j),tempUse),[2:2:10])/sum(hist(unitBar.nSpkEach(UIDsort{r,c}(j),tempUse),[2:2:10]));
        tempHistNon = [tempHistNon; hist(unitBar.nSpkEach(UIDsort{r,c}(j),tempUse),[2:2:10])];
    end
    cumSpk{r,c} = tempHistNon; % GETTING A MASSIVE # IN BIN WITH 2 SPIKES - whyyy (naybe this is fine)
    subplot(size(regKey,2),2,i); imagesc(tempHist, [0 0.2]); hold on; title(strcat(regKey(1,r),' ',modKey(1,c)));
    xticks([1 2 3 4 5]);
    xticklabels({num2str(2) num2str(4) num2str(6) num2str(8) num2str(10)});
    if r == size(regKey,1)
        ylabel('Unit');
        xlabel('# spikes');
    end
    i = i+1;
    tempHist = [];
    tempHistNon = [];
    for j = 1:length(UIDsort{r,c+1})
        tempUse = unitBar.nSpkEach(UIDsort{r,c+1}(j),:) > 0;
        tempHist(j,:) = hist(unitBar.nSpkEach(UIDsort{r,c+1}(j),tempUse),[2:2:10])/sum(hist(unitBar.nSpkEach(UIDsort{r,c+1}(j),tempUse),[2:2:10]));
        tempHistNon = [tempHistNon; hist(unitBar.nSpkEach(UIDsort{r,c+1}(j),tempUse),[2:2:10])]; %non-normalized
    end
    cumSpk{r,c+1} = tempHistNon;
    subplot(size(regKey,2),2,i); imagesc(tempHist, [0 0.2]); hold on; title(strcat(regKey(1,r),' ',modKey(1,c+1)));
    xticks([1 2 3 4 5]);
    xticklabels({num2str(2) num2str(4) num2str(6) num2str(8) num2str(10)});
    i = i+1;
end
hold off;
saveas(gcf,strcat(specPath, saveSchem, '.numSpkSC.png'));
barProp.cumSpk = cumSpk;
cumMet.cumSpk = cumSpk;
if ~showPlt
    close all
end

%% [2.5] FR per unit (per event)
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
saveas(gcf,strcat(specPath, saveSchem, '.avgFRun.png'));
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
saveas(gcf,strcat(specPath, saveSchem, '.avgFRevt.png'));
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
saveas(gcf,strcat(specPath, saveSchem, '.avgFRevtDur.png'));

% Firing rate imageSC
figure('Position', get(0, 'Screensize'));
hold on;
title('Spread of FR per event per cell type');
c = 2;
i = 1;
cumFR = cell(size(regKey,2),size(modKey,2));
for r = 1:size(regKey,2)
    tempHist = [];
    tempHistNon = [];
    for j = 1:length(UIDsort{r,c})
        tempUse = [];
        tempUse = unitBar.FReach(UIDsort{r,c}(j))>0;
        tempHist(j,:) = hist(unitBar.FReach(UIDsort{r,c}(j),tempUse),[2:2:20])/sum(hist(unitBar.FReach(UIDsort{r,c}(j),tempUse),[2:2:20]));
        tempHistNon = [tempHistNon; hist(unitBar.FReach(UIDsort{r,c}(j),tempUse),[2:2:20])];
    end
    cumFR{r,c} = tempHistNon;
    subplot(size(regKey,2),2,i); imagesc(tempHist, [0 0.05]); hold on; title(strcat(regKey(1,r),' ',modKey(1,c)));
    xticks([1:10]);
    xticklabels({num2str(2) num2str(4) num2str(6) num2str(8) num2str(10)...
        num2str(12) num2str(14) num2str(16) num2str(18) num2str(20)});
    if r == size(resCell,1)
        ylabel('Unit');
        xlabel('Firing Rate');
    end
    i = i+1;
    tempHist = [];
    tempHistNon = [];
    for j = 1:length(UIDsort{r,c+1})
        tempUse = [];
        tempUse = unitBar.FReach(UIDsort{r,c+1}(j))>0;
        tempHist(j,:) = hist(unitBar.FReach(UIDsort{r,c+1}(j),tempUse),[2:2:20])/sum(hist(unitBar.FReach(UIDsort{r,c+1}(j),tempUse),[2:2:20]));
        tempHistNon = [tempHistNon; hist(unitBar.FReach(UIDsort{r,c+1}(j),tempUse),[2:2:20])];
    end
    cumFR{r,c+1} = tempHistNon;
    subplot(size(regKey,2),2,i); imagesc(tempHist, [0 0.05]); hold on; title(strcat(regKey(1,r),' ',modKey(1,c+1)));
    xticks([1:10]);
    xticklabels({num2str(2) num2str(4) num2str(6) num2str(8) num2str(10)...
        num2str(12) num2str(14) num2str(16) num2str(18) num2str(20)});
    i = i+1;
end
hold off;
saveas(gcf,strcat(specPath, saveSchem, '.FRsc.png'));
barProp.cumFR = cumFR;
cumMet.cumFR = cumFR;
if ~showPlt
    close all
end

%% [2.6] Population SAVE
saveStructs(barProp,specPath,basename)
save(strcat('Z:\home\Lindsay\Barrage\CumMet\',animName,'.',basename,'.cumMet.mat'),'cumMet', '-v7.3');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

function saveStructs(structure,specPath,basename)
    unpackStruct(structure);
    clear structure FR %SOME THINGS ARE NOT BEING SAVED BECAUSE THEY ARE TOO BIG, BUT THAT'S OKAY FOR NOW
    save(strcat(specPath,basename,'.props.mat'));
end

function unpackStruct(structure)
    fn = fieldnames(structure);
    for i = 1:numel(fn)
        fni = string(fn(i));
        field = structure.(fni);
        assignin('caller', fni, field);
    end
end