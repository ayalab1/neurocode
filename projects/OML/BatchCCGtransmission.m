function [pTE,pTI,pAll,pCA3,regionCode,varargin,basepath] = BatchCCGtransmission(basepath,varargin)

binSize = 4/10000; duration = 0.12; % cell excplorer ce_MonoSynConvClick parameters
[spikes,regionID,regionNames,spikesCell] = GetAyaSpikes(basepath);
cell_metrics = getStruct(basepath,'cell_metrics');
MergePoints = getStruct(basepath,'MergePoints');
SleepStateEpisodes = getStruct(basepath,'SleepStateEpisodes');
sleep = ConsolidateIntervals(MergePoints.timestamps(cellfun(@(x) any(strfind(lower(x),'sleep')),MergePoints.foldernames),:),'epsilon',1);
sws = SleepStateEpisodes.ints.NREMepisode;
pyr = cellfun(@(x) ~isempty(strfind(x,'Pyramidal')), cell_metrics.putativeCellType)';
if isempty(varargin)
    cheese = false;
    load(fullfile(basepath,[basenameFromBasepath(basepath) '.placeFields.cellinfo']),'placeFieldStats');
    isplacecell = [cellfun(@(x) isnan(x{1}.x(1)),placeFieldStats.mapStats) cellfun(@(x) isnan(x{2}.x(1)),placeFieldStats.mapStats)];
else
    cheese = true;
    condition = varargin{1};
    if condition==0, isplacecell = [zeros(size(pyr)) ones(size(pyr))]; else isplacecell = [ones(size(pyr)) zeros(size(pyr))]; end
end

category = isplacecell(:,1)*1 + isplacecell(:,2)*2+1;  % 1 = not a place cell, 2 = place cell only for stim direction, 3 = place cell only for off direction, 4 = place cell for both directions

CA1 = find(cellfun(@(x) any(strfind(lower(x),'ca1')),regionNames));
CA3 = find(cellfun(@(x) any(strfind(lower(x),'ca3')),regionNames));
Unknown = find(cellfun(@(x) any(strfind(lower(x),'unknown')),regionNames));
regionCode = ismember(regionID,CA1)*1 + ismember(regionID,Unknown)*2 + ismember(regionID,CA3)*3;

pTE = []; pTI = [];
if cheese, excConditions = 2; else excConditions = 0:2; end
try
    for exc = excConditions
        if exc == 0 % Inhibitory connections
            list = cell_metrics.putativeConnections.inhibitory;
        elseif exc==1 % Excitatory connections
            list = cell_metrics.putativeConnections.excitatory;
        else % all cell pairs (connected or not)
            nUnits = length(spikesCell);
            list = [repelem((1:nUnits)',nUnits) repmat((1:nUnits)',nUnits,1)];
        end
        if isempty(list), list = zeros(0,2); end

        pTransmission = nan(size(list,1),8);
        for period = 1:2
            interval = Restrict(sws,sleep(period,:));
            restricted = Restrict(spikes,interval);
            [ccg,t] = CCG(restricted(:,1),restricted(:,2),'binSize',binSize,'duration',duration);
            for i = 1:size(list,1)
                try
                ccgI = ccg(:,list(i,1),list(i,2));
                spikesI = sum(restricted(:,2)==list(i,1));
                [trans,prob,prob_uncor,pred] = ce_GetTransProb(ccgI,  spikesI,  binSize,  0.020);
                pTransmission(i,period) = trans;
                catch
                    pTransmission(i,period)  = nan;
                end
            end
        end
        pTransmission(:,3:4) = list;
        pTransmission(:,5) = pyr(list(:,2)); % is the target a pyramidal cell (1) or a pyramidal neuron (0)
        pTransmission(:,6) = category(list(:,1)); % is the source cell a place cell (no=1, on direction = 2, off direction = 3, both = 4);
        pTransmission(:,7:8) = regionCode(list);
        if exc==0, pTI = pTransmission; elseif exc==1, pTE = pTransmission; else, pAll = pTransmission; end
    end
catch
    keyboard
end

spikes0 = spikes;
%% response to CA3

% add an "pooled CA3" unit
pooledCA3 = spikes(ismember(spikes0(:,2),find(regionCode==3 & pyr))); pooledCA3(:,2) = nUnits+1;

regionCode(nUnits+1) = 3;
spikes = sortrows([spikes0; pooledCA3]);

list = [ones(nUnits,1)*nUnits+1 (1:nUnits)'];
pTransmission = nan(size(list,1),8);
for period = 1:2
    interval = Restrict(sws,sleep(period,:));
    restricted = Restrict(spikes,interval);
    [ccg,t] = CCG(restricted(:,1),restricted(:,2),'binSize',binSize,'duration',duration);
    for i = 1:size(list,1)
        try
            ccgI = ccg(:,list(i,1),list(i,2));
            spikesI = sum(restricted(:,2)==list(i,1));
            [trans,prob,prob_uncor,pred] = ce_GetTransProb(ccgI,  spikesI,  binSize,  0.020);
            pTransmission(i,period) = trans;
        catch
            pTransmission(i,period)  = nan;
        end
    end
end

pTransmission(:,3:4) = list;
pTransmission(:,5) = pyr(list(:,2)); % is the target a pyramidal cell (1) or a pyramidal neuron (0)
pTransmission(:,6) = category(list(:,2)); % is the target cell a place cell (no=1, on direction = 2, off direction = 3, both = 4);
pTransmission(:,7:8) = regionCode(list);

pCA3 = pTransmission;


