function PlotOrderedRipplePETH(varargin)

% parse inputs
p = inputParser;
p.PartialMatching = true; % allows the use of abbreviations
addParameter(p, 'basepath', pwd, @(x) any([isfolder(x), iscell(x)]));
addParameter(p, 'triggerEvent', 'ripples', @ischar);
addParameter(p, 'sortingEvent', 'deltaWaves', @ischar);
addParameter(p, 'brainRegion', {'PFC', 'MEC', 'EC', 'ILA', 'PL','Cortex'}, @(x) iscell(x) || iscchar(x));
addParameter(p, 'gapDuration', [1 4], @isnumeric);
addParameter(p, 'durations', [-1 1]*0.5, @isnumeric);
addParameter(p, 'nBins', 201, @isnumeric);
addParameter(p, 'nQuantiles', 100, @isnumeric);
addParameter(p, 'smooth', 2, @isnumeric);
parse(p, varargin{:})
basepath = p.Results.basepath;
triggerEventName = p.Results.triggerEvent;
sortingEventName = p.Results.sortingEvent;
regions = p.Results.brainRegion;
durations = p.Results.durations;
nBins = p.Results.nBins;
nQuantiles = p.Results.nQuantiles;
smooth = p.Results.smooth;


[parentFolder,basename] = fileparts(basepath);
[~,projectName] = fileparts(parentFolder);
sessionID = [projectName '_' basename];

triggerStruct = getStruct(basepath,triggerEventName);
sortingStruct = getStruct(basepath,sortingEventName);

[~,regionCell,regionNames,spikesCell,~,~] = GetAyaSpikes(basepath,false);
okRegion = find(contains(regionNames,regions));
if isempty(okRegion)
    error(['No cells from the desired brain regions found.']);
end
spikes = sortrows(Group(spikesCell{ismember(regionCell,okRegion)})); 

try
    sortingEvents = sortingStruct.peaks;
catch
    sortingEvents = sortingStruct.timestamps(:,1);
end
sortingEvents = sortingEvents(:);
triggerEvents = triggerStruct.timestamps(:,1);


[h,ht] = PETH(spikes(:,1),triggerEvents,'durations',durations,'nBins',nBins);
[h_sorting,ht] = PETH(sortingEvents,triggerEvents,'durations',durations,'nBins',nBins);

rt = triggerEvents - sortingEvents(FindClosest(sortingEvents,triggerEvents));
ok = InIntervals(rt,durations);
[~,m] = max(Shrink(sortby(h_sorting(ok,:),rt(ok)),floor(sum(ok)/nQuantiles),1),[],2);
PlotColorMap(Smooth(Shrink(sortby(h(ok,:),rt(ok)),floor(sum(ok)/nQuantiles),1),smooth),'x',ht);
hold all
plot(ht(m),1:length(m),'w','linewidth',2);
PlotHVLines(0,'v','k--','linewidth',2);
xlabel(['time from ' triggerEventName ' (s)']); ylabel([triggerEventName ' ID ordered by ' sortingEventName ' proximity']);
set(gca,'box','off','TickDir','out','fontsize',15);
title(strrep([projectName '-' basename],'_','-'));
drawnow























