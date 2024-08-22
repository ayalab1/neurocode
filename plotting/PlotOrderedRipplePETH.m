function PlotOrderedRipplePETH(varargin)

%PlotOrderedRipplePETH - Compute an event-triggered PETH where ripples are
% sorted according to their timing relative to a third event (e.g. delta waves).
%
% This function is useful to check, for example, the cortical responses to 
% ripples and how they depend on ripple-delta coupling. 
% [mat,t,m] = PETH(data,events,<options>)
%
%  USAGE
%
%    PlotOrderedRipplePETH(<options>)
%
%    <options>      optional list of property-value pairs (see table below)
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     basepath          directory where spiking and events will be pulled from
%     triggerEvent      name for trigger event file which will be the center of
%                       the PETH (default = 'ripples')
%     triggerSelection  function handle in case you want to select for specific
%                       trigger events. This is a function that will be applied
%                       on the event structure, for example: @(x) x.amplitude>100
%                       or @(x) x.duration<0.12 (default = no restriction).
%     sortingEvent      name of the event file which will be used to sort
%                       the individual rows of the PETH (default = deltaWaves)
%     brainRegion       a cell containing the acceptable brainRegion tags to 
%                       be pulled from cell_metrics. Only spikes from units with the 
%                       provided tags will be used to compute the PETHs.
%                       (default = {'PFC', 'MEC', 'EC', 'ILA', 'PL','Cortex'})
%     durations         durations before and after each synchronizing event
%                       (in s) (default = [-0.5 0.5])
%     nBins             number of time bins around the events (default = 201)
%     smooth            standard deviation for Gaussian kernel (default = 2 bins)
%                       applied to the PETH
%     nQuantiles        number of rows for the ordered PETH. The triggerEvents
%                       are divided into nQuantiles groups according to their
%                       timing around the sorting events, and each row of the 
%                       plotted color map will represent the mean spiking
%                       response profile to the triggerEvents of each group
%    =========================================================================
%
%  OUTPUT
%  There is no output. The function plots a figure.
%
%  EXAMPLES
%
%  % observe a bad cortical response to ripples that is entirely due to delta wave timing:
%  cd Y:\OJRproject\OJR42\day12; figure; PlotOrderedRipplePETH;
% 
%  % observe a good cortical response to ripples regardless of delta wave timing:
%  cd Y:\OJRproject\OJR34\day9; figure; PlotOrderedRipplePETH;
%
% observe hippocampal responses to all opto-pulses as a function of ripple timing:
% figure; PlotOrderedRipplePETH('brainRegion',{'CA1','HPC','HPF'},...
%       'sortingEvent','ripples','triggerEvent','pulses','durations',[-1 1]*0.2);
%
% observe hippocampal responses to opto-pulses shorter than 120ms as a function of ripple timing:
% figure; PlotOrderedRipplePETH('brainRegion',{'CA1','HPC','HPF'},...
%       'sortingEvent','ripples','triggerEvent','pulses','durations',[-1 1]*0.2,...
%        'triggerSelection',@(x) x.duration<0.12);
%
%  SEE
%
%    See also PETH, JointPETH, PlotJoingPETH.
%
% Copyright (C) 2024 by Ralitsa Todorova
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.


% parse inputs
p = inputParser;
p.PartialMatching = true; % allows the use of abbreviations
addParameter(p, 'basepath', pwd, @(x) any([isfolder(x), iscell(x)]));
addParameter(p, 'triggerEvent', 'ripples', @ischar);
addParameter(p, 'triggerSelection', @(x) true(size(x.timestamps(:,1))), @isfun);
addParameter(p, 'sortingEvent', 'deltaWaves', @ischar);
addParameter(p, 'brainRegion', {'PFC', 'MEC', 'EC', 'ILA', 'PL','Cortex'}, @(x) iscell(x) || ischar(x));
addParameter(p, 'gapDuration', [1 4], @isnumeric);
addParameter(p, 'durations', [-1 1]*0.5, @isnumeric);
addParameter(p, 'nBins', 201, @isnumeric);
addParameter(p, 'nQuantiles', 100, @isnumeric);
addParameter(p, 'smooth', 2, @isnumeric);
parse(p, varargin{:})
basepath = p.Results.basepath;
triggerEventName = p.Results.triggerEvent;
triggerSelection = p.Results.triggerSelection;
sortingEventName = p.Results.sortingEvent;
regions = p.Results.brainRegion;
durations = p.Results.durations;
nBins = p.Results.nBins;
nQuantiles = p.Results.nQuantiles;
smooth = p.Results.smooth;


[parentFolder,basename] = fileparts(basepath);
[~,projectName] = fileparts(parentFolder);

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
triggerEvents = triggerStruct.timestamps(triggerSelection(triggerStruct),1);

[h,~] = PETH(spikes(:,1),triggerEvents,'durations',durations,'nBins',nBins);
[h_sorting,ht] = PETH(sortingEvents,triggerEvents,'durations',durations,'nBins',nBins);

rt = triggerEvents - sortingEvents(FindClosest(sortingEvents,triggerEvents)); % rt = relative time
ok = InIntervals(rt,durations);
[~,m] = max(Shrink(sortby(h_sorting(ok,:),rt(ok)),floor(sum(ok)/nQuantiles),1),[],2);
subplot(2,1,1);
PlotColorMap(Smooth(Shrink(sortby(h(ok,:),rt(ok)),floor(sum(ok)/nQuantiles),1),smooth),'x',ht);
hold all
plot(ht(m),1:length(m),'w','linewidth',2);
PlotHVLines(0,'v','k--','linewidth',2);
xlabel(['time from ' triggerEventName ' (s)']); ylabel([triggerEventName ' ID ordered by ' sortingEventName ' proximity']);
set(gca,'box','off','TickDir','out','fontsize',12);
title([strrep([projectName '-' basename],'_','-') ': ' strjoin(regionNames(okRegion),', ')]);
drawnow

subplot(2,1,2);
semplot(ht,h/diff(ht(1:2)),'k',1);
set(gca,'box','off','TickDir','out','fontsize',12);
ylabel('mean multiunit firing rate (Hz)');
PlotHVLines(0,'v','k--','linewidth',2);


% % Alternative code, correcting for the effect of the sorting events:
% subplot(2,1,2); cla
% [h1,ht] = PETH(spikes(:,1),sortingEvents,'durations',durations,'nBins',nBins);
% [h2,ht] = PETH(triggerEvents,sortingEvents,'durations',durations,'nBins',nBins);
% [~,~,d] = JointPETH(h1/diff(ht(1:2)),h2/diff(ht(1:2)),1);
% corrected = CircularShift(d,ceil(size(d,2)/2)-(1:size(d,2)))';
% semplot(ht,corrected); set(gca,'box','off','TickDir','out','fontsize',12);
% ylabel('corrected response (Hz)'); PlotHVLines(0,'v','k--','linewidth',2);





















