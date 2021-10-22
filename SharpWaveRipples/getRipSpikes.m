function [spkEventTimes] = getRipSpikes(varargin)
%
%    [spkEventTimes] = getRipSpikes(varargin)
% Saves spike times of all units inside given events in different ways:
%   1. Absolute and relative time of spikes by unit and by event or
%      nonevent period
%   2. Absolute and relative time of spikes by unit
%   3. Absolute and relative time of spikes by event
%   4. Absolute and relative time of spikes outside events
%
% INPUTS
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'basepath'    full path where session is located (default pwd)
%     'events'      It can be either the following options:
%                   1. Buzcode structure for specific events (ripples, UDStates, ...)
%                      By default it will load ripples (output from bz_DetectSWR).
%                      Specifically, if not provided, it loads this ripple 
%                      structure from 'basepath' (if provided), or from current
%                      folder (if not). The input structure must have the 
%                      following field:
%                        .timestamps: Nx2 matrix with starting and ending times 
%                                     (in sec) of each event.
%                   2. A Nx2 matrix with starting and ending times (in sec)
%                      of each event, just like the .timestamps field.
%                       (N: number of events)
%     'spikes'      buzcode ripple structure (from bz_GetSpikes). 
%                   If not provided, it loads it from 'basepath' (if provided),
%                   or from current folder (if not)
%     'UIDs'        A Mx1 boolean matrix with 1s for units to be considered
%                   and 0s for units to be discarded.(M: number of units)
%     'padding'     additional time before and after ripple end to still search for spikes. 
%                   (default is 0.05 sec)
%     'savePath'    Definition of the path to save the output structure
%                   (default: pwd)
%     'saveMat'   	Saves file, logical (default: true) 
%
%    =========================================================================
%
% OUTPUTS
%
% spkEventTimes structure with the following fields:
%
%    =========================================================================
%     Properties        Values
%    -------------------------------------------------------------------------
%     .padding          Value of padding used
%     .EventDuration    Nx1 double. Contains the duration of each event
%     .NonEventDur      (N+1)x1 double. Contains the duration of each 
%                       interval between events.
%	  .UnitEventAbs     MxN cell matrix. In each cell, absolute times of
%                       spikes for that particular unit and event
%	  .UnitEventRel     MxN cell matrix. In each cell, relative times of
%                       spikes (relative to the starting time of events)  
%                       for that particular unit and event
%     .UnitNonEventAbs  Mx(N+1) cell matrix. In each cell, absolute times
%                       of spikes for that particular unit and intermittent
%                       period between events.
%     .UnitNonEventRel  Mx(N+1) cell matrix. In each cell, relative times
%                       of spikes for that particular unit and intermittent
%                       period between events.
%	  .UnitAbs          1xM cell matrix. In each cell, absolute times of
%                       spikes for that particular unit across all events
%	  .UnitRel          1xM cell matrix. In each cell, relative times of
%                       spikes for that particular unit across all events
%     .UnitNonAbs       1xM cell matrix. In each cell, absolute times of
%                       spikes for that particular unit across all
%                       non-event periods.
%     .UnitNonRel       1xM cell matrix. In each cell, relative times of
%                       spikes for that particular unit across all
%                       non-event periods.
%	  .EventAbs         2xN cell matrix. In the first row, absolute times 
%                       of spikes for that particular event across all 
%                       units. In the second row, the UID associated to the 
%                       above spike
%	  .EventRel         2xN cell matrix. In the first row, relative times 
%                       of spikes for that particular event across all 
%                       units. In the second row, the UID associated to the
%                       above spike
%     .NonEventAbs      2x(N+1) cell matrix. In the first row, absolute
%                       times of spikes for that particular period between
%                       events, across all units. In the second row, the
%                       UID associated to the above spike. 
%     .NonEventRel      2x(N+1) cell matrix. In the first row, relative
%                       times of spikes for that particular period between
%                       events, across all units. In the second row, the
%                       UID associated to the above spike.
%
%    =========================================================================
%
%    Antonio FR, 2017

% Parse inputs 
p = inputParser;
addParameter(p,'basepath',pwd,@isstr);
addParameter(p,'events',[], @(x) isnumeric(x) || isstruct(x));
addParameter(p,'spikes',{},@isstruct);
addParameter(p,'UIDs',[],@islogical);
addParameter(p,'padding',0.05,@isnumeric);
addParameter(p,'savePath',pwd,@isstr);
addParameter(p,'saveMat', true, @islogical);

parse(p,varargin{:});
basepath = p.Results.basepath;
events = p.Results.events;
spikes = p.Results.spikes;
UIDs = p.Results.UIDs;
padding = p.Results.padding;
savePath = p.Results.savePath;
saveMat = p.Results.saveMat;

% Get session info
basename = basenameFromBasepath(basepath);
load([basepath filesep basename '.session.mat']);
sesEpoch = session.epochs{end};
sesEnd = sesEpoch.stopTime;

% Default events, UIDs, and spikes
if isempty(spikes)
    spikes = load([basepath filesep basename '.spikes.cellinfo.mat']);
    spikes = spikes.spikes;
end
if isempty(UIDs)
    UIDs = ones(size(spikes.UID));
end
if isempty(events)
    events = load([basepath filesep basename '.ripples.events.mat']);
    events = events.ripples;
end

% Starting and ending timestamps
if isnumeric(events)
    timestamps = events;
elseif isstruct(events)
    timestamps = events.timestamps;
else
    warning('Events must be either a Nx2 vector or a bz event structure!');
end

%% Get spikes for each unit and each ripple

% We will save spike times of different units in different ways:
spkEventTimes = {};
spkEventTimes.padding = padding;

% 1. Absolute and relative time of spikes by unit and by event (or nonevent
%    period)
for unit = 1:length(spikes.UID)
    if UIDs(unit)
        for event = 1:size(timestamps,1)
            % Start and end of ripple
            tini(event) = timestamps(event,1) - padding;
            tend(event) = timestamps(event,2) + padding;
            spkEventTimes.EventDuration(event,1) = tend(event)-tini(event);
            % Spikes of this unit within this ripple interval
            tsUnitEvent = spikes.times{unit};
            tsUnitEvent = tsUnitEvent(tsUnitEvent>=tini(event) & tsUnitEvent<=tend(event));
            % Absolute time of spikes by unit and by ripple
            spkEventTimes.UnitEventAbs{unit,event} = tsUnitEvent';
            % Relative time of spikes by unit and by ripple to ripple start
            spkEventTimes.UnitEventRel{unit,event} = tsUnitEvent' - tini(event);
        end
        tini(size(timestamps,1)+1) = sesEnd;
        tend = [0 tend];
        for event = 1:length(tini)
            spkUnit = spikes.times{unit};
            spkEventTimes.NonEventDur(event,1) = tini(event)-tend(event);
            tsUnitNonEvent = spkUnit(spkUnit<tini(event) & spkUnit>tend(event));
            spkEventTimes.UnitNonEventAbs{unit,event} = tsUnitNonEvent';
            spkEventTimes.UnitNonEventRel{unit,event} = tsUnitNonEvent' - tend(event);
        end
    end
end

% 2. Absolute and relative time of spikes by unit
for unit = 1:length(spikes.UID)
    if UIDs(unit)
        spkEventTimes.UnitAbs{unit} = cell2mat(spkEventTimes.UnitEventAbs(unit,:));
        spkEventTimes.UnitRel{unit} = cell2mat(spkEventTimes.UnitEventRel(unit,:));
        spkEventTimes.UnitNonAbs{unit} = cell2mat(spkEventTimes.UnitNonEventAbs(unit,:));
        spkEventTimes.UnitNonRel{unit} = cell2mat(spkEventTimes.UnitNonEventRel(unit,:));
    end
end

% 3. Absolute and relative time of spikes by ripple
for event = 1:size(timestamps,1)
    spkEventTimes.EventAbs{event} = [];
    spkEventTimes.EventRel{event} = [];
    for unit = 1:length(spikes.UID)
        if UIDs(unit)
            spkEventTimes.EventAbs{event} = [ spkEventTimes.EventAbs{event}, ...
                                        [cell2mat(spkEventTimes.UnitEventAbs(unit,event)); ...
                                         cell2mat(spkEventTimes.UnitEventAbs(unit,event))*0+spikes.UID(unit)] ];
            spkEventTimes.EventRel{event} = [ spkEventTimes.EventRel{event}, ...
                                         [cell2mat(spkEventTimes.UnitEventRel(unit,event)); ...
                                          cell2mat(spkEventTimes.UnitEventRel(unit,event))*0+spikes.UID(unit)] ];
        end
    end
    spkEventTimes.EventAbs{event} = sortrows(spkEventTimes.EventAbs{event}')';
    spkEventTimes.EventRel{event} = sortrows(spkEventTimes.EventRel{event}')';
end

% 4. Absolute and relative time of spikes by nonevent periods
for event = 1:size(spkEventTimes.NonEventDur,1)
    spkEventTimes.NonEventAbs{event} = [];
    spkEventTimes.NonEventRel{event} = [];
    for unit = 1:length(spikes.UID)
        if UIDs(unit)
            spkEventTimes.NonEventAbs{event} = [ spkEventTimes.NonEventAbs{event}, ...
                                        [cell2mat(spkEventTimes.UnitNonEventAbs(unit,event)); ...
                                         cell2mat(spkEventTimes.UnitNonEventAbs(unit,event))*0+spikes.UID(unit)] ];
            spkEventTimes.EventRel{event} = [ spkEventTimes.NonEventRel{event}, ...
                                         [cell2mat(spkEventTimes.UnitNonEventRel(unit,event)); ...
                                          cell2mat(spkEventTimes.UnitNonEventRel(unit,event))*0+spikes.UID(unit)] ];
        end
    end
    spkEventTimes.NonEventAbs{event} = sortrows(spkEventTimes.NonEventAbs{event}')';
    spkEventTimes.NonEventRel{event} = sortrows(spkEventTimes.NonEventRel{event}')';
end

% 5. Save
if saveMat
   save(strcat(savePath,'\',basename,'.spkEventTimes.mat'),'spkEventTimes');
end

end