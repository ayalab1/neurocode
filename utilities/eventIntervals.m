function eventInt = eventIntervals(events,intervals,collapse)
%eventIntervals - Separate given  event struct by given time intervals
%   events -  event structure containing pulses
%   intervals - set of time intervals
%   collapse - collpase individual struct for each interval into one
%   Requires: first field of events struct must be timestamps
%   Jacob Kerr, 2020.
%   AntonioFR, 2021

if ~isstruct(events)
    eventsL = load(events);
    loadFields = fieldnames(events);
    events = getfield(eventsL,loadFields);
end

loadFields = fieldnames(events);
fields = fieldnames(events);

if ~isnumeric(intervals)
    load(intervals);
end

for i = 1:size(intervals,1)
    tmp = InIntervals(getfield(events,fields{1,1}),intervals(i,:)); % this function needs to be compiled
    timestampField = getfield(events,fields{1,1});
    eventInt{i}.timestamps = timestampField(tmp==1,:);
    for j = 2:size(fields)
        if size(getfield(events, fields{j,1}),1) == size(getfield(events,fields{1,1}),1)
            tempArr = getfield(events, fields{j,1});
            eventInt{i}.(fields{j,1}) = tempArr(tmp==1,:);
        elseif size(getfield(events, fields{j,1}),1) == size(getfield(events,fields{1,1}),1) &&...
                size(getfield(events, fields{j,1}),2) > 1
            newTmp = InIntervals(getfield(events,fields{j,1}),intervals(i,:));
            tempArr = getfield(events, fields{j,1});
            eventInt{i}.(fields{j,1}) = tempArr(newTmp==1,:);
        else
            eventInt{i}.(fields{j,1}) = getfield(events, fields{j,1});
        end
    end
    %clear tmp1;
end

if collapse
    for i = 1:numel(eventInt)
        if isempty(eventInt{i}.timestamps)
            eventInt{i} = [];
        end
    end
    eventInt(cellfun('isempty', eventInt)) = [];
    eventT = cell2mat(eventInt);
    eventInt = collapseStruct(eventT,'match');
end
end


