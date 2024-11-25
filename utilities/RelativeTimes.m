function [rt,intervalID] = RelativeTimes(t,intervals,values)

%Relative times - compute the relative times of samples within intervals
% For example, a sample with a timestamp in the middle of an interval
% would have relative time 0.5. A sample that falls outside of the provided
% intervals would return as NaN. Multiple-column intervals are supported.
%
%   USAGE
%
%    [r,w] = RelativeTimes(t,intervals,values)
%
%    t              timestamps that will be converted in relative time
%    intervals      the [start stop] timestamps of the intervals that will
%                   serve as basis for computing the relative time. Multiple
%                   columns are allowed (for example [start peak stop] 
%                   timestamps)
%    values         the values assigned to each column in intervals
%                   (default = 0:size(intervals,2)-1), so values coinciding
%                   with intervals(:,1) would have a value of 0 and values
%                   coinciding with intervals(:,end) would have a value of
%                   size(intervals,2)-1.
%
%   OUTPUT
%
%     rt            relative times (one for each value of t)
%     intervalID    for each value of t, the row of 'intervals' in which
%                   it was detected
%
%   EXAMPLE
%
% if we have a spindle power signal [t spindlePower], then:
% rt = RelativeTimes(t, spindles(:,[1 2 3]),0:2)
% b = ceil(rt*20); % bins
% m = Accumulate(b(~isnan(b)),spindlePower(~isnan(b)),'mode','mean');
% plot(linspace(-1,1,3),m); 
% will plot the average spindle power in relative spindle time, with the
% first half covering the power between the spindle start to the spindle
% peak, and the second half covering the power between the spindle peak to
% the spindle end.
%
% (C) 2016 by Ralitsa Todorova
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.


if nargin<3
    values = 0:size(intervals,2)-1;
end

[rt,intervalID] = deal(nan(length(t),size(intervals,2)-1));
for i=1:size(intervals,2)-1
    [rt(:,i),intervalID(:,i)] = RelativeTime(t,intervals(:,[i i+1]));
end

% In case a timestamp is present in more than one event, we will select the
% event that is closest to the timestamp, and ignore the timestamps's position relative to
% other events. (In the above example, if a timestamp is 0.2 seconds after a spindle end, 
% but also 0.1 seconds before a spindle start, the 0.2 value will be ignored).
[~,columnsToKeep] = max(abs(0.5-rt),[],2);
indicesToKeep = sub2ind(size(rt),(1:length(t))',columnsToKeep);

for i=1:size(intervals,2)-1
    rt(:,i) = rt(:,i)*(values(i+1)-values(i)) + values(i);
end

rt = rt(indicesToKeep); intervalID = intervalID(indicesToKeep);

% ------------------------------- Helper function -------------------------------

function [rt, intervalID] = RelativeTime(t, intervals)

t = t(:,1); %only timestamps are taken

if any(diff(intervals(:,1))<0) % we need to order the intervals so that they are monotonically increasing
    sorting = true;
    [~,order] = sort(intervals(:,1));
    intervals = intervals(order,:);
else
    sorting = false;
end

rt = nan(size(t));

if size(intervals,2)==3
    intervalID = zeros(length(t), size(intervals,2)-1);
    for i=1:size(intervals,2)-1
        theseIntervals = intervals(:,i:i+1);
        [inIntervals,intervalID(:,i)] = InIntervals(t, theseIntervals);
        rt(inIntervals) = RelativeTime(t(inIntervals), theseIntervals)+i-1;
    end
    return
end
[inIntervals, intervalID] = InIntervals(t, intervals);
rt(inIntervals,1) = (t(inIntervals) - intervals(intervalID(inIntervals),1))./...
    (intervals(intervalID(inIntervals),end) - intervals(intervalID(inIntervals),1));

if sorting
    intervalID(intervalID>0) = order(intervalID(intervalID>0));
end

