function truncated = TruncateIntervals(intervals, duration)

%TruncateIntervals
% 
% Sometimes it makes sense to limit the duration of the intervals you are
% analysing: for example, taking the first hour of slow wave sleep to
% score reactivation/replay. TruncateIntervals will truncate the inputted
% intervals for the desired duration, cutting away parts of the intervals
% that go after the limit.
%
%  USAGE
%
%    truncated = TruncateIntervals(intervals, duration)
%
%    intervals      list of (start,stop) pairs
%    duration       duration of the truncated intervals
%
%  NOTE
%
%    For more advanced time restriction of samples, use <a href="matlab:help InIntervals">InIntervals</a>.
%
%    *Also note that the idx output is not implemented for when sep=true
%
%  EXAMPLE
%
%    postSleep = Restrict(SleepStateEpisodes.ints.NREMepisode, [trials(end) Inf]); % get all SWS periods which follow the task
%    postSleep = TruncateIntervals(postSleep, 3600); % truncate the intervals to only the first hour of continuous SWS
%
%
% Copyright (C) 2023 by Ralitsa Todorova
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

if sum(diff(intervals,[],2))<duration % if the intervals are already shorter than the desired duration
    truncated = intervals; % there is nothing to truncate
    return
end

limit = Unshift(duration,intervals); % transform the limit in absolute time
truncated = intervals(intervals(:,1)<limit,:); % take only intervals that start before the limit is reached
if truncated(end,2)>limit, % if the end of the last interval goes over the limit
    truncated(end,2) = limit; % cut it short  
end