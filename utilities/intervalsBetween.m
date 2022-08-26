function [nonEvents] = intervalsBetween(start,stop,sr,events)
%intervalsBetween - Pull times between event periods
%
%  This function is for pulling the intervals between event times as if
%  they are individual events (exporting start and stop times for these
%  periods).
%
%  USAGE
%
%    % Identify non-ripple times during sleep
%    [nonRippleTimes] = intervalsBetween(start_preSleep, stop_preSleep, ripples_times)
%
%  INPUTS
%    start       Start time (integer, seconds) of your overarching period
%    stop        Stop time (integer, seconds) of your overarching period
%    sr          Sampling rate of the timestamps (Hz)
%    events      Nx2 matrix of ordered start and stop times (seconds) for N
%                events of interest. 
% 
%  OUTPUTS
%    nonEvents      Mx2 matrix of start and stop times (seconds) for M 
%                   periods between events during [start stop] period
%
%  EXAMPLES
%
%    intervalsBetween(0,10,0.1,[1,2; 5,7]) returns [0,0.9; 2.1,4.9; 7.1,10]
%
%  SEE
%
%    ...

% Lindsay Karaba, 2022
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

nonEvents = [];
nEct = 1;
t_int = 1/sr;
if isstring(stop)
   stop = str2double(stop); 
end

startEvt = find(events(:,1) >= start,1,'first');
stopEvt = find(events(:,2) <= stop,1,'last');

if ~(events(startEvt,1)==start)
    nonEvents(nEct,1) = start;
    nonEvents(nEct,2) = events(startEvt,1)-t_int;
    nEct = nEct+1;
end

for i = startEvt:stopEvt-1
    nonEvents(nEct,1) = events(i,2)+t_int;
    nonEvents(nEct,2) = events(i+1,1)-t_int;
    nEct = nEct+1;
end

if ~(events(stopEvt,2)==stop)
    nonEvents(nEct,1) = events(stopEvt,2)+t_int;
    nonEvents(nEct,2) = stop;
end
    
end