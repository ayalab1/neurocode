function [EV, REV] = ExplainedVariance(spikes, preIntervals, intervals, postIntervals)

% This returns the explained variance as reported by NcNaughton 1998.
% 'spikes' is the [timestamp neuronID] matrix, and 'preIntervals', 'intervals', and
% 'postIntervals' are the [start stop] time bins that describe pre-sleep, behaviour, 
% and post-sleep, respectively. These will be used to compute the correlation matrices
% within each epoch. It is recommended that the [start stop] bins be ~50 - 200 ms of duration,
% although larger bin sizes have been reported (Louie and Wilson, 2001)
%
% EXAMPLE:
% sleepIntervals1 = [0:0.1:552; 0.1:0.1:552.1]';
% behIntervals = [552.1:0.1:1102; 552.2:0.1:1102.1]';
% sleepIntervals2 = [1102.1:0.1:1520; 1102.2:0.1:1520.1]';
% EV = ExplainedVariance(spikes, sleepIntervals1, behIntervals, sleepIntervals2]
%
% Copyright (C) 2017 by Ralitsa Todorova
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

rMaze = CorrInIntervals(spikes,intervals);
r1 = CorrInIntervals(spikes,preIntervals);
r2 = CorrInIntervals(spikes,postIntervals);

numerator = nancorr(rMaze(:),r2(:)) - nancorr(rMaze(:),r1(:))*nancorr(r1(:),r2(:));
denominator = sqrt((1- nancorr(rMaze(:),r1(:))^2).*(1- nancorr(r1(:),r2(:))^2));
EV = (numerator/denominator)^2;

if nargout>1, REV = ExplainedVariance(spikes, postIntervals, intervals, preIntervals); end