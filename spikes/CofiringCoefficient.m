function [r,p] = CofiringCoefficient(spikes,intervals1,intervals2)

%CofiringCoefficient - Compute the cofiring coefficient between two sets of intervals.
%
% Compute the cofiring coefficient of unit pairs between two sets of intervals
% (such as theta cycles and ripples). First, count the spikes emitted by each
% unit in each interval of the first set. Compute Pearson's correlation between
% the spike counts for each cell pair, and store the results in a vector. Then
% repeat this for the second set of intervals. The cofiring coefficient is
% defined as Pearson's correlation between the two resulting vectors. Thus,
% a high value indicates that cell pair correlations are similar between the two
% sets of intervals.
%
%  USAGE
%
%    [r,p] = CofiringCoefficient(spikes,intervals1,intervals2)
%
%  INPUT
%
%    spikes         two-column matrix of timestamps and unit IDs,
%                   provided by <a href="matlab:help GetSpikeTimes">GetSpikeTimes</a> using 'output' = 'numbered'
%    intervals1     two-column matrix of start and end times
%    intervals2     two-column matrix of start and end times
%
%  OUTPUT
%
%    r              Pearson's r between the neuronal pair correlations in
%                   intervals1 and in  intervals2 (cofiring coefficient).
%    p              p-value for Pearson's test
%

% Copyright (C) 2014-2018 by Ralitsa Todorova and Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

if nargin < 3,
    error('Incorrect number of parameters (type ''help <a href="matlab:help CofiringCoefficient">CofiringCoefficient</a>'' for details).');
end

if isempty(spikes) || ~isdmatrix(spikes,'@2') || ~isivector(spikes(:,2)),
    error('Incorrect spikes (type ''help <a href="matlab:help CofiringCoefficient">CofiringCoefficient</a>'' for details).');
end
if isempty(intervals1) || ~isdmatrix(intervals1,'@2'),
    error('Incorrect intervals1 (type ''help <a href="matlab:help CofiringCoefficient">CofiringCoefficient</a>'' for details).');
end
if isempty(intervals2) || ~isdmatrix(intervals2,'@2'),
    error('Incorrect intervals2 (type ''help <a href="matlab:help CofiringCoefficient">CofiringCoefficient</a>'' for details).');
end

% Compute Pearson's r between pairs of units, resulting in a column vector, where each line is a pair of units.
rPairs1 = reshape(CorrInIntervals(spikes,intervals1),[],1);
rPairs2 = reshape(CorrInIntervals(spikes,intervals2),[],1);

% Calculate Pearson's correlation between the pair correlations.
% (If a unit did not fire within intervals1 or intervals2, all pairs with this unit must be ignored,
% otherwise r would be NaN).
ok = ~isnan(rPairs1) & ~isnan(rPairs2); 
if sum(ok) == 0,
    r = NaN;
    return
end
if nargout > 1,
    [r,p] = corr(rPairs1(ok),rPairs2(ok));
else
    r = corr(rPairs1(ok),rPairs2(ok));
end

end

% ------------------------------- Helper function -------------------------------

function C = CorrInIntervals(spikes,intervals);

% Computes Pearson's correlations for the spikes of each pair of units within a set of intervals.
% The spikes emitted by each unit are counted in each interval, then the two spike counts are correlated.
% (Cii is set to NaN).

n = max(spikes(:,2)); % number of units


% List the interval number and unitID for each spike.
[inIntervals,intervalIDs] = InIntervals(spikes(:,1),intervals);
spikeList = [intervalIDs(inIntervals) spikes(inIntervals,2)];

% Generate the spike count matrix for each unit (column) and each interval (row).
spikeCountMatrix = Accumulate(spikeList,1,[n size(intervals,1)]);

% Compute the correlation matrix C. Because C is symmetrical (Cij=Cji), keep only upper triangular elements
% in order to remove duplicates, so that each pair is counted only once.
C = corr(spikeCountMatrix);
C(logical(tril(ones(size(C))))) = NaN;

end
