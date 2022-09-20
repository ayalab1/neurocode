function [r,p] = CofiringCoefficient(spikes,intervals1,intervals2)

% [CofiringCoefficient] - [Compute the cofiring coefficient between two
% sets of intervals]
%
% [Compute the cofiring coefficient of unit pairs between two sets of intervals
% (such as theta cycles and ripples). First, count the spikes emitted by each
% unit in each interval of the first set. Compute Pearson's correlation between
% the spike counts for each cell pair, and store the results in a vector. Then
% repeat this for the second set of intervals. The cofiring coefficient is
% defined as Pearson's correlation between the two resulting vectors. Thus,
% a high value indicates that cell pair correlations are similar between the two
% sets of intervals.]
%
%  USAGE
%
%    [r,p] = CofiringCoefficient(spikes,intervals1,intervals2)
%
%  INPUT
%
%    [spikes]         [two-column matrix of timestamps and unit IDs,
%                      prvided by <a href="matlab:help
%                      GetSpikeTimes">GetSpikeTimes</a> using 'output' = 'numbered']
%    [intervals1]     [two-column matrix of start and end times]
%    [intervals2]     [two-column matrix of start and end times]
%
%  OUTPUT
%
%    [r]              [Pearson's r between the neuronal pair correlations in
%                      intervals1 and in  intervals2 (cofiring coefficient)]
%    [p]              [p-value for Pearson's test]
%
%  EXAMPLES
%
%  SEE
%
%   Dependencies: CorrInIntervals, InIntervals, Accumulate
%
% [Ralitsa Todorova and Michaël Zugaro] [2014-2021]
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


% ------------------------------- Helper function -------------------------------

function [corrMatrix,pMatrix,psr] = CorrInIntervals(spikes, intervals, varargin)

% Computes Pearson's correlation for each pair of neuron's spikes within
% the intervals provided. For each interval, the spikes emitted by each
% unit are counted, and if units M and N fire in a correlated manner within the
% provided windows, corr(M,N) would be high. corr(k,k), (which would always be
% equal to 1) is set to NaN by default.
%
%
% Copyright (C) 2017 by Ralitsa Todorova
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

output = 'triangular';
intervalIDs = (1:length(intervals))';

for i = 1:2:length(varargin),
    if ~ischar(varargin{i}),
        error(['Parameter ' num2str(i+2) ' is not a property (type ''help <a href="matlab:help CorrinIntervals">CorrinIntervals</a>'' for details).']);
    end
    switch(lower(varargin{i})),
        case 'output',
            output = varargin{i+1};
            if ~isstring(output,'full','triangular'),
                error('Incorrect value for property ''output'' (type ''help <a href="matlab:help CorrinIntervals">CorrinIntervals</a>'' for details).');
            end
        case 'ids', %interval IDs (in case of within-interval gaps)
            intervalIDs = varargin{i+1};
            if ~isvector(intervalIDs),
                error('Incorrect value for property ''ids'' (type ''help <a href="matlab:help CorrinIntervals">CorrinIntervals</a>'' for details).');
            end
            
        otherwise,
            error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help CorrinIntervals">CorrinIntervals</a>'' for details).']);
    end
end

n = max(spikes(:,2));
id = spikes(:,2);

psr = zeros(size(intervals,1),max(id));
for i=1:n,
    psr(:,i) = CountInIntervals(spikes(id==i),intervals);
end

if ~isempty(psr)
    [corrMatrix pMatrix] = corr(psr);
else
    [corrMatrix pMatrix] = deal(nan(size(psr,2)));
end

% Compute the correlation matrix C. Because C is symmetrical (Cij=Cji), keep only upper triangular elements
% in order to remove duplicates ((M,N)=(N,M) and diagonal), so that each pair is counted only once.

if strcmpi(output,'triangular'),
    corrMatrix(logical(tril(ones(size(corrMatrix))))) = NaN;
    pMatrix(logical(tril(ones(size(pMatrix))))) = NaN;
elseif strcmpi(output,'full'),
%     corrMatrix(logical(eye(size(corrMatrix)))) = NaN;
end


