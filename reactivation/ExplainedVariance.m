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


% ------------------------------- Helper functions -------------------------------

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