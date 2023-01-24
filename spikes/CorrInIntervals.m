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