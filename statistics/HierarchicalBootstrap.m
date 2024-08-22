function [p,means] = HierarchicalBootstrap(data,varargin)

% HierarchicalBootstrap:
%
% The function implements a hierarchical bootstrap as described in
% Saravanan et al. (2021) https://doi.org/10.1101/819334 with the addition
% of a within-subject design bootstrap, where data from the same e.g.
% animal is sampled together for both conditions (rather than independently
% for each condition, which would ignore the paired design).
% The function will give you a p-value testing if the values of condition
% 1 are systematically higher than the values corresponding to condition 2,
% after controlling for all the in-between nesting variables.
% The output "means" is the bootstrapped average of the values of conditions
% 1 and 2 for each of the iterations.
%
%  USAGE
%
%    [p,means] = HierarchicalBootstrap(data)
%
%    data           The expected format is long format, with the first 'nColumns'
%                   columns being the values to be tested, and subsequent
%                   columns indicate nested groupings of the data, from
%                   deepest to highest (see EXAMPLE)
%    nColumns       Number of columns in which data is present (default = 1).
%                   Other columns are assumed to contain grouping values.
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'nIterations' how many iterations to perform, resampling from the data
%                   (default = 1000)
%     'fun'         the function handle of the function computing the average
%                   value per condition (default = @mean)
%    =========================================================================
%
%  OUTPUT
%
%    p              p-value of the test that the values of condition 1 are
%                   higher than the values of condition 2 after controlling
%                   for the other nesting variables
%    means          bootstrapped results for the mean of each condition (one for
%                   each iteration)
%
%
% EXAMPLE
%
% data = [fieldSize, cellID, session, animal, condition]
% [p,means] = HierarchicalBootstrap(data);
% "p<0.05" indicates that the fieldSizes of condition 1 are significantly
% higher than the fieldSizes of condition 2, after controlling for
% cellID, session, and animal.
% pValues = mean(means<=0); % is data from the each group significantly higher than zero
%
% Copyright (C) 2023-2024 by Ralitsa Todorova
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

if iscell(data)
    % find out how deep the cell nesting goes:
    % Somehow group them into a matrix:
    error('Grouping of nested data in cells not yet implemented. Please provide matrix format');
end

% if any(~ismember(data(:,end),[1 2]))
%     error('The last column should provide the grouping variable (1 or 2)');
% end

nIterations = 1000;
average = @nanmean; % change to nanmedian for testing if the group medians are different
nColumns = 1;

% Parse parameters
for i = 1:2:length(varargin),
    if ~ischar(varargin{i}),
        builtin('error',['Parameter ' num2str(i+2) ' is not a property (type ''help <a href="matlab:help HierarchicalBootstrap">HierarchicalBootstrap</a>'' for details).']);
    end
    switch(lower(varargin{i})),
        case {'niterations','nshuffles','n'}
            nIterations = varargin{i+1};
            if ~(isscalar(nIterations) && (nIterations>0))
                builtin('error',['Incorrect value for property ''' varargin{i} ''' (type ''help <a href="matlab:help HierarchicalBootstrap">HierarchicalBootstrap</a>'' for details).']);
            end
        case {'fun','function','handle','average'}
        	  average = varargin{i+1};
        	  if ~isa(average,'function_handle')
        	      builtin('error',['Incorrect value for property ''' varargin{i} ''' (type ''help <a href="matlab:help HierarchicalBootstrap">HierarchicalBootstrap</a>'' for details).']);
              end
        case {'ncolumns'}
            nColumns = varargin{i+1};
            if ~(isscalar(nColumns) && (nColumns>0) && nColumns<size(data,2)-1)
                builtin('error',['Incorrect value for property ''' varargin{i} ''' (type ''help <a href="matlab:help HierarchicalBootstrap">HierarchicalBootstrap</a>'' for details).']);
            end
        otherwise,
            builtin('error',['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help HierarchicalBootstrap">HierarchicalBootstrap</a>'' for details).']);
    end
end

% Remove empty groups (sessions, animals):
for k=2:size(data,2)-1
    [~,~,data(:,k)] = unique(data(:,k));
end
nGroups = max(data(:,end));

means = nan(nIterations,nGroups);

level=size(data,2)-1;
while level>nColumns
    % Now, consider the lewest level of grouping (e.g. at the level of the cell)
    cells0{level} = cell(max(data(:,level)),1);
    for i=1:length(cells0{level}),
        ok = data(:,level)==i;
        cells0{level}{i} = data(ok,:);
    end
    cellConditions0{level} = Accumulate(data(:,[level end]));
    level = level-1;
end

for k=1:nIterations,
    level=size(data,2)-1;
    while level>nColumns
        cells = cells0{level};
        cellConditions = cellConditions0{level};
        % Separate between and within-subject designs:
        if all(any(cellConditions==0,2)) % between-subject case
            % sample for each condition independetly
            [~,cellCondition] = max(cellConditions,[],2);
            idx = []; for j=1:nGroups, indices = find(cellCondition==j); n = length(indices); idx = [idx; indices(ceil(rand(n,1)*n))]; end
        else
            % sample on this level, regardless of condition
            n = length(cells); idx = ceil(rand(n,1)*n);
        end
        cells = cells(idx);
        if level==2 % Last level; it is time to sample with replacement the values WITHIN these cells:
            for i=1:length(cells), n=size(cells{i},1); idx = ceil(rand(n,1)*n); cells{i} = cells{i}(idx,:); end
        end
        resampled = cat(1,cells{:}); % now this level is done
        level = level-1;
    end

    % now that only the "conditions" level remains, get the mean
    for j=1:nGroups,
        try
            means(k,j) = average(resampled(resampled(:,end)==j,1:nColumns));
        catch
            means(k,j) = nan;
        end
    end
end

p = mean(means(:,1)<means(:,2));

