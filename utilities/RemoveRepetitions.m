function [clean,cleanIndices,discarded] = RemoveRepetitions(vector, varargin)

%RemoveRepetitions - Removes repetitions from a vector or a matrix.
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'mode'        'all' (default) - removes all repetitions, keeping only
%                   the first instance of the value (lie 'unique' but without 
%                   changing the order)
%                   'consecutive' - removes a value only if it is an
%                   consecutive repetition of the value preceding it.
%                   'matrix' - removes repetitions within each row - these
%                   will be substituted by NaNs.
%
%   EXAMPLES:
%   [clean,cleanIndices,discarded] = RemoveRepetitions([1 1 3 4 5 5 5 2 1]),
%   clean = [1 3 4 5 2];
%   cleanIndices = [1 3 4 5 8];
%   discarded = [0 1 0 0 0 1 1 0 1];
%
%   [clean,cleanIndices,discarded] = RemoveRepetitions([1 1 3; 4 5 5; 5 2 1], 'mode', 'matrix'),
%   clean = [1 NaN 3; 4 5 NaN; 5 2 1];
%   cleanIndices = [1 1; 1 3; 2 1; 2 2; 3 1; 3 2; 3 3];
%   discarded = [0 1 0; 0 0 1; 0 0 0];
%
%   [clean,cleanIndices,discarded] = RemoveRepetitions([1 1 3 4 5 5 5 2 1], 'mode', 'consecutive'),
%   clean = [1 3 4 5 2 1];
%   cleanIndices = [1 3 4 5 8 9];
%   discarded = [0 1 0 0 0 1 1 0 0];
%
% (C) 2016 by Ralitsa Todorova
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

mode = 'all';

% Options
for i = 1:2:length(varargin),
    if ~ischar(varargin{i}),
        error(['Parameter ' num2str(i+1) ' is not a property (type ''help <a href="matlab:help RemoveRepetitions">RemoveRepetitions</a>'' for details).']);
    end
    switch(lower(varargin{i})),
        case 'mode',
            mode = lower(varargin{i+1});
            if ~isastring(mode,'consecutive','all','matrix')
                error('Incorrect value for property ''mode'' (type ''help <a href="matlab:help RemoveRepetitions">RemoveRepetitions</a>'' for details).');
            end
        otherwise,
            error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help RemoveRepetitions">RemoveRepetitions</a>'' for details).']);
    end
end

if size(vector,2)>size(vector,1) && ~strcmp(mode, 'matrix'),
    turn = true;
    vector = vector(:);
else
    turn = false;
end

if strcmp(mode, 'consecutive'),
    discarded = [false; diff(vector)==0];
    clean = vector(~discarded);
    cleanIndices = find(~discarded);
elseif strcmp(mode, 'matrix'),
    % Each nan is unique, causing problems. Make them zeros to avoid such problems:
    discarded = false(size(vector));
    nans = isnan(vector);
    vector(isnan(vector)) = 0; % should be a safe value
    [count, value] = hist(vector', unique(vector));
    [elementIndex, sequenceIndex] = find(count>1);
    for i=unique(sequenceIndex)',
        repeatedValues = value(elementIndex(sequenceIndex==i));
        repeatedValues(repeatedValues==0) = []; %remove repeats due to NaNs
        for j=1:length(repeatedValues),
            
            thismany = sum(vector(i,:)==repeatedValues(j));
            discarded(i, find(vector(i,:)==repeatedValues(j),round(thismany-1),'last')) = 1; % letting the first one stay, the last n-1 will be discarded
        end
    end
    [x y] = find(~discarded);
    cleanIndices = sortrows([x y]);
    clean = vector;
    clean(discarded) = NaN;
    % Put the NaNs back in:
    clean(nans) = NaN;
else
    [count, value] = hist(vector, unique(vector));
    repeatedValues = value(count > 1);
    cleanIndices = find(~ismember(vector,repeatedValues));
    for i=1:length(repeatedValues),
        cleanIndices = [cleanIndices; find(vector==repeatedValues(i),1,'first')];
    end
    cleanIndices = sort(cleanIndices);
    discarded = ~Unfind(cleanIndices,length(vector));
    clean = vector(cleanIndices);
end


if turn,
    clean = clean';
    discarded = discarded';
end