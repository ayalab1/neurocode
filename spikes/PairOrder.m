function [score] = PairOrder(targetSequence, sequenceList, varargin)

%Apply the method used by Gupta et al (Redish) 2010
% sequenceList can be a matrix with NaNs, for example:
%
%       1   3   4   2   1   6   NaN
%       NaN 2   NaN 1   2   3   5
%       4   1   1   5   2   6   6
% where each row is a sequence. 'score' will have a score (row) for each of
% the rows in sequenceList
%
% EXAMPLE:
% [~,peakBin] = max(firingMaps,[],2);
% [~,correctOrder] = sort(peakBin); % the order of place fields (as
%                       ordered by the location of their peak)
% sequenceList = GetSequenceList(spikes, thetaCycles)
% thetaScores = PairOrder(correctOrder,sequenceList)
%
%
%  INPUTS
%    targetSequence     sequence with intended order
%    sequenceList       sequence list actually recorder. See above
%    <options>   optional list of property-value pairs (see table below)
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'circular'        if data is circular 'on.' If not 'off' (Default)
%     'normalize'       normalize score

%    =========================================================================
% 
%  OUTPUTS
%      score     Description of input 1
%   
%
% SEE ALSO
%
% Copyright (C) 2016-2022 by Ralitsa Todorova
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.
%
%
%------------------------------------------------------------------------


if ~isvector(targetSequence),
    error('targetSequence should be a vector');
end

if numel(unique(targetSequence)) < numel(targetSequence),
    %     warning('Target sequence involves repetitions. Only the first appearance of each member will be counted.');
end

normalise = 'off';
circular = 'off';

for i = 1:2:length(varargin),
    if ~ischar(varargin{i}),
        error(['Parameter ' num2str(i+2) ' is not a property (type ''help <a href="matlab:help PairOrder">PairOrder</a>'' for details).']);
    end
    switch(lower(varargin{i})),
        case 'circular',
            circular = varargin{i+1};
            if ~isastring(circular,'on','off'),
                error('Incorrect value for property ''circular'' (type ''help <a href="matlab:help PairOrder">PairOrder</a>'' for details).');
            end
        case {'normalise','normalize'},
            normalise = varargin{i+1};
            if ~isastring(normalise,'on','off'),
                error('Incorrect value for property ''normalise'' (type ''help <a href="matlab:help PairOrder">PairOrder</a>'' for details).');
            end
        otherwise,
            error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help PairOrder">PairOrder</a>'' for details).']);
    end
end


%% Decode sequence so that it is in order (i.e. rename the neuron IDs so that they reflect the order observed in the target sequence)
% E.g. if we're comparing it will [5 1 4 9], neuron 5 needs to be named 1, etc).

% Create a dictionary-like vector 'code'. The j-th element of this vector is the position of neuron 'j' in the target sequence
code = nan(1,max([max(sequenceList(:)) max(targetSequence)]));
code(targetSequence) = Renumber(targetSequence);
code(end+1) = NaN; % add an extra element for nans
seqNoNans = sequenceList;
seqNoNans(isnan(sequenceList)) = length(code); % the NaN element just added

% The working sequence with the proper order:
m = code(seqNoNans);

% if data is circular, subtract the mean from each sequence
if strcmp(circular,'on'),
    m = 2*pi*m/max(m(:));
    m = wrap(m - repmat(circ_mean(m,[],2),1,size(m,2)));
end

score = zeros(size(sequenceList(:,1)));

% We go through each column, where the first column contains the first spike of the sequence
for i=1:size(m,2),
    % We compare the i-th spike id with the ids that come after it
    score = score+nansum(sign(m-repmat(m(:,i),1,size(m,2))),2);
    m(:,i) = nan; % we then remove the id we have just evaluated, so we only count each pair once
end


if strcmp(normalise,'off'),
    return
end

% We normalise by the maximal possible score, which is dependent on the sequence length
n = sum(~isnan(sequenceList),2);
% (n + n-1 + n-2 ... + 2 + 1 i.e. (n*(n-1)/2));
maxscore = n.*(n-1)/2;

% This estimate is too high if a repetition was present
[~,column] = meshgrid(1:size(sequenceList,2),1:size(sequenceList,1)); % sequuence ID (1 1 1 1 ... 1 1; 2 2 2 2 ... 2 2; 3.....)
occurrences = Accumulate([reshape(column(~isnan(sequenceList)),[],1) reshape(sequenceList(~isnan(sequenceList)),[],1)],1); %How many times each neuron fired in each sequence
% subtract the points that each repetition deducts.
% EXAMPLE:
% [1 2 3 4 5] gives 10 points.
% [1 2 2 4 5] gives 9 as we lose the point for 2<3: 2*(2-1)/2 = 1 point
% [1 2 2 2 5] gives 7 as we lose the points for 2<3<4: 3*(3-1)/2 = 3 points
scoreToSubtract =  occurrences.*(occurrences-1)./2;
maxscore = maxscore-sum(scoreToSubtract,2);

score = score./maxscore;

end

