function [pieces ids] = SplitIntervals(intervals,varargin)

%SplitIntervals - split intervals into even pieces. 
% Specify either 'pieceSize' or 'nPieces'.
% EXAMPLE 1 (with pieceSize): 
% SplitIntervals([2 5; 6 10], 'pieceSize', 2) (i.e. split these intervals into 2-second pieces)
% gives intervals = [2 4; 6 8; 8 10] and ids = [1; 2; 2]
% EXAMPLE 2 (with nPieces): 
% SplitIntervals([2 5; 6 10], 'nPieces', 2) (i.e. split each of these intervals into 2)
% gives intervals = [2 3.5; 3.5 5; 6 8; 8 0] and ids = [1; 1; 2; 2]
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'pieceSize'   the required size of the split intervals
%     'nPieces'     alternatively, one can provide the required number of pieces.
%                   Thus the piece size is dictated by the overall duration of 
%                   the intervals, divided by nPieces.
%     'mode'        in case an interval is not equally divisible into pieces of
%                   the desired pieceSize, the function can behave in three
%                   possible ways: 'discard', 'keep', and 'round'. If the last 
%                   piece at the end of an interval is shorter than "pieceSize", 
%                   it will be discarded in all cases ('discard'), kept in all
%                   cases ('keep'), or kept only if it is longer than half the
%                   desired duration ('round'). (default = 'discard'). 
%     'extend'      When the function is called with a desired 'pieceSize' in
%                   mode 'keep' or 'half', this option will round up the 
%                   end of the last piece so that it matches the duration of
%                   all the other pieces (pieceSize). (default = false). 
%    ===========================================================================
%
% OUTPUTS
%
%    pieces         resulting intervals (split)
%    ids            for each final piece, the index of the original interval
%                   to which it belongs
%    =========================================================================
% Hint: to get overlapping windows, simply use
% sortrows([SplitIntervals(intervals, window); SplitIntervals(intervals+window/2, window)])
%
% EXAMPLE 1 (with nPieces): 
% SplitIntervals([2 5; 6 10], 'nPieces', 2) (i.e. split each of these intervals into 2)
% gives intervals = [2 3.5; 3.5 5; 6 8; 8 0] and ids = [1; 1; 2; 2]
%
% EXAMPLE 2 (with pieceSize): 
% SplitIntervals([2 5; 6 10], 'pieceSize', 2) (i.e. split these intervals into 2-second pieces)
% gives intervals = [2 4; 6 8; 8 10] and ids = [1; 2; 2]
% SplitIntervals([0 2.1],'pieceSize',1,'mode','discard') generates [0 1; 1 2], discarding the
% bin [2 2.1] which is too short (desired duration is 1 second).
% SplitIntervals([0 2.1],'pieceSize',1,'mode','keep') generates [0 1; 1 2; 2 2.1]
% SplitIntervals([0 2.1],'pieceSize',1,'mode','keep','extend',true) generates [0 1; 1 2; 2 3]
% extending the last piece so that it is also 1 second long.
%
% Copyright (C) 2019-2023 Ralitsa Todorova & (C) 2023 by Federica Lareno Faccini
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

pieceSize = 0.02;
mode = 'discard';
extend = false;
nPieces = [];

for i = 1:2:length(varargin)
    if ~ischar(varargin{i})
        error(['Parameter ' num2str(i+2) ' is not a property (type ''help <a href="matlab:help SplitIntervals">SplitIntervals</a>'' for details).']);
    end
    switch(lower(varargin{i}))
        case 'piecesize'
            pieceSize = varargin{i+1};
            if ~isvector(pieceSize) || length(pieceSize) ~= 1
                error('Incorrect value for property ''pieceSize'' (type ''help <a href="matlab:help SplitIntervals">SplitIntervals</a>'' for details).');
            end
        case 'mode'
            mode = lower(varargin{i+1});
            if ~isastring(lower(mode),'discard','keep','round')
                error('Incorrect value for property ''mode'' (type ''help <a href="matlab:help SplitIntervals">SplitIntervals</a>'' for details).');
            end
        case 'npieces'
            nPieces = varargin{i+1};
            if ~isvector(nPieces) || length(nPieces) ~= 1
                error('Incorrect value for property ''nPieces'' (type ''help <a href="matlab:help SplitIntervals">SplitIntervals</a>'' for details).');
            end
        case 'extend'
            extend = varargin{i+1};
            if isastring(lower(extend),'on','off')
                extend = strcmpi(extend,'on'); % transform to logical
            end
            if length(extend)==1
                if ~islogical(extend), extend = extend>0; end % allow for 0/1 usage and transform it to true/false
            else
                error('Incorrect value for property ''discard'' (type ''help <a href="matlab:help SplitIntervals">SplitIntervals</a>'' for details).');
            end
        otherwise
            error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help SplitIntervals">SplitIntervals</a>'' for details).']);
    end
end

%% Split each interval into we equal pieces (nPieces)

if ~isempty(nPieces)
    matrix = nan(size(intervals,1),nPieces+1);
    matrix(:,[1 end]) = intervals;
    for i=2:(nPieces)
        matrix(:,i) = matrix(:,1)+(i-1)*(matrix(:,end)-matrix(:,1))/nPieces;
    end
    pieces = [reshape(matrix(:,1:end-1)',[],1) reshape(matrix(:,2:end)',[],1)];
    ids = ceil((1:size(pieces,1))'/nPieces);
    return
end

%% Split intervals into pieces of equal durations (pieceSize)

d = diff(intervals,[],2);
switch mode
    case 'discard'
        piecesPerInterval = floor(d./pieceSize);
    case 'keep'
        piecesPerInterval = ceil(d./pieceSize);
    case 'round'
        piecesPerInterval = round(d./pieceSize);
end
if ~any(piecesPerInterval>0), pieces = zeros(0,2); ids = zeros(0,1); return; end
indicesNotEmpty = find(piecesPerInterval>0); try indicesNotEmpty(end+1) = indicesNotEmpty(end); end
firstPiecePosition = CumSum([1;piecesPerInterval(1:end-1)]);
% create the pieces by filling in the interval starts in their respective positions
pieces = zeros(sum(piecesPerInterval),1);
pieces(firstPiecePosition,1) = intervals(:,1);
reset = pieces>0; reset(1)=1;
time2add = ones(size(pieces))*pieceSize;
time2add = CumSum(time2add,reset) - pieceSize; %reset at new bin starts and make sure baseline is 0
ids = CumSum(reset);
if pieceSize<1, time2add = round(time2add/pieceSize)*pieceSize; end %fix underflow errers: time2add should be divisible by pieceSize
pieces = CumSum(pieces,reset) + time2add; %repeat first pieces for the duration of the interval
pieces(:,2) = pieces(:,1)+pieceSize;
if any(strcmpi(mode,{'keep','round'})) && ~extend
    overflow = pieces(:,2)>intervals(ids,2);
    pieces(overflow,2) = intervals(ids(overflow),2);
end

try
    ids = indicesNotEmpty(ids);
catch
    warning('There is an issue with the indices. Check results manually');
end


end


