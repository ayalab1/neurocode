function [pieces ids] = SplitIntervals(intervals,varargin)

% Splits intervals into even pieces. There are two different ways to call this
% function: by either providing 'pieceSize' or 'nPieces'.
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

pieceSize = 0.02;
nPieces = [];

for i = 1:2:length(varargin),
    if ~ischar(varargin{i}),
        error(['Parameter ' num2str(i+2) ' is not a property (type ''help <a href="matlab:help SplitIntervals">SplitIntervals</a>'' for details).']);
    end
    switch(lower(varargin{i})),
        case 'piecesize',
            pieceSize = varargin{i+1};
            if ~isvector(pieceSize) || length(pieceSize) ~= 1,
                error('Incorrect value for property ''pieceSize'' (type ''help <a href="matlab:help SplitIntervals">SplitIntervals</a>'' for details).');
            end
        case 'npieces',
            nPieces = varargin{i+1};
            if ~isvector(nPieces) || length(nPieces) ~= 1,
                error('Incorrect value for property ''nPieces'' (type ''help <a href="matlab:help SplitIntervals">SplitIntervals</a>'' for details).');
            end
        otherwise,
            error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help SplitIntervals">SplitIntervals</a>'' for details).']);
    end
end

%% SPLIT EACH INTERVAL INTO EQUAL PIECES (nPieces)
if ~isempty(nPieces),
    matrix = nan(size(intervals,1),nPieces+1);
    matrix(:,[1 end]) = intervals;
    for i=2:(nPieces),
        matrix(:,i) = matrix(:,1)+(i-1)*(matrix(:,end)-matrix(:,1))/nPieces;
    end
    pieces = [reshape(matrix(:,1:end-1)',[],1) reshape(matrix(:,2:end)',[],1)];
    ids = ceil((1:size(pieces,1))'/nPieces);
    return
end

%% SPLIT INTERVALS INTO PIECES OF EQUAL ABSULUTE LENGTHS (pieceSize)
d = diff(intervals,[],2);
piecesPerInterval = floor(d./pieceSize);
firstPiecePosition = CumSum([1;piecesPerInterval(1:end-1)]);
% create the pieces by filling in the interval starts in their respective positions
pieces = zeros(sum(piecesPerInterval),1);
pieces(firstPiecePosition,1) = intervals(:,1);
reset = pieces>0; reset(1)=1;
time2add = ones(size(pieces))*pieceSize;
time2add = CumSum(time2add,reset) - pieceSize; %reset at new bin starts and make sure baseline is 0
if pieceSize<1, time2add = round(time2add/pieceSize)*pieceSize; end %fix underflow errers: time2add should be divisible by pieceSize
pieces = CumSum(pieces,reset) + time2add; %repeat first pieces for the duration of the interval
pieces(:,2) = pieces(:,1)+pieceSize;
ids = CumSum(reset);

end


