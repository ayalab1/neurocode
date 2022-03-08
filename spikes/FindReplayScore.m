function [r,p,st,sp,rShuffled,aShuffled,bShuffled,c,cShuffled,jump,jumpShuffled,maxJump,maxJumpShuffled] = FindReplayScore(matrix,varargin)

%FindReplayScore 
%
%
% USAGE
%
%    [r,p,a,b,rShuffled,c,cShuffled,jump,jumpShuffled,maxJump,
%    maxJumpShuffled] = FindReplayScore(matrix,threshold,<options>);
%
%    matrix         
%    threshold      
%    
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'threshold'   considered distance from the line (default = 15 bins)
%     'nShuffles'   default = 500
%     'jumps'       compute jump distance (sum(jump)/length(windows)/length(track))
%                   and max jump distance (Silva et al.,2015) (default = 'off')
%     'wcorr'       compute weigthed correlation (cf Silva et al.,2015)
%                   (type ''help <a href="matlab:help
%                   WeightedCorrCirc">WeightedCorrCirc</a>'' for details)
%     'circular'    for circular-linear data (default = 'on')
%    =========================================================================
%
%   OUTPUT
%
%     r                 replay score of matrix
%     p                 p-value of replay score
%     st                start position bin of the fitted line
%     sp                stop position bin of the fitted line
%     rShuffled         replay scores of shuffled matrices
%     stShuffled        st of shuffled matrices
%     spShuffled        sp of shuffled matrices
%     c                 circular weigthed correlation of matrix
%     cShuffled         circular weigthed correlations of shuffled matrices
%     jump              jump value of matrix  
%     jumpShuffled      jump values of shuffled matrices
%     maxJump           max jump value of matrix
%     maxJumpShuffled   max jump values of shuffled matrices
%
% Raly & C??line ;)

% Defaults values
nShuffles = [];
p = [];
rShuffled = [];
aShuffled = [];
bShuffled = [];
cShuffled = [];
maxJumpShuffled = [];
jumpShuffled = [];
jump = [];
maxJump = [];
jumps = 'off';
wcorr = 'off';
circular = 'on';
threshold = 15;

if nargin < 1,
	error('Incorrect number of parameters (type ''help <a href="matlab:help FindReplayScore">FindReplayScore</a>'' for details).');
end

% Parse parameters
for i = 1:2:length(varargin),
	if ~ischar(varargin{i}),
		builtin('error',['Parameter ' num2str(i+2) ' is not a property (type ''help <a href="matlab:help FindReplayScore">FindReplayScore</a>'' for details).']);
	end
	switch(lower(varargin{i})),
		case 'threshold',
			threshold = varargin{i+1};
            if ~isdscalar(threshold,'>0'),
				builtin('error','Incorrect value for property ''threshold'' (type ''help <a href="matlab:help FindReplayScore">FindReplayScore</a>'' for details).');
            end
        case 'nshuffles',
			nShuffles = varargin{i+1};
            if ~isdscalar(nShuffles,'>=0'),
				builtin('error','Incorrect value for property ''nShuffles'' (type ''help <a href="matlab:help FindReplayScore">FindReplayScore</a>'' for details).');
            end
        case 'jumps',
            jumps = varargin{i+1};
            if ~isastring(jumps,'on','off'),
                builtin('error','Incorrect value for property ''jumps'' (type ''help <a href="matlab:help FindReplayScore">FindReplayScore</a>'' for details).');
            end
        case 'circular',
            circular = varargin{i+1};
            if ~isastring(circular,'on','off'),
                builtin('error','Incorrect value for property ''circular'' (type ''help <a href="matlab:help FindReplayScore">FindReplayScore</a>'' for details).');
            end
        case 'wcorr',
            wcorr = varargin{i+1};
            if ~isastring(wcorr,'on','off'),
                builtin('error','Incorrect value for property ''wcorr'' (type ''help <a href="matlab:help FindReplayScore">FindReplayScore</a>'' for details).');
            end
        otherwise,
			builtin('error',['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help FindReplayScore">FindReplayScore</a>'' for details).']);
	end
end

% Prepare some default values
nBinsY = size(matrix,1); nBinsX = size(matrix,2); 
if isempty(nShuffles),
    nShuffles = 500;
end

%% Get the matrix of sums 
% 'sums' is the score matrix for each possible (x,y) point. 
% The sum of the total score of a replay event would be the sum of the scores of all of its points
% where there are (nBinsX) points, each with its preferred X orientation (among nBinsY choices), 
% where all the points should fit on a line.

matrix = reshape(matrix,nBinsY,[]);

[x,y]= meshgrid((-threshold:1:threshold),1:nBinsY);
indices = mod(x+y-1,nBinsY)+1;
matrixID = zeros(1,1,size(matrix,2));matrixID(:) = 1:size(matrix,2);
indX = reshape(repmat(indices,1,size(matrix,2)),[nBinsY threshold*2+1 size(matrix,2)]);
ind = sub2ind(size(matrix),indX,repmat(matrixID,size(indices)));
sums = squeeze(sum(matrix(ind),2));

%% Get magic indices describing all possible lines of the sums matrix:

% a and b are the start and end angles (from 1 to nBins) of the linear fit
if strcmp(circular,'on'),
    a = reshape(repmat((1:nBinsY),nBinsY*2-1,1),[],1); % the line can start anywhere
    b = repmat(-(nBinsY-1):(nBinsY-1),1,nBinsY)' + a; % the line has to finish not farther than nBinsY-1 below or above a
else
    a = reshape(repmat((1:nBinsY),nBinsY,1),[],1); % the line can start anywhere from 1 to nBinsY
    b = reshape(repmat((1:nBinsY)',1,nBinsY),[],1); % the line can end anywhere from 1 to nBinsY
end
% make matrix of x-elements 'x'.
% the first column of x is a, and its last column is b
x = [a zeros(length(a),nBinsX-2) b];
% the middle columns go progressively from a to b
for i=2:(nBinsX-1),
    x(:,i) = x(:,1)+(i-1)*(x(:,end)-x(:,1))/(nBinsX-1);
end
% x is rounded so that it designates bin indices
x = mod(round(x)-1,nBinsY)+1;
% y-coordinates
y = repmat(1:nBinsX,size(x,1),1);
indices = sub2ind(size(sums),x,y);

scores = nanmean(sums(indices),2);
[r,ind] = max(scores);
st = a(ind); sp = b(ind);

if strcmp(wcorr,'on'),
    c = WeightedCorrCirc(matrix);
end

if strcmp(jumps,'on'),
    jumpShuffled = nan(nShuffles,1);
    maxJumpShuffled = nan(nShuffles,1);
    [~,goodWindows] = find(~isnan(matrix(1,:)));
    neighbour = diff(goodWindows);
    if any(neighbour==1),
        neighbour(neighbour~=1) = 0;
        [~,wherePmax] = max(matrix(:,goodWindows));
        d = diff(wherePmax);
        delta = mod(d(logical(neighbour)),nBinsY);
        delta(delta>nBinsY/2) = delta(delta>nBinsY/2) - nBinsY;
        maxJump = max(abs(delta));
        jump = sum(abs(delta))/length(delta);
        jumps = 'ok';
    else
        jump = nan;
        maxJump = nan;
    end
end

%% Shuffle to get a p-value

if nShuffles>0,
    rShuffled = nan(nShuffles,1);
    cShuffled = nan(nShuffles,1);
    maxJumpShuffled = nan(nShuffles,1);
    jumpShuffled = nan(nShuffles,1);
    aShuffled = nan(nShuffles,1);
    bShuffled = nan(nShuffles,1);
    for i=1:nShuffles,
        
        shift = round(rand(1,nBinsX)*(nBinsY-1));
        mockSums = CircularShift(sums,shift);
        [rShuffled(i,1),ind] = max(nanmean(mockSums(indices),2));
        aShuffled(i,1) = a(ind); bShuffled(i,1) = b(ind);

        if strcmp(wcorr,'on'),
            mockMatrix = CircularShift(matrix,shift);
            cShuffled(i,1) = WeightedCorrCirc(mockMatrix);
        end
        
        if strcmp(jumps,'ok'),
            mockMatrix = CircularShift(matrix,shift);
            [~,wherePmax] = max(mockMatrix(:,goodWindows));
            d = diff(wherePmax);
            delta = mod(d(logical(neighbour)),nBinsY);
            delta(delta>nBinsY/2) = delta(delta>nBinsY/2) - nBinsY;
            maxJumpShuffled(i,1) = max(abs(delta));
            jumpShuffled(i,1) = sum(abs(delta))/length(delta);
        end
        
    end
    p = sum(rShuffled>r)/nShuffles;    
end


% INDICES = indices;
% 
% [x,y]= meshgrid((-threshold:1:threshold),1:nBinsY);
% indices = mod(x+y-1,nBinsY)+1;
% matrixID = zeros(1,1,size(matrix,2));matrixID(:) = 1:size(matrix,2);
% indX = reshape(repmat(indices,1,size(matrix,2)),[nBinsY threshold*2+1 size(matrix,2)]);
% ind = sub2ind(size(matrix),indX,repmat(matrixID,size(indices)));
%     
% rShuffled2 = nan(nShuffles,1);
% for i=1:nShuffles,        
%     shift = round(rand(1,nBinsX)*(nBinsY-1));
%     mockMatrix = CircularShift(matrix,shift);
%     mockSums = squeeze(sum(mockMatrix(ind),2));        
%     rShuffled2(i,1) = max(nanmean(mockSums(INDICES),2));
% end

% figure;subplot(2,1,1);hist(rShuffled,50);
% subplot(2,1,2);hist(rShuffled2,50);
    
