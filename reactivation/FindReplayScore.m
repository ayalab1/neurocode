function [r,p,st,sp,rShuffled,aShuffled,bShuffled,c,cShuffled,jump,jumpShuffled,maxJump,maxJumpShuffled,quadrantScore] = FindReplayScore(matrix,varargin)

% FindReplayScore
%
%
% USAGE
%
%    [r,p,a,b,rShuffled,c,cShuffled,jump,jumpShuffled,maxJump, maxJumpShuffled,quadrantScore] = FindReplayScore(matrix,threshold,<options>);
%
%  INPUT
%
%    matrix     probability matrix for a specific event ("estimations" output of ReconstructPosition)       
%    threshold  distance away from the fitted line (in bins) to count towards the score (see Davidson et al. 2009)
%    
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'threshold'   considered distance from the line (default = 15 bins)
%     'nShuffles'   default = 500
%     'jumps'       compute jump distance and max jump distance (Silva et al.,2015)
%                   (sum(jump)/length(windows)/length(track))
%                    (default = 'off')
%     'wcorr'       compute weigthed correlation (cf Silva et al.,2015)
%                   (type ''help <a href="matlab:help
%                   WeightedCorrCirc">WeightedCorrCirc</a>'' for details)
%     'circular'    for circular-linear data (default = 'on')
%     'shuffle'     either 'column' (spatial) or 'temporal' (default = 'column');
%    =========================================================================
%
%   OUTPUT
%
%     r               replay score of matrix
%     p               p-value of replay score
%     st              start position bin of the fitted line
%     sp              stop position bin of the fitted line
%     rShuffled       replay scores of shuffled matrices
%     stShuffled      st of shuffled matrices
%     spShuffled      sp of shuffled matrices
%     c               circular weigthed correlation of matrix
%     cShuffled       circular weigthed correlations of shuffled
%                        matrices
%     jump            jump value of matrix  
%     jumpShuffled    jump values of shuffled matrices
%     maxJump         max jump value of matrix
%     maxJumpShuffled max jump values of shuffled matrices
%     quadrantScore   quadrant score (as described by Feng et al 2016)
%
%  SEE ALSO
%
% Copyright (C) 2018-2023 by Ralitsa Todorova and Celine Drieu
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.
%
%------------------------------------------------------------------------

% Defaults values
nShuffles = 500;
p = nan;
jump = nan;
maxJump = nan;
r = nan;
quadrantScore = nan;
st = nan;
sp = nan;
c = nan;
shuffle = 'column';
if nargout>7, wcorr = 'on'; else, wcorr = 'off'; end
if nargout>9, jumps = 'on'; else, jumps = 'off'; end
circular = 'off';
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
        case 'shuffle',
            shuffle = lower(varargin{i+1});
            if ~isastring(shuffle,'column','temporal'),
                builtin('error','Incorrect value for property ''shuffle'' (type ''help <a href="matlab:help FindReplayScore">FindReplayScore</a>'' for details).');
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

% Default values
nBinsY = size(matrix,1); nBinsX = size(matrix,2); 
[rShuffled, aShuffled, bShuffled, cShuffled, maxJumpShuffled, jumpShuffled] = deal(nan(1,nShuffles));

if isempty(matrix)
    return
end

%% Get the matrix of sums 
% 'sums' is the score matrix for each possible (x,y) point. 
% The sum of the total score of a replay event would be the sum of the scores of all of its points
% where there are (nBinsX) points, each with its preferred X orientation (among nBinsY choices), 
% where all the points should fit on a line.

matrix = reshape(matrix,nBinsY,[]);

[x,y]= meshgrid((-threshold:1:threshold),1:nBinsY);
if strcmp(circular,'on'),
    indices = mod(x+y-1,nBinsY)+1;
    matrixID = zeros(1,1,size(matrix,2));matrixID(:) = 1:size(matrix,2);
    indX = reshape(repmat(indices,1,size(matrix,2)),[nBinsY threshold*2+1 size(matrix,2)]);
    ind = sub2ind(size(matrix),indX,repmat(matrixID,size(indices)));
    sums = squeeze(sum(matrix(ind),2));
else
%     matrix0  = [matrix; nan(size(matrix(1,:)))]; % append nans for the cases where the linear fit goes outside of the matrix
    matrix0  = [matrix; median(matrix)]; % Silva appended the median probability instead of nans
    indices = (x+y); indices(indices<1 | indices>nBinsY) = nBinsY+1;
    matrixID = zeros(1,1,size(matrix,2));matrixID(:) = 1:size(matrix,2);
    indX = reshape(repmat(indices,1,size(matrix,2)),[nBinsY threshold*2+1 size(matrix,2)]);
    ind = sub2ind(size(matrix0),indX,repmat(matrixID,size(indices)));
    sums = squeeze(mean(matrix0(ind),2))*size(x,2);
    sums(sums>1) = 1; % make sure if peak is at the end of the matrix, the presence of nans doesn't overly inflate the possible score
end

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
x = bsxfun(@times,linspace(0,1,nBinsX),(b-a))+a;
if strcmp(circular,'on'),
    % x is rounded so that it designates bin indices
    x = mod(round(x)-1,nBinsY)+1;
else
    x = round(x);
end
indices = bsxfun(@plus, x, (0:nBinsX-1)*nBinsY);

scores = mean(sums(indices),2);
[r,ind] = max(scores);
st = a(ind); sp = b(ind);

if strcmp(wcorr,'on'),
    if strcmp(circular,'on'),
        c = WeightedCorrCirc(matrix);
    else
        c = WeightedCorr(matrix);
    end
end

if strcmp(jumps,'on'),
    jumpShuffled = nan(1,nShuffles);
    maxJumpShuffled = nan(1,nShuffles);
    [~,goodWindows] = find(~isnan(matrix(1,:)) & (max(matrix)>min(matrix) + 0.0000001));
    neighbour = diff(goodWindows);
    if any(neighbour==1),
        neighbour(neighbour~=1) = 0;
        [~,wherePmax] = max(matrix(:,goodWindows));
        d = diff(wherePmax);
        if strcmp(circular,'on'),
            delta = mod(d(logical(neighbour)),nBinsY);
            delta(delta>nBinsY/2) = delta(delta>nBinsY/2) - nBinsY;
        else
            delta = d;
        end
        maxJump = max(abs(delta));
        jump = sum(abs(delta))/length(delta);
        jumps = 'ok';
    else
        jump = nan;
        maxJump = nan;
    end
end

[nBins,nPieces] = size(matrix);
spatialBinID = (1:nBins)'/nBins; temporalBinID = linspace(0,1,nPieces);
underestimating = spatialBinID>0.375 & spatialBinID<0.5; overestimating = spatialBinID>0.5 & spatialBinID<0.625;
earlyPhase = temporalBinID>=0.2 & temporalBinID<0.5; latePhase = temporalBinID>0.5 & temporalBinID<=0.8;
Qok = false(size(matrix)); Qok(underestimating,earlyPhase) = true; Qok(overestimating,latePhase) = true;
Qcontrol = false(nBins,nPieces); Qcontrol(underestimating,latePhase) = true; Qcontrol(overestimating,earlyPhase) = true;
score = nanmean(nanmean(matrix(Qok)));
score(2) = nanmean(nanmean(matrix(Qcontrol)));
quadrantScore = (score(1)-score(2))./sum(score,2);

%% Shuffle to get a p-value

if nShuffles==0
    return
end
rShuffled = nan(1,nShuffles);
cShuffled = nan(1,nShuffles);
maxJumpShuffled = nan(1,nShuffles);
jumpShuffled = nan(1,nShuffles);
aShuffled = nan(1,nShuffles);
bShuffled = nan(1,nShuffles);
if strcmp(shuffle,'column')
    for i=1:nShuffles,
        shift = round(rand(1,nBinsX)*(nBinsY-1));
        mockSums = CircularShift(sums,shift);
        [rShuffled(i),ind] = max(mean(mockSums(indices),2));
        aShuffled(i) = a(ind); bShuffled(i) = b(ind);
        if strcmp(wcorr,'on') || strcmp(jumps,'ok')
            mockMatrix = CircularShift(matrix,shift);
            if strcmp(wcorr,'on')
                if strcmp(circular,'on'),
                    cShuffled(i) = WeightedCorrCirc(mockMatrix);
                else
                    cShuffled(i) = WeightedCorr(mockMatrix);
                end
            end
            if strcmp(jumps,'ok')
                [~,wherePmax] = max(mockMatrix(:,goodWindows));
                d = diff(wherePmax);
                if strcmp(circular,'on'),
                    delta = mod(d(logical(neighbour)),nBinsY);
                    delta(delta>nBinsY/2) = delta(delta>nBinsY/2) - nBinsY;
                else
                    delta = d;
                end
                maxJumpShuffled(i) = max(abs(delta));
                jumpShuffled(i) = sum(abs(delta))/length(delta);
            end
        end
    end
else % temporal shuffle
    for i=1:nShuffles,
        [~,order] = sort(rand(1,nBinsX));
        mockSums = sums(:,order);
        [rShuffled(i),ind] = max(mean(mockSums(indices),2));
        aShuffled(i) = a(ind); bShuffled(i) = b(ind);
        if strcmp(wcorr,'on') || strcmp(jumps,'ok')
            mockMatrix = matrix(:,order);
            if strcmp(wcorr,'on')
                if strcmp(circular,'on'),
                    cShuffled(i) = WeightedCorrCirc(mockMatrix);
                else
                    cShuffled(i) = WeightedCorr(mockMatrix);
                end
            end
            if strcmp(jumps,'ok')
                [~,wherePmax] = max(mockMatrix(:,ismember(order,goodWindows)));
                d = diff(wherePmax);
                if strcmp(circular,'on'),
                    delta = mod(d(logical(neighbour)),nBinsY);
                    delta(delta>nBinsY/2) = delta(delta>nBinsY/2) - nBinsY;
                else
                    delta = d;
                end
                maxJumpShuffled(i) = max(abs(delta));
                jumpShuffled(i) = sum(abs(delta))/length(delta);
            end
        end
    end
end

p = sum(rShuffled>=r)/nShuffles;

% ------------------------------- Helper function -------------------------------

function rho = WeightedCorrCirc(weights,x,alpha)

% Provide a matrix of weights, and this function will check the 
% correlation between the X and Y dimensions of the matrix.
% You can provide the X-values and the (angular) Y-values for 
% each column and row.
%
% Copyright (C) 2017 Céline Drieu
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

if exist('x','var') && ~isempty('x'),
    if isvector(x), x = repmat(x(:)',size(weights,1),1); end
else
    [x,~] = meshgrid(1:size(weights,2),1:size(weights,1));
end
if exist('alpha','var') && ~isempty('alpha'),
    if isvector(alpha), alpha = repmat(alpha(:),1,size(weights,2)); end
else
    [~,alpha] = meshgrid(1:size(weights,2),1:size(weights,1));
    alpha = 2*pi*alpha/size(weights,1);
end


rxs = WeightedCorr(weights,x,sin(alpha));
rxc = WeightedCorr(weights,x,cos(alpha));
rcs = WeightedCorr(weights,sin(alpha),cos(alpha));

% compute angular-linear correlation (equ. 27.47)
rho = sqrt((rxc^2 + rxs^2 - 2*rxc*rxs*rcs)/(1-rcs^2));


function c = WeightedCorr(weights,x,y)

% Provide a matrix of weights, and this function will check the 
% correlation between the X and Y dimensions of the matrix.
% You can provide the X-values and the Y-values for each column
% and row.
%
% Copyright (C) 2017 Céline Drieu
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.


if exist('x','var') && ~isempty('x'),
    if isvector(x), x = repmat(x(:)',size(weights,1),1); end
else
    [x,~] = meshgrid(1:size(weights,2),1:size(weights,1));
end
if exist('y','var') && ~isempty('y'),
    if isvector(y), y = repmat(y(:),1,size(weights,2)); end
else
    [~,y] = meshgrid(1:size(weights,2),1:size(weights,1));
end

x = x(:);
y = y(:);
w = weights(:);

mX = nansum(w.*x)./nansum(w);
mY = nansum(w.*y)./nansum(w);


covXY = nansum(w.*(x-mX).*(y-mY))./nansum(w(:));
covXX = nansum(w.*(x-mX).*(x-mX))./nansum(w(:));
covYY = nansum(w.*(y-mY).*(y-mY))./nansum(w(:));

c = covXY ./ sqrt(covXX.*covYY);
