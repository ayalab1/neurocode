function [map,stats] = Map(v,z,varargin)

% Chez raly
%Map - Map z on (x,y) where x, y and z are time-varying variables (samples).
%
%  Compute a continuous map, where one time-varying variable z is represented
%  as a function of one or two time-varying variables x and y. The variable z
%  can either be a point process (typically, a list of spike timestamps) or a
%  continuous measure (e.g. the instantaneous velocity of the animal, the
%  spectral power of an LFP channel in a given frequency band, the coherence
%  between two oscillating LFP channels, etc.) Typical examples of x and y
%  include spatial coordinates and angular directions.
%
%  An occupancy map is also computed.
%
%  USAGE
%
%    map = Map([t1 x y],[t2 z],<options>)
%
%    t1             timestamps for x and y
%    x              x values in [0,1]
%    y              optional y values in [0,1]
%    t2             timestamps for z
%    z              optional z values
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'smooth'      smoothing size in bins (0 = no smoothing, default = 2)
%     'nBins'       number of horizontal and vertical bins (default = [50 50])
%     'minTime'     minimum time spent in each bin (in s, default = 0)
%     'mode'        'interpolate' to interpolate missing points (< minTime)
%                   (default for a continuous z) or 'discard' to discard them
%                   (default for a point process)
%     'maxDistance' maximal distance for interpolation (default = 5)
%     'maxGap'      z values recorded during time gaps between successive (x,y)
%                   samples exceeding this threshold (e.g. undetects) will not
%                   be interpolated; also, such long gaps in (x,y) sampling
%                   will be clipped to 'maxGap' to compute the occupancy map
%                   (default = 0.100 s)
%     'type'        three letters (one for X, one for Y and one for Z) indi-
%                   cating which coordinates are linear ('l') and which are
%                   circular ('c') - for 1D data, only two letters are used
%                   (default 'lll')
%     'resolution'  [maxX maxY] positions resolution (default = [0 1])
%    =========================================================================
%
%  OUTPUT
%
%    map.x          x bins
%    map.y          y bins
%    map.z          average map (z continuous)
%                   or rate map (z point process)
%    map.count      count map (z point process)
%    map.time       occupancy map (in s)
%
%  NOTES
%
%    x values are arranged in columns and y values in rows in all output matrices
%    (e.g. 'map.z').
%
%  SEE
%
%    See also MapStats, FiringMap, PlotColorMap, Accumulate.

% Copyright (C) 2002-2014 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Check number of parameters
if nargin < 2 | mod(length(varargin),2) ~= 0,
    error('Incorrect number of parameters (type ''help <a href="matlab:help Map">Map</a>'' for details).');
end

% Check parameter sizes
if size(v,2) < 2,
    error('Parameter ''[t x y]'' should have at least 2 columns (type ''help <a href="matlab:help Map">Map</a>'' for details).');
end
if (size(z,2) < 1 || size(z,2) > 2) && ~isempty(z),
    error('Parameter ''z'' should have 1 or 2 columns (type ''help <a href="matlab:help Map">Map</a>'' for details).');
end

% Default values
maxGap = 0.1;
map.x = [];
map.y = [];
map.count = [];
map.time = [];
map.z = [];
smooth = 2;
nBins = 50;
minTime = 0;
type = 'lll';
mode = 'discard';
maxDistance = 5;
resolution = [];

if isempty(v) || size(v,1) < 2, return; end

% Some info about x, y and z
pointProcess = (isempty(z) | size(z,2) == 1);
if ~pointProcess, mode = 'interpolate'; end
t = v(:,1);
x = v(:,2);
if size(v,2) >= 3,
    y = v(:,3);
else
    y = [];
end

% Parse parameter list
for i = 1:2:length(varargin),
    if ~ischar(varargin{i}),
        error(['Parameter ' num2str(i+2) ' is not a property (type ''help <a href="matlab:help Map">Map</a>'' for details).']);
    end
    switch(lower(varargin{i})),
        
        case 'smooth',
            smooth = varargin{i+1};
            if ~isdvector(smooth,'>=0') || length(smooth) > 2,
                error('Incorrect value for property ''smooth'' (type ''help <a href="matlab:help Map">Map</a>'' for details).');
            end
            
        case 'nbins',
            nBins = varargin{i+1};
            if ~isivector(nBins,'>0') || length(nBins) > 2,
                error('Incorrect value for property ''nBins'' (type ''help <a href="matlab:help Map">Map</a>'' for details).');
            end
            
        case 'mintime',
            minTime = varargin{i+1};
            if ~isdscalar(minTime,'>=0'),
                error('Incorrect value for property ''minTime'' (type ''help <a href="matlab:help Map">Map</a>'' for details).');
            end
            
        case 'maxgap',
            maxGap = varargin{i+1};
            if ~isdscalar(maxGap,'>=0'),
                error('Incorrect value for property ''maxGap'' (type ''help <a href="matlab:help Map">Map</a>'' for details).');
            end
            
        case 'maxdistance',
            maxDistance = varargin{i+1};
            if ~isdscalar(maxDistance,'>=0'),
                error('Incorrect value for property ''maxDistance'' (type ''help <a href="matlab:help Map">Map</a>'' for details).');
            end
            
        case 'mode',
            mode = lower(varargin{i+1});
            if ~isastring(mode,'interpolate','discard'),
                error('Incorrect value for property ''mode'' (type ''help <a href="matlab:help Map">Map</a>'' for details).');
            end
            
        case 'resolution',
            resolution = lower(varargin{i+1});
            if  ~isdvector(resolution,'>=0') || length(resolution) > 2
                error('Incorrect value for property ''resolution'' (type ''help <a href="matlab:help Map">Map</a>'' for details).');
            end
            
        case 'type',
            type = lower(varargin{i+1});
            if (isempty(y) && ~isastring(type,'cc','cl','lc','ll')) || (~isempty(y) && ~isastring(type,'ccl','cll','lcl','lll','ccc','clc','lcc','llc')),
                error('Incorrect value for property ''type'' (type ''help <a href="matlab:help Map">Map</a>'' for details).');
            end
            
        otherwise,
            error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help Map">Map</a>'' for details).']);
            
    end
end

% Make sure x and y are normalized
if max(x) > 1 || min(x) < 0 || ~isempty(resolution),
    if isempty(resolution),
        xLimits = [min(x) max(x)];
    else
        xLimits = [0 resolution(1)]; 
    end
    x = (x-xLimits(1))/diff(xLimits);
else
    xLimits = [0 1];
end
if ~isempty(y),
    if max(y) > 1 || min(y) < 0 || ~isempty(resolution),
        if isempty(resolution),
            yLimits = [min(y) max(y)];
        else
            yLimits = [0 resolution(2)]; 
        end
        y = (y-yLimits(1))/diff(yLimits);
    else
        yLimits = [0 1];
    end
end

% Number of bins for x and y
nBinsX = nBins(1);
if length(nBins) == 1,
    nBinsY = nBinsX;
    nBins(2) = nBins;
else
    nBinsY = nBins(2);
end

% Bin x and y
x = Bin(x,[0 1],nBinsX);
if ~isempty(y),
    y = Bin(y,[0 1],nBinsY);
end

% Duration for each (X,Y) sample (clipped to maxGap)
dt = diff(t);dt(end+1)=dt(end);dt(dt>maxGap) = maxGap;

if pointProcess,
    % Count occurrences for each (x,y) timestamp
    n = CountInIntervals(z,[t t+dt]);
else
    % Interpolate z at (x,y) timestamps
    [z,discarded] = Interpolate(z,t,'maxGap',maxGap);
    if isempty(z), return; end
    if strcmp(type(end),'c'),
        range = isradians(z(:,2));
        z(:,2) = exp(j*z(:,2));
    end
%   make sure z and t have the same number of rows (remove any discarded values)
    if any(discarded),
        x(discarded,:) = []; dt(discarded,:) = [];
        if ~isempty(y),y(discarded) = []; end
    end
    n = 1;
end

% Computations
if isempty(y),
    % 1D (only x)
    map.x = linspace(xLimits(1),xLimits(2),nBinsX);
    map.count = Accumulate(x,n,nBinsX);
    map.time = Accumulate(x,dt,nBinsX);
    valid = map.time > minTime;
    map.count = Smooth(Interpolate1(map.x,map.count,valid,mode,maxDistance),smooth,'type',type(1))';
    map.time = Smooth(Interpolate1(map.x,map.time,valid,mode,maxDistance),smooth,'type',type(1))';
    if pointProcess,
        map.z = map.count./(map.time+eps);
    else
        map.z = Accumulate(x,z(:,2),'size',nBinsX,'mode','mean');
        map.z = Smooth(Interpolate1(map.x',map.z,valid,mode,maxDistance),smooth,'type',type(1))';
    end
else
    % 2D (x and y)
    map.x = linspace(xLimits(1),xLimits(2),nBinsX);
    map.y = linspace(yLimits(1),yLimits(2),nBinsY);
    map.count = Accumulate([x y],n,nBins);
    map.time = Accumulate([x y],dt,nBins);
    valid = map.time > minTime;
    map.count = Smooth(Interpolate2(map.x,map.y,map.count,valid,mode,maxDistance),smooth,'type',type(1:2))';
	map.time = Smooth(Interpolate2(map.x,map.y,map.time,valid,mode,maxDistance),smooth,'type',type(1:2))';
    if pointProcess,
        map.z = map.count./(map.time+eps);
    else
        map.z = Accumulate([x y],z(:,2),'size',nBins);
        map.z = Interpolate2(map.x,map.y,map.z,valid,mode,maxDistance)';
        nans = isnan(map.z);
        map.z = nansmooth(map.z,smooth,'type',type(1:2));
        map.z(nans) = nan;
        map.z = map.z./(map.count+eps);
    end
end

% Circular z
if ~exist('range','var'),range=1;end
if strcmp(type(end),'c'), map.z = wrap(angle(map.z),range); end


% Interpolate or discard regions with insufficient sampling
if strcmp(mode,'discard'),
    map.z(map.time<=minTime) = 0;
end

% ------------------------------- Helper functions -------------------------------

% Interpolate if required (1D)
function yint = Interpolate1(x,y,valid,mode,maxDistance)

if strcmp(mode,'discard'),
	yint = y;
else
	yint = interp1(x(valid),y(valid),x);
end

% Interpolate if required (2D)
function zint = Interpolate2(x,y,z,valid,mode,maxDistance)

if strcmp(mode,'discard'),
	% In discard mode, do nothing
	zint = z;
else
	% In interpolation mode, interpolate missing points (where time < minTime) using other points
	d = DistanceTransform(valid);
    xx = repmat(x',1,length(y));
    yy = repmat(y,length(x),1);
	if exist('scatteredInterpolant') == 2,
		F = scatteredInterpolant(xx(d==0),yy(d==0),z(d==0));
        zint = F(xx,yy);
    else
        if any(imag(z(:))),
            Freal = TriScatteredInterp(xx(d==0),yy(d==0),real(z(d==0)));
            zintReal = Freal(xx,yy);
            Fimaginary = TriScatteredInterp(xx(d==0),yy(d==0),imag(z(d==0)));
            zintImaginary = Fimaginary(xx,yy);
            zint = zintReal + 1i.*zintImaginary;
        else
            F = TriScatteredInterp(xx(d==0),yy(d==0),z(d==0));
            zint = F(xx,yy);
        end
    end
	% (do not interpolate missing points too distant from valid points)
	zint(d>maxDistance) = z(d>maxDistance);
	zint(isnan(zint)) = z(isnan(zint));
end

