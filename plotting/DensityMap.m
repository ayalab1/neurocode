function varargout = DensityMap(x, y, varargin)

%DensityMap - colormap equivalent of the dotplot produced by PlotXY 
%
% This funtion produces a colormap equivalent of the dotplot producesd by
% PlotXY. Provide two vectors x and y (of equal lengths) to get a map with
% values of x varying in different columns (abscissa) and values of y
% increasing down the rows (ordinate). The values in the matrix correspond
% to the probability of encountering the corresponding x and y values
% together.
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'cutoffs'     lower and upper cutoff values ([] = autoscale, default)
%     'output'      'full' to get the full structure 'map', or 'partial' 
%                   (default) for just a count (density) map. 
%     'nBins'       number of bins (rows and columns) to use for the map
%     'smooth'      default smoothing is 2
%     'show'        default = 'on';
%
% Copyright (C) 2016 by Ralitsa Todorova
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.


% Default values:
cutoffs = [];
output = 'partial';
nBins = 50;
show = 'on';
smooth = 2;
type = 'lll';

if length(x)<2,
    warning('Color map cannot be produced on the basis of the single value (vector is required). Returning value.')
    map = x;
    return
end

for i = 1:2:length(varargin),
    if ~ischar(varargin{i}),
        error(['Parameter ' num2str(i+2) ' is not a property (type ''help <a href="matlab:help DensityMap">DensityMap</a>'' for details).']);
    end
    switch(lower(varargin{i})),
        case 'output',
            output = varargin{i+1};
            if ~isstring(output,'full','partial'),
                error('Incorrect value for property ''output'' (type ''help <a href="matlab:help DensityMap">DensityMap</a>'' for details).');
            end
        case 'cutoffs',
            cutoffs = varargin{i+1};
            if ~isvector(cutoffs) || length(cutoffs) ~= 2,
                error('Incorrect value for property ''cutoffs'' (type ''help <a href="matlab:help DensityMap">DensityMap</a>'' for details).');
            end
        case 'nbins',
            nBins = varargin{i+1};
            if ~isvector(nBins) || length(nBins) > 2,
                error('Incorrect value for property ''nBins'' (type ''help <a href="matlab:help DensityMap">DensityMap</a>'' for details).');
            end
        case 'show',
            show = varargin{i+1};
            if ~isastring(show,'on','off'),
                error('Incorrect value for property ''show'' (type ''help <a href="matlab:help DensityMap">DensityMap</a>'' for details).');
            end
        case 'smooth',
	  smooth = varargin{i+1};
	  if ~isvector(smooth) || length(smooth) > 2,
	      error('Incorrect value for property ''smooth'' (type ''help <a href="matlab:help DensityMap">DensityMap</a>'' for details).');
	  end
        case 'type',
	  type = lower(varargin{i+1});
	  if (isempty(y) && ~isastring(type,'cc','cl','lc','ll')) || (~isempty(y) && ~isastring(type,'ccl','cll','lcl','lll','ccc','clc','lcc','llc')),
	      error('Incorrect value for property ''type'' (type ''help <a href="matlab:help Map">DensityMap</a>'' for details).');
	  end
        otherwise,
	  error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help DensityMap">DensityMap</a>'' for details).']);
    end
end


if size(x,2) == 2 %if index column is provided
    t(:,1) = x(:,1);
    x = x(:,2);
else
    t(:,1) = 1:length(x);
end


if ~exist('y', 'var'),
    y=x;
    x = t;
end

if length(y)~= length(t)
    error('Lengths of x and y need to be the same');
end

a=x; b=y;
% Transforming to conform with the criteria for 'Map'
% a = (x-min(x))/(max(x)-min(x));
% b = (y-min(y))/(max(y)-min(y));

% Get rid of nans if any:
ok = abs(a)<Inf & abs(b)<Inf; a = a(ok); b = b(ok); t = t(ok);

map = Map0([t a b], t, 'smooth', smooth, 'nBins', nBins,'type',type);

% Reversing the transform

% map.x = map.x*(max(x) - min(x)) + min(x);
% map.y = map.y*(max(y) - min(y)) + min(y);


if strcmpi(show,'on')
    PlotColorMap(map.count, 'x', map.x, 'y', map.y, 'cutoffs', cutoffs);
end



if strcmpi(output,'partial')
    x = map.x;
    y = map.y;
    map = map.count;
end

% if nargout==0,
%     map = 'you forgot the semicolon ;)';
% end
if nargout>0
    varargout{1} = map;
    varargout{2} = x;
    varargout{3} = y;
end

