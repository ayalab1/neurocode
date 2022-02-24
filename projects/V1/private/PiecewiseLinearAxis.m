function a = PiecewiseLinearAxis(x,axis,varargin)

%PiecewiseLinearAxis - Label axis using piecewise linear values.
%
%  USAGE
%
%    a = PiecewiseLinearAxis(x,axis,<options>)
%
%    x              axis values
%    axis           optional axis: 'x' or 'y' (default = 'x')
%    <options>      options for function <a href="matlab:help plot">plot</a>

% Copyright (C) 2015 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Check number of parameters
if nargin < 1,
  error('Incorrect number of parameters (type ''help <a href="matlab:help PiecewiseLinearAxis">PiecewiseLinearAxis</a>'' for details).');
end

% Check parameters
if ~isdvector(x,'<'),
  error('Incorrect values (type ''help <a href="matlab:help PiecewiseLinearAxis">PiecewiseLinearAxis</a>'' for details).');
end
if nargin < 2,
	axis = 'x';
else
	axis = lower(axis);
end
if ~isastring(axis,'x','y'),
	varargin = {axis,varargin{:}};
	axis = 'x';
end
if isempty(varargin)
	varargin = {'k:'};
end

% Find transitions (gaps) between linear pieces, and show them graphically
dx = diff(x);
binSize = median(dx);
gaps = [find(dx>2*binSize)+1];
a = gca;
set(a,[axis 'lim'],[min(x)-binSize/2 max(x)+binSize/2]);
if isempty(gaps), return; end

% Because Matlab uses linear axes, we will need a reference vector X starting at x(1), ending at x(end),
% and of the same length as x. This is how Matlab really views the values on the axis.
X = linspace(x(1),x(end),length(x))';

% Graphically indicate gap locations
if strcmp(axis,'x'),
	PlotHVLines(X(gaps),varargin{:});
else
	PlotHVLines(X(gaps),'h',varargin{:});
end

% We will ask Matlab to set optimal ticks for each linear piece
% To this end, we create temporary, invisible axes of the same physical size (in the figure window)
% as each piece, set its limits, and get tick information from there; we then transfer these back
% to the original axis
P = get(a,'position');
gaps = [1;gaps;length(x)];
ticks = [];
tickLabels = [];
for i = 2:length(gaps),
	% Determine location of start/stop of current piece, in Matlab 'coordinates' (i.e. using X)
	start = X(gaps(i-1));
	stop = X(gaps(i));
	% Create temporary, invisible axes, and set their physical width (in the figure window) so that
	% it matches that of current linear piece
	if strcmp(axis,'x'),
		factor = P(3)/(x(end)-x(1))
		p = [0 0.25 (stop-start)*factor 0.5];
	else
		factor = P(4)/(x(end)-x(1));
		p = [0.25 0 0.5 (stop-start)*factor];
	end
	tmp = axes('position',p,'visible','off');
	% Set limits as the beginning and end of current piece
	set(tmp,[axis 'lim'],[x(gaps(i-1)) x(gaps(i)-1)]);
	% Get locations of ticks computed by Matlab
	l = str2num(get(tmp,[axis 'ticklabel']));
	% Translate into original axis coordinates (X)
	t = interp1(x,X,l);
	% Store and move to next piece
	tickLabels = [tickLabels;l];
	ticks = [ticks;t];
	delete(tmp);
end
% Set tick locations and labels
set(a,[axis 'tick'],ticks,[axis 'TickLabel'],tickLabels);
