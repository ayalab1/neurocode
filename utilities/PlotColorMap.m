function PlotColorMap(data,dimm,varargin)

%PlotColorMap - Plot a color map.
%
%  Plot a color map (e.g. the firing field of a place cell).
%
%  USAGE
%
%    PlotColorMap(data,dimm,<options>)
%
%    data           data
%    dimm           optional luminance map
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'newfig'	    'on'/'off' (default = 'on')
%     'x'           abscissae
%     'y'           ordinates
%     'threshold'   dimm values below this limit are zeroed (default = 0.01)
%     'cutoffs'     lower and upper cutoff values ([] = autoscale, default)
%     'hgamma'      gamma-like correction for hue (1 = no correction, default)
%     'colorspace'  'HSV' or 'RGB' (default = 'RGB')
%     'bar'         draw a color bar (default = 'off')
%     'type'        either 'linear' or 'circular' (default 'linear')
%     'ydir'        either 'normal' (default) or 'reverse' (useful when the
%                   x and y coordinates correspond to spatial positions,
%                   as video cameras measure y in reverse direction)
%    =========================================================================
%
%  NOTE
%
%    The luminance map is used to dimm the color map. A single scalar value
%    is interpreted as a constant luminance map. If this parameter is not
%    provided, normal equiluminance is assumed (i.e. scalar value of 1).
%
%  EXAMPLE
%
%    fm = FiringMap(positions,spikes);      % firing map for a place cell
%    figure;PlotColorMap(fm.rate,fm.time);  % plot, dimming with occupancy map
%
%  SEE
%
%    See also FiringMap, PhaseMap, MTSpectrogram, PlotShortTimeCCG.

% Copyright (C) 2004-2012 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Default values
newfig='off';
cutoffs = [];
hgamma = 1;
gamma = 1;
hg = 0;
colorspace = 'rgb';
threshold = 0.01;
drawBar = 0;
type = 'linear';
[y,x] = size(data);
x = 1:x;y = 1:y;
ydir = 'normal';

if nargin < 1,
	error('Incorrect number of parameters (type ''help <a href="matlab:help PlotColorMap">PlotColorMap</a>'' for details).');
end
if nargin == 1,
	dimm = 1;
end
if isa(dimm,'char'),
	varargin = {dimm varargin{:}};
	dimm = 1;
end

% Parse parameter list
for i = 1:2:length(varargin),
	if ~ischar(varargin{i}),
		error(['Parameter ' num2str(i+2) ' is not a property (type ''help <a href="matlab:help PlotColorMap">PlotColorMap</a>'' for details).']);
	end
	switch(lower(varargin{i})),

		case 'threshold',
			threshold = varargin{i+1};
			if ~isdscalar(threshold,'>=0'),
				error('Incorrect value for property ''threshold'' (type ''help <a href="matlab:help PlotColorMap">PlotColorMap</a>'' for details).');
			end
		case 'x',
			x = varargin{i+1};
			if ~isdvector(x),
				error('Incorrect value for property ''x'' (type ''help <a href="matlab:help PlotColorMap">PlotColorMap</a>'' for details).');
			end
		case 'y',
			y = varargin{i+1};
			if ~isdvector(y),
				error('Incorrect value for property ''y'' (type ''help <a href="matlab:help PlotColorMap">PlotColorMap</a>'' for details).');
			end
		case 'cutoffs',
			cutoffs = varargin{i+1};
			if ~isdvector(cutoffs,'#2','<'),
				error('Incorrect value for property ''cutoffs'' (type ''help <a href="matlab:help PlotColorMap">PlotColorMap</a>'' for details).');
			end
		case 'hgamma',
			hg = 1;
			hgamma = varargin{i+1};
			if ~isdscalar(hgamma,'>=0'),
				error('Incorrect value for property ''hgamma'' (type ''help <a href="matlab:help PlotColorMap">PlotColorMap</a>'' for details).');
			end
		case 'gamma',
			gamma = varargin{i+1};
			if ~isdscalar(gamma,'>=0'),
				error('Incorrect value for deprecated property ''gamma'' (type ''help <a href="matlab:help PlotColorMap">PlotColorMap</a>'' for details).');
			end
		case 'colorspace',
			colorspace = lower(varargin{i+1});
			if ~isstring_FMAT(colorspace,'hsv','rgb'),
				error('Incorrect value for property ''colorspace'' (type ''help <a href="matlab:help PlotColorMap">PlotColorMap</a>'' for details).');
			end
		case 'bar',
			drawBar = lower(varargin{i+1});
			if ~isstring_FMAT(drawBar,'on','off'),
				error('Incorrect value for property ''bar'' (type ''help <a href="matlab:help PlotColorMap">PlotColorMap</a>'' for details).');
			end
		case 'type',
			type = lower(varargin{i+1});
			if ~isstring_FMAT(type,'linear','circular'),
				error('Incorrect value for property ''type'' (type ''help <a href="matlab:help PlotColorMap">PlotColorMap</a>'' for details).');
			end
		case 'ydir',
			ydir = lower(varargin{i+1});
			if ~isstring_FMAT(ydir,'normal','reverse'),
				error('Incorrect value for property ''ydir'' (type ''help <a href="matlab:help PlotColorMap">PlotColorMap</a>'' for details).');
			end
		case 'newfig'
			newfig = lower(varargin{i+1});
			if ~isstring_FMAT(newfig,'on','off'),
				error('Incorrect value for property ''newfig'' (type ''help <a href="matlab:help PlotColorMap">PlotColorMap</a>'' for details).');
			end
		otherwise,
			error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help PlotColorMap">PlotColorMap</a>'' for details).']);

	end
end

if hg == 0,
	hgamma = 1/gamma;
end

if ~isempty(cutoffs),
	m = cutoffs(1);
	M = cutoffs(2);
else
	m = min(min(data));
	M = max(prctile(data,90));
end
if m == M, M = m+1; end
if isnan(m), m = 0; M = 1; end

if length(dimm) == 1,
	dimm = dimm*ones(size(data));
end

if strcmp(newfig,'on')
  f = figure;
end

p = imagesc(x,y,data,[m M]);
set(gca,'color',[1 1 1]);

if any(dimm~=1),
 	alpha(p,1./(1+threshold./(dimm+eps)));
end
set(gca,'ydir',ydir,'tickdir','out','box','off');
if ~isempty(x) && length(x) ~= 1,
	binX = (max(x)-min(x))/length(x);
	set(gca,'xlim',[min(x)-binX/2 max(x)+binX/2]);
end
if ~isempty(y) && length(y) ~= 1,
	binY = (max(y)-min(y))/length(y);
	set(gca,'ylim',[min(y)-binY/2 max(y)+binY/2]);
end

colormap(Bright(100,'hgamma',hgamma,'type',type));

if strcmp(drawBar,'on'),
	b = colorbar('vert');
% 	set(b,'xtick',[],'tickdir','out','box','off','ycolor','k');
	set(gcf,'currentaxes',gca);
end