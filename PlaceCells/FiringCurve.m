function [curve,stats] = FiringCurve(samples,spikes,varargin)

%FiringCurve - Compute firing curve (e.g. for a head direction cell).
%
%  Compute firing curve (e.g. for a head direction cell, or a place or grid cell
%  on a linear track), as well as occupancy and spike count curves. Optionnally,
%  curve statistics can also be computed, including in-field peak and mean firing
%  rates, firing field width, etc.
%
%  USAGE
%
%    [curve,stats] = FiringCurve(samples,spikes,<options>)
%
%    samples        angular or linear <a href="matlab:help samples">samples</a>
%    spikes         spike timestamps
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'smooth'      smoothing size in bins (0 = no smoothing, default = 2)
%     'nBins'       number of bins (default = 50)
%     'minTime'     minimum time spent in each bin (in s, default = 0)
%     'maxGap'      z values recorded during time gaps between successive (x,y)
%                   samples exceeding this threshold (e.g. undetects) will not
%                   be interpolated; also, such long gaps in (x,y) sampling
%                   will be clipped to 'maxGap' to compute the occupancy map
%                   (default = 0.100 s)
%     'type'        'linear' for linear data, 'circular' for angular data
%                   (default 'linear')
%     'threshold'   values above threshold*peak belong to the field
%                   (default = 0.2)
%     'minSize'     fields smaller than this size are considered spurious
%                   and ignored (default = 10)
%     'minPeak'     peaks smaller than this size are considered spurious
%                   and ignored (default = 1)
%     'mode'        'interpolate' to interpolate missing points (< minTime),
%                   or 'discard' to discard them (default)
%     'maxDistance' maximal distance for interpolation (default = 5)
%     'verbose'     display processing information (default = 'off')
%    =========================================================================
%
%  OUTPUT
%
%    curve.x             x bins
%    curve.rate          average firing rate curve (in Hz)
%    curve.count         spike count curve
%    curve.time          occupancy curve (in s)
%
%    stats.x             location of the peak rate (in bins)
%    stats.peak          in-field peak rate
%    stats.mean          in-field mean value
%    stats.size          field size (in bins)
%    stats.field         field (1 = bin in field, 0 = bin not in field)
%    stats.fieldX        field x boundaries (in bins)
%    stats.specificity   spatial specificity (Skaggs et al., 1993)
%
%    For circular data:
%
%    stats.m             mean angle
%    stats.mode          distribution mode, prefered angle
%    stats.r             mean resultant length
%    stats.k             von Mises concentration
%
%  NOTE
%
%    This function is provided for convenience. It simply calls <a href="matlab:help Map">Map</a> and <a href="matlab:help MapStats">MapStats</a>
%    using the same parameters. The outputs are the same except for curve.z which
%    is replaced by curve.rate.
%
%  SEE
%
%    See also Map, MapStats, FiringMap, PlotXY.

% Copyright (C) 2005-2016 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Default values
type = 'linear';

% Check number of parameters
if nargin < 2,
  error('Incorrect number of parameters (type ''help <a href="matlab:help FiringCurve">FiringCurve</a>'' for details).');
end

if size(samples,2) ~= 2,
  error('Parameter ''samples'' is not a Nx2 matrix (type ''help <a href="matlab:help FiringCurve">FiringCurve</a>'' for details).');
end

im = 1;argsm = {};
is = 1;argss = {};
% Parse parameter list
for i = 1:2:length(varargin),
	if ~ischar(varargin{i}),
		error(['Parameter ' num2str(i+2) ' is not a property (type ''help <a href="matlab:help FiringCurve">FiringCurve</a>'' for details).']);
	end
	switch(lower(varargin{i})),
		case 'type',
			argss{is} = 'type';
			argss{is+1} = varargin{i+1}(1);
			is = is+2;
			argsm{im} = 'type';
			argsm{im+1} = [varargin{i+1}(1) 'l'];
			im = im+2;
		case {'threshold','minsize','minpeak','verbose','debug','localmax'},
			argss{is} = varargin{i};
			argss{is+1} = varargin{i+1};
			is = is+2;
		otherwise,
			argsm{im} = varargin{i};
			argsm{im+1} = varargin{i+1};
			im = im+2;
  end
end

curve = Map(samples,spikes,argsm{:});
if nargout == 2,
	stats = MapStats(curve,argss{:});
end
curve.rate = curve.z;
curve = rmfield(curve,'z');
