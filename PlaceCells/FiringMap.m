function [map,stats] = FiringMap(positions,spikes,varargin)

%FiringMap - Compute firing map (e.g. for a place cell, head direction cell, etc.)
%
%  Compute firing map (e.g. for a place cell), as well as occupancy and spike
%  count maps. Optionally, field statistics can also be computed, including
%  in-field peak and mean firing rates, firing field size, etc.
%
%  USAGE
%
%    [map,stats] = FiringMap(positions,spikes,<options>)
%
%    positions      position <a href="matlab:help samples">samples</a> in 
%                   1D or 2 dimensional data.
%    spikes         spike timestamps
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'smooth'      smoothing size in bins (0 = no smoothing, default = 2)
%     'nBins'       number of horizontal and vertical bins (default = [50 50])
%     'minTime'     minimum time spent in each bin (in s, default = 0)
%     'mode'        'interpolate' to interpolate missing points (< minTime),
%                   or 'discard' to discard them (default)
%     'maxDistance' maximal distance for interpolation (default = 5)
%     'maxGap'      z values recorded during time gaps between successive (x,y)
%                   samples exceeding this threshold (e.g. undetects) will not
%                   be interpolated; also, such long gaps in (x,y) sampling
%                   will be clipped to 'maxGap' to compute the occupancy map
%                   (default = 0.100 s)
%     'type'        two letters (one for X and one for Y) indicating which
%                   coordinates are linear ('l') and which are circular ('c')
%                   (default 'll')
%     'threshold'   values above threshold*peak belong to the field
%                   (default = 0.2)
%     'minSize'     fields smaller than this size are considered spurious
%                   and ignored (default = 100)
%     'minPeak'     peaks smaller than this size are considered spurious
%                   and ignored (default = 1)
%     'nShuffles'   number of shuffles to compute the p-value for specificity
%                   (default = 0, no p-value computed)
%     'verbose'     display processing information (default = 'off')
%    =========================================================================
%
%  OUTPUT
%
%    map.x               x bins
%    map.y               y bins
%    map.rate            average firing rate map (in Hz)
%    map.count           spike count map
%    map.time            occupancy map (in s)
%
%    stats.x             abscissa of the peak rate (in bins)
%    stats.y             ordinate of the peak rate (in bins)
%    stats.peak          in-field peak rate
%    stats.mean          in-field mean value
%    stats.size          field size (in bins)
%    stats.field         field (1 = bin in field, 0 = bin not in field)
%    stats.fieldX        field x boundaries (in bins)
%    stats.fieldY        field y boundaries (in bins)
%    stats.specificity   spatial specificity (Skaggs et al., 1993)
%    stats.p             spatial specificity p-value (computed over shuffles)
%
%  NOTE
%
%    This function is provided for convenience. It simply calls <a href="matlab:help Map">Map</a> and <a href="matlab:help MapStats">MapStats</a>
%    using the same parameters. The outputs are the same except for map.z which
%    is replaced by map.rate.
%
%  SEE
%
%    See also Map, MapStats, PlotColorMap.

% Copyright (C) 2023 Ralitsa Todorova, 2004-2013 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Check number of parameters
if nargin < 2,
  error('Incorrect number of parameters (type ''help <a href="matlab:help FiringMap">FiringMap</a>'' for details).');
end

nDim = size(positions,2)-1;
if ~ismember(nDim,[1 2])
  error('Parameter ''positions'' is not a Nx2 or Nx3 matrix (type ''help <a href="matlab:help FiringMap">FiringMap</a>'' for details).');
end

nShuffles = 0;
im = 1;argsm = {};
is = 1;argss = {};
% Parse parameter list
for i = 1:2:length(varargin),
    if ~ischar(varargin{i}),
        error(['Parameter ' num2str(i+2) ' is not a property (type ''help <a href="matlab:help FiringMap">FiringMap</a>'' for details).']);
    end
    switch(lower(varargin{i})),
        case {'nshuffles'}, % argument for this function only
            nShuffles = varargin{i+1};
        case {'mintime'}, % arguments to be given to both Map and MapStats
            argss{is} = varargin{i};
            argss{is+1} = varargin{i+1};
            is = is+2;
            argsm{im} = varargin{i};
            argsm{im+1} = varargin{i+1};
            im = im+2;
        case 'type', % special case to be given to both Map and MapStats
            argss{is} = 'type';
            argss{is+1} = varargin{i+1}(1:nDim);
            is = is+2;
            argsm{im} = 'type';
            argsm{im+1} = [varargin{i+1}(1:nDim) 'l'];
            im = im+2;
        case {'threshold','minsize','minpeak','verbose','debug','localmax'}, % arguments to be given to MapStats only
            argss{is} = varargin{i};
            argss{is+1} = varargin{i+1};
            is = is+2;
        otherwise, % arguments to be given to Map only
            argsm{im} = varargin{i};
            argsm{im+1} = varargin{i+1};
            im = im+2;
    end
end

map = Map(positions,spikes,argsm{:});
if nargout == 2,
    stats = MapStats(map,argss{:});
    if nShuffles>0
        maxTime = max(positions(:,1)); specificity = nan(nShuffles,1);
        for i=1:nShuffles
            shifted = sortrows([rem(positions(:,1)+rand(1)*maxTime,maxTime) positions(:,2:end)]);
            shuffled = Map(shifted,spikes,argsm{:});
            s = MapStats(shuffled,argss{:});
            specificity(i) = s.specificity;
        end
        stats.p = 1-sum(specificity<stats.specificity)./sum(~isnan(specificity));
    else
        stats.p = nan;
    end
end
map.rate = map.z;
map = rmfield(map,'z');
