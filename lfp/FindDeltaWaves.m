function delta = FindDeltaWaves(filtered,varargin)

%FindDeltaWaves - Find cortical delta waves (1-6Hz waves).
%
%  USAGE
%
%    delta = FindDeltaWaves(filtered,<options>)
%
%    filtered       delta-band filtered LFP <a href="matlab:help samples">samples</a> (one channel). This must
%                   be restricted to slow wave sleep periods for the algorithm
%                   to perform best.
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'thresholds'  thresholds for z-scored minimum peak and trough amplitudes
%                   (default = [1 2 0 1.5], see NOTE below)
%     'durations'   min and max wave durations in ms (default = [150 500])
%    =========================================================================
%
%  OUTPUT
%
%    delta          for each delta wave, times and z-scored amplitudes of
%                   the beginning, peak and trough of the wave, in a Nx6
%                   matrix [start_t peak_t end_t start_z peak_z end_z]
%
%  NOTE
%
%    To be selected, candidate delta waves must fulfill one of two amplitude
%    conditions. The peak and trough must both exceed a threshold, but one can
%    compensate for the other, i.e. a large peak requires a smaller trough and
%    vice versa.
%
%    More precisely, there are two thresholds for the peak, p and P (p<P), and
%    two for the trough, t and T (t<T). The amplitude must fulfill one of the
%    following conditions:
%
%           peak > p and trough > T
%      or   peak > P and trough > t
%
%    Thresholds are given in SDs, and the above conditions relate to absolute
%    values.
%
%  SEE
%
%    See also FilterLFP, FindRipples.

% Copyright (C) 2012-2017 Michaël Zugaro, 2012-2015 Nicolas Maingret, 
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Default values
highPeak = 2 ; % Threshold for filtered signal (number of SDs)
lowPeak = 1;
highTrough = 1.5;
lowTrough = 0;
minDuration = 150; % min time between successive zero crossings (in ms)
maxDuration = 450; % max time between successive zero crossings (in ms)

% Check number of parameters
if nargin < 1 | mod(length(varargin),2) ~= 0,
  error('Incorrect number of parameters (type ''help <a href="matlab:help FindDeltaWaves">FindDeltaWaves</a>'' for details).');
end

% Check parameter sizes
if ~isdmatrix(filtered,'@2'),
	error('Parameter ''filtered'' is not a Nx2 matrix (type ''help <a href="matlab:help FindDeltaWaves">FindDeltaWaves</a>'' for details).');
end

% Parse parameter list
for i = 1:2:length(varargin),
	if ~ischar(varargin{i}),
		error(['Parameter ' num2str(i+2) ' is not a property (type ''help <a href="matlab:help FindDeltaWaves">FindDeltaWaves</a>'' for details).']);
	end
	switch(lower(varargin{i})),
		case {'thresholds','amplitudes'},
			thresholds = varargin{i+1};
			if ~isdvector(thresholds,'#4'),
				error('Incorrect value for property ''thresholds'' (type ''help <a href="matlab:help FindDeltaWaves">FindDeltaWaves</a>'' for details).');
			end
			lowPeak = thresholds(1);
			highPeak = thresholds(2);
			lowTrough = thresholds(3);
			highTrough = thresholds(4);
			if lowPeak > highPeak || lowTrough > highTrough,
				error('Inconsistent amplitude thresholds (type ''help <a href="matlab:help FindDeltaWaves">FindDeltaWaves</a>'' for details).');
			end
		case 'durations',
			durations = varargin{i+1};
			if ~isdvector(durations,'#2','<','>0'),
				error('Incorrect value for property ''durations'' (type ''help <a href="matlab:help FindDeltaWaves">FindDeltaWaves</a>'' for details).');
			end
			if durations(2) < 1,
				warning('Delta wave min and max durations are less than 1 ms, assuming seconds.');
				durations = durations * 1000;
			end
			minDuration = durations(1); maxDuration = durations(2);
		otherwise,
			error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help FindDeltaWaves">FindDeltaWaves</a>'' for details).']);
	end
end

% Find local minima and maxima corresponding to beginning, peak and end of delta waves
% This is done by finding zero crossings of the (z-scored) derivative of the signal

% Differentiate, filter and z-score signal
z = filtered;
z(:,2) = [diff(filtered(:,2));0];
z = bz_Filter(z,'passband',[0 6],'order',8,'FMAlegacy',true);
z(:,2) = zscore(z(:,2));

% Find positions (in # samples) of zero crossings
[up,down] = ZeroCrossings(z);
down = find(down);
up = find(up);
if down(1) < up(1), down(1) = []; end

% List positions (in # samples) of successive up,down,up crossings in an Nx3 matrix
n = length(up);
where = [up(1:n-1) down(1:n-1) up(2:n)];

% List positions but also z-scored amplitudes in an Nx6 matrix (positions then amplitudes)
z = filtered;
z(:,2) = zscore(filtered(:,2));
delta = z(where,:);
delta = reshape(delta,size(where,1),6);

% Discard waves that are too long or too short
duration = delta(:,3) - delta(:,1);
delta(duration<minDuration/1000|duration>maxDuration/1000,:) = [];

% Threshold z-scored peak and trough amplitudes
peak = delta(:,5);
trough = delta(:,6);
case1 = peak > highPeak & trough <= -lowTrough;
case2 = peak >= lowPeak & trough < -highTrough;
delta = delta(case1|case2,:);
