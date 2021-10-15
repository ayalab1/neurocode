function spindles = FindSpindles(filtered,varargin)

%FindSpindles - Find thalamo-cortical spindles (9-17Hz oscillations).
%
%  USAGE
%
%    spindles = FindSpindles(filtered,<options>)
%
%    filtered       spindle-band filtered LFP <a href="matlab:help samples">samples</a> (one channel). This must
%                   be restricted to slow wave sleep periods for the algorithm
%                   to perform best.
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'threshold'   minimum z-scored envelope tail (default = 2.5)
%     'peak'        minimum z-scored envelope peak (default = 5)
%     'durations'   min and max spindle durations in ms (default = [400 3500])
%    =========================================================================
%
%  OUTPUT
%
%    spindles       for each spindle, start, peak and stop times, and z-scored
%                   peak amplitude, in a Nx4 matrix [start peak_t end peak_z]
%

% Copyright (C) 2012-2016 Michaël Zugaro, 2012-2015 Nicolas Maingret, 2016 Ralitsa Todorova
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Defaults
minDistance = 0.5;
threshold = 2.5;
minPeak = 5 ;
minDuration = 0.4;
maxDuration = 3.5;

% Constants
minBoutDuration = 0.050;
window = 0.320;

% Check number of parameters
if nargin < 1 | mod(length(varargin),2) ~= 0,
  error('Incorrect number of parameters (type ''help <a href="matlab:help FindSpindles">FindSpindles</a>'' for details).');
end

% Check parameter sizes
if ~isdmatrix(filtered,'@2'),
	error('Parameter ''filtered'' is not a Nx2 matrix (type ''help <a href="matlab:help FindSpindles">FindSpindles</a>'' for details).');
end

% Parse parameter list
for i = 1:2:length(varargin),
	if ~ischar(varargin{i}),
		error(['Parameter ' num2str(i+2) ' is not a property (type ''help <a href="matlab:help FindSpindles">FindSpindles</a>'' for details).']);
	end
	switch(lower(varargin{i})),
		case 'threshold',
			threshold = varargin{i+1};
			if ~isdscalar(threshold,'>0'),
				error('Incorrect value for property ''threshold'' (type ''help <a href="matlab:help FindSpindles">FindSpindles</a>'' for details).');
			end
		case 'peak',
			minPeak = varargin{i+1};
			if ~isdscalar(minPeak,'>0'),
				error('Incorrect value for property ''peak'' (type ''help <a href="matlab:help FindSpindles">FindSpindles</a>'' for details).');
			end
		case 'durations',
			durations = varargin{i+1};
			if ~isdvector(durations,'#2','>0'),
				error('Incorrect value for property ''durations'' (type ''help <a href="matlab:help FindSpindles">FindSpindles</a>'' for details).');
			end
			minDuration = durations(1)/1000;
			maxDuration = durations(2)/1000;
		otherwise,
			error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help FindSpindles">FindSpindles</a>'' for details).']);
	end
end


% Compute low-pass filtered squared envelope of z-scored signal
time = filtered(:,1);
frequency = 1/median(diff(time));
windowLength = round(frequency*window);
window = ones(windowLength,1)/windowLength;
envelope = abs(hilbert(zscore(filtered(:,2)))).^2;
envelope = filtfilt(window,1,envelope);

% Find start/stop as indices of threshold crossings, and discard incomplete pairs
crossings = envelope > threshold;
crossings = diff(crossings);
start = find(crossings>0);
stop = find(crossings<0);
if time(stop(1)) < time(start(1)), stop(1) = []; end
if time(start(end)) > time(stop(end)), start(end) = []; end

% Discard short events
tooShort = time(stop) - time(start) < minBoutDuration;
start(tooShort) = [];
stop(tooShort) = [];

% Merge events when they are too close
while true,
	tooClose = find(time(start(2:end)) - time(stop(1:end-1)) < minDistance);
	if isempty(tooClose), break; end
	start(tooClose+1) = [];
	stop(tooClose) = [];
end

% Discard events when envelope peak is too small
% For each timestamp of the LFP, get the id of the spindle
[~,id] = InIntervals(time,time([start stop]));
% Find peak amplitudes and indices
% (add 1 so that Accumulate does not complain when id=0, i.e. time points that do not belong to any candidate interval)
[peak_z,peak_i] = Accumulate(id+1,envelope,'mode','max');
% (now get rid of id=0, i.e. time points that do not belong to any candidate interval)
peak_z(1) = []; peak_i(1) = []; 
% Keep only spindles with large enough amplitude
tooSmall = peak_z < minPeak;
start(tooSmall) = [];
stop(tooSmall) = [];
peak_z(tooSmall) = [];
peak_i(tooSmall) = [];

% Replace indices with timestamps
start = time(start);
stop = time(stop);
peak_t = time(peak_i);

% Spindle start may be inaccurate due to leading delta wave, we need to correct for this
% 1) Update threshold to 1/3 peak (if this is higher than previous threshold)
newThrehold = peak_z/3;
% 2) Select spindles that need correction
indices = find(newThrehold>threshold); 
% 3) Find last threshold crossing before peak (vectorized code)
% Candidate intervals for the new spindle start = between the previous spindle end and current spindle peak
if indices(1) > 1, candidateIntervals = [stop(indices-1) peak_t(indices)];
else candidateIntervals = [0 peak_t(indices(1)); stop(indices(2:end)-1) peak_t(indices(2:end))]; end
% For each timestamp, get the id of the interval in which it falls
[in,id] = InIntervals(time,candidateIntervals);
belowThreshold = false(size(id));
belowThreshold(in) = envelope(in) < newThrehold(indices(id(in)));
dt = mode(diff(time));
start(indices) = Accumulate(id(in&belowThreshold),time(in&belowThreshold),'mode','max')+dt;


% Merge events when they are too close
n = length(peak_z);
while true,
	tooClose = find(start(2:end) - stop(1:end-1) < minDistance);
	if isempty(tooClose), break; end
	start(tooClose+1) = [];
	stop(tooClose) = [];
	[peak_z,i] = max([peak_z(tooClose) peak_z(tooClose+1)],[],2);
	peak_t = peak_t(tooClose+i-1);
end

% Discard events that are too long or too short
duration = stop-start;
bad = abs(duration)<minDuration | abs(duration)>maxDuration;
start(bad) = [];
stop(bad) = [];
peak_z(bad) = [];
peak_t(bad) = [];

spindles = [start peak_t stop peak_z];
spindles = unique(spindles,'rows');
