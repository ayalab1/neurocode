function [coherence,phase,f,sc,sp] = MTCoherence(lfp1,lfp2,varargin)

%MTCoherence - Compute LFP coherence by multi-taper estimation.
%
%  The coherence is computed as the average coherogram, using the <a href="http://www.chronux.org">chronux</a> toolbox.
%
%  USAGE
%
%    [coherence,phase,f,sc,sp] = MTCoherence(lfp1,lfp2,<options>)
%
%    lfp1,lfp2      unfiltered LFP <a href="matlab:help samples">samples</a> (one channel each).
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'frequency'   sampling rate (in Hz) (default = from timestamps if
%                   available, otherwise 1250Hz)
%     'range'       frequency range (in Hz) (default = all)
%     'tapers'      relative resolution and order of the tapers [NW K]
%                   (default = [3 5])
%     'pad'         FFT padding (see help for <a href="matlab:help coherencyc">coherencyc</a>) (default = 0)
%     'show'        plot results (default = 'off')
%    =========================================================================
%
%  NOTES
%
%    The LFPs can be provided either as time stamped matrices (lists of time-voltage
%    pairs), or as a voltage vectors - in which case the frequency must be specified.
%
%  OUTPUT
%
%    coherence      coherence magnitude
%    phase          coherence phase
%    f              frequency bins
%    sc             error on magnitude
%    sp             error on phase
%
%  DEPENDENCIES
%
%    This function requires the <a href="http://www.chronux.org">chronux</a> toolbox.
%
%  SEE
%
%    See also MTCoherogram, MTSpectrum, MTSpectrogram, PlotMean.

% Copyright (C) 2010-2014 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Make sure chronux is installed and functional
CheckChronux('coherencyc');

% Defaults
f = 1250;
frequency = [];
range = [];
show = 'off';
tapers = [3 5];
pad = 0;

% Check number of parameters
if nargin < 2 | mod(length(varargin),2) ~= 0,
  error('Incorrect number of parameters (type ''help <a href="matlab:help MTCoherence">MTCoherence</a>'' for details).');
end

% Check parameter sizes
if size(lfp1,2) ~= 1 && size(lfp1,2) ~= 2,
	error('Parameter ''lfp1'' is not a vector or a Nx2 matrix (type ''help <a href="matlab:help MTCoherence">MTCoherence</a>'' for details).');
end
if size(lfp2,2) ~= 1 && size(lfp2,2) ~= 2,
	error('Parameter ''lfp2'' is not a vector or a Nx2 matrix (type ''help <a href="matlab:help MTCoherence">MTCoherence</a>'' for details).');
end

% Parse parameter list
v = {};
for i = 1:2:length(varargin),
	if ~ischar(varargin{i}),
		error(['Parameter ' num2str(i+2) ' is not a property (type ''help <a href="matlab:help MTCoherence">MTCoherence</a>'' for details).']);
	end
	switch(lower(varargin{i})),
		case 'frequency',
			frequency = varargin{i+1};
			if ~isdscalar(frequency,'>0'),
				error('Incorrect value for property ''frequency'' (type ''help <a href="matlab:help MTCoherence">MTCoherence</a>'' for details).');
			end
		case 'range',
			range = varargin{i+1};
			if ~isdvector(range,'#2','<','>=0'),
				error('Incorrect value for property ''range'' (type ''help <a href="matlab:help MTCoherence">MTCoherence</a>'' for details).');
			end
		case 'tapers',
			tapers = varargin{i+1};
			if ~isivector(tapers,'#2','>0'),
				error('Incorrect value for property ''tapers'' (type ''help <a href="matlab:help MTCoherence">MTCoherence</a>'' for details).');
			end
		case 'pad',
			pad = varargin{i+1};
			if ~isiscalar(pad,'>-1'),
				error('Incorrect value for property ''pad'' (type ''help <a href="matlab:help MTCoherence">MTCoherence</a>'' for details).');
			end
		case 'show',
			show = varargin{i+1};
			if ~isastring(show,'on','off'),
				error('Incorrect value for property ''show'' (type ''help <a href="matlab:help MTCoherence">MTCoherence</a>'' for details).');
			end
		otherwise,
			error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help MTCoherence">MTCoherence</a>'' for details).']);
	end
	if ~strcmp(varargin{i},'show'), v = {v{:},varargin{i:i+1}}; end
end

% Compute spectrogram
[coherogram,phase,~,f] = MTCoherogram(lfp1,lfp2,v{:});

% Compute coherence
coherogram = coherogram';
coherence = mean(coherogram);
sc = sqrt(var(coherogram));
phase = phase';
[phase,sp] = CircularConfidenceIntervals(phase);

% Plot coherence
if strcmp(lower(show),'on'),
	figure;
	subplot(2,1,1);
	PlotMean(f,coherence,coherence-sc,coherence+sc,':');
	xlabel('Frequency (Hz)');
	ylabel('Coherence');
	title('Coherence Amplitude');
	subplot(2,1,2);
	PlotMean(f,phase*180/pi,sp(1,:)*180/pi,sp(2,:)*180/pi,':');
	xlabel('Frequency (Hz)');
	ylabel('Phase (°)');
	title('Coherence Phase');
end
