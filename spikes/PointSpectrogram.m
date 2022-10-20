function [spectrogram,t,f] = PointSpectrogram(times,varargin)

%PointSpectrogram - Compute point process spectrogram by multi-taper estimation
%
% [This is an optimized and simplified implementation inspired by the Chronux toolbox.]
%
% USAGE
%
%  [[spectrogram,t,f] = MTPointSpectrogram(times,<options>)]
%
%  INPUT
%
%  [times]     [event times, e.g. spike times]
%  <options>   optional list of property-value pairs (see table below)
%
%   == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == = 
%   Properties  Values
%  -------------------------------------------------------------------------
%   ['frequency'] [sampling rate (in Hz) (default = 20000Hz)]
%   ['range']     [frequency range (in Hz) (default = all)]
%   ['window']    [duration (in s) of the time window (default = 5)]
%   ['overlap']   [overlap between successive windows (default = window/2)]
%   ['step']      [step between successive windows (default = window/2)]
%   ['tapers']    [relative resolution and order of the tapers [NW K]
%                 (default = [3 5])]
%   ['pad']       [FFT padding (see help for <a href = "matlab:help
%                 mtspecgrampt">mtspecgrampt</a>) (default = 0)]
%   ['show']      [plot spectrogram (default = 'off')]
%   ['cutoffs']   [cutoff values for color plot (default = [0 13])]
%   == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == = 
%
% NOTES
%
%  [The time displacement between successive short time spectra can be supplied
%  either as a 'step' (explicit time difference) or as an 'overlap' (between
%  successive time windows)]
%
% OUTPUT
%
%  [spectrogram] [time-frequency matrix]
%  [t]           [time bins]
%  [f]           [frequency bins]
%
% SEE
%
%  See also SpectrogramBands, PlotColorMap.
%
% Copyright (C) [2020-2022] by Ralitsa Todorova, [2002-2012] by Michaël Zugaro, and chronux
% 
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Defaults
Fs = 20000;
window = 5;
range = [1 20];
overlap = [];
step = [];
show = 'off';
cutoffs = [0 13];
tapers = [3 5];
pad = 0;

% Check number of parameters
if nargin < 1 | mod(length(varargin),2) ~= 0,
 error('Incorrect number of parameters (type ''help <a href = "matlab:help MTPointSpectrogram">MTPointSpectrogram</a>'' for details).');
end

% Check parameter sizes
if ~isdvector(times),
	error('Parameter ''times'' is not a vector (type ''help <a href = "matlab:help MTPointSpectrogram">MTPointSpectrogram</a>'' for details).');
end

% Parse parameter list
for i = 1:2:length(varargin),
	if ~ischar(varargin{i}),
		error(['Parameter ' num2str(i+2) ' is not a property (type ''help <a href = "matlab:help MTPointSpectrogram">MTPointSpectrogram</a>'' for details).']);
	end
	switch(lower(varargin{i})),
		case 'frequency',
			Fs = varargin{i+1};
			if ~isdscalar(Fs,'>0'),
				error('Incorrect value for property ''frequency'' (type ''help <a href = "matlab:help MTPointSpectrogram">MTPointSpectrogram</a>'' for details).');
			end
		case 'range',
			range = varargin{i+1};
			if ~isdvector(range,'#2','>= 0','<'),
				error('Incorrect value for property ''range'' (type ''help <a href = "matlab:help MTPointSpectrogram">MTPointSpectrogram</a>'' for details).');
			end
		case 'window',
			window = varargin{i+1};
			if ~isdscalar(window,'>0'),
				error('Incorrect value for property ''window'' (type ''help <a href = "matlab:help MTPointSpectrogram">MTPointSpectrogram</a>'' for details).');
			end
		case 'overlap',
			overlap = varargin{i+1};
			if ~isdscalar(overlap,'>0'),
				error('Incorrect value for property ''overlap'' (type ''help <a href = "matlab:help MTPointSpectrogram">MTPointSpectrogram</a>'' for details).');
			end
		case 'step',
			step = varargin{i+1};
			if ~isdscalar(step,'>0'),
				error('Incorrect value for property ''step'' (type ''help <a href = "matlab:help MTPointSpectrogram">MTPointSpectrogram</a>'' for details).');
			end
		case 'tapers',
			tapers = varargin{i+1};
			if ~isivector(tapers,'#2','>0'),
				error('Incorrect value for property ''tapers'' (type ''help <a href = "matlab:help MTPointSpectrogram">MTPointSpectrogram</a>'' for details).');
			end
		case 'pad',
			pad = varargin{i+1};
			if ~isdscalar(pad,'>-1'),
				error('Incorrect value for property ''pad'' (type ''help <a href = "matlab:help MTPointSpectrogram">MTPointSpectrogram</a>'' for details).');
			end
		case 'show',
			show = varargin{i+1};
			if ~isstring(show,'on','off'),
				error('Incorrect value for property ''show'' (type ''help <a href = "matlab:help FindRipples">FindRipples</a>'' for details).');
			end
		case 'cutoffs',
			cutoffs = varargin{i+1};
			if ~isdvector(cutoffs,'#2','>= 0','<'),
				error('Incorrect value for property ''cutoffs'' (type ''help <a href = "matlab:help MTPointSpectrogram">MTPointSpectrogram</a>'' for details).');
      end
		otherwise,
			error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href = "matlab:help MTPointSpectrogram">MTPointSpectrogram</a>'' for details).']);
	end
end

% Determine step/overlap
if isempty(step),
	if isempty(overlap),
		step = window/2;
	end
else
	if ~isempty(overlap) && overlap ~= window-step,
		error('Incompatible ''step'' and ''overlap'' parameters (type ''help <a href = "matlab:help MTPointSpectrogram">MTPointSpectrogram</a>'' for details).');
	end
end

% Compute and plot spectrogram
t = (min(times)+window/2:step:max(times)-window/2)';
nSamplesPerWindow = round(Fs*window); % number of samples in window
nfft = max(2^(nextpow2(nSamplesPerWindow)+pad),nSamplesPerWindow);
fAll = linspace(0,Fs,nfft)'; 
ok = fAll >= range(1) & fAll<= range(end);
Nf = sum(ok);
tapers = dpss(nSamplesPerWindow,tapers(1),tapers(2))*sqrt(Fs);
nWindows = length(t);
spectrogram = zeros(Nf,nWindows);
H = fft(tapers,nfft,1); % fft of tapers
% restrict fft of tapers to required frequencies
f = fAll(ok); H = H(ok,:);
w = 2*pi*f; % angular frequencies at which ft is to be evaluated
for i = 1:nWindows,
    timegrid = linspace(t(i)-window/2,t(i)+window/2,nSamplesPerWindow);
    data = Restrict(times,timegrid([1 end]));
    if isempty(data), continue; end
    data_proj = interp1(timegrid',tapers,data);
    exponential = exp(-1i*w*(data-timegrid(1))');
    J = exponential*data_proj-H*length(data)/length(timegrid);
    spectrogram(:,i) = squeeze(mean(conj(J).*J,2));
end
if strcmpi(show,'on'),
    logTransformed = log(spectrogram);
    PlotColorMap(logTransformed,1,'x',t,'y',f,'cutoffs',cutoffs);	
    xlabel('Time (s)');	ylabel('Frequency (Hz)');title('Power Spectrogram');
end