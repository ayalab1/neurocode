function [spectra,f] = PointSpectra(timesCell,varargin)

%PointSpectrogra - Compute point process spectrogram of multiple units using multi-taper estimation
%
% This is an optimized and simplified implementation inspired by the Chronux toolbox.
%
% USAGE
%
%  [spectra,f] = PointSpectra(times,<options>)
%
%  INPUT
%
%  times        a cell array containing the event times of each event type, e.g. spike times of each cell
%  <options>    optional list of property-value pairs (see table below)
%
%   == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == = 
%   Properties  Values
%  -------------------------------------------------------------------------
%   'frequency' sampling rate (in Hz) (default = 20000Hz)
%   'range'     frequency range (in Hz) (default = all)
%   'window'    duration (in s) of the time window (default = 5)
%   'overlap'   overlap between successive windows (default = window/2)
%   'step'      step between successive windows (default = window/2)
%   'tapers'    relative resolution and order of the tapers [NW K] (default = [3 5])
%   'pad'       FFT padding (see help for <a href = "matlab:help mtspecgrampt">mtspecgrampt</a>) (default = 0)
%   'show'      plot spectrogram (default = 'off')
%   'cutoffs'   cutoff values for color plot (default = [0 13])
%   == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == = 
%
% NOTES
%
%  The time displacement between successive short time spectra can be supplied
%  either as a 'step' (explicit time difference) or as an 'overlap' (between
%  successive time windows)
%
% OUTPUT
%
%  spectra     unit-frequency matrix (each row is the spectrum of an individual unit)
%  f           frequency bins (correspond to the columns of "spectra")
%
% SEE
%
%  See also PointSpectrogram, RelativePointSpectra
%
% Copyright (C) 2020-2023 by Ralitsa Todorova, 2002-2012 by Michaël Zugaro, and chronux
% 
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Defaults
Fs = 1250;
window = 5;
range = [1 20];
overlap = [];
step = [];
show = 'off';
cutoffs = [0 13];
tapers0 = [3 5];
pad = 0;

% Check number of parameters
if nargin < 1 || mod(length(varargin),2) ~= 0,
 error('Incorrect number of parameters (type ''help <a href = "matlab:help MTPointSpectrogram">MTPointSpectrogram</a>'' for details).');
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
			tapers0 = varargin{i+1};
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
times = cat(1,timesCell{:});
timesRange = [min(times) max(times)];
t = mean(timesRange);
window = floor(diff(timesRange));
nSamplesPerWindow = round(Fs*window); % number of samples in window
nfft = max(2^(nextpow2(nSamplesPerWindow)+pad),nSamplesPerWindow);
fAll = linspace(0,Fs,nfft)'; 
ok = fAll >= range(1) & fAll<= range(end);
Nf = sum(ok);
tapers = dpss(nSamplesPerWindow,tapers0(1),tapers0(2))*sqrt(Fs);
spectra = zeros(Nf,length(timesCell));
H = fft(tapers,nfft,1); % fft of tapers
% restrict fft of tapers to required frequencies
f = fAll(ok); H = H(ok,:);
w = 2*pi*f; % angular frequencies at which ft is to be evaluated
timegrid = linspace(timesRange(1),timesRange(2),nSamplesPerWindow);
for i = 1:length(timesCell),
    times = timesCell{i};
    data = Restrict(times,timegrid([1 end]));
    if isempty(data), continue; end
    data_proj = interp1(timegrid',tapers,data);
    exponential = exp(-1i*w*(data-timegrid(1))');
    J = exponential*data_proj-H*length(data)/length(timegrid);
    spectra(:,i) = squeeze(mean(conj(J).*J,2));
end
if strcmpi(show,'on'),
    logTransformed = log(spectrogram);
    PlotColorMap(logTransformed,1,'x',t,'y',f,'cutoffs',cutoffs);	
    xlabel('Time (s)');	ylabel('Frequency (Hz)');title('Power Spectrogram');
end