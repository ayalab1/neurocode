function [spectrogram,t,f] = PointSpectrogram(times,varargin)

%PointSpectrogram - Compute point process spectrogram by multi-taper estimation.
%
% The spectrogram is computed using the <a href = "http://www.chronux.org">chronux</a> toolbox.
%
% USAGE
%
%  [spectrogram,t,f] = MTPointSpectrogram(times,<options>)
%
%  times     event times, e.g. spike times
%  <options>   optional list of property-value pairs (see table below)
%
%   == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == = 
%   Properties  Values
%  -------------------------------------------------------------------------
%   'frequency'  sampling rate (in Hz) (default = 20000Hz)
%   'range'    frequency range (in Hz) (default = all)
%   'window'   duration (in s) of the time window (default = 5)
%   'overlap'   overlap between successive windows (default = window/2)
%   'step'    step between successive windows (default = window/2)
%   'tapers'   relative resolution and order of the tapers [NW K]
%          (default = [3 5])
%   'pad'     FFT padding (see help for <a href = "matlab:help mtspecgrampt">mtspecgrampt</a>) (default = 0)
%   'show'    plot spectrogram (default = 'off')
%   'cutoffs'   cutoff values for color plot (default = [0 13])
%   == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == = 
%
% NOTES
%
%  The time displacement between successive short time spectra can be supplied
%  either as a 'step' (explicit time difference) or as an 'overlap' (between
%  successive time windows).
%
% OUTPUT
%
%  spectrogram  time-frequency matrix
%  t       time bins
%  f       frequency bins
%
% DEPENDENCIES
%
%  This function requires the <a href = "http://www.chronux.org">chronux</a> toolbox.
%
% SEE
%
%  See also SpectrogramBands, PlotColorMap.
% I put all the functions calculating a point spectrogram in this function 
% in order to avoid multiple calculations of the fft of tapers, thus dramatically
% increasing calculation speed. This was just for optimisation purposes, all 
% programs were developed as referenced (i.e. by Michaël Zugaro and Chronux).

% Copyright (C) 2004,2012 by Michaël Zugaro, chronux, and 2020 Ralitsa Todorova
% (optimization)
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Ralitsa Todorova put all the functions calculating a point spectrogram in this function 
% in order to avoid multiple calculations of the fft of tapers, thus dramatically
% increasing calculation speed. This was just for optimisation purposes, all 
% programs were developed as referenced (i.e. by Michaël Zugaro and Chronux).

% Defaults
Fs = 20000;
window = 5;
range = [];
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
f = getfgrid(Fs,nfft,range); Nf = length(f);
tapers = dpss(nSamplesPerWindow,tapers(1),tapers(2))*sqrt(Fs);
nWindows = length(t);
spectrogram = zeros(Nf,nWindows);
H = fft(tapers,nfft,1); % fft of tapers
f = (0:Fs/nfft:Fs)'; % all possible frequencies
% restrict fft of tapers to required frequencies
indexF = find(f>=range(1) & f<=range(end)); 
f = f(indexF); H = H(indexF,:);
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

%% [spectrogram,t,f] = mtspecgrampt(times,[window step],parameters); % Chronux function
%   function [S,t,f,R,Serr] = mtspecgrampt(data,movingwin,parameters,fscorr)
% Multi-taper time-frequency spectrum - point process times
%
% Usage:
%
% [S,t,f,R,Serr] = mtspecgrampt(data,movingwin,parameters,fscorr)
% Input:
%    data    (structure array of spike times with dimension channels/trials;
%          also accepts 1d array of spike times) -- required
%    movingwin     (in the form [window,winstep] i.e length of moving
%                         window and step size.
%
%    parameters: structure with fields tapers, pad, Fs, range, err, trialave
%    - optional
%      tapers : precalculated tapers from dpss or in the one of the following
%          forms:
%          (1) A numeric vector [TW K] where TW is the
%            time-bandwidth product and K is the number of
%            tapers to be used (less than or equal to
%            2TW-1).
%          (2) A numeric vector [W T p] where W is the
%            bandwidth, T is the duration of the data and p
%            is an integer such that 2TW-p tapers are used. In
%            this form there is no default i.e. to specify
%            the bandwidth, you have to specify T and p as
%            well. Note that the units of W and T have to be
%            consistent: if W is in Hz, T must be in seconds
%            and vice versa. Note that these units must also
%            be consistent with the units of parameters.Fs: W can
%            be in Hz if and only if parameters.Fs is in Hz.
%            The default is to use form 1 with TW = 3 and K = 5
%           Note that T has to be equal to movingwin(1).
%
%	    pad		  (padding factor for the FFT) - optional (can take values -1,0,1,2...).
%          -1 corresponds to no padding, 0 corresponds to padding
%          to the next highest power of 2 etc.
%			   	 e.g. For N = 500, if PAD = -1, we do not pad; if PAD = 0, we pad the FFT
%			   	 to 512 points, if pad = 1, we pad to 1024 points etc.
%			   	 Defaults to 0.
%      Fs  (sampling frequency) - optional. Default 1.
%      range  (frequency band to be used in the calculation in the form
%                  [fmin fmax])- optional.
%                  Default all frequencies between 0 and Fs/2
%      err (error calculation [1 p] - Theoretical error bars; [2 p] - Jackknife error bars
%                  [0 p] or 0 - no error bars) - optional. Default 0.
%      trialave (average over trials/channels when 1, don't average when 0) - optional. Default 0
%    fscorr  (finite size corrections, 0 (don't use finite size corrections) or
%        1 (use finite size corrections) - optional
%        (available only for spikes). Defaults 0.
%
% Output:
%    S    (spectrogram with dimensions time x frequency x channels/trials if trialave = 0;
%        dimensions time x frequency if trialave = 1)
%    t    (times)
%    f    (frequencies)
%
%    Serr  (error bars) - only if err(1)>= 1
%
%% function [S,f,R,Serr] = mtspectrumpt(data,parameters,fscorr,t)
% Multi-taper spectrum - point process times
%
% Usage:
%
% [S,f,R,Serr] = mtspectrumpt(data,parameters,fscorr,t)
% Input:
%    data    (structure array of spike times with dimension channels/trials;
%          also accepts 1d array of spike times) -- required
%    parameters: structure with fields tapers, pad, Fs, range, err, trialave
%    - optional
%      tapers : precalculated tapers from dpss or in the one of the following
%          forms:
%          (1) A numeric vector [TW K] where TW is the
%            time-bandwidth product and K is the number of
%            tapers to be used (less than or equal to
%            2TW-1).
%          (2) A numeric vector [W T p] where W is the
%            bandwidth, T is the duration of the data and p
%            is an integer such that 2TW-p tapers are used. In
%            this form there is no default i.e. to specify
%            the bandwidth, you have to specify T and p as
%            well. Note that the units of W and T have to be
%            consistent: if W is in Hz, T must be in seconds
%            and vice versa. Note that these units must also
%            be consistent with the units of parameters.Fs: W can
%            be in Hz if and only if parameters.Fs is in Hz.
%            The default is to use form 1 with TW = 3 and K = 5
%
%	    pad		  (padding factor for the FFT) - optional (can take values -1,0,1,2...).
%          -1 corresponds to no padding, 0 corresponds to padding
%          to the next highest power of 2 etc.
%			   	 e.g. For N = 500, if PAD = -1, we do not pad; if PAD = 0, we pad the FFT
%			   	 to 512 points, if pad = 1, we pad to 1024 points etc.
%			   	 Defaults to 0.
%      Fs  (sampling frequency) - optional. Default 1.
%      range  (frequency band to be used in the calculation in the form
%                  [fmin fmax])- optional.
%                  Default all frequencies between 0 and Fs/2
%      err (error calculation [1 p] - Theoretical error bars; [2 p] - Jackknife error bars
%                  [0 p] or 0 - no error bars) - optional. Default 0.
%      trialave (average over channels/trials when 1, don't average when 0) - optional. Default 0
%    fscorr  (finite size corrections, 0 (don't use finite size corrections) or
%        1 (use finite size corrections) - optional
%        (available only for spikes). Defaults 0.
%    t    (time grid over which the tapers are to be calculated:
%           this argument is useful when calling the spectrum
%           calculation routine from a moving window spectrogram
%           calculation routine). If left empty, the spike times
%           are used to define the grid.
% Output:
%    S    (spectrum with dimensions frequency x channels/trials if trialave = 0;
%        dimension frequency if trialave = 1)
%    f    (frequencies)
%    R    (rate)
%    Serr  (error bars) - only if err(1)>= 1
%
%% function [J,Msp,Nsp] = mtfftpt(data,tapers,nfft,t,f,findx)
% Multi-taper fourier transform for point process given as times
%
% Usage:
% [J,Msp,Nsp] = mtfftpt (data,tapers,nfft,t,f,findx) - all arguments required
% Input: 
%    data    (struct array of times with dimension channels/trials; 
%          also takes in 1d array of spike times as a column vector) 
%    tapers   (precalculated tapers from dpss) 
%    nfft    (length of padded data) 
%    t      (time points at which tapers are calculated)
%    f      (frequencies of evaluation)
%    findx    (index corresponding to frequencies f) 
% Output:
%    J (fft in form frequency index x taper index x channels/trials)
%    Msp (number of spikes per sample in each channel)
%    Nsp (number of spikes in each channel)

function [f,findx]=getfgrid(Fs,nfft,fpass)
% Helper function that gets the frequency grid associated with a given fft based computation
% Called by spectral estimation routines to generate the frequency axes 
% Usage: [f,findx]=getfgrid(Fs,nfft,fpass)
% Inputs:
% Fs        (sampling frequency associated with the data)-required
% nfft      (number of points in fft)-required
% fpass     (band of frequencies at which the fft is being calculated [fmin fmax] in Hz)-required
% Outputs:
% f         (frequencies)
% findx     (index of the frequencies in the full frequency grid). e.g.: If
% Fs=1000, and nfft=1048, an fft calculation generates 512 frequencies
% between 0 and 500 (i.e. Fs/2) Hz. Now if fpass=[0 100], findx will
% contain the indices in the frequency grid corresponding to frequencies <
% 100 Hz. In the case fpass=[0 500], findx=[1 512].
if nargin < 3; error('Need all arguments'); end;
df=Fs/nfft;
f=0:df:Fs; % all possible frequencies
f=f(1:nfft);
if length(fpass)~=1;
   findx=find(f>=fpass(1) & f<=fpass(end));
else
   [fmin,findx]=min(abs(f-fpass));
   clear fmin
end;
f=f(findx);
