function [spectrogram,t,f] = WaveletSpectrogram(lfp,varargin)

%WaveletSpectrogram - Compute LFP wavelet spectrogram
%
%  The spectrogram is computed using the <a href="http://paos.colorado.edu/research/wavelets/">wavelet toolbox by C. Torrence and G. Compo</a>.
%
%  USAGE
%
%    [spectrogram,t,f] = WaveletSpectrogram(lfp,<options>)
%
%    lfp            unfiltered LFP <a href="matlab:help samples">samples</a> (one channel).
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'range'       frequency range (in Hz) (default = all)
%     'resolution'  desired frequency resolution (sub-octaves per octave, default = 0.05)
%     'step'        desired step between successive windows (default = window/2):
%                   output will be downsampled to this value
%     'interpolate' interpolate to create equidistant frequency bins instead
%                   of bins in power of 2 (default = false)
%    =========================================================================
%
%  OUTPUT
%
%    spectrogram    time-frequency matrix
%    t              time bins
%    f              frequency bins
%
%  SEE
%    See "http://paos.colorado.edu/research/wavelets/"
%    See also MTSpectrogram, MTSpectrum, SpectrogramBands, MTCoherence, MTCoherogram, PlotColorMap
%
% Copyright (C) 2018-2022 by Ralitsa Todorova
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Check number of parameters
if nargin < 1 || mod(length(varargin),2) ~= 0,
    error(['Incorrect number of parameters (type ''help <a href="matlab:help WaveletSpectrogram">WaveletSpectrogram</a>'' for details).']);
end

% Check parameter sizes
if size(lfp,2) ~= 2,
    error(['Parameter ''lfp'' is not a Nx2 matrix (type ''help <a href="matlab:help WaveletSpectrogram">WaveletSpectrogram</a>'' for details).']);
end

% Defaults
range = [];
pad = 1;      % pad the time series with zeroes (recommended)
resolution = 0.05;    % this will do (1/resolution) sub-octaves per octave
step = [];
interpolate = false;

% Parse parameter list
for i = 1:2:length(varargin),
    if ~ischar(varargin{i}),
        error(['Parameter ' num2str(i+2) ' is not a property (type ''help <a href="matlab:help WaveletSpectrogram">WaveletSpectrogram</a>'' for details).']);
    end
    switch(lower(varargin{i})),
        case 'range',
            range = varargin{i+1};
            if ~isdvector(range,'#2','>=0','<'),
                error(['Incorrect value for property ''' varargin{i} ''' (type ''help <a href="matlab:help WaveletSpectrogram">WaveletSpectrogram</a>'' for details).']);
            end
        case 'step',
            step = varargin{i+1};
            if ~isdscalar(step,'>0'),
                error(['Incorrect value for property ''' varargin{i} ''' (type ''help <a href="matlab:help MTSpectrogram">MTSpectrogram</a>'' for details).']);
            end
        case 'resolution',
            resolution = varargin{i+1};
            if ~isdvector(resolution,'#1','>0'),
                error(['Incorrect value for property ''' varargin{i} ''' (type ''help <a href="matlab:help WaveletSpectrogram">WaveletSpectrogram</a>'' for details).']);
            end
        case 'interpolate',
            interpolate = varargin{i+1};
            if isastring(interpolate,'on','off'), if isastring(interpolate,'on'), interpolate = true; else interpolate = false; end; end
            if ~islogical(interpolate)
                error(['Incorrect value for property ''' varargin{i} ''' (type ''help <a href="matlab:help WaveletSpectrogram">WaveletSpectrogram</a>'' for details).']);
            end
        otherwise,
            error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help WaveletSpectrogram">WaveletSpectrogram</a>'' for details).']);
    end
end

t = lfp(:,1);
signal = lfp(:,2);
dt = mode(diff(t));

if isempty(range),
    % Provide the maximum computable range possible
    s0 = 2*dt;    % maximum allowable frequency at least twice the sampling frequency
    j1 = ceil(10/resolution);    % up to 10 powers-of-two with dj (=resolution) sub-octaves each
else
    s0 = 1/range(2);
    j1 = ceil(log2(range(2)./range(1))/resolution);
end

% Wavelet transform:
[wave,~,scale,~] = wavelet(signal,dt,pad,resolution,s0,j1,'Morlet');
power = flipud((abs(wave)).^2);        % compute wavelet power spectrum
f0 = 1./flipud(scale(:));

% Downsample to desired temporal resolution
if ~isempty(step),
    t0 = t;
    t = (t0(1):step:t(end))';
    power = interp1(t0,power',t)';
end

if interpolate
    % Get equally spaced frequencies (rather than logscale)
    f = (min(f0):min(diff(f0)):max(f0))';
    spectrogram = interp1(f0,power,f);
else
    spectrogram = power;
    f = f0;
end

% ------------------------------- Helper function -------------------------------

function [wave,period,scale,coi] = wavelet(Y,dt,pad,dj,s0,J1,mother,param)

if (nargin < 8), param = -1; end
if (nargin < 7), mother = -1; end
if (nargin < 6), J1 = -1; end
if (nargin < 5), s0 = -1; end
if (nargin < 4), dj = -1; end
if (nargin < 3), pad = 0; end
if (nargin < 2)
	error('Must input a vector Y and sampling time DT')
end

n1 = length(Y);

if (s0 == -1), s0=2*dt; end
if (dj == -1), dj = 1./4.; end
if (J1 == -1), J1=fix((log(n1*dt/s0)/log(2))/dj); end
if (mother == -1), mother = 'MORLET'; end

%....construct time series to analyze, pad if necessary
x(1:n1) = Y - mean(Y);
if (pad == 1)
	base2 = fix(log(n1)/log(2) + 0.4999);   % power of 2 nearest to N
	x = [x,zeros(1,2^(base2+1)-n1)];
end
n = length(x);

%....construct wavenumber array used in transform [Eqn(5)]
k = [1:fix(n/2)];
k = k.*((2.*pi)/(n*dt));
k = [0., k, -k(fix((n-1)/2):-1:1)];

%....compute FFT of the (padded) time series
f = fft(x);    % [Eqn(3)]

%....construct SCALE array & empty PERIOD & WAVE arrays
scale = s0*2.^((0:J1)*dj);
period = scale;
wave = zeros(J1+1,n);  % define the wavelet array
wave = wave + i*wave;  % make it complex

% loop through all scales and compute transform
for a1 = 1:J1+1
	[daughter,fourier_factor,coi,dofmin]=wave_bases(mother,k,scale(a1),param);
	wave(a1,:) = ifft(f.*daughter);  % wavelet transform[Eqn(4)]
end

period = fourier_factor*scale;
coi = coi*dt*[1E-5,1:((n1+1)/2-1),fliplr((1:(n1/2-1))),1E-5];  % COI [Sec.3g]
wave = wave(:,1:n1);  % get rid of padding before returning

function [daughter,fourier_factor,coi,dofmin] = wave_bases(mother,k,scale,param)

mother = upper(mother);
n = length(k);

if (strcmp(mother,'MORLET'))  %-----------------------------------  Morlet
	if (param == -1), param = 6.; end
	k0 = param;
	expnt = -(scale.*k - k0).^2/2.*(k > 0.);
	norm = sqrt(scale*k(2))*(pi^(-0.25))*sqrt(n);    % total energy=N   [Eqn(7)]
	daughter = norm*exp(expnt);
	daughter = daughter.*(k > 0.);     % Heaviside step function
	fourier_factor = (4*pi)/(k0 + sqrt(2 + k0^2)); % Scale-->Fourier [Sec.3h]
	coi = fourier_factor/sqrt(2);                  % Cone-of-influence [Sec.3g]
	dofmin = 2;                                    % Degrees of freedom
elseif (strcmp(mother,'PAUL'))  %--------------------------------  Paul
	if (param == -1), param = 4.; end
	m = param;
	expnt = -(scale.*k).*(k > 0.);
	norm = sqrt(scale*k(2))*(2^m/sqrt(m*prod(2:(2*m-1))))*sqrt(n);
	daughter = norm*((scale.*k).^m).*exp(expnt);
	daughter = daughter.*(k > 0.);     % Heaviside step function
	fourier_factor = 4*pi/(2*m+1);
	coi = fourier_factor*sqrt(2);
	dofmin = 2;
elseif (strcmp(mother,'DOG'))  %--------------------------------  DOG
	if (param == -1), param = 2.; end
	m = param;
	expnt = -(scale.*k).^2 ./ 2.0;
	norm = sqrt(scale*k(2)/gamma(m+0.5))*sqrt(n);
	daughter = -norm*(i^m)*((scale.*k).^m).*exp(expnt);
	fourier_factor = 2*pi*sqrt(2./(2*m+1));
	coi = fourier_factor/sqrt(2);
	dofmin = 1;
else
	error('Mother must be one of MORLET,PAUL,DOG')
end





