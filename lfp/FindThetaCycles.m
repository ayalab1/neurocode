function [peaktopeak,troughs, amplitude] = FindThetaCycles(lfp,varargin)

%FindThetaCycles - Find intervals that qualify as theta cycles from the lfp signal

% lfp should be in the [timestamps signal] format. Using CleanLFP before
% calling this function is recommended. The [start stop] intervals as well
% as the troughs take theta asymmetry into account.
%
% It's recommended to provide lfp restricted to a behavior session (excluding
% sleep sessions). The function will take into account theta asymmetry to find
% the exact peak and trough timestamps as described by Belluscio et al (2012).
%
% USAGE
%    [peaktopeak,troughs] = FindThetaCycles(lfp,<options>)
%
% INPUTS
%    lfp                unfiltered LFP (one channel) to use, in [timestamps signal] format
%    <options>          optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties        Values
%    -------------------------------------------------------------------------
%     'minDuration'     minimum theta cycle duration (default = 0.1 seconds)
%     'maxDuration'     maximum theta cycle duration (default = 0.2 seconds)
%     'artefactThreshold' threshold to use to pass on to CleanLFP to exclude
%                       artifacts when computing theta amplitude (default = 5)
%     'baseline'        interval(s) of the behavior session, excluding sleep sessions,
%                       (there is no need to restrict to running epochs) provided in [start stop]
%                       matrix format. These intervals will be used to estimate
%                       the expected theta amplitude and compute amplitude thresholds.
%                       (default = [0 Inf]);
%     'marginAmplitude' the minimum theta amplitude (in sd-s). If the theta amplitude
%                       during a cycle is lower than this amplitude, the cycle will
%                       be discarded (default = -1). Note that this should be low
%                       because during a behavioral epoch, theta oscillations are
%                       expected (as the animal is running) in the majority of the
%                       session.
%    =========================================================================
%
% OUTPUT
%
%    peaktopeak      theta cycle [start stop] timestamps (interval between two peaks)
%    troughs         the timestamp of the trough corresponding to the theta cycles in peaktopeak
%
% SEE ALSO
%
%    See also auto_theta_cycles, CleanLFP
%
% Copyright (C) 2018-2022 by Ralitsa Todorova
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

%% Parsing parameters
broadPassband = [1 80]; % we work with broad bandpass filtered (1-80Hz like in Belluscio et al (2012) signal
p = inputParser;
addParameter(p,'maxDuration',0.2,@isnumeric);
addParameter(p,'minDuration',0.1,@isnumeric);
addParameter(p,'artefactThreshold',5,@isnumeric); % thresholds for CleanLFP to discount artefacts
addParameter(p,'amplitudeStdThreshold',-1,@isnumeric); % if at any point in the cycle, amplitude falls below this many sds below the mean, the cycle will be discarded
addParameter(p,'baseline',[0 Inf],@isnumeric); % the whole provided period will be used to compute the expected theta amplitude

parse(p,varargin{:})
maxDuration = p.Results.maxDuration;
minDuration = p.Results.minDuration;
artefactThreshold = p.Results.artefactThreshold;
amplitudeStdThreshold = p.Results.amplitudeStdThreshold;
baseline = p.Results.baseline;

%% Compute reference theta cycles
thetaFiltered = FilterLFP(lfp, 'passband', [4 12]);
[~, amplitude, ~] = Phase(thetaFiltered);
isBaseline = InIntervals(amplitude,baseline);
amplitudeThreshold = nanmean(amplitude(isBaseline,2)) + amplitudeStdThreshold*nanstd(amplitude(isBaseline,2));

troughs = SineWavePeaks(thetaFiltered,'mode','troughs');
troughtotrough = [troughs(1:end-1) troughs(2:end)];

%% Shift peaks to avoid the bias of trying to fit a sine onto an asymmetric wave

filtered = FilterLFP(lfp, 'passband', broadPassband);
t = filtered(:,1);

% the real peak is the lowest point between two troughs
maxima = FindLocalMaxima(filtered(:,2));
[~,w] = InIntervals(filtered(maxima),troughtotrough);
ok = w>0; maxima = maxima(ok); w = w(ok);
[~,keepers] = Accumulate(w,filtered(maxima,2),'mode','max');
maxima = maxima(keepers(~isnan(keepers)));
peaks = t(maxima);

peaktopeak = [peaks(1:end-1) peaks(2:end)];

% the real trough is the lowest point between two peaks
minima = FindLocalMinima(filtered(:,2));
[~,w] = InIntervals(filtered(minima),peaktopeak);
ok = w>0; minima = minima(ok); w = w(ok);
[~,keepers] = Accumulate(w,filtered(minima,2),'mode','min');
minima = minima(keepers(~isnan(keepers)));
troughs = t(minima);

% keep track of the theta cycles to keep
ok = true(length(peaktopeak),1);
% remove cycles that are too long or too short
badsize = diff(peaktopeak,[],2) > maxDuration | diff(peaktopeak,[],2) < minDuration;
ok(badsize) = false;
% apply minimum amplitude threshold
% the derivative threshold is set to Inf because theta is a slow signal
% and fast artefacts captured by the derivative are not relevant
[~,bad,~] = CleanLFP(lfp,'thresholds',[artefactThreshold Inf],'manual',false);
amplitude(bad,2) = nan;
nottheta = amplitude(~(amplitude(:,2)>amplitudeThreshold),1);
% intervals containing moments of low amplitude theta are not theta cycles
ok = CountInIntervals(nottheta, peaktopeak)==0 & ok;

troughs = troughs(ok);
peaktopeak = peaktopeak(ok,:);
amplitude = amplitude(ok,2);
end