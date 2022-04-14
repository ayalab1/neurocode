function [peaktopeak troughs] = FindThetaCycles(lfp)

% FindThetaCycles: find intervals that qualify as theta cycles from the lfp signal
% 
% lfp should be in the [timestamps signal] format. Using CleanLFP before
% calling this function is recommended. The [start stop] intervals as well
% as the troughs take theta asymmetry into account.
%
% Copyright (C) 2018 by Ralitsa Todorova
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

thetaFiltered = FilterLFP(lfp, 'passband', 'theta');
[~, amplitude, ~] = Phase(thetaFiltered);

troughs = SineWavePeaks(thetaFiltered,'mode','troughs');
troughtotrough = [troughs(1:end-1) troughs(2:end)];

%% Shift peaks to avoid the bias of trying to fit a sine onto an asymmetric wave

% we work with broad bandpass filtered (1-80Hz like in Belluscio et al (2012) signal
filtered = FilterLFP(lfp, 'passband', [1 80]);
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
badsize = diff(peaktopeak,[],2) > 0.2 | diff(peaktopeak,[],2) < 0.1;
ok(badsize) = false;
% if at any point in the cycle, amplitude falls below 1 sds below the mean
nottheta = amplitude(zscore(amplitude(:,2))<-1,1); 
[l, nottheta] = InIntervals(nottheta, peaktopeak); % this is not a theta cycle
ok(unique(nottheta(l)),:) = false;

troughs = troughs(ok);
peaktopeak = peaktopeak(ok,:);
