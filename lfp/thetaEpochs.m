function [SleepState] = thetaEpochs(basepath, varargin)
% [SleepState] = thetaEpochs(basepath)
% Find theta and non-theta epochs using the theta ratio determined by the SleepScoreMaster
%
% Note:
%   - theta and non-theta are only determined during wake state
%   - non-theta is determined by high emg and low theta, so most likely chewing or grooming

% F.Sharif 2020
%    =========================================================================
%  USAGE
%
%INPUT
%   [basePath]         [directory: '/whatevetPath/baseName/']
%   [manual_thres]      [bool: to manually adjust thresholds]
%
%    =========================================================================
%
%OUTPUT
%   Appends the the data to the output of SleepScoreMaster (SleepState.mat)
%   =========================================================================

% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

p = inputParser;
addParameter(p, 'manual_thres', false, @islogical);
parse(p, varargin{:});
manual_thres = p.Results.manual_thres;

basename = basenameFromBasepath(basepath);
sleepstatepath = fullfile(basepath, [basename, '.SleepState.states.mat']);
load(sleepstatepath, 'SleepState')

thratio = SleepState.detectorinfo.detectionparms.SleepScoreMetrics.thratio;
ththr = SleepState.detectorinfo.detectionparms.SleepScoreMetrics.histsandthreshs.THthresh;
emg = SleepState.detectorinfo.detectionparms.SleepScoreMetrics.EMG;
emgthr = SleepState.detectorinfo.detectionparms.SleepScoreMetrics.histsandthreshs.EMGthresh;

% determine high theta and also high emg
Theta_NDX = thratio > ththr & emg > emgthr;

A(:, 1) = SleepState.idx.states;
A(Theta_NDX, 2) = 7;

% non theta is determined by wake state and low theta
Non_thetaNDX = A(:, 1) == 1 & A(:, 2) == 0;
A(Non_thetaNDX, 2) = 9;

if manual_thres
    figure();
    ts = SleepState.detectorinfo.detectionparms.SleepScoreMetrics.t_clus;
    theta_intervals = ts(findIntervals(A(:, 2) == 7));

    ax1 = subplot(2, 1, 1);
    plot(ts, thratio)
    yline(ththr, 'color', 'w')
    axis tight
    ylabel('theta ratio')
    PlotIntervals(theta_intervals)

    ax2 = subplot(2, 1, 2);
    plot(ts, emg)
    yline(emgthr, 'color', 'w')
    axis tight
    ylabel('emg')
    xlabel('time (s)')
    PlotIntervals(theta_intervals)

    darkBackground(gcf, [0.1, 0.1, 0.1], [0.7, 0.7, 0.7])
    linkaxes([ax1, ax2], 'x');
    disp('Feel free to change "ththr" and/or emgthr then type "dbcont" to continue.')
    keyboard;

    % determine high theta and also high emg
    Theta_NDX = thratio > ththr & emg > emgthr;

    A(Theta_NDX, 2) = 7;
end

% store results
SleepState.idx.theta_epochs.states = A(:, 2);
SleepState.idx.theta_epochs.timestamps = SleepState.idx.timestamps;
SleepState.idx.theta_epochs.statenames{7} = 'THETA';
SleepState.idx.theta_epochs.statenames{9} = 'nonTHETA';
[INT] = bz_IDXtoINT(SleepState.idx.theta_epochs);
SleepState.ints.THETA = INT.THETAstate;
SleepState.ints.nonTHETA = INT.nonTHETAstate;

% save back to original structure
save(sleepstatepath, 'SleepState');
end