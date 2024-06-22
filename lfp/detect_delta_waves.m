function detect_delta_waves(varargin)
% detect_delta_waves: Detect delta waves in LFP data.
%
% This function is a wrapper for the FindDeltaWaves function, which
% identifies delta waves in local field potential (LFP) data. It processes
% the LFP data, selects an appropriate channel based on delta power, and
% performs delta wave detection on the clean LFP data.
%
%  USAGE
%
%    detect_delta_waves('ParameterName', ParameterValue, ...)
%
%    =========================================================================
%     ParameterName    ParameterValue
%    -------------------------------------------------------------------------
%     'basepath'       a string specifying the directory path for the session
%                      data (default = current directory)
%     'channel'        a channel or a list of channels specifying the channel
%                      to use to detect delta waves (1-indexing expected)
%     'threshold'      a number specifying the minimum peak-to-trough amplitude
%                      (default = 3) in z-units.
%     'brainRegion'    a string or a cell of multiple strings specifying the
%                      brain region labels to analyze (default = {'PFC', 'IL',
%                      'PL','Cortex','MEC', 'EC'}. This will only be used to
%                      select a channel if no channel is provided by the user.
%                      In addition, if 'verify_firing' is true, spiking data
%                      from the selected regions will be used do detect the
%                      cortical firing rate suppression during delta waves
%                      to adjust the delta wave amplitude threshold.
%     'manual_pick_channel'    a boolean specifying whether to open Neuroscope2
%                      for the user to select a channel manually if no
%                      channel if no channel is provided by the user (default = false)
%     'use_sleep_score_delta_channel'   a boolean specifying whether to use
%                      the channel previously used for sleep scoring if no
%                      channel is provided by the user (default = false). Note
%                      that this may not be a cortical channel.
%     'EMG_threshold'  a number specifying the threshold for EMG (Electromyography)
%                      data. Periods below this threshold immobility periods
%                      will be used to detect delta waves (default = 0.6).
%                      Note: prividing "Inf" will waive the immobility requirement
%     'NREM_restrict'  a boolean specifying whether to restrict to previously
%                      detected NREM periods (default = true)
%     'passband'       the frequency passband to filter the signal when detecting
%                      delta waves (default = [1 6]) Hz.
%     'showfig'        a boolean specifying whether to show delta/firing
%                      psth (default = false)
%     'verify_firing'  a boolean specifying whether to increase the threshold
%                      to make sure cortical firing rate is suppressed
%                      during the detected events (default = false). This will
%                      load spiking data from the brain regions privided in
%                      'brainRegion'
%     'ignore_intervals'  set of intervals to ignore (e.g opto stim intervals)
%                      (default = empty)
%    =========================================================================
%
%  EXAMPLE
%   detect_delta_waves('basepath', 'Y:\OJRproject\OJR57\day12', 'brainRegion', 'PFC');
%
%  SEE
%
%    Requires the toolbox "neurocode"
% Example:
%   detect_delta_waves('basepath', '/path/to/data', 'brainRegion', 'PFC');
%
% Depends on: neurocode
%
% Copyright (C) 2023-2024 by Ryan Harvey, Ralitsa Todorova
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% parse inputs
p = inputParser;
p.PartialMatching = true; % allows the use of abbreviations (e.g. "EMG_thres" instead of "EMG_threshold"
addParameter(p, 'basepath', pwd, @(x) any([isfolder(x), iscell(x)]));
addParameter(p, 'channel', [], @isnumeric);
addParameter(p, 'peak_to_trough_ratio', 3, @isnumeric); % legacy name for "threshold"
addParameter(p, 'threshold', 3, @isnumeric);
addParameter(p, 'brainRegion', {'PFC', 'MEC', 'EC', 'ILA', 'PL', 'Cortex'}, @(x) any(iscell(x), iscchar(x)));
addParameter(p, 'manual_pick_channel', false, @islogical);
addParameter(p, 'use_sleep_score_delta_channel', false, @islogical);
addParameter(p, 'EMG_threshold', 0.6, @isnumeric);
addParameter(p, 'NREM_restrict', true, @islogical);
addParameter(p, 'passband', [1, 6], @isnumeric);
addParameter(p, 'showfig', false, @islogical);
addParameter(p, 'verify_firing', true, @islogical);
addParameter(p, 'ignore_intervals', [], @isnumeric);

parse(p, varargin{:})
basepath = p.Results.basepath;


% if single basepath, package in cell array
if ~iscell(basepath)
    basepath = {basepath};
end
% iterate over basepaths
for i = 1:length(basepath)
    run(basepath{i}, p)
end
end

function run(basepath, p)
% local function to carry out detection

brainRegion = p.Results.brainRegion;
EMG_threshold = p.Results.EMG_threshold;
threshold = p.Results.peak_to_trough_ratio; % legacy name for "threshold"
threshold = p.Results.threshold;
passband = p.Results.passband;
channel = p.Results.channel;
NREM_restrict = p.Results.NREM_restrict;
manual_pick_channel = p.Results.manual_pick_channel;
use_sleep_score_delta_channel = p.Results.use_sleep_score_delta_channel;
showfig = p.Results.showfig;
verify_firing = p.Results.verify_firing;
ignore_intervals = p.Results.ignore_intervals;

disp(basepath)

basename = basenameFromBasepath(basepath);

% check if file aready exists
event_file = fullfile(basepath, [basename, '.deltaWaves.events.mat']);
if exist(event_file, "file")
    % verify file contains minimum fields (timestamps, peaks)
    load(event_file, 'deltaWaves')
    if isfield(deltaWaves, 'timestamps') && isfield(deltaWaves, 'peaks')
        if verify_firing
            % if newer detector was used, firing was already verified
            try
                if deltaWaves.detectorinfo.detectionparms.verify_firing
                    return
                end
            catch
                warning("older file type... verifying spiking")
            end
            [~, st] = importSpikes('basepath', basepath, 'brainRegion', brainRegion);
            if st.isempty
                warning("no cortical cells, impossible to verify firing")
                return
            end
            [delta_psth, ts] = PETH(st.spikes, deltaWaves.peaks);
            delta_psth_avg = mean(delta_psth);
            if mean(delta_psth_avg) < mean(delta_psth_avg(ts > -0.1 & ts < 0.1))
                warning("Firing rate does not dip around delta, change paramaters")
                figure;
                subplot(2, 1, 1)
                plot(ts, delta_psth_avg)
                subplot(2, 1, 2)
                [~, idx] = sort(mean(delta_psth(:, ts > -0.1 & ts < 0.1), 2));
                PlotColorMap(Shrink(sortby(delta_psth, idx), 72, 1), 'x', ts);
                colormap("parula")
                keyboard
            end
        end
        return
    else
        warning('existing file is out of date, and will be overwritten')
    end
end

% load session meta data
load(fullfile(basepath, [basename, '.session.mat']), 'session')

if NREM_restrict
    load(fullfile(basepath, [basename, '.SleepState.states.mat']), 'SleepState')
    lfp_restrict = SleepState.ints.NREMstate;
    if isempty(lfp_restrict)
        warning("no NREM sleep, skipping session")
        return
    end
else
    lfp_restrict = [-inf, inf];
end

% checks for manually added delta and if doesn't exist will open neuroscope
if manual_pick_channel
    if isfield(session.channelTags, 'delta')
        channel = session.channelTags.delta.channels;
        warning('using delta channel that was already tagged')
    else
        disp("Pick delta channel from neuroscope")
        disp("add channel tag as 'delta'")
        NeuroScope2('basepath', basepath)
        % load session meta data
        load(fullfile(basepath, [basename, '.session.mat']), 'session')
        if isfield(session.channelTags, 'delta')
            channel = session.channelTags.delta.channels;
        else
            warning("you did not pick a channel, or did not name it 'delta'")
        end
    end
end

if use_sleep_score_delta_channel
    channel = SleepState.detectorinfo.detectionparms.SleepScoreMetrics.SWchanID + 1;
end

if isempty(channel) || length(channel) > 1

    if isempty(channel)
        % find nearest region name
        regions = fields(session.brainRegions);
        region_idx = contains(regions, brainRegion);
        if isempty(region_idx)
            warning("no channels for brainRegion")
            return
        end
        region_name = regions(region_idx);

        % find channels for this brain region
        channels = [];
        for region_i = 1:length(region_name)
            channels = [channels, session.brainRegions.(region_name{region_i}).channels];
        end

        % find and remove bad channels
        channels = channels(~ismember(channels, session.channelTags.Bad.channels));

        if isempty(channels)
            warning("no channels for brainRegion")
            return
        end
    else
        channels = channel;
    end


    [lfp, timestamps] = load_lfp(basepath, session, channels, lfp_restrict);

    % get delta power to choose channel
    try
        pBand = bandpower(double(lfp), ...
            session.extracellular.srLfp, passband);

        pTot = bandpower(double(lfp), ...
            session.extracellular.srLfp, ...
            [1, (session.extracellular.srLfp / 2) - 1]);
    catch
        for c = 1:size(lfp, 2)
            pBand(c) = bandpower(double(lfp(:, c)), ...
                session.extracellular.srLfp, passband);

            pTot(c) = bandpower(double(lfp(:, c)), ...
                session.extracellular.srLfp, ...
                [1, (session.extracellular.srLfp / 2) - 1]);
        end
    end
    % find max theta power, normalized by wide band
    [~, c_idx] = max(pBand./pTot);

    % only leave channel and load lfp
    channel = channels(:, c_idx);
    lfp = lfp(:, c_idx);
else
    [lfp, timestamps] = load_lfp(basepath, session, channel, lfp_restrict);
end

% remove artifacts from lfp
try
    [clean_lfp, ~, ~] = CleanLFP([timestamps, double(lfp)], ...
        'thresholds', [6, 10], 'manual', false);
catch
    % if CleanLFP fails, use original... CleanLFP needs testing
    clean_lfp = [timestamps, double(lfp)];
end
% remove lfp struct because we now have clean lfp
clear lfp

if EMG_threshold < Inf
    % keep waves during low emg
    EMG = getStruct(basepath, 'EMG');
    immobility = EMG.timestamps(FindInterval(EMG.data < EMG_threshold));
    immobility(diff(immobility, [], 2) < 1, :) = [];
else
    immobility = [-inf, inf];
end

% detect delta waves
deltas0 = FindDeltaWaves(Restrict(clean_lfp, immobility));

% restict to deltas that are above the amplitude threshold
deltas = deltas0(deltas0(:, 5)-deltas0(:, 6) > threshold, :);


% remove deltas that intersect with ignore intervals
if ~isempty(ignore_intervals)
    keep_intervals = ~IntervalsIntersect(deltas(:, [1, 3]), ignore_intervals);
    deltas = deltas(keep_intervals, :);
end

if verify_firing
    % Verify that spiking decreases at the peak of the delta
    [~, st] = importSpikes('basepath', basepath, 'brainRegion', brainRegion);
    if ~st.isempty
        [delta_psth, ts] = PETH(st.spikes, deltas(:, 2));
        % smooth delta psth
        delta_psth_smooth = delta_psth;
        delta_psth_smooth = nanzscore(Smooth(delta_psth_smooth, [0, 1]), [], 2);
        % extract delta response
        response = mean(delta_psth_smooth(:, InIntervals(ts, [-1, 1]*0.05)), 2);
        % sort responses and locate when value passes threshold
        ordered = sortrows([deltas(:, 5) - deltas(:, 6), response], -1);
        threshold = ordered(find(Smooth(ordered(:, 2), 100) <= -0.5, 1, 'last'), 1);
        % only keep deltas with sufficient response
        deltas = deltas(deltas(:, 5)-deltas(:, 6) > threshold, :);
        % verify that the dip in firing around delta
        delta_psth_avg = mean(delta_psth);
        if mean(delta_psth_avg) < mean(delta_psth_avg(ts > -0.1 & ts < 0.1))
            warning("Firing rate does not dip around delta, change paramaters")
            figure;
            subplot(2, 1, 1)
            plot(ts, delta_psth_avg)
            subplot(2, 1, 2)
            PlotColorMap(Shrink(sortby(delta_psth, -(deltas(:, 5) - deltas(:, 6))), 72, 1), 'x', ts);
            colormap("parula")
            keyboard
        end
        % Optionally, view the firing around the detected events
        if showfig
            figure;
            subplot(2, 1, 1)
            plot(ts, delta_psth_avg)
            subplot(2, 1, 2)
            PlotColorMap(Shrink(sortby(delta_psth, -(deltas(:, 5) - deltas(:, 6))), 72, 1), 'x', ts);
            colormap("parula")
        end
    else
        verify_firing = false;
        warning("no cortical cells, impossible to verify firing");
    end
end

% store to cell explorer event file
deltaWaves.timestamps = deltas(:, [1, 3]);
deltaWaves.peaks = deltas(:, 2);
deltaWaves.amplitude = deltas(:, 5) - deltas(:, 6);
deltaWaves.amplitudeUnits = 'zscore';
deltaWaves.values = deltas(:, 4:end); % the signal value of the [start, peak, trough] of each delta wave
deltaWaves.eventID = [];
deltaWaves.eventIDlabels = [];
deltaWaves.eventIDbinary = false;
deltaWaves.center = median(deltaWaves.timestamps, 2);
deltaWaves.duration = deltaWaves.timestamps(:, 2) - deltaWaves.timestamps(:, 1);
deltaWaves.detectorinfo.detectorname = 'detect_delta_waves.m';
deltaWaves.detectorinfo.detectiondate = datetime("today");
deltaWaves.detectorinfo.detectionintervals = [clean_lfp(1, 1), clean_lfp(end, 1)];
deltaWaves.detectorinfo.detectionparms = p.Results;
deltaWaves.detectorinfo.detectionparms.verify_firing = verify_firing;
deltaWaves.detectorinfo.detectionchannel = channel - 1;
deltaWaves.detectorinfo.detectionchannel1 = channel;
deltaWaves.detectorinfo.threshold = threshold;
save(event_file, 'deltaWaves');
end

function [lfp, restricted_timestamps] = load_lfp(basepath, session, channels, restrict_intervals)

fileName = fullfile(basepath, [basenameFromBasepath(basepath), '.lfp']);
if ~exist(fileName, 'file')
    fileName = fullfile(basepath, [basenameFromBasepath(basepath), '.eeg']);
end
nChannels = session.extracellular.nChannels;

filenamestruct = dir(fileName);
% determine number of bytes per sample
dataTypeNBytes = numel(typecast(cast(0, 'int16'), 'uint8'));
% Number of samples per channel
nSamp = filenamestruct.bytes / (nChannels * dataTypeNBytes);
% memory map file
mmf = memmapfile(fileName, 'Format', {'int16', [nChannels, nSamp], 'data'}, ...
    'Writable', false);

timestamps = (1:nSamp)' / session.extracellular.srLfp;

[restricted_timestamps, idx] = Restrict(timestamps, restrict_intervals);

lfp = mmf.Data.data(channels, idx)';
end
