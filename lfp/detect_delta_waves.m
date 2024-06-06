function detect_delta_waves(varargin)
% detect_delta_waves: Detect delta waves in LFP data.
%
% This function is a wrapper for the FindDeltaWaves function, which
% identifies delta waves in local field potential (LFP) data. It processes
% the LFP data, selects an appropriate channel based on delta power, and
% performs delta wave detection on the clean LFP data.
%
% Usage:
%   detect_delta_waves('ParameterName', ParameterValue, ...)
%
% Parameters:
%   - 'basepath' (string): The base directory path for the data. Default is
%       the current working directory.
%   - 'brainRegion' (string): The brain region to analyze (e.g., 'PFC').
%       Default is 'PFC', 'MEC', 'EC'.
%   - 'EMG_thres' (double): The threshold for EMG (Electromyography) data to
%       identify immobility periods. Default is 0.6.
%   - 'peak_to_trough_ratio' (double): The minimum peak-to-trough ratio
%       required for a delta wave to be considered. Default is 3.
%   - 'passband' (numeric array): The frequency passband for delta power
%       calculation. Default is [1, 6].
%   - 'channel' (int): The delta channel. Default is empty
%   - 'NREM_restrict' (logical): restrict to detected NREM epochs
%   - 'manual_pick_channel' (logical): will open neuroscope for you to pick channel
%   - 'use_sleep_score_delta_channel' (logical): use detected delta channel
%       from sleep score. Note/warning, may not be cortical channel.
%   - 'showfig' (logical): option to show delta/firing psth
%   - 'verify_firing' (logical): verify that firing rate is suppressed in
%       delta wave
%   - 'ignore_intervals' (numeric): set of intervals to ignore. Use case:
%       opto stim intervals to ignore
%
% Example:
%   detect_delta_waves('basepath', '/path/to/data', 'brainRegion', 'PFC');
%
% Depends on: neurocode
%
% Author: Ryan H
% Date: 2023

% parse inputs
p = inputParser;
addParameter(p, 'basepath', pwd, @(x) any([isfolder(x), iscell(x)]));
addParameter(p, 'brainRegion', {'PFC', 'MEC', 'EC', 'ILA', 'PL','Cortex'}, @(x) any(iscell(x), iscchar(x)));
addParameter(p, 'EMG_thres', 0.6, @isdouble);
addParameter(p, 'peak_to_trough_ratio', 3, @isdouble);
addParameter(p, 'passband', [1, 6], @isnumeric);
addParameter(p, 'channel', [], @isnumeric);
addParameter(p, 'NREM_restrict', true, @islogical);
addParameter(p, 'manual_pick_channel', false, @islogical);
addParameter(p, 'use_sleep_score_delta_channel', false, @islogical);
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
EMG_thres = p.Results.EMG_thres;
peak_to_trough_ratio = p.Results.peak_to_trough_ratio;
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

% detect delta waves
deltas0 = FindDeltaWaves(clean_lfp);

% keep waves during low emg
EMG = getStruct(basepath, 'EMG');
immobility = EMG.timestamps(FindInterval(EMG.data < EMG_thres));
immobility(diff(immobility, [], 2) < 1, :) = [];
deltas0 = Restrict(deltas0, immobility);

% restict to deltas that are above peak_to_trough_ratio
deltas = deltas0(deltas0(:, 5)-deltas0(:, 6) > peak_to_trough_ratio, :);

% remove deltas that intersect with ignore intervals
if ~isempty(ignore_intervals)
    keep_intervals = ~IntervalsIntersect(deltas(:, [1, 3]), ignore_intervals);
    deltas = deltas(keep_intervals, :);
end

if verify_firing
    % Verify that spiking decreases at the peak of the delta
    try
        [~, st] = importSpikes('basepath', basepath, 'brainRegion', brainRegion);
        if ~st.isempty
            [delta_psth, ts] = PETH(st.spikes, deltas(:, 2));
            z = delta_psth; z = nanzscore(Smooth(z,[0 1]),[],2);
            response = mean(z(:,InIntervals(ts,[-1 1]*0.05)),2);
            ordered = sortrows([deltas(:,5)-deltas(:,6) response],-1);
            peak_to_trough_ratio = ordered(find(Smooth(ordered(:,2),100)<=-0.5,1,'last'),1);
            deltas = deltas(deltas(:,5)-deltas(:,6) > peak_to_trough_ratio, :);

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
            warning("no cortical cells, impossible to verify firing")
        end
    catch
        warning("no cortical cells, impossible to verify firing")
    end
end

% store to cell explorer event file
deltaWaves.timestamps = deltas(:, [1, 3]);
deltaWaves.peaks = deltas(:, 2);
deltaWaves.amplitude = deltas(:, 5)-deltas(:, 6);
deltaWaves.amplitudeUnits = 'zscore';
deltaWaves.eventID = [];
deltaWaves.eventIDlabels = [];
deltaWaves.eventIDbinary = false;
deltaWaves.center = median(deltaWaves.timestamps, 2);
deltaWaves.duration = deltaWaves.timestamps(:, 2) - deltaWaves.timestamps(:, 1);
deltaWaves.detectorinfo.detectorname = 'detect_delta_waves.m';
deltaWaves.detectorinfo.detectiondate = datetime("today");
deltaWaves.detectorinfo.detectionintervals = [clean_lfp(1, 1), clean_lfp(end, 1)];
deltaWaves.detectorinfo.detectionparms = p.Results;
deltaWaves.detectorinfo.detectionchannel = channel - 1;
deltaWaves.detectorinfo.detectionchannel1 = channel;
deltaWaves.detectorinfo.peak_to_trough_ratio = peak_to_trough_ratio;

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