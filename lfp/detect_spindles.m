function detect_spindles(varargin)
% detect_spindles: Detect spindles in Local Field Potential (LFP) data.
%
% This function detects spindles in LFP data by preprocessing the signal,
% selecting appropriate channels based on spindle power, and running the
% detection algorithm. It uses multiple parameters to customize the process,
% such as the brain region to analyze, EMG thresholds, and whether to restrict
% detection to specific intervals (e.g., NREM sleep or immobility).
%
% Usage:
%   detect_spindles('ParameterName', ParameterValue, ...)
%
% Parameters:
%   'basepath'             (string, default: current directory)
%       Path to the session directory containing the data.
%
%   'channel'              (numeric array, default: [])
%       Channel(s) to use for spindle detection (1-indexed). If empty, a
%       channel is automatically selected based on spindle power.
%
%   'peak_threshold'       (numeric, default: 5)
%       Minimum z-score of the peak amplitude for spindle detection.
%
%   'tail_threshold'       (numeric, default: 2.5)
%       Minimum z-score for the start and end of a spindle.
%
%   'brainRegion'          (string, char, or cell array of strings, default: {'PFC', 'MEC', 'LEC', 'EC', 'ILA', 'PL', 'Cortex', 'CTX'})
%       Specifies brain regions for channel selection. If no channel is
%       provided, the function selects channels from these regions.
%
%   'manual_pick_channel'  (logical, default: false)
%       If true, opens Neuroscope2 for manual channel selection if no channel
%       is provided. The user must tag the selected channel as 'spindle'.
%
%   'EMG_threshold'        (numeric, default: 0.6)
%       Threshold for EMG (electromyography) data. Periods with EMG below
%       this value are considered immobility periods for spindle detection.
%       Set to Inf to disable the immobility restriction.
%
%   'NREM_restrict'        (logical, default: true)
%       If true, restrict spindle detection to previously detected NREM
%       sleep intervals.
%
%   'passband'             (numeric array, default: [9, 17])
%       Frequency range (in Hz) used for filtering the signal when detecting
%       spindles.
%
%   'ignore_intervals'     (numeric array, default: [])
%       Intervals to ignore during spindle detection (e.g., optogenetic
%       stimulation periods).
%
%   'forceDetect'          (logical, default: false)
%       If true, forces spindle detection even if a detection file already
%       exists. Overwrites outdated detection files.
%
%   'save_filename'        (string, default: 'spindles')
%       Filename for saving spindle detection results.
%
% Outputs:
%   Results are saved to a file in the session directory in CellExplorer
%   event format. The output file includes spindle start/stop times, peak
%   amplitudes, durations, and metadata about the detection process.
%
% Example:
%   detect_spindles('basepath', '/path/to/data', 'brainRegion', 'PFC');
%
% Notes:
% - If no channel is provided, the function automatically selects a channel
%   based on spindle power within the specified brain regions.
% - If 'manual_pick_channel' is true, Neuroscope2 is launched to allow the
%   user to select a channel.
% - EMG and NREM restrictions help ensure that detected spindles occur during
%   immobility or specific sleep states.
%
% Dependencies:
% - Neuroscope2 (optional for manual channel selection).
% - CellExplorer event file format for saving results.


% parse inputs
p = inputParser;
p.PartialMatching = true; % allows the use of abbreviations (e.g. "EMG_thres" instead of "EMG_threshold"
addParameter(p, 'basepath', pwd, @(x) any([isfolder(x), iscell(x)]));
addParameter(p, 'channel', [], @isnumeric);
addParameter(p, 'peak_threshold', 5, @isnumeric);
addParameter(p, 'tail_threshold', 2.5, @isnumeric);
addParameter(p, 'brainRegion', {'PFC', 'MEC', 'LEC', 'EC', 'ILA', 'PL', 'Cortex', 'CTX'}, ...
    @(x) ischar(x) || isstring(x) || (iscell(x) && all(cellfun(@(y) ischar(y) || isstring(y), x))));
addParameter(p, 'manual_pick_channel', false, @islogical);
addParameter(p, 'EMG_threshold', 0.6, @isnumeric);
addParameter(p, 'NREM_restrict', true, @islogical);
addParameter(p, 'passband', [9, 17], @isnumeric);
addParameter(p, 'ignore_intervals', [], @isnumeric);
addParameter(p, 'forceDetect', false, @islogical);
addParameter(p, 'save_filename', 'spindles', @ischar);

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
peak_threshold = p.Results.peak_threshold;
tail_threshold = p.Results.tail_threshold;
passband = p.Results.passband;
channel = p.Results.channel;
NREM_restrict = p.Results.NREM_restrict;
manual_pick_channel = p.Results.manual_pick_channel;
ignore_intervals = p.Results.ignore_intervals;
forceDetect = p.Results.forceDetect;
save_filename = p.Results.save_filename;

disp(basepath)

basename = basenameFromBasepath(basepath);

% check if file aready exists
event_file = fullfile(basepath, [basename, '.', save_filename, '.events.mat']);
if exist(event_file, "file") && ~forceDetect
    % verify file contains minimum fields (timestamps, peaks)
    load(event_file, 'spindle')
    if isfield(spindle, 'timestamps') && isfield(spindle, 'peaks')
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

% checks for manually added spindle and if doesn't exist will open neuroscope
if manual_pick_channel
    if isfield(session.channelTags, 'spindle')
        channel = session.channelTags.spindle.channels;
        warning('using spindle channel that was already tagged')
    else
        disp("Pick spindle channel from neuroscope")
        disp("add channel tag as 'spindle'")
        NeuroScope2('basepath', basepath)
        % load session meta data
        load(fullfile(basepath, [basename, '.session.mat']), 'session')
        if isfield(session.channelTags, 'spindle')
            channel = session.channelTags.spindle.channels;
        else
            warning("you did not pick a channel, or did not name it 'spindle'")
        end
    end
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

    % get spindle power to choose channel
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
    % reshape interval if needed
    if size(immobility, 2) == 1
        immobility = immobility';
    end
else
    immobility = [-inf, inf];
end

% filter in spindle band
filtered = FilterLFP(clean_lfp, 'passband', passband);

% detect spindles, Nx4 matrix [start peak_t end peak_z]
spindles_ = FindSpindles(Restrict(filtered, immobility), ...
    'threshold', tail_threshold, 'peak', peak_threshold);

% put into table (like pandas but much worse)
spindles = array2table(spindles_, 'VariableNames', {'start', 'peak_t', 'stop', 'peak_z'});

% remove spindles that intersect with ignore intervals
if ~isempty(ignore_intervals)
    keep_intervals = ~IntervalsIntersect(spindles{:, {'start', 'stop'}}, ignore_intervals);
    spindles = spindles(keep_intervals, :);
end

% store to cell explorer event file
spindles_ = struct();
spindles_.timestamps = spindles{:, {'start', 'stop'}};
spindles_.peaks = spindles.peak_t;
spindles_.amplitude = spindles.peak_z;
spindles_.amplitudeUnits = 'zscore';
spindles_.eventID = [];
spindles_.eventIDlabels = [];
spindles_.eventIDbinary = false;
spindles_.center = median(spindles_.timestamps, 2);
spindles_.duration = spindles_.timestamps(:, 2) - spindles_.timestamps(:, 1);
spindles_.detectorinfo.detectorname = 'detect_spindles.m';
spindles_.detectorinfo.detectiondate = datetime("today");
spindles_.detectorinfo.detectionintervals = [clean_lfp(1, 1), clean_lfp(end, 1)];
spindles_.detectorinfo.detectionparms = p.Results;
try
    spindles_.detectorinfo.detectionchannel1 = channel;
    spindles_.detectorinfo.detectionchannel = channel - 1;
catch

end
spindles_.detectorinfo.threshold = [tail_threshold, peak_threshold];

spindles = spindles_;
save(event_file, 'spindles');
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
