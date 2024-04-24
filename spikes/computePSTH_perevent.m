function [PSTH] = computePSTH_perevent(event, spikes, varargin)
% COMPUTEPSTH_PEREVENT Calculate PSTH per event
%
%   [PSTH, index_abs] = COMPUTEPSTH_PEREVENT(event, spikes, varargin) calculates
%   the PSTH for each event using the given spikes and event data.
%
%   Required input arguments:
%   - event: struct containing event timestamps
%   - spikes: struct containing spike timestamps
%
%   Optional input arguments (specified as parameter-value pairs):
%   - binCount: number of bins for the PSTH (default: 100)
%   - alignment: alignment of time ('onset', 'center', 'peaks', 'offset') (default: 'onset')
%   - binDistribution: distribution of bins around events (default: [0.4, 0.2, 0.4])
%   - duration: duration of PSTH (default: 0.15 seconds)
%   - smoothing: Gaussian smoothing window size (default: 5)
%   - percentile: percentile of event duration distribution (default: 99)
%   - plots: flag to show plots (default: true)
%   - eventName: title used for plots
%   - maxWindow: maximum window size in seconds (default: 10)
%   - zscorePlot: flag to plot z-scored response (default: true)
%
%   Output:
%   - PSTH: struct containing PSTH data
%   - index_abs: absolute indices of units for plotting

% Parse input arguments
p = inputParser;
addParameter(p, 'binCount', 100, @isnumeric);
addParameter(p, 'alignment', 'onset', @ischar);
addParameter(p, 'binDistribution', [0.4, 0.2, 0.4], @isnumeric);
addParameter(p, 'duration', 0.15, @isnumeric);
addParameter(p, 'smoothing', 5, @isnumeric);
addParameter(p, 'percentile', 99, @isnumeric);
addParameter(p, 'plots', true, @islogical);
addParameter(p, 'eventName', '', @ischar);
addParameter(p, 'maxWindow', 10, @isnumeric);
addParameter(p, 'zscorePlot', true, @islogical);
parse(p, varargin{:});

binCount = p.Results.binCount;
alignment = p.Results.alignment;
binDistribution = p.Results.binDistribution;
duration = p.Results.duration;
smoothing = p.Results.smoothing;
percentile = p.Results.percentile;
plots = p.Results.plots;
eventName = p.Results.eventName;
maxWindow = p.Results.maxWindow;
zscorePlot = p.Results.zscorePlot;

% If no duration is given, an optimal duration is determined
if duration == 0
    durations = diff(event.timestamps');
    stim_duration = prctile(sort(durations), percentile);
    duration = min(max(round(stim_duration * 1000), 50) / 1000, maxWindow);
end

binSize = max(round(duration / binCount * 1000), 1) / 1000;

% Determine event alignment
switch alignment
    case 'onset'
        event_times = event.timestamps(:, 1);
        padding = binDistribution(1) / binDistribution(2) * duration;
        binsToKeep = int64(ceil(padding / binSize):ceil((duration * 2 + padding) / binSize));
    case 'center'
        event_times = mean(event.timestamps, 2);
        padding = 0;
        binsToKeep = 1:duration * 2 / binSize;
    case 'offset'
        event_times = event.timestamps(:, 2);
        padding = binDistribution(3) / binDistribution(2) * duration;
        binsToKeep = int64(ceil(padding / binSize):ceil((duration * 2 + padding) / binSize));
    case 'peaks'
        event_times = event.peaks;
        padding = 0;
        binsToKeep = 1:duration * 2 / binSize;
end

disp(['  ', num2str(length(event_times)), '  events, duration set to: ', num2str(duration), ' sec, aligned to ', alignment, ', with binsize: ', num2str(binSize)])

% Determining the bins interval for metrics
binsPre = 1:floor(binDistribution(1) * length(binsToKeep));
binsEvents = floor(binDistribution(1) * length(binsToKeep)) + 1:floor((binDistribution(1) + binDistribution(2)) * length(binsToKeep));
binsPost = floor((binDistribution(1) + binDistribution(2)) * length(binsToKeep)) + 1:length(binsToKeep);

% Initialize PSTH_out to store PSTH for each event
PSTH_out = nan(numel(binsToKeep), numel(event_times));

% Iterate over each event
for i = 1:numel(event_times)
    disp(['Processing event ', num2str(i), ' out of ', num2str(numel(event_times))]);
    % Concatenate spike times with the current event time
    spike_times = [vertcat(spikes.times{:}); event_times(i)];
    % Assign cluster indices to spikes and events
    spike_cluster_index = [ones(size(vertcat(spikes.times{:}))); 2 * ones(size(event_times))];
    % Compute CCG for the concatenated spike times
    [ccg, time] = CCG(spike_times, spike_cluster_index, 'binSize', binSize, 'duration', (duration + padding) * 2);
    % Store the normalized PSTH for the current event
    PSTH_out(:, i) = ccg(binsToKeep + 1, 2, 1) ./ numel(event_times) / binSize;
end

% Extract the time vector
time = time(binsToKeep + 1);

% Plot the PSTH
if plots
    figure;
    plot(time, mean(PSTH_out, 2), 'LineWidth', 2);
    xlabel('Time (s)');
    ylabel('Mean PSTH');
    title('Population PSTH');
end
% Initialize cell array to store spike times for each event
spike_times_in_event = cell(numel(event.timestamps), 1);

% Extract spike times for each event
for i = 1:numel(event.timestamps)
    % Extract spike times for the i-th event
    spike_times_in_event{i} = spikes.times{i};
end

% Plot the raster
figure;
hold on;
for i = 1:numel(spike_times_in_event)
    plot(spike_times_in_event{i}, i * ones(size(spike_times_in_event{i})), '.', 'MarkerSize', 10);
end
xlabel('Time (s)');
ylabel('Event Number');
title('Raster Plot');
% Assign PSTH_out to output variable PSTH
PSTH.PSTH_out = PSTH_out;
PSTH.time = time;

end

