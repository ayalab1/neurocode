function [mean_psth, sem_psth, psth_window, raster] = computePSTH_perevent(event_matrix, spikes, varargin)
    p = inputParser;

    % Define input parameters
    addParameter(p, 'binCount', 100, @isnumeric);        % Number of bins
    addParameter(p, 'duration', 0.15, @isnumeric);       % Duration of PSTH (half window, in seconds)
    addParameter(p, 'plots', true, @islogical);          % Show plots?
    addParameter(p, 'eventName', '', @ischar);           % Event name for plots

    parse(p, varargin{:});

    % Retrieve parsed parameters
    binCount = p.Results.binCount;
    duration = p.Results.duration;
    plots = p.Results.plots;
    eventName = p.Results.eventName;

    % Calculate PSTH for each event
    event_count = size(event_matrix, 1);
    psth_matrix = zeros(event_count, binCount);
    raster = cell(event_count, 1);

    % Define the PSTH time window
    psth_window = linspace(-duration/2, duration/2, binCount);

    for event_idx = 1:event_count
        event_time = event_matrix(event_idx, :);
        spikes_in_event = spikes.times{1}(spikes.times{1} >= event_time(1) & spikes.times{1} <= event_time(2));
        spike_counts = histcounts(spikes_in_event, binCount, 'BinLimits', [-duration/2, duration/2]);
        psth_matrix(event_idx, :) = spike_counts / (duration / binCount);
        raster{event_idx} = spikes_in_event;
    end

    % Compute mean PSTH and SEM across all events
    mean_psth = mean(psth_matrix, 1);
    sem_psth = std(psth_matrix, 1) / sqrt(event_count);

    % Plot mean PSTH with SEM (optional)
    if plots
        figure;
        subplot(2, 1, 1);
        plot(psth_window, mean_psth, 'LineWidth', 2);
        hold on;
        shadedErrorBar(psth_window, mean_psth, sem_psth, '-b');
        xlabel('Time (s)');
        ylabel('Firing Rate (spikes/s)');
        title(['Mean PSTH for ' eventName]);
        legend('Mean', 'SEM');

        % Plot raster for each event
        subplot(2, 1, 2);
        for event_idx = 1:event_count
            scatter(raster{event_idx}, ones(size(raster{event_idx})) * event_idx, 'k', '.');
            hold on;
        end
        ylim([0.5, event_count + 0.5]);
        xlabel('Time (s)');
        ylabel('Event');
        title('Raster Plot');
    end
end
