function add_event_psth_cell_explorer(basepath)
% add_event_psth_cell_explorer: adds event psth after cell_metrics creation
%
% Here is a function to add the psths for your events to cell_metrics after
% you have already created your cell_metrics file. Normally, this is done during
% pre processing, but maybe you forgot to detect ripples or maybe you have
% a new type of event that you want to view in cell explorer


% set default params
parameters.excludeMetrics = {'none'};
parameters.metrics = 'all';
parameters.ignoreEventTypes = {'MergePoints'};
parameters.showFigures = false;
parameters.metricsToExcludeManipulationIntervals = {'waveform_metrics', ...
    'PCA_features', 'acg_metrics', 'monoSynaptic_connections', ...
    'theta_metrics', 'spatial_metrics', 'event_metrics', 'psth_metrics'};

% load session, cell metrics, and spikes
basename = basenameFromBasepath(basepath);

load(fullfile(basepath, [basename, '.session.mat']), 'session');
preferences = preferences_ProcessCellMetrics(session);

load(fullfile(basepath, [basename, '.cell_metrics.cellinfo.mat']), 'cell_metrics');

load(fullfile(basepath, [basename, '.spikes.cellinfo.mat']), 'spikes');

if ~isfield(spikes, 'sr')
    spikes.sr = session.extracellular.sr;
end

% generate psths for existing event files
if any(contains(parameters.metrics, {'event_metrics', 'all'})) && ~any(contains(parameters.excludeMetrics, {'event_metrics'}))
    field2remove = {'rippleCorrelogram', 'events', 'rippleModulationIndex', 'ripplePeakDelay'};
    test = isfield(cell_metrics, field2remove);
    cell_metrics = rmfield(cell_metrics, field2remove(test));

    eventFiles = dir(fullfile(basepath, [basename, '.*.events.mat']));
    eventFiles = {eventFiles.name};
    eventFiles(contains(eventFiles, parameters.ignoreEventTypes)) = [];

    if ~isempty(eventFiles)
        dispLog('Event metrics', basename)
        for i = 1:length(eventFiles)
            eventName = strsplit(eventFiles{i}, '.');
            eventName = eventName{end-2};
            eventOut = load(fullfile(basepath, eventFiles{i}));
            disp(['  Importing ', eventName]);
            if strcmp(fieldnames(eventOut), eventName)
                if isfield(preferences.psth, eventName) && isstruct(preferences.psth.(eventName))
                    psth_parameters = preferences.psth.(eventName);
                else
                    psth_parameters = preferences.psth;
                end

                PSTH = calc_PSTH(eventOut.(eventName), ...
                    spikes, ...
                    'binCount', psth_parameters.binCount, ...
                    'alignment', psth_parameters.alignment, ...
                    'binDistribution', psth_parameters.binDistribution, ...
                    'duration', psth_parameters.duration, ...
                    'smoothing', psth_parameters.smoothing, ...
                    'percentile', psth_parameters.percentile, ...
                    'intervals', psth_parameters.intervals, ...
                    'eventName', eventName, ...
                    'plots', parameters.showFigures);

                if size(PSTH.responsecurve, 2) == cell_metrics.general.cellCount
                    cell_metrics.events.(eventName) = num2cell(PSTH.responsecurve, 1);
                    cell_metrics.general.events.(eventName).event_file = eventFiles{i};
                    cell_metrics.general.events.(eventName).x_bins = PSTH.time * 1000;
                    cell_metrics.general.events.(eventName).x_label = 'Time (ms)';
                    cell_metrics.general.events.(eventName).alignment = PSTH.alignment;

                    cell_metrics.([eventName, '_modulationIndex']) = PSTH.modulationIndex;
                    cell_metrics.([eventName, '_modulationRatio']) = PSTH.modulationRatio;
                    cell_metrics.([eventName, '_modulationPeakResponseTime']) = PSTH.modulationPeakResponseTime;
                    cell_metrics.([eventName, '_modulationSignificanceLevel']) = PSTH.modulationSignificanceLevel;
                end
            end
        end
    end
end

% save results back to cell_metrics file
save(fullfile(basepath, [basename, '.cell_metrics.cellinfo.mat']), 'cell_metrics');

end