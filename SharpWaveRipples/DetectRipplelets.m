function ripplelets = DetectRipplelets(varargin)
    % DetectRipple Process ripple events and merge them based on inter-ripple interval.
    %   This function loads ripple events, merges ripples with inter-ripple interval
    %   less than the specified threshold, and saves the processed ripple data.
    %
    %   Parameters:
    %       'basepath' - Basepath of session file (default: pwd)
    %       'single_iri' - Inter-ripple interval threshold in seconds (default: 0.2)
    %       'savedata' - Flag to save the processed data (default: 1)
    %       'plot_wavelet' - Flag to plot the wavelet (default: false)
    %       'lfp_channel' - LFP channel number (required if plot_wavelet is true)
    %
    %
    %   Code by Wenbo Tang, adapted to Neurocode by Heathlarsson 2024
    %   Code adapted from Yamamoto and Tonegawa, 2017 Neuron
    %
    
    % Input parser
    p = inputParser;
    addParameter(p, 'basepath', pwd, @ischar);
    addParameter(p, 'single_iri', 0.2, @isnumeric);
    addParameter(p, 'savedata', 1, @isnumeric);
    addParameter(p, 'plot_wavelet', false, @islogical);
    addParameter(p, 'lfp_channel', [], @(x) isnumeric(x) && isscalar(x));

    parse(p, varargin{:});

    % Assign parsed input to variables, using defaults where necessary
    basepath = p.Results.basepath;
    single_iri = p.Results.single_iri;
    savedata = p.Results.savedata;
    plot_wavelet = p.Results.plot_wavelet;
    lfp_channel = p.Results.lfp_channel;

    if plot_wavelet && isempty(lfp_channel)
        error('lfp_channel must be provided if plot_wavelet is true');
    end

    % Load ripples
    basename = basenameFromBasepath(basepath);
    load([basename, '.ripples.events.mat']); % Load ripple events
    start = ripples.timestamps(:, 1);
    stop = ripples.timestamps(:, 1) + ripples.duration;

    % Merge the ripples with inter-ripple interval < single_iri
    while true
        single_idx = find(start(2:end) - start(1:end-1) < single_iri);
        if isempty(single_idx), break; end
        start(single_idx + 1) = [];
        stop(single_idx) = [];
    end

    % Process merged ripples
    ripcount = zeros(length(stop), 1);
    peak = zeros(length(stop), 1);
    amplitude = zeros(length(stop), 1);
    if isfield(ripples, 'peaksSw')
        amplitudeSw = zeros(length(stop), 1);
        peakSw = zeros(length(stop), 1);
    end

    for i = 1:length(stop)
        current_riptime = [start(i), stop(i)];
        ripple_idx = find(ripples.timestamps(:, 1) >= current_riptime(1) & ...
            (ripples.timestamps(:, 1) + ripples.duration) <= current_riptime(2));
        ripcount(i) = length(ripple_idx);
        peak(i) = nanmean(ripples.peaks(ripple_idx));
        amplitude(i) = max(ripples.amplitude(ripple_idx));
        if isfield(ripples, 'peaksSw')
            amplitudeSw(i) = max(ripples.amplitudeSw(ripple_idx));
            peakSw(i) = nanmean(ripples.peaksSw(ripple_idx));
        end
    end

    % Save results in a struct
    ripplelets.timestamps = [start, stop];
    ripplelets.duration = stop - start;
    ripplelets.peaks = peak';
    ripplelets.ripcount = ripcount';
    ripplelets.amplitude = amplitude';
    if isfield(ripples, 'peaksSw')
        ripplelets.peakSw = peakSw';
        ripplelets.amplitudeSw = amplitudeSw';
    end
    ripplelets.detectorinfo = ripples.detectorinfo;

    % Save the processed data if required
    if savedata
        save(fullfile(basepath, 'ripplelets.events.mat'), 'ripplelets');
    end

    % Plotting wavelet if requested
    if plot_wavelet
        if ~exist(fullfile(basepath, 'Ripple_Profile'), 'dir')
            mkdir(fullfile(basepath, 'Ripple_Profile'));
        end
        lfpRip = getLFP(lfp_channel, 'basepath', basepath);
        try 
            [wavAvg, lfpAvg] = eventWavelet(lfpRip, ripples.peaks(1:500), 'twin', [0.1 0.1]);
        catch
            [wavAvg, lfpAvg] = eventWavelet(lfpRip, ripples.peaks, 'twin', [0.1 0.1]);
        end
        saveas(gcf, fullfile(basepath, 'Ripple_Profile', 'swrWaveletSample.png'));
    end
end