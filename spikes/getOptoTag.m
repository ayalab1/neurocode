function optoTag = getOptoTag(varargin)
% getOptoTag: optogenetic response + plotting (digital OR analog pulses) + SALT (Kepecs-style)
%
% Main features
%   • Digital pulses (getDigitalIn) OR Analog pulses (basename.pulses.events.mat)
%   • Restrict to MergePoints epochs (rows of MergePoints.timestamps)
%   • Raster aligned to pulse start with pulse underlines + start/stop ticks
%   • PSTH (Hz) plotted as its own panel (same width as raster)
%   • CellExplorer waveform + narrow ACG per unit (from day*.cell_metrics.cellinfo.mat)
%   • Population heatmaps (Hz and Z) in UID order, bone colormap
%   • SALT (Kepecs-style): first-spike latency within 10 ms bins
%       - Test bins aligned to stim onset OR poststim (after pulse end)
%       - Baseline bins are random (within same epoch windows) excluding:
%           [pulseStart - preStimExclude, pulseEnd + postStimExclude]
%       - Uses modified Jensen–Shannon divergence against baseline-derived null
%   • Per-unit SALT diagnostics:
%       - Latency histogram (baseline vs stim) same bins, different colors
%       - Divergence histogram (null vs stim-baseline) same bins, different colors
%   • One global scatter: waveform correlation r vs -log10(SALT p_emp)
%
% Example usage:
%   optoTag = getOptoTag('digitalCh',14,'epochs',[5 6 7 8],'winSize',0.015,...
%                        'doSALT',true,'saltMode','stim','showFigures',false);
%
%   optoTag = getOptoTag('useAnalog',true,'analogCh',1,'epochs',[5 6],...
%                        'winSize',0.02,'doSALT',true,'saltMode','poststim',...
%                        'postOffsetStart',0.010,'showFigures',false,'tag','sleep');

%% ---------------------- Parse inputs ----------------------
p = inputParser;

% Pulse source
addParameter(p,'useAnalog',false,@islogical);
addParameter(p,'analogCh',[],@isnumeric);          % filter pulses.analogChannel
addParameter(p,'eventGroupID',[],@isnumeric);      % filter pulses.eventGroupID
addParameter(p,'digitalCh',[],@isnumeric);         % required if useAnalog=false

% Core
addParameter(p,'spikes',[],@isstruct);
addParameter(p,'basepath',pwd,@ischar);
addParameter(p,'binSize',0.001,@isnumeric);
addParameter(p,'winSize',0.5,@isnumeric);          % scalar or [left right]
addParameter(p,'epochs',[],@isnumeric);            % MergePoints rows

% Plot/save controls
addParameter(p,'plotFigures',true,@islogical);     % create per-neuron plots
addParameter(p,'showFigures',true,@islogical);     % show or keep invisible
addParameter(p,'tag','',@ischar);                  % opto_tag_<tag>
addParameter(p,'savePath','',@ischar);             % if empty, save in basepath
addParameter(p,'saveMat',true,@islogical);

% Cosmetics
addParameter(p,'pulseAlpha',0.22,@isnumeric);      % transparency for pulse underlines/ticks

% SALT controls
addParameter(p,'doSALT',true,@islogical);
addParameter(p,'saltMode','stim',@ischar);         % 'stim' or 'poststim'
addParameter(p,'saltWin',0.010,@isnumeric);        % 10 ms window
addParameter(p,'saltBin',0.0005,@isnumeric);       % 0.5 ms bins
addParameter(p,'saltNumBaseline',300,@isnumeric);
addParameter(p,'preStimExclude',0.010,@isnumeric); % 10 ms BEFORE pulse start excluded
addParameter(p,'postStimExclude',0.100,@isnumeric);% 100 ms AFTER pulse end excluded
addParameter(p,'postOffsetStart',0,@isnumeric);    % for poststim: start after pulse end

parse(p,varargin{:});

useAnalog        = p.Results.useAnalog;
analogCh         = p.Results.analogCh;
eventGroupID     = p.Results.eventGroupID;
digitalCh        = p.Results.digitalCh;

spikes           = p.Results.spikes;
basepath         = p.Results.basepath;
binSize          = p.Results.binSize;
winSize          = p.Results.winSize;
epochs           = p.Results.epochs;

plotFigures      = p.Results.plotFigures;
showFigures      = p.Results.showFigures;
tag              = p.Results.tag;
savePath         = p.Results.savePath;
saveMat          = p.Results.saveMat;

pulseAlpha       = p.Results.pulseAlpha;

doSALT           = p.Results.doSALT;
saltMode         = lower(p.Results.saltMode);
saltWin          = p.Results.saltWin;
saltBin          = p.Results.saltBin;
saltNumBaseline  = p.Results.saltNumBaseline;
preStimExclude   = p.Results.preStimExclude;
postStimExclude  = p.Results.postStimExclude;
postOffsetStart  = p.Results.postOffsetStart;

if ~(strcmp(saltMode,'stim') || strcmp(saltMode,'poststim'))
    error('saltMode must be ''stim'' or ''poststim''.');
end

prevPath = pwd;
cd(basepath);

%% ---------------------- winSize handling ----------------------
if numel(winSize) == 2
    winL = abs(winSize(1));
    winR = abs(winSize(2));
else
    winL = winSize;
    winR = winSize;
end

%% ---------------------- Load session + MergePoints ----------------------
basename = basenameFromBasepath(basepath);

sessionFile     = fullfile(basepath,[basename '.session.mat']);
mergepointsFile = fullfile(basepath,[basename '.MergePoints.events.mat']);

if exist(sessionFile,'file')
    load(sessionFile,'session');
else
    error('Missing session file: %s', sessionFile);
end

if exist(mergepointsFile,'file')
    load(mergepointsFile,'MergePoints');
else
    error('Missing MergePoints file: %s', mergepointsFile);
end

%% ---------------------- Epoch windows (validated) ----------------------
nMP = size(MergePoints.timestamps,1);
if isempty(epochs)
    epochWindows = MergePoints.timestamps;
else
    epochs = epochs(:)'; % row vector
    bad = epochs < 1 | epochs > nMP | ~isfinite(epochs) | mod(epochs,1)~=0;
    if any(bad)
        if isfield(MergePoints,'foldernames')
            disp('MergePoints.foldernames (row order):');
            disp(MergePoints.foldernames(:));
        end
        error('Requested epochs [%s] but MergePoints.timestamps has %d rows (valid: 1..%d).', ...
            num2str(epochs), nMP, nMP);
    end
    epochWindows = MergePoints.timestamps(epochs,:);
end
nEpochs = size(epochWindows,1);

%% ---------------------- Output folder ----------------------
if isempty(savePath)
    rootSave = basepath;
else
    rootSave = savePath;
end

if isempty(tag)
    outdir = fullfile(rootSave,'opto_tag');
else
    outdir = fullfile(rootSave, ['opto_tag_' tag]);
end
if ~exist(outdir,'dir'), mkdir(outdir); end

%% ---------------------- Load spikes ----------------------
if isempty(spikes)
    try
        spikes = importSpikes('basepath', basepath);
    catch
        spikes = loadSpikes;
    end
end
nCells = length(spikes.UID);

for iC = 1:nCells
    if ~isempty(spikes.times{iC})
        spikes.times{iC} = sort(spikes.times{iC}(:));
    end
end

%% ---------------------- Load CellExplorer cell_metrics (optional) ----------------------
hasCellMetrics = false;
cell_metrics = [];
cmFile = dir('*cell_metrics*.mat');
if ~isempty(cmFile)
    try
        load(cmFile(1).name,'cell_metrics');
        hasCellMetrics = true;
        fprintf('Loaded cell_metrics: %s\n', cmFile(1).name);
    catch
        warning('Found cell_metrics but could not load it.');
    end
end

%% ---------------------- Load pulses (analog OR digital) ----------------------
pulses.timestamps = [];
pulses.channel    = [];
pulses.source     = ternary(useAnalog,'analog','digital');

if useAnalog
    pulsesFile = fullfile(basepath, [basename '.pulses.events.mat']);
    if ~exist(pulsesFile,'file')
        error('useAnalog=true but file not found: %s', pulsesFile);
    end
    S = load(pulsesFile,'pulses');
    if ~isfield(S,'pulses') || ~isfield(S.pulses,'timestamps')
        error('Pulses file did not contain pulses.timestamps');
    end
    pA = S.pulses;

    ts = pA.timestamps;
    if isempty(ts), error('pulses.timestamps is empty.'); end
    keep = true(size(ts,1),1);

    if ~isempty(analogCh)
        if ~isfield(pA,'analogChannel')
            error('Requested analogCh but pulses.analogChannel missing.');
        end
        keep = keep & ismember(pA.analogChannel(:), analogCh(:));
    end
    if ~isempty(eventGroupID)
        if ~isfield(pA,'eventGroupID')
            error('Requested eventGroupID but pulses.eventGroupID missing.');
        end
        keep = keep & ismember(pA.eventGroupID(:), eventGroupID(:));
    end

    ts = ts(keep,:);
    if isempty(ts), error('No analog pulses remain after filtering.'); end

    pulses.timestamps = ts;
    if isfield(pA,'analogChannel')
        pulses.channel = pA.analogChannel(keep);
    else
        pulses.channel = ones(size(ts,1),1);
    end

else
    if isempty(digitalCh)
        error('Must provide digitalCh when useAnalog=false.');
    end

    fs = session.extracellular.sr;
    digitalIn = getDigitalIn('all','fs',fs,'periodLag',.01);

    for ii = 1:length(digitalCh)
        ts = digitalIn.ints{digitalCh(ii)};
        if isempty(ts), continue; end
        if size(ts,1)==2 && size(ts,2)~=2
            ts = ts';
        end
        pulses.timestamps = [pulses.timestamps; ts];
        pulses.channel    = [pulses.channel; repmat(digitalCh(ii), size(ts,1), 1)];
    end

    if isempty(pulses.timestamps)
        error('No digital pulses found for requested digitalCh.');
    end
end

% sort pulses by start
[~, ord] = sort(pulses.timestamps(:,1));
pulses.timestamps = pulses.timestamps(ord,:);
pulses.channel    = pulses.channel(ord,:);

%% ---------------------- Restrict pulses to epochs ----------------------
[~, idx] = Restrict(pulses.timestamps, epochWindows);
pulses.timestamps = pulses.timestamps(idx,:);
pulses.channel    = pulses.channel(idx,:);

if isempty(pulses.timestamps)
    warning('No pulses remain after epoch restriction.');
    optoTag = struct();
    cd(prevPath);
    return
end

pulseStarts = pulses.timestamps(:,1);
pulseEnds   = pulses.timestamps(:,2);
pulseDur    = pulseEnds - pulseStarts;
pulseDurationMed  = median(pulseDur);
pulseDurationMean = mean(pulseDur);
nPulses = numel(pulseStarts);

%% ---------------------- PSTH timebase ----------------------
t = -winL : binSize : winR;
nBins = numel(t);
edgesPSTH = [t - binSize/2, t(end) + binSize/2]; % length nBins+1

%% ---------------------- Allocate outputs ----------------------
optoTag.response     = nan(nCells,nBins);
optoTag.responseZ    = nan(nCells,nBins);
optoTag.pulseRate    = nan(nCells,1);
optoTag.baselineRate = nan(nCells,1);
optoTag.isTagged     = nan(nCells,1);

% SALT outputs
optoTag.SALT = struct();
optoTag.SALT.doSALT          = doSALT;
optoTag.SALT.mode            = saltMode;
optoTag.SALT.saltWin         = saltWin;
optoTag.SALT.saltBin         = saltBin;
optoTag.SALT.nBaselineSets   = saltNumBaseline;
optoTag.SALT.preExclude      = preStimExclude;
optoTag.SALT.postExclude     = postStimExclude;
optoTag.SALT.postOffsetStart = postOffsetStart;

optoTag.SALT.p_emp   = nan(nCells,1);
optoTag.SALT.p_ks    = nan(nCells,1);
optoTag.SALT.D_jsd   = nan(nCells,1);
optoTag.SALT.R_stim  = nan(nCells,1);
optoTag.SALT.R_base  = nan(nCells,1);
optoTag.SALT.dR      = nan(nCells,1);
optoTag.SALT.nTrials = nan(nCells,1);
optoTag.SALT.lat_ms  = nan(nCells,1);
optoTag.SALT.jit_ms  = nan(nCells,1);
optoTag.SALT.note    = strings(nCells,1);

% SALT diagnostics (cell arrays)
optoTag.SALT.D_post_each = cell(nCells,1); % stim-baseline distances
optoTag.SALT.D_null      = cell(nCells,1); % baseline-baseline distances
optoTag.SALT.postLat_ms  = cell(nCells,1); % latencies, includes "no spike" at win+bin
optoTag.SALT.baseLat_ms  = cell(nCells,1); % one baseline sample (for plotting)

% SALT distributions stored per neuron
optoTag.SALT.Q_dist      = cell(nCells,1); % stim distribution
optoTag.SALT.P_ex_dist   = cell(nCells,1); % one baseline distribution
optoTag.SALT.edges_ms    = [];             % shared
optoTag.SALT.noSpike_ms  = [];

%% ---------------------- Compute PSTH per neuron (pooled pulses) ----------------------
for iCell = 1:nCells
    stAll = Restrict(spikes.times{iCell}, epochWindows);

    counts = zeros(1,nBins);
    for pIdx = 1:nPulses
        rel = stAll - pulseStarts(pIdx);
        rel = rel(rel >= -winL & rel <= winR);
        if ~isempty(rel)
            counts = counts + histcounts(rel, edgesPSTH);
        end
    end

    psth = counts / (nPulses * binSize); % Hz
    optoTag.response(iCell,:) = psth;

    baseIdx = t < 0;
    mu = mean(psth(baseIdx));
    sd = std(psth(baseIdx));
    if sd == 0 || isnan(sd)
        optoTag.responseZ(iCell,:) = nan(size(psth));
    else
        optoTag.responseZ(iCell,:) = (psth - mu) / sd;
    end

    durClip = min(pulseDurationMed, winR);
    duringIdx = (t >= 0) & (t <= durClip);
    beforeIdx = (t >= -min(pulseDurationMed, winL)) & (t < 0);

    optoTag.pulseRate(iCell)    = mean(psth(duringIdx));
    optoTag.baselineRate(iCell) = mean(psth(beforeIdx));
    optoTag.isTagged(iCell)     = mean(optoTag.responseZ(iCell,duringIdx), 'omitnan') > 2;
end

%% ---------------------- SALT: baseline pool + test starts ----------------------
baseStartPool = buildBaselineStartPoolInEpoch(epochWindows, pulses.timestamps, saltWin, preStimExclude, postStimExclude);

if strcmp(saltMode,'stim')
    testT0 = pulseStarts;
else
    testT0 = pulseEnds + postOffsetStart;
end

%% ---------------------- SALT per neuron ----------------------
if doSALT
    for iCell = 1:nCells
        stAll = Restrict(spikes.times{iCell}, epochWindows);

        salt = runSALT_kepecsSameEpoch_DIAG(stAll, testT0, baseStartPool, ...
            'saltWin',saltWin,'saltBin',saltBin,'nBaselineSets',saltNumBaseline);

        optoTag.SALT.p_emp(iCell)   = salt.p_emp;
        optoTag.SALT.p_ks(iCell)    = salt.p_ks;
        optoTag.SALT.D_jsd(iCell)   = salt.D_jsd_median;
        optoTag.SALT.R_stim(iCell)  = salt.R_post;
        optoTag.SALT.R_base(iCell)  = salt.R_base;
        optoTag.SALT.dR(iCell)      = salt.R_post - salt.R_base;
        optoTag.SALT.nTrials(iCell) = salt.nTrials;
        optoTag.SALT.lat_ms(iCell)  = salt.lat_median_ms;
        optoTag.SALT.jit_ms(iCell)  = salt.lat_sd_ms;
        optoTag.SALT.note(iCell)    = string(salt.note);

        % diagnostics
        optoTag.SALT.D_post_each{iCell} = salt.D_post_each;
        optoTag.SALT.D_null{iCell}      = salt.D_null;
        optoTag.SALT.postLat_ms{iCell}  = salt.postLat_ms;
        optoTag.SALT.baseLat_ms{iCell}  = salt.baseLat_ex_ms;

        % distributions
        optoTag.SALT.Q_dist{iCell}    = salt.Q_dist;
        optoTag.SALT.P_ex_dist{iCell} = salt.P_ex_dist;
        if isempty(optoTag.SALT.edges_ms)
            optoTag.SALT.edges_ms   = salt.edges_ms;
            optoTag.SALT.noSpike_ms = salt.noSpike_ms;
        end
    end
end

%% ---------------------- Waveform correlation r (each unit vs population mean) ----------------------
optoTag.waveform_r = nan(nCells,1);
if hasCellMetrics
    wfCell = cell(nCells,1);
    maxLen = 0;
    for iCell = 1:nCells
        [~, wf] = getWaveformForCell(cell_metrics, iCell);
        wfCell{iCell} = wf;
        maxLen = max(maxLen, numel(wf));
    end

    if maxLen > 0
        wfMat = nan(nCells, maxLen);
        good = false(nCells,1);
        for iCell = 1:nCells
            wf = wfCell{iCell};
            if ~isempty(wf)
                wf = wf(:)';
                wfMat(iCell,1:numel(wf)) = wf;
                good(iCell) = true;
            end
        end

        if any(good)
            wfRef = mean(wfMat(good,:), 1, 'omitnan');
            for iCell = 1:nCells
                if good(iCell)
                    r = corr(wfMat(iCell,:)', wfRef(:), 'rows','complete');
                    optoTag.waveform_r(iCell) = r;
                end
            end
        end
    end
end

%% ---------------------- Population heatmaps (epoch-averaged; NO SORTING) ----------------------
epochPSTH = nan(nCells, nBins, nEpochs);
pulseCountsPerEpoch = zeros(nEpochs,1);

for e = 1:nEpochs
    epWin = epochWindows(e,:);
    [~, pIdxE] = Restrict(pulses.timestamps, epWin);
    if isempty(pIdxE), continue; end

    pStartsE = pulses.timestamps(pIdxE,1);
    nPulsesE = numel(pStartsE);
    pulseCountsPerEpoch(e) = nPulsesE;

    for iCell = 1:nCells
        stE = Restrict(spikes.times{iCell}, epWin);

        countsE = zeros(1,nBins);
        for p = 1:nPulsesE
            rel = stE - pStartsE(p);
            rel = rel(rel >= -winL & rel <= winR);
            if ~isempty(rel)
                countsE = countsE + histcounts(rel, edgesPSTH);
            end
        end
        epochPSTH(iCell,:,e) = countsE / (nPulsesE * binSize);
    end
end

validEpochs = pulseCountsPerEpoch > 0;
if any(validEpochs)
    popPSTH = mean(epochPSTH(:,:,validEpochs), 3, 'omitnan');
else
    popPSTH = nan(nCells,nBins);
end

baseIdx = t < 0;
muB = mean(popPSTH(:,baseIdx), 2, 'omitnan');
sdB = std(popPSTH(:,baseIdx), 0, 2, 'omitnan');
sdB(sdB==0 | isnan(sdB)) = NaN;
popPSTH_Z = (popPSTH - muB) ./ sdB;

figVis = ternary(showFigures,'on','off');
figH = figure('Visible',figVis,'Color','w','Position',[150 150 1400 600]);

ax1 = subplot(1,2,1);
imagesc(ax1, t, 1:nCells, popPSTH);
axis(ax1,'tight'); set(ax1,'YDir','normal','TickDir','out');
xlabel(ax1,'Time (s)'); ylabel(ax1,'Neuron index (UID order)');
title(ax1,'Epoch-averaged PSTH (Hz)'); xline(ax1,0,'k--');
colorbar(ax1); colormap(ax1,'bone');

ax2 = subplot(1,2,2);
imagesc(ax2, t, 1:nCells, popPSTH_Z);
axis(ax2,'tight'); set(ax2,'YDir','normal','TickDir','out');
xlabel(ax2,'Time (s)'); ylabel(ax2,'Neuron index (UID order)');
title(ax2,'Epoch-averaged PSTH (Z; baseline t<0)'); xline(ax2,0,'k--');
colorbar(ax2); colormap(ax2,'bone');

saveas(figH, fullfile(outdir,'PopulationHeatmaps.png'));
close(figH);

%% ---------------------- Population scatter: waveform r vs SALT p_emp ----------------------
if doSALT
    figS = figure('Visible',figVis,'Color','w','Position',[200 200 950 420]);
    x = optoTag.waveform_r;
    y = optoTag.SALT.p_emp;

    yy = -log10(y);
    scatter(x, yy, 28, 'filled'); grid on;
    xlabel('Waveform correlation r (vs population mean)');
    ylabel('-log10(SALT p_{emp})');
    title('Waveform similarity vs SALT significance');
    yline(-log10(0.01),'k--'); % Kepecs threshold often shown
    saveas(figS, fullfile(outdir,'WaveformCorr_vs_SALT.png'));
    close(figS);
end

%% ---------------------- Save outputs ----------------------
optoTag.timestamps      = t;
optoTag.pulses          = pulses;
optoTag.winL            = winL;
optoTag.winR            = winR;
optoTag.binSize         = binSize;
optoTag.winSize         = winSize;
optoTag.epochs          = epochs;

optoTag.pulseDurations  = pulseDur;
optoTag.pulseDur_median = pulseDurationMed;
optoTag.pulseDur_mean   = pulseDurationMean;
optoTag.pulsesPerEpoch  = pulseCountsPerEpoch;

optoTag.basepath        = basepath;
optoTag.basename        = basename;
optoTag.outdir          = outdir;

if saveMat
    save(fullfile(outdir,'OptoTag.cellinfo.mat'),'optoTag');
end

%% ---------------------- Per-neuron wide figures ----------------------
if plotFigures
    pulseRGB = [0.2 0.4 1];

    for iCell = 1:nCells
        st = Restrict(spikes.times{iCell}, epochWindows);

        % raster points
        rast_x = [];
        rast_y = [];
        for pIdx = 1:nPulses
            rel = st - pulseStarts(pIdx);
            rel = rel(rel >= -winL & rel <= winR);
            if ~isempty(rel)
                rast_x = [rast_x; rel(:)];
                rast_y = [rast_y; pIdx*ones(numel(rel),1)];
            end
        end

        % PSTH
        psth = optoTag.response(iCell,:);
        psth_sm = smooth(psth, 11);

        % WIDE layout: keep adding panels to the right (no stacked column)
        fig = figure('Visible',figVis,'Color','w','Position',[40 120 3200 520]);

        ax_raster = axes('Position',[0.03 0.18 0.16 0.72]);
        ax_psth   = axes('Position',[0.21 0.18 0.16 0.72]);
        ax_wf     = axes('Position',[0.39 0.18 0.10 0.72]);
        ax_acg    = axes('Position',[0.51 0.18 0.10 0.72]);
        ax_lat    = axes('Position',[0.63 0.18 0.16 0.72]); % Kepecs middle: latency hist
        ax_dist   = axes('Position',[0.81 0.18 0.16 0.72]); % Kepecs bottom: divergence hist

        % -------- Raster with pulse underlines + start/stop ticks --------
        hold(ax_raster,'on');

        for pIdx = 1:nPulses
            relEnd = pulseEnds(pIdx) - pulseStarts(pIdx);

            % underline
            h1 = plot(ax_raster, [0 relEnd], [pIdx pIdx], 'LineWidth', 2);
            setLineAlpha(h1, pulseRGB, pulseAlpha);

            % tick at start
            h2 = plot(ax_raster, [0 0], [pIdx-0.25 pIdx+0.25], 'LineWidth', 1);
            setLineAlpha(h2, pulseRGB, pulseAlpha);

            % tick at end
            h3 = plot(ax_raster, [relEnd relEnd], [pIdx-0.25 pIdx+0.25], 'LineWidth', 1);
            setLineAlpha(h3, pulseRGB, pulseAlpha);
        end

        plot(ax_raster, rast_x, rast_y, 'k.', 'MarkerSize', 5);

        xlabel(ax_raster,'Time (s)');
        ylabel(ax_raster,'Pulse #');
        title(ax_raster, sprintf('CellIdx %d (UID %d) — %s', iCell, spikes.UID(iCell), ternary(useAnalog,'Analog','Digital')));
        xlim(ax_raster, [-winL winR]);
        ylim(ax_raster, [0 nPulses+1]);
        set(ax_raster,'TickDir','out','FontSize',12);
        box(ax_raster,'off');

        % -------- PSTH --------
        hold(ax_psth,'on');
        plot(ax_psth, t, psth_sm, 'k', 'LineWidth', 1.5);
        xline(ax_psth,0,'k--');

        % pulse duration bar at bottom
        durClip = min(pulseDurationMed, winR);
        yl = ylim(ax_psth);
        hbar = plot(ax_psth, [0 durClip], [yl(1) yl(1)], 'LineWidth', 5);
        setLineAlpha(hbar, pulseRGB, min(0.45, pulseAlpha+0.15));

        xlabel(ax_psth,'Time (s)');
        ylabel(ax_psth,'Firing rate (Hz)');
        title(ax_psth,'PSTH');
        xlim(ax_psth, [-winL winR]);
        set(ax_psth,'TickDir','out','FontSize',12);
        box(ax_psth,'off');

        % SALT text in top-right of PSTH
        if doSALT
            saltStr = sprintf([ ...
                'SALT (%s)\n' ...
                'p = %.3g\n' ...
                'D = %.4f\n' ...
                'Rstim = %.2f  Rbase = %.2f\n' ...
                'n = %d\n' ...
                'lat = %.2f \\pm %.2f ms'], ...
                saltMode, ...
                optoTag.SALT.p_emp(iCell), ...
                optoTag.SALT.D_jsd(iCell), ...
                optoTag.SALT.R_stim(iCell), optoTag.SALT.R_base(iCell), ...
                optoTag.SALT.nTrials(iCell), ...
                optoTag.SALT.lat_ms(iCell), ...
                optoTag.SALT.jit_ms(iCell));

            text(ax_psth, 0.98, 0.98, saltStr, ...
                'Units','normalized', ...
                'HorizontalAlignment','right', ...
                'VerticalAlignment','top', ...
                'FontSize',10, ...
                'Interpreter','none', ...
                'BackgroundColor','w', ...
                'Margin',6);
        end

        % -------- Waveform --------
        cla(ax_wf); hold(ax_wf,'on'); box(ax_wf,'off');
        if hasCellMetrics
            [wf_t_ms, wf] = getWaveformForCell(cell_metrics, iCell);
            if ~isempty(wf) && ~isempty(wf_t_ms)
                plot(ax_wf, wf_t_ms, wf, 'k','LineWidth',1.5);
                xlabel(ax_wf,'Time (ms)'); ylabel(ax_wf,'Amp');
                title(ax_wf,'Waveform');
            else
                text(ax_wf,0.5,0.5,'No waveform','HorizontalAlignment','center');
                axis(ax_wf,'off');
            end
        else
            text(ax_wf,0.5,0.5,'No cell\_metrics','HorizontalAlignment','center');
            axis(ax_wf,'off');
        end

        % -------- ACG (CellExplorer narrow) --------
        cla(ax_acg); hold(ax_acg,'on'); box(ax_acg,'off');
        if hasCellMetrics && isfield(cell_metrics,'acg') && isfield(cell_metrics.acg,'narrow') ...
                && size(cell_metrics.acg.narrow,2) >= iCell && ~isempty(cell_metrics.acg.narrow(:,iCell))
            acg_n = cell_metrics.acg.narrow(:,iCell);
            area(ax_acg, acg_n, 'FaceColor','k','EdgeColor','none');
            title(ax_acg,'ACG (narrow)');
            xlabel(ax_acg,'Bin'); ylabel(ax_acg,'Count');
        else
            text(ax_acg,0.5,0.5,'No ACG','HorizontalAlignment','center');
            axis(ax_acg,'off');
        end

        % -------- SALT latency histogram (baseline vs stim) --------
        cla(ax_lat); hold(ax_lat,'on'); box(ax_lat,'off');
        if doSALT && ~isempty(optoTag.SALT.postLat_ms{iCell}) && ~isempty(optoTag.SALT.baseLat_ms{iCell})

            stimLat = optoTag.SALT.postLat_ms{iCell}(:);
            baseLat = optoTag.SALT.baseLat_ms{iCell}(:);

            win_ms = 1000*saltWin;
            % convert "no spike" sentinel (>win_ms) to exactly win_ms for display
            stimLat(stimLat > win_ms) = win_ms;
            baseLat(baseLat > win_ms) = win_ms;

            edgesL = 0:(1000*saltBin):win_ms;
            if edgesL(end) < win_ms, edgesL = [edgesL win_ms]; end

            histogram(ax_lat, baseLat, edgesL, 'EdgeColor','none', 'FaceAlpha',0.35, 'FaceColor',[0.5 0.5 0.5]);
            histogram(ax_lat, stimLat, edgesL, 'EdgeColor','none', 'FaceAlpha',0.35, 'FaceColor',pulseRGB);

            xlim(ax_lat,[0 win_ms]);
            xlabel(ax_lat,'First-spike latency (ms)');
            ylabel(ax_lat,'Count');
            title(ax_lat,'Latency hist (baseline vs test)');
            legend(ax_lat, {'baseline','test'}, 'Location','best', 'Box','off');
        else
            text(ax_lat,0.5,0.5,'No SALT latencies','HorizontalAlignment','center');
            axis(ax_lat,'off');
        end

        % -------- SALT divergence histogram (null vs stim-baseline) --------
        cla(ax_dist); hold(ax_dist,'on'); box(ax_dist,'off');
        if doSALT && ~isempty(optoTag.SALT.D_null{iCell}) && ~isempty(optoTag.SALT.D_post_each{iCell})
            Dnull = optoTag.SALT.D_null{iCell}(:);
            Dsb   = optoTag.SALT.D_post_each{iCell}(:);

            allD = [Dnull; Dsb];
            allD = allD(~isnan(allD));
            if isempty(allD)
                text(ax_dist,0.5,0.5,'No SALT distances','HorizontalAlignment','center');
                axis(ax_dist,'off');
            else
                nbins = 30;
                edgesD = linspace(min(allD), max(allD), nbins+1);

                histogram(ax_dist, Dnull, edgesD, 'EdgeColor','none', 'FaceAlpha',0.35, 'FaceColor',[0.5 0.5 0.5]);
                histogram(ax_dist, Dsb,   edgesD, 'EdgeColor','none', 'FaceAlpha',0.35, 'FaceColor',pulseRGB);

                ylD = ylim(ax_dist);
                plot(ax_dist, [optoTag.SALT.D_jsd(iCell) optoTag.SALT.D_jsd(iCell)], ylD, 'k-', 'LineWidth', 1.5);

                title(ax_dist,'Divergence hist');
                xlabel(ax_dist,'JSD'); ylabel(ax_dist,'Count');
                legend(ax_dist, {'baseline-baseline (null)','test-baseline','D (median)'}, 'Location','best', 'Box','off');
            end
        else
            text(ax_dist,0.5,0.5,'No SALT distances','HorizontalAlignment','center');
            axis(ax_dist,'off');
        end

        saveas(fig, fullfile(outdir, sprintf('UID%03d_optoTag_WIDE.png', spikes.UID(iCell))));
        close(fig);
    end
end

cd(prevPath);
end

%% ========================= Helpers =========================

function out = ternary(cond, a, b)
if cond, out = a; else, out = b; end
end

function setLineAlpha(h, rgb, a)
% Try RGBA if supported; otherwise lighten the color to mimic transparency.
try
    h.Color = [rgb(:)' a];
catch
    h.Color = rgb(:)' + (1-rgb(:)')*0.55;
end
end

function [wf_t_ms, wf] = getWaveformForCell(cell_metrics, iCell)
wf_t_ms = [];
wf = [];

if ~isfield(cell_metrics,'waveforms') || ~isfield(cell_metrics.waveforms,'time')
    return
end

t = cell_metrics.waveforms.time;
if iscell(t)
    if numel(t) >= iCell && ~isempty(t{iCell})
        t = t{iCell};
    else
        t = t{1};
    end
end
wf_t_ms = 1000 * t(:)';

W = [];
if isfield(cell_metrics.waveforms,'raw') && numel(cell_metrics.waveforms.raw) >= iCell && ~isempty(cell_metrics.waveforms.raw{iCell})
    W = cell_metrics.waveforms.raw{iCell};
elseif isfield(cell_metrics.waveforms,'filt_zscored') && numel(cell_metrics.waveforms.filt_zscored) >= iCell && ~isempty(cell_metrics.waveforms.filt_zscored{iCell})
    W = cell_metrics.waveforms.filt_zscored{iCell};
end

if isempty(W), return; end

if ismatrix(W) && size(W,1) > 1 && size(W,2) > 1
    wf = mean(W,1);
else
    wf = W(:)';
end
end

function baseStarts = buildBaselineStartPoolInEpoch(epochWindows, pulseIntervals, winLen, preEx, postEx)
% Baseline start pool:
%   all t0 in epochWindows (1 ms grid) where [t0,t0+winLen] fits,
%   excluding overlap with expanded stim windows [pulseStart-preEx, pulseEnd+postEx].
dt = 0.001; % 1 ms grid
baseStarts = [];

for k = 1:size(epochWindows,1)
    a = epochWindows(k,1);
    b = epochWindows(k,2);
    if (b - a) <= winLen
        continue
    end
    t0 = (a:dt:(b-winLen))';
    baseStarts = [baseStarts; t0];
end

if isempty(baseStarts) || isempty(pulseIntervals)
    return
end

badInts = [pulseIntervals(:,1) - preEx, pulseIntervals(:,2) + postEx];

s = baseStarts;
e = baseStarts + winLen;
keep = true(size(baseStarts));

for i = 1:size(badInts,1)
    keep = keep & ~(s < badInts(i,2) & e > badInts(i,1));
end

baseStarts = baseStarts(keep);
end

function salt = runSALT_kepecsSameEpoch_DIAG(spikeTimes, testT0, baseStartPool, varargin)
% Kepecs-style SALT within same-epoch baseline pool + diagnostics
p = inputParser;
addParameter(p,'saltWin',0.010,@isnumeric);
addParameter(p,'saltBin',0.0005,@isnumeric);
addParameter(p,'nBaselineSets',300,@isnumeric);
parse(p,varargin{:});

saltWin = p.Results.saltWin;
saltBin = p.Results.saltBin;
nB      = p.Results.nBaselineSets;

spikeTimes = spikeTimes(:);
testT0     = sort(testT0(:));
nTrials    = numel(testT0);

salt = emptySalt();
salt.note = '';

if nTrials < 5
    salt.note = 'Too few trials for SALT.';
    return
end
if numel(baseStartPool) < nTrials
    salt.note = sprintf('Baseline pool too small (%d starts) for %d trials.', numel(baseStartPool), nTrials);
    return
end

% latency histogram edges + explicit "no spike" bin
edges = 0:saltBin:saltWin;
if edges(end) < saltWin, edges = [edges saltWin]; end
noSpikeVal = saltWin + saltBin; % sentinel (seconds)

% ---- TEST latencies ----
postLat = firstSpikeLatencies_sortedEvents(spikeTimes, testT0, saltWin); % seconds, NaN allowed
R_post  = mean(~isnan(postLat));

postVals = postLat; postVals(isnan(postVals)) = noSpikeVal; % seconds with sentinel
Q = latencyDistWithNoSpike(postVals, edges, noSpikeVal);

% ---- One representative baseline draw for plotting baseline latency hist ----
idxEx  = randsample(numel(baseStartPool), nTrials, true);
baseT0 = sort(baseStartPool(idxEx));
baseLat_ex = firstSpikeLatencies_sortedEvents(spikeTimes, baseT0, saltWin);
R_base = mean(~isnan(baseLat_ex));
baseVals_ex = baseLat_ex; baseVals_ex(isnan(baseVals_ex)) = noSpikeVal;
P_ex = latencyDistWithNoSpike(baseVals_ex, edges, noSpikeVal);

% ---- Many baseline sets ----
P = nan(nB, numel(Q));
D_post_each = nan(nB,1);

for b = 1:nB
    idx = randsample(numel(baseStartPool), nTrials, true);
    baseT0b = sort(baseStartPool(idx)); % critical for fast scan
    baseLat = firstSpikeLatencies_sortedEvents(spikeTimes, baseT0b, saltWin);
    baseVals = baseLat; baseVals(isnan(baseVals)) = noSpikeVal;

    P(b,:) = latencyDistWithNoSpike(baseVals, edges, noSpikeVal);
    D_post_each(b) = JSDiv(Q, P(b,:));
end

D_post = median(D_post_each,'omitnan');

% ---- Null distances: baseline vs baseline ----
Dnull = nan(nB,1);
for b = 1:nB
    b2 = randi(nB);
    while b2 == b
        b2 = randi(nB);
    end
    Dnull(b) = JSDiv(P(b,:), P(b2,:));
end

p_emp = (1 + sum(Dnull >= D_post,'omitnan')) / (1 + sum(~isnan(Dnull)));

% Optional KS (not Kepecs) vs representative baseline draw
[~, p_ks] = kstest2(postVals, baseVals_ex);

succ = postLat(~isnan(postLat));
lat_med_ms = 1000*median(succ,'omitnan');
lat_sd_ms  = 1000*std(succ,0,'omitnan');

salt.nTrials       = nTrials;
salt.R_post        = R_post;
salt.R_base        = R_base;
salt.D_jsd_median  = D_post;
salt.p_emp         = p_emp;
salt.p_ks          = p_ks;
salt.lat_median_ms = lat_med_ms;
salt.lat_sd_ms     = lat_sd_ms;

% diagnostics
salt.D_post_each   = D_post_each;
salt.D_null        = Dnull;
salt.postLat_ms    = 1000*postVals;      % includes no-spike sentinel at (saltWin+saltBin)*1000
salt.baseLat_ex_ms = 1000*baseVals_ex;   % includes no-spike sentinel

% distributions saved
salt.Q_dist     = Q;
salt.P_ex_dist  = P_ex;
salt.edges_ms   = 1000*edges;
salt.noSpike_ms = 1000*noSpikeVal;
end

function lat = firstSpikeLatencies_sortedEvents(spikeTimes, eventTimes, win)
% Fast single-pass; eventTimes should be sorted.
lat = nan(numel(eventTimes),1);
if isempty(spikeTimes) || isempty(eventTimes), return; end

spikeTimes = spikeTimes(:);
eventTimes = sort(eventTimes(:));

si = 1; nS = numel(spikeTimes);
for i = 1:numel(eventTimes)
    t0 = eventTimes(i);
    while si <= nS && spikeTimes(si) < t0
        si = si + 1;
    end
    if si > nS, break; end
    dt = spikeTimes(si) - t0;
    if dt >= 0 && dt <= win
        lat(i) = dt;
    end
end
end

function p = latencyDistWithNoSpike(vals, edges, noSpikeVal)
c = histcounts(vals(vals~=noSpikeVal), edges);
c_no = sum(vals==noSpikeVal);
p = [c(:)' c_no];
eps0 = 1e-12;
p = p + eps0;
p = p ./ sum(p);
end

function D = JSDiv(P,Q)
P = P(:)'; Q = Q(:)';
P = P./sum(P); Q = Q./sum(Q);
M = 0.5*(P+Q);
D = 0.5*KLDiv(P,M) + 0.5*KLDiv(Q,M);
end

function D = KLDiv(P,Q)
D = sum(P .* log2(P ./ Q));
end

function s = emptySalt()
s = struct('nTrials',0,'R_post',NaN,'R_base',NaN,'D_jsd_median',NaN,'p_emp',NaN,'p_ks',NaN,...
    'lat_median_ms',NaN,'lat_sd_ms',NaN,'note','');
end