function [EMGFromLFP] = getEMGFromLFP(basepath, varargin)
% Based on Erik Schomburg's work and code.  Grabs channels and calculates
% their correlations in the 300-600Hz band over sliding windows of 0.5sec.

% Channels are automatically selected and are a combination of first and last channels
% on each shank.  This is based on the xml formatting standard that channel ordering goes
% from superficial to deep for each channel/spike group.
%
% Special channels should be 0-indexed, per neuroscope convention
% Requires .lfp/lfp and .xml.  Assumes each spikegroup in the .xml
% represents a "shank"
%
% Mean pairwise correlations are calculated for each time point.

%   =========================================================================

%USAGE
% [EMGCorr] = bz_EMGCorrFromLFP(basePath)
%
% INPUTS
%       basePath            - string combination of basepath and basename of recording
%                             example: '/animal/recording/recording01'
%
%   (Optional)
%       'restrict'          - interval of time (relative to recording) to sleep score
%                            default = [0 inf]
%       'specialChannels'   - vector of 'special' channels that you DO want to use for EMGCorr calc (will be added to those auto-selected by spike group)
%       'rejectChannels'    - vector of 'bad' channels that you DO NOT want to use for EMGCorr calc
%       'restrictChannels'  - use only these channels (Neuroscope numbers)
%       'saveMat'           - true/false - default:true
%       'saveLocation'      - default: basePath
%       'overwrite'         - true/false - overwrite saved EMGFromLFP.LFP.mat
%                             default: false
%       'samplingFrequency' - desired sampling rate for EMG output. default:2
%       'noPrompts'     (default: false) prevents prompts about saving/adding metadata
%       'fromDat'           -uses the .dat file instead of .lfp (default:false)
%       'maxAutoChannels'   - maximum number of spike groups sampled when
%                             channels are picked automatically (one channel
%                             per group). Recordings with more usable spike
%                             groups use this many groups spread evenly
%                             across the spike-group list, because the cost
%                             grows with the number of channel pairs.
%                             Recordings with at most this many usable
%                             groups are unaffected.
%                             specialChannels are always added on top.
%                             Use inf to disable the cap. default: 16
%
%
% OUTPUTS
%
%       EMGFromLFP              - struct of the LFP datatype. saved at
%                               basePath/baseName.EMGFromLFP.LFP.mat
%          .timestamps          - timestamps (in seconds) that match .data samples
%          .data                - correlation data
%          .channels            - channel #'s used for analysis
%          .detectorName        - string name of function used
%          .samplingFrequency   - 1 / sampling rate of EMGCorr data
%          .detectorinfo        - how channels were selected (method,
%                                 maxAutoChannels, spike groups used, ...)
%
%   =========================================================================

% Erik Schomburg, Brendon Watson, Dan Levenstein, David Tingley, 2017
% Updated: Rachel Swanson 5/2017

%% Buzcode name of the EMGCorr.LFP.mat file
recordingname = basenameFromBasepath(basepath);
matfilename = fullfile(basepath, [recordingname, '.EMGFromLFP.LFP.mat']);

%% xmlPameters
p = inputParser;
addParameter(p, 'restrict', [0, inf], @isnumeric)
addParameter(p, 'specialChannels', [], @isnumeric)
addParameter(p, 'rejectChannels', [], @isnumeric)
addParameter(p, 'restrictChannels', [], @isnumeric)
addParameter(p, 'saveMat', true, @islogical)
addParameter(p, 'saveLocation', '', @isstr)
addParameter(p, 'overwrite', false, @islogical)
addParameter(p, 'samplingFrequency', 2, @isnumeric)
addParameter(p, 'noPrompts', false, @islogical);
addParameter(p, 'fromDat', false, @islogical);
addParameter(p, 'chInfo', [], @isstruct);
addParameter(p, 'maxAutoChannels', 16, @(x) isnumeric(x) && isscalar(x) && x >= 2);
parse(p, varargin{:})

restrict = p.Results.restrict;
specialChannels = p.Results.specialChannels;
rejectChannels = p.Results.rejectChannels;
restrictChannels = p.Results.restrictChannels;
saveMat = p.Results.saveMat;
overwrite = p.Results.overwrite;
samplingFrequency = p.Results.samplingFrequency;
noPrompts = p.Results.noPrompts;
fromDat = p.Results.fromDat;
chInfo = p.Results.chInfo;
maxAutoChannels = p.Results.maxAutoChannels;

if ~isempty(p.Results.saveLocation)
    matfilename = fullfile(p.Results.saveLocation, [recordingname, '.EMGFromLFP.LFP.mat']);
end

%% Check if EMGCorr has already been calculated for this recording
%If the EMGCorr file already exists, load and return with EMGCorr in hand
if exist(matfilename, 'file') && ~overwrite
    disp('EMGFromLFP Correlation already calculated - loading from EMGFromLFP.LFP.mat')
    load(matfilename)
    if exist('EMGCorr', 'var') %for backcompatability
        EMGFromLFP = EMGCorr;
    end
    if ~exist('EMGFromLFP', 'var')
        disp([matfilename, ' does not contain a variable called EMGFromLFP'])
    end
    return
end
disp('Calculating EMGFromLFP from High Frequency LFP Correlation')

load(fullfile(basepath, [recordingname, '.session.mat']), 'session')
nChannels = session.extracellular.nChannels;
SpkGrps = session.extracellular.spikeGroups.channels;
Fs = session.extracellular.srLfp;
if exist(fullfile(basepath, [recordingname, '.lfp']), 'file')
    lfpFile = checkFile('basepath', basepath, 'fileTypes', {'.lfp', '.eeg'});
    lfpFile = [basepath, filesep, lfpFile(1).name];
elseif exist(fullfile(basepath, [recordingname, '.eeg']), 'file')
    lfpFile = [basepath, filesep, recordingname, '.eeg'];
end
if fromDat
    datFile = checkFile('basepath', basepath, 'fileTypes', {'.dat'});
    datFile = [basepath, filesep, datFile(1).name];
    datFs = session.extracellular.sr;
end

%%

binScootS = 1 ./ samplingFrequency;
binScootSamps = round(Fs*binScootS); % must be integer, or error on line 190

%% Pick channels and to analyze
% get spike groups,
% pick every other one... unless specialshanks, in which case pick non-adjacent
%This is potentially dangerous in combination with rejectChannels... i.e.
%what if you pick every other shank but then the ones you pick are all
%reject because noisy shank.

% xcorrs_chs is a list of channels that will be loaded
% spkgrpstouse is a list of spike groups to find channels from

nUsableSpikeGroups = NaN;
spkgrpsused = [];
if ~isempty(restrictChannels) % If restrict channel case:
    xcorr_chs = restrictChannels;
    selectionMethod = 'restrictChannels';
else
    % get list of spike groups (aka shanks) that should be used
    usablechannels = [];
    spkgrpstouse = [];

    for gidx = 1:length(SpkGrps)
        usableshankchannels{gidx} = setdiff(SpkGrps{gidx}, rejectChannels);
        usablechannels = cat(2, usablechannels, usableshankchannels{gidx});
        if ~isempty(usableshankchannels{gidx})
            spkgrpstouse = cat(2, spkgrpstouse, gidx);
        end
    end

    % okay, lets limit usableshankchannels with spkgrpstouse. Thus, empty shanks will be excluded.
    usableshankchannels = usableshankchannels(spkgrpstouse);
    nUsableSpikeGroups = length(spkgrpstouse);

    % The number of channel pairs, and with it the run time, grows
    % quadratically with the number of spike groups. Keep at most
    % maxAutoChannels groups, spread evenly over the spike-group list.
    % Recordings with <= maxAutoChannels usable groups take the unchanged
    % path below.
    if nUsableSpikeGroups > maxAutoChannels
        keepgroups = unique(round(linspace(1, nUsableSpikeGroups, floor(maxAutoChannels))));
        usableshankchannels = usableshankchannels(keepgroups);
        spkgrpstouse = spkgrpstouse(keepgroups);
        selectionMethod = 'onePerSpikeGroupCapped';
    elseif length(usableshankchannels) > 1
        selectionMethod = 'onePerSpikeGroup';
    else
        selectionMethod = 'fiveFromSingleSpikeGroup';
    end
    spkgrpsused = spkgrpstouse;

    % get number of channels per shank. If one valid shank, choose 5 channels
    if length(usableshankchannels) > 1
        n = 1;
    else
        n = 5;
    end

    % get list of channels (1 from each good spike group)
    xcorr_chs = [];
    for gidx = 1:length(usableshankchannels)

        %grab random channel on each shank
        if ~isempty(usableshankchannels{gidx})
            if n == 1
                % sample with replacement for case of multi-shank
                randChfromShank = usableshankchannels{gidx}(randi(length(usableshankchannels{gidx}), n, 1));
            elseif n > 1
                % sample without replacement for case of 1 shank
                randChfromShank = usableshankchannels{gidx}(randperm(length(usableshankchannels{gidx}), n));
            end
            xcorr_chs = [xcorr_chs, randChfromShank];

        end
    end
    xcorr_chs = unique([xcorr_chs, specialChannels]);
end

%% Read and filter channel
% read channels
% old: bz_sessionInfo is 0 indexed (neuroscope) channels,
% but Loadbinary.m takes 1-indexed channel #'s
% new: should use bz_GetLFP

switch fromDat
    case false
        lfp = loadBinary(lfpFile, 'nChannels', nChannels, 'channels', xcorr_chs, ...
            'start', restrict(1), 'duration', diff(restrict), 'frequency', Fs); %read and convert to mV

        if isempty(lfp)
            error('LFP file empty. Must have valid binary file to compute EMG.')
        end

    case true
        lfp = loadBinary(datFile, 'nChannels', nChannels, 'channels', xcorr_chs, ...
            'start', restrict(1), 'duration', diff(restrict), 'frequency', datFs, ...
            'downsample', datFs./Fs); %read and convert to mV
        if isempty(lfp)
            error('DAT file empty. Must have valid binary file to compute EMG.')
        end
end

% Filter first in high frequency band to remove low-freq physiologically
% correlated LFPs (e.g., theta, delta, SPWs, etc.)

maxfreqband = floor(max([625, Fs / 2]));
% xcorr_freqband = [275 300 600 625]; % Hz
xcorr_freqband = [275, 300, maxfreqband - 25, maxfreqband]; % Hz
lfp = filtsig_in(lfp, Fs, xcorr_freqband);

%% xcorr 'strength' is the summed correlation coefficients between channel
% pairs for a sliding window of 25 ms
xcorr_window_samps = round(binScootS*Fs);
xcorr_window_inds = -xcorr_window_samps:xcorr_window_samps; %+- that number of ms in samples

% Pearson correlation of every channel pair in each window, summed over
% pairs. Windows are processed in blocks of bins: each channel's windows are
% centered and scaled once per block, so a pair costs one elementwise
% product instead of a per-bin corr() call. Same formula as corr (NaN for a
% flat window), results agree with the former per-bin loop to rounding error.
timestamps = (1 + xcorr_window_inds(end)):binScootSamps:(size(lfp, 1) - xcorr_window_inds(end));
numbins = length(timestamps);
EMGCorr = zeros(numbins, 1);
nCh = length(xcorr_chs);
nPairs = nCh * (nCh - 1) / 2;
nWin = length(xcorr_window_inds);
% bins per block, keeping the block of windowed data near 256 MB
binsPerBlock = max(1, floor(256e6/(8 * nWin * max(nCh, 1))));
for binindstart = 1:binsPerBlock:numbins
    binindend = min(numbins, binindstart+binsPerBlock-1);
    nb = binindend - binindstart + 1;
    sampleinds = xcorr_window_inds(:) + timestamps(binindstart:binindend);
    z = reshape(double(lfp(sampleinds(:), :)), nWin, nb, nCh);
    z = z - sum(z, 1) / nWin;
    z = z ./ sqrt(sum(z.^2, 1));
    for j = 1:nCh
        for k = (j + 1):nCh
            r = sum(z(:, :, j).*z(:, :, k), 1);
            r(r > 1) = 1; % guard against rounding, as corr does
            r(r < -1) = -1;
            EMGCorr(binindstart:binindend) = EMGCorr(binindstart:binindend) + r(:);
        end
    end
end

EMGCorr = EMGCorr / (length(xcorr_chs) * (length(xcorr_chs) - 1) / 2); % normalize

EMGFromLFP.timestamps = timestamps' ./ Fs;
EMGFromLFP.data = EMGCorr;
EMGFromLFP.channels = xcorr_chs;
EMGFromLFP.detectorName = 'getEMGFromLFP';
EMGFromLFP.samplingFrequency = samplingFrequency;
EMGFromLFP.detectorinfo.channelSelection = selectionMethod;
EMGFromLFP.detectorinfo.maxAutoChannels = maxAutoChannels;
EMGFromLFP.detectorinfo.nUsableSpikeGroups = nUsableSpikeGroups;
EMGFromLFP.detectorinfo.spikeGroupsUsed = spkgrpsused;
EMGFromLFP.detectorinfo.specialChannels = specialChannels;
EMGFromLFP.detectorinfo.rejectChannels = rejectChannels;
EMGFromLFP.detectorinfo.nChannelPairs = nPairs;

if ~any(~isnan(EMGFromLFP.data))
    keyboard
end
if saveMat
    %Save in buzcodeformat
    save(matfilename, 'EMGFromLFP');
end


function [filt_sig, Filt] = filtsig_in(sig, Fs, filtband_or_Filt)
% [filt_sig, Filt] = filtsig(sig, dt_ms, filtband_or_Filt)
%
% Created by: Erik Schomburg, 2011

if isnumeric(filtband_or_Filt)
    h = fdesign.bandpass(filtband_or_Filt(1), filtband_or_Filt(2), filtband_or_Filt(3), filtband_or_Filt(4), ...
        60, 1, 60, Fs);
    Filt = design(h, 'butter', 'MatchExactly', 'passband');
else
    Filt = filtband_or_Filt;
end

if ~isempty(sig)
    if iscell(sig)
        filt_sig = cell(size(sig));
        for i = 1:length(sig(:))
            filt_sig{i} = filter(Filt, sig{i});
            filt_sig{i} = filter(Filt, filt_sig{i}(end:-1:1));
            filt_sig{i} = filt_sig{i}(end:-1:1);
        end
    elseif ((size(sig, 1) > 1) && (size(sig, 2) > 1))
        filt_sig = zeros(size(sig));
        for i = 1:size(filt_sig, 2)
            filt_sig(:, i) = filter(Filt, sig(:, i));
            filt_sig(:, i) = filter(Filt, filt_sig(end:-1:1, i));
            filt_sig(:, i) = filt_sig(end:-1:1, i);
        end
    else
        filt_sig = filter(Filt, sig);
        filt_sig = filter(Filt, filt_sig(end:-1:1));
        filt_sig = filt_sig(end:-1:1);
    end
else
    filt_sig = [];
end
