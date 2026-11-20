function tests = test_getEMGFromLFP
% Create a test suite for the getEMGFromLFP function
tests = functiontests(localfunctions);
end

function setup(testCase)
% Setup function to create a temporary directory and mock files
testCase.TestData.tempDir = tempname;
mkdir(testCase.TestData.tempDir);

% Derive the recording name from the basepath
[~, recordingname] = fileparts(testCase.TestData.tempDir);

% Create a mock session.mat file with the correct name
session.extracellular.nChannels = 32;
session.extracellular.spikeGroups.channels = {1:8, 9:16, 17:24, 25:32};
session.extracellular.srLfp = 1250; % Sampling rate for LFP
session.extracellular.sr = 20000; % Sampling rate for dat file

% Save the session.mat file with the correct name
save(fullfile(testCase.TestData.tempDir, [recordingname, '.session.mat']), 'session');

% Create a mock .lfp file
nChannels = session.extracellular.nChannels;
Fs = session.extracellular.srLfp;
duration = 10; % 10 seconds of data
data = randn(duration * Fs, nChannels); % Random LFP data
fid = fopen(fullfile(testCase.TestData.tempDir, [recordingname, '.lfp']), 'w');
fwrite(fid, data', 'int16');
fclose(fid);
end

function test_basic_functionality(testCase)
% Test basic functionality of getEMGFromLFP
basepath = testCase.TestData.tempDir;

% Call the function
EMGFromLFP = getEMGFromLFP(basepath, 'restrict', [0, 5], 'samplingFrequency', 2);

% Verify the output structure
verifyTrue(testCase, isfield(EMGFromLFP, 'timestamps'));
verifyTrue(testCase, isfield(EMGFromLFP, 'data'));
verifyTrue(testCase, isfield(EMGFromLFP, 'channels'));
verifyTrue(testCase, isfield(EMGFromLFP, 'detectorName'));
verifyTrue(testCase, isfield(EMGFromLFP, 'samplingFrequency'));

% Verify the dimensions of the output
verifyEqual(testCase, length(EMGFromLFP.timestamps), length(EMGFromLFP.data));
verifyEqual(testCase, EMGFromLFP.samplingFrequency, 2);

% Verify that the data is not all NaNs
verifyFalse(testCase, all(isnan(EMGFromLFP.data)));
end

function test_with_special_channels(testCase)
% Test the function with special channels specified
basepath = testCase.TestData.tempDir;

% Call the function with special channels
specialChannels = [1, 2, 3];
EMGFromLFP = getEMGFromLFP(basepath, 'restrict', [0, 5], ...
    'specialChannels', specialChannels, 'samplingFrequency', 2);

% Verify that the special channels are included in the output
verifyTrue(testCase, all(ismember(specialChannels, EMGFromLFP.channels)));
end

function test_with_reject_channels(testCase)
% Test the function with reject channels specified
basepath = testCase.TestData.tempDir;

% Call the function with reject channels
rejectChannels = [1, 2, 3];
EMGFromLFP = getEMGFromLFP(basepath, 'restrict', [0, 5], ...
    'rejectChannels', rejectChannels, 'samplingFrequency', 2);

% Verify that the reject channels are not included in the output
verifyFalse(testCase, any(ismember(rejectChannels, EMGFromLFP.channels)));
end

function test_with_restrict_channels(testCase)
% Test the function with restrict channels specified
basepath = testCase.TestData.tempDir;

% Call the function with restrict channels
restrictChannels = [1, 2, 3];
EMGFromLFP = getEMGFromLFP(basepath, 'restrict', [0, 5], ...
    'restrictChannels', restrictChannels, 'samplingFrequency', 2);

% Verify that only the restrict channels are included in the output
verifyEqual(testCase, EMGFromLFP.channels, restrictChannels);
end

function test_with_fromDat_option(testCase)
% Test the function with the fromDat option
basepath = testCase.TestData.tempDir;

% Create a mock .dat file
nChannels = 32;
Fs = 20000; % Sampling rate for dat file
duration = 10; % 10 seconds of data
data = randn(duration * Fs, nChannels); % Random data
fid = fopen(fullfile(testCase.TestData.tempDir, 'test.dat'), 'w');
fwrite(fid, data', 'int16');
fclose(fid);

% Call the function with fromDat option
EMGFromLFP = getEMGFromLFP(basepath, 'restrict', [0, 5], ...
    'fromDat', true, 'samplingFrequency', 2);

% Verify the output structure
verifyTrue(testCase, isfield(EMGFromLFP, 'timestamps'));
verifyTrue(testCase, isfield(EMGFromLFP, 'data'));
verifyTrue(testCase, isfield(EMGFromLFP, 'channels'));
verifyTrue(testCase, isfield(EMGFromLFP, 'detectorName'));
verifyTrue(testCase, isfield(EMGFromLFP, 'samplingFrequency'));

% Verify the dimensions of the output
verifyEqual(testCase, length(EMGFromLFP.timestamps), length(EMGFromLFP.data));
verifyEqual(testCase, EMGFromLFP.samplingFrequency, 2);

% Verify that the data is not all NaNs
verifyFalse(testCase, all(isnan(EMGFromLFP.data)));
end

function test_few_groups_selection_and_values_match_legacy(testCase)
% Layouts with <= maxAutoChannels spike groups must give the same
% channels as the former implementation for the same random seed, and the
% same EMG values up to floating-point rounding.
basepath = testCase.TestData.tempDir;
cases = {{}, {'rejectChannels', [1, 2, 9]}, {'specialChannels', [5, 30]}};
for c = 1:numel(cases)
    opts = cases{c};
    rng(7 + c);
    newEMG = getEMGFromLFP(basepath, 'restrict', [0, 8], 'saveMat', false, ...
        'overwrite', true, 'noPrompts', true, opts{:});
    rng(7 + c);
    [legacyChannels, legacyData] = legacyEMGFromLFP(basepath, [0, 8], opts{:});
    verifyEqual(testCase, newEMG.channels, legacyChannels);
    verifyEqual(testCase, newEMG.data, legacyData, 'AbsTol', 1e-12);
    verifyEqual(testCase, newEMG.detectorinfo.channelSelection, 'onePerSpikeGroup');
end
end

function test_single_shank_selection_matches_legacy(testCase)
basepath = makeSession(testCase, {1:16}, 16);
rng(3);
newEMG = getEMGFromLFP(basepath, 'restrict', [0, 8], 'saveMat', false, ...
    'overwrite', true, 'noPrompts', true);
rng(3);
[legacyChannels, legacyData] = legacyEMGFromLFP(basepath, [0, 8]);
verifyEqual(testCase, newEMG.channels, legacyChannels);
verifyEqual(testCase, numel(newEMG.channels), 5);
verifyEqual(testCase, newEMG.data, legacyData, 'AbsTol', 1e-12);
verifyEqual(testCase, newEMG.detectorinfo.channelSelection, 'fiveFromSingleSpikeGroup');
end

function test_many_spike_groups_are_capped(testCase)
% Many small spike groups: 42 groups of 6 channels.
nGroups = 42;
groups = arrayfun(@(g) (g - 1)*6+(1:6), 1:nGroups, 'UniformOutput', false);
basepath = makeSession(testCase, groups, nGroups*6);

EMG = getEMGFromLFP(basepath, 'restrict', [0, 4], 'saveMat', false, ...
    'overwrite', true, 'noPrompts', true);
verifyEqual(testCase, numel(EMG.channels), 16);
verifyEqual(testCase, EMG.detectorinfo.channelSelection, 'onePerSpikeGroupCapped');
verifyEqual(testCase, EMG.detectorinfo.nUsableSpikeGroups, nGroups);
expectedGroups = unique(round(linspace(1, nGroups, 16)));
verifyEqual(testCase, EMG.detectorinfo.spikeGroupsUsed, expectedGroups);
groupOfChannel = ceil(EMG.channels/6);
verifyEqual(testCase, sort(groupOfChannel), expectedGroups);
verifyEqual(testCase, EMG.detectorinfo.nChannelPairs, 16*15/2);

% special channels are added on top of the cap
EMG = getEMGFromLFP(basepath, 'restrict', [0, 4], 'saveMat', false, ...
    'overwrite', true, 'noPrompts', true, 'specialChannels', [2, 3]);
verifyTrue(testCase, all(ismember([2, 3], EMG.channels)));
verifyLessThanOrEqual(testCase, numel(EMG.channels), 18);

% the cap can be raised or disabled
EMG = getEMGFromLFP(basepath, 'restrict', [0, 4], 'saveMat', false, ...
    'overwrite', true, 'noPrompts', true, 'maxAutoChannels', inf);
verifyEqual(testCase, numel(EMG.channels), nGroups);
verifyEqual(testCase, EMG.detectorinfo.channelSelection, 'onePerSpikeGroup');
end

% ---------------------------------------------------------------------
function basepath = makeSession(testCase, groups, nChannels)
basepath = tempname;
mkdir(basepath);
testCase.addTeardown(@() rmdir(basepath, 's'));
[~, recordingname] = fileparts(basepath);
session.extracellular.nChannels = nChannels;
session.extracellular.spikeGroups.channels = groups;
session.extracellular.srLfp = 1250;
session.extracellular.sr = 20000;
save(fullfile(basepath, [recordingname, '.session.mat']), 'session');
data = round(1000*randn(10*1250, nChannels));
fid = fopen(fullfile(basepath, [recordingname, '.lfp']), 'w');
fwrite(fid, data', 'int16');
fclose(fid);
end

function [xcorr_chs, EMGCorr] = legacyEMGFromLFP(basepath, restrict, varargin)
% Frozen copy of the channel selection, filtering and per-bin corr() loop
% of getEMGFromLFP before the many-spike-group change (neurocode 732cbbd).
p = inputParser;
addParameter(p, 'specialChannels', []);
addParameter(p, 'rejectChannels', []);
parse(p, varargin{:});
specialChannels = p.Results.specialChannels;
rejectChannels = p.Results.rejectChannels;

[~, recordingname] = fileparts(basepath);
load(fullfile(basepath, [recordingname, '.session.mat']), 'session');
nChannels = session.extracellular.nChannels;
SpkGrps = session.extracellular.spikeGroups.channels;
Fs = session.extracellular.srLfp;
binScootS = 1 ./ 2;
binScootSamps = round(Fs*binScootS);
corrChunkSz = 20;

usablechannels = [];
spkgrpstouse = [];
for gidx = 1:length(SpkGrps)
    usableshankchannels{gidx} = setdiff(SpkGrps{gidx}, rejectChannels); %#ok<AGROW>
    usablechannels = cat(2, usablechannels, usableshankchannels{gidx});
    if ~isempty(usableshankchannels{gidx})
        spkgrpstouse = cat(2, spkgrpstouse, gidx);
    end
end
usableshankchannels = usableshankchannels(spkgrpstouse);
if length(usableshankchannels) > 1
    n = 1;
else
    n = 5;
end
xcorr_chs = [];
for gidx = 1:length(usableshankchannels)
    if ~isempty(usableshankchannels{gidx})
        if n == 1
            randChfromShank = usableshankchannels{gidx}(randi(length(usableshankchannels{gidx}), n, 1));
        elseif n > 1
            randChfromShank = usableshankchannels{gidx}(randperm(length(usableshankchannels{gidx}), n));
        end
        xcorr_chs = [xcorr_chs, randChfromShank]; %#ok<AGROW>
    end
end
xcorr_chs = unique([xcorr_chs, specialChannels]);

lfp = loadBinary(fullfile(basepath, [recordingname, '.lfp']), 'nChannels', nChannels, ...
    'channels', xcorr_chs, 'start', restrict(1), 'duration', diff(restrict), 'frequency', Fs);
maxfreqband = floor(max([625, Fs / 2]));
h = fdesign.bandpass(275, 300, maxfreqband-25, maxfreqband, 60, 1, 60, Fs);
Filt = design(h, 'butter', 'MatchExactly', 'passband');
filt_sig = zeros(size(lfp));
for i = 1:size(filt_sig, 2)
    filt_sig(:, i) = filter(Filt, lfp(:, i));
    filt_sig(:, i) = filter(Filt, filt_sig(end:-1:1, i));
    filt_sig(:, i) = filt_sig(end:-1:1, i);
end
lfp = filt_sig;

xcorr_window_samps = round(binScootS*Fs);
xcorr_window_inds = -xcorr_window_samps:xcorr_window_samps;
timestamps = (1 + xcorr_window_inds(end)):binScootSamps:(size(lfp, 1) - xcorr_window_inds(end));
EMGCorr = zeros(length(timestamps), 1);
for j = 1:(length(xcorr_chs))
    for k = (j + 1):length(xcorr_chs)
        c1 = [];
        c2 = [];
        binind = 0;
        binindstart = 1;
        for i = timestamps
            binind = binind + 1;
            c1 = cat(2, c1, lfp(i + xcorr_window_inds, j));
            c2 = cat(2, c2, lfp(i + xcorr_window_inds, k));
            if size(c1, 2) == corrChunkSz || i == timestamps(end)
                tmp = diag(corr(c1, c2));
                EMGCorr(binindstart:binind) = EMGCorr(binindstart:binind) + tmp;
                c1 = [];
                c2 = [];
                binindstart = binind + 1;
            end
        end
    end
end
EMGCorr = EMGCorr / (length(xcorr_chs) * (length(xcorr_chs) - 1) / 2);
end
