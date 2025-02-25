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
EMGFromLFP = getEMGFromLFP(basepath, 'restrict', [0, 5],...
    'specialChannels', specialChannels, 'samplingFrequency', 2);

% Verify that the special channels are included in the output
verifyTrue(testCase, all(ismember(specialChannels, EMGFromLFP.channels)));
end

function test_with_reject_channels(testCase)
% Test the function with reject channels specified
basepath = testCase.TestData.tempDir;

% Call the function with reject channels
rejectChannels = [1, 2, 3];
EMGFromLFP = getEMGFromLFP(basepath, 'restrict', [0, 5],...
    'rejectChannels', rejectChannels, 'samplingFrequency', 2);

% Verify that the reject channels are not included in the output
verifyFalse(testCase, any(ismember(rejectChannels, EMGFromLFP.channels)));
end

function test_with_restrict_channels(testCase)
% Test the function with restrict channels specified
basepath = testCase.TestData.tempDir;

% Call the function with restrict channels
restrictChannels = [1, 2, 3];
EMGFromLFP = getEMGFromLFP(basepath, 'restrict', [0, 5],...
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
EMGFromLFP = getEMGFromLFP(basepath, 'restrict', [0, 5],...
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