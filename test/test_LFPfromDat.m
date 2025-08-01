function tests = test_LFPfromDat
% test_LFPfromDat A test suite for the LFPfromDat function.
%   This script uses the MATLAB unit testing framework to verify
%   the functionality of LFPfromDat.m.

tests = functiontests(localfunctions);

end

function testLFPFileCreationAndSize(testCase)
% testLFPFileCreationAndSize Verifies that LFPfromDat creates a file
% with the correct size after downsampling.

%% --- Setup Phase ---
% Create a temporary directory for the test data
testDir = fullfile(pwd, 'LFPfromDat_test_data');
if ~exist(testDir, 'dir')
    mkdir(testDir);
end

% Define dummy session parameters
session.general.name = 'test_session';
session.extracellular.nChannels = 16;
session.extracellular.sr = 30000;    % Input sample rate (Hz)
session.extracellular.srLfp = 1250;  % Output LFP sample rate (Hz)

% Create a dummy .dat file
originalDuration = 10; % seconds of raw data
nSamples = originalDuration * session.extracellular.sr;
datFileName = fullfile(testDir, [session.general.name, '.dat']);
fid = fopen(datFileName, 'w');
fwrite(fid, randi([-500, 500], [session.extracellular.nChannels, nSamples], 'int16'), 'int16');
fclose(fid);

% Save dummy session info, as LFPfromDat needs it
sessionFileName = fullfile(testDir, [session.general.name, '.session.mat']);
save(sessionFileName, 'session');

% Add helper getSession function to the path for the test
addpath(testDir);

% NOTE: LFPfromDat relies on a getSession.m function.
% For this test to be self-contained, we'll create a mock getSession.m
% that just loads the file we just created.
createMockGetSession(testDir, sessionFileName);

%% --- Run the Function ---
% Execute LFPfromDat with the test data
LFPfromDat(testDir);

%% --- Verification Phase ---
lfpFileName = fullfile(testDir, [session.general.name, '.lfp']);
lfpFileInfo = dir(lfpFileName);

% Check 1: Verify that the .lfp file was created
verifyTrue(testCase, isfile(lfpFileName), 'LFP file was not created.');

% Check 2: Verify the size of the created .lfp file
sampleRatio = session.extracellular.sr / session.extracellular.srLfp;
expectedSamples = nSamples / sampleRatio;
expectedSize = expectedSamples * session.extracellular.nChannels * 2; % 2 bytes per int16
% Allow for a small tolerance in file size due to potential header/footer differences
% or minor rounding errors, though it should be exact in this case.
tolerance = 0; % The file size should be exact with this setup
verifyEqual(testCase, lfpFileInfo.bytes, expectedSize, 'AbsTol', tolerance, ...
    'File size mismatch for the generated .lfp file.');

%% --- Teardown Phase ---
% Remove all generated files and the directory
delete(datFileName);
delete(sessionFileName);
delete(lfpFileName);
delete(fullfile(testDir, 'getSession.m')); % Clean up mock function
rmpath(testDir);
rmdir(testDir);
end

function testCustomParametersAndSize(testCase)
% testCustomParametersAndSize Verifies LFPfromDat with custom parameters.

%% --- Setup Phase ---
testDir = fullfile(pwd, 'LFPfromDat_custom_test_data');
if ~exist(testDir, 'dir')
    mkdir(testDir);
end

session.general.name = 'test_custom_session';
session.extracellular.nChannels = 32;
session.extracellular.sr = 40000;
session.extracellular.srLfp = 2000; % Custom LFP sample rate

originalDuration = 5; % seconds
nSamples = originalDuration * session.extracellular.sr;
datFileName = fullfile(testDir, [session.general.name, '.dat']);
fid = fopen(datFileName, 'w');
fwrite(fid, randi([-1000, 1000], [session.extracellular.nChannels, nSamples], 'int16'), 'int16');
fclose(fid);

sessionFileName = fullfile(testDir, [session.general.name, '.session.mat']);
save(sessionFileName, 'session');

addpath(testDir);
createMockGetSession(testDir, sessionFileName);

%% --- Run the Function ---
LFPfromDat(testDir, 'outFs', 2000, 'lopass', 800, 'useGPU', false);

%% --- Verification Phase ---
lfpFileName = fullfile(testDir, [session.general.name, '.lfp']);
lfpFileInfo = dir(lfpFileName);

verifyTrue(testCase, isfile(lfpFileName), 'LFP file was not created with custom parameters.');

sampleRatio = session.extracellular.sr / 2000; % Using custom outFs
expectedSamples = nSamples / sampleRatio;
expectedSize = expectedSamples * session.extracellular.nChannels * 2;

verifyEqual(testCase, lfpFileInfo.bytes, expectedSize, 'AbsTol', 0, ...
    'File size mismatch for the custom parameter test.');

%% --- Teardown Phase ---
delete(datFileName);
delete(sessionFileName);
delete(lfpFileName);
delete(fullfile(testDir, 'getSession.m'));
rmpath(testDir);
rmdir(testDir);
end

function testFileAlreadyExists(testCase)
% testFileAlreadyExists Verifies that LFPfromDat skips processing if the .lfp
% file already exists.

%% --- Setup Phase ---
testDir = fullfile(pwd, 'LFPfromDat_exist_test_data');
if ~exist(testDir, 'dir')
    mkdir(testDir);
end
session.general.name = 'test_exist_session';
session.extracellular.nChannels = 1;
session.extracellular.sr = 30000;
session.extracellular.srLfp = 1000;
datFileName = fullfile(testDir, [session.general.name, '.dat']);
fid = fopen(datFileName, 'w');
fwrite(fid, randi([-1, 1], [1, 30000], 'int16'), 'int16'); % 1 second of data
fclose(fid);
sessionFileName = fullfile(testDir, [session.general.name, '.session.mat']);
save(sessionFileName, 'session');
lfpFileName = fullfile(testDir, [session.general.name, '.lfp']);
% Create a dummy, small LFP file
fid = fopen(lfpFileName, 'w');
fwrite(fid, zeros(1, 100, 'int16'), 'int16');
fclose(fid);

initialSize = dir(lfpFileName).bytes;

addpath(testDir);
createMockGetSession(testDir, sessionFileName);

%% --- Run the Function ---
% This should not overwrite the file or change its size
LFPfromDat(testDir);

%% --- Verification Phase ---
finalSize = dir(lfpFileName).bytes;

% Check that the file was not overwritten
verifyEqual(testCase, finalSize, initialSize, 'The existing LFP file was unexpectedly modified.');

%% --- Teardown Phase ---
delete(datFileName);
delete(sessionFileName);
delete(lfpFileName);
delete(fullfile(testDir, 'getSession.m'));
rmpath(testDir);
rmdir(testDir);
end


function createMockGetSession(testDir, sessionFileName)
% Helper function to create a mock getSession.m file for testing
fid = fopen(fullfile(testDir, 'getSession.m'), 'w');
fprintf(fid, 'function session = getSession(varargin)\n');
fprintf(fid, '    p = inputParser;\n');
fprintf(fid, '    addParameter(p, ''basepath'', [], @isfolder);\n');
fprintf(fid, '    parse(p, varargin{:});\n');
fprintf(fid, '    basepath = p.Results.basepath;\n');
fprintf(fid, '    sessionFile = dir(fullfile(basepath, ''*.session.mat''));\n');
fprintf(fid, '    if ~isempty(sessionFile)\n');
fprintf(fid, '        load(fullfile(basepath, sessionFile(1).name), ''session'');\n');
fprintf(fid, '    else\n');
fprintf(fid, '        error(''Mock session file not found'');\n');
fprintf(fid, '    end\n');
fprintf(fid, 'end\n');
fclose(fid);
end