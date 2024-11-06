function tests = test_cut_dat
tests = functiontests(localfunctions);
end

function testFileRenamingAndProcessing(testCase)
% Setup: Create mock .dat files to simulate ephys data files
mockFiles = {'amplifier.dat', 'auxiliary.dat', 'digitalin.dat'};
nChannels = 32; % Updated channel count
sampleRate = 20000; % Hz
originalDuration = 2; % seconds
cutDuration = 1; % seconds (to test trimming)
nSamples = originalDuration * sampleRate * nChannels; % Total samples for 2 seconds

for i = 1:length(mockFiles)
    fid = fopen(mockFiles{i}, 'w');
    fwrite(fid, ones(1, nSamples, 'int16'), 'int16'); % Write 2 seconds of sample data
    fclose(fid);
end

% Run the function with maxTime set to 1 second
cut_dat(nChannels, cutDuration, 'sample_rate', sampleRate);

% Verify: Check that _original.dat files were created and cut correctly
for i = 1:length(mockFiles)
    originalFile = [mockFiles{i}(1:end - 4), '_original.dat'];
    newFile = mockFiles{i};

    % Check existence of original and trimmed files
    verifyTrue(testCase, isfile(originalFile), ['Missing file: ', originalFile]);
    verifyTrue(testCase, isfile(newFile), ['Processed file not found: ', newFile]);

    % Check that the new file has data corresponding to 1 second
    newFileInfo = dir(newFile);
    expectedSize = cutDuration * sampleRate * nChannels * 2; % 1 sec of int16 data with 32 channels
    % Set a tolerance of 64 bytes
    tolerance = 64;
    verifyTrue(testCase, abs(newFileInfo.bytes-expectedSize) <= tolerance, ...
        sprintf('File size mismatch for: %s (within tolerance of %d bytes)', newFile, tolerance));
end

% Cleanup: Remove generated files
for i = 1:length(mockFiles)
    delete(mockFiles{i});
    delete([mockFiles{i}(1:end - 4), '_original.dat']);
end
end
