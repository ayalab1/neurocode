function tests = test_read_Intan_RHD2000_file
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
root = fileparts(fileparts(mfilename('fullpath')));
testCase.TestData.originalPath = path;
testCase.TestData.folder = tempname;
mkdir(testCase.TestData.folder);
% Exercise the public reader without displaying its file picker.
mock = fullfile(testCase.TestData.folder, 'uigetfile.m');
fid = fopen(mock, 'w');
fprintf(fid, ['function [file, folder, index] = uigetfile(varargin)\n' ...
    'global NEUROCODE_INTAN_TEST_FILE\n' ...
    'if isempty(NEUROCODE_INTAN_TEST_FILE)\n' ...
    'file = 0; folder = 0; index = 0; return;\nend\n' ...
    '[folder, name, extension] = fileparts(NEUROCODE_INTAN_TEST_FILE);\n' ...
    'file = [name extension]; folder = [folder filesep]; index = 1;\nend\n']);
fclose(fid);
addpath(fullfile(root, 'preProcessing', 'intan'));
addpath(testCase.TestData.folder, '-begin');
end

function teardownOnce(testCase)
path(testCase.TestData.originalPath);
clear uigetfile;
clear global NEUROCODE_INTAN_TEST_FILE;
rmdir(testCase.TestData.folder, 's');
end

function setup(testCase)
global NEUROCODE_INTAN_TEST_FILE;
testCase.TestData.output = tempname(testCase.TestData.folder);
mkdir(testCase.TestData.output);
NEUROCODE_INTAN_TEST_FILE = fullfile(testCase.TestData.output, 'test.rhd');
end

function testRawExportsAcrossVersions(testCase)
for version = {[1 0], [1 1], [1 2], [1 3], [2 0], [3 0]}
    data = writeFixture(testCase, version{1}, 0, false);
    read_Intan_RHD2000_file('streamToDat', true);
    verifyEqual(testCase, readDat(testCase, 'amplifier', 'int16'), ...
        reshape(int16(double(data.amplifier) - 32768), [], 1));
    verifyEqual(testCase, readDat(testCase, 'auxiliary', 'uint16'), ...
        reshape(repelem(data.auxiliary, 1, 4), [], 1));
    verifyEqual(testCase, readDat(testCase, 'analogin', 'uint16'), data.analog(:));
    verifyEqual(testCase, readDat(testCase, 'digitalin', 'uint16'), data.digital(:));
    verifyEqual(testCase, readDat(testCase, 'time', 'int32'), data.time(:));
end
end

function testNotchMatchesInMemoryAcrossBlocks(testCase)
for version = {[1 3], [2 0], [3 0]}
    for mode = [1 2]
        writeFixture(testCase, version{1}, mode, false);
        read_Intan_RHD2000_file;
        expected = evalin('base', 'int16(amplifier_data(:))');
        read_Intan_RHD2000_file('streamToDat', true);
        verifyEqual(testCase, readDat(testCase, 'amplifier', 'int16'), expected);
    end
end
end

function testNoSignalChannels(testCase)
data = writeFixture(testCase, [2 0], 0, true);
read_Intan_RHD2000_file('streamToDat', true);
for name = {'amplifier', 'auxiliary', 'analogin', 'digitalin'}
    verifyEmpty(testCase, readDat(testCase, name{1}, 'uint16'));
end
verifyEqual(testCase, readDat(testCase, 'time', 'int32'), data.time(:));
end

function testIncompleteBlockRejectedBeforeExport(testCase)
writeFixture(testCase, [2 0], 0, false);
global NEUROCODE_INTAN_TEST_FILE;
fid = fopen(NEUROCODE_INTAN_TEST_FILE, 'a');
fwrite(fid, 1, 'uint8');
fclose(fid);
opened = openFileIds;
verifyError(testCase, @() read_Intan_RHD2000_file('streamToDat', true), ...
    'Intan:IncompleteDataBlock');
verifyEqual(testCase, openFileIds, opened);
verifyEmpty(testCase, dir(fullfile(testCase.TestData.output, '*.dat')));
end

function testOutputOpenFailureClosesFiles(testCase)
writeFixture(testCase, [2 0], 0, false);
mkdir(fullfile(testCase.TestData.output, 'analogin.dat'));
opened = openFileIds;
verifyError(testCase, @() read_Intan_RHD2000_file('streamToDat', true), ...
    'Intan:OutputOpenFailed');
verifyEqual(testCase, openFileIds, opened);
end

function testTimestampGapAcrossBlocks(testCase)
data = writeFixture(testCase, [2 0], 0, false);
global NEUROCODE_INTAN_TEST_FILE;
fid = fopen(NEUROCODE_INTAN_TEST_FILE, 'r+', 'ieee-le');
fseek(fid, data.headerBytes + data.blockBytes, 'bof');
fwrite(fid, 1000, 'int32');
fclose(fid);
output = evalc('read_Intan_RHD2000_file(''streamToDat'', true)');
verifySubstring(testCase, output, '2 gaps in timestamp data found');
end

function testUnsignedTimestampBits(testCase)
data = writeFixture(testCase, [1 0], 0, false);
global NEUROCODE_INTAN_TEST_FILE;
fid = fopen(NEUROCODE_INTAN_TEST_FILE, 'r+', 'ieee-le');
fseek(fid, data.headerBytes, 'bof');
fwrite(fid, uint32(2^31 + (0:59)), 'uint32');
fclose(fid);
read_Intan_RHD2000_file('streamToDat', true);
ticks = readDat(testCase, 'time', 'uint32');
verifyEqual(testCase, ticks(1:60), uint32(2^31 + (0:59)'));
end

function testDefaultStillPopulatesWorkspace(testCase)
data = writeFixture(testCase, [2 0], 0, false);
read_Intan_RHD2000_file;
verifyEqual(testCase, evalin('base', 'amplifier_data'), double(data.amplifier) - 32768);
verifyEqual(testCase, evalin('base', 'aux_input_data'), 37.4e-6 * double(data.auxiliary));
verifyEqual(testCase, evalin('base', 't_amplifier'), double(data.time) / 1000);
verifyEmpty(testCase, dir(fullfile(testCase.TestData.output, '*.dat')));
end

function testHeaderOnlyDoesNotExport(testCase)
data = writeFixture(testCase, [2 0], 0, false);
global NEUROCODE_INTAN_TEST_FILE;
fid = fopen(NEUROCODE_INTAN_TEST_FILE, 'r');
header = fread(fid, data.headerBytes, 'uint8=>uint8');
fclose(fid);
fid = fopen(NEUROCODE_INTAN_TEST_FILE, 'w');
fwrite(fid, header, 'uint8');
fclose(fid);
read_Intan_RHD2000_file('streamToDat', true);
verifyEmpty(testCase, dir(fullfile(testCase.TestData.output, '*.dat')));
end

function testCancelledPicker(testCase)
global NEUROCODE_INTAN_TEST_FILE;
NEUROCODE_INTAN_TEST_FILE = '';
read_Intan_RHD2000_file('streamToDat', true);
verifyEmpty(testCase, dir(fullfile(testCase.TestData.output, '*.dat')));
end

function data = writeFixture(testCase, version, notchMode, noChannels)
global NEUROCODE_INTAN_TEST_FILE;
fid = fopen(NEUROCODE_INTAN_TEST_FILE, 'w', 'ieee-le');
cleanup = onCleanup(@() fclose(fid));
fwrite(fid, hex2dec('c6912702'), 'uint32');
fwrite(fid, version, 'int16');
fwrite(fid, 1000, 'single');
fwrite(fid, 1, 'int16');
fwrite(fid, [1 1 400 1 1 400], 'single');
fwrite(fid, notchMode, 'int16');
fwrite(fid, [1000 1000], 'single');
for i = 1:3
    writeString(fid, '');
end
hasTemp = version(1) > 1 || version(2) >= 1;
if hasTemp
    fwrite(fid, ~noChannels, 'int16');
end
if version(1) > 1 || version(2) >= 3
    fwrite(fid, 13, 'int16');
end
if version(1) > 1
    writeString(fid, '');
end
types = [0 0 1 1 2 3 3 4 4 5];
orders = [0 1 0 1 0 0 1 1 15 0];
if noChannels
    fwrite(fid, 0, 'int16');
else
    fwrite(fid, 1, 'int16');
    writeString(fid, 'Port A');
    writeString(fid, 'A');
    fwrite(fid, [1 numel(types) 2], 'int16');
    for i = 1:numel(types)
        writeString(fid, sprintf('A-%03d', i));
        writeString(fid, sprintf('A-%03d', i));
        fwrite(fid, [orders(i) i-1 types(i) 1 0 0 0 0 0 0], 'int16');
        fwrite(fid, [1000 0], 'single');
    end
end
data.headerBytes = ftell(fid);
samples = 128;
if version(1) == 1
    samples = 60;
end
data.amplifier = uint16(reshape(mod(0:6*samples-1, 65536), 2, []));
data.amplifier(:, 1:2) = uint16([0 32768; 65535 32767]);
data.auxiliary = uint16(reshape(1000:1000+6*samples/4-1, 2, []));
data.analog = uint16(reshape(30000:30000+6*samples-1, 2, []));
data.digital = uint16(mod(0:3*samples-1, 2)*2 + mod(floor((0:3*samples-1)/3), 2)*32768);
data.time = int32(-10:3*samples-11);
if version(1) == 1 && version(2) < 2
    data.time = int32(0:3*samples-1);
end
for block = 1:3
    start = ftell(fid);
    indices = (block-1)*samples + (1:samples);
    fwrite(fid, data.time(indices), 'int32');
    if ~noChannels
        fwrite(fid, data.amplifier(:, indices)', 'uint16');
        auxIndices = (block-1)*samples/4 + (1:samples/4);
        fwrite(fid, data.auxiliary(:, auxIndices)', 'uint16');
        fwrite(fid, 50000, 'uint16'); % supply voltage
        if hasTemp
            fwrite(fid, -100, 'int16');
        end
        fwrite(fid, data.analog(:, indices)', 'uint16');
        fwrite(fid, data.digital(indices), 'uint16');
        fwrite(fid, 1234*ones(1, samples), 'uint16'); % digital output
    end
    data.blockBytes = ftell(fid) - start;
end
end

function writeString(fid, value)
fwrite(fid, 2*numel(value), 'uint32');
fwrite(fid, uint16(value), 'uint16');
end

function ids = openFileIds
if exist('openedFiles', 'builtin') || exist('openedFiles', 'file')
    ids = openedFiles;
else
    ids = fopen('all');
end
end

function data = readDat(testCase, name, precision)
fid = fopen(fullfile(testCase.TestData.output, [name '.dat']), 'r', 'ieee-le');
cleanup = onCleanup(@() fclose(fid));
data = fread(fid, Inf, [precision '=>' precision]);
end
