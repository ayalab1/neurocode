function LFPfromDat(basepath, varargin)
%LFPfromDat Extract Local Field Potential (LFP) from wideband .dat file
%   This function performs extraction of LFP signals from raw wideband
%   neural data (.dat file) through filtering and downsampling.
%
%   Syntax:
%   LFPfromDat(basepath)
%   LFPfromDat(basepath, 'ParameterName', ParameterValue, ...)
%
%   Inputs:
%   basepath - Path to directory containing the .dat file (default: current directory)
%
%   Optional Name-Value Parameters:
%   'datFile'     - Name of .dat file (default: [basename '.dat'])
%   'outFs'       - Output sampling rate (default: session.extracellular.srLfp)
%   'lopass'      - Low-pass cutoff frequency (default: 450 Hz)
%   'useGPU'      - Enable GPU acceleration (default: false)
%   'inFs'        - Input sampling rate (default: session.extracellular.sr)
%   'localDir'    - Temporary directory for processing (default: basepath)
%   'session'     - Session metadata structure (default: loaded from basepath)
%
%   Outputs:
%   Creates a .lfp file in the basepath directory with downsampled, filtered data
%
%   Processing Pipeline:
%   1. Reads raw data in chunks to minimize memory usage
%   2. Applies optimized sinc filter (low-pass)
%   3. Downsamples to target rate
%   4. Writes results to binary file
%
%   Example:
%   % Basic usage
%   LFPfromDat('/path/to/data')
%
%   % Custom parameters
%   LFPfromDat('/path/to/data', 'lopass', 300, 'outFs', 1250, 'useGPU', true)
%
%   Dependencies:
%   - iosr.dsp package
%   - Parallel Computing Toolbox (for GPU support)
%   - getSession.m (from CellExplorer)
%
%   See also: fast_sinc_filter_matrix, getSession
%
%   Note: Will skip processing if .lfp file already exists

%% Import
import iosr.dsp.*

%% Input handling
if ~exist('basepath', 'var')
    basepath = pwd;
end

p = inputParser;
addParameter(p, 'datFile', [], @ischar);
addParameter(p, 'outFs', [], @isnumeric);
addParameter(p, 'lopass', 450, @isnumeric);
addParameter(p, 'useGPU', false, @islogical);
addParameter(p, 'inFs', [], @isnumeric);
addParameter(p, 'localDir', [], @isfolder);
addParameter(p, 'session', [], @isstruct);

parse(p, varargin{:})
datFile = p.Results.datFile;
outFs = p.Results.outFs;
lopass = p.Results.lopass;
useGPU = p.Results.useGPU;
inFs = p.Results.inFs;
localDir = p.Results.localDir;
session = p.Results.session;

if isempty(session)
    session = getSession('basepath', basepath);
end
basename = session.general.name;

if isempty(datFile)
    datFile = [basename, '.dat'];
elseif ~strcmp(datFile(end-3:end), '.dat')
    datFile = [datFile, '.dat'];
end

% GPU setup - simplified
if useGPU && gpuDeviceCount < 1
    warning('No GPU device found, continuing without..')
    useGPU = false;
end

if useGPU
    g = gpuDevice(1);
end

sizeInBytes = 2;

%% Housekeeping
fInfo = checkFile('basepath', basepath, 'filename', datFile, 'searchSubdirs', false);
fdat = fInfo.name;

if isempty(inFs)
    inFs = session.extracellular.sr;
end

nbChan = session.extracellular.nChannels;

if isempty(outFs)
    outFs = session.extracellular.srLfp;
end

if lopass > outFs / 2
    warning('low pass cutoff beyond Nyquist')
end

ratio = lopass / (inFs / 2);
sampleRatio = (inFs / outFs);

% Output file
if ~isempty(localDir)
    flfp = fullfile(localDir, [basename, '.lfp']);
else
    flfp = fullfile(basepath, [basename, '.lfp']);
end

%% chunking
chunksize = 1e5; % Start with original size
if mod(chunksize, sampleRatio) ~= 0
    chunksize = chunksize + sampleRatio - mod(chunksize, sampleRatio);
end

% Buffer size
ntbuff = 525;
if mod(ntbuff, sampleRatio) ~= 0
    ntbuff = ntbuff + sampleRatio - mod(ntbuff, sampleRatio);
end

nBytes = fInfo.bytes;
nbChunks = floor(nBytes/(nbChan * sizeInBytes * chunksize)) - 1;

%% Check if file exists
if exist([basepath, filesep, basename, '.lfp'], 'file') || ...
        exist([basepath, filesep, basename, '.eeg'], 'file')
    fprintf('LFP file already exists \n')
    return
end

%% Main processing
fidI = fopen(fdat, 'r');
fidout = fopen(flfp, 'a');

fprintf('Starting LFP extraction...\n')
tic;

% Pre-allocate output array once
DATA = zeros(nbChan, chunksize/sampleRatio, 'int16');

for ibatch = 1:nbChunks
    if mod(ibatch, 10) == 0
        if ibatch ~= 10
            fprintf(repmat('\b', [1, length([num2str(round(100*(ibatch - 10)/nbChunks)), ' percent complete'])]))
        end
        fprintf('%d percent complete', round(100*ibatch/nbChunks));
    end

    % Read data
    if ibatch > 1
        fseek(fidI, ((ibatch - 1) * (nbChan * sizeInBytes * chunksize))-(nbChan * sizeInBytes * ntbuff), 'bof');
        dat = fread(fidI, nbChan*(chunksize + 2 * ntbuff), 'int16');
        try
            dat = reshape(dat, [nbChan, (chunksize + 2 * ntbuff)]);
        catch
            fidI = fopen(fdat, 'r');
            fseek(fidI, ((ibatch - 1) * (nbChan * sizeInBytes * chunksize))-(nbChan * sizeInBytes * ntbuff), 'bof');
            dat = fread(fidI, nbChan*(chunksize + 2 * ntbuff), 'int16');
            dat = reshape(dat, [nbChan, (chunksize + 2 * ntbuff)]);
        end
    else
        dat = fread(fidI, nbChan*(chunksize + ntbuff), 'int16');
        dat = reshape(dat, [nbChan, (chunksize + ntbuff)]);
    end

    % vectorized filtering
    dat = double(dat);

    if useGPU
        dat = gpuArray(dat);
        filtered_data = fast_sinc_filter_matrix(dat, ratio, useGPU);

        % Downsample
        if ibatch == 1
            DATA = gather(int16(real(filtered_data(:, sampleRatio:sampleRatio:end-ntbuff))));
        else
            DATA = gather(int16(real(filtered_data(:, ntbuff+sampleRatio:sampleRatio:end-ntbuff))));
        end
    else
        filtered_data = fast_sinc_filter_matrix(dat, ratio, useGPU);

        % Downsample
        if ibatch == 1
            DATA = int16(real(filtered_data(:, sampleRatio:sampleRatio:end-ntbuff)));
        else
            DATA = int16(real(filtered_data(:, ntbuff+sampleRatio:sampleRatio:end-ntbuff)));
        end
    end

    % Write to file
    fwrite(fidout, DATA(:), 'int16');
end

% Handle remainder
remainder = nBytes / (sizeInBytes * nbChan) - nbChunks * chunksize;
if remainder > 0
    fseek(fidI, ((nbChunks) * (nbChan * sizeInBytes * chunksize))-(nbChan * sizeInBytes * ntbuff), 'bof');
    dat = fread(fidI, nbChan*(remainder + ntbuff), 'int16');
    dat = reshape(dat, [nbChan, (remainder + ntbuff)]);

    dat = double(dat);

    if useGPU
        dat = gpuArray(dat);
        filtered_data = fast_sinc_filter_matrix(dat, ratio, useGPU);
        DATA = gather(int16(real(filtered_data(:, ntbuff+sampleRatio:sampleRatio:end))));
    else
        filtered_data = fast_sinc_filter_matrix(dat, ratio, useGPU);
        DATA = int16(real(filtered_data(:, ntbuff+sampleRatio:sampleRatio:end)));
    end

    fwrite(fidout, DATA(:), 'int16');
end

fclose(fidI);
fclose(fidout);

if useGPU
    reset(g);
    gpuDevice([]);
end

elapsed = toc;
fprintf('\nLFP extraction completed in %.2f seconds\n', elapsed);

if ~isempty(localDir)
    movefile(flfp, fullfile(basepath, [basename, '.lfp']));
end

disp('LFP file created')
end

function y = fast_sinc_filter_matrix(x, ratio, useGPU)
%FAST_SINC_FILTER_MATRIX Optimized matrix-based sinc filtering
% This version uses vectorized operations but avoids FFT overhead for moderate sizes

[nChan, nSamples] = size(x);

% Create sinc kernel
N = 525; % Use original kernel size
n = -floor(N/2):floor(N/2);
B = sinc(ratio*n) .* ratio;

if useGPU
    B = gpuArray(B);
end

% Use matrix-based convolution for better performance
if nSamples < 50000
    % For smaller chunks, use direct convolution (faster)
    y = zeros(nChan, nSamples);
    if useGPU
        y = gpuArray(y);
    end

    for i = 1:nChan
        y(i, :) = conv(x(i, :), B, 'same');
    end
else
    % For larger chunks, use frequency domain but with optimizations
    kernelLen = length(B);
    fftLen = 2^nextpow2(nSamples+kernelLen-1);

    % Pre-compute filter FFT
    B_fft = fft(B, fftLen);

    % Process all channels at once
    y = zeros(nChan, nSamples);
    if useGPU
        y = gpuArray(y);
        B_fft = gpuArray(B_fft);
    end

    for i = 1:nChan
        % Zero-pad input
        x_padded = [x(i, :), zeros(1, fftLen-nSamples)];
        if useGPU
            x_padded = gpuArray(x_padded);
        end

        % FFT-based convolution
        X_fft = fft(x_padded);
        Y_fft = X_fft .* B_fft;
        y_full = ifft(Y_fft);

        % Extract same-sized output (matching 'same' option)
        delay = floor(kernelLen/2);
        y(i, :) = y_full(delay+1:delay+nSamples);
    end
end
end