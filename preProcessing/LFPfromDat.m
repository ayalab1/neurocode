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
%   'filterOrder' - Order of the sinc filter (controls sharpness and kernel length) (default: 525, but dynamically adjusted)
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
%   - Parallel Computing Toolbox (for GPU support)
%   - getSession.m (from CellExplorer)
%
%
%   Note: Will skip processing if .lfp file already exists

%% Input handling
if ~exist('basepath', 'var')
    basepath = pwd;
end

p = inputParser;
addParameter(p, 'datFile', [], @ischar);
addParameter(p, 'outFs', [], @isnumeric);
addParameter(p, 'lopass', 450, @isnumeric);
addParameter(p, 'useGPU', true, @islogical);
addParameter(p, 'inFs', [], @isnumeric);
addParameter(p, 'localDir', [], @isfolder);
addParameter(p, 'session', [], @isstruct);
addParameter(p, 'filterOrder', 525, @isnumeric); % New parameter for filter order

parse(p, varargin{:})
datFile = p.Results.datFile;
outFs = p.Results.outFs;
lopass = p.Results.lopass;
useGPU = p.Results.useGPU;
inFs = p.Results.inFs;
localDir = p.Results.localDir;
session = p.Results.session;
filterOrder = p.Results.filterOrder;

if isempty(session)
    try
        session = getSession('basepath', basepath);
    catch ME
        error('LFPfromDat:SessionLoadError', 'Could not load session metadata. Please ensure getSession.m is available and a .session.mat file exists in the basepath. Error: %s', ME.message);
    end
end

basename = session.general.name;

if isempty(datFile)
    datFile = [basename, '.dat'];
elseif ~strcmp(datFile(end-3:end), '.dat')
    datFile = [datFile, '.dat'];
end

% GPU setup
if useGPU && gpuDeviceCount < 1
    warning('No GPU device found, continuing without..')
    useGPU = false;
end

if useGPU
    g = gpuDevice(1);
end

sizeInBytes = 2; % int16 is 2 bytes

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

% Handle lopass > outFs / 2 more robustly.
% If lopass is above Nyquist, set it to Nyquist and warn.
if lopass >= outFs / 2
    warning('LFPfromDat:NyquistWarning', 'Low-pass cutoff (%.1f Hz) is at or beyond the Nyquist frequency (%.1f Hz) for the output sampling rate (%.1f Hz). Setting lopass to Nyquist - epsilon.', lopass, outFs/2, outFs);
    lopass = (outFs / 2) - 1; % Set to just below Nyquist to avoid issues
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

% Dynamically adjust filter order (N) and buffer size (ntbuff)
% A common heuristic for filter order is to relate it to the sampling rate
% and cutoff frequency, or to use a fixed number of cycles.
% Here, we use the 'filterOrder' parameter but also ensure it's reasonable
% relative to the `sampleRatio`.
% The buffer size `ntbuff` should be at least half the filter length to correctly handle
% edge effects. A safe approach is to make it slightly larger than filter length.
N = max(filterOrder, round(3*sampleRatio)); % Ensure minimum filter order based on downsample ratio
if mod(N, 2) == 0 % Ensure N is odd for sinc filter symmetry around 0
    N = N + 1;
end
ntbuff = N + round(sampleRatio); % Buffer size larger than filter length
if mod(ntbuff, sampleRatio) ~= 0
    ntbuff = ntbuff + sampleRatio - mod(ntbuff, sampleRatio);
end

nBytes = fInfo.bytes;
nbChunks = floor(nBytes/(nbChan * sizeInBytes * chunksize)); % Adjusted calculation for last chunk

%% Check if file exists
if exist([basepath, filesep, basename, '.lfp'], 'file') || ...
        exist([basepath, filesep, basename, '.eeg'], 'file')
    fprintf('LFP file already exists \n')
    return
end

%% Main processing
fidI = fopen(fdat, 'r');
if fidI == -1
    error('LFPfromDat:FileOpenError', 'Could not open input .dat file: %s', fdat);
end

fidout = fopen(flfp, 'w'); % Change 'a' to 'w' to ensure a fresh file
if fidout == -1
    error('LFPfromDat:FileOpenError', 'Could not create output .lfp file: %s', flfp);
end

fprintf('Starting LFP extraction...\n')
tic;

% Loop through chunks
for ibatch = 0:nbChunks - 1 % Start from 0 to simplify fseek logic
    if mod(ibatch, 10) == 0
        if ibatch ~= 0
            fprintf(repmat('\b', [1, length([num2str(round(100*(ibatch - 10)/nbChunks)), ' percent complete'])]))
        end
        fprintf('%d percent complete', round(100*ibatch/nbChunks));
    end

    % Calculate read position and size
    startByte = (ibatch * nbChan * sizeInBytes * chunksize) - (nbChan * sizeInBytes * ntbuff);
    if ibatch == 0
        startByte = 0; % No pre-buffer for the first chunk
        readSamples = nbChan * (chunksize + ntbuff);
    else
        readSamples = nbChan * (chunksize + 2 * ntbuff);
    end

    fseek(fidI, startByte, 'bof');
    dat = fread(fidI, readSamples, 'int16');

    % Handle partial read at the end if the last chunk is smaller
    if numel(dat) < readSamples
        warning('LFPfromDat:PartialRead', 'Partial read for chunk %d. Data might be truncated.', ibatch);
        % Adjust reshape size based on actual read data
        dat = reshape(dat, [nbChan, numel(dat) / nbChan]);
    else
        dat = reshape(dat, [nbChan, readSamples / nbChan]);
    end

    dat = double(dat);

    % Filter
    if useGPU
        dat = gpuArray(dat);
        filtered_data = fast_sinc_filter_matrix(dat, ratio, N, useGPU); % Pass N
    else
        filtered_data = fast_sinc_filter_matrix(dat, ratio, N, useGPU); % Pass N
    end

    % Downsample and extract relevant portion (handle buffer)
    if ibatch == 0
        % For the first chunk, take from sampleRatio:end-ntbuff
        processed_chunk = filtered_data(:, sampleRatio:sampleRatio:end-ntbuff);
    else
        % For subsequent chunks, take from ntbuff+sampleRatio:end-ntbuff
        processed_chunk = filtered_data(:, ntbuff+sampleRatio:sampleRatio:end-ntbuff);
    end

    if useGPU
        processed_chunk = gather(int16(real(processed_chunk)));
    else
        processed_chunk = int16(real(processed_chunk));
    end

    % Write to file
    fwrite(fidout, processed_chunk(:), 'int16');
end

% Handle the very last chunk (remainder) separately
% The previous loop structure often leads to double-counting or missing data
% at the end. A simpler approach is to handle it after the main loop.
% Calculate actual samples remaining after full chunks
totalSamples = nBytes / (sizeInBytes * nbChan);
processedSamples = nbChunks * chunksize;
remainder = totalSamples - processedSamples;

if remainder > 0
    fprintf('\nProcessing remainder chunk...\n');
    % Seek to the start of the remainder, with appropriate pre-buffer
    startByte = (nbChunks * nbChan * sizeInBytes * chunksize) - (nbChan * sizeInBytes * ntbuff);
    fseek(fidI, startByte, 'bof');

    readSamples = nbChan * (remainder + ntbuff);
    dat = fread(fidI, readSamples, 'int16');
    if numel(dat) < readSamples
        warning('LFPfromDat:PartialReadEnd', 'Partial read for final remainder chunk. Data might be truncated.');
        dat = reshape(dat, [nbChan, numel(dat) / nbChan]);
    else
        dat = reshape(dat, [nbChan, readSamples / nbChan]);
    end

    dat = double(dat);

    if useGPU
        dat = gpuArray(dat);
        filtered_data = fast_sinc_filter_matrix(dat, ratio, N, useGPU);
    else
        filtered_data = fast_sinc_filter_matrix(dat, ratio, N, useGPU);
    end

    % For remainder, take from ntbuff+sampleRatio to the actual end
    processed_remainder = filtered_data(:, ntbuff+sampleRatio:sampleRatio:end);

    if useGPU
        processed_remainder = gather(int16(real(processed_remainder)));
    else
        processed_remainder = int16(real(processed_remainder));
    end

    fwrite(fidout, processed_remainder(:), 'int16');
end


fclose(fidI);
fclose(fidout);

if useGPU
    reset(g);
end

elapsed = toc;
fprintf('\nLFP extraction completed in %.2f seconds\n', elapsed);

if ~isempty(localDir) && ~strcmp(localDir, basepath) % Only move if localDir is different from basepath
    try
        movefile(flfp, fullfile(basepath, [basename, '.lfp']));
        disp('LFP file moved to basepath.')
    catch ME
        warning('LFPfromDat:MoveFileError', 'Could not move LFP file from %s to %s. Error: %s', flfp, fullfile(basepath, [basename, '.lfp']), ME.message);
    end
end

disp('LFP file created')
end

function y = fast_sinc_filter_matrix(x, ratio, N, useGPU)
%FAST_SINC_FILTER_MATRIX Optimized matrix-based sinc filtering
% This version uses vectorized operations but avoids FFT overhead for moderate sizes

[nChan, nSamples] = size(x);

% Create sinc kernel - Fix 1: N is now passed as an argument
n = -floor(N/2):floor(N/2);
B = sinc(ratio*n) .* ratio;

if useGPU
    B = gpuArray(B);
end


% For larger chunks, use frequency domain
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
    % For 'same' convolution, the delay is floor(kernelLen/2)
    delay = floor(kernelLen/2);
    y(i, :) = y_full(delay+1:delay+nSamples);
end
end
