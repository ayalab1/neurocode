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
% Corrected chunking and filtering logic to fix file size mismatch bug
% A chunk size that's a multiple of the sample ratio is used for clean downsampling.
chunksize = 1e5;
if mod(chunksize, sampleRatio) ~= 0
    chunksize = chunksize - mod(chunksize, sampleRatio);
    if chunksize < sampleRatio
        chunksize = sampleRatio;
    end
end

% Dynamically adjust filter order (N) and buffer size (ntbuff)
N = max(filterOrder, round(3*sampleRatio));
if mod(N, 2) == 0
    N = N + 1;
end
ntbuff = N;

totalSamples = fInfo.bytes / (nbChan * sizeInBytes);

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

fidout = fopen(flfp, 'w');
if fidout == -1
    error('LFPfromDat:FileOpenError', 'Could not create output .lfp file: %s', flfp);
end

fprintf('Starting LFP extraction...\n')
tic;

processedSamples = 0;
pre_buffer = []; % No pre-buffer for the first chunk

WaitMessage = parfor_wait(ceil(totalSamples/chunksize),...
    'ReportInterval', 20);

while processedSamples < totalSamples
    % Determine how many samples to read this iteration
    remainingSamples = totalSamples - processedSamples;
    readSamples = min(chunksize, remainingSamples);

    % Read chunk data
    fseek(fidI, (processedSamples - size(pre_buffer, 2))*nbChan*sizeInBytes, 'bof');
    current_chunk = fread(fidI, [nbChan, readSamples + size(pre_buffer, 2)], 'int16');

    % If this is the last chunk, it might read fewer samples than requested
    if size(current_chunk, 2) < readSamples + size(pre_buffer, 2)
        current_chunk = reshape(current_chunk, [nbChan, numel(current_chunk) / nbChan]);
    end

    dataToFilter = double(current_chunk);

    % Filter the data
    if useGPU
        dataToFilter = gpuArray(dataToFilter);
        filtered_data = fast_sinc_filter_matrix(dataToFilter, ratio, N, useGPU);
    else
        filtered_data = fast_sinc_filter_matrix(dataToFilter, ratio, N, useGPU);
    end

    % Extract the downsampled portion corresponding to the *original* chunk
    % not including the pre-buffer
    if isempty(pre_buffer)
        downsampled_chunk = filtered_data(:, sampleRatio:sampleRatio:readSamples);
    else
        downsampled_chunk = filtered_data(:, size(pre_buffer, 2)+sampleRatio:sampleRatio:size(pre_buffer, 2)+readSamples);
    end

    % Update the pre-buffer for the next iteration
    pre_buffer = double(current_chunk(:, end-ntbuff+1:end));

    if useGPU
        downsampled_chunk = gather(int16(real(downsampled_chunk)));
    else
        downsampled_chunk = int16(real(downsampled_chunk));
    end

    % Write to file
    fwrite(fidout, downsampled_chunk(:), 'int16');
    processedSamples = processedSamples + readSamples;

    % if mod(processedSamples, chunksize*10) == 0
    %     fprintf('%d percent complete\n', round(100*processedSamples/totalSamples));
    % end
    WaitMessage.Send;
end
WaitMessage.Destroy;

fclose(fidI);
fclose(fidout);

if useGPU
    reset(g);
end

elapsed = toc;
fprintf('\nLFP extraction completed in %.2f seconds\n', elapsed);

if ~isempty(localDir) && ~strcmp(localDir, basepath)
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
