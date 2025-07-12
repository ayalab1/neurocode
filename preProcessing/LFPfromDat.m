function LFPfromDat(basepath, varargin)
% LFPFROMDAT Extract Local Field Potential (LFP) from wideband .dat file
%   [perform lowpass (2 X output Fs) sinc filter on wideband data
%   subsample the filtered data and save as a new flat binary
%   basename must have basename.dat and basename.xml
%   basepath is the full path for basename.dat]
%
% INPUTS
% [basePath]    [path where the recording files are located]
%   (options)
%    =========================================================================
%     Properties    Values
%    ------------------------------------------------------------------------
%       ['datFile']     [specify which file you'd like to compute lfp from]
%       ['outFs']       [(default: 1250) downsampled frequency of the .lfp output file]
%       ['lopass']      [(default: 450) low pass filter frequency]
%       ['useGPU']      [(default: false) whether to use GPU for processing]
%       ['inFs']        [input sampling frequency]
%       ['localDir']    [temporary directory for processing]
%       ['session']     [session structure]
%       ['maxMemoryGB'] [(default: auto) maximum memory to use in GB]
%       ['useMemMap']   [(default: auto) use memory mapping for very large files]
%
% OUTPUT
%   [Creates file:   basePath/baseName.lfp]
%
% [SMckenzie, BWatson, DLevenstein, kathrynmcclain] [2018-2022]
% Optimized version [2025]

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
addParameter(p, 'maxMemoryGB', [], @isnumeric);
addParameter(p, 'useMemMap', [], @islogical);

parse(p, varargin{:})
datFile = p.Results.datFile;
outFs = p.Results.outFs;
lopass = p.Results.lopass;
useGPU = p.Results.useGPU;
inFs = p.Results.inFs;
localDir = p.Results.localDir;
session = p.Results.session;
maxMemoryGB = p.Results.maxMemoryGB;
useMemMap = p.Results.useMemMap;

if isempty(session)
    session = getSession('basepath', basepath);
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
    fprintf('Using GPU: %s\n', g.Name);
end

sizeInBytes = 2;

%% Get file info and metadata
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
    warning('Low pass cutoff beyond Nyquist frequency')
end

ratio = lopass / (inFs / 2);
sampleRatio = (inFs / outFs);

% Output file
if ~isempty(localDir)
    flfp = fullfile(localDir, [basename, '.lfp']);
else
    flfp = fullfile(basepath, [basename, '.lfp']);
end

%% Check if file already exists
if exist([basepath, '\', basename, '.lfp'], 'file') || exist([basepath, '\', basename, '.eeg'], 'file')
    fprintf('LFP file already exists \n')
    return
end

%% Optimize chunk size based on available memory
nBytes = fInfo.bytes;
totalSamples = nBytes / (nbChan * sizeInBytes);

% Get available memory
if isempty(maxMemoryGB)
    try
        memInfo = memory;
        availableMemory = memInfo.MaxPossibleArrayBytes;
    catch
        availableMemory = 8e9; % Default to 8GB if memory info unavailable
    end
else
    availableMemory = maxMemoryGB * 1e9;
end

% Calculate optimal chunk size
% Account for: input data (int16), filtered data (double), output data (int16)
% Safety factor of 4 for intermediate calculations
memoryPerSample = nbChan * (2 + 8 + 2) * 4; % bytes per sample
maxChunkSize = floor(availableMemory/memoryPerSample);
maxChunkSize = min(maxChunkSize, 2e6); % Cap at 2M samples for reasonable progress updates

% Ensure chunk size is multiple of sample ratio
chunksize = maxChunkSize;
if mod(chunksize, sampleRatio) ~= 0
    chunksize = chunksize - mod(chunksize, sampleRatio);
end

fprintf('Optimized chunk size: %d samples (%.1f MB per chunk)\n', chunksize, chunksize*nbChan*sizeInBytes/1e6);

%% Set up memory mapping for very large files
if isempty(useMemMap)
    useMemMap = nBytes > 50e9; % Auto-enable for files > 50GB
end

% Buffer size for filter
ntbuff = 525; % default filter size in iosr toolbox
if mod(ntbuff, sampleRatio) ~= 0
    ntbuff = ntbuff + sampleRatio - mod(ntbuff, sampleRatio);
end

nbChunks = floor(totalSamples/chunksize);
if nbChunks == 0
    nbChunks = 1;
    chunksize = totalSamples;
end

%% Pre-allocate arrays for maximum efficiency
fprintf('Pre-allocating arrays...\n');
maxOutputSize = ceil(chunksize/sampleRatio);

% Pre-allocate main processing arrays - corrected dimensions
if useGPU
    try
        % For transposed data: [samples x channels]
        dat_gpu = gpuArray(zeros(chunksize+2*ntbuff, nbChan));
        filtered_gpu = gpuArray(zeros(chunksize+2*ntbuff, nbChan));
        fprintf('GPU arrays pre-allocated\n');
    catch ME
        warning(ME.message);
        useGPU = false;
    end
end

% CPU arrays not needed for pre-allocation since we create them on-the-fly

%% Set up file I/O
fprintf('Setting up optimized file I/O...\n');

if useMemMap && ~useGPU
    % Memory-mapped file access for very large files
    fprintf('Using memory-mapped file access for large dataset\n');
    try
        mmap = memmapfile(fdat, 'Format', {'int16', [nbChan, totalSamples], 'data'});
        useMemMapFile = true;
    catch
        warning('Memory mapping failed, using standard file I/O');
        useMemMapFile = false;
        fidI = fopen(fdat, 'r');
    end
else
    useMemMapFile = false;
    fidI = fopen(fdat, 'r');
end

fidout = fopen(flfp, 'w'); % Use 'w' instead of 'a' for better performance

%% Main processing loop - OPTIMIZED
fprintf('Starting optimized LFP extraction...\n');
tic;

for ibatch = 1:nbChunks

    % Progress reporting
    if mod(ibatch, max(1, floor(nbChunks/20))) == 0
        elapsed = toc;
        eta = elapsed * (nbChunks - ibatch) / ibatch;
        fprintf('Progress: %d/%d (%.1f%%) - ETA: %.1f min\n', ...
            ibatch, nbChunks, 100*ibatch/nbChunks, eta/60);
    end

    %% Read data chunk
    if useMemMapFile
        % Memory-mapped reading
        startIdx = (ibatch - 1) * chunksize + 1;
        if ibatch == 1
            endIdx = min(startIdx+chunksize+ntbuff-1, totalSamples);
            dat = double(mmap.Data.data(:, startIdx:endIdx));
        else
            startIdx = startIdx - ntbuff;
            endIdx = min(startIdx+chunksize+2*ntbuff-1, totalSamples);
            dat = double(mmap.Data.data(:, startIdx:endIdx));
        end
    else
        % Standard file reading
        if ibatch > 1
            fseek(fidI, ((ibatch - 1) * (nbChan * sizeInBytes * chunksize))-(nbChan * sizeInBytes * ntbuff), 'bof');
            dat = fread(fidI, nbChan*(chunksize + 2 * ntbuff), 'int16');
        else
            dat = fread(fidI, nbChan*(chunksize + ntbuff), 'int16');
        end

        try
            if ibatch > 1
                dat = reshape(dat, [nbChan, (chunksize + 2 * ntbuff)]);
            else
                dat = reshape(dat, [nbChan, (chunksize + ntbuff)]);
            end
            dat = double(dat);
        catch
            warning('Reshape failed at chunk %d, trying file reload...', ibatch);
            if ~useMemMapFile
                fclose(fidI);
                fidI = fopen(fdat, 'r');
                if ibatch > 1
                    fseek(fidI, ((ibatch - 1) * (nbChan * sizeInBytes * chunksize))-(nbChan * sizeInBytes * ntbuff), 'bof');
                    dat = fread(fidI, nbChan*(chunksize + 2 * ntbuff), 'int16');
                    dat = reshape(dat, [nbChan, (chunksize + 2 * ntbuff)]);
                else
                    dat = fread(fidI, nbChan*(chunksize + ntbuff), 'int16');
                    dat = reshape(dat, [nbChan, (chunksize + ntbuff)]);
                end
                dat = double(dat);
            end
        end
    end

    %% OPTIMIZED FILTERING - Transpose for better performance
    % sincFilter processes columns, so transpose to make time dimension first
    dat_transposed = dat'; % Now [samples x channels]

    if useGPU
        % GPU processing
        dat_gpu(1:size(dat_transposed, 1), 1:size(dat_transposed, 2)) = dat_transposed;

        % Apply filter - now each channel is a column, processed in parallel
        filtered_gpu(1:size(dat_transposed, 1), 1:size(dat_transposed, 2)) =...
            iosr.dsp.sincFilter(dat_gpu(1:size(dat_transposed, 1), 1:size(dat_transposed, 2)), ratio);

        % Transpose back to [channels x samples]
        filtered_data = filtered_gpu(1:size(dat_transposed, 1), 1:size(dat_transposed, 2))';

        % Downsample
        if ibatch == 1
            DATA = int16(real(filtered_data(:, sampleRatio:sampleRatio:end-ntbuff)));
        else
            DATA = int16(real(filtered_data(:, ntbuff+sampleRatio:sampleRatio:end-ntbuff)));
        end

        % Transfer back to CPU for writing
        DATA = gather(DATA);

    else
        % CPU processing - transpose for column-wise processing
        filtered_transposed = iosr.dsp.sincFilter(dat_transposed, ratio);
        filtered_data = filtered_transposed'; % Back to [channels x samples]

        % Downsample
        if ibatch == 1
            DATA = int16(real(filtered_data(:, sampleRatio:sampleRatio:end-ntbuff)));
        else
            DATA = int16(real(filtered_data(:, ntbuff+sampleRatio:sampleRatio:end-ntbuff)));
        end
    end

    %% Write data
    fwrite(fidout, DATA, 'int16');
end

%% Process remainder
remainder = totalSamples - nbChunks * chunksize;
if remainder > 0
    fprintf('Processing remainder: %d samples\n', remainder);

    if useMemMapFile
        startIdx = nbChunks * chunksize - ntbuff + 1;
        endIdx = totalSamples;
        dat = double(mmap.Data.data(:, startIdx:endIdx));
    else
        fseek(fidI, ((nbChunks - 1) * (nbChan * sizeInBytes * chunksize))-(nbChan * sizeInBytes * ntbuff), 'bof');
        dat = fread(fidI, nbChan*(remainder + ntbuff), 'int16');
        dat = reshape(dat, [nbChan, (remainder + ntbuff)]);
        dat = double(dat);
    end

    % Process remainder with corrected orientation
    if useGPU
        dat_transposed = dat';
        dat_gpu(1:size(dat_transposed, 1), 1:size(dat_transposed, 2)) = dat_transposed;
        filtered_gpu(1:size(dat_transposed, 1), 1:size(dat_transposed, 2)) =...
            iosr.dsp.sincFilter(dat_gpu(1:size(dat_transposed, 1), 1:size(dat_transposed, 2)), ratio);
        filtered_data = filtered_gpu(1:size(dat_transposed, 1), 1:size(dat_transposed, 2))';
        output_size = floor(remainder/sampleRatio);
        DATA = gather(int16(real(filtered_data(:, ntbuff+sampleRatio:sampleRatio:end))));
        DATA = DATA(:, 1:output_size);
    else
        dat_transposed = dat';
        filtered_transposed = iosr.dsp.sincFilter(dat_transposed, ratio);
        filtered_data = filtered_transposed';
        output_size = floor(remainder/sampleRatio);
        DATA = int16(real(filtered_data(:, ntbuff+sampleRatio:sampleRatio:end)));
        DATA = DATA(:, 1:output_size);
    end

    fwrite(fidout, DATA, 'int16');
end

%% Cleanup
if useMemMapFile
    clear mmap;
else
    fclose(fidI);
end
fclose(fidout);

if useGPU
    reset(g);
    gpuDevice([]);
end

elapsed = toc;
fprintf('LFP extraction completed in %.1f minutes\n', elapsed/60);
fprintf('Processing speed: %.1f MB/min\n', (nBytes / 1e6)/(elapsed / 60));

% Move file if using local directory
if ~isempty(localDir)
    movefile(flfp, fullfile(basepath, [basename, '.lfp']));
end

fprintf('LFP file created: %s\n', fullfile(basepath, [basename, '.lfp']));

end