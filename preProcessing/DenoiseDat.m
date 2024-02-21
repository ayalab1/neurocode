function DenoiseDat(datFile,session,varargin)

% DenoiseDat - remove the first PCA ('noise') component from your .dat file
% 
% WARNING: THIS FUNCTION WILL MODIFY YOUR .DAT FILE! Make sure original files are backed up!
%
% This function computes PCA (each shank is treated separately) over a
% short (preferably noisy) sample of the data (from the period "baseline", 
% user-selected to include noisy data). Subsequently, the part of the 
% signal that is explained by this "noisy" component is removed, resulting
% in a denoised ('cleaned') .dat file, which should improve spike sorting.
% Note: the selection of a "baseline" period containing the problematic noise
% (physical or electrical noise) is crucial for the denoising to work. 
%
% Note that this manipulation will also remove any real signal shared among
% channels, such as important oscillations. It is therefore recommended
% to NOT apply this de-noising to the original dat files as this will 
% permanently modify them and delete this signal. Instead, it is recommended
% that you preserve your original files in their subsession folders, and 
% you apply this denoising step to the .dat file only after the LFP file
% has been created (the denoised .dat file will be relatively flat, missing
% a lot of oscillation information of interest for the .lfp). 
%
%  USAGE
%
%    DenoiseDat(session, datFile, <options>)
%
%  INPUTS
%  
%    datFile        full path towards the file to be denoised. 
%                   WARNING: THIS FILE WILL BE MODIFIED! 
%                   Make sure original files are backed up!
%    session        Cell Explorer format variable containing parameters about 
%                   the recording (sampling rate, shank configuration, number 
%                   of channels)
%
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties         Values
%    -------------------------------------------------------------------------
%   'SSD_path'           path to SSD. If provided, .dat file will be copied
%                        and edited there for faster performance, and 
%                        afterwards copied back, overwriting the original file.
%   'baseline'           interval (in seconds) that will be used to compute
%                        the PCA. It is recommended that this includes a 
%                        noisy period (default = whole recording)
%   'sampleDuration'     duration (in seconds) of the random sample that 
%                        will be selected from the "baseline" period to 
%                        compute the PCA (default = 20). Note that longer 
%                        periods increase memory demands. 
%   'secondsPerChunk'    duration (in seconds) of a chunk of data (default = 10). 
%                        After computing the noise component, noise is 
%                        progressively removed from the .dat file in chunks. 
%                        Note that larger chunks increase memory demands. 
%   'rejectChannels'     channels (base1) that should be excluded from denoising
%   'verbose'            Display progress messages (default = true)
%    =========================================================================
%
%  EXAMPLE 
%  
%    LFPfromDat; % Produce LFP before denoising the .dat file to preserve slow oscillations affecting all channels (e.g. delta waves)
%    DenoiseDat(session,datFile,'SSD','D:\'); % remove first PCA 'noise' component from .dat file to improve spike sorting 
%    kilosortFolder = KiloSortWrapper('SSD_path','D:\','rejectchannels',excludeChannels,'datFilename',SSD_file); % spike sort
%
% Ralitsa Todorova and Weiwei Chen 2024
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.
%-------------------------------------------------------------------------


p = inputParser;
addParameter(p, 'SSD_path', [], @isfolder);
addParameter(p, 'baseline', [], @isdmatrix); % interval to use as basis for PCA
addParameter(p, 'sampleDuration', 20, @isdscalar); % cap duration of the noise random sample to be used to compute PCA
addParameter(p, 'secondsPerChunk', 10, @isdscalar);
addParameter(p, 'rejectChannels', [], @isivector);
addParameter(p, 'verbose', true, @islogical);
parse(p, varargin{:});

SSD_path = p.Results.SSD_path;
baseline = p.Results.baseline;
sampleDuration = p.Results.sampleDuration;
secondsPerChunk = p.Results.secondsPerChunk;
rejectChannels = p.Results.rejectChannels; %!
verbose = p.Results.verbose;

%% Load information from the session variable
nChannels = session.extracellular.nChannels;
samplingRate = session.extracellular.sr;
shanks = session.extracellular.spikeGroups.channels; nShanks = length(shanks);
for i=1:nShanks, shanks{i}(ismember(shanks{i},rejectChannels)) = []; end % remove any channels from this process
basepath = session.general.basePath;
basename = basenameFromBasepath(basepath);
if ~exist(datFile,'file'), error([datFile ' file does not exist.']); end

if isempty(SSD_path), SSD = false; else, SSD = true; end
if SSD, SSD_file = fullfile(SSD_path,[basename '.dat']);
    copyfile(datFile,SSD_file);
    if verbose, disp([datestr(clock) ': dat file copied for ' basepath '. Finding noisy periods...']); end
    m = memmapfile(SSD_file, 'Format','int16','Writable',true);
else % work directly on the .dat file
    m = memmapfile(datFile, 'Format','int16','Writable',true);
end

data = reshape(m.data,nChannels,[]);
nSamples = size(data,2);

% Select timestamp indices that will be taken as baseline to compute PCA:
if isempty(baseline), baseline = [0 nSamples/samplingRate]; end % use the whole recording as a "baseline"
baselineSamples = round(baseline*samplingRate);
if sum(diff(baselineSamples,[],2))>samplingRate*sampleDuration % baseline is larger than the required number of indices
    % select a random sample of indices within the baseline period
    idx = sort(randperm(round(sum(diff(baselineSamples,[],2))),round(samplingRate*sampleDuration)))'; % select random indices among all possible indices
    idx = round(Unshift(idx,baselineSamples)); % shift them to the corresponding indices within the baseline
else
    if isdvector(baselineSamples)
        idx = baselineSamples(1):baselineSamples(2);
    else
        idx = linspaceVector(baselineSamples(:,1),baselineSamples(:,2)); % all possible indices within the baseline period
    end
end

if verbose, disp([datestr(clock) ': Computing PCA over noise periods in ' basepath '...']); end
% Compute the PCA shank by shank
V = cell(nShanks,1); vs = cell(nShanks,1);
baselineData = data(:,idx);
for i=1:nShanks
    channels = shanks{i};
    shankData = baselineData(channels,:);
    c = corr(double(shankData'));
    [v,~] = eig(c);
    V{i} = zeros(nChannels,1); V{i}(channels) = v(:,1);
    vs{i} = v(:,1);
end
if verbose, disp([datestr(clock) ': PCA computed. Removing noise components in ' datFile ' file...']); end
%%
% Modify the .dat file in chunks:
chunkSize = ceil(samplingRate*secondsPerChunk);
indicesChunk1 = (1:chunkSize*nChannels)'; % these are the .dat file indices for the first chunk
nChunks = ceil(nSamples/chunkSize);
for chunk = 1:nChunks
    indicesChunk = indicesChunk1+(chunkSize*nChannels)*(chunk-1);
    dataChunk = double(reshape(m.Data(indicesChunk),nChannels,[])); % Consider adding a constant if you're re-running denoising!
    noise = zeros(size(dataChunk),'int16');
    for i=1:nShanks
        projection = dataChunk'*V{i};
        % project back to channels:
        noise(shanks{i},:) = int16(bsxfun(@times,projection',vs{i}));
    end
    m.Data(indicesChunk) =   m.Data(indicesChunk)-noise(:);
end

if SSD,
    if verbose, disp([datestr(clock) ': Noise removed. Copying back local dat file to ' basepath '...']); end
    copyfile(SSD_file,datFile);
else
    disp([datestr(clock) ': Noise removed!']);
end

clear m





















