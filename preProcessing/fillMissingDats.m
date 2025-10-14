function fillMissingDats(varargin)
% [fillMissingDats(varagin)]
%
% [fills in missing dats and related files (analogin, digital in etc.) that
% are missing]
%
%
%  INPUTS
%  [parser]  [input parser for differing options, see below]
%   =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     ['basepath']  [basepath of missing dats. should be multiple
%                   basepaths]
%    =======================================================================
%
%
%  OUTPUTS
%     NA
%
%
%  SEE ALSO
%
% [HeathLarsson] [2021-2022]
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

p = inputParser;
otherdattypes = {'analogin'; 'digitalin'; 'auxiliary'; 'time'; 'supply'};
isFileType = @(x) sum(strcmp(x, otherdattypes)) == 1;
addParameter(p, 'basepath', pwd, @isfolder)
addParameter(p, 'fileType', [], isFileType)
addParameter(p, 'behaviorOnly', false, @islogical);

parse(p, varargin{:})
basepath = p.Results.basepath;
fileType = p.Results.fileType;
behaviorOnly = p.Results.behaviorOnly;

% get file types and data types
files_table = table();
files_table.files = {'amplifier', 'auxiliary', 'digitalin', 'digitalout', 'analogin', 'time', 'supply'}';
files_table.data_type = {'int16', 'uint16', 'uint16', 'uint16', 'uint16', 'int32', 'uint16'}';

%% Get session info
session = getSession('basepath', basepath); % Peter's sessionInfo
ampNch = session.extracellular.nChannels;

%% check location of each file and what to concat
typeFiles = dir([basepath, filesep, '*', filesep, fileType, '.dat']);

% For behavior-only sessions, use digitalin.dat as reference instead of amplifier.dat
if behaviorOnly || ampNch == 0
    disp('Behavior-only mode: using digitalin.dat as reference');
    ampFiles = dir([basepath, filesep, '*', filesep, 'digitalin.dat']);
    if isempty(ampFiles)
        warning('No digitalin.dat files found in subsessions. Cannot fill missing files.');
        return
    end
else
    % Standard mode: use amplifier.dat or continuous.dat
    ampFiles = dir([basepath, filesep, '*', filesep, 'amplifier.dat']);
    contFiles = dir([basepath, filesep, '**', filesep, 'continuous.dat']); %check for openEphys
    ampFiles = cat(1, ampFiles, contFiles);
end

typeFolders = {typeFiles.folder};
ampFolders = {ampFiles.folder};

fillInds = find(~ismember(ampFolders, typeFolders)); %index of folders that need fill

if isempty(fillInds)
    disp(['All subsessions already have ', fileType, '.dat files. No filling needed.']);
    return
end

% calculate number of channels for fill file
refTypeInd = find(ismember(typeFolders, ampFolders), 1);

if isempty(refTypeInd)
    warning(['Cannot find a subsession with both ', fileType, '.dat and reference file. Cannot determine channel count.']);
    return
end

refTypeSize = typeFiles(refTypeInd).bytes;
refAmpInd = find(strcmp(ampFolders, typeFiles(refTypeInd).folder), 1);

if isempty(refAmpInd)
    warning('Cannot find matching reference file. Cannot fill missing files.');
    return
end

refAmpSize = ampFiles(refAmpInd).bytes;

% Calculate number of channels based on file type
if behaviorOnly || ampNch == 0
    % For behavior-only, use digitalin as reference (16 channels, uint16)
    digitalNch = 16; % Standard Intan digital input channels
    refNch = digitalNch;
    refBytesPerSample = 2; % uint16
    typeNch = refTypeSize * refNch / refAmpSize;
else
    % Standard mode with amplifier channels
    typeNch = refTypeSize * ampNch / refAmpSize;
end

%% Fill in missing files
for ii = 1:length(fillInds)
    fillIdx = fillInds(ii);
    localAmpSize = ampFiles(fillIdx).bytes;
    
    % Calculate number of samples based on reference file
    if behaviorOnly || ampNch == 0
        % For behavior-only, reference is digitalin.dat (16 channels, uint16, 2 bytes)
        nPoints = localAmpSize / (refNch * refBytesPerSample);
    else
        % Standard mode with amplifier
        nPoints = localAmpSize / (ampNch * 2);
    end

    zeroData = zeros(typeNch, nPoints);

    filepath = [ampFiles(fillIdx).folder, filesep, fileType, '.dat'];
    disp(['Creating missing file: ', filepath]);
    fid = fopen(filepath, 'w');
    fwrite(fid, zeroData, files_table.data_type{contains(files_table.files, fileType)});
    fclose(fid);
end

disp(['Successfully filled ', num2str(length(fillInds)), ' missing ', fileType, '.dat files']);
