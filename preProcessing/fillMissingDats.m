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

parse(p, varargin{:})
basepath = p.Results.basepath;
fileType = p.Results.fileType;

% get file types and data types
files_table = table();
files_table.files = {'amplifier', 'auxiliary', 'digitalin', 'digitalout', 'analogin', 'time', 'supply'}';
files_table.data_type = {'int16', 'uint16', 'uint16', 'uint16', 'uint16', 'int32', 'uint16'}';

%% Get session info
session = getSession('basepath', basepath); % Peter's sessionInfo
basename = session.general.name;
ampNch = session.extracellular.nChannels;

%% check location of each file and what to concat
typeFiles = dir([basepath, filesep, '*', filesep, fileType, '.dat']);
ampFiles = dir([basepath, filesep, '*', filesep, 'amplifier.dat']);
typeFolders = {typeFiles.folder};
ampFolders = {ampFiles.folder};

fillInds = find(~ismember(ampFolders, typeFolders)); %index of folders that need fill

% calculate number of channels to for fill file
refTypeInd = find(ismember(typeFolders, ampFolders), 1);
refTypeSize = typeFiles(refTypeInd).bytes;
refAmpInd = find(strcmp(ampFolders, typeFiles(refTypeInd).folder));
refAmpSize = ampFiles(refAmpInd).bytes;
typeNch = refTypeSize * ampNch / refAmpSize;

%% Fill in missing files
for ii = 1:length(fillInds)
    fillIdx = fillInds(ii);
    localAmpSize = ampFiles(fillIdx).bytes;
    nPoints = localAmpSize / (ampNch * 2);

    zeroData = zeros(typeNch, nPoints);

    filepath = [ampFiles(fillIdx).folder, filesep, fileType, '.dat'];
    fid = fopen(filepath, 'w');
    fwrite(fid, zeroData, files_table.data_type{contains(files_table.files, fileType)});
    fclose(fid);
end
