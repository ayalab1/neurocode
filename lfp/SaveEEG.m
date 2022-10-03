function SaveEEG(varargin)

% This function will make a new lfp file (.eeg) in which the new reference
% is the mean of all other electrodes. 

%Using the mean of all channels as a reference is a common practice in EEG studies ("average reference"), which
% is why it may be fitting to name that file .eeg and keep the original
% .lfp file for other analyses (such as estimating the EMG based on 
% cross-shank correlations, which would not work as well on this new 
% re-referenced lfp file). 
% Note: This function requires a [basename].lfp file and a '[basename].session.mat'
% file (to get the number of channels) in the [basepath] folder.
%
% USAGE
%
%    SaveEEG(<options>)
%
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties        Values
%    -------------------------------------------------------------------------
%     'basepath'        the folder where the lfp file is stored (default = 
%                       current directory)
%     'rejectChannels'  a list of channels (1-indexing, so add 1 to the number
%                       in neuroscope) which will be ignored when computing
%                       the mean
%    =========================================================================
%
%
% Copyright (C) 2022 Ralitsa Todorova
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

p = inputParser;
addParameter(p,'basepath',pwd,@isfolder);
addParameter(p,'rejectChannels',[]);
parse(p,varargin{:}); % add the optional inputs to replace the defaults
% assign all varriabler from "p.Results"
fields = fieldnames(p.Results);
cellfun(@(x,y) assignin('caller', x, y), fields, struct2cell(p.Results)); 
basename = basenameFromBasepath(basepath);
load(fullfile(basepath,[basename '.session.mat']),'session');
nChannels = session.extracellular.nChannels;
okChannels = find(~ismember((1:nChannels)',rejectChannels(:)));
eegFile = fullfile(basepath,[basename '.eeg']); % This is where we will store the normalized LFP file
if ~exist(eegFile,'file')
    lfpFile = fullfile(basepath,[basename '.lfp']);
    copyfile(lfpFile,eegFile);
    file = memmapfile(eegFile,'Format','int16','Writable',true);
    data = reshape(file.Data,nChannels,[]);
    m = int16(mean(data(okChannels,:)));
    newData = data;
    newData(okChannels,:) = bsxfun(@minus,data(okChannels,:),m);
    file.data = newData(:);
    clear file
else
    disp(['EEG file exists already. Exiting...']);
end



