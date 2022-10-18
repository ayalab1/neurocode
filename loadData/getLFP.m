function [lfp,info] = getLFP(varargin)
%[ getLFP] - Get local field potentials.
%
%  Load local field potentials from disk.
%
%    =========================================================================
%  USAGE
%
%    [lfp] = getLFP(channels,<options>)
%
%  INPUTS
%
%    channels(required) -must be first input, numeric row vector
%                        of channels to load or use keyword 'all' for all
%                        channels. Input is 1-based indexing.
%  Name-value paired inputs:
%    'basepath'           - folder in which .lfp file will be found (default
%                           is pwd)
%                           folder should follow buzcode standard:
%                           whateverPath/baseName
%                           and contain file baseName.lfp
%    'basename'           -base file name to load
%    'intervals'          -list of time intervals [0 10; 20 30] to read from
%                           the LFP file (default is [0 inf])
%    'downsample'         -factor to downsample the LFP (i.e. 'downsample',5
%                           will load a 1250Hz .lfp file at 250Hz)
%    'noPrompts'          -logical (default) to supress any user prompts
%    'fromDat'            -option to load directly from .dat file (default:false)
%
%    =========================================================================
%  OUTPUT
%
%    lfp              [Nt x  1 + Nd] matrix of the LFP data, where the first
%                     column is the timestamps (in seconds) of the lfp data
%                     and each subsequent column is the LFP data for each
%                     channel.
%    info             structure with metadata describing the lfp
%    .channels        [Nd X 1] vector of channel ID's
%    .samplingRate    LFP sampling rate [default = 1250]
%    .intervals       [N x 2] matrix of start/stop times of LFP interval
%
%
%  EXAMPLES
%
%    % channel ID 5 (= # 6), from 0 to 120 seconds
%    lfp = getLFP(5,'restrict',[0 120]);
%    % same, plus from 240.2 to 265.23 seconds
%    lfp = getLFP(5,'restrict',[0 120;240.2 265.23]);
%    % multiple channels
%    lfp = getLFP([1 2 3 4 10 17],'restrict',[0 120]);
%    % channel # 3 (= ID 2), from 0 to 120 seconds
%    lfp = getLFP(3,'restrict',[0 120],'select','number');

%    =========================================================================
% Copyright (C) 2004-2011 by Michaël Zugaro, 2017 David Tingley,
% 2020 kathryn mcclain, 2022 Ralitsa Todorova & Laura Berkowitz
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.
%
% TODO
% expand channel selection options (i.e. region or spikegroup)
%% Parse the inputs!

channelsValidation = @(x) isnumeric(x) || strcmp(x,'all');

% parse args
p = inputParser;
addRequired(p,'channels',channelsValidation)
addParameter(p,'intervals',[],@isnumeric)
addParameter(p,'restrict',[],@isnumeric)
addParameter(p,'basepath',pwd,@isfolder);
addParameter(p,'downsample',1,@isnumeric);
addParameter(p,'noPrompts',false,@islogical);
addParameter(p,'fromDat',false,@islogical);

parse(p,varargin{:})
channels = p.Results.channels;
downsamplefactor = p.Results.downsample;
basepath = p.Results.basepath;
noPrompts = p.Results.noPrompts;
fromDat = p.Results.fromDat;

basename = basenameFromBasepath(basepath);

% doing this so you can use either 'intervals' or 'restrict' as parameters to do the same thing
intervals = p.Results.intervals;
restrict = p.Results.restrict;
if isempty(intervals) && isempty(restrict) % both empty
    intervals = [0 Inf];
elseif isempty(intervals) && ~isempty(restrict) % intervals empty, restrict isn't
    intervals = restrict;
end

%% let's check that there is an appropriate LFP file
if isempty(basename)
    %disp('No basename given, so we look for a *lfp/*eeg file...')
    switch fromDat
        case false
            d = dir([basepath filesep '*lfp']);
        case true
            d = dir([basepath filesep '*dat']);
    end
    if length(d) > 1 % we assume one .lfp file or this should break
        error('there is more than one .lfp file in this directory?');
    elseif length(d) == 0
        d = dir([basepath filesep '*eeg']);
        if isempty(d)
            error('could not find an lfp/eeg file..')
        end
        structure.Filename = d.name;
        basename = strsplit(structure.Filename,'.');
        if length(basename) > 2
            base = [];
            for i=1:length(basename)-1
                base = [base basename{i} '.'];
            end
            basename = base(1:end-1);  % this is an fugly hack to make things work with Kenji's naming system...
        else
            basename = basename{1};
        end
    end
else
    switch fromDat
        case false
            d = dir([basepath filesep basename '.lfp']);
        case true
            d = dir([basepath filesep basename '.dat']);
    end
    
    if length(d) > 1 % we assume one .lfp file or this should break
        error('there is more than one .lfp file in this directory?');
    elseif length(d) == 0
        d = dir([basepath filesep basename '.eeg']);
        if isempty(d)
            error('could not find an lfp/eeg file..')
        end
    end
    structure.Filename = d.name;
end

%% things we can parse from sessionInfo or xml file

session = getSession('basepath',basepath);

if fromDat
    samplingRate = session.extracellular.sr;
else
    samplingRate = session.extracellular.srLfp;
end

samplingRateLFP_out = samplingRate./downsamplefactor;
structure.samplingRate = samplingRateLFP_out;

if mod(samplingRateLFP_out,1)~=0
    error('samplingRate/downsamplefactor must be an integer')
end


%% Channel load options
%Right now this assumes that all means channels 0:nunchannels-1 (neuroscope
%indexing), we could also add options for this to be select region or spike
%group from the xml...
if strcmp(channels,'all')
    chInfo = hackInfo('basepath',basepath);
    channels = chInfo.one.channels;
else
    %Put in something here to collapse into X-Y for consecutive channels...
    display(['Loading Channels ',num2str(channels),' (1-indexing)'])
end

%% get the data
disp('loading LFP file...')
nIntervals = size(intervals,1);
% returns lfp/bz format
for i = 1:nIntervals
    structure(i).duration = (intervals(i,2)-intervals(i,1));
    structure(i).interval = [intervals(i,1) intervals(i,2)];
    
    % Load data and put into struct
    % we assume 0-indexing like neuroscope, but loadBinary uses 1-indexing to
    % load....
    structure(i).data = loadBinary([basepath filesep structure.Filename],...
        'duration',double(structure(i).duration),...
        'frequency',samplingRate,'nchannels',session.extracellular.nChannels,...
        'start',double(structure(i).interval(1)),'channels',channels,...
        'downsample',downsamplefactor);
    structure(i).timestamps = (1:length(structure(i).data))'/samplingRateLFP_out;
    % check if duration is inf, and reset to actual duration...
    if structure(i).interval(2) == inf
        structure(i).interval(2) = length(structure(i).timestamps)/structure(i).samplingRate;
    end
    if structure(i).interval(1)>0 % shift the timestamps accordingly
        add = floor(intervals(i,1)*samplingRateLFP_out)/samplingRateLFP_out; % in practice, the interval starts at the nearest lfp timestamp
        structure(i).timestamps = structure(i).timestamps + add;
        structure(i).timestamps = structure(i).timestamps - 1/samplingRateLFP_out; % when using intervals the lfp actually starts 0s away from the first available sample
    end
    
    % Get regions from session or anatomical_map
    if isfield(session,'brainRegions')
        [anatomical_map,channel_map] = get_maps(session);
        anatomical_map = get_map_from_session(session,anatomical_map,channel_map);
        info.region = anatomical_map';
    elseif isfile(fullfile(basepath,'anatomical_map.csv'))
        [anatomical_map,~] = get_maps(session);
        [anatomical_map,~] = get_anatomical_map_csv(basepath,anatomical_map);
        info.region = anatomical_map';
    end
    
end

% save the data into a single matrix (concatenate individual intervals)
lfp = [cat(1,structure.timestamps) double(cat(1,structure.data))];
info.intervals = cat(1,structure.interval);
info.channels = channels;
info.samplingRate = samplingRateLFP_out;
end


function anatomical_map = get_map_from_session(session,anatomical_map,channel_map)
if isfield(session,'brainRegions')
    regions = fields(session.brainRegions);
    for i = 1:length(regions)
        region_idx = ismember(channel_map,...
            session.brainRegions.(regions{i}).channels);
    end
end
end

function [anatomical_map,channel_map] = get_maps(session)
max_channels = max(cellfun('length',session.extracellular.electrodeGroups.channels));
anatomical_map = cell(max_channels,session.extracellular.nElectrodeGroups);
channel_map = nan(size(anatomical_map));

for i = 1:session.extracellular.nElectrodeGroups
    n_ch = length(session.extracellular.electrodeGroups.channels{i});
    anatomical_map(1:n_ch,i) = repmat({'Unknown'},1,n_ch);
    channel_map(1:n_ch,i) = session.extracellular.electrodeGroups.channels{i};
end
end

function [anatomical_map,pull_from_session] = get_anatomical_map_csv(basepath,anatomical_map)
pull_from_session = false;
filename = fullfile(basepath,'anatomical_map.csv');
if ~exist(filename,'file')
    warning('no .anatomical_map.csv... ')
    disp('will try to pull from basename.session')
    disp('you can check anatomical_map and rerun')
    
    writetable(cell2table(anatomical_map),...
        fullfile(basepath,'anatomical_map.csv'),...
        'WriteVariableNames',0)
    pull_from_session = true;
    return
end
anatomical_map = table2cell(readtable(filename,'ReadVariableNames',false));
end





