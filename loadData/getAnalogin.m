function analogin = getAnalogin(varargin)
% getLFP - Get downsample .analogin file after running analogfromDat.
%
%  Load .analogin from disk. No longer dependent on
%  FMAT/SetCurrentSession.
%
%  Function is intended to load .analogin in .lfp style to enable
%   extraction of analog signals from stimulation (opto/stim fiber) that
%   can be used as lfp type file for lfp type analysis.
%
%  USAGE
%
%    [analogin] = getAnalogin(channels,<options>)
%
%  INPUTS
%
%    channels(required) -must be first input, numeric
%                        list of channels to load (use keyword 'all' for all)
%                        channID is 1-indexing
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
%  OUTPUT
%
%    lfp             struct of lfp data. Can be a single struct or an array
%                    of structs for different intervals.  lfp(1), lfp(2),
%                    etc for intervals(1,:), intervals(2,:), etc
%    .data           [Nt x Nd] matrix of the LFP data
%    .timestamps     [Nt x 1] vector of timestamps to match LFP data
%    .interval       [1 x 2] vector of start/stop times of LFP interval
%    .channels       [Nd X 1] vector of channel ID's
%    .samplingRate   LFP sampling rate [default = 1250]
%    .duration       duration, in seconds, of LFP interval
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

% [HLarsson] [2023-2024]
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.
%
% NOTES
% -'select' option has been removed, it allowed switching between 0 and 1
%   indexing.  This should no longer be necessary with .lfp.mat structs
%
% TODO
% add saveMat input
% expand channel selection options (i.e. region or spikegroup)
% add forcereload
%
% SEE ALSO
% this function was written off of getLFP but altered for .analogin
% 
%% Parse the inputs!

channelsValidation = @(x) isnumeric(x) || strcmp(x,'all');

% parse args
p = inputParser;
addRequired(p,'channels',channelsValidation)
addParameter(p,'basename','',@isstr)
addParameter(p,'intervals',[],@isnumeric)
addParameter(p,'restrict',[],@isnumeric)
addParameter(p,'basepath',pwd,@isfolder);
addParameter(p,'downsample',1,@isnumeric);
addParameter(p,'saveMat',false,@islogical);
addParameter(p,'forceReload',false,@islogical);
addParameter(p,'noPrompts',false,@islogical);
addParameter(p,'fromDat',false,@islogical);

parse(p,varargin{:})
channels = p.Results.channels;
downsamplefactor = p.Results.downsample;
basepath = p.Results.basepath;
basename = p.Results.basename; if isempty(basename), basename = basenameFromBasepath(basepath); end
noPrompts = p.Results.noPrompts;
fromDat = p.Results.fromDat;

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
    disp('No basename given, so we look for a .analogin file...')
    switch fromDat
        case false
            d = dir([basepath filesep '*analogin']);
        case true
            d = dir(['analogin.dat']);
    end
    if length(d) > 1 % we assume one .lfp file or this should break
        error('there is more than one .analogin file in this directory?');
    elseif isempty(d)
        d = dir([basepath filesep '*analogin']);
        if isempty(d)
            error('could not find an .analogin file...')
        end
    end
    analogin.Filename = d.name;
    basename = strsplit(analogin.Filename,'.');
    if length(basename) > 2
        base = [];
        for i=1:length(basename)-1
            base = [base basename{i} '.'];
        end
        basename = base(1:end-1);
    else
        basename = basename{1};
    end
    
else
    switch fromDat
        case false
            d = dir([basepath filesep basename '.analogin']);
        case true
            d = dir('analogin.dat');
    end
    
    if length(d) > 1 % we assume one .lfp file or this should break
        error('there is more than one .analogin file in this directory?');
    elseif isempty(d)
        d = dir([basepath filesep basename '.analogin']);
        if isempty(d)
            error('could not find an .analogin file..')
        end
    end
    analogin.Filename = d.name;
end

%% things we can parse from sessionInfo or xml file

session = getSession('basepath',basepath);

if fromDat
    samplingRate = session.extracellular.sr;
else
    samplingRate = session.extracellular.srLfp;
end

samplingRateLFP_out = samplingRate./downsamplefactor;

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
    disp(['Loading Channels ',num2str(channels),' (1-indexing)'])
end

%% get the data
disp('loading LFP file...')
nIntervals = size(intervals,1);
% returns lfp/bz format
for i = 1:nIntervals
    analogin(i).duration = (intervals(i,2)-intervals(i,1));
    analogin(i).interval = [intervals(i,1) intervals(i,2)];
    
    % Load data and put into struct
    % we assume 0-indexing like neuroscope, but loadBinary uses 1-indexing to
    % load....
    analogin(i).data = loadBinary([basepath filesep analogin.Filename],...
        'duration',double(analogin(i).duration),...
        'frequency',samplingRate,'nchannels',session.extracellular.nChannels,...
        'start',double(analogin(i).interval(1)),'channels',channels,...
        'downsample',downsamplefactor);
    analogin(i).timestamps = (1:length(analogin(i).data))'/samplingRateLFP_out;
    analogin(i).channels = channels;
    analogin(i).samplingRate = samplingRateLFP_out;
    % check if duration is inf, and reset to actual duration...
    if analogin(i).interval(2) == inf
        analogin(i).interval(2) = length(analogin(i).timestamps)/analogin(i).samplingRate;
        analogin(i).duration = (analogin(i).interval(i,2)-analogin(i).interval(i,1));
    end
    if analogin(i).interval(1)>0 % shift the timestamps accordingly
        % in practice, the interval starts at the nearest lfp timestamp
        add = floor(intervals(i,1)*samplingRateLFP_out)/samplingRateLFP_out; 
        analogin(i).timestamps = analogin(i).timestamps + add;
        % when using intervals the lfp actually starts 0s away from the first available sample
        analogin(i).timestamps = analogin(i).timestamps - 1/samplingRateLFP_out; 
    end
    
    % Get regions from session or anatomical_map
    if isfield(session,'brainRegions')
        [anatomical_map,channel_map] = get_maps(session);
        anatomical_map = get_map_from_session(session,anatomical_map,channel_map);
        analogin(i).region = get_region(channels, anatomical_map,channel_map);
    elseif isfile(fullfile(basepath,'anatomical_map.csv'))
        [anatomical_map,channel_map] = get_maps(session);
        [anatomical_map,~] = get_anatomical_map_csv(basepath,anatomical_map);
        analogin(i).region = get_region(channels, anatomical_map,channel_map);
    else
        disp('No brain regions associated with channels found. Saving ''Unkown''')
        [anatomical_map,channel_map] = get_maps(session);
        analogin(i).region = get_region(channels, anatomical_map,channel_map);
    end
end
end

function region = get_region(channels, anatomical_map,channel_map)

% reshape so are both row vectors for comparison with channels
channel_map = channel_map(:)'; % make into row vector
anatomical_map = anatomical_map(:)';
% index region for channels
region = anatomical_map(ismember(channel_map(:)',channels));

end

function anatomical_map = get_map_from_session(session,anatomical_map,channel_map)
if isfield(session,'brainRegions')
    regions = fields(session.brainRegions);
    for i = 1:length(regions)
        region_idx = ismember(channel_map,...
            session.brainRegions.(regions{i}).channels);
        anatomical_map(region_idx) = {regions{i}};
        
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