function [lfp] = getLFP(varargin)
%[ getLFP] - Get local field potentials.
%
%  Load local field potentials from disk. No longer dependent on
%  FMAT/SetCurrentSession.
%
%    =========================================================================
%  USAGE
%
%    [lfp] = getLFP(channels,<options>)
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
%    =========================================================================
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

%    =========================================================================
% Copyright (C) 2004-2011 by Michaël Zugaro
% editied by David Tingley, 2017
% kathryn mcclain 2020
% Ralitsa Todorova 2022
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
%% Parse the inputs!

channelsValidation = @(x) isnumeric(x) || strcmp(x,'all');

% parse args
p = inputParser;
addRequired(p,'channels',channelsValidation)
addParameter(p,'basename','',@isstr)
addParameter(p,'intervals',[],@isnumeric)
addParameter(p,'restrict',[],@isnumeric)
addParameter(p,'basepath',pwd,@isstr);
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
   end
   lfp.Filename = d.name;
   basename = strsplit(lfp.Filename,'.');
   if length(basename) > 2
       base = [];
       for i=1:length(basename)-1
          base = [base basename{i} '.'];
       end
       basename = base(1:end-1);  % this is an fugly hack to make things work with Kenji's naming system...
   else
       basename = basename{1};
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
   lfp.Filename = d.name;   
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
    display(['Loading Channels ',num2str(channels),' (1-indexing)'])
end

%% get the data
disp('loading LFP file...')
nIntervals = size(intervals,1);
% returns lfp/bz format
for i = 1:nIntervals
    lfp(i).duration = (intervals(i,2)-intervals(i,1));
    lfp(i).interval = [intervals(i,1) intervals(i,2)];

    % Load data and put into struct
    % we assume 0-indexing like neuroscope, but loadBinary uses 1-indexing to
    % load....
    lfp(i).data = loadBinary([basepath filesep lfp.Filename],...
        'duration',double(lfp(i).duration),...
                  'frequency',samplingRate,'nchannels',session.extracellular.nChannels,...
                  'start',double(lfp(i).interval(1)),'channels',channels,...
                  'downsample',downsamplefactor);
    lfp(i).timestamps = (1:length(lfp(i).data))'/samplingRateLFP_out;
    lfp(i).channels = channels;
    lfp(i).samplingRate = samplingRateLFP_out;
    % check if duration is inf, and reset to actual duration...
    if lfp(i).interval(2) == inf
        lfp(i).interval(2) = length(lfp(i).timestamps)/lfp(i).samplingRate;
        lfp(i).duration = (lfp(i).interval(i,2)-lfp(i).interval(i,1));
    end
    if lfp(i).interval(1)>0 % shift the timestamps accordingly
        add = floor(intervals(i,1)*samplingRateLFP_out)/samplingRateLFP_out; % in practice, the interval starts at the nearest lfp timestamp
        lfp(i).timestamps = lfp(i).timestamps + add; 
        lfp(i).timestamps = lfp(i).timestamps - 1/samplingRateLFP_out; % when using intervals the lfp actually starts 0s away from the first available sample
    end
    if isfield(session,'brainRegions') && isfield(session,'channels')
        [~,~,regionidx] = intersect(lfp(i).channels,session.channels,'stable');
        lfp(i).region = session.brainRegions(regionidx); % match region order to channel order..
    end
end
