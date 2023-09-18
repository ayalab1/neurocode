function analogInp = computeAnalogInputs(varargin)
%[computeAnalogInputs] - extracts the analog pulses for given channels

%    =========================================================================
%  USAGE
%
%INPUT
%   [analogCh]         [channel number. Before 1-indexing]
%                       
%=========================================================================
%     Properties    Values
%    ------------------------------------------------------------------------
%     ['savMat']       [option to save as .mat output. Default false]
%     ['downsample']   [Whether to downsample analogin. Default false]
%
%     ['downsampfreq'] [If downsample is true, then what frequency 
%                       to downsample to. Default is 1250]
%     ['intervals']    [Intervals to extract in form of [start stop]
%                           Default is [0 inf]]
%     ['restrict']     [Similar for intervals, for FMA preference of
%                           Restrict]
%
%  OUTPUT
%   [Creates file:   basePath/baseName.analogin]
%
%   [Hlearsson] [2023-2024] -[added intervals, downsampling equal to getLFP
%                               and output structure equal to getLFP]
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.
%%  Process inputs
p = inputParser;

addParameter(p,'analogCh', [], @isnumeric)
addParameter(p,'filename',[],@isstring)
addParameter(p,'saveMat',false,@islogical)
addParameter(p, 'downsample', false, @logical)
addParameter(p,'downsampfreq', 1250, @isnumeric)
addParameter(p,'intervals',[],@isnumeric)
addParameter(p,'restrict',[],@isnumeric)

parse(p, varargin{:});
analogCh = p.Results.analogCh;
saveMat = p.Results.saveMat;
downsample = p.Results.downsample;
downsampfreq = p.Results.downsampfreq;

% Checking analog input file
filename = p.Results.filename;
if isempty(filename)
    filename = 'analogin.dat';
end
if ~exist(filename,'file')
    error('Cannot find analog input file');
end
fileinfo = dir(filename);
if fileinfo.bytes == 0
    error('Analog input file empty');
end
% doing this so you can use either 'intervals' or 'restrict' as parameters to do the same thing
intervals = p.Results.intervals;
restrict = p.Results.restrict;
if isempty(intervals) && isempty(restrict) % both empty
    intervals = [0 Inf];
elseif isempty(intervals) && ~isempty(restrict) % intervals empty, restrict isn't
    intervals = restrict;
end
%% Get info from info.rhd file

% Finding file
infofile = ls('*.rhd');
if isempty(infofile) % look in subdirectories
    d = dir(pwd);
    for i = 1:length(d)
        subdir = d(i).name;
        if isfolder(subdir)
            subfile = ls(strcat(subdir,'\*.rhd'));
            if ~isempty(subfile)
                infofile = strcat(subdir,filesep,subfile);
                break
            end
        end
    end
end
if isempty(infofile)
    error('Cannot find info.rhd file')
end

% Getting relevant info
intaninfo = read_Intan_Info_Wrapper(infofile);
n_active_channels = length(intaninfo.board_adc_channels);
for i = 1:n_active_channels
    active_channels(i) = intaninfo.board_adc_channels(i).native_order;
end
if isempty(analogCh)
    wantedInds = 1:length(active_channels);
elseif sum(ismember(analogCh,active_channels))<length(analogCh)
    error('Analog channel requested not active during recording')
else
    wantedInds = zeros(length(analogCh),1);
    for i = 1:length(analogCh)
        wantedInds(i) = find(active_channels == analogCh(i));
    end
end

%% Read the analog data
num_samples = fileinfo.bytes/(n_active_channels*2);
fid = fopen(filename,'r');
data = fread(fid, [n_active_channels, num_samples], 'uint16');
fclose(fid);
data = data(wantedInds,:);


%% Downsampling - note this downsampling is using loadBinary to ensure compatibility with lfp - Hlarsson 2023
fs_analog = intaninfo.frequency_parameters.board_adc_sample_rate;
if downsample
    nIntervals = size(intervals,1);
    disp('loading Analogin file and downsampling...');
    downsampleFactor = round(fs_analog/downsampfreq);
    for i = 1:nIntervals
        analogInp = struct();
        analogInp.Filename = 'analogin.dat';
        analogInp(i).duration = (intervals(i,2)-intervals(i,1));
        analogInp(i).interval = [intervals(i,1) intervals(i,2)];
        
        % Load data and put into struct
        % we assume 0-indexing like neuroscope, but loadBinary uses 1-indexing to
        % load....
        active_channels = active_channels +1; %binarize it
        analogInp(i).data = loadBinary('analogin.dat',...
            'duration',double(analogInp(i).duration),...
            'frequency',fs_analog,'nchannels',n_active_channels,...
            'start',double(analogInp(i).interval(1)),'channels',active_channels,...
            'downsample',downsampleFactor);
        analogInp.data = (data'-6800)*5/(59000-6800);
        analogInp(i).timestamps = (1:length(analogInp(i).data))'/downsampfreq;
        analogInp(i).channels = active_channels;
        analogInp(i).samplingRate = downsampfreq;
        % check if duration is inf, and reset to actual duration...
        if analogInp(i).interval(2) == inf
            analogInp(i).duration = [];
            analogInp(i).interval(2) = length(analogInp(i).timestamps)/analogInp(i).samplingRate;
            analogInp(i).duration = (analogInp(i).interval(i,2)-analogInp(i).interval(i,1));
        end
        if analogInp(i).interval(1)>0 % shift the timestamps accordingly
            % in practice, the interval starts at the nearest lfp timestamp
            add = floor(intervals(i,1)*downsampfreq)/downsampfreq; 
            analogInp(i).timestamps = analogInp(i).timestamps + add;
            % when using intervals the lfp actually starts 0s away from the first available sample
            analogInp(i).timestamps = analogInp(i).timestamps - 1/downsampfreq; 
        end
        
        % Assign region as analog to be able to integrate with metadata
        % neccessities
        analogInp(i).region ='analog';
    
    end
end

%% Output
% set up downsampling if necessary - added above downsampling HLarsson 2023

if ~downsample
        nIntervals = size(intervals,1);
    disp('loading Analogin file and downsampling...');
    downsampleFactor = 1;
    for i = 1:nIntervals
        analogInp = struct();
        analogInp.Filename = 'analogin.dat';
        analogInp(i).duration = (intervals(i,2)-intervals(i,1));
        analogInp(i).interval = [intervals(i,1) intervals(i,2)];
        
        % Load data and put into struct
        % we assume 0-indexing like neuroscope, but loadBinary uses 1-indexing to
        % load....
        active_channels = active_channels +1; %binarize it
        analogInp(i).data = loadBinary('analogin.dat',...
            'duration',double(analogInp(i).duration),...
            'frequency',fs_analog,'nchannels',n_active_channels,...
            'start',double(analogInp(i).interval(1)),'channels',active_channels,...
            'downsample',downsampleFactor);
        analogInp.data = (data'-6800)*5/(59000-6800);
        analogInp(i).timestamps = (1:length(analogInp(i).data))'/downsampfreq;
        analogInp(i).channels = active_channels;
        analogInp(i).samplingRate = downsampfreq;
        % check if duration is inf, and reset to actual duration...
        if analogInp(i).interval(2) == inf
            analogInp(i).duration = [];
            analogInp(i).interval(2) = length(analogInp(i).timestamps)/analogInp(i).samplingRate;
            analogInp(i).duration = (analogInp(i).interval(i,2)-analogInp(i).interval(i,1));
        end
        if analogInp(i).interval(1)>0 % shift the timestamps accordingly
            % in practice, the interval starts at the nearest lfp timestamp
            add = floor(intervals(i,1)*downsampfreq)/downsampfreq; 
            analogInp(i).timestamps = analogInp(i).timestamps + add;
            % when using intervals the lfp actually starts 0s away from the first available sample
            analogInp(i).timestamps = analogInp(i).timestamps - 1/downsampfreq; 
        end
        
        % Assign region as analog to be able to integrate with metadata
        % neccessities
        analogInp(i).region ='analog';
    end
end

if saveMat
    %savepath = [fileinfo.folder, 'analogInput.mat'];
    savepath = fileinfo.folder;
    basename = basenameFromBasepath(savepath);
    save(fullfile(savepath, [basename '.analogInput.behavior.mat']),'analogInp','-v7.3')
end
    
end
