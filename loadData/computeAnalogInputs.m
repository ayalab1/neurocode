function analogInp = computeAnalogInputs(varargin)
%[computeAnalogInputs] - extracts the analog pulses for given channels

%    =========================================================================
%  USAGE
%
%INPUT
%   [analogCh]         [channel number]
%   (options)
%   [ 'saveMat' ]      (default: prompt)
%   ['fs' ]   DAQ (intan) samlping 
%                       

%    =========================================================================

%OUTPUT
%  analogInp     
%   =========================================================================
% 

% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.
%%  Process inputs
p = inputParser;

addParameter(p,'analogCh', [], @isnumeric)
addParameter(p,'fs',20000,@isnumeric)
addParameter(p,'filename',[],@isstring)
addParameter(p,'saveMat',false,@islogical)

parse(p, varargin{:});
analogCh = p.Results.analogCh;
fs = p.Results.fs;
saveMat = p.Results.saveMat;

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

%% Put together output data structure
try % temporary solution
fs_analog = intaninfo.frequency_parameters.board_adc_sample_rate;
fs_lfp = intaninfo.frequency_parameters.amplifier_sample_rate;
catch
fs_analog = intaninfo.supply_voltage_channels.board_adc_sample_rate;
fs_lfp = intaninfo.supply_voltage_channels.amplifier_sample_rate;
end

% set up downsampling if necessary
if fs_analog ~= fs_lfp
    error('analog sampling rate not equal to lfp sampling rate')
end

analogInp = struct();
analogInp.timestamps = [0 : 1/fs_analog : (num_samples-1)/fs_analog]';
analogInp.data = (data'-6800)*5/(59000-6800); % convert to voltage- values taken from data
analogInp.samplingRate = fs_analog;
analogInp.inputChannels = active_channels;

if saveMat
    %savepath = [fileinfo.folder, 'analogInput.mat'];
    savepath = fileinfo.folder;
    basename = basenameFromBasepath(savepath);
    save(fullfile(savepath, [basename '.analogInput.behavior.mat']),'analogInp','-v7.3')
end
    
end
