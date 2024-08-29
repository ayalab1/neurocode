function auto_theta_cycles(varargin)
% auto_theta_cycles: automatically select theta channel and calc cycles
%
% auto_theta_cycles: selects deep ca1 channels and locates the channel that
% maximizes theta power (pow 4-12 / pow 1-nyquist). Then it finds theta
% cycles using FindThetaCycles.m. It will save basename.thetacycles.events.mat
%
%
%  USAGE
%
%    auto_theta_cycles(varargin)
%
%  INPUT
% =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%    <options>      optional list of property-value pairs (see table below)
%     'basepath'  basepath to session (can be cell array of many) (default = pwd)
%     'passband'  frequency of theta band (default = [4,12])
%     'maximize_theta_power' whether to find the channel that maximizes theta power (default=true)
%     'run_parallel'  to run multiple in parallel (default = false)
%     'overwrite'  to run overwrite existing file (default = false)
% =========================================================================
%
%  OUTPUT
%       basename.thetacycles.events.mat
%
%   Dependencies - basenameFromBasepath,getLFP,FindThetaCycles
%   - assumes
%       *accurate deep/sup classification '.deepSuperficialfromRipple.channelinfo.mat'
%       *channel brain region labels '.session.mat'
%
% Ryan Harvey 2022
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

p = inputParser;
addParameter(p,'basepath',pwd)
addParameter(p,'passband',[4,12])
addParameter(p,'maximize_theta_power',true)
addParameter(p,'run_parallel',false)
addParameter(p,'overwrite',false)

parse(p,varargin{:})
basepath = p.Results.basepath;
passband = p.Results.passband;
maximize_theta_power = p.Results.maximize_theta_power;
run_parallel = p.Results.run_parallel;
overwrite = p.Results.overwrite;

if ~iscell(basepath)
    basepath = {basepath};
end

% iterate over basepaths
if run_parallel && length(basepath)>1
    parfor i = 1:length(basepath)
        run(basepath{i},passband,maximize_theta_power,overwrite)
    end
else
    for i = 1:length(basepath)
        run(basepath{i},passband,maximize_theta_power,overwrite)
    end
end
end

function run(basepath,passband,maximize_theta_power,overwrite)
disp(basepath)
basename = basenameFromBasepath(basepath);

% pass if file already exists
if exist(fullfile(basepath,[basename,'.thetacycles.events.mat']),'file') && ~overwrite
    return
end

% find deep ca1 lfp channel
[lfp,channel] = get_deep_ca1_lfp(basepath,passband,maximize_theta_power);
if isempty(lfp)
    disp('no ca1 lfp')
    return
end

% find theta cycles, pass cleaned signal in using clean lfp
[peaktopeak, troughs, amplitude] = FindThetaCycles(...
    CleanLFP([lfp.timestamps,double(lfp.data)],'thresholds',[8 Inf])...
    );

% package output
thetacycles.timestamps = peaktopeak;
thetacycles.peaks = troughs;
thetacycles.amplitude = amplitude;
thetacycles.amplitudeUnits = 'z-units';
thetacycles.eventID = [];
thetacycles.eventIDlabels = [];
thetacycles.center = median(peaktopeak,2);
thetacycles.duration = peaktopeak(:,2) - peaktopeak(:,1);
thetacycles.detectorinfo.method = 'auto_theta_cycles';
thetacycles.detectorinfo.theta_channel = channel;

% save to basepath
save(fullfile(basepath,[basename,'.thetacycles.events.mat']),'thetacycles')
end

function [lfp,channel] = get_deep_ca1_lfp(basepath,passband,maximize_theta_power)
% get_deep_ca1_lfp: locates a deep ca1 channel that maximizes theta power

lfp = [];
channel=[];

basename = basenameFromBasepath(basepath);

load(fullfile(basepath,[basename,'.session.mat']))

if ~exist(fullfile(basepath,[basename,'.deepSuperficialfromRipple.channelinfo.mat']),'file')
    classification_DeepSuperficial(session,'basepath',basepath,'sample_ripples',true);
end
load(fullfile(basepath,[basename,'.deepSuperficialfromRipple.channelinfo.mat']))

% find deep ca1 channels to check
try
    deep_channels = deepSuperficialfromRipple.channel(contains(deepSuperficialfromRipple.channelClass,'Deep'));
catch % different spelling
    deep_channels = deepSuperficialfromRipple.channels(contains(deepSuperficialfromRipple.channelClass,'Deep'));
end
if isempty(deep_channels)
    return
end

% check for CA1 channels
deep_channels_temp = identify_good_deep_hpc_channel(session,deep_channels,'CA1');

% if no CA1, check CA2
if isempty(deep_channels_temp)
    deep_channels_temp = identify_good_deep_hpc_channel(session,deep_channels,'CA2');
end
% if no CA1 or CA2, return
if isempty(deep_channels_temp)
    return
else
    deep_channels = deep_channels_temp;
end

% load deep channels
[r,c] = size(deep_channels);
if r>c
    deep_channels = deep_channels';
end

if maximize_theta_power
    % try to load downsampled to same time
    try
        lfp = getLFP(deep_channels,'basepath',basepath,'downsample',10,...
            'basename',basename);

        % if sample rate cannot be factored by 10, load entire file
    catch
        lfp = getLFP(deep_channels,'basepath',basepath,...
            'basename',basename);
    end

    % get theta power to choose channel
    try
        pBand = bandpower(single(lfp.data),...
            lfp.samplingRate,passband);

        pTot = bandpower(single(lfp.data),...
            lfp.samplingRate,...
            [1,(lfp.samplingRate/2)-1]);
    catch
        for c = 1:size(lfp,2)
            pBand(c) = bandpower(single(lfp.data(:,c)),...
                lfp.samplingRate,passband);

            pTot(c) = bandpower(single(lfp.data(:,c)),...
                lfp.samplingRate,...
                [1,(lfp.samplingRate/2)-1]);
        end
    end
    % find max theta power, normalized by wide band
    [~,c_idx] = max(pBand./pTot);

    % only leave theta channel
    lfp = getLFP(lfp.channels(:,c_idx),'basepath',basepath,...
        'basename',basename);

else
    lfp = getLFP(randsample(deep_channels,1),'basepath',basepath,...
        'basename',basename);
end
channel = lfp.channels;
end

function deep_channels = identify_good_deep_hpc_channel(session,deep_channels,region)

% search for ca1 channels
brainRegions = fields(session.brainRegions);
ca1_layers = brainRegions(contains(brainRegions,region));

if isempty(ca1_layers)
    deep_channels = [];
    return
end
ca1_channels = [];
for region = ca1_layers'
    ca1_channels = [ca1_channels;session.brainRegions.(region{1}).channels'];
end
ca1_channels = unique(ca1_channels);

% find channels that are deep and ca1
deep_channels = deep_channels(ismember(deep_channels,ca1_channels))';

% kick out bad channels
try
    deep_channels = deep_channels(~ismember(deep_channels,session.channelTags.Bad.channels));
catch
end

if isempty(deep_channels)
    deep_channels = [];
    return
end
end