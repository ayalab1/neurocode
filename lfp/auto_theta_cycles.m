function auto_theta_cycles(varargin)
% auto_theta_cycles: automatically select theta channel and calc cycles
%
% auto_theta_cycles: selects deep ca1 channels and locates the channel that
% maximizes theta power (pow 6-12 / pow 1-nyquist). Then it finds theta
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
%     ['basepath']  basepath to session (can be cell array of many) (default = pwd)
%     ['passband']  frequency of theta band (default = [6,12])
%     ['run_parallel']  to run multiple in parallel (default = false)
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
addParameter(p,'passband',[6,12])
addParameter(p,'run_parallel',false)

parse(p,varargin{:})
basepath = p.Results.basepath;
passband = p.Results.passband;
run_parallel = p.Results.run_parallel;

if ~iscell(basepath)
    basepath = {basepath};
end

% iterate over basepaths
if run_parallel && length(basepath)>1
    parfor i = 1:length(basepath)
        run(basepath{i},passband)
    end
else
    for i = 1:length(basepath)
        run(basepath{i},passband)
    end
end
end

function run(basepath,passband)

basename = basenameFromBasepath(basepath);

% % pass if file already exists
% if exist(fullfile(basepath,[basename,'.thetacycles.events.mat']),'file')
%     return
% end

% find deep ca1 lfp channel 
lfp = get_deep_ca1_lfp(basepath,passband);
if isempty(lfp)
    disp('no ca1 lfp')
    return
end

% find theta cycles
[peaktopeak, troughs] = FindThetaCycles([lfp.timestamps,double(lfp.data)]);

% package output
thetacycles.timestamps = peaktopeak;
thetacycles.peaks = troughs;
thetacycles.amplitude = [];
thetacycles.amplitudeUnits = [];
thetacycles.eventID = [];
thetacycles.eventIDlabels = [];
thetacycles.center = median(peaktopeak,2);
thetacycles.duration = peaktopeak(:,2) - peaktopeak(:,1);
thetacycles.detectorinfo.method = 'auto_theta_cycles';
thetacycles.detectorinfo.theta_channel = lfp.channels;

% save to basepath
save(fullfile(basepath,[basename,'.thetacycles.events.mat']),'thetacycles')
end

function lfp = get_deep_ca1_lfp(basepath,passband)
% get_deep_ca1_lfp: locates a deep ca1 channel that maximizes theta power

basename = basenameFromBasepath(basepath);

load(fullfile(basepath,[basename,'.session.mat']))
load(fullfile(basepath,[basename,'.deepSuperficialfromRipple.channelinfo.mat']))

% find deep ca1 channels to check
try
    deep_channels = deepSuperficialfromRipple.channel(contains(deepSuperficialfromRipple.channelClass,'Deep'));
catch
    deep_channels = deepSuperficialfromRipple.channels(contains(deepSuperficialfromRipple.channelClass,'Deep'));
end
if isempty(deep_channels)
    lfp = [];
    return
end

% search for ca1 channels
brainRegions = fields(session.brainRegions);
ca1_layers = brainRegions(contains(brainRegions,'CA1'));
if isempty(ca1_layers)
    lfp = [];
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
deep_channels = deep_channels(~ismember(deep_channels,session.channelTags.Bad.channels));

if isempty(deep_channels)
    lfp = [];
    return
end

% load deep channels
[r,c] = size(deep_channels);
if r>c
    deep_channels = deep_channels';
end
lfp = getLFP(deep_channels,'basepath',basepath,...
    'basename',basename,'downsample',10);

% get theta power to choose channel
try
    pBand = bandpower(double(lfp.data),...
        lfp.samplingRate,passband);
    
    pTot = bandpower(double(lfp.data),...
        lfp.samplingRate,...
        [1,(lfp.samplingRate/2)-1]);
catch
    for c = 1:size(lfp.data,2)
        pBand(c) = bandpower(double(lfp.data(:,c)),...
            lfp.samplingRate,passband);
        
        pTot(c) = bandpower(double(lfp.data(:,c)),...
            lfp.samplingRate,...
            [1,(lfp.samplingRate/2)-1]);
    end
end
% find max theta power, normalized by wide band
[~,c_idx] = max(pBand./pTot);

% only leave theta channel
lfp = getLFP(lfp.channels(:,c_idx),'basepath',basepath,...
    'basename',basename,'downsample',2);
end