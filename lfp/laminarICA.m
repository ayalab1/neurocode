function [ica] = laminarICA(varargin)
%    [ica] = doICA(varargin)
%
% Performs Independent Component Analysis (ICA) decomposition

%performs ICA of a laminar LFP profile using the logistic infomax ICA algorithm of Bell & Sejnowski (1995) with
% the natural gradient feature of Amari, Cichocki & Yang.
% The resulting indenpendent components (ICs)correspond to different
% synaptic sources of the LFPs (e.g. CA3 and entorhinal inputs for CA1
% LFPs). For more details about the method implementation see Fernandez-Ruiz,JNeurosci, 2012 and
% Schomburg et al., Neuron, 2014
%
% This function is still work in progress
%
%
% INPUTS
% <optional>
% basepath      Default pwd
% lfp           A matrix of the data organized by [timestamps x Nchannels]
%                   first output of getLFP.m
% region_tag    String input to inidcate a brain region (i.e. CA1 for all
%               CA1 channels). When indexing channels, contains is used so if specific
%               channels are desired be specific (e.g. CA1sp for pyramidal or CA1sr for
%               radiatum etc).
% chinfo        A structure containing
%               If not provided, runs getLFP('all') on basepath
%               IMPORTANT: lfp provided must be in the correct order of
%               channels! If not provided, this function will by default
%               load the channels in the correct order from the first
%               shank.
% shankNum      Specifies which shank to load channels from if lfp is not provided.
%               Default 1
% brain_state   Fieldname from SleepState.states for specified state (i.e.
%               'THETA' for all theta states or 'REM' for all rapid eye
%               movement states. Will load lfp using intervals form
%               SleepState.ints.(brain_state).
% passband      Prefiltering passband interval, default [30 200]
% nICs          Number of independent components to extract. Default 8.
% saveMat       Save results, default true.
% force         Force analysis (disable loading option if already computed,
%                   default false)
% plotWeights   Default true
% plotCFC       Will also calculate theta-gamma CFC for each IC. Default true.
% regionChan    Optional, if a .hippocampalLayers.channelinfo.mat file
%               exists, or if input is provided, it will label these
%               channels on the ica output. Default order is or, pyr, rad,
%               slm
% chanRange     Optional. Index of channels to specify whether you only want a subset of
%               channels from the shank to be used.
%
% TO DO: INCLUDE IMPORTANTS IMPUTS TO runica.m as additional arguments!
%
% OUTPUT
% ica           a buzcode structure with the following fields:
% .activations  independent components (this is the data matrix to use for
%               analysis)
% .data         independent components organized by explanatory variance.
% .timestamps
% .topo         Weight matrix
% .sphere       data sphering matrix (chans,chans)
% .weights      ICA weight matrix (comps,chans)
% .meanvar      Explained variance.
%
% Copyright (C) 2021 by Antonio FR. (using previous code from Manu V, Ipshita Z)
% This functions is a  wrapper for the EEGLAB function 'runica.mat' from Scott
% Makeig (CNL/The Salk Institute, La Jolla, 1996-). Copyright (C) 2004-2011
%
% Updated 10/2022 by Laura Berkowitz
%
% To do:
%       - Add option for applying ICA for each shank if desired.
%       - Remove potentially depreceiated aspects of code see comments for
%         details.
%       - Add option to combine epochs e.g. Theta states during Wake states
%         etc.
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

%% Parse inputs
p = inputParser;
addParameter(p,'basepath',pwd,@ischar);
addParameter(p,'lfp',[],@isstruct);
addParameter(p,'region_tag',[],@iscell);
addParameter(p,'brain_state',[],@ischar);
addParameter(p,'shankNum',1,@isnumeric)
addParameter(p,'saveMat',true,@islogical);
addParameter(p,'force',true,@islogical);
addParameter(p,'passband',[30 200],@isnumeric)
addParameter(p,'nICs',15,@isnumeric)
addParameter(p,'plotWeights',true,@islogical)
addParameter(p,'plotCFC',true,@islogical)
addParameter(p,'regionChan',[],@isnumeric) % possibly depreciated
addParameter(p,'chanRange',[], @isnumeric) % maybe we should not use this and instead use regions (i.e. CA1 for all CA1 channels). Cause channel range could change from day to day or subject to subject.
addParameter(p,'thetaChannel',[],@isnumeric) % This is not used. Reason for this input? -LB 10/21/2022

parse(p,varargin{:});
basepath = p.Results.basepath;
lfp = p.Results.lfp;
region_tag = p.Results.region_tag;
brain_state = p.Results.brain_state;
saveMat = p.Results.saveMat;
force = p.Results.force;
passband = p.Results.passband;
nICs = p.Results.nICs;
shankNum = p.Results.shankNum;
regionChan = p.Results.regionChan;
plotWeights = p.Results.plotWeights;
plotCFC = p.Results.plotCFC;
chanRange = p.Results.chanRange;
thetaChannel = p.Results.thetaChannel;


cd(basepath); % will change this once I figure out which downstream functions require PWD

% check if file already created
if ~isempty(brain_state)
    % check if there exists an ICA file for the indicated brain state
    targetFile = dir([basepath,filesep,'*.ica*',brain_state,'.channelInfo.mat']);
else
    targetFile = dir([basepath,filesep,'*.ica.channelInfo.mat']);
end

% Check if created
if ~isempty(targetFile) && ~force
    disp('ICA already computed! Loading file.');
    load(targetFile.name);
    return
end

% get basename from basepath, and load session
basename = basenameFromBasepath(basepath);

%% Get the channels and regions desired for the ICA.
[channelOrder,region,regionChan] =  select_channels(basepath,chanRange,regionChan,region_tag,shankNum);

%% preprocess lfp for ica (load lfp and filter with passband)
lfp =  preprocess_lfp(channelOrder,basepath,passband,brain_state,lfp);
lfp.region = region;
lfp.timestamps = lfp.timestamps';

% Run runica (saves file to basepath)
ica =  main(lfp,nICs);
clear lfp

if saveMat
    save_ica(ica,basepath,brain_state,region_tag);
end


%% Plotting below
% Plot the components
if plotWeights
    % legacy plot
    plot_weights(ica,regionChan,basepath,nICs)
    
    % plots individual weights
    plot_weights_individual(ica,basepath,nICs)
    
end

%% cross frequency coupling
if plotCFC
    if exist('hippocampalLayers','var')
        channel = hippocampalLayers.pyramidal;
        % search for region of interest in session.
    elseif any(contains(ica.region,'CA1slm'))
        channel = ica.channels(contains(ica.region,'CA1slm'));
        channel = channel(floor(length(channel)/2)); % choose middle channel
        % Else search for tag in anatomical_map
    else
        %Pick the middle channel from the shank for LFP
        channel = channelOrder(floor(length(channelOrder)/2));
    end
    
    % creates plot and saves
    plot_CFC(ica,channel,ica.channels,basepath)
    
end

end

% main functions
function lfp =  preprocess_lfp(channelOrder,basepath,passband,brain_state,lfp)
basename = basenameFromBasepath(basepath);

% grab interval for epoching if indicated
if ~isempty(brain_state)
    % load brain states
    load(fullfile(basepath,[basename,'.SleepState.states.mat']))
    % get intervals of when target states occur
    intervals = SleepState.ints.(brain_state);
else
    intervals = [0 Inf]; % default for getLFP
end

% Load lfp
if isempty(lfp)
    lfp = load_lfp(basepath,channelOrder);
    % index epochs of interest
    [in_idx,~,~] = InIntervals(lfp.timestamps,intervals);
    lfp.data =  lfp.data(in_idx,:);
    lfp.timestamps =  lfp.timestamps(in_idx');
    lfp.in_idx = in_idx; % save this so plotting examples can be indexed using this reference
    lfp.intervals = intervals;
    lfp.channels = channelOrder;
end

% It is recommmend to filter LFPs in the band of interest (e.g. gamma)
disp('Filtering...');
lfp.data = bz_Filter(lfp.data,'passband',passband,'order',4,'filter','butter');
end

function [channelOrder,region,regionChan] =  select_channels(basepath,chanRange,regionChan,region_tag,shankNum)
% select_channels handles input arguments to 1) select desired
% channels (chanRange, regionChan, shankNum,

% get basename from basepath, and load session
basename = basenameFromBasepath(basepath);
session = loadSession(basepath,basename);

% Get channel map and regions from session
[channelOrder, region] = get_channels(session);

% If shankNum is populated, then use value to index channelOrder and region
if ~isempty(shankNum)
    channelOrder = channelOrder{shankNum};
    region = region{shankNum};
    region = region'; % transpose so same dim as channelOrder
end

% Apply channel range if indicated
if ~isempty(chanRange)
    % if its multiple shanks, channelOrder will be a cell array of channels
    if iscell(channelOrder)
        
        for shank = 1:length(channelOrder)
            if isempty(channelOrder{shank})
                continue
            end
            channelOrder{shank} = channelOrder{shank}(chanRange);
            region{shank} = region{shank}(chanRange);
        end
        
    else
        % extract from double
        channelOrder = channelOrder(chanRange);
        region= region(chanRange);
    end
    
elseif ~isempty(region_tag)
    if iscell(channelOrder)
        
        for shank = 1:length(channelOrder)
            if isempty(channelOrder{shank})
                continue
            end
            channelOrder{shank} = channelOrder{shank}(contains(region{shank},region_tag));
            region{shank} = region{shank}(contains(region{shank},region_tag));
        end
        
    else
        % extract from double
        channelOrder = channelOrder(contains(region,region_tag));
        region= region(contains(region,region_tag));
    end
end


% find region info (Seems depreceiated? We use anatomical_map from channel_mapping.m to update session.brainRegion. Consider removing LB 10/22)
if ~isempty(regionChan)
    regFile = dir([basepath,filesep,'*.hippocampalLayers.channelinfo.mat']);
    if ~isempty(regFile)
        load(regFile.name);
        regionChan = hippocampalLayers.all;
    end
end

end

function ica =  main(lfp,nICs,varargin)
% performs ica on filtered lfp using runica and saves output as structure
% to basepath;
%
% inputs:
% lrate - learning rate <<1 (heuristic)
% saveMat - logical indicator to save ica structure to basepath (default is
%           true)
% brain_state - character that is added to the file name and indicates ICA as a function of
%               brain state;
% output:
% structure containig outputs of ICA

p = inputParser;
addParameter(p,'lrate',1.0000e-03,@isnumeric);

parse(p,varargin{:});
lrate = p.Results.lrate;

% Perform ICA
[weights,sphere,meanvar,bias,signs,lrates,data] = runica(lfp.data','lrate',lrate,'pca',nICs);

% Normalize sphere
sphere = sphere./norm(sphere);

% Generate mixing and unmixing matrix
unmixing = weights*sphere;
if (size(unmixing,1)==size(unmixing,2)) && rank(unmixing)==size(unmixing,1)
    mixing = inv(unmixing);
else
    mixing = pinv(unmixing);
end

activations = unmixing * lfp.data';

% Populate output structure
ica.activations = activations;
ica.data = data';
ica.timestamps = lfp.timestamps;
ica.sphere = sphere;
ica.weights = weights;
ica.meanvar = meanvar;
ica.topo = mixing;
ica.unmixing = unmixing;
ica.samplingRate = lfp.samplingRate;
ica.channels = lfp.channels;
ica.region = lfp.region;
ica.in_idx = lfp.in_idx;

end

% helper functions
function save_ica(ica,basepath,brain_state,region_tag)
basename = basenameFromBasepath(basepath);

    disp('Saving results...');
    if ~isempty(brain_state) | ~isempty(region_tag)
        save([basepath,filesep,basename '.ica_',brain_state,'_',region_tag{:},'.channelInfo.mat'],'ica');
    else
        save([basepath,filesep,basename '.ica.channelInfo.mat'],'ica');
    end

end

function lfp = load_lfp(basepath,channelOrder)

basename = basenameFromBasepath(basepath);
session = loadSession(basepath,basename);

% find lfp file
lfp_path = dir(fullfile(basepath,'*.lfp'));
lfp_path = fullfile(lfp_path.folder,lfp_path.name);

% Load lfp
lfp.data = loadBinary(lfp_path,...
    'frequency',session.extracellular.srLfp,...
    'nchannels',session.extracellular.nChannels,...
    'start',0,...
    'duration',Inf,...
    'channels',channelOrder);

% lfp sample rate
Fs = session.extracellular.srLfp;

% prepare structure for below fuctions
lfp.samplingRate = Fs;
lfp.timestamps = (0:length(lfp.data)-1)/Fs;
lfp.channels = channelOrder; 

end

function [channelOrder, region] = get_channels(session)

% use get maps to find brainRegion and channel_map
[anatomical_map,channel_map] = get_maps(session);
anatomical_map = get_map_from_session(session,anatomical_map,channel_map); % updates regions from session and is same dims as channel_map

% loop through shanks and get regions and channels
for shankNum = 1:size(channel_map,2)
    % get channels from session and region from anatomical_map per
    % shank
    temp_channels  =  session.extracellular.electrodeGroups.channels{1,shankNum};
    bad_idx = ismember(temp_channels,session.channelTags.Bad.channels(:));
    temp_region = anatomical_map(find(~bad_idx),shankNum);
    % remove bad channels
    temp_channels(bad_idx) = [];
    % save
    channelOrder{shankNum} = temp_channels;
    region{shankNum} = temp_region;
end

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

% plotting functions

function plot_weights_individual(ica,basepath,nICs)
fig = figure;
fig.Position = [1638,-107,1505,695];
fig.Color = [1 1 1];
for i = 1:nICs
    subplot(1,nICs,i)
    plot(ica.topo(:,i),repmat([1:numel(ica.channels)]'',1,1),'LineWidth',2);
    hold on;
    plot(ones(numel(ica.channels)),repmat([1:numel(ica.channels)]',1,1),'Color','r')
    if i == 1
        set(gca,'YDir','reverse','ytick',repmat([1:numel(ica.channels)]',1,1),'yticklabel',ica.region)
    else
        set(gca,'YDir','reverse')
    end
    hold on
    ylim([1 size(ica.weights,2)])
end

if ~isfolder('ICA')
    mkdir('ICA')
end
saveas(gcf,[basepath,filesep,'ICA\individual_weights.png']);
end

function plot_weights(ica,regionChan,basepath,nICs)

fig = figure;
fig.Color = [1 1 1];
plot(ica.topo,repmat([1:numel(ica.channels)]',1,nICs),'LineWidth',2);
set(gca,'YDir','reverse')
hold on
legend({num2str(ica.channels)})
ylim([1 size(ica.weights,2)])
if ~isempty(regionChan)
    for ii = 1:length(regionChan)
        idx = find(ica.channels == regionChan(ii));
        if ~isempty(idx)
            line([-300 300],[idx idx],'Color',[0.5 0.5 0.5],'LineStyle','--');
        end
    end
end
if ~isfolder('ICA')
    mkdir('ICA')
end
saveas(gcf,[basepath,filesep,'ICA\Weights.png']);

end

function plot_CFC(ica,pyrCh,channelOrder,basepath,varargin)
% Plot the comodulagram for theta from CA1pyr and IC
%
% Currently chooses theta from REMstate as epochs tend to be shorter and theta is clean. 
%
% To-do: 
%   - create input to choose awake Theta state LB 10/2022
%
% inputs:
%   ica - output of ICA analysis
%   channel - channel of CA1 pyramidal layer
%   basepath - path to session that contains SleepStates.states.mat
%
%  (optional)
%   phaserange: vector of frequency range of interest
%               lowerbound:step:upperbound (default 5:0.5:12)
%   amprange:   vecto of amplitude range (default 30:5:200)
%
%   outputs:
%       saves figure of coupling across ICA components

p = inputParser;
addParameter(p,'phaserange',5:0.5:12,@isnumeric);
addParameter(p,'amprange', 30:5:200,@isnumeric);

parse(p,varargin{:});
phaserange = p.Results.phaserange;
amprange = p.Results.amprange;

basename = basenameFromBasepath(basepath);

% Load sleep state to pull interval with theta epoch
load(fullfile(basepath,[basename,'.SleepState.states.mat']))
interval = SleepState.ints.THETA;

% [lfpTheta,lfpICA] = getLFP(pyrCh,'basepath',basepath,'interval',interval); % Only take 1000 seconds worth of data to keep the computation quick
% load lfp (all lfp)
lfp = load_lfp(basepath,channelOrder);

% Keep pyramidal channel only and index same index as ica
lfp.data = lfp.data(ica.in_idx,find(lfp.channels == pyrCh));
lfp.timestamps = lfp.timestamps(ica.in_idx');
lfp.channels = pyrCh; 

% choose the biggest THETA epoch
[~,idx] = max(SleepState.ints.THETA(:,2) - SleepState.ints.THETA(:,1));
interval = interval(idx,:);

% limit LFP to epoch and find index for data in interval for ICA 
[status,~,~] = InIntervals(lfp.timestamps,interval);
lfp.data =  lfp.data(status,:);
lfp.timestamps =  lfp.timestamps(status');
in = find(ica.timestamps >= interval(1) & ica.timestamps <= interval(2));

% if interval is greater than 10 seconds, limit to 10 seconds 
% if size(lfp.data,1) > lfp.samplingRate*10 
%     %limit lfp 
%     lfp.data =  lfp.data(1:lfp.samplingRate*10,:);
%     lfp.timestamps =  lfp.timestamps(1:lfp.samplingRate*10);
%     %limit to first second of interval
%     in = in(1:lfp.samplingRate*10);
% end

%% For each ICA, run CFC

% Make lfp structure, with first channel as the Pyr lfp, and remaining
% channels as the activations for each ICA components.
lfp.data(:,2:size(ica.weights,1)+1) = ica.activations(:,in')';
lfp.channels = 1:size(lfp.data,2);

% Create the plot
[comodulogram] = CFCPhaseAmp(lfp,phaserange,amprange,'phaseCh',1,'ampCh',2:(size(ica.weights,1)+1));

% save to basepath
set(gcf,'Position',[100 100 1400 300])
saveas(gcf,[basepath,filesep,'ICA\CFC.png']);

end