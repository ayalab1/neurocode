function [spikeT] = importSpikes(varargin)
% [spikes] = importSpikes(varargin)
% Import into the workspace spiketimes for selected units taking as input a
% pre-computed buzcode spikes.cellinfo structure 

% INPUTS
%    basepath        -path to recording where spikes.cellinfo is. Default: pwd
%    UID             -vector subset of UID's to load. Default: all.
%    channel         -specific channels (neuroscope indexing) to use. Default: all

% These inputs requiere CellExplorer cell_metrics pre calculated. NOT IMPLEMENTED YET
%    brainRegion     -string region ID to load neurons from specific region
%    cellType        -cell type to load
%    sleepState      -string sleep state to keep spikes falling within the
%                       specified sleep state

% OUTPUT
%    spikeT - cellinfo struct with the following fields
%          .UID            -unique identifier for each neuron in a recording
%          .times          -cell array of timestamps (seconds) for each neuron
%          .sessionName             

% TODO:
% Add funcionality to load specific cell types or regions
% Do we want more output fields?

% AntonioFR, 8/20; Lindsay K, 11/21

%% Parse inputs 

p = inputParser;
addParameter(p,'basepath',pwd,@isstr);
addParameter(p,'brainRegion','',@isstr); 
addParameter(p,'cellType','',@isstr); 
addParameter(p,'sleepState','',@isstr);
addParameter(p,'UID',[],@isvector);
addParameter(p,'spikes',[],@isstruct);  
addParameter(p,'session',[],@isstruct);  
addParameter(p,'channel',[],@isnumeric);

parse(p,varargin{:})
basepath = p.Results.basepath;
region = p.Results.brainRegion;
type = p.Results.cellType;
state = p.Results.sleepState;
UID = p.Results.UID;
spikes = p.Results.spikes;
session = p.Results.session;
channel = p.Results.channel;

%% Load spikes 
basename = basenameFromBasepath(basepath);

if isempty(spikes) && exist(fullfile(basepath,[basename,'.spikes.cellinfo.mat'])) 
    load(fullfile(basepath,[basename,'.spikes.cellinfo.mat']))
end

%% Output structure

spikeT.UID = spikes.UID;
spikeT.times = spikes.times;

if ~isempty(UID)
   spikeT.UID = spikes.UID(UID);
   spikeT.times = spikes.times(UID);
elseif ~isempty(channel)
    setUn = [];
    for i = 1:length(channel)
        setUn = [setUn find(spikes.maxWaveformCh == channel(i))];
    end
    setUn = sort(setUn);
    spikeT.UID = spikes.UID(setUn);
    spikeT.times = spikes.times(setUn);
end

if ~isempty(region)
    load(fullfile(basepath,[basename,'.cell_metrics.cellinfo.mat']));
    tempUID = []; tempTimes = [];
    tempUID = spikes.UID(strcmp(cell_metrics.brainRegion, region));
    tempTimes = spikes.times(strcmp(cell_metrics.brainRegion, region));
    [~,~,useInd] = intersect(tempUID, spikeT.UID);
    spikeT.UID = spikeT.UID(useInd);
    spikeT.times = spikeT.times(useInd);
end

if ~isempty(type)
    load(fullfile(basepath,[basename,'.cell_metrics.cellinfo.mat']));
    tempUID = []; tempTimes = [];
    tempUID = spikes.UID(strcmp(cell_metrics.putativeCellType, type));
    tempTimes = spikes.times(strcmp(cell_metrics.putativeCellType, type));
    [~,~,useInd] = intersect(tempUID, spikeT.UID);
    spikeT.UID = spikeT.UID(useInd);
    spikeT.times = spikeT.times(useInd);
end

if ~isempty(state)
    load(fullfile(basepath,[basename,'.SleepState.states.mat']));
    if isfield(SleepState.ints, state)
        intervals = eval(strcat('SleepState.ints.',state));
        tempUID = spikeT.UID; tempTimes = spikeT.times; keepct = 1;
        clear spikeT.UID spikeT.times
        for i = 1:length(tempUID)
            tempSpk = [];
            tempSpk = Restrict(tempTimes{i},intervals);
            if ~isempty(tempSpk)
                spikeT.times{keepct} = tempSpk;
                spikeT.UID(keepct) = tempUID(i);
                keepct = keepct+1;
            end
        end
    else
        disp(fieldnames(SleepState.ints));
        error('Incorrect sleep state entered. Please choose from the above instead');
    end
end        
end