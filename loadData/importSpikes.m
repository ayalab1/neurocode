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

% OUTPUT
%    spikeT - cellinfo struct with the following fields
%          .UID            -unique identifier for each neuron in a recording
%          .times          -cell array of timestamps (seconds) for each neuron
%          .sessionName             

% TODO:
% Add funcionality to load specific cell types or regions
% Do we want more output fields?

% AntonioFR, 8/20

%% Parse inputs 

p = inputParser;
addParameter(p,'basepath',pwd,@isstr);
addParameter(p,'brainRegion','',@isstr); 
addParameter(p,'cellType','',@isstr); 
addParameter(p,'UID',[],@isvector);
addParameter(p,'spikes',[],@isstruct);  
addParameter(p,'session',[],@isstruct);  
addParameter(p,'channel',[],@isnumeric);

parse(p,varargin{:})
basepath = p.Results.basepath;
region = p.Results.brainRegion;
type = p.Results.cellType;
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

if isempty(channel)
    if ~isempty(region) && ~isempty(type)
       load(fullfile(basepath,[basename,'.cell_metrics.cellinfo.mat']));
       spikeT.UID = []; spikeT.times = [];
       for i = 1:numel(cell_metrics.brainRegion)
           if contains(cell_metrics.brainRegion{i},region) && ...
              strcmp(cell_metrics.putativeCellType{i},type)  
              spikeT.UID = cat(2,spikeT.UID,spikes.UID(i));
              spikeT.times = cat(2,spikeT.times,spikes.times(i));
           end
           spikeT.brainRegion = region;
           spikeT.cellTypes = type;
       end
    elseif ~isempty(region) && isempty(type)
       load(fullfile(basepath,[basename,'.cell_metrics.cellinfo.mat']));
       spikeT.UID = []; spikeT.times = [];
       for i = 1:numel(cell_metrics.brainRegion)
           if contains(cell_metrics.brainRegion{i},region)  
              spikeT.UID = cat(2,spikeT.UID,spikes.UID(i));
              spikeT.times = cat(2,spikeT.times,spikes.times(i));
           end
           spikeT.brainRegion = region;
       end    
    elseif isempty(region) && ~isempty(type)
       load(fullfile(basepath,[basename,'.cell_metrics.cellinfo.mat']));
       spikeT.UID = []; spikeT.times = [];
       for i = 1:numel(cell_metrics.brainRegion)
           if strcmp(cell_metrics.putativeCellType{i},type)  
              spikeT.UID = cat(2,spikeT.UID,spikes.UID(i));
              spikeT.times = cat(2,spikeT.times,spikes.times(i));
           end
            spikeT.cellTypes = type;
       end    
    end
else
    if ~isempty(region) && ~isempty(type)
       load(fullfile(basepath,[basename,'.cell_metrics.cellinfo.mat']));
       spikeT.UID = []; spikeT.times = [];
       for i = 1:numel(cell_metrics.brainRegion)
           if contains(cell_metrics.brainRegion{i},region) && ...
              strcmp(cell_metrics.putativeCellType{i},type)
                if ~isempty(find(setUn==i))
                  spikeT.UID = cat(2,spikeT.UID,spikes.UID(i));
                  spikeT.times = cat(2,spikeT.times,spikes.times(i));
                end
           end
           spikeT.brainRegion = region;
           spikeT.cellTypes = type;
       end
    elseif ~isempty(region) && isempty(type)
       load(fullfile(basepath,[basename,'.cell_metrics.cellinfo.mat']));
       spikeT.UID = []; spikeT.times = [];
       for i = 1:numel(cell_metrics.brainRegion)
           if contains(cell_metrics.brainRegion{i},region)  
               if ~isempty(find(setUn==i))
                  spikeT.UID = cat(2,spikeT.UID,spikes.UID(i));
                  spikeT.times = cat(2,spikeT.times,spikes.times(i));
               end
           end
           spikeT.brainRegion = region;
       end    
    elseif isempty(region) && ~isempty(type)
       load(fullfile(basepath,[basename,'.cell_metrics.cellinfo.mat']));
       spikeT.UID = []; spikeT.times = [];
       for i = 1:numel(cell_metrics.brainRegion)
           if strcmp(cell_metrics.putativeCellType{i},type)  
               if ~isempty(find(setUn==i))
                  spikeT.UID = cat(2,spikeT.UID,spikes.UID(i));
                  spikeT.times = cat(2,spikeT.times,spikes.times(i));
               end
           end
            spikeT.cellTypes = type;
       end    
    end
end

