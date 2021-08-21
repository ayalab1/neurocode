function [spikeT] = importSpikes(varargin)
% [spikes] = importSpikes(varargin)
% Import into the workspace spiketimes for selected units taking as input a
% pre-computed buzcode spikes.cellinfo structure 

% INPUTS
%    basepath        -path to recording where spikes.cellinfo is. Default: pwd
%    UID             -vector subset of UID's to load. Default: all.

% These inputs requiere CellExplorer cell_metrics pre calculated. NOT IMPLEMENTED YET
%    region          -string region ID to load neurons from specific region
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
addParameter(p,'region','',@isstr); 
addParameter(p,'UID',[],@isvector);
addParameter(p,'spikes',[],@isstruct); % Load existing spikes structure 
addParameter(p,'session',[],@isstruct); % A buzsaki lab session struct

parse(p,varargin{:})
basepath = p.Results.basepath;
region = p.Results.region;
UID = p.Results.UID;
spikes = p.Results.spikes;
session = p.Results.session;

%% Load spikes 
cd(basepath);
basename = bz_BasenameFromBasepath(pwd);

if isempty(spikes) && exist(fullfile(basepath,[basename,'.spikes.cellinfo.mat'])) 
    load(fullfile(basepath,[basename,'.spikes.cellinfo.mat']))
end

%% Output structure
if ~isempty(UID)
   spikeT.UID = spikes.UID(UID);
   spikeT.times = spikes.times(UID);
end
   spikeT.sessionName = spikes.sessionName;

end

