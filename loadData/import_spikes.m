function spikes = import_spikes(varargin)
% import_spikes: loads spikes and metadata from cell metrics
% note: this is the successor to importSpikes.m
%
% Inputs:
%   basepath: optional, path to recording session
%   brainRegion: optional, string, char, or cell array
%   putativeCellType: optional, string, char, or cell array
%
% Outputs:
%   spikes
%       times: spike times for each neuron
%       UID: unique unit ID
%       brainRegion: brain region of neuron
%       putativeCellType: cell type of neuron
%
% Ryan H

p = inputParser;
addParameter(p, 'basepath', pwd, @isfolder);
addParameter(p, 'brainRegion', '', @(x) ischar(x) || isstring(x) || iscell(x));
addParameter(p, 'putativeCellType', '', @(x) ischar(x) || isstring(x));

parse(p, varargin{:})
basepath = p.Results.basepath;
brainRegion = p.Results.brainRegion;
putativeCellType = p.Results.putativeCellType;

basename = basenameFromBasepath(basepath);

load([basepath, filesep, basename, '.cell_metrics.cellinfo.mat'], 'cell_metrics');

spikes = cell_metrics.spikes;
spikes.UID = cell_metrics.UID;
spikes.brainRegion = cell_metrics.brainRegion;
spikes.putativeCellType = cell_metrics.putativeCellType;

% restrict to brainRegion
if ~isempty(brainRegion)
    valid_uid = cell_metrics.UID(contains(cell_metrics.brainRegion, brainRegion));
    keep_idx = ismember(spikes.UID, valid_uid);
    spikes = restrict_cells(spikes, keep_idx);
end

% restrict to putativeCellType
if ~isempty(putativeCellType)
    valid_uid = cell_metrics.UID(contains(cell_metrics.putativeCellType, putativeCellType));
    keep_idx = ismember(spikes.UID, valid_uid);
    spikes = restrict_cells(spikes, keep_idx);
end

% remove cells labeled bad
if isfield(cell_metrics, 'tags')
    if isfield(cell_metrics.tags, 'Bad')
        keep_idx = ~ismember(spikes.UID, cell_metrics.tags.Bad);
        spikes = restrict_cells(spikes, keep_idx);
    end
end

% verify outputs lengths
n_cells = length(spikes.times);
assert(length(spikes.UID) == n_cells)
assert(length(spikes.brainRegion) == n_cells)
assert(length(spikes.putativeCellType) == n_cells)

end

function spikes = restrict_cells(spikes, keep_idx)
spikes.times = spikes.times(keep_idx);
spikes.UID = spikes.UID(keep_idx);
spikes.brainRegion = spikes.brainRegion(keep_idx);
spikes.putativeCellType = spikes.putativeCellType(keep_idx);
end