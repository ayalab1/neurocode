function [spikes, varargout] = importSpikes(varargin)
% importSpikes: loads spikes and metadata from cell metrics
% note: this is the successor to importSpikes.m
%
% Inputs:
%   basepath: optional, path to recording session
%   brainRegion: optional, string, char, or cell array
%   cellType: optional, string, char, or cell array
%   UID: optional, numeric array
%   state: optional, string or char, sleep state to restrict to
%
% Outputs:
%   spikes
%       times: spike times for each neuron
%       UID: unique unit ID
%       brainRegion: brain region of neuron
%       cellType: cell type of neuron
%       n_spikes: number of spikes for each neuron
%       support: time support of spikes
%  varargout:
%   spikeArray: SpikeArray object, optional output
%
%
% Examples:
% spikes = importSpikes('basepath', 'Z:\Data\AYAold\AYA7\day19', ...
%     'brainRegion', 'CA1', 'cellType', 'Pyr')
%
% spikes = importSpikes('basepath', 'Z:\Data\AYAold\AYA7\day19', ...
%     'brainRegion', {'CA1', 'CA3'}, 'cellType', 'Pyr')
%
% spikes = importSpikes('basepath', 'Z:\Data\AYAold\AYA7\day19', ...
%     'brainRegion', {'CA1', 'CA3'}, 'cellType', {'Pyr', 'Int'})
%
% spikes = importSpikes('basepath', 'Z:\Data\AYAold\AYA7\day19', ...
%     'brainRegion', {'CA1', 'CA3'}, 'cellType', {'Pyr', 'Int'}, ...
%     'state', 'NREMstate')
%
% spikes = importSpikes('basepath', 'Z:\Data\AYAold\AYA7\day19', ...
%     'UID', [1, 2, 3])
%
% [spikes, spike_array] = importSpikes('basepath', 'Z:\Data\AYAold\AYA7\day19', ...
%     'brainRegion', 'CA1', 'cellType', 'Pyr')
% spike_array
% <SpikeArray: 15 units> of length 06:23:07.270 hours
%
% Ryan H

% parse inputs
p = inputParser;
addParameter(p, 'basepath', pwd, @isfolder);
addParameter(p, 'brainRegion', '', @(x) ischar(x) || isstring(x) || iscell(x));
addParameter(p, 'cellType', '', @(x) ischar(x) || isstring(x) || iscell(x));
addParameter(p, 'UID', [], @(x) isnumeric(x));
addParameter(p, 'state', '', @(x) ischar(x) || isstring(x));

parse(p, varargin{:})
basepath = p.Results.basepath;
brainRegion = p.Results.brainRegion;
cellType = p.Results.cellType;
UID = p.Results.UID;
state = p.Results.state;

% load cell metrics
basename = basenameFromBasepath(basepath);
if ~exist([basepath, filesep, basename, '.cell_metrics.cellinfo.mat'], 'file')
    warning('no cell_metrics file')
    return
end
load([basepath, filesep, basename, '.cell_metrics.cellinfo.mat'], 'cell_metrics');

% create spikes struct
spikes = cell_metrics.spikes;
spikes.UID = cell_metrics.UID;
spikes.brainRegion = cell_metrics.brainRegion;
spikes.cellType = cell_metrics.putativeCellType;

% add time support by getting the min and max times
spikes.support = [min(cellfun(@min, spikes.times)), max(cellfun(@max, spikes.times))];

% restrict to brainRegion
if ~isempty(brainRegion)
    valid_uid = spikes.UID(contains(spikes.brainRegion, brainRegion));
    keep_idx = ismember(spikes.UID, valid_uid);
    spikes = restrict_cells(spikes, keep_idx);
end

% restrict to cellType
if ~isempty(cellType)
    valid_uid = spikes.UID(contains(spikes.cellType, cellType));
    keep_idx = ismember(spikes.UID, valid_uid);
    spikes = restrict_cells(spikes, keep_idx);
end

% restrict to UID
if ~isempty(UID)
    keep_idx = ismember(spikes.UID, UID);
    spikes = restrict_cells(spikes, keep_idx);
end

% remove cells labeled bad
if isfield(cell_metrics, 'tags')
    if isfield(cell_metrics.tags, 'Bad')
        keep_idx = ~ismember(spikes.UID, cell_metrics.tags.Bad);
        spikes = restrict_cells(spikes, keep_idx);
    end
end

% restrict to state
if ~isempty(state)
    load(fullfile(basepath, [basename, '.SleepState.states.mat']), 'SleepState');
    if any(contains(state, fields(SleepState.ints)))
        spikes.times = cellfun(@(x) Restrict(x, SleepState.ints.(state)), ...
            spikes.times, 'UniformOutput', false);
    else
        disp(fieldnames(SleepState.ints));
        error('Incorrect sleep state entered. Please choose from the above instead');
    end
end

% number of spikes per cell
spikes.n_spikes = cellfun(@length, spikes.times);

% verify outputs lengths
n_cells = length(spikes.times);
assert(length(spikes.UID) == n_cells)
assert(length(spikes.brainRegion) == n_cells)
assert(length(spikes.cellType) == n_cells)
assert(length(spikes.n_spikes) == n_cells)

% output SpikeArray if requested, otherwise just return spikes struct
if nargout == 2
    if isempty(spikes.times)
        varargout{1} = SpikeArray();
    else
        varargout{1} = SpikeArray(spikes.times);
    end
end
end

function spikes = restrict_cells(spikes, keep_idx)
% restricts spikes struct to keep_idx
spikes.times = spikes.times(keep_idx);
spikes.UID = spikes.UID(keep_idx);
spikes.brainRegion = spikes.brainRegion(keep_idx);
spikes.cellType = spikes.cellType(keep_idx);
end