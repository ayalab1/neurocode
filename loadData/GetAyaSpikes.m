function [spikes,regionID,regionNames,spikesCell,order] = GetAyaSpikes(basepath)

% Loads the spikes of the session in "basepath" in a matrix [timestamp id] format
% This is a simple function loading data saved in CellExplorer format into
% a structure-free matrix and vector format, which FMAToolbox users might find useful.
%
%  USAGE
%
%    [spikes,regionID,regionNames] = GetAyaSpikes(basepath);
%
%    basepath       The full path to the folder where the session files are
%                   stored.
%
%  OUTPUT
%
%    spikes         list of [t,ID] couples, where the IDs are ordered according 
%                   to brain region
%    regionID       list of brain region identifiers matched for each ID of 
%                   "spikes"
%    regionNames    list of the names of all the brain regions as labelled in 
%                   cell_metrics.brainRegion
%
%  EXAMPLE
% 
%  [spikes,regionID,regionNames] = GetAyaSpikes('N:\V1test\AO52\day3');
% The first spike of the session was fired at spikes(1,1) seconds (e.g. 0.1578 seconds). 
% This spike was fired by unit spikes(1,2) (e.g. unit 4). This unit corresponds to the brain 
% region labelled as regionNames{regionID(spikes(1,2))}. Here, regionID(4) was 3, so the third 
% brain region (alphabetical order), corresponding to regionNames{3}, which was 'SomatosensoryCx'.
%
% Copyright (C) 2022 by Ralitsa Todorova
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

basename = basenameFromBasepath(basepath);
if exist(fullfile(basepath,[basename '.cell_metrics.cellinfo.mat']),'file')
    load(fullfile(basepath,[basename '.cell_metrics.cellinfo.mat']),'cell_metrics');
    regions = cell_metrics.brainRegion;
    regionNames = unique(regions);
    regionCell = zeros(length(regions),1);
    for i=1:length(regionNames)
        regionCell(strcmp(regions,regionNames{i})) = i;
    end
    spikesCell = cell_metrics.spikes.times';
    % reorder based on region:
    [regionID,order] = sort(regionCell);
    spikesCell = spikesCell(order);
    % make a second ID column
    for u=1:length(spikesCell)
        spikesCell{u,1}(:,2) = u; 
    end
    spikes = sortrows(cell2mat(spikesCell));
    spikesCell = cell_metrics.spikes.times(order)';
else
    spikes = zeros(0,2); regionID = zeros(0,1); regionNames = {};
end











