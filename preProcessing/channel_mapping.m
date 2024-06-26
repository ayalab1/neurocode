function channel_mapping(varargin)
%
%
% [channel_mapping: updates basename.session with brain regions per channel
% and creates a adjustable .csv with brain region labels]
%
% [In the csv, each column is a shank and each row is a level on your probe
% corresponding to each channel. This can be adjusted as you drive down,
% simply by copying the .csv from the last session and making minor edits]
%
% Usage:
%       [1. To be ran after basename.session is created to set up initial
%       anatomical map for an animal. This will set up a unlabeled .csv.
%       You are then to label each channel and run this function again to
%       add that data to basename.session. IMPORTANT NOTE: Label the .csv 
%       just as you see the channels in Neuroscope Multi-column view. 
%
%       2. Run gui_session(session) and update brain
%       region manually, then run this function to be consistent.
%
%       3. To be run to integrate older data.
%       Ex: channel_mapping('pull_from_cell_metrics',true)
%       This approach will look within cell_metrics.cellinfo for the
%       location where each unit fired max. Manual adjustment of the
%       anatomical map csv may be need after as many channels will be
%       unknown if the cell counts are low]
%
%  INPUTS
%    [parser]      [input parser, see below]
%    <options>   optional list of property-value pairs (see table below)
%    =========================================================================
%     Properties            Values
%    -------------------------------------------------------------------------
% ['basepath']                   [basepath paramter. Default pwd]             
% ['fig'],true)                  [simple debugging/summary figs. Default true]
% ['save_csv']                   [save output to basepath. Default true]
% ['pull_from_cell_metrics']     [to populate map from cell_metrics. Default
%                                 False]
% ['force_cell_metric_overwrite'][will bypass the warning if .session
%                                 already has regions. Default false]
% ['save_session']               [save session to basepath. Default true]
% ['save_cell_metrics'],true) %  [save updated cell metrics to basepath.
%                                 Default true]
% ['session']                    [session.mat file]
% ['show_gui_session']           [Show gui session. Default false]
%    ========================================================================
%   
%  OUTPUT
%   NA
%
%
% NOTE 
%  [* if you make a manual edit in the .csv, make sure to run this
%       function again to propagate that change to basename.session] 
%
% TODO: Add functionality to optionally pull regions from basename.session
%           to populate the .csv after .csv has been populated. One may 
%           wish to update either basename.session or .csv. 
%
%
% [Ryan H] [2021-2022]
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.
%
%------------------------------------------------------------------------


p = inputParser;
addParameter(p,'basepath',pwd) % path to folder
addParameter(p,'fig',true) % simple debugging/summary figs
addParameter(p,'save_csv',true) % save output to basepath
addParameter(p,'pull_from_cell_metrics',false) % to populate map from cell_metrics
addParameter(p,'force_cell_metric_overwrite',false) % will bypass the warning if .session already has regions
addParameter(p,'save_session',true) % save session to basepath
addParameter(p,'save_cell_metrics',true) % save updated cell metrics to basepath
addParameter(p,'session',[])
addParameter(p,'show_gui_session',false)

parse(p,varargin{:})
basepath = p.Results.basepath;
fig = p.Results.fig;
save_csv = p.Results.save_csv;
pull_from_cell_metrics = p.Results.pull_from_cell_metrics;
force_cell_metric_overwrite = p.Results.force_cell_metric_overwrite;
save_session = p.Results.save_session;
save_cell_metrics = p.Results.save_cell_metrics;
session = p.Results.session;
show_gui_session = p.Results.show_gui_session;

% get the basename and load your basename.session if none was provided
basename = basenameFromBasepath(basepath);
if isempty(session)
    load(fullfile(basepath,[basename,'.session.mat']))
end

% generate anatomical and channel map from .session
% anatomical map will be blank at this point
[anatomical_map,channel_map] = get_maps(session);

% get the number of unknown channels
initial_unknown = get_unknown_count(anatomical_map);

% populate from cell metrics
% here we use previously labeled units to label channels
% this is done to integrate older data that has already been labeled
if pull_from_cell_metrics
    anatomical_map = get_region_from_cell_metrics(basepath,...
        basename,...
        channel_map,...
        anatomical_map,...
        session,...
        force_cell_metric_overwrite);
end

% pull from csv that has been already been generated
pull_from_session = false;
if ~pull_from_cell_metrics
    [anatomical_map,pull_from_session] = get_anatomical_map_csv(basepath,...
        anatomical_map);
end

% check if session has already been updated and generate map from that
if pull_from_session
    anatomical_map = get_map_from_session(session,anatomical_map,channel_map);
end

% check if anatomical_map has been populated
if get_unknown_count(anatomical_map) == initial_unknown
    disp('regions not found... edit your .csv')
    return
end

% update session with anatomical_map data
session = get_updated_session(anatomical_map,channel_map,session);

% save .session
if save_session
    save(fullfile(basepath,[basename,'.session.mat']),'session')
end

% Write the table to a CSV file
if save_csv
    writetable(cell2table(anatomical_map),...
        fullfile(basepath,'anatomical_map.csv'),...
        'WriteVariableNames',0)
end

% update cell labels in cell metrics
if save_cell_metrics &&...
        exist(fullfile(basepath,[basename,'.cell_metrics.cellinfo.mat']),'file')
    disp('updating units in cell metrics')
    load(fullfile(basepath,[basename,'.cell_metrics.cellinfo.mat']))
    chListBrainRegions = findBrainRegion(session);
    for j = 1:cell_metrics.general.cellCount
        if isfield(session,'brainRegions') && ~isempty(session.brainRegions)
            try
                cell_metrics.brainRegion{j} = chListBrainRegions{cell_metrics.maxWaveformCh1(j)};
            catch
                cell_metrics.brainRegion{j} = 'Unknown';
            end
        end
    end
    save(fullfile(basepath,[basename,'.cell_metrics.cellinfo.mat']),'cell_metrics')
end
if save_cell_metrics &&...
        exist(fullfile(basepath,[basename,'.unsorted.cell_metrics.cellinfo.mat']),'file')
    disp('updating units in the unsorted cell metrics')
    load(fullfile(basepath,[basename,'.unsorted.cell_metrics.cellinfo.mat']))
    chListBrainRegions = findBrainRegion(session);
    for j = 1:cell_metrics.general.cellCount
        if isfield(session,'brainRegions') && ~isempty(session.brainRegions)
            try
                cell_metrics.brainRegion{j} = chListBrainRegions{cell_metrics.maxWaveformCh1(j)};
            catch
                cell_metrics.brainRegion{j} = 'Unknown';
            end
        end
    end
    save(fullfile(basepath,[basename,'.unsorted.cell_metrics.cellinfo.mat']),'cell_metrics')
end

if fig
    try 
        chanCoords = getStruct(basepath,'chanCoords'); 
        generateChannelMap1(session,anatomical_map,channel_map,chanCoords); 
    catch
        generateChannelMap1(session,anatomical_map,channel_map); 
    end
    exportgraphics(gcf,fullfile(basepath,'anatomical_map.png'),'Resolution',150)
end

if show_gui_session
    gui_session(session)
end
end

function initial_unknown = get_unknown_count(anatomical_map)
anatomical_map_vec = anatomical_map(:);
anatomical_map_vec = anatomical_map_vec(~cellfun('isempty',anatomical_map_vec));
initial_unknown = sum(contains(anatomical_map_vec,'Unknown'));
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

function anatomical_map = get_region_from_cell_metrics(basepath,...
                                                        basename,...
                                                        channel_map,...
                                                        anatomical_map,...
                                                        session,...
                                                        force_cell_metric_overwrite)
if ~force_cell_metric_overwrite                                                    
    if isfield(session,'brainRegions')
        warning('Careful, basename.session already has brainRegions field, you may by overwriting')
        disp(session.brainRegions)
        continue_ = input('type 1 to continue:   ');
    else
        continue_ = 1;
    end
    if continue_ ~= 1
        return
    end
end
load(fullfile(basepath,[basename,'.cell_metrics.cellinfo.mat']))
for i = 1:length(cell_metrics.maxWaveformCh1)
    anatomical_map(channel_map == cell_metrics.maxWaveformCh1(i)) =...
        cell_metrics.brainRegion(i);
end
end

function [anatomical_map,pull_from_session] = get_anatomical_map_csv(basepath,anatomical_map)
pull_from_session = false;
filename = fullfile(basepath,'anatomical_map.csv');
if ~exist(filename,'file')
    warning('no .anatomical_map.csv... ')
    disp('will try to pull from basename.session')
    disp('you can check anatomical_map and rerun')
    
    writetable(cell2table(anatomical_map),...
        fullfile(basepath,'anatomical_map.csv'),...
        'WriteVariableNames',0)
    pull_from_session = true;
    return
end
anatomical_map = table2cell(readtable(filename,'ReadVariableNames',false));
end

function session = get_updated_session(anatomical_map,channel_map,session)
% mat to vectors
anatomical_map_vec = anatomical_map(:);
channel_map_vec = channel_map(:);

% remove non existant channels
nonempty_idx = channel_map(:) >= 0;
anatomical_map_vec = anatomical_map_vec(nonempty_idx);
channel_map_vec = channel_map_vec(nonempty_idx);

% iter through each region and add to session
regions = unique(anatomical_map_vec);
regions = regions(~cellfun('isempty',regions));

% clear out brain regions in order to populate 
session.brainRegions = [];

for region = regions'
    session.brainRegions.(region{1}).channels =...
        channel_map_vec(ismember(anatomical_map_vec,region))';
    % add electrode groups
    [~,c]=find(ismember(channel_map,session.brainRegions.(region{1}).channels));
    session.brainRegions.(region{1}).electrodeGroups = unique(c)';
end
end

function generateChannelMap1(session,anatomical_map,channel_map,chanCoords)

channel_map_vec = channel_map(:);
anatomical_map_vec = anatomical_map(:);
anatomical_map_vec = anatomical_map_vec(~isnan(channel_map_vec));
channel_map_vec = channel_map_vec(~isnan(channel_map_vec));
anatomical_map_vec(contains(anatomical_map_vec,'Unknown')) = {''};

for i = 1:length(anatomical_map_vec)
    label{i} = [anatomical_map_vec{i},' ',num2str(channel_map_vec(i))];
end

% check if chanCoords exist and generate using defaults
% also, save updated session file with generated coords
if ~exist('chanCoords','var') && ~isfield(session.extracellular,'chanCoords')
    session.extracellular.chanCoords = generateChanCoords(session);
    basepath = session.general.basePath;
    if exist(basepath,'dir')
        save(fullfile(basepath,[basenameFromBasepath(basepath),'.session.mat']),'session')
    end
end

% attempt to pull channel coords from session
if isfield(session.extracellular,'chanCoords')
    chanCoords.x = session.extracellular.chanCoords.x;
    chanCoords.y = session.extracellular.chanCoords.y;
end

% if channel coords are empty, try to create them here
if ~exist('chanCoords','var') || isempty(chanCoords.x)
    try
       chanMap = generateChannelMap(session,'reorder',false);
    catch
       chanMap = createChannelMap(session,'reorder',true);
    end
    chanCoords.x = chanMap.xcoords(:);
    chanCoords.y = chanMap.ycoords(:);
end
% chanCoords.source = chanMap.source;
% chanCoords.layout = chanMap.layout;
% chanCoords.shankSpacing = chanMap.shankSpacing;

x_range = range(chanCoords.x);
y_range = range(chanCoords.y);
if x_range > y_range
    fig_width = 1600;
    fig_height = ceil(fig_width*y_range/x_range)+200;
else
    fig_height = 1000;
    fig_width = ceil(fig_height*x_range/y_range)+200;
end
fig1 = figure('Name','Channel map','position',[5,5,fig_width,fig_height]);
movegui(fig1,'center')
plot(chanCoords.x,chanCoords.y,'.k'), hold on

xlim([min(chanCoords.x) - 50 ,max(chanCoords.x) + 50])
ylim([min(chanCoords.y) - 50 ,max(chanCoords.y) + 50])

[~,idx] = sort(channel_map_vec);
text(chanCoords.x,chanCoords.y,anatomical_map_vec(idx),...
    'VerticalAlignment','bottom','HorizontalAlignment','center','fontsize',8);
title({' ','Channel map',' '}), xlabel('X (um)'), ylabel('Y (um)')

end
