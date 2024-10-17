function cell_metrics = multiCellMetrics(basepath, ifSave, ifOpen)

% [multiCellMetrics- to visualize multiple sessions in cellExplorer GUI]
%
%  [Helper function to load and visualize multiple sessions within the
%  cellExplorer GUI. After you have multiple sessions spike sorted and 
%  processed, you can open them all in cell explorer]

%
%  USAGE
%
%    [Load all cell_metrics from a project
%    multiCellMetrics('X:\data\Barrage');]
%
%    [Load all cell_metrics from an animal
%    multiCellMetrics('X:\data\Barrage\NN2');]
%
%    [Load a custom list of sessions (provide basepaths and basenames)
%    multiCellMetrics(customBasepaths, customBasenames);]
%
%  INPUTS
%  [basepath]       [If asking the function to search for all cell_metrics
%                   folders within a given path, provide a character string
%                   leading to the folder of interest (ie project or animal
%                   folder). Otherwise, if interested in loading custom
%                   sessions, provide an Nx1 cell of basepaths of interest
%                   (characters within each cell)]
%   [ifSave]       [Option to save the created dataframe to the relevant
%                   folder. Default is false]  
% 
%  OUTPUTS
%    N/A
%
%
%  SEE ALSO
%
%    [PreprocessSpikes]
%
%
% [Ryan Harvey] [2021-2022]
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

if nargin < 1
    error('Please input session basepaths or a parent folder path');
elseif nargin < 2
    ifSave = 0;
    ifOpen = 1;
elseif nargin < 3
    ifOpen = 1;
end

if isstring(basepath)||ischar(basepath)
    if contains(basepath,'.csv')
        useSess = readtable(basepath);
        for i = 1:size(useSess,1)
           basepath_use{i} = useSess.basepath{i};
           basename_use{i} = useSess.basename{i};
        end
    else
        cd(basepath);
        data_path = basepath;

        % look for all the cell_metrics.cellinfo.mat files
        files = dir([data_path,'\**\*.cell_metrics.cellinfo.mat']);

        % pull out basepaths and basenames
        for i = 1:length(files)
            basepath_use{i} = files(i).folder;
            basename_use{i} = basenameFromBasepath(files(i).folder);
        end
    end
elseif iscell(basepath)
    basepath_use = basepath;
    for i = 1:length(basepath_use)
        basename_use{i} = basenameFromBasepath(basepath{i});
    end
else
    error('Incorrect basepath input. Please input either a character basepath or a cell array of custom basepaths');
end

if ifSave
    % write data frame of basepaths to project folder
    % using this .csv, you will not need to again search for all the cell_metrics
    % files, which can take some time if you have many sessions
    df = table();
    df.basepath = basepath_use';
    df.basename = basename_use';
    writetable(df,['sessions_',date,'.csv']);
end

% load all cell metrics
cell_metrics = loadCellMetricsBatch('basepaths',basepath_use,'basenames',basename_use);

% pull up gui to inspect all units in your project
if ifOpen
    cell_metrics = CellExplorer('metrics',cell_metrics);
end
end