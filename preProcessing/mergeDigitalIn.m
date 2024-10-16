function digitalIn = mergeDigitalIn(varargin)

% This function merges the digitalIn.events.mat files inside subsession 
% folders and saves the merged digitalIn.events.mat into the basepath folder.
% Subsession timestamps (taken from MergePoints) are added to each subsession's
% digitialIn structure. 
% 
%  USAGE
%
%    digitalIn = mergeDigitalIn(<options>)
%
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'basepath'    path to the folder containing all session folders (default = 
%                   current directory)
%     'forceDetect' a boolean indicating whether to overwrite an existing
%                   digitalIn.events.mat file in the basepath (defeault = false)
%    =========================================================================
%
%  OUTPUT
%     Operates on files in specified folder.  No output variable
%
% Copyright (C) 2024 by Ralitsa Todorova
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

p = inputParser;
addParameter(p,'basepath',pwd,@isfolder); % by default, current folder
addParameter(p,'forceDetect',false,@islogical); 

parse(p,varargin{:});
basepath = p.Results.basepath;
forceDetect = p.Results.forceDetect;

filename = fullfile(basepath,'digitalIn.events.mat');
if ~forceDetect && exist(filename,'file')
    disp([filename ' already exists. Aborting...']);
    return
end

% Get folder list from MergePoints:
MergePoints = getStruct(basepath,'MergePoints');
folderList = MergePoints.foldernames';
add = MergePoints.timestamps(:,1);
toMerge = false(size(folderList,1),1); digsToMerge = cell(size(folderList,1),1); nEventTypes = zeros(size(folderList,1),1);
for subfolderID = 1:size(folderList,1)
    if exist(fullfile(basepath,folderList{subfolderID},'digitalIn.events.mat'),'file')
        load(fullfile(basepath,folderList{subfolderID},'digitalIn.events.mat'),'digitalIn');
        % add the subsession start to all event timestamps
        fields = fieldnames(digitalIn);
        nEventTypes(subfolderID,1) = length(digitalIn.(fields{1}));
        for j=1:length(fields)
            field = fields{j};
            if ~strcmp(field,'dur') % don't add timestamps to "durations" 
                digitalIn.(field) = cellfun(@(x) x+add(subfolderID),digitalIn.(field),'UniformOutput',false);
            end
        end
        digsToMerge{subfolderID} = digitalIn;
        toMerge(subfolderID) = true;
    end
end

% start with a structure with the maximum number of event types
digitalIn = digsToMerge{find(nEventTypes==max(nEventTypes),1)};
fields = fieldnames(digitalIn);
for j=1:length(fields)
    field = fields{j};
    for k=1:max(nEventTypes)
        cellsToMerge = cellfun(@(x) x.(field){k},digsToMerge(nEventTypes>=k),'UniformOutput',false);
        if any(strcmp(field,{'ints','dur'}))
            digitalIn.(field){k} = cat(2,cellsToMerge{:});
        else
            digitalIn.(field){k} = cat(1,cellsToMerge{:});
        end
    end
end









%% Save merged events

save(filename,'digitalIn');














