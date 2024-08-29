function basepathList = updateBasepathList(parentDirectory,varargin)

%updateBasepathList - get a list of valid basepaths within the queried 
% directory. Folders which contain a 'session.mat' file will be returned
% in a cell and saved in a 'basepathList.mat' file for later loading.
%
% Copyright (C) 2024 Ralitsa Todorova
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

p = inputParser;
addParameter(p,'saveDir',parentDirectory,@(x) isempty(x) || isfolder(x)); % by default, queried directory

% addParameter(p,'pullData',[],@isdir); To do...
parse(p,varargin{:});
saveDir = p.Results.saveDir;

basepathList = {};
if iscell(parentDirectory) % Call the function separately for each cell
    for i=1:length(parentDirectory)
        output = updateBasepathList(parentDirectory{i},varargin{:});
        if ~isempty(output)
            basepathList = cat(1,basepathList,output);
        end
    end
    return
end

list = dir(parentDirectory);


% Is this a basepath?
% If a "session" file exists, this is a basepath. Don't go into subfolders
basename = basenameFromBasepath(parentDirectory);
if exist(fullfile(parentDirectory,[basename '.session.mat']),'file')
    basepathList = {parentDirectory};
    return
end

% If no "session" file exists, this is not a basepath, and to look for basepaths,
% we need to look into deeper subfolders.

list = list(~cellfun(@(x) ismember(x(1),{'.','#'}),{list.name})); % remove subfolders that start with a dot ('.'). This includes the '.phy' folder
subdirectories = list([list.isdir]);
subdirectories = cellfun(@(x,y) fullfile(x,y),{subdirectories.folder}',{subdirectories.name}','UniformOutput',false);
for i=1:length(subdirectories)
    output = updateBasepathList(subdirectories{i},'saveDir',[]); % Don't save output during these iterations
    if ~isempty(output)
        basepathList = cat(1,basepathList,output);
    end
end

if ~isempty(saveDir)
    save(fullfile(saveDir,'basepathList'),'basepathList');
end

