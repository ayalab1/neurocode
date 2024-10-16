function basepathList = getBasepathList(parentDirectory)

%Load a list of basepaths within the queried directory which
% have been previously saved with "updateBasepathList".
%
% EXAMPLES:
% % get a list of basepaths within a project folder:
% basepathList = getBasepathList('Y:\OJRproject'); 
% % get a list of basepaths within the whole ayadata3:
% basepathList = getBasepathList('X:\'); 
% % get a list of basepaths within any of the current servers 
% % (aleph, aleph2, ayadata1-4, ayadataB1-B4)
% basepathList = {'K:\','L:\','T:\','U:\','V:\','W:\','X:\','Y:\','Z:\'});
%
% Copyright (C) 2024 Ralitsa Todorova
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

basepathList = {};
if iscell(parentDirectory) % Call the function separately for each cell
    for i=1:length(parentDirectory)
        output = getBasepathList(parentDirectory{i});
        if ~isempty(output)
            basepathList = cat(1,basepathList,output);
        end
    end
else
    if exist(fullfile(parentDirectory,'basepathList.mat'),'file')
        load(fullfile(parentDirectory,'basepathList.mat'),'basepathList');
    else
        error(['File ' fullfile(parentDirectory,'basepathList.mat') ' not found. Please call updateBasepathList to create it.']);
    end
end
