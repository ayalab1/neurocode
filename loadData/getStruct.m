function struct = getStruct(basepath,extension,fieldName)

% Load custom data extensions associated with a given session.
% This is a helper function similar to CellExplorer's loadStruct.
% If multiple .mat files match your query, the function will return
% the first one (in alphabetical order).
%
%  USAGE
%
%    struct = getStruct(basepath,extension)
%
%    basepath          full path where the session is located (default = pwd)
%    extension         extension to look for in the filename of the .mat file
%                      (default = 'ripples')
%
%  EXAMPLES
%
%    ripples = getStruct(basepath,'ripples');  % load the ripples structure
%    MergePoints = getStruct(basepath,'Merge');  % load the MergePoints structure
%
% Copyright (C) 2022 by Ralitsa Todorova
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

list = dir(basepath);
filenames = cat(1,{list.name})'; filenames(cellfun(@length,filenames)<3) = [];
ismat = cellfun(@(x) ~isempty(strfind(x,'.mat')),filenames);
ok = cellfun(@(x) ~isempty(strfind(x,extension)),filenames) & ismat;
index = find(ok,1);
if isempty(index)
    error(['No ''' extension ''' file found in ' basepath]);
    return
end
if nargin<3,fieldName = extension; end
filename = filenames{index};
[parts] = strsplit(filename,'.');

state = warning; state = state(1).state; warning('off'); % supress warnings for the following code
try
%     fieldName = parts{end-2}; % parts{end} is mat
    struct = load(fullfile(basepath,filename),fieldName);
    struct = struct.(fieldName);
catch
    try fieldName = parts{end-1};
        struct = load(fullfile(basepath,filename),fieldName);
        struct = struct.(fieldName);
    catch fieldName = parts{end-2};
        struct = load(fullfile(basepath,filename),fieldName);
        struct = struct.(fieldName);
    end
end
warning(state); % return to previous warning state
