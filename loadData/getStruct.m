function [struct,filename] = getStruct(basepath,extension,fieldName)

%getStruct - Load custom data extensions associated with a given session.

% This is a helper function similar to CellExplorer's loadStruct.
% If multiple .mat files match your query, the function will return
% the first one (in alphabetical order).
%
%    =========================================================================
%  USAGE
%
%%INPUT
%    basepath          full path where the session is located (default = pwd)
%    extension         extension to look for in the filename of the .mat file
%                      (default = 'ripples')
%
%    =========================================================================

%OUTPUT
%   struct         

%  EXAMPLES
%
%    ripples = getStruct(basepath,'ripples');  % load the ripples structure
%    MergePoints = getStruct(basepath,'Merge');  % load the MergePoints structure
%
%    =========================================================================
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
indices = find(ok); 
if isempty(indices)
    disp(['No ''' extension ''' file found in ' basepath]);
    error(['No ''' extension ''' file found in ' basepath]);
    return
end
if nargin<3,fieldName = extension; end
if length(indices)>1,
    candidates = filenames(indices);
    containsBasename = cellfun(@(x) contains(x,basenameFromBasepath(basepath)),candidates);
    if any(containsBasename), candidates(~containsBasename) = []; end
    % pick the shortest one of the remaining filenames (to avoid taking e.g. ripples.events.old.mat rather than ripples.events.mat)
    [~,index] = sort(cellfun(@length,candidates));
    filename = candidates{index(1)};
    if length(indices)<5 % display the filenames
        str = filenames{indices(1)}; for i=2:length(indices), str = [str '; ' filenames{indices(i)}]; end
        disp(['Multiple (' num2str(length(indices)) ') ' extension ' files exist in ' basepath ':']);
        disp(str);
        disp(['Choosing file "' filename '".']);
%         cd(basepath); keyboard
    else
        disp(['Multiple (' num2str(length(indices)) ') ' extension ' files exist in ' basepath '. Choosing file "' filename '".']);
    end
    
else
    filename = filenames{indices};
end

[parts] = strsplit(filename,'.');

state = warning; state = state(1).state; warning('off'); % supress warnings for the following code
if ~exist(fullfile(basepath,filename),'file')
    error([fullfile(basepath,filename) ' does not exist.']);
end

vars = whos('-file',fullfile(basepath,filename)); % list the variables actually contained in the file
if length(vars)==1 % if there is a single variable, load that variable regardless of what it's called inside the file
    struct = load(fullfile(basepath,filename),vars.name);
    struct = struct.(vars.name);
else % try to find the variable specified in "fieldName"
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
end
warning(state); % return to previous warning state
