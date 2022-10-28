function uiimport(varargin)

% This function overrides the original "uiimport" function
% It asks the user to confirm if they want to open a large (>1GB file), 
% and cancels the action if "No" is clicked.
% This prevents accidental openings of .dat file, which can freeze matlab.
%
% Copyright (C) 2022 by Ralitsa Todorova
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.


currentFolder = pwd;
paths = which('-all','uiimport'); 
builtInPath = fileparts(paths{2});
cd(builtInPath);
handle = str2func('uiimport');
cd(currentFolder);

q = dir(varargin{1});
q(1).bytes
if q(1).bytes/1024/1024/1024>1 % more than 1GB file
    button = questdlg(...
        ['Are you sure you want to open this ' num2str(q(1).bytes/1024/1024/1024) ' GB file?'], ...
        'MATLAB','Yes','No','No');
    switch button
        case 'No',
            return
    end
end

handle(varargin{:});