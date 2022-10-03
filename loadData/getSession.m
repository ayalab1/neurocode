function session = getSession(varargin)

%[getSession] -  finds session.mat struct for given basepath

%    =========================================================================
%  USAGE
%
%INPUT
%   [basePath]         [directory: '/whatevetPath/baseName/']

%    =========================================================================

%OUTPUT
%   session         .session.mat file for the basepath 
%   =========================================================================


% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

p = inputParser;
addParameter(p,'basepath',pwd,@isfolder)

parse(p,varargin{:});
basepath = p.Results.basepath;

sesfile = checkFile('basepath',basepath,'fileType','.session.mat');

% weird error running some files, temporary fix
if size(sesfile, 1) > 1
    sesfile = sesfile(1,:);
    warning('sesfile is too large, confirm that file chosen is correct');
    disp(['folder: ' sesfile.folder]);
    disp(['file: ' sesfile.name]);
end

load([sesfile.folder,filesep,sesfile.name],'session');

end
