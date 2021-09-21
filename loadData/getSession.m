function session = getSession(varargin)

p = inputParser;
addParameter(p,'basepath',pwd,@isstr)

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