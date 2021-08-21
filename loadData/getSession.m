function session = getSession(varargin)

p = inputParser;
addParameter(p,'basepath',pwd,@isstr)

parse(p,varargin{:});
basepath = p.Results.basepath;

sesfile = checkFile('basepath',basepath,'fileType','.session.mat');

load([sesfile.folder,filesep,sesfile.name],'session');

end