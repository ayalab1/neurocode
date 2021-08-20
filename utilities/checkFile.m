function file = checkFile(varargin)

% Checks that file exists and returns proper basepath and filename
% intpus (optional):
%   filename (optional) - name of file, MUST have this or filetype
%   basepath (optional) - path to file, default current working directory
%   fileType (optional) - .whatever filetype (only needed if not part of
%                           filename)
%
% outputs:
%   file - matlab file struct
%
% kathryn mcclain 2020

p = inputParser;
addParameter(p,'filename',[],@ischar)
addParameter(p,'basepath',pwd,@ischar)
addParameter(p,'fileType',[],@ischar)
addParameter(p,'searchSubdirs',true,@islogical)

parse(p,varargin{:})
filename = p.Results.filename;
basepath = p.Results.basepath;
fileType = p.Results.fileType;
searchSubdirs = p.Results.searchSubdirs;

if isempty(filename) && isempty(fileType)
    error('Need a filename or a file type...gotta gimme something')
end

if ~isempty(fileType) && ~strcmp(fileType(1),'.')
    fileType = ['.',fileType];
end

if isempty(filename)
    file = dir([basepath,filesep,'*',fileType]);
    if searchSubdirs
        subFile = dir([basepath,filesep,'*',filesep,'*',fileType]);
        file = [file; subFile];
    end
    if isempty(file)
        error(['Cant find file of type: ',fileType,' in this location: ',basepath])
    end
else
    if ~strcmp(filename(end-(length(fileType)-1):end),fileType)
        filename = [filename, fileType];
    end
    file = dir([basepath,filesep,filename]);
    if searchSubdirs
        subFile = dir([basepath,filesep,'*',filesep,filename]);
        file = [file; subFile];
    end
    if isempty(file)
        error(['Cant find file: ',filename,' in this location: ',basepath])
    end
end
end


