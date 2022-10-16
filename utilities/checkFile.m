function file = checkFile(varargin)

% checkfile - Checks that file exists and returns proper basepath and filename
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
addParameter(p,'fileTypes',[],@iscell)
addParameter(p,'searchSubdirs',true,@islogical)
addParameter(p,'searchSuperdirs',true,@islogical)

parse(p,varargin{:})
filename = p.Results.filename;
basepath = p.Results.basepath;
fileType = p.Results.fileType;
fileTypes = p.Results.fileTypes;
searchSubdirs = p.Results.searchSubdirs;
searchSuperdirs = p.Results.searchSuperdirs;

if ~isempty(fileType) && ~isempty(fileTypes)
    error('Why dont you put everything inside fileTypes...')
end

if isempty(filename) && isempty(fileType) && isempty(fileTypes)
    error('Need a filename or a file type...gotta gimme something')
end

if ~isempty(fileType) && ~strcmp(fileType(1),'.')
    fileType = ['.',fileType];
end

if ~isempty(fileTypes)
    for i=1:length(fileTypes)
         if ~strcmp(fileTypes{i}(1),'.')
            fileTypes{i} = ['.',fileTypes{i}];
         end
    end
end

if isempty(filename) && ~isempty(fileType)
    file = dir([basepath,filesep,'*',fileType]);
    if isempty(file) && searchSubdirs
        subFile = dir([basepath,filesep,'*',filesep,'*',fileType]);
        ignore = cellfun(@(x) strcmp(x,'settings.xml'),{subFile(1:end).name})'; % ignore "settings.xml" files
        subFile(ignore) = [];
        file = [file; subFile];
        if ~isempty(file)
            warning(['Found file of type: ',fileType,' in subdirectory']);
        end
    end
    if isempty(file) && searchSuperdirs
        newdir = fileparts(basepath);
        superFile = dir([newdir,filesep,'*',fileType]);
        file = [file; superFile];
        if ~isempty(file)
            warning(['Found file of type: ',fileType,' in parent directory']);
        end
    end
    if isempty(file)
        error(['Cant find file of type: ',fileType,' in this location: ',basepath])
    end
elseif isempty(filename) && ~isempty(fileTypes)
    for i =1:length(fileTypes)
        file_tmp{i} = dir([basepath,filesep,'*',fileTypes{i}]);
            if searchSubdirs
                subFile = dir([basepath,filesep,'*',filesep,'*',fileTypes{i}]);
                file_tmp{i} = [file_tmp{i}; subFile];
            end
        if ~isempty(file_tmp{i})
            file = [file_tmp{i}; subFile];
            return
        end
    end
else
    if ~strcmp(filename(end-(length(fileType)-1):end),fileType)
        filename = [filename, fileType];
    end
    file = dir([basepath,filesep,filename]);
    if searchSubdirs
        subFile = dir([basepath,filesep,'*',filesep,filename]);
        file = [file; subFile];
        if ~isempty(file)
            warning(['Found file of type: ',fileType,' in subdirectory']);
        end
    end
    if isempty(file) && searchSuperdirs
        mydir  = pwd;
        idcs   = strfind(mydir,'\');
        newdir = mydir(1:idcs(end)-1);
        superFile = dir([newdir,filesep,'*',fileType]);
        file = [file; superFile];
        if ~isempty(file)
            warning(['Found file of type: ',fileType,' in parent directory']);
        end
    end
    if isempty(file)
        error(['Cant find file: ',filename,' in this location or near: ',basepath])
    end
end
end
