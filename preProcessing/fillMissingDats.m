function fillMissingDats(varargin)

p = inputParser;
otherdattypes = {'analogin';'digitalin';'auxiliary';'time';'supply'};
isFileType = @(x) sum(strcmp(x,otherdattypes))==1;
addParameter(p,'basepath',cd,@isstr)
addParameter(p,'fileType',[],isFileType)

parse(p,varargin{:})
basepath = p.Results.basepath;
fileType = p.Results.fileType;

%% Get session info
session = getSession('basepath',basepath); % Peter's sessionInfo
basename = session.general.name;
ampNch = session.extracellular.nChannels;

%% check location of each file and what to concat
typeFiles = dir([basepath,filesep,'*',filesep,fileType,'.dat']);
ampFiles = dir([basepath,filesep,'*',filesep,'amplifier.dat']);
typeFolders = {typeFiles.folder};
ampFolders = {ampFiles.folder};

fillInds = find(~ismember(ampFolders,typeFolders)); %index of folders that need fill

% calculate number of channels to for fill file
refTypeInd = find(ismember(typeFolders,ampFolders),1);
refTypeSize = typeFiles(refTypeInd).bytes;
refAmpInd = find(strcmp(ampFolders,typeFiles(refTypeInd).folder));
refAmpSize = ampFiles(refAmpInd).bytes;
typeNch = refTypeSize*ampNch/refAmpSize;

%% Fill in missing files
for ii = 1:length(fillInds)
    fillIdx = fillInds(ii);
    localAmpSize = ampFiles(fillIdx).bytes;
    nPoints = localAmpSize/(ampNch*2);
    
    zeroData = zeros(typeNch,nPoints);
    
    filepath = [ampFiles(fillIdx).folder,filesep,fileType,'.dat'];
    fid = fopen(filepath,'w');
    fwrite(fid,zeroData,'int16');
    fclose(fid);
end
