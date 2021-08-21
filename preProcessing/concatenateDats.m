function ConcatenateDats(varargin)
% bz_ConcatenateDats - Concatenate raw .dat files found in a session folder
% - for intan type recordings 
% 
% ALGORITHM OUTLINE: looks for .dat files in a folder (or in subfolders) to
% concatenate together.  The concatenation happens via system commands 
% ("cat" command for linux/mac, "copy" command if windows/pc).  Uses
% different assumptions to find and recognize relevant .dats depending on
% the acquisition system.  
% 
% REQUIREMENTS: Assumes you are in or pointed to a directory containing 
% subdirectories for various recording files from a single session. *It is 
% assumed that an earlier-acquired data file/folder will have a name that
% is sorted alphanumerically earlier.  Alphanumeric sorting order is
% assumed to be the recording temporal sequence.
% Works with acquisition systems: Intan  - 
%   1) intan: wherein subfolders are inside the session folder.  Each
%   subfolder contains simultaneously-recorded .dat files recorded for a
%   continuous period of time.  Start/stop recording commands each create a
%   new folder.  *It is assumed that the alphanumeric sorting of these 
%   folders corresponds with their sequence in acquisiton time.*  
%   These folders contain
%       - info.rhd files with metadata about the recording. 
%       - amplifier.dat - int16 file with usually neural data from the
%           headstage
%       - auxiliary.dat (optional) - uint16 file from auxiliary channels on
%           the headstage - often accelerometer
%       - analogin.dat (optional) - uint16 file from analogin channels on 
%           main board 
%       - digitalin.dat (optional) - uint16 file recording all 16 digital 
%           channels on main board 
%       - time.dat - int32 file giving recording sample indicies (e.g. 
%           0,1,2,3...) for each sample recorded in other channels
%       - supply.dat - uint16 file showing voltage supplied to (?preamp?)
%   
%
%  USAGE
%
%    ConcatenateDats(basepath,deleteoriginaldatsbool,sortFiles)
%
%  INPUTS
%
%    basepath          computer path to session folder.  Defaults to
%                      current folder if no input given
%    deleteoriginaldatsbool  - boolean denoting whether to delete (1) or
%                              not delete (0) original .dats after
%                              concatenation.  Default = 0. Not recommended.
%    sortFiles               - boolean denoting whether to sort files according 
%                              to time of recording (1) or
%                              not (0) and thus sort them alphabetically 
%                              Default = 0.
%
%  OUTPUT
%     Operates on files in specified folder.  No output variable
%
%  EXAMPLES
%      Can be called directly or via bz_PreprocessExtracellEphysSession.m
%
% Brendon Watson 2017
% Antonio FR, 2018
% kathryn mcclain 2020


%% Handling inputs
% basic session name and and path
p = inputParser;
addParameter(p,'basepath',cd,@isstr)

parse(p,varargin{:})
basepath = p.Results.basepath;

%% Get session info
session = getSession('basepath',basepath); 
basename = session.general.name;

%% If the dats are already merged quit
if exist(fullfile(basepath,[basename,'.dat']),'file')
    disp('.dat already exists in session directory, not merging subdats')
    return
end

%% check location of each file and what to concat
files = dir([basepath,filesep,'*\*.dat']);
datFolders = {files.folder};
datTypes = {files.name};
subfolders = unique(datFolders);
%sort subfolders
subTimes = cell2mat(cellfun(@(X) str2num(X(end-5:end)),subfolders,'UniformOutput',false));
[~,subOrder] = sort(subTimes);
subfolders = subfolders(subOrder);

types = unique(datTypes);
nSubfolders = length(subfolders);
nTypes = length(types);

fileCheck = zeros(nTypes,nSubfolders);
for subIdx = 1:nSubfolders
    for typeIdx = 1:nTypes
        match = strcmp(datFolders,subfolders(subIdx))&strcmp(datTypes,types(typeIdx));
        fileCheck(typeIdx,subIdx) = find(match);
    end
end

typesToConcat = find(sum(fileCheck~=0,2)==nSubfolders);
if isempty(typesToConcat)
    disp('Nothing to concatenate')
end

%% Concatenate
clear sizecheck
for i = 1:length(typesToConcat)
    type = types{typesToConcat(i)};
    if strcmp(type,'amplifier.dat')
        newPath = [basepath,filesep,basename,'.dat'];
    else
        newPath = [basepath,filesep,type];
    end
    
    oldPaths = cell(nSubfolders,1);
    for j = 1:nSubfolders
        oldPaths{j} = [subfolders{j},filesep,type];
    end
    
    if isunix
        cs = strjoin(oldPaths);
        catstring = ['! cat ', cs, ' > ',newPath];
    elseif ispc%As of 4/9/2017 - never tested
        cs = strjoin(oldPaths, ' + ');
        catstring = ['! copy /b ', cs, ' ',newPath];
    end

    eval(catstring)%execute concatenation

    % Check that size of resultant .dat is equal to the sum of the components
    newFile = dir(newPath);
    newSize = newFile.bytes;
    oldSize = sum([files(fileCheck(typesToConcat(i),:)).bytes]);
    if newSize==oldSize
        disp([type ' concatenated and size checked'])
        sizecheck.(type(1:end-4))= true;
    else
        error('New .dat size not right')
    end
    
end

%% Calculate merge points
timeInd = find(strcmp(types,'time.dat'));

if isempty(timeInd) || sum(fileCheck(timeInd,:)~=0)~=nSubfolders
    disp('Not enough time files so merge points no go')
else
    
    firstLastTimePoints = zeros(nSubfolders,2);
    transitiontimes_samp = zeros(nSubfolders,1);
    for i = 1:nSubfolders
        relFile = files(fileCheck(timeInd,i));
        f = fopen([relFile.folder,filesep,relFile.name],'r');
        % Determine total number of samples in file
        fileStart = ftell(f);
        
        %Read the first time point
        firsttimepoint = fread(f,1,'int32');
        status = fseek(f,-4,'eof'); %int32 = 4 bytes
        lasttimepoint = fread(f,1,'int32');
        fileStop = ftell(f);
        
        firstLastTimePoints(i,:) = [firsttimepoint lasttimepoint];
        numsamples(i) = fileStop./4;
        if i==1
            transitiontimes_samp = firstLastTimePoints(i,:);
        else
            transitiontimes_samp(i,:) = firstLastTimePoints(i,:)+transitiontimes_samp(i-1,2)+1;
        end
    end
    
    transitiontimes_sec = transitiontimes_samp./session.extracellular.sr; %convert to seconds
    
    % Make the events.mat file that saves all the merge information
    
    eventsfilename = fullfile(basepath,[basename,'.MergePoints.events.mat']);
    
    MergePoints.timestamps = transitiontimes_sec;
    MergePoints.timestamps_samples = transitiontimes_samp;
    MergePoints.firstlasttimpoints_samples = firstLastTimePoints;
    MergePoints.foldernames = subfolders;
    MergePoints.filesmerged = files;
    MergePoints.sizecheck = sizecheck;
    MergePoints.detectorinfo.detectorname = 'concatenateDats';
    MergePoints.detectorinfo.detectiondate = datestr(now,'yyyy-mm-dd');
    
    %Saving SleepStates
    save(eventsfilename,'MergePoints');
end



