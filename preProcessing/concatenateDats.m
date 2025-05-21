function concatenateDats(varargin)
% Concatenate raw .dat files found in a session folder
% for Intan AND/OR openEphys recordings
%
%% INPUTS
% basepath                Basepath for experiment. It contains all session
%                         folders. If not provided takes pwd.
% fillMissingDatFiles     Logical option to fill any dat file types which
%                         may not be present in every subsession. Default
%                         is false.
% sortFiles               Logical option to sort files by their Intan or
%                         OpenEphys timestamp. Setting to false will
%                         default to altSort ordering (below). If altSort
%                         is empty, files will be sorted alphabetically.
% altSort                 Numerical array of indices ordering your 
%                         subsession dat files. If, for example, you have 
%                         subsession folders containing dat files labeled 
%                         as FolderA; FolderB; FolderC; and want the order 
%                         to be concatenated as "C, A, B", input altSort as 
%                         [2, 3, 1]; Default is false, which sorts by
%                         date/time (YYMMDD_HHMMSS for Intan, 
%                         YYYY-MM-DD_HH-MM-SS for OpenEphys).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;
addParameter(p, 'basepath', pwd, @isfolder); % by default, current folder
addParameter(p, 'fillMissingDatFiles', false, @islogical);
addParameter(p, 'sortFiles', true, @islogical);
addParameter(p, 'altSort', [], @isnumeric);

parse(p, varargin{:});

basepath = p.Results.basepath;
fillMissingDatFiles = p.Results.fillMissingDatFiles;
sortFiles = p.Results.sortFiles;
altSort = p.Results.altSort;

if sortFiles&&(~isempty(altSort))
    error('sortFiles cannot be empty while altSort provides an order. Please choose to either sort by time (sortFiles=true) or designate a manual order of concatenation (altSort)');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
basename =  basenameFromBasepath(basepath);
[datpaths, recordingnames] = acqID(basepath, sortFiles, altSort);

otherdattypes = {'analogin'; 'digitalin'; 'auxiliary'; 'time'; 'supply'};
toFill = zeros(size(otherdattypes));
fileBase = cell(size(datpaths));

if fillMissingDatFiles
    % check if any of the folders to concatenate contain these files
    for i = 1:size(datpaths,2)
        curDat = datpaths{i};
        cutoff = find(curDat=='\', 1, 'last');
        curDat = curDat(2:cutoff); %remove formatting
        fileBase{i} = curDat;
        for j = 1:size(otherdattypes,1)
            if exist([curDat otherdattypes{j} '.dat'])
                toFill(j) = true;
            end
        end
    end
else
    datCount = zeros(size(datpaths));
    for i = 1:size(datpaths,2)
        curDat = datpaths{i};
        cutoff = find(curDat=='\', 1, 'last');
        curDat = curDat(2:cutoff); %remove formatting
        fileBase{i} = curDat;
        for j = 1:size(otherdattypes,1)
            if exist([curDat otherdattypes{j} '.dat'])
                datCount(j) = datCount(j)+1;
            end
        end
    end
    for j = 1:length(datCount)
        if datCount(j)<length(datpaths)
            toFill(j) = false;
            disp([otherdattypes{j} ' present in ' num2str(datCount(j))...
                '/' num2str(length(datpaths)) ' subfolders, will not concatenate.']);
            %If receiving this message, can make 'fillMissingDatFiles' true
            %to fill the subfolders missing the files with zero'd data
        else
            toFill(j) = true;
        end
    end
end
fileTypes = otherdattypes(logical(toFill));
for i = 1:size(fileTypes,1)
    newpaths.(fileTypes{i}) = fullfile(basepath, [fileTypes{i}, '.dat']);
end

%% Concatenate amplifier/continuous files
newpaths.amplifier = fullfile(basepath, [basename, '.dat']);
if ispc
    if size(datpaths,2) > 1
        for didx = 1:size(datpaths,2) - 1
            datpathsplus.amplifier{didx} = [datpaths{didx}, ' +'];
            t = dir(datpaths{didx}(2:end-2));
            datsizes.amplifier(didx) = t.bytes;
            for j = 1:size(fileTypes,1)
                datpathsplus.(fileTypes{j}){didx} = cat(2,'"',fileBase{didx},fileTypes{j},'.dat "', ' +');
                t = dir(datpathsplus.(fileTypes{j}){didx}(2:end-4));
                datsizes.(fileTypes{j}){didx} = t.bytes;
            end
        end
        t = dir(datpaths{end}(2:end-2));
        datsizes.amplifier(end+1) = t.bytes;
        %Last file string shouldn't end with '+'
        datpathsplus.amplifier{length(datpaths)} = [datpaths{length(datpaths)}];
        for j = 1:size(fileTypes,1)
            datpathsplus.(fileTypes{j}){length(datpaths)} = cat(2,'"',fileBase{didx},fileTypes{j},'.dat "');
            t = dir(datpathsplus.(fileTypes{j}){end}(2:end-2));
            datsizes.(fileTypes{j}){end+1} = t.bytes;
        end
    else
        datpathsplus.amplifier = datpaths;
        t = dir(datpaths{1}(2:end-2));
        datsizes.amplifier = t.bytes;
        for j = 1:size(fileTypes,1)
            datpathsplus.(fileTypes{j}) = cat(2,'"',fileBase{didx},fileTypes{j},'.dat "');
            t = dir(datpathsplus.(fileTypes{j}){1}(2:end-2));
            datsizes.(fileTypes{j}){1} = t.bytes;
        end
    end
    cs = strjoin(datpathsplus.amplifier);
    catstring.amplifier = ['! copy /b ', cs, ' ', newpaths.amplifier];
    for j = 1:size(fileTypes,1)
        cs = strjoin(datpathsplus.(fileTypes{j}));
        catstring.(fileTypes{j}) = ['! copy /b ', cs, ' ', newpaths.(fileTypes{j})];
    end
else
    error('This function cannot yet handle non-pc commands');
    %if isunix
    %cs = strjoin(datpaths.amplifier);
    %catstring = ['! cat ', cs, ' > ', newdatpath];
end

if ~exist(newpaths.amplifier, "file")
    disp('Concatenating Amplifier Dats... be patient');
    eval(catstring.amplifier) %execute concatention
else
    disp("Amplifier Dats already Concatenated")
end

for j = 1:size(fileTypes,1)
    if ~exist(newpaths.(fileTypes{j}), "file")
        disp(['Concatenating ' (fileTypes{j}) ' Dats... be patient']);
        eval(catstring.(fileTypes{j})) %execute concatention
    else
        disp([(fileTypes{j}) ' Dats already Concatenated']);
    end
end

t = dir(newpaths.amplifier);
if t.bytes ~= sum(datsizes.amplifier)
    error('New .dat size not right.  Exiting')
    return
else
    sizecheck.amplifier = true;
    disp('Primary .dats concatenated and size checked')
end

%% Create MergePoints file based on concatenation
sessionInfo = LoadXml(fullfile(basepath,[basename, '.xml']));
nSamp = [];
for didx = 1:length(datpaths)
   t = dir(datpaths{didx}(2:end-2)); 
   dataTypeNBytes = numel(typecast(cast(0, 'int16'), 'uint8'));
   nSamp(didx) = t.bytes/(sessionInfo.nChannels*dataTypeNBytes);
end

cumsum_nSamp = cumsum([nSamp]);
starts = [0,cumsum_nSamp(1:end-1)];
transitiontimes_samp = [starts',cumsum_nSamp'];
transitiontimes_sec = transitiontimes_samp./sessionInfo.SampleRate;

firstlasttimepoints = [zeros(length(nSamp),1),nSamp'];

MergePoints.timestamps = transitiontimes_sec;
MergePoints.timestamps_samples = transitiontimes_samp;
MergePoints.firstlasttimpoints_samples = firstlasttimepoints;
MergePoints.foldernames = recordingnames;
MergePoints.detectorinfo.detectorname = 'concatenateDatsNew';
MergePoints.detectorinfo.detectiondate = datestr(now,'yyyy-mm-dd');

save(fullfile(basepath,[basename,'.MergePoints.events.mat']),'MergePoints');
end