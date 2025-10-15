function concatenateDats(varargin)
% Concatenate raw .dat files found in a session folder
% for Intan AND/OR openEphys recordings
%
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
% ignoreFolders           Folder names that contain dat folders which
%                         should be ignored. Input should be a list of
%                         strings. Most often, this applies to a 'backup'
%                         folder containing original copies of the data.
%                         Example input may look like: ["backup",
%                         "ignore"].
% behaviorOnly            Set to true for behavior-only sessions (no
%                         amplifier.dat). Reads metadata from info.rhd or
%                         settings.xml. Skips LFP, Kilosort, and SleepScore.
%                         Default false

p = inputParser;
addParameter(p, 'basepath', pwd, @isfolder); % by default, current folder
addParameter(p, 'fillMissingDatFiles', false, @islogical);
addParameter(p, 'sortFiles', true, @islogical);
addParameter(p, 'altSort', [], @isnumeric);
addParameter(p, 'ignoreFolders', "", @isstring);
addParameter(p, 'behaviorOnly', false, @islogical);

parse(p, varargin{:});

basepath = p.Results.basepath;
fillMissingDatFiles = p.Results.fillMissingDatFiles;
sortFiles = p.Results.sortFiles;
altSort = p.Results.altSort;
ignoreFolders = p.Results.ignoreFolders;
behaviorOnly = p.Results.behaviorOnly;

if sortFiles && (~isempty(altSort))
    error('sortFiles cannot be empty while altSort provides an order. Please choose to either sort by time (sortFiles=true) or designate a manual order of concatenation (altSort)');
end

basename = basenameFromBasepath(basepath);
[datpaths, recordingnames] = acqID(basepath, sortFiles, altSort, ignoreFolders, behaviorOnly);
if isempty(datpaths)
    disp('no subsessions detected, exiting concatenation');
    return
else
    disp('Concatenating subsessions in the following order: ');
    for i = 1:size(datpaths, 2)
        disp(datpaths{i}(2:end - 2)); %don't print added formating for concatentation
    end
end

if behaviorOnly
    % In behavior-only mode, digitalin is the primary file (replaces amplifier)
    % so we don't include it in otherdattypes
    otherdattypes = {'analogin'; 'auxiliary'; 'time'; 'supply'};
else
    % Standard mode: amplifier is primary, digitalin is in otherdattypes
    otherdattypes = {'analogin'; 'digitalin'; 'auxiliary'; 'time'; 'supply'};
end
toFill = zeros(size(otherdattypes));
fileBase = cell(size(datpaths));

if fillMissingDatFiles
    % check if any of the folders to concatenate contain these files
    for i = 1:size(datpaths, 2)
        curDat = datpaths{i};
        cutoff = find(curDat == filesep, 1, 'last');
        curDat = curDat(2:cutoff); %remove formatting
        fileBase{i} = curDat;
        for j = 1:size(otherdattypes, 1)
            if exist([curDat, otherdattypes{j}, '.dat'], 'file')
                toFill(j) = true;
            end
        end
    end
    for ii = 1:length(toFill)
        if toFill(ii) == 1
            fillMissingDats('basepath', basepath, 'fileType', otherdattypes{ii}, 'behaviorOnly', behaviorOnly);
        end
    end
else
    datCount = zeros(size(otherdattypes));
    for i = 1:size(datpaths, 2)
        curDat = datpaths{i};
        cutoff = find(curDat == filesep, 1, 'last');
        curDat = curDat(2:cutoff); %remove formatting
        fileBase{i} = curDat;
        for j = 1:size(otherdattypes, 1)
            if exist([curDat, otherdattypes{j}, '.dat'], 'file')
                datCount(j) = datCount(j) + 1;
            end
        end
    end
    for j = 1:length(datCount)
        if datCount(j) < length(datpaths)
            toFill(j) = false;
            disp([otherdattypes{j}, ' present in ', num2str(datCount(j)), ...
                filesep, num2str(length(datpaths)), ' subfolders, will not concatenate.']);
            %If receiving this message, can make 'fillMissingDatFiles' true
            %to fill the subfolders missing the files with zero'd data
        else
            toFill(j) = true;
        end
    end
end
fileTypes = otherdattypes(logical(toFill));
for i = 1:size(fileTypes, 1)
    newpaths.(fileTypes{i}) = fullfile(basepath, [fileTypes{i}, '.dat']);
end

%% Concatenate amplifier/continuous files
newpaths.amplifier = fullfile(basepath, [basename, '.dat']);
if ispc
    if size(datpaths, 2) > 1
        for didx = 1:size(datpaths, 2) - 1
            datpathsplus.amplifier{didx} = [datpaths{didx}, ' +'];
            t = dir(datpaths{didx}(2:end - 2));
            datsizes.amplifier(didx) = t.bytes;
            for j = 1:size(fileTypes, 1)
                datpathsplus.(fileTypes{j}){didx} = cat(2, '"', fileBase{didx}, fileTypes{j}, '.dat "', ' +');
                t = dir(datpathsplus.(fileTypes{j}){didx}(2:end - 4));
                datsizes.(fileTypes{j}){didx} = t.bytes;
            end
        end
        t = dir(datpaths{end}(2:end - 2));
        datsizes.amplifier(end+1) = t.bytes;
        %Last file string shouldn't end with '+'
        datpathsplus.amplifier{length(datpaths)} = [datpaths{length(datpaths)}];
        for j = 1:size(fileTypes, 1)
            datpathsplus.(fileTypes{j}){length(datpaths)} = cat(2, '"', fileBase{end}, fileTypes{j}, '.dat "');
            t = dir(datpathsplus.(fileTypes{j}){end}(2:end - 2));
            datsizes.(fileTypes{j}){end, +1} = t.bytes;
        end
    else
        datpathsplus.amplifier = datpaths;
        t = dir(datpaths{1}(2:end - 2));
        datsizes.amplifier = t.bytes;
        for j = 1:size(fileTypes, 1)
            datpathsplus.(fileTypes{j}){1} = cat(2, '"', fileBase{1}, fileTypes{j}, '.dat "');
            t = dir(datpathsplus.(fileTypes{j}){1}(2:end - 2));
            datsizes.(fileTypes{j}){1} = t.bytes;
        end
    end
    cs = strjoin(datpathsplus.amplifier);
    catstring.amplifier = ['! copy /b ', cs, ' ', newpaths.amplifier];
    for j = 1:size(fileTypes, 1)
        cs = strjoin(datpathsplus.(fileTypes{j}));
        catstring.(fileTypes{j}) = ['! copy /b ', cs, ' ', newpaths.(fileTypes{j})];
    end
else
    % Linux/Unix: build input lists and concatenate with cat
    nSubs = size(datpaths, 2);
    datsizes.amplifier = zeros(1, nSubs);
    datpathsplus.amplifier = cell(1, nSubs);

    for didx = 1:nSubs
        % Unquote original entry (stored like: "full/path/file.dat ")
        amp_in = datpaths{didx}(2:end-2);

        % Resolve to exactly one file and record size
        info = dir(amp_in);                          % may be scalar, empty, or multiple [web:39]
        if isempty(info)
            error('File not found: %s', amp_in);
        end
        info = info(~[info.isdir]);                  % drop directories [web:39]
        [~,targetName,targetExt] = fileparts(amp_in);% parse target name/ext [web:43]
        match = strcmp({info.name}, [targetName targetExt]); % exact match
        if any(match)
            info = info(find(match,1,'first'));
        elseif numel(info) == 1
            % Single non-dir entry; accept as best-effort
            info = info(1);
        else
            error('Path did not resolve to a single file: %s', amp_in);
        end

        datpathsplus.amplifier{didx} = ['"', amp_in, '"'];
        datsizes.amplifier(didx) = info.bytes;

        % Build per-type inputs
        for j = 1:size(fileTypes, 1)
            in_path = [fileBase{didx}, fileTypes{j}, '.dat'];
            datpathsplus.(fileTypes{j}){didx} = ['"', in_path, '"'];
            % Optional bookkeeping (not used in final size check for amplifier)
            info2 = dir(in_path);                    % may be empty if not selected for concat [web:39]
            if ~isempty(info2), info2 = info2(~[info2.isdir]); end
            if ~isempty(info2)
                datsizes.(fileTypes{j}){didx} = info2(1).bytes; %#ok<AGROW>
            else
                datsizes.(fileTypes{j}){didx} = 0; %#ok<AGROW>
            end
        end
    end

    % Compose shell commands with proper quoting and redirection
    cs = strjoin(datpathsplus.amplifier, ' ');
    catstring.amplifier = ['! cat ', cs, ' > "', newpaths.amplifier, '"'];

    for j = 1:size(fileTypes, 1)
        cs = strjoin(datpathsplus.(fileTypes{j}), ' ');
        catstring.(fileTypes{j}) = ['! cat ', cs, ' > "', newpaths.(fileTypes{j}), '"'];
    end
end

if ~exist(newpaths.amplifier, "file")
    disp('Concatenating Amplifier Dats... be patient');
    eval(catstring.amplifier) %execute concatention
else
    disp("Amplifier Dats already Concatenated")
end

for j = 1:size(fileTypes, 1)
    if ~exist(newpaths.(fileTypes{j}), "file")
        disp(['Concatenating ', (fileTypes{j}), ' Dats... be patient']);
        eval(catstring.(fileTypes{j})) %execute concatention
    else
        disp([(fileTypes{j}), ' Dats already Concatenated']);
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
% For behavior-only sessions, read metadata from info.rhd or session file
% instead of requiring XML file
if behaviorOnly
    % Try to load session file to get sampling rate
    sessionFile = fullfile(basepath, [basename, '.session.mat']);
    if exist(sessionFile, 'file')
        load(sessionFile, 'session');
        sessionInfo.SampleRate = session.extracellular.sr;
        sessionInfo.nChannels = 1; % For digitalin.dat, we use bytes directly
    else
        % Try to read from info.rhd
        rhdFile = dir(fullfile(basepath, '**', '*.rhd'));
        if ~isempty(rhdFile)
            try
                [~, ~, ~, ~, ~, ~, frequency_parameters, ~] = ...
                    read_Intan_RHD2000_file_bz('basepath', rhdFile(1).folder);
                sessionInfo.SampleRate = frequency_parameters.board_dig_in_sample_rate;
                sessionInfo.nChannels = 1; % For digitalin.dat
            catch
                warning('Could not read sampling rate from info.rhd. Using default 20000 Hz');
                sessionInfo.SampleRate = 20000;
                sessionInfo.nChannels = 1;
            end
        else
            warning('No session file or info.rhd found. Using default sampling rate 20000 Hz');
            sessionInfo.SampleRate = 20000;
            sessionInfo.nChannels = 1;
        end
    end
    
    % For behavior-only, calculate samples from digitalin.dat (uint16, 2 bytes per sample)
    nSamp = [];
    for didx = 1:length(datpaths)
        % Find digitalin.dat in this subsession folder
        subsessionPath = datpaths{didx}(2:end - 2); % Remove quotes and spaces
        [parentPath, ~, ~] = fileparts(subsessionPath);
        digitalinPath = fullfile(parentPath, 'digitalin.dat');
        
        if exist(digitalinPath, 'file')
            t = dir(digitalinPath);
            % digitalin.dat: uint16, 2 bytes per sample
            nSamp(didx) = t.bytes / 2;
        else
            error(['digitalin.dat not found in ', parentPath]);
        end
    end
else
    % Standard mode: use XML file
    sessionInfo = LoadXml(fullfile(basepath, [basename, '.xml']));
    nSamp = [];
    for didx = 1:length(datpaths)
        t = dir(datpaths{didx}(2:end - 2));
        dataTypeNBytes = numel(typecast(cast(0, 'int16'), 'uint8'));
        nSamp(didx) = t.bytes / (sessionInfo.nChannels * dataTypeNBytes);
    end
end

cumsum_nSamp = cumsum(nSamp);
starts = [0, cumsum_nSamp(1:end-1)];
transitiontimes_samp = [starts', cumsum_nSamp'];
transitiontimes_sec = transitiontimes_samp ./ sessionInfo.SampleRate;

firstlasttimepoints = [zeros(length(nSamp), 1), nSamp'];

MergePoints.timestamps = transitiontimes_sec;
MergePoints.timestamps_samples = transitiontimes_samp;
MergePoints.firstlasttimpoints_samples = firstlasttimepoints;
MergePoints.foldernames = recordingnames;
MergePoints.detectorinfo.detectorname = 'concatenateDatsNew';
MergePoints.detectorinfo.detectiondate = datestr(now, 'yyyy-mm-dd');

save(fullfile(basepath, [basename, '.MergePoints.events.mat']), 'MergePoints');
end