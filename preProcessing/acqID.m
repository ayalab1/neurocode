function [datpaths, recordingnames] = acqID(basepath, sortFiles, altSort, ignoreFolders)

%% [datpaths, recordingnames] = acqID(basepath, sortFiles, altSort)

% Function to acquire the ID of dat files (Intan or openEphys) within the
% given session folder (basepath) and sort according to the indicated
% paradigm (sortFiles, altSort). The function returns the path to and names
% of the found dat files in sorted order. This script is to be used in
% conjunction with concatenateDats.m

%% INPUTS
% basepath                Basepath for experiment. It contains all session
%                         folders. If not provided takes pwd.
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

%% OUTPUTS
% datpaths                Returns cell array of ordered datpaths in "PATH "
%                         format, ready for concatenation.
% recordingnames          Returns cell array of ordered subsession folder
%                         names.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if sortFiles && (~isempty(altSort))
    error('sortFiles cannot be empty while altSort provides an order. Please choose to either sort by time (sortFiles=true) or designate a manual order of concatenation (altSort)');
end
if ignoreFolders == ""
    ignoreFolders = []; %set to empty if default
end

d = dir(basepath);
datpaths = {};
recordingnames = {};
expNum = [];
recNum = [];

listing = dir("**");
allFolders = struct2table(listing);
useIDX = find(contains(allFolders.name, 'continuous.dat') | contains(allFolders.name, 'amplifier.dat'));
folderNames = allFolders(useIDX, :).folder;

%this assumes that your recording folders with the date are within your
%current basepath directory
removeID = [];
for i = 1:size(useIDX, 1)
    checkPath = cat(2, '"', allFolders(useIDX(i), :).folder{1}, filesep, allFolders(useIDX(i), :).name{1}, ' "');
    usePath = true;
    for f = 1:length(ignoreFolders)
        if contains(checkPath, ignoreFolders(f))
            usePath = false;
        end
    end
    if usePath
        datpaths{i} = checkPath;
        folderPath = allFolders(useIDX(i), :).folder{1};

        % Remove basepath from folderPath
        relPath = strrep(folderPath, basepath, '');
        relPath = regexprep(relPath, ['^' filesep], ''); % Remove leading separator

        % Split relative path into parts
        parts = split(relPath, filesep);

        if numel(parts) == 1
            % Intan format: only one folder after basepath
            recordingnames{i} = parts{1};
            expNum(i) = '1';
            recNum(i) = '1';
        else
            recordingnames{i} = relParts{1};
            expIdx = find(startsWith(relParts, 'experiment'), 1);
            recIdx = find(startsWith(relParts, 'recording'), 1);
            expNum(i) = ternary(~isempty(expIdx), extractAfter(relParts{expIdx}, 'experiment'), '1');
            recNum(i) = ternary(~isempty(recIdx), extractAfter(relParts{recIdx}, 'recording'), '1');
        end
    else
        fprintf('.dat file found nested in a folder labeled %s . Skipping: \n', ignoreFolders(f));
        disp(checkPath);
        if i ~= numel(useIDX)  % won't need to remove if a real slot isn't filled after
            removeID(end+1) = i;
        end
    end
end
if ~isempty(removeID)
    datpaths(removeID) = [];
    recordingnames(removeID) = [];
    expNum(removeID) = [];
    recNum(removeID) = [];
end

% datpaths and recordingnames are ordered alphabetically by default. If
% not maintaining the default order (ie sorting by time or manual order
% rather than by folder name), this loop will be entered
if sortFiles %if we sort by time
    names2sort_intan = cellfun(@(x) regexpi(x, '(?<=_)\d{6}', 'match'), recordingnames, 'UniformOutput', false);
    dateStart_OE = cellfun(@(x) regexpi(x, '\d{4}\-\d{2}\-\d{2}\_\d{2}\-\d{2}\-\d{2}', 'match'), recordingnames, 'UniformOutput', false);
    ordering = nan(1, size(names2sort_intan, 2));
    for i = 1:size(names2sort_intan, 2)
        if isempty(names2sort_intan{i})
            %assume open ephys
            find_date_start = strfind(recordingnames{i}, dateStart_OE{i});
            ordering(i) = str2num([recordingnames{i}(find_date_start + [2:3, 5:6, 8:9, 11:12, 14:15, 17:18]), expNum(i), recNum(i)]);
        else
            ordering(i) = str2num([names2sort_intan{i}{1}, names2sort_intan{i}{2}, expNum(i), recNum(i)]);
        end
    end

    [~, order] = sort(ordering);
    datpaths = datpaths(order);
    recordingnames = recordingnames(order);
elseif ~sortFiles && ~isempty(altSort) %if we sort by manual ordering
    datpaths = datpaths(altSort);
    recordingnames = recordingnames(altSort);
end
end