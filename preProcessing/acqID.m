function [datpaths, recordingnames] = acqID(basepath, sortFiles, altSort)

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

%% OUTPUTS
% datpaths                Returns cell array of ordered datpaths in "PATH "
%                         format, ready for concatenation.
% recordingnames          Returns cell array of ordered subsession folder
%                         names.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if sortFiles && (~isempty(altSort))
    error('sortFiles cannot be empty while altSort provides an order. Please choose to either sort by time (sortFiles=true) or designate a manual order of concatenation (altSort)');
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
for i = 1:size(useIDX, 1)
    datpaths{i} = cat(2, '"', allFolders(useIDX(i), :).folder{1}, '\', allFolders(useIDX(i), :).name{1}, ' "');
    startIDX = size(basepath, 2) + 2;
    seps = find(allFolders(useIDX(i), :).folder{1} == '\');
    seps(seps < startIDX) = [];
    if isempty(seps)
        recordingnames{i} = allFolders(useIDX(i), :).folder{1}(startIDX:end); %intan
        expNum(i) = '1';
        recNum(i) = '1';
    else
        recordingnames{i} = allFolders(useIDX(i), :).folder{1}(startIDX:seps(1) - 1); %openEphys
        expIDX = strfind(allFolders(useIDX(i), :).folder{1}, "\experiment");
        expNum(i) = (allFolders(useIDX(i), :).folder{1}(expIDX + 11:seps(find(seps > expIDX, 1, 'first')) - 1));
        recIDX = strfind(allFolders(useIDX(i), :).folder{1}, "\recording");
        recNum(i) = (allFolders(useIDX(i), :).folder{1}(recIDX + 10:seps(find(seps > recIDX, 1, 'first')) - 1));
    end
end

% datpaths and recordingnames are ordered alphabetically by default. If
% not maintaining the default order (ie sorting by time or manual order
% rather than by folder name), this loop will be entered
if sortFiles %if we sort by time
    names2sort_intan = cellfun(@(x) regexpi(x, '(?<=_)\d{6}', 'match'), recordingnames, 'UniformOutput', false);
    ordering = nan(1, size(names2sort_intan, 2));
    for i = 1:size(names2sort_intan, 2)
        if isempty(names2sort_intan{i})
            %assume open ephys
            ordering(i) = str2num([recordingnames{i}([3:4, 6:7, 9:10, 12:13, 15:16, 18:19]), expNum(i), recNum(i)]);
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