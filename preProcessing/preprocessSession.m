function preprocessSession(varargin)

% preprocessSession(varargin)

%   Master function to run the basic pre-processing pipeline for an
%   individual sessions. Is based on sessionsPipeline.m but in this case
%   works on an individual session basis no in a folfer with multiple ones
%
%
% INPUTS
%   input parser  inputs as opiton list of values - see below
%   <options>       optional list of property-value pairs (see table below)
%  =========================================================================
%   Properties    Values
%  -------------------------------------------------------------------------
% basepath                Basepath for experiment. It contains all session
%                         folders. If not provided takes pwd
% analogChannels          List of analog channels with pulses to be detected (it
%                         supports Intan Buzsaki Edition)
% digitalChannels         List of digital channels with pulses to be detected (it
%                         supports Intan Buzsaki Edition)
% forceSum                Force make folder summary (overwrite, if necessary).
%                         Default false
% cleanArtifacts          Remove artifacts from dat file. By default, if there is
%                         analogEv in folder, is true
% stateScore              Run automatic brain state detection with SleepScoreMaster.
%                         Default true
% LFPbeforeKilo           Option to generate LFP and run LFP-based
%                         functions before KiloSort. Default true
% spikeSort               Run automatic spike sorting using Kilosort. Default true
% sortFiles               Sort subsessions with the date and timestamp in the end
%                         of the folder name (ignore alphabetical order)
% removeNoise             Remove first (noise) ICA component from .dat file before
%                         spike sorting
% cleanRez                Run automatic noise detection of the Kilosort results
%                        (these will be pre-labelled as noise in phy). Default true
% getPos                  get tracking positions. Default false
% runSummary              run summary analysis using AnalysisBatchScript.
%                         Defualt false
% pullData                Path for raw data. Look for not analized session to
%                         copy to the main folder basepath. To do...
% path_to_dlc_bat_file    path to your dlc bat file to analyze your videos
%                         (see neurocode\behavior\dlc for example files)
% nKilosortRuns           Number of desired Kilosort runs (default = 1). The
%                         function will break down the shanks into "nKilosortRuns"
%                         groups for each run
% kiloShankSplit          Shanks delineation for multi-kilosort run. If you
%                         have a four shank probe and want to sort shanks 1
%                         and 2 together, skip 3, and sort 4 on it's own,
%                         input as: [1 1 0 2]. Shanks 1 and 2 will be the
%                         first kilosort run, and shank 4 will be the
%                         second. Shank 3 is set to 0 so it will not be
%                         run. Default [].
% sortFiles               Logical option to sort files by their Intan or
%                         OpenEphys timestamp. Setting to false will
%                         default to altSort ordering (below). If altSort
%                         is empty, files will be sorted alphabetically.
% altSort                 Numerical array of indices ordering your
%                         subsession dat files. If, for example, you have
%                         subsession folders containing dat files labeled
%                         as FolderA; FolderB; FolderC; and want the order
%                         to be concatenated as "C, A, B", input altSort as
%                         [2, 3, 1]; Default is false, which, when
%                         sortFiles is true, sorts by date/time
%                         (YYMMDD_HHMMSS for Intan, YYYY-MM-DD_HH-MM-SS for
%                         OpenEphys).
% ignoreFolders           Folder names that contain dat folders which
%                         should be ignored. Input should be a list of
%                         strings. Most often, this applies to a 'backup'
%                         folder containing original copies of the data.
%                         Example input may look like: ["backup",
%                         "ignore"].
%
%  OUTPUTS
%    N/A
%
%  NOTES
%   TODO:
%   - Improve auto-clustering routine
%   - Test cleaning and removing artifacts routines
%  SEE ALSO
%
% Copyright (C) 2020-2023 by AntonioFR, 2021 Azahara Oliva,
%               2022-2025 Lindsay Karaba, 2022-2024 Heath Larson,
%               2022-2024 Ryan Harvey, 2022-2024 Ralitsa Todorova
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Set up parameters and parse inputs

p = inputParser;
addParameter(p, 'basepath', pwd, @isfolder); % by default, current folder
addParameter(p, 'fillMissingDatFiles', false, @islogical);
addParameter(p, 'analogInputs', false, @islogical);
addParameter(p, 'analogChannels', [], @isnumeric);
addParameter(p, 'digitalInputs', false, @islogical);
addParameter(p, 'digitalChannels', [], @isnumeric);
addParameter(p, 'getAcceleration', false, @islogical);
addParameter(p, 'cleanArtifacts', false, @islogical);
addParameter(p, 'stateScore', true, @islogical);
addParameter(p, 'LFPbeforeKilo', true, @islogical);
addParameter(p, 'spikeSort', true, @islogical);
addParameter(p, 'cleanRez', true, @islogical);
addParameter(p, 'getPos', false, @islogical);
addParameter(p, 'removeNoise', false, @islogical); % denoising method removing the first PCA component
addParameter(p, 'runSummary', false, @islogical);
addParameter(p, 'SSD_path', 'D:\KiloSort', @ischar); % Path to SSD disk. Make it empty to disable SSD
addParameter(p, 'path_to_dlc_bat_file', '', @isfile)
addParameter(p, 'nKilosortRuns', 1, @isnumeric);
addParameter(p, 'kiloShankSplit', [], @isnumeric);
addParameter(p, 'sortFiles', true, @islogical);
addParameter(p, 'altSort', [], @isnumeric);
addParameter(p, 'ignoreFolders', "", @isstring);
addParameter(p, 'SWChannels', 0, @isnumeric);
addParameter(p, 'ThetaChannels', 0, @isnumeric);
addParameter(p, 'clean_rez_params', { ...
    'mahalThreshold', inf, ...
    'minNumberOfSpikes', 100, ...
    'multiTrough', false, ...
    'isi', false, ...
    'singleBin', false, ...
    'global', false, ...
    }, @iscell);

% addParameter(p,'pullData',[],@isdir); To do...
parse(p, varargin{:});

basepath = p.Results.basepath;
fillMissingDatFiles = p.Results.fillMissingDatFiles;
analogInputs = p.Results.analogInputs;
analogChannels = p.Results.analogChannels;
digitalInputs = p.Results.digitalInputs;
digitalChannels = p.Results.digitalChannels;
getAcceleration = p.Results.getAcceleration;
cleanArtifacts = p.Results.cleanArtifacts;
stateScore = p.Results.stateScore;
LFPbeforeKilo = p.Results.LFPbeforeKilo;
spikeSort = p.Results.spikeSort;
cleanRez = p.Results.cleanRez;
getPos = p.Results.getPos;
removeNoise = p.Results.removeNoise;
runSummary = p.Results.runSummary;
SSD_path = p.Results.SSD_path;
path_to_dlc_bat_file = p.Results.path_to_dlc_bat_file;
nKilosortRuns = p.Results.nKilosortRuns;
kiloShankSplit = p.Results.kiloShankSplit;
sortFiles = p.Results.sortFiles;
altSort = p.Results.altSort;
ignoreFolders = p.Results.ignoreFolders;
clean_rez_params = p.Results.clean_rez_params;
SWChannels = p.Results.SWChannels;
ThetaChannels = p.Results.ThetaChannels;

if removeNoise
    removeNoise = false;
    warning('removeNoise not implemented')
end

% Check for active breakpoints
dbstatus_ = dbstatus('-completenames');
if ~isempty(dbstatus_)
    warning('Active breakpoints detected! Consider removing breakpoints (dbclear all)');
end

if ~exist(basepath, 'dir')
    error('path provided does not exist')
end
cd(basepath)

%% Pull meta data

% Get session names
if strcmp(basepath(end), filesep)
    basepath = basepath(1:end-1);
end
[~, basename] = fileparts(basepath);

% Get xml file in order
if ~exist([basepath, '\', basename, '.xml'], 'file')
    xmlFile = checkFile('fileType', '.xml', 'searchSubdirs', true);
    xmlFile = xmlFile(1);
    if ~(strcmp(xmlFile.folder, basepath) && strcmp(xmlFile.name(1:end-4), basename))
        copyfile([xmlFile.folder, filesep, xmlFile.name], [basepath, filesep, basename, '.xml'])
    end
end

% Check info.rhd
% (assumes this will be the same across subsessions)
rhdFile = dir("**/*.rhd");
if ~isempty(rhdFile)
    rhdFile = rhdFile(1);
    if ~(strcmp(rhdFile.folder, basepath) && strcmp(rhdFile.name(1:end-4), basename))
        copyfile([rhdFile.folder, filesep, rhdFile.name], [basepath, filesep, basename, '.rhd'])
    end
else
    disp('No rhd file found. This will be the case if only processing openEphys files. Skipping step');
end

%% Make SessionInfo
% Manually ID bad channels at this point. automating it would be good
session = sessionTemplate(basepath, 'showGUI', false);
save(fullfile(basepath, [basename, '.session.mat']), 'session');

%% Concatenate sessions
disp('Concatenate session folders...');
concatenateDats('basepath', basepath, 'fillMissingDatFiles', fillMissingDatFiles, ...
    'sortFiles', sortFiles, 'altSort', altSort, 'ignoreFolders', ignoreFolders);

%% run again to add epochs from basename.MergePoints.mat
session = sessionTemplate(basepath, 'showGUI', false);
save(fullfile(basepath, [basename, '.session.mat']), 'session');

%% Process additional inputs - CHECK FOR OUR LAB

% Analog inputs
% check the two different fucntions for delaing with analog inputs and proably rename them
if analogInputs
    if ~isempty(analogChannels)
        analogInp = computeAnalogInputs('analogCh', analogChannels, 'saveMat', true);
    else
        analogInp = computeAnalogInputs('analogCh', [], 'saveMat', true);
    end

    % analog pulses ...
    [pulses] = getAnalogPulses('samplingRate', session.extracellular.sr);
end

% Digital inputs
if digitalInputs
    if ~isempty(digitalChannels)
        % need to change to only include specified channels
        getDigitalIn('all', 'fs', session.extracellular.sr, 'digUse', digitalChannels);
    else
        getDigitalIn('all', 'fs', session.extracellular.sr);
    end
end

% Auxilary input
if getAcceleration
    computeIntanAccel('saveMat', true);
end

%% Clean data  - CHECK FOR OUR LAB
% Remove stimulation artifacts
if cleanArtifacts && analogInputs
    [pulses] = getAnalogPulses(analogInp, 'analogCh', analogChannels);
    cleanPulses(pulses.ints{1}(:));
end

if LFPbeforeKilo
    %Make LFP
    runLFP(basepath, basename, session);

    % Get brain states
    % an automatic way of flaging bad channels is needed
    if stateScore
        if ~exist('pulses', 'var')
            pulses = [];
        end
        runStateScore(basepath, pulses, session, SWChannels, ThetaChannels);
    end

    % remove noise from data for cleaner spike sorting
    if removeNoise
        runRemoveNoise(basepath, basename, session);
    end

    %% Kilosort concatenated sessions - Needs to be changed to probes, not shanks HLR 01/05/2023
    if spikeSort
        runKiloSort(nKilosortRuns, kiloShankSplit, session, SSD_path, clean_rez_params, cleanRez);
    end
else
    % remove noise from data for cleaner spike sorting
    if removeNoise
        runRemoveNoise(basepath, basename, session);
    end

    % Kilosort concatenated sessions - Needs to be changed to probes, not shanks HLR 01/05/2023
    if spikeSort
        runKiloSort(nKilosortRuns, kiloShankSplit, session, SSD_path, clean_rez_params, cleanRez);
    end

    % Make LFP
    runLFP(basepath, basename, session);

    % Get brain states
    % an automatic way of flaging bad channels is needed
    if stateScore
        if ~exist('pulses', 'var')
            pulses = [];
        end
        runStateScore(basepath, pulses, session, SWChannels, ThetaChannels);
    end
end

%% Get tracking positions
if getPos
    % check for pre existing deeplab cut
    if exist(path_to_dlc_bat_file, 'file')
        system(path_to_dlc_bat_file, '-echo')
    end
    % put tracking into standard format
    general_behavior_file('basepath', basepath)
end

%% Summary -  NOT WELL IMPLEMENTED YET
if runSummary
    if forceSum || (~isempty(dir('*Kilosort*')) && (isempty(dir('summ*')))) % is kilosorted but no summ
        disp('running summary analysis...');
        sessionSummary;
    end
end

%% logging
% log params used
results = p.Results;
save(fullfile(basepath, 'preprocessSession_params.mat'), 'results')

% log script
% saves a text file of the current code used
targetFile = fullfile(basepath, 'preprocessSession.log');
copyfile(which('preprocessSession.m'), targetFile);

end

function runLFP(basepath, basename, session)
try
    try
        LFPfromDat(basepath, 'outFs', 1250, 'useGPU', true);
    catch e
        fprintf(1, 'The identifier was:\n%s', e.identifier);
        fprintf(1, 'There was an error! The message was:\n%s', e.message);
        if (exist([basepath, '\', basename, '.lfp'], "file") ~= 0)
            fclose([basepath, '\', basename, '.lfp']); %if the above run failed after starting the file
            delete([basepath, '\', basename, '.lfp']);
        end
        LFPfromDat(basepath, 'outFs', 1250, 'useGPU', true);
    end
catch e
    fprintf(1, 'The identifier was:\n%s', e.identifier);
    fprintf(1, 'There was an error! The message was:\n%s', e.message);
    try
        warning('LFPfromDat failed, trying ResampleBinary')
        ResampleBinary([basepath, '\', basename, '.dat'], ...
            [basepath, '\', basename, '.lfp'], session.extracellular.nChannels, 1, 16);
    catch e
        warning('LFP file could not be generated, moving on');
        fprintf(1, 'The identifier was:\n%s', e.identifier);
        fprintf(1, 'There was an error! The message was:\n%s', e.message);
    end
end
% 'useGPU'=true gives an error if CellExplorer in the path. Need to test if
% it is possible to remove the copy of iosr toolbox from CellExplorer -
% seems to be fixed? as of 9/22
end

function runStateScore(basepath, pulses, session, SWChannels, ThetaChannels)
try
    if ~isempty(pulses)
        SleepScoreMaster(basepath, 'noPrompts', true, ...
            'ignoretime', pulses.intsPeriods, ...
            'rejectChannels', session.channelTags.Bad.channels, ...
            'SWChannels', SWChannels, ...
            'ThetaChannels', ThetaChannels);
    else
        SleepScoreMaster(basepath, 'noPrompts', true, ...
            'rejectChannels', session.channelTags.Bad.channels, ...
            'SWChannels', SWChannels, ...
            'ThetaChannels', ThetaChannels);
    end
catch e
    warning('Problem with SleepScore scoring... unable to calculate');
    fprintf(1, 'The identifier was:\n%s', e.identifier);
    fprintf(1, 'There was an error! The message was:\n%s', e.message);
end
end

function runRemoveNoise(basepath, basename, session)

try
    EMGFromLFP = getStruct(basepath, 'EMGFromLFP');
catch e
    fprintf(1, 'The identifier was:\n%s', e.identifier);
    fprintf(1, 'There was an error! The message was:\n%s', e.message);

    EMGFromLFP = getEMGFromLFP(basepath, 'noPrompts', true, 'saveMat', true);
end

baseline = EMGFromLFP.timestamps(FindInterval(EMGFromLFP.data > quantile(EMGFromLFP.data, 0.99))); % select the period of top 1% EMG activity as the denoising baseline
DenoiseDat(fullfile(basepath, [basename, '.dat']), session, 'baseline', baseline);
end


function runKiloSort(nKilosortRuns, kiloShankSplit, session, SSD_path, clean_rez_params, cleanRez)
if nKilosortRuns > 1 % if more than one Kilosort cycle desired, break the shanks down into the desired number of kilosort runs
    shanks = session.extracellular.spikeGroups.channels;
    if isempty(kiloShankSplit)
        kilosortGroup = ceil(((1:length(shanks)) / nKilosortRuns));
    else
        kilosortGroup = kiloShankSplit;
        nKilosortRuns = max(kiloShankSplit);
    end
    for i = 1:nKilosortRuns
        channels = cat(2, shanks{kilosortGroup == i});
        excludeChannels = find(~ismember((1:session.extracellular.nChannels), channels));
        excludeChannels = cat(2, excludeChannels, session.channelTags.Bad.channels);
        excludeChannels = unique(excludeChannels);
        if (length(excludeChannels) == session.extracellular.nChannels)
            warning(['Run number ', num2str(i), ' excluded, moving on']);
        else
            kilosortFolder = KiloSortWrapper('SSD_path', SSD_path, 'rejectchannels', excludeChannels);
            if cleanRez
                load(fullfile(kilosortFolder, 'rez.mat'), 'rez');
                CleanRez(rez, 'savepath', kilosortFolder, clean_rez_params{:});
            end
        end
    end
else

    %% single sort
    kilosortFolder = KiloSortWrapper('SSD_path', SSD_path, ...
        'rejectchannels', session.channelTags.Bad.channels); % 'NT',20*1024 for long sessions when RAM is overloaded
    if cleanRez
        load(fullfile(kilosortFolder, 'rez.mat'), 'rez');
        CleanRez(rez, 'savepath', kilosortFolder, clean_rez_params{:});
    end

end
end