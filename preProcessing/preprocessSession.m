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
% [sortFiles]             Concatenate .dat files based on intan time or
%                         alphabetically (default = true, based on time).
%                         Alternatively, enter false for alphabetical sort.
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
%               2022-2023 Lindsay Karaba, 2022-2024 Heath Larson,
%               2022-2024 Ryan Harvey, 2022-2024 Ralitsa Todorova
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Set up parameters and parse inputs

p = inputParser;
addParameter(p, 'basepath', pwd, @isfolder); % by default, current folder
addParameter(p, 'fillMissingDatFiles', false, @islogical);
addParameter(p, 'fillTypes', [], @iscellstr);
addParameter(p, 'analogInputs', false, @islogical);
addParameter(p, 'analogChannels', [], @isnumeric);
addParameter(p, 'digitalInputs', false, @islogical);
addParameter(p, 'digitalChannels', [], @isnumeric);
addParameter(p, 'getAcceleration', false, @islogical);
addParameter(p, 'cleanArtifacts', false, @islogical);
addParameter(p, 'stateScore', true, @islogical);
addParameter(p, 'spikeSort', true, @islogical);
addParameter(p, 'cleanRez', true, @islogical);
addParameter(p, 'getPos', false, @islogical);
addParameter(p, 'removeNoise', false, @islogical); % denoising method removing the first PCA component
addParameter(p, 'runSummary', false, @islogical);
addParameter(p, 'SSD_path', 'D:\KiloSort', @ischar) % Path to SSD disk. Make it empty to disable SSD
addParameter(p, 'path_to_dlc_bat_file', '', @isfile)
addParameter(p, 'nKilosortRuns', 1, @isnumeric);
addParameter(p, 'sortFiles', true, @islogical);
addParameter(p, 'clean_rez_params', { ...
    'mahalThreshold', 12, ...
    'minNumberOfSpikes', 20, ...
    'multiTrough', true, ...
    'isi', true, ...
    'singleBin', true, ...
    'global', true, ...
    }, @iscell);


% addParameter(p,'pullData',[],@isdir); To do...
parse(p, varargin{:});

basepath = p.Results.basepath;
fillMissingDatFiles = p.Results.fillMissingDatFiles;
fillTypes = p.Results.fillTypes;
analogInputs = p.Results.analogInputs;
analogChannels = p.Results.analogChannels;
digitalInputs = p.Results.digitalInputs;
digitalChannels = p.Results.digitalChannels;
getAcceleration = p.Results.getAcceleration;
cleanArtifacts = p.Results.cleanArtifacts;
stateScore = p.Results.stateScore;
spikeSort = p.Results.spikeSort;
cleanRez = p.Results.cleanRez;
getPos = p.Results.getPos;
removeNoise = p.Results.removeNoise;
runSummary = p.Results.runSummary;
SSD_path = p.Results.SSD_path;
path_to_dlc_bat_file = p.Results.path_to_dlc_bat_file;
nKilosortRuns = p.Results.nKilosortRuns;
sortFiles = p.Results.sortFiles;
clean_rez_params = p.Results.clean_rez_params;


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
if ~exist([basepath '\' basename '.xml'], 'file')
    xmlFile = checkFile('fileType', '.xml', 'searchSubdirs', true);
    xmlFile = xmlFile(1);
    if ~(strcmp(xmlFile.folder, basepath) && strcmp(xmlFile.name(1:end-4), basename))
        copyfile([xmlFile.folder, filesep, xmlFile.name], [basepath, filesep, basename, '.xml'])
    end
end

% Check info.rhd
% (assumes this will be the same across subsessions)
rhdFile = checkFile('fileType', '.rhd', 'searchSubdirs', true);
rhdFile = rhdFile(1);
if ~(strcmp(rhdFile.folder, basepath) && strcmp(rhdFile.name(1:end-4), basename))
    copyfile([rhdFile.folder, filesep, rhdFile.name], [basepath, filesep, basename, '.rhd'])
end

%% Make SessionInfo
% Manually ID bad channels at this point. automating it would be good
session = sessionTemplate(basepath, 'showGUI', false);
save(fullfile(basepath, [basename, '.session.mat']), 'session');

%% Fill missing dat files of zeros
if fillMissingDatFiles
    if isempty(fillTypes)
        fillTypes = {'analogin'; 'digitalin'; 'auxiliary'; 'time'; 'supply'};
    end
    for ii = 1:length(fillTypes)
        fillMissingDats('basepath', basepath, 'fileType', fillTypes{ii});
    end
end

%% Concatenate sessions
disp('Concatenate session folders...');
concatenateDats(basepath, sortFiles);

%% run again to add epochs from basename.MergePoints.m
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

%% Make LFP
try
    try
        LFPfromDat(basepath, 'outFs', 1250, 'useGPU', false);
    catch e
        fprintf(1, 'The identifier was:\n%s', e.identifier);
        fprintf(1, 'There was an error! The message was:\n%s', e.message);
        if (exist([basepath, '\', basename, '.lfp'], "file") ~= 0)
            fclose([basepath, '\', basename, '.lfp']); %if the above run failed after starting the file
            delete([basepath, '\', basename, '.lfp']);
        end
        LFPfromDat(basepath, 'outFs', 1250, 'useGPU', false);
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

%% Clean data  - CHECK FOR OUR LAB
% Remove stimulation artifacts
if cleanArtifacts && analogInputs
    [pulses] = getAnalogPulses(analogInp, 'analogCh', analogChannels);
    cleanPulses(pulses.ints{1}(:));
end

%% Get brain states
% an automatic way of flaging bad channels is needed
if stateScore
    try
        if exist('pulses', 'var')
            SleepScoreMaster(basepath, 'noPrompts', true, 'ignoretime', pulses.intsPeriods, 'rejectChannels', session.channelTags.Bad.channels); % try to sleep score
        else
            SleepScoreMaster(basepath, 'noPrompts', true, 'rejectChannels', session.channelTags.Bad.channels); % takes lfp in base 0
        end
    catch e
        warning('Problem with SleepScore scoring... unable to calculate');
        fprintf(1, 'The identifier was:\n%s', e.identifier);
        fprintf(1, 'There was an error! The message was:\n%s', e.message);
    end
end


% remove noise from data for cleaner spike sorting
if removeNoise
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

%% Kilosort concatenated sessions - Needs to be changed to probes, not shanks HLR 01/05/2023
if spikeSort
    if nKilosortRuns > 1 % if more than one Kilosort cycle desired, break the shanks down into the desired number of kilosort runs
        shanks = session.extracellular.spikeGroups.channels;
        kilosortGroup = ceil(((1:length(shanks)) / nKilosortRuns));
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
        %     PhyAutoClustering(kilosortFolder);
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
