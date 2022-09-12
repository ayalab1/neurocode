function  preprocessSession(varargin)

%         preprocessSession(varargin)

%   Master function to run the basic pre-processing pipeline for an
%   individual sessions. Is based on sessionsPipeline.m but in this case
%   works on an individual session basis no in a folfer with multiple ones.
%

% INPUTS
%   <options>       optional list of property-value pairs (see table below)
%   basepath        - Basepath for experiment. It contains all session
%                       folders. If not provided takes pwd.
%   analogCh       - List of analog channels with pulses to be detected (it support Intan Buzsaki Edition).
%   forceSum       - Force make folder summary (overwrite, if necessary). Default false.
%   cleanArtifacts - Remove artifacts from dat file. By default, if there is analogEv in folder, is true.
%   stateScore     - Run automatic brain state detection with SleepScoreMaster. Default true.
%   spikeSort      - Run automatic spike sorting using Kilosort. Default true.
%   cleanRez       - Run automatic noise detection of the Kilosort results (these will be pre-labelled as noise in phy). Default true.
%   getPos         - get tracking positions. Default false.
%   runSummary     - run summary analysis using AnalysisBatchScrip. Default false.
%   pullData       - Path for raw data. Look for not analized session to copy to the main folder basepath. To do...
%   path_to_dlc_bat_file - path to your dlc bat file to analyze your videos (see neurocode\behavior\dlc for example files)
%   multiKilosort  - Nx1 cell of channels to exclude in each of N separate
%                    KiloSort runs. Within each cell, there should be a 
%                    numeric array of channels to exclude. Default {}. 
%                    Useful in the event of long recordings with multiple 
%                    electrodes or many recording sites. 
%
%  HISTORY:
%   AntonioFR, 5/20

%  TO DO:
%   - Improve auto-clustering routine
%   - Test cleaning and removing artifacts routines


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Set up parameters and parse inputs

p = inputParser;
addParameter(p,'basepath',pwd,@isfolder); % by default, current folder
addParameter(p,'fillMissingDatFiles',false,@islogical);
addParameter(p,'fillTypes',[],@iscellstr);
addParameter(p,'analogInputs',false,@islogical);
addParameter(p,'analogChannels',[],@isnumeric);
addParameter(p,'digitalInputs',false,@islogical);
addParameter(p,'digitalChannels',[],@isnumeric);
addParameter(p,'getAcceleration',false,@islogical);
addParameter(p,'cleanArtifacts',false,@islogical);
addParameter(p,'stateScore',true,@islogical);
addParameter(p,'spikeSort',true,@islogical);
addParameter(p,'cleanRez',true,@islogical);
addParameter(p,'getPos',false,@islogical);
addParameter(p,'removeNoise',false,@islogical); % raly: noise removal is bad, it removes periods 20ms after (because of the filter shifting) a peak in high gamma. See ayadata1\home\raly\Documents\notes\script_NoiseRemoval_bad.m for details.
addParameter(p,'runSummary',false,@islogical);
addParameter(p,'SSD_path','D:\KiloSort',@ischar)    % Path to SSD disk. Make it empty to disable SSD
addParameter(p,'path_to_dlc_bat_file','',@isfile) 
addParameter(p,'multiKilosort',{},@iscell);

% addParameter(p,'pullData',[],@isdir); To do...
parse(p,varargin{:});

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
multiKilosort = p.Results.multiKilosort;

if ~exist(basepath,'dir')
    error('path provided does not exist')
end
cd(basepath)

%% Pull meta data

% Get session names
if strcmp(basepath(end),filesep)
    basepath = basepath(1:end-1);
end
[~,basename] = fileparts(basepath);

% Get xml file in order
xmlFile = checkFile('fileType','.xml','searchSubdirs',true);
xmlFile = xmlFile(1);
if ~(strcmp(xmlFile.folder,basepath)&&strcmp(xmlFile.name(1:end-4),basename))
    copyfile([xmlFile.folder,filesep,xmlFile.name],[basepath,filesep,basename,'.xml'])
end

% Check info.rhd
% (assumes this will be the same across subsessions)
rhdFile = checkFile('fileType','.rhd','searchSubdirs',true);
rhdFile = rhdFile(1);
if ~(strcmp(rhdFile.folder,basepath)&&strcmp(rhdFile.name(1:end-4),basename))
    copyfile([rhdFile.folder,filesep,rhdFile.name],[basepath,filesep,basename,'.rhd'])
end

%% Make SessionInfo
% Manually ID bad channels at this point. automating it would be good
session = sessionTemplate(basepath,'showGUI',false);
save(fullfile(basepath,[basename, '.session.mat']),'session');

%% Fill missing dat files of zeros
if fillMissingDatFiles
    if isempty(fillTypes)
        fillTypes = {'analogin';'digitalin';'auxiliary';'time';'supply'};
    end
    for ii = 1:length(fillTypes)
        fillMissingDats('basepath',basepath,'fileType',fillTypes{ii});
    end
end
%% Concatenate sessions
disp('Concatenate session folders...');
concatenateDats(basepath,1);

%% run again to add epochs from basename.MergePoints.m
session = sessionTemplate(basepath,'showGUI',false);
save(fullfile(basepath,[basename, '.session.mat']),'session');

%% Process additional inputs - CHECK FOR OUR LAB

% Analog inputs
% check the two different fucntions for delaing with analog inputs and proably rename them
if analogInputs
    if  ~isempty(analogChannels)
        analogInp = computeAnalogInputs('analogCh',analogChannels,'saveMat',true,'fs',session.extracellular.sr);
    else
        analogInp = computeAnalogInputs('analogCh',[],'saveMat',true,'fs',session.extracellular.sr);
    end
    
    % analog pulses ...
    [pulses] = getAnalogPulses('samplingRate',session.extracellular.sr);
end

% Digital inputs
if digitalInputs
    if ~isempty(digitalChannels)
        % need to change to only include specified channels
        digitalInp = getDigitalIn('all','fs',session.extracellular.sr);
    else
        digitalInp = getDigitalIn('all','fs',session.extracellular.sr);
    end
end

% Auxilary input
if getAcceleration
    accel = computeIntanAccel('saveMat',true);
end

%% Make LFP
try
    LFPfromDat(basepath,'outFs',1250,'useGPU',false);
catch
    try
        warning('LFPfromDat failed, trying ResampleBinary')
        ResampleBinary([basepath '\' basename '.dat'],[basepath '\' basename '.lfp'],session.extracellular.nChannels,1,16);
    catch
        warning('LFP file could not be generated, moving on');
    end
end

% 'useGPU'=true gives an error if CellExplorer in the path. Need to test if
% it is possible to remove the copy of iosr toolbox from CellExplorer

%% Clean data  - CHECK FOR OUR LAB
% Remove stimulation artifacts
if cleanArtifacts && analogInputs
    [pulses] = getAnalogPulses(analogInp,'analogCh',analogChannels);
    cleanPulses(pulses.ints{1}(:));
end

% remove noise from data for cleaner spike sorting
if removeNoise
    NoiseRemoval(basepath); % not very well tested yet
end

%% Get brain states
% an automatic way of flaging bad channels is needed
if stateScore
    try
        if exist('pulses','var')
            SleepScoreMaster(basepath,'noPrompts',true,'rejectChannels',session.channelTags.Bad.channels,'ignoretime',pulses.intsPeriods); % try to sleep score
            thetaEpochs(basepath);
        else
            SleepScoreMaster(basepath,'noPrompts',true,'rejectChannels',session.channelTags.Bad.channels); % takes lfp in base 0
            thetaEpochs(basepath);
        end
    catch
        warning('Problem with SleepScore scoring... unable to calculate');
    end
end

%% Kilosort concatenated sessions
if spikeSort
    if length(multiKilosort)==0
        kilosortFolder = KiloSortWrapper('SSD_path',SSD_path,'rejectchannels',session.channelTags.Bad.channels);
        %     PhyAutoClustering(kilosortFolder);
        if cleanRez
            load(fullfile(kilosortFolder,'rez.mat'),'rez');
            CleanRez(rez,'savepath',kilosortFolder);
        end
    else
        for i = 1:length(multiKilosort)
            excludeKilo = [];
            excludeKilo = session.channelTags.Bad.channels;
            excludeKilo = [excludeKilo multiKilosort{i}];
            kilosortFolder = KiloSortWrapper('SSD_path',SSD_path,'rejectchannels',excludeKilo);
            if cleanRez
                load(fullfile(kilosortFolder,'rez.mat'),'rez');
                CleanRez(rez,'savepath',kilosortFolder);
            end
        end
    end
end

%% run dlc
if exist(path_to_dlc_bat_file,'file')
    system(path_to_dlc_bat_file,'-echo')
end

%% Get tracking positions - TO FIX
if getPos
    % put tracking into standard format
    general_behavior_file('basepath',basepath)
end

%% Summary -  NOT WELL IMPLEMENTED YET
if runSummary
    if forceSum || (~isempty(dir('*Kilosort*')) && (isempty(dir('summ*')))) % is kilosorted but no summ
        disp('running summary analysis...');
        sessionSummary;
    end
end
end