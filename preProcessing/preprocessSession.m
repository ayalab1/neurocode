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
%   getPos         - get tracking positions. Default true. 
%   runSummary     - run summary analysis using AnalysisBatchScrip. Default false.
%   pullData       - Path for raw data. Look for not analized session to copy to the main folder basepath. To do...
%
%  HISTORY: 
%   AntonioFR, 5/20

%  TO DO:
%   - Improve auto-clustering routine 

% write file to keep track of 

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
addParameter(p,'stateScore',false,@islogical);
addParameter(p,'spikeSort',true,@islogical);
addParameter(p,'getPos',true,@islogical);
addParameter(p,'removeNoise',true,@islogical);
addParameter(p,'runSummary',false,@islogical);

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
getPos = p.Results.getPos;
removeNoise = p.Results.removeNoise;
runSummary = p.Results.runSummary;

if ~exist('basepath')
    error('path provided does not exist')
end
cd(basepath);

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

%% Make SessionInfo
% ID bad channels at this point. automating it would be good

session = sessionTemplate(pwd,'showGUI',true); %
save([basename '.session.mat'],'session');

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

cd(basepath);

disp('Concatenate session folders...');
concatenateDats();
% NEED TO FIX MergePoints.foldernames or revert to original ConcatenateDats

%% Process additional inputs

% Analog input
if analogInputs
    if  ~isempty(analogChannels)
        analogInp = computeAnalogInputs('analogCh',analogChannels,'saveMat',true);
    else
        analogInp = computeAnalogInputs('analogCh',[],'saveMat',true); 
    end
end
% analog pulses ... 

% Digital input
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

LFPfromDat(pwd,'outFs',1250,'useGPU',true);
% 'useGPU'=true gives an error if CellExplorer in the path. Need to test if
% it is possible to remove the copy of iosr toolbox from CellExplorer

%% Clean data
% Remove stimulation artifacts
if cleanArtifacts && analogInputs
    [pulses] = getAnalogPulses(analogInp,'analogCh',analogChannels);
    cleanPulses(pulses.ints{1}(:));
end

% remove noise from data for cleaner spike sorting
if removeNoise
    NoiseRemoval(pwd); % not very well tested yet
end

%% Get brain states
% an automatic way of flaging bad channels is needed 
if stateScore 
    try 
        if exist('pulses','var')
            SleepScoreMaster(pwd,'noPrompts',true,'ignoretime',pulses.intsPeriods); % try to sleep score
            thetaEpochs(pwd); 
        else
            SleepScoreMaster(pwd,'noPrompts',true); % takes lfp in base 0
            thetaEpochs(pwd); 
        end
    catch
        disp('Problem with SleepScore skyping...');
    end
end

%% Kilosort concatenated sessions
if spikeSort
    kilosortFolder = KiloSortWrapper('SSD_path','E:');
    PhyAutoClustering(strcat(basepath,filesep,kilosortFolder));
end
%% Get tracking positions 
if getPos
    getSessionTracking('optitrack',true);
end

%% Summary -  NOT WELL IMPLEMENTED YET
if runSummary
    if forceSum || (~isempty(dir('*Kilosort*')) && (isempty(dir('summ*')))) % is kilosorted but no summ
        disp('running summary analysis...');
        sessionSummary;         
    end
end
end

