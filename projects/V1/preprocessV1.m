% function  preprocessV1(varargin)
% 
% %         preprocessSession(varargin)
% 
% %   Master function to run the basic pre-processing pipeline for an
% %   individual sessions. Is based on sessionsPipeline.m but in this case
% %   works on an individual session basis no in a folfer with multiple ones.
% % 
% 
% % INPUTS
% %   <options>       optional list of property-value pairs (see table below)
% %   basepath        - Basepath for experiment. It contains all session
% %                       folders. If not provided takes pwd.
% %   analogCh       - List of analog channels with pulses to be detected (it support Intan Buzsaki Edition).
% %   forceSum       - Force make folder summary (overwrite, if necessary). Default false.
% %   cleanArtifacts - Remove artifacts from dat file. By default, if there is analogEv in folder, is true.
% %   stateScore     - Run automatic brain state detection with SleepScoreMaster. Default true.
% %   spikeSort      - Run automatic spike sorting using Kilosort. Default true.
% %   getPos         - get tracking positions. Default true. 
% %   runSummary     - run summary analysis using AnalysisBatchScrip. Default false.
% %   pullData       - Path for raw data. Look for not analized session to copy to the main folder basepath. To do...
% %
% %  HISTORY: 
% %   AntonioFR, 5/20
% 
% %  TO DO:
% %   - Improve auto-clustering routine 
% %   - Test cleaning and removing artifacts routines
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %% Set up parameters and parse inputs
% 
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
addParameter(p,'getPos',false,@islogical);
addParameter(p,'removeNoise',true,@islogical);
addParameter(p,'runSummary',false,@islogical);
addParameter(p,'SSD_path','D:\KiloSort',@ischar)    % Path to SSD disk. Make it empty to disable SSD

% addParameter(p,'pullData',[],@isdir); To do... 
% parse(p,varargin{:});
parse(p);
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
SSD_path = p.Results.SSD_path;

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
rhdFile = checkFile('fileType','.rhd','searchSubdirs',true);
rhdFile = rhdFile(1);
if ~(strcmp(rhdFile.folder,basepath)&&strcmp(rhdFile.name(1:end-4),basename))
    copyfile([rhdFile.folder,filesep,rhdFile.name],[basepath,filesep,basename,'.rhd'])
end

%% Make SessionInfo
% Manually ID bad channels at this point. automating it would be good

session = sessionTemplate(pwd,'showGUI',false); % show GUI only after concatenating data and getting MergePoints
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
concatenateDats(pwd,0,0);

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
try
    digitalInp = getDigitalIn('all','fs',session.extracellular.sr);
catch
    display('Error with digital inputs. This step was skipped.');
end
% Auxilary input
if getAcceleration
    accel = computeIntanAccel('saveMat',true); 
end

%% Make LFP

LFPfromDat(pwd,'outFs',1250,'useGPU',false);

% 'useGPU'=true ; gives an error if CellExplorer in the path. Need to test if
% it is possible to remove the copy of iosr toolbox from CellExplorer

%% Clean data  - CHECK FOR OUR LAB
% Remove stimulation artifacts
if cleanArtifacts && analogInputs
    [pulses] = getAnalogPulses(analogInp,'analogCh',analogChannels);
    cleanPulses(pulses.ints{1}(:));
end

% remove noise from data for cleaner spike sorting
% if removeNoise
%     NoiseRemoval(pwd); % not very well tested yet: Raly: this is bad
% end

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
    kilosortFolder = KiloSortWrapper('SSD_path',SSD_path);
    load(fullfile(kilosortFolder,'rez.mat'),'rez');
    CleanRez(rez,'savepath',kilosortFolder);
    %     PhyAutoClustering(kilosortFolder);
end
% move datfile to kilosort folder
% movefile(fullfile(basepath,[basename '.dat']),fullfile(kilosortFolder,[basename '.dat']));

input('Indicate when you are done with phy spike sorting')
% movefile(fullfile(kilosortFolder,[basename '.dat']),fullfile(basepath,[basename '.dat']));
cd(basepath);
session = sessionTemplate(pwd,'showGUI',true,'remove_folder_date',true);
f = dir('Kilosort*');

spikes = loadSpikes('session',session,'clusteringpath',[f.folder filesep f.name]);
cell_metrics = ProcessCellMetrics('session',session,'manualAdjustMonoSyn',false,'excludeMetrics',{'deepSuperficial'});
channel_mapping
close all
% copyfile(fullfile(kilosortFolder,[basename '.session.mat']),fullfile(basepath,[basename '.session.mat']));
% copyfile(fullfile(kilosortFolder,[basename '.cell_metrics.cellinfo.mat']),fullfile(basepath,[basename '.cell_metrics.cellinfo.mat']));
% copyfile(fullfile(kilosortFolder,[basename '.spikes.cellinfo.mat']),fullfile(basepath,[basename '.spikes.cellinfo.mat']));
% copyfile(fullfile(kilosortFolder,[basename '.mono_res.cellinfo.mat']),fullfile(basepath,[basename '.mono_res.cellinfo.mat']));
% copyfile(fullfile(kilosortFolder,[basename '.noiseLevel.channelInfo.mat']),fullfile(basepath,[basename '.noiseLevel.channelInfo.mat']));
% copyfile(fullfile(kilosortFolder,[basename '.chanCoords.channelInfo.mat']),fullfile(basepath,[basename '.chanCoords.channelInfo.mat']));

%% Sync stimuli timestamps
SyncStimuliTimestamps


%% Get tracking positions - TO FIX
if getPos
    getSessionTracking('optitrack',false);
end

%% Summary -  NOT WELL IMPLEMENTED YET
if runSummary
    if forceSum || (~isempty(dir('*Kilosort*')) && (isempty(dir('summ*')))) % is kilosorted but no summ
        disp('running summary analysis...');
        sessionSummary;         
    end
end


