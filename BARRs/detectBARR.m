%% Simple BARR detection
function [HSE] = detectBARR(varargin)

p = inputParser;
addParameter(p, 'basepath', pwd, @isfolder);
addParameter(p, 'force', false, @islogical);
addParameter(p, 'nSigma', 3, @isnumeric);
addParameter(p, 'EMGThresh', 0.8, @isnumeric);
addParameter(p, 'Hz', 30, @isnumeric);
addParameter(p, 'ft', 0.3, @isnumeric);
addParameter(p, 'numEvt', 5, @isnumeric);
addParameter(p, 'unMin', 2, @isnumeric);
addParameter(p, 'spkNum', 5, @isnumeric);
addParameter(p, 'spkHz', 100, @isnumeric);
addParameter(p, 'unMax', 0, @isnumeric);
addParameter(p, 'pareDur', 0.2, @isnumeric);
addParameter(p, 'stim', 0, @isnumeric);
addParameter(p, 'remRip', false, @islogical);

parse(p, varargin{:});

basepath = p.Results.basepath;
force = p.Results.force;
nSigma = p.Results.nSigma;
EMGThresh = p.Results.EMGThresh;
Hz = p.Results.Hz;
ft = p.Results.ft;
numEvt = p.Results.numEvt;
unMin = p.Results.unMin;
spkNum = p.Results.spkNum;
spkHz = p.Results.spkHz;
unMax = p.Results.unMax;
pareDur = p.Results.pareDur;
stim = p.Results.stim;
remRip = p.Results.remRip;

%% Check if BARRs have already been detected
basename = basenameFromBasepath(basepath);
animal = animalFromBasepath(basepath);

if exist([basepath '\Barrage_Files\' basename '.HSE.mat'])&&~force
    disp('BARRs already detected! Loading file...');
    load([basepath '\Barrage_Files\' basename '.HSE.mat']);
    return
end

%% Prep for detection
savePath = convertStringsToChars(strcat(basepath, '\Barrage_Files\', basename, '.'));
if ~exist([basepath '\Barrage_Files'])
    mkdir('Barrage_Files');
end

pullSpikes([basepath '\Barrage_Files']); %Get region/cell type spike files

if ~exist(strcat(basepath,'\Barrage_Files\',basename,'.CA2pyr.cellinfo.mat'))
    disp('No CA2 pyramidal cells detected, exiting');
    return
end

if ~exist(strcat(basepath,'\',basename,'.SleepState.states.mat'))
    disp('SleepScore not detected, calculating');
    load([basepath '\' basename '.session.mat']);
    SleepState = SleepScoreMaster(basepath,'rejectChannels',session.channelTags.Bad.channels);
    clear session
else
    load(strcat(basepath,'\',basename,'.SleepState.states.mat'));
end

%% Detect good candidate units
unitsForDetection(Hz, ft, numEvt, savePath);
loadPath = strcat(savePath,'useSpk.cellinfo.mat');
note_all = "DetectUn";

load([savePath 'useSpk.UIDkeep.mat']);
load([savePath 'useSpk.cellinfo.mat']); %load in the spikes that we've picked and run to save

if isempty(spikes.UID)
    load([savePath 'CA2pyr.cellinfo.mat']);
    UIDkeep = spikes.UID;
    save([savePath 'useSpk.UIDkeep.mat'],'UIDkeep');
    warning('unitsForDetection did not return any units, defaulting to all CA2pyr');
    note_all = "Failed DetectUn";
end

%% Run BARR detection
HSE = find_HSE_BARR('spikes',spikes,'nSigma',nSigma,'binSz',0.005,'tSmooth',0.02,...
                'tSepMax',0.005,'mindur',0.05,'maxdur',10,'lastmin',0.05,'EMGThresh',EMGThresh,...
                'Notes',note_all,'sstd',-1*(nSigma-0.5),'estd',(nSigma-0.5),...
                'recordMetrics',true,'remRip',remRip);

HSE = pareBARRs(HSE, spikes, savePath, unMin, spkNum, pareDur, spkHz, stim, unMax);

%% Save NeuroScope2 file
BARR_N2('basepath',basepath);

%% Run analysis script
BARR_PSTH(savePath,"NREM");
end