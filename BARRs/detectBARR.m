function [HSE] = detectBARR(varargin)
%% Simple BARR detection
%
% This script is a simple detector for CA2 Barrages. The detector works by
% identifying CA2 pyramidal cells which are most likely to barrage (long
% periods of high firing rate), then detecting increases in this
% subpopulation's firing rate. There is an additional cleaning step that
% refines which barrages are kept. 
%
%%%%%%%%%%%%%%%%%%%%
%%% REQUIREMENTS %%%
%%%%%%%%%%%%%%%%%%%%
% Region and cell type must be updated in cell_metrics for more accurate
% detection. CA2 pyramidal cells must be present in order for the detection
% to run. 
%
% Sleep State should also be run and verified ahead of time. 
%
%%%%%%%%%%%%%%
%%% INPUTS %%%
%%%%%%%%%%%%%%
% 
% basepath:     Full path where session is located. Default: pwd
% force:        Logical option to force redetection of BARRs. Default:
%               false
% nSigma:       Number of standard deviations that the subpopulation firing
%               rate must pass for detection. For most sessions this should 
%               be between 2-3. Default: 3
% EMGThresh:    EMG threshold for removing EMG/noisy events. Default: 0.8
% Hz:           Firing rate minimum (in Hz) for BARR unit detection. For 
%               most sessions this should be at least 20-25Hz. Default: 30
%               Hz
% ft:           Time over which cells much fire at "Hz" in order to be
%               flagged as a potential BARR unit, in seconds. For most 
%               sessions, this should be 0.3-0.35. Default: 0.3
% numEvt:       Number of times a potential BARR unit must fire at "Hz" for
%               "ft" seconds to be kept. Because of the restrictive
%               threshold, a small number (10 and under) generally works
%               well. Default: 5
% unMin:        Number of BARR units which must participate in a BARR in
%               order for it to be kept. For most sessions, 2-3 works.
%               Default: 2
% spkNum:       Number of spikes per BARR unit needed in order to be
%               considered a participant in a particular BARR. For most
%               sessions, 4-7 spikes works best. Default: 5
% spkHz:        Firing rate needed for a single BARR unit to be allowed to
%               drive a BARR. Note that this does not mean it's the only 
%               cell firing, just that it's the only flagged cell firing. 
%               Default: 100
% unMax:        Maximum number of flagged units which are allowed to
%               participate in a BARR for it to be kept. This may help
%               prevent misdetection of SWRs. When set to 0, it is ignored.
%               Default: 0
% pareDur:      Minimum duration of BARRs kept, in seconds. This should
%               generally be kept at or below 0.2. Default: 0.2
% stim:         Logical option to indicate whether or not the session
%               includes stimulation (ie ripple generation). If so, 
%               barrages overlapping with SWRs will be removed. Default:
%               false
% remRip:       Logical option to remove BARRs which overlap with ripples.
%               Default: false
%
%%%%%%%%%%%%%%%
%%% OUTPUTS %%%
%%%%%%%%%%%%%%%
% 
% HSE:              BARR structure including...
%   timestamps:     Nx2 timestamps (s) for ALL detections before refinement
%   peaks:          Nx1 peaks (s) for ALL detections before refinement
%   amplitudes:     1xN amplitudes for ALL detections before refinement
%   amplitudeUnits: Unit in which "amplitudes" is stored
%   duration:       1xN duration (s) for ALL detections before refinement
%   center:         1xN center timestamp (s) for ALL detections before
%                   refinement
%   detectorinfo:   Variables used for detection
%   keep:           Mx1 IDXs of refined BARRs, regardless of sleepState. 
%   NREM:           Px1 IDXs of refined BARRs only in NREM sleep
%   nonRip:         Qx1 IDXs of refined BARRs, excluding BARRs overlapping
%                   with SWRs
%   fin:            Qx1 final IDXs related to HSE.timestamps (including
%                   nonRip, which is not always optimal). 
% 
%%% NOTE: 
%   BARR timestamps should generally be indexed as:
%        HSE.timestamps(HSE.keep(HSE.NREM),:)
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

pullSpikes('basepath', basepath, 'savePath', [basepath '\Barrage_Files'], 'force', true); %Get region/cell type spike files

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

HSE = pareBARRs(basepath, HSE, spikes, savePath, unMin, spkNum, pareDur, spkHz, stim, unMax);

%% Save NeuroScope2 file
BARR_N2(basepath);

%% Run analysis script
BARR_PSTH(savePath,"NREM");
end