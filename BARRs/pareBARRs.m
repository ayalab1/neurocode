%% Pare events
function HSE = pareBARRs(basepath, HSE, spikes, savePath, numCell, numSpk, pareDur, singleHz, stim, maxCell)
%% Refine BARR detection
%
% This script pares back BARR detections based on restrictions provided as
% inputs. 
%
%%%%%%%%%%%%%%
%%% INPUTS %%%
%%%%%%%%%%%%%%
% 
%% REQUIRED %%
% basepath:     Location of HSE structure.
% HSE:          BARR structure created by find_HSE_BARR.m
% spikes:       Spikes structure for BARR identified pyramidal cells.
% savePath:     Prefix for HSE file to be saved, ie. [basepath
%               '\Barrage_Files\' basename '.'].
% numCell:      Number of BARR units which must participate in a BARR in
%               order for it to be kept. For most sessions, 2-3 works.
% numSpk:       Number of spikes per BARR unit needed in order to be
%               considered a participant in a particular BARR. For most
%               sessions, 4-7 spikes works best.
% pareDur:      Minimum duration of BARRs kept, in seconds. This should
%               generally be kept at or below 0.2.
% singleHz:     Firing rate needed for a single BARR unit to be allowed to
%               drive a BARR. Note that this does not mean it's the only 
%               cell firing, just that it's the only flagged cell firing.
%% OPTIONAL %%
% stim:         Logical option to indicate whether or not the session
%               includes stimulation (ie ripple generation). If so, 
%               barrages overlapping with SWRs will be removed. Default:
%               false   
% maxCell:      Maximum number of flagged units which are allowed to
%               participate in a BARR for it to be kept. This may help
%               prevent misdetection of SWRs. When set to 0, it is ignored.
%               Default: 0
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
if nargin < 8
    stim = 0;
    maxCell = 0;
elseif nargin <9
    maxCell = 0;
end

if ~contains(savePath, basepath)
    warning('savePath is not within basepath, is this correct?');
    keyboard
end

basename = basenameFromBasepath(basepath);
numEvt = size(HSE.timestamps,1);
tempInd = zeros(numEvt,1);
tempPass = zeros(numEvt,length(spikes.UID));
tempSing = zeros(numEvt,1);
for e = 1:numEvt
    if ((HSE.timestamps(e,2)-HSE.timestamps(e,1))>=pareDur)||(pareDur==0)
        for u = 1:length(spikes.UID)
            tempPass(e,u) = length(find((spikes.times{u}>=HSE.timestamps(e,1))&(spikes.times{u}<=HSE.timestamps(e,2))));
            if tempPass(e,u) >= numSpk
                tempInd(e) = tempInd(e) + 1;
            end
            if ((tempPass(e,u)/HSE.duration(e)) >= singleHz)&&(tempPass(e,u)>3)
                tempSing(e) = 1;
            end
        end
    end
end
tempKeepMin = find(tempInd >= numCell);
if maxCell == 0
    tempKeep = tempKeepMin;
else
    tempKeepMax = find(tempInd <= maxCell);
    tempKeep = intersect(tempKeepMin, tempKeepMax);
end
tempKeep2 = find(tempSing==1);
HSE.keep = sort(unique(cat(1,tempKeep,tempKeep2)));
tempKeep3 = [];
if stim==1
    load([basepath '/' basename '.ripples.events.mat']);
    for j = 1:length(HSE.keep)
        met = [];
        met = (ripples.peaks >= (HSE.timestamps(HSE.keep(j),1)))&(ripples.peaks <= (HSE.timestamps(HSE.keep(j),2)));
        if sum(met)==0
            tempKeep3 = [tempKeep3; HSE.keep(j)];
        end
    end
    HSE.keep = tempKeep3;
    HSE.note = "HSE.keep has removed barrages overlapping with ripples";
end


load(strcat(basepath,filesep,basename,'.SleepState.states.mat'));
if ~isempty(SleepState.ints.NREMstate)
    HSEnREM = eventIntervals(HSE,SleepState.ints.NREMstate,1);
    if ~isempty(HSEnREM)
        [~,HSE.NREM] = intersect(HSE.peaks(HSE.keep), HSEnREM.peaks);
    else
        HSE.NREM = [];
    end
    save([savePath 'HSE.mat'], 'HSE');
end
load([basepath '/' basename '.ripples.events.mat']);
HSE.nonRip = [];
for j = 1:length(HSE.NREM)
    met = [];
    met = (ripples.peaks >= (HSE.timestamps(HSE.keep(HSE.NREM(j)),1)-0.1))&(ripples.peaks <= (HSE.timestamps(HSE.keep(HSE.NREM(j)),2)+0.1));
    if sum(met)==0
        HSE.nonRip = [HSE.nonRip; j];
    end
end
HSE.fin = HSE.keep(HSE.NREM(HSE.nonRip)); 
disp(['Final num: ' num2str(length(HSE.fin))]);
disp(['Num removed: ' num2str(length(HSE.NREM) - length(HSE.fin))]);
save([savePath 'HSE.mat'], 'HSE');
end