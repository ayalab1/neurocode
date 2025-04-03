%% Pare events
function HSE = pareBARRs(HSE, spikes, savePath, numCell, numSpk, pareDur, singleHz, stim, maxCell)
if nargin < 8
    stim = 0;
    maxCell = 0;
elseif nargin <9
    maxCell = 0;
end
basepath = pwd;
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


load(strcat(basepath,'\',basename,'.SleepState.states.mat'));
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