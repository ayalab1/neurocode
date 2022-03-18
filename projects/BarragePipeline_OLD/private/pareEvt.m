%% Pare events
function HSE = pareEvt(HSE, spikes, numCell, numSpk, singleHz, maxCell)
if nargin < 6
    maxCell = 0;
end
basepath = pwd;
basename = basenameFromBasepath(basepath);
numEvt = size(HSE.timestamps,1);
tempInd = zeros(numEvt,1);
tempPass = zeros(numEvt,length(spikes.UID));
tempSing = zeros(numEvt,1);
for e = 1:numEvt
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
tempKeepMin = find(tempInd >= numCell);
if maxCell == 0
    tempKeep = tempKeepMin;
else
    tempKeepMax = find(tempInd <= maxCell);
    tempKeep = intersect(tempKeepMin, tempKeepMax);
end
tempKeep2 = find(tempSing==1);
HSE.keep = sort(unique(cat(1,tempKeep,tempKeep2)));
load(strcat(basepath,'\',basename,'.SleepState.states.mat'));
if ~isempty(SleepState.ints.NREMstate)
    HSEnREM = eventIntervals(HSE,SleepState.ints.NREMstate,1);
    if ~isempty(HSEnREM)
        [~,HSE.NREM] = intersect(HSE.peaks(HSE.keep), HSEnREM.peaks);
    else
        HSE.NREM = [];
    end
    save([basepath '\Barrage_Files\' basename '.HSE.mat'], 'HSE');
end
end