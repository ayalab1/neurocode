%% Pare events
function HSE = pareEvt(HSE, spikes, numCell, numSpk, singleHz)
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
tempKeep = find(tempInd >= numCell);
tempKeep2 = find(tempSing==1);
HSE.keep = sort(unique(cat(1,tempKeep,tempKeep2)));
save([basepath '\Barrage_Files\' basename '.HSE.mat'], 'HSE');
end