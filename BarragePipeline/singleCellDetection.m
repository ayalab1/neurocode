%% Set paths
basepath = pwd;
basename = basenameFromBasepath(basepath);

if ~exist(strcat(basepath,'\','Barrage_Files'))
    mkdir('Barrage_Files');
end

savePath = strcat(basepath, '\Barrage_Files\', basename, '.');
load([basename '.cell_metrics.cellinfo.mat']);

%% Produce spike structure
spikes = [];
spikes = importSpikes('cellType',"Pyramidal Cell",'sleepState',"NREMstate");
save([savePath 'NREMpyr.cellinfo.mat'], 'spikes');
%% Flag burst events
burstThresh = 0.095;
burstEvts = cell(size(spikes.times,2),1);
burstSz = cell(size(spikes.times,2),1);
avgBurstSz = NaN(size(spikes.times,2),1);
for i = 1:size(spikes.times,2)
    tempISI = diff(spikes.times{i});
    evtstart = spikes.times{i};
    evtstart = evtstart(1:length(evtstart)-1);
    evtstop = spikes.times{i};
    evtstop = evtstop(2:end);
    evtdur = evtstop-evtstart;
    evtpeak = evtstart + (evtdur/2);
    evtamp = zeros(length(evtstart),1);
    flagConc = (evtdur <= burstThresh);
    [start,stop,~,~,numCat] = CatCon(evtstart,evtstop,evtpeak,evtamp,flagConc);
    keep_burst = (numCat>=5); % needs to have at least 5 spikes
    avgBurstSz(i) = mean(numCat(logical(keep_burst)));
    burstSz{i} = numCat(logical(keep_burst));
    clear evtstart evtstop evtdur evtpeak evtamp flagConc keep_burst samples
end

%% Get spike count info
thresh = 7;
normSpkThresh = 0.4;
flag = NaN(size(spikes.times,2),1);
checking = NaN(size(spikes.times,2),4);

for i = 1:size(spikes.times,2)
   resTemp = burstSz{i};
   checking(i,1) = sum(resTemp>thresh);
   checking(i,2) = (sum(resTemp>thresh)/length(resTemp));
   flag(i) = (sum(resTemp>thresh)/length(resTemp)) > normSpkThresh;
   checking(i,3) = flag(i);
end

checking(:,4) = spikes.UID';

%% Pull the appropriate UIDs to run through HSE_b
UIDkeep = spikes.UID(logical(flag));
keep.UID = spikes.UID(logical(flag));
for i = 1:length(keep.UID)
    keep.times{i} = spikes.times{spikes.UID(i)};
end
spikes = [];
spikes = keep;
% spikes = importSpikes('UID', UIDkeep);
save([savePath 'brstDt.cellinfo.mat'], 'spikes');
save([savePath 'brstDt.UIDkeep.mat'],'UIDkeep');