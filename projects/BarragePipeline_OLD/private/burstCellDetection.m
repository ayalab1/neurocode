function checking = burstCellDetection(normSpkThresh, loadPath)
%% Set paths
basepath = pwd;
basename = basenameFromBasepath(basepath);

if ~exist(strcat(basepath,'\','Barrage_Files'))
    mkdir('Barrage_Files');
end

savePath = strcat(basepath, '\Barrage_Files\', basename, '.');
load([basename '.cell_metrics.cellinfo.mat']);

%% Produce spike structure
if (nargin <2)||isempty(loadPath)
    spikes = [];
    spikes = importSpikes('cellType',"Pyramidal Cell",'sleepState',"NREMstate");
    % spikes = importSpikes('cellType',"Pyramidal Cell",'brainRegion',"CA2");
    save([savePath 'NREMpyr.cellinfo.mat'], 'spikes');
else
    load(loadPath);
end
%% Flag burst events
burstThresh = 2*(1/20); %maybe look at this in terms of Hz
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
    keep_burst = (numCat>=3); % needs to have at least 5 spikes
    avgBurstSz(i) = mean(numCat(logical(keep_burst)));
    burstSz{i} = numCat(logical(keep_burst));
    clear evtstart evtstop evtdur evtpeak evtamp flagConc keep_burst samples
end

%% Get spike count info

% Maybe track the percentile of number of bursts? So cells in the top
% whatever threshold of burst lengths are kept? Count is vulnerable to
% length of the recording, and fraction of bursts seems to be vulnerable to
% the burstiness properties of the cell itself... Also maybe # spikes that
% are a part of a longer burst?

% useAvgBurstSz = avgBurstSz(~isnan(avgBurstSz));
% thresh = mean(useAvgBurstSz)+(std(useAvgBurstSz));
thresh = 3;
% normSpkThresh = 0.25;
flag = NaN(size(spikes.times,2),1);
checking = NaN(size(spikes.times,2),4);
% useThese = [10,19,23,25,31,47,68,72,77,81,83,85,91,92,93,96,97,101,103,113,125];
for i = 1:size(spikes.times,2)
   resTemp = burstSz{i};
   checking(i,1) = sum(resTemp>thresh);
   checking(i,2) = (sum(resTemp>thresh)/length(resTemp));
   if sum(resTemp>thresh) >= 1
       flag(i) = (sum(resTemp>thresh)/length(resTemp)) > normSpkThresh;
   else
       flag(i) = 0;
   end
%     flag(i) = logical(~isempty(find(useThese==spikes.UID(i),1)));
   checking(i,3) = flag(i);
end

checking(:,4) = spikes.UID';

%% Pull the appropriate UIDs to run through HSE_b
UIDkeep = spikes.UID(logical(flag));
keep.UID = spikes.UID(logical(flag));
keep.times = cell(1,sum(flag));
fc = 1;
for i = 1:length(flag)
    if logical(flag(i))
        keep.times{fc} = spikes.times{i};
        fc = fc+1;
    end
end
spikes = [];
spikes = keep;
save([savePath 'useSpk.cellinfo.mat'], 'spikes');
save([savePath 'useSpk.UIDkeep.mat'],'UIDkeep');
end