function [] = unitsForDetection()
% Pull units based on their firing rate characteristics
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
ki = 1;
for i = 1:length(spikes.UID)
%     if (~contains("CA1",cell_metrics.brainRegion(spikes.UID(i))))&&(~contains("EC",cell_metrics.brainRegion(spikes.UID(i))))&&(~contains("DG",cell_metrics.brainRegion(spikes.UID(i))))
    if contains("CA2",cell_metrics.brainRegion(spikes.UID(i)))||contains("CA3",cell_metrics.brainRegion(spikes.UID(i)))
        keep.UID(ki) = spikes.UID(i);
        keep.times{ki} = spikes.times{i};
        ki = ki+1;
    end
end
spikes = []; spikes = keep;
% spikes = importSpikes('cellType',"Pyramidal Cell",'brainRegion',"CA2");
save([savePath 'NREMpyr.cellinfo.mat'], 'spikes');
%% Get cell by cell firing rate
tSmooth = 0.05; binsz = 0.01;
flagging = nan(length(spikes.UID),4);
flagging(:,1) = spikes.UID';
parfor u = 1:length(spikes.UID)
    [unFR{u},~,~,ts{u}] = spkRtHist(spikes.times{u}, 'tSmooth', tSmooth, 'binsz', binsz, 'ifz', false);
    unFR{u} = unFR{u}/binsz; %can add the condition that the length of these high firing events must be
    % at least t bins?
    evtstart{u} = ts{u};
    evtstop{u} = ts{u};
    evtstop{u} = [evtstop{u}(2:end) (evtstop{u}(end)+(evtstart{u}(2)-evtstart{u}(1)))];
    evtdur{u} = evtstop{u}-evtstart{u};
    evtpeak{u} = evtstart{u} + (evtdur{u}/2);
    evtamp{u} = zeros(length(evtstart{u}),1);
    flagConc{u} = (unFR{u} >= 20);
    [start{u},stop{u},~,~,numCat{u}] = CatCon(evtstart{u},evtstop{u},evtpeak{u},evtamp{u},flagConc{u});
    samples{u} = Restrict(spikes.times{u}, [start{u}' stop{u}']);
end
clear evtstart evtstop evtdur evtpeak evtamp flagConc
% thresh = mean(flagging(:,2))+(1.2*std(flagging(:,2)));
% flagging(:,3)=flagging(:,2)>=thresh;

%% Pull the appropriate UIDs to run through HSE_b
for u = 1:length(spikes.UID)
    flagging(u,2) = length(find(unFR{u} >= 20)); 
    flagging(u,3) = length(find(numCat{u} >= (0.3/binsz))); %300 ms
    if mean(flagging(:,3)>10)
        flagging(u,4) = flagging(u,3)>=5; %this is the thresh we should play with maybe?
    else
        flagging(u,4) = flagging(u,3)>0;
    end
%     br(u,1) = cell_metrics.brainRegion(spikes.UID(u));
end
clear samples
flag = logical(flagging(:,4));
UIDkeep = spikes.UID(flag);
% keptBR = br(flag);
keep.UID = spikes.UID(flag);
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
save([savePath 'brstDt.cellinfo.mat'], 'spikes');
save([savePath 'brstDt.UIDkeep.mat'],'UIDkeep');
end