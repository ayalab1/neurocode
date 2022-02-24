function [] = unitsForDetection(Hz, ft, numEvt)
% Pull units based on their firing rate characteristics
% INPUTS:
% Hz -> firing rate that each neuron must hit to be considered 
% ft -> amount of time the neuron must fire at the above firing rate to be
%       kept, in seconds. 
%%%
if nargin < 1
    Hz = 20;
    ft = 0.3;
elseif nargin < 2
    ft = 0.3;
end
if ft > 1
    warning('Make sure your firing time is in seconds!');
end
%% Set paths
basepath = pwd;
basename = basenameFromBasepath(basepath);

if ~exist(strcat(basepath,'\','Barrage_Files'))
    mkdir('Barrage_Files');
end

savePath = strcat(basepath, '\Barrage_Files\', basename, '.');
load([basename '.cell_metrics.cellinfo.mat']);

%% Produce spike structure
load([savePath 'CA2pyr.cellinfo.mat']);
% ki = 1;
% for i = 1:length(spikes.UID)
% %     if (~contains("CA1",cell_metrics.brainRegion(spikes.UID(i))))&&(~contains("EC",cell_metrics.brainRegion(spikes.UID(i))))&&(~contains("DG",cell_metrics.brainRegion(spikes.UID(i))))
%     if contains("CA2",cell_metrics.brainRegion(spikes.UID(i)))%||contains("CA3",cell_metrics.brainRegion(spikes.UID(i)))
%         keep.UID(ki) = spikes.UID(i);
%         keep.times{ki} = spikes.times{i};
%         ki = ki+1;
%     end
% end
% spikes = []; spikes = keep;
% spikes = importSpikes('cellType',"Pyramidal Cell",'brainRegion',"CA2");
% save([savePath 'NREMpyr.cellinfo.mat'], 'spikes');
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
    flagConc{u} = (unFR{u} >= Hz); %CHANGED FROM 20
    [start{u},stop{u},~,~,numCat{u}] = CatCon(evtstart{u},evtstop{u},evtpeak{u},evtamp{u},flagConc{u});
    samples{u} = Restrict(spikes.times{u}, [start{u}' stop{u}']);
end
clear evtstart evtstop evtdur evtpeak evtamp flagConc
% thresh = mean(flagging(:,2))+(1.2*std(flagging(:,2)));
% flagging(:,3)=flagging(:,2)>=thresh;

%% Pull the appropriate UIDs to run through HSE_b
for u = 1:length(spikes.UID)
    flagging(u,2) = length(find(unFR{u} >= Hz)); %CHANGED FROM 20
    flagging(u,3) = length(find(numCat{u} >= (ft/binsz))); %300 ms, find #times cell fires at 15Hz for at least 300ms
    if ~isnan(flagging(u,3)) %if any events (I was doing mean(flagging(:,3)>10) for some reason??? Why????)
        flagging(u,4) = flagging(u,3)>=numEvt; %..this is the thresh we should play with maybe?
    else
        flagging(u,4) = 0; %was flagging(u,4) = flagging(u,3)>0; (WHY)
    end
%     br(u,1) = cell_metrics.brainRegion(spikes.UID(u));
end
clear samples
flag = logical(flagging(:,4));
% flag = ismember(spikes.UID, [90 70 35 6 18 81 10 2 121 100]);
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
save([savePath 'useSpk.cellinfo.mat'], 'spikes');
save([savePath 'useSpk.UIDkeep.mat'],'UIDkeep');
end