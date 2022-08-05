function [curves1,curves2,curves11,curves12,curves21,curves22,basepath] = temp(basepath)

%% basepath = 'M:\Data\Can\OLM21\day10';
if ~exist('basepath','var') || isempty(basepath), basepath = pwd; end
basename = basenameFromBasepath(basepath);
[parentFolder,dayName] = fileparts(basepath);
[~,projectName] = fileparts(parentFolder);
sessionID = [projectName '_' dayName];

load(fullfile(basepath,[basename '.animal.behavior.mat']),'behavior');
load(fullfile(basepath,[basename '.cell_metrics.cellinfo.mat']),'cell_metrics');
% load(fullfile(basepath,[basename '.ripples.events.mat']),'ripples');

stim = behavior.stimON;
run = behavior.run;

%%

% Make sure our each run epoch is confined within the same trial: no run epoch should begin in one trial and continue over the next trial
run = SubtractIntervals(run,bsxfun(@plus,sort(behavior.trials(:)),[-1 1]*0.00001));
run(diff(run,[],2)<0.5,:) = []; % remove run epochs lasting for less than 0.5s

% load spikes
spikesCell = cell_metrics.spikes.times';
nSpikes = cellfun(@(x) sum(x(:,1)>run(1) & x(:,1)<run(end)), spikesCell);
pyr = cellfun(@(x) ~isempty(strfind(x,'Pyramidal')), cell_metrics.putativeCellType)';
spikes = Group(spikesCell{pyr & nSpikes>200});
regions = cell_metrics.brainRegion;
regionNames = unique(regions);
regionCell = zeros(length(regions),1);
for i=1:length(regionNames)
    regionCell(strcmp(regions,regionNames{i})) = i;
end

nonrun = SubtractIntervals([0 Inf],run); nonstim = SubtractIntervals([0 Inf],stim);

trials = behavior.trials;
trialKind = 2-behavior.trialID(:,2);
onTrials = SubtractIntervals(trials(trialKind==1,:),nonrun);
offTrials = SubtractIntervals(trials(trialKind==2,:),nonrun);
pos1 = behavior.positionTrials{1}; pos1 = pos1(~any(isnan(pos1),2),:);
pos2 = behavior.positionTrials{2}; pos2 = pos2(~any(isnan(pos2),2),:);

tt = pos1(:,1); backward = tt(FindInterval(diff(pos1(:,2))<0)); 
stimIntervals = SubtractIntervals(onTrials,ConsolidateIntervals([backward; nonstim; nonrun])); % only stim run periods
tt = pos2(:,1); backward = tt(FindInterval(diff(pos2(:,2))<0)); 
outsideOfRange = tt(FindInterval(~InIntervals(pos2(:,2),quantile(pos1(InIntervals(pos1,stim),2),[0.05 0.95]))));
nonstimIntervals = SubtractIntervals(offTrials,ConsolidateIntervals([stim; backward; outsideOfRange; nonrun])); % make sure nonstimIntervals only take place during offTrials

pos = sortrows([Restrict(pos1,stimIntervals) ; Restrict([pos2(:,1) 1+pos2(:,2)],nonstimIntervals)]); pos(:,2) = pos(:,2)./2;

% compute firing curves
clear curves curves1 curves2 curves11 curves12 curves21 curves22
for i=1:max(spikes(:,2))
    map = FiringCurve(Restrict(pos1,onTrials,'shift','on'),Restrict(spikes(spikes(:,2)==i),onTrials,'shift','on'));
    curves1(i,:) = map.rate;
    map = FiringCurve(Restrict(pos2,offTrials,'shift','on'),Restrict(spikes(spikes(:,2)==i),offTrials,'shift','on'));
    curves2(i,:) = map.rate;

    map = FiringCurve(Restrict(pos1,onTrials(1:2:end,:),'shift','on'),Restrict(spikes(spikes(:,2)==i),onTrials(1:2:end,:),'shift','on'));
    curves11(i,:) = map.rate;
    map = FiringCurve(Restrict(pos1,onTrials(2:2:end,:),'shift','on'),Restrict(spikes(spikes(:,2)==i),onTrials(2:2:end,:),'shift','on'));
    curves12(i,:) = map.rate;
    map = FiringCurve(Restrict(pos2,offTrials(1:2:end,:),'shift','on'),Restrict(spikes(spikes(:,2)==i),offTrials(1:2:end,:),'shift','on'));
    curves21(i,:) = map.rate;
    map = FiringCurve(Restrict(pos2,offTrials(2:2:end,:),'shift','on'),Restrict(spikes(spikes(:,2)==i),offTrials(2:2:end,:),'shift','on'));
    curves22(i,:) = map.rate;
end



% i=i+1;
% 
% j=1;
% 
% clf
% subplot(2,2,3); plot(firingMaps.rateMaps{i}{j});
% try PlotHVLines(peakPosition(i,j),'v','r--'); end
% subplot(2,2,2); plot(firingMaps.occupancy{i}{j});
% try PlotHVLines(peakPosition(i,j),'v','r--'); end
% subplot(2,2,1); plot(firingMaps.countMaps{i}{j});
% try PlotHVLines(peakPosition(i,j),'v','r--'); end
% title(i)
% 
% %%
% load('day19.placeFields.cellinfo.mat')
% load('day19.ratemap.firingRateMap.mat')
% load('day19.animal.behavior.mat')
% load('day19.cell_metrics.cellinfo.mat')
% spikesCell = cell_metrics.spikes.times';
% spikes = Group(spikesCell{:,1});
% run = behavior.run;
% basepath = pwd
% [parentFolder,dayName] = fileparts(basepath);
% [~,animalName] = fileparts(parentFolder);
% switch animalName
% case 'OML18'
% channel = 34;
% case 'OML19'
% channel = 34;
% case 'OLM21'
% channel = 22;
% case 'OML22'
% channel = 163;
% otherwise
% keyboard
% end
% lfp = getLFP(channel+1);
% theta = bz_Filter(lfp,'passband',[5 15]);
% phases = [theta.timestamps, theta.phase];
% phases = Phase([theta.timestamps theta.data],spikes(:,1));
% 
% %%
% xy = [behavior.position.x(:),behavior.position.y(:)];
% trials = behavior.trials; trialKind = 2-behavior.trialID(:,2);
% rotated = RotateCoordinates(xy,0.64,[0 0]); 
% [~,outliersX] = RemoveOutliers(rotated(:,1),10);
% [~,outliersY] = RemoveOutliers(rotated(:,2),10);
% outliers = outliersX | outliersY;
% rotated(outliers,:) = nan;  
% center = [mean([min(rotated(:,1)) max(rotated(:,1))]),...
%     mean([min(rotated(:,2)) max(rotated(:,2))])];
% % center
% rotated = bsxfun(@minus,rotated,center);
% % PlotXY(rotated);
% 
% figure(1);
% ppos = [interp1(behavior.timestamps(:),-rotated(:,1),spikes(:,1)) interp1(behavior.timestamps(:),rotated(:,1),spikes(:,1))];
% % ppos = [interp1(behavior.timestamps(:),behavior.position.linearized(:),spikes(:,1)) interp1(behavior.timestamps(:),behavior.position.linearized(:),spikes(:,1))];
% % ppos = [interp1(behavior.positionTrials{1}(:,1),behavior.positionTrials{1}(:,2),spikes(:,1)) interp1(behavior.positionTrials{2}(:,1),behavior.positionTrials{2}(:,2),spikes(:,1))];
% ppos = [interp1(behavior.positionTrials{1}(:,1),behavior.positionTrials{1}(:,2),spikes(:,1)) interp1(behavior.positionTrials{2}(:,1),behavior.positionTrials{2}(:,2),spikes(:,1))];
% % 
% % i=1;
% % i = i+1; % 23?
% % i= 42; %26; 42
% 
% clf
% for j=1:2
%     subplot(2,1,j)
%     ok = spikes(:,2)==i & InIntervals(spikes(:,1),run) & InIntervals(spikes(:,1),trials(trialKind==j,:));
%     ok = find(ok);
%     [~,w] = InIntervals(spikes(ok,1),trials(trialKind==j,:));
%     RasterPlot([ppos(ok,j) w]);
% %     thatthat = [ppos(ok) w thatthat];
%     plot(repmat(ppos(ok,j),2,1),[phases(ok,2);phases(ok,2)+2*pi],'.')
% %     DensityMap(repmat(ppos(ok,j),2,1),[phases(ok,2);phases(ok,2)+2*pi],'smooth',1,'nBins',[100 100]);
%     ylim([0 2*pi*2]);
%     title(num2str(i));
% %     title([num2str(i) ': ' num2str(peakPosition(i,j))]);
% %     try PlotHVLines(peakPosition(i,j)/100,'v','r--'); end
% end
% 
% %%
% 
% clf
% 
