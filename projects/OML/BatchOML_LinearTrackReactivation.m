function [m,z,post,templateKind,basepath] = BatchOML_LinearTrackReactivation(basepath,conditionType,taskType)

% batch = StartBatch(@BatchOML_LinearTrackReactivation,'OMLproject.batch');

%%
if ~exist('basepath','var') || isempty(basepath), basepath = pwd; end
basename = basenameFromBasepath(basepath);
[parentFolder,dayName] = fileparts(basepath);
[~,projectName] = fileparts(parentFolder);
sessionID = [projectName '_' dayName];

disp([datestr(clock) ': Starting session ' basepath '...']);

load(fullfile(basepath,[basename '.animal.behavior.mat']),'behavior');
load(fullfile(basepath,[basename '.cell_metrics.cellinfo.mat']),'cell_metrics');
% load(fullfile(basepath,[basename '.ripples.events.mat']),'ripples');

stim = behavior.stimON;
run = behavior.run;

load(fullfile(basepath,[basename '.MergePoints.events.mat']),'MergePoints');
sleep = ConsolidateIntervals(MergePoints.timestamps(cellfun(@(x) any(strfind(lower(x),'sleep')),MergePoints.foldernames),:));
try
    load(fullfile(basepath,[basename '.SleepState.states.mat']));
catch
    SleepScoreMaster(basepath,'noPrompts',true,'rejectchannels',[]);
    load(fullfile(basepath,[basename '.SleepState.states.mat']));
end

sws = SleepState.ints.NREMstate;
preSleep = SubtractIntervals(sleep(1,:), SubtractIntervals([0 Inf],sws));
postSleep = SubtractIntervals(sleep(2:end,:), SubtractIntervals([0 Inf],sws));

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
pos1 = behavior.positionTrials{1};
pos2 = behavior.positionTrials{2};

tt = pos1(:,1); backward = tt(FindInterval(diff(pos1(:,2))<0)); 
stimIntervals = SubtractIntervals(onTrials,ConsolidateIntervals([backward; nonstim; nonrun])); % only stim run periods
tt = pos2(:,1); backward = tt(FindInterval(diff(pos2(:,2))<0)); 
outsideOfRange = tt(FindInterval(~InIntervals(pos2(:,2),quantile(pos1(InIntervals(pos1,stim),2),[0.05 0.95]))));
nonstimIntervals = SubtractIntervals(offTrials,ConsolidateIntervals([stim; backward; outsideOfRange; nonrun])); % make sure nonstimIntervals only take place during offTrials

% make sure the order is presleep -> task -> postsleep
preSleep = SubtractIntervals(preSleep,[trials(1) Inf]); postSleep = SubtractIntervals(postSleep,[0 trials(end)]);
bins = Bins(0,sleep(end),0.05,0.01); 
bins = Restrict(bins,[preSleep; postSleep]);
t = mean(bins,2);
pre = InIntervals(t,preSleep);
post = InIntervals(t,postSleep);

[templates,correlations,weights] = ActivityTemplates(Restrict(spikes,stimIntervals,'shift','on'),'binsize',timescale);
% re = ReactivationStrength(spikes,templates);
nTemplates = size(templates,3);
re = ReactivationStrength(spikes,templates,'bins',bins);
[templates,correlations,weights] = ActivityTemplates(Restrict(spikes,nonstimIntervals,'shift','on'),'binsize',timescale);
re2 = ReactivationStrength(spikes,templates,'bins',bins);
re = [re re2(:,2:end)];
nTemplates(2) = size(templates,3);

m = zeros(sum(nTemplates),2);
zre = nanzscore(re);
for k=1:sum(nTemplates)
    m(k,:) = [nanmean(re(pre,k+1)) nanmean(re(post,k+1))];
    z(k,:) = [nanmean(zre(pre,k+1)) nanmean(zre(post,k+1))];
end

templateKind = repelem([1;2],nTemplates);


% rng(0); bursts = FindBursts(Restrict(spikes,sws,'shift','on'),'thresholds',[0 3],'smooth',0.01);
% bursts(:,1:3) = reshape(Unshift(reshape(bursts(:,1:3),[],1),sws),[],3);
% duration = diff(bursts(:,[1 3]),[],2);
% bursts(duration>1,:) = [];
% 
% pre = InIntervals(bursts,preSleep);
% post = InIntervals(bursts,postSleep);
% 
% %%
% bins = Bins(0,sleep(end),0.05,0.001);
% [in,w] = InIntervals(mean(bins,2),bursts(:,[1 3]));
% re = ReactivationStrength(spikes,templates,'bins',bins(in,:));
% m = zeros(size(bursts,1),nTemplates);
% for k=1:nTemplates
%     m(:,k) = Accumulate(w(in),re(:,k+1),'mode','mean','size',[size(bursts,1)]);
% end
% 

















