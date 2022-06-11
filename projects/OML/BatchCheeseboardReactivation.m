function [m,z,basepath,conditionType,taskType] = BatchCheeseboardReactivation(basepath,conditionType,taskType)

% batch = StartBatch(@BatchCheeseboardReactivation,'OMLcheese.batch');

%%
timescale = 0.02;
basename = basenameFromBasepath(basepath);
[parentFolder,dayName] = fileparts(basepath);
[~,projectName] = fileparts(parentFolder);
sessionID = [projectName '_' dayName];

disp([datestr(clock) ': starting session ' basepath '...']);

cd(basepath);
MergePoints = getStruct(basepath,'MergePoints');
% ripples = getStruct(basepath,'ripples')
sleep = ConsolidateIntervals(MergePoints.timestamps(cellfun(@(x) any(strfind(lower(x),'sleep')),MergePoints.foldernames),:));
load(fullfile(basepath,'cheesbTrialsDayOffset.mat'),'cheesbTrialsDay');
trials = [cheesbTrialsDay.trials.start cheesbTrialsDay.trials.end];
postprobe = [cheesbTrialsDay.postProbe.start cheesbTrialsDay.postProbe.end];
postprobe = trials;

cell_metrics = getStruct(basepath,'cell_metrics');
pyr = cellfun(@(x) ~isempty(strfind(x,'Pyramidal')), cell_metrics.putativeCellType)';
spikesCell = cell_metrics.spikes.times';
nSpikes = cellfun(@(x) sum(InIntervals(x,trials)), spikesCell);
spikes = Group(spikesCell{pyr & nSpikes>200});

try
    load(fullfile(basepath,[basename '.SleepState.states.mat']));
catch
    SleepScoreMaster(basepath,'noPrompts',true,'rejectchannels',[]);
    load(fullfile(basepath,[basename '.SleepState.states.mat']));
end

sws = SleepState.ints.NREMstate;
preSleep = SubtractIntervals(sleep(1,:), SubtractIntervals([0 Inf],sws));
postSleep = SubtractIntervals(sleep(2:end,:), SubtractIntervals([0 Inf],sws));

% make sure the order is presleep -> task -> postsleep
preSleep = SubtractIntervals(preSleep,[trials(1) Inf]); postSleep = SubtractIntervals(postSleep,[0 trials(end)]);
[templates,correlations,weights] = ActivityTemplates(Restrict(spikes,postprobe,'shift','on'),'binsize',timescale);
% re = ReactivationStrength(spikes,templates);
nTemplates = size(templates,3);
bins = Bins(0,sleep(end),timescale,min([0.01 timescale])); 
% bins = Restrict(bins,[preSleep; trials; postSleep]);
t = mean(bins,2);
pre = InIntervals(t,preSleep);
post = InIntervals(t,postSleep);
task = InIntervals(t,trials);
re = ReactivationStrength(spikes,templates,'bins',bins);
m = zeros(nTemplates,3); z = m;
zre = nanzscore(re);
for k=1:nTemplates
    m(k,:) = [nanmean(re(pre,k+1)) nanmean(re(task,k+1))  nanmean(re(post,k+1))];
    z(k,:) = [nanmean(zre(pre,k+1)) nanmean(zre(task,k+1))  nanmean(re(post,k+1))];
end



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

















