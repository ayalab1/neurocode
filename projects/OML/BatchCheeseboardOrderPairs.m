function [p,ab,p0,ab0,basepath,timescale,conditionType,taskType] = BatchCheeseboardOrderPairs(basepath,conditionType,taskType)

% batch = StartBatch(@BatchCheeseboardReactivation,'OMLcheese.batch');

%%
rng(0);

if ~exist('basepath','var') || isempty(basepath), basepath = pwd; end

basename = basenameFromBasepath(basepath);
[parentFolder,dayName] = fileparts(basepath);
[~,projectName] = fileparts(parentFolder);
sessionID = [projectName '_' dayName];
timescale = 0.02;

disp([datestr(clock) ': starting session ' basepath '...']);

cd(basepath);
MergePoints = getStruct(basepath,'MergePoints');
% ripples = getStruct(basepath,'ripples')
sleep = ConsolidateIntervals(MergePoints.timestamps(cellfun(@(x) any(strfind(lower(x),'sleep')),MergePoints.foldernames),:));
load(fullfile(basepath,'cheesbTrialsDayOffset.mat'),'cheesbTrialsDay');
trials = [cheesbTrialsDay.trials.start cheesbTrialsDay.trials.end];
postprobe = [cheesbTrialsDay.postProbe.start cheesbTrialsDay.postProbe.end];

cell_metrics = getStruct(basepath,'cell_metrics');
pyr = cellfun(@(x) ~isempty(strfind(x,'Pyramidal')), cell_metrics.putativeCellType)';
spikesCell = cell_metrics.spikes.times';
nSpikes = cellfun(@(x) sum(InIntervals(x,trials)), spikesCell);
spikes = Group(spikesCell{pyr & nSpikes>200});

% postprobe = [cheesbTrialsDay.postProbe.start cheesbTrialsDay.postProbe.end];

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
spikes = Restrict(spikes,ConsolidateIntervals(sortrows([preSleep; trials; postSleep; postprobe])));

nUnits = max(spikes(:,2));
indices = combnk(1:nUnits,2);
pairs = false(length(indices),nUnits); % lArray stands for "logical array", where each line is an assembly and each column is a neuron,
% so [1 1 0 1 0 0;...] would mean that neurons 1, 2 and 4 together form the first assembly
pairs(sub2ind(size(pairs),(1:size(pairs,1))',indices(:,1)))=true;
pairs(sub2ind(size(pairs),(1:size(pairs,1))',indices(:,2)))=true;

activationIntervals = SSA_Activations(spikes,pairs,timescale);
for j=1:size(activationIntervals,1), s = spikes(ismember(spikes(:,2),find(pairs(j,:))),:);
    code = cumsum(pairs(j,:))';
    activationIntervals{j,1}(:,3) = code(s(FindClosest(s(:,1),activationIntervals{j,1}(:,1)),2)); 
end
activations = cellfun(@(x) [mean(x(:,1:2),2) x(:,3)],activationIntervals,'UniformOutput',0); activations = Group(activations{:});
task = InIntervals(activations,postprobe);
ab = Accumulate(activations(:,[3 2]),task);
flip = find(ab(:,2)>ab(:,1));
activations(ismember(activations(:,end),flip),2) = 3-activations(ismember(activations(:,end),flip),2);

pre = InIntervals(activations,preSleep);
post = InIntervals(activations,postSleep);
task = InIntervals(activations,trials);

ab = [Accumulate(activations(:,[3 2]),pre) Accumulate(activations(:,[3 2]),task) Accumulate(activations(:,[3 2]),post)];
p = [ab(:,1)./(ab(:,1)+ab(:,2)) ab(:,3)./(ab(:,3)+ab(:,4)) ab(:,5)./(ab(:,5)+ab(:,6))];


%% Same, for shuffled spikes
rng(0);
shuffled = spikes; 
shuffled(InIntervals(spikes(:,1),preSleep),2) = Scramble(shuffled(InIntervals(spikes(:,1),preSleep),2));
shuffled(InIntervals(spikes(:,1),trials),2) = Scramble(shuffled(InIntervals(spikes(:,1),trials),2));
shuffled(InIntervals(spikes(:,1),postSleep),2) = Scramble(shuffled(InIntervals(spikes(:,1),postSleep),2));
activationIntervals = SSA_Activations(shuffled,pairs,timescale);
for j=1:size(activationIntervals,1), s = shuffled(ismember(shuffled(:,2),find(pairs(j,:))),:);
    code = cumsum(pairs(j,:))';
    activationIntervals{j,1}(:,3) = code(s(FindClosest(s(:,1),activationIntervals{j,1}(:,1)),2)); 
end
activations = cellfun(@(x) [mean(x(:,1:2),2) x(:,3)],activationIntervals,'UniformOutput',0); activations = Group(activations{:});
task = InIntervals(activations,postprobe); ab0 = Accumulate(activations(:,[3 2]),task); flip = find(ab0(:,2)>ab0(:,1));
activations(ismember(activations(:,end),flip),2) = 3-activations(ismember(activations(:,end),flip),2);

pre = InIntervals(activations,preSleep);
post = InIntervals(activations,postSleep);
task = InIntervals(activations,trials);

ab0 = [Accumulate(activations(:,[3 2]),pre) Accumulate(activations(:,[3 2]),task) Accumulate(activations(:,[3 2]),post)];
p0 = [ab0(:,1)./(ab0(:,1)+ab0(:,2)) ab0(:,3)./(ab0(:,3)+ab0(:,4)) ab0(:,5)./(ab0(:,5)+ab0(:,6))];











