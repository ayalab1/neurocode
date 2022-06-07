function [errorsOn,errorsOff,outputsOn,outputsOff,outputsBoth,pre,post,durations,basepath] = BatchCanReplay(basepath)

%% basepath = 'M:\Data\Can\OLM21\day10';
if ~exist('basepath','var') || isempty(basepath), basepath = pwd; end
basename = basenameFromBasepath(basepath);
[parentFolder,dayName] = fileparts(basepath);
[~,projectName] = fileparts(parentFolder);
sessionID = [projectName '_' dayName];

name = [sessionID '_BatchCanReplay_outputs.mat'];
try
    load(['M:\home\raly\Documents\code\new\tempFiles\' name],'errorsOn','errorsOff','outputsOn','outputsOff','outputsBoth','pre','post','durations','basepath');
    return
end
disp(['Starting session ' basepath '...']);


load(fullfile(basepath,[basename '.animal.behavior.mat']),'behavior');
load(fullfile(basepath,[basename '.cell_metrics.cellinfo.mat']),'cell_metrics');
% load(fullfile(basepath,[basename '.ripples.events.mat']),'ripples'); 
try
    stim = behavior.stimON;
catch
    % convert from old formats (for OML18 and OML19)
    load(fullfile(basepath,'stimInt.mat'), 'stimInt');
    stim = stimInt.On;
    try 
        load(fullfile(basepath,[basename '_Pos.mat']), 'Pos');
        takeColumn = (diff(sum(Pos.ON(:,2:3)>0),[],2)>0)+2;
        posOn = Pos.ON(:,[1 takeColumn]);
        posOn(:,2) = posOn(:,2)-min(posOn(:,2)); posOn(:,2) = posOn(:,2)/max(posOn(:,2)); % normalize (separately as this was badly done for some directions)
        takeColumn = (diff(sum(Pos.Off(:,2:3)>0),[],2)>0)+2;
        posOff = Pos.Off(:,[1 takeColumn]);
        posOff(:,2) = posOff(:,2)-min(posOff(:,2)); posOff(:,2) = posOff(:,2)/max(posOff(:,2)); % normalize (separately as this was badly done for some directions)
        posOff(:,2) = 1-posOff(:,2); % flip because the code below expects it
    catch  % OML18 days 5 6 7 don't have either structure, so load posTrials in those...
        load(fullfile(basepath,['posTrials.mat']), 'posTrials');
        posOn = posTrials{2}; posOff = posTrials{1};
    end
    q = sortrows([posOn; posOff]);
    behavior.position.x = q(:,2)';
    behavior.position.y = ones(size(q(:,2)))';
    behavior.positionTrials{1} = posOn;
    behavior.positionTrials{2} = posOff;
    isBoundary = zscore(diff(posOn(:,1)))>1;
    trialsOn = [posOn([true; isBoundary],1) posOn([isBoundary;true],1)];
    isBoundary = zscore(diff(posOff(:,1)))>1;
    trialsOff = [posOff([true; isBoundary],1) posOff([isBoundary;true],1)];
    trials = Group(trialsOn,trialsOff);
    trials(:,3) = 2-trials(:,3);
    behavior.trials = trials(:,1:2);
    behavior.trialID = trials(:,[3 3]);
end

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

try
    run = behavior.run;
catch
    % Old way to define running:
    ok = ~isnan(behavior.position.x(:)) & ~isnan(behavior.position.y(:));
    interpolated = interp1(behavior.timestamps(ok)',[behavior.timestamps(ok)' behavior.position.x(ok)',behavior.position.y(ok)'],behavior.timestamps');
    interpolated(isnan(interpolated(:,1)),:) = [];
    speed = LinearVelocity(interpolated,5); t = speed(:,1);
    allPositions = sortrows([behavior.positionTrials{1}; behavior.positionTrials{2}]);
    [up,down] = ZeroCrossings([allPositions(:,1),allPositions(:,2)-0.5]);
    midpoints = allPositions(up|down,1);
    topSpeed = interp1(speed(:,1),speed(:,2),midpoints);
    threshold = nanmedian(topSpeed)/10; % the threshold is 10% of the peak speed
    run = t(FindInterval(speed(:,2)>threshold));
    run = ConsolidateIntervals(run,'epsilon',0.01);
    [in,w] = InIntervals(behavior.timestamps(:),run);
    peak = Accumulate(w(in),behavior.speed(in)','mode','max');

    % remove outliers (data in between sessions gives outlier speeds)
    [~,isOutlier] = RemoveOutliers(peak);
    % remove run epochs that don't reach the speed threshold
    run(peak<0.1 | isOutlier,:) = [];
    run(IntervalsIntersect(run,sleep),:) = [];
end


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

% Make sure the saved linearized positions go from 0 to 1
minimum = min([min(behavior.positionTrials{1}(:,2)) min(behavior.positionTrials{2}(:,2))]);
maximum = max([max(behavior.positionTrials{1}(:,2)) max(behavior.positionTrials{2}(:,2))]);
behavior.positionTrials{1}(:,2) = (behavior.positionTrials{1}(:,2) - minimum)./(maximum-minimum);
behavior.positionTrials{2}(:,2) = (behavior.positionTrials{2}(:,2) - minimum)./(maximum-minimum);
pos0 = sortrows([behavior.positionTrials{1}; behavior.positionTrials{2}(:,1) 2-behavior.positionTrials{2}(:,2)]);
if sum(diff(pos0(:,2))<0) > sum(diff(pos0(:,2))>0) % backwards, needs to be flipped!
    pos0(:,2) = 2-pos0(:,2);
end
pos1 = Restrict(pos0(pos0(:,2)<1,:),run);
pos2 = Restrict(pos0(pos0(:,2)>1,:),run); pos2(:,2) = pos2(:,2)-1;
try
    onTrials = behavior.trials(behavior.trialID(:,2)==1,:); % Can promises that behavior.trialID(:,2)==1 is always stim and 0 is always nonstim trials
    offTrials = behavior.trials(behavior.trialID(:,2)==0,:);
catch
    onTrials = behavior.trials(IntervalsIntersect(behavior.trials,stim),:);
    offTrials = behavior.trials(~IntervalsIntersect(behavior.trials,stim),:);
end

% the first direction is stim (for labels in later analyses)
if sum(InIntervals(pos2,stim)) > sum(InIntervals(pos1,stim)) % otherwise, swap them
    this = pos1; pos1 = pos2; pos2 = this;
end

stimIntervals = SubtractIntervals(run,SubtractIntervals([0 Inf],stim)); % only stim run periods
nonstimIntervals = SubtractIntervals(run,stim); 
nonstimIntervals = SubtractIntervals(nonstimIntervals,SubtractIntervals([0 Inf],offTrials)); % make sure nonstimIntervals only take place during offTrials
pos = sortrows([Restrict(pos1,stimIntervals) ; Restrict([pos2(:,1) 1+pos2(:,2)],nonstimIntervals)]);
pos(:,2) = pos(:,2)./2;

% compute firing curves
clear curves curves1 curves2
for i=1:max(spikes(:,2))
    map = FiringCurve(Restrict(pos1,stimIntervals,'shift','on'),Restrict(spikes(spikes(:,2)==i),stimIntervals,'shift','on'));
    curves1(i,:) = map.rate;
    map = FiringCurve(Restrict(pos2,nonstimIntervals,'shift','on'),Restrict(spikes(spikes(:,2)==i),nonstimIntervals,'shift','on'));
    curves2(i,:) = map.rate;
    map = FiringCurve(Restrict(pos,run,'shift','on'),Restrict(spikes(spikes(:,2)==i),run,'shift','on'));
    curves(i,:) = map.rate;
end

clf
[~,m] = max(Smooth(curves,[0 2],'type','cc'),[],2);
PlotColorMap(sortby(nanzscore([curves],[],2),m))
PlotColorMap(sortby(nanzscore([curves1 curves2 curves],[],2),m))

SaveFig(fullfile('M:\home\raly\results\OML\fields',[sessionID '-Place_field_coverage']))

%%
% figure; 
clf
try load(fullfile(basepath,'thetaCyclesTask.mat'),'cycles');
catch load(fullfile(basepath,'thetaCyclesTrack.mat'),'cycles');
end
bad = IntervalsIntersect(cycles, SubtractIntervals([0 Inf],run));
cycles(bad,:) = [];

% r = ripples.timestamps(:,1:2);
[splitCycles] = SplitIntervals(cycles,'nPieces',6);
id = repmat((1:6)',length(cycles),1);
in = repelem(InIntervals(cycles,stimIntervals),6,1);

[estimations,actual,errors,average] = ReconstructPosition(pos1,spikes,splitCycles(in,:),'training',stimIntervals,'id',id(in));
errorsOn = errors;
on = nanmean(reshape(errors,size(errors,1),6,[]),3);

ok = repelem(InIntervals(cycles,nonstimIntervals),6,1);
[estimations,actual,errors,average] = ReconstructPosition(pos2,spikes,splitCycles(ok,:),'training',nonstimIntervals,'id',id(ok));
errorsOff = errors;
off = nanmean(reshape(errors,size(errors,1),6,[]),3);
% [estimations,actual,errors,average] = ReconstructPosition(pos,spikes,splitCycles(in,:),'training',run,'id',id(in));
% on = nanmean(reshape(errors,size(errors,1),6,[]),3);
% 
% [estimations,actual,errors,average] = ReconstructPosition(pos,spikes,splitCycles(~in,:),'training',run,'id',id(~in));
% off = nanmean(reshape(errors,size(errors,1),6,[]),3);

subplot(1,2,1);
average = on;
[score,p,a,b,~,~,~,c,cShuffled] = FindReplayScore(average,'circular','off','wcorr','on');
PlotColorMap(repmat(average,1,2));
hold on; plot([1 6 6.5 [1 6]+6],[a b nan a b],'k--','linewidth',2);
plot([1 6 6.5 [1 6]+6],[a b nan a b]-15,'k','linewidth',1); plot([1 6 6.5 [1 6]+6],[a b nan a b]+15,'k','linewidth',1);
set(gca,'ytick',100,'yticklabel','0','xtick','');
PlotHVLines(100,'h','w--','linewidth',2);
ylabel('(decoded position) - (current position)');
title(strrep([sessionID ' stim ON'],'_','-'));

subplot(1,2,2);
average = off;
[score,p,a,b,~,~,~,c,cShuffled] = FindReplayScore(average,'circular','off','wcorr','on');
PlotColorMap(repmat(average,1,2));
hold on;
plot([1 6 6.5 [1 6]+6],[a b nan a b],'k--','linewidth',2);
plot([1 6 6.5 [1 6]+6],[a b nan a b]-15,'k','linewidth',1); plot([1 6 6.5 [1 6]+6],[a b nan a b]+15,'k','linewidth',1);
set(gca,'ytick',100,'yticklabel','0','xtick','');
PlotHVLines(100,'h','w--','linewidth',2);
ylabel('(decoded position) - (current position)');
title(strrep([sessionID ' stim OFF'],'_','-'));

% ylim([50.5 150.5]);

clims
SaveFig(fullfile('M:\home\raly\results\OML\theta',[sessionID '-theta-sequences-On-Off']))

%% Check if spikes are enough to reconstruct current position

clf
pieceSize = 0.5; step = 0.01;
x = [0 6]+15;
subplot(1,2,1);
% on
% intervals = stimIntervals;
intervals = SubtractIntervals(onTrials,SubtractIntervals([0 Inf],run));
bins = Bins(0,sum(diff(intervals,[],2)),pieceSize,step);
bins(:,1) = Unshift(bins(:,1),intervals);
bins(:,2) = Unshift(bins(:,2),intervals);
bins(diff(bins,[],2)>pieceSize*1.0001,:) = [];
[estimations,actual,errors] = ReconstructPosition(pos1,spikes,bins,'training',intervals);
relative_t = (1:size(estimations,2))'*step;
imagesc(estimations,'x',relative_t);
set(gca,'YDir','normal','tickdir','out');
hold all
plot(relative_t,actual*200,'k.-')
clim([0 0.05])
xlim(x);
ColorMap(gca,[1 1 1],[1 0 0]);
xlabel('time (s) (time is restricted to on-trials)');
title(['on trials: ' num2str(pieceSize*1000) 'ms widnows and ' num2str(step*1000) 'ms step.']);
% cycles_relative_t = relative_t(FindClosest(mean(bins,2),cycles));
% PlotHVLines(Restrict(cycles_relative_t,xlim),'v','k--');

% x = [19 25];
subplot(1,2,2);
% off
intervals = nonstimIntervals;
bins = Bins(0,sum(diff(intervals,[],2)),pieceSize,step);
bins(:,1) = Unshift(bins(:,1),intervals);
bins(:,2) = Unshift(bins(:,2),intervals);
bins(diff(bins,[],2)>pieceSize*1.0001,:) = [];
[estimations,actual,errors] = ReconstructPosition(pos2,spikes,bins,'training',intervals);
relative_t = (1:size(estimations,2))'*step;

imagesc(estimations,'x',relative_t);
set(gca,'YDir','normal','tickdir','out');
hold all
plot(relative_t,actual*200,'k.-')
clim([0 0.05])
xlim(x);
ColorMap(gca,[1 1 1],[0 0 1]);
xlabel('time (s) (time is restricted to off-trials)');
title(['off trials: ' num2str(pieceSize*1000) 'ms widnows and ' num2str(step*1000) 'ms step.']);
% cycles_relative_t = relative_t(FindClosest(mean(bins,2),cycles));
% PlotHVLines(Restrict(cycles_relative_t,xlim),'v','k--');

SaveFig(fullfile('M:\home\raly\results\OML\reconstructionQuality',[sessionID '-General-Reconstruction-Quality']))

%% Replay
clf
% bursts = FindBursts(spikes,'thresholds',[0 3],'smooth',0.01);
bursts = FindBursts(Restrict(spikes,sws,'shift','on'),'thresholds',[0 3],'smooth',0.01);
bursts(:) = Unshift(bursts(:),sws);
duration = diff(bursts(:,[1 3]),[],2);
bursts(duration>1,:) = [];

pre = InIntervals(bursts,preSleep);
post = InIntervals(bursts,postSleep);

r = bursts(:,[1 3]); leeway = 0.01;
durations = diff(r,[],2);
intervals = r; intervals0 = intervals;
intervals(1:end-1,2) = min([intervals0(1:end-1,2)+leeway mean([intervals0(1:end-1,2) intervals0(2:end,1)],2)],[],2); intervals(end,2) = intervals0(end,2) + leeway;
intervals(2:end,1) = max([intervals0(2:end,1)-leeway mean([intervals0(2:end,1) intervals0(1:end-1,2)],2)],[],2); intervals(1,1) = intervals0(1,1) - leeway;
[rWindows,rID] = SplitIntervals(intervals,'pieceSize',0.02);

% % First, off
% nonstimIntervals = SubtractIntervals(run,stim);
% [rEstimations] = ReconstructPosition(pos2,spikes,rWindows,'training',nonstimIntervals);
% empty = (max(rEstimations)==min(rEstimations))';
% [score,pValue,seqStartStop,stops] = deal(nan(size(r,1),1)); seqStartStop(:,2) = nan;
% rEstimationsCell = cell(size(r,1),1);
% for i=1:size(r,1)
%     rEstimationsCell{i,1} = rEstimations(:,rID==i);
% end
% parfor i=1:length(r)
%     if sum(rID==i & ~empty)>=3
%         [score(i),pValue(i),seqStartStop(i),stops(i)] = FindReplayScore(rEstimationsCell{i},'circular','off');
%     end
%     if rem(i,100)==0, display([num2str(i) '/' num2str(size(r,1))]); end
% end
% seqStartStop(:,2) = stops;
% % save(['M:\home\raly\Documents\code\new\tempFiles\' name],'estimations','actual','errors','average','rEstimations','score','pValue','seqStartStop','bursts');

% First, off
name = [sessionID '_OML_off_Cell.mat'];
[rEstimations] = ReconstructPosition(pos2,spikes,rWindows,'training',nonstimIntervals);
empty = (max(rEstimations)==min(rEstimations))';
rEstimationsCell = cell(size(r,1),1);
for i=1:size(r,1)
    rEstimationsCell{i,1} = rEstimations(:,rID==i);
end

tic
try
    load(['M:\home\raly\Documents\code\new\tempFiles\' name],'outputs');
    if size(outputs,1)~=size(rEstimationsCell,1)
        disp(['Saved replay scores don''t fit the detected events. Recomputing...']);
        error(['Saved replay scores don''t fit the detected events. Recomputing...']);
    end
catch
%     error('!')
    outputs = cell(size(r,1),13);
    parfor i=1:size(r,1)
        if sum(rID==i & ~empty)>=3
            outputs(i,:) = FindReplayScoreCell(rEstimationsCell{i},'circular','off');
        end
        if rem(i,100)==0, display([num2str(i) '/' num2str(size(r,1))]); end
    end
    save(['M:\home\raly\Documents\code\new\tempFiles\' name],'bursts','rEstimationsCell','outputs','intervals','pos2');
end
toc
outputsOff = outputs; rEstimationsCellOff = rEstimationsCell;

% Then, on
name = [sessionID '_OML_on_Cell.mat'];
[rEstimations] = ReconstructPosition(pos1,spikes,rWindows,'training',stimIntervals);
tic
rEstimationsCell = cell(size(r,1),1);
for i=1:size(r,1)
    rEstimationsCell{i,1} = rEstimations(:,rID==i);
end
try
    load(['M:\home\raly\Documents\code\new\tempFiles\' name],'outputs');
    if size(outputs,1)~=size(rEstimationsCell,1)
        disp(['Saved replay scores don''t fit the detected events. Recomputing...']);
        error(['Saved replay scores don''t fit the detected events. Recomputing...']);
    end
catch
%     error('!')
    outputs = cell(size(r,1),13);
    parfor i=1:size(r,1)
        if sum(rID==i & ~empty)>=3
            outputs(i,:) = FindReplayScoreCell(rEstimationsCell{i},'circular','off');
        end
        if rem(i,100)==0, display([num2str(i) '/' num2str(size(r,1))]); end
    end
    save(['M:\home\raly\Documents\code\new\tempFiles\' name],'bursts','rEstimationsCell','outputs','intervals','pos1');
end
toc
outputsOn = outputs; rEstimationsCellOn = rEstimationsCell;

% Together (merged positions, on is from 0 to 0.5, off is from 0.5 to 1)
name = [sessionID '_OML_Cell.mat'];
[rEstimations] = ReconstructPosition(pos,spikes,rWindows,'training',run);
tic
rEstimationsCell = cell(size(r,1),1);
for i=1:size(r,1)
    rEstimationsCell{i,1} = rEstimations(:,rID==i);
end
try
    load(['M:\home\raly\Documents\code\new\tempFiles\' name],'outputs');
    if size(outputs,1)~=size(rEstimationsCell,1)
        disp(['Saved replay scores don''t fit the detected events. Recomputing...']);
        error(['Saved replay scores don''t fit the detected events. Recomputing...']);
    end
catch
%     error('!')
    outputs = cell(size(r,1),13);
    parfor i=1:size(r,1)
        if sum(rID==i & ~empty)>=3
            outputs(i,:) = FindReplayScoreCell(rEstimationsCell{i},'circular','off','threshold',10);
        end
        if rem(i,100)==0, display([num2str(i) '/' num2str(size(r,1))]); end
    end
    save(['M:\home\raly\Documents\code\new\tempFiles\' name],'bursts','rEstimationsCell','outputs','intervals','pos');
end
toc
outputsBoth = outputs; rEstimationsCellBoth = rEstimationsCell;

try
    name = [sessionID '_BatchCanReplay_outputs.mat'];
    save(['M:\home\raly\Documents\code\new\tempFiles\' name],'errorsOn','errorsOff','outputsOn','outputsOff','outputsBoth','pre','post','durations','basepath');
    disp(['Saved the results for session ' basepath '...']);
end

% outsputs are [r,p,st,sp,rShuffled,aShuffled,bShuffled,c,cShuffled,jump,jumpShuffled,maxJump,maxJumpShuffled]

%%
nBins = [101 25];
bad = cellfun(@isempty,outputs(:,1));
outputs = outputsOff;
smooth = [0 1];
clf
reOn = false(size(outputsOn,1),1);
reOff = false(size(outputsOn,1),1);
for i=1:size(outputsOn,1), try reOn(i,1) = outputsOn{i,1}>outputsOff{i,1}; end; end
for i=1:size(outputsOn,1), try reOff(i,1) = outputsOn{i,1}<outputsOff{i,1}; end; end
conditionNames = {strrep([sessionID ' ON'],'_','-'),strrep([sessionID ' OFF'],'_','-')};
periodNames = {'pre-task sleep','post-task sws'};
for condition = 1:2
    if condition==1
        outputs = outputsOn;
        bad = cellfun(@isempty,outputs(:,1));
        bad = bad | ~(reOn);
    else
        outputs = outputsOff;
        bad = cellfun(@isempty,outputs(:,1));
        bad = bad | ~(reOff);
    end
    for k=1:2
        if k==1
            ok = ~bad & pre;
        else
            ok = ~bad & post;
        end
        subplot(2,2,k+(condition-1)*2);
        scores = cell2mat(outputs(ok,1));
        nShuffles = size(outputs{find(ok,1),5},1);
        shuffleID = repmat((1:nShuffles)',size(scores,1),1);
        shuffledScores = cell2mat(outputs(ok,5));
        duration = diff(r(ok,:),[],2);
        slopes = (cell2mat(outputs(ok,4))-cell2mat(outputs(ok,3)))/size(rEstimations,1)*4./duration; % in m/s
        shuffledSlopes = (cell2mat(outputs(ok,7))-cell2mat(outputs(ok,6)))/size(rEstimations,1)*4./repelem(duration,nShuffles);
        normalizedShuffledSlopes = (shuffledSlopes+100)./200; % from -100 to 100 m/s
        normalizedSlopes = (slopes+100)./200; % from -100 to 100 m/s

        % slopes = cell2mat(outputs(ok,4))-cell2mat(outputs(ok,3));
        % shuffledSlopes = cell2mat(outputs(ok,7))-cell2mat(outputs(ok,6));
        % normalizedShuffledSlopes = (shuffledSlopes+size(rEstimations,1))/(size(rEstimations,1)*2);
        % normalizedSlopes = (slopes+size(rEstimations,1))/(size(rEstimations,1)*2);

        mapsShuffled = nan([nBins([2 1]) nShuffles]);
        for i=1:nShuffles
            mapsShuffled(:,:,i) = DensityMap(shuffledSlopes(shuffleID==i),shuffledScores(shuffleID==i)*0.8,'nBins',nBins,'smooth',0,'show','off');
        end
        map = DensityMap(normalizedSlopes,scores*0.8,'nBins',nBins,'smooth',0,'show','off');
        mShuffled = Smooth(nanmean(mapsShuffled,3),0); stdShuffled = Smooth(nanstd(mapsShuffled,[],3),0);
        zmap = (map-mShuffled)./stdShuffled; zmap(zmap==Inf) = max(zmap(zmap~=Inf));
        smoothed = zmap; smoothed(isnan(zmap)) = 0; 
        PlotColorMap(Smooth(smoothed,smooth),'x',linspace(-100,100,nBins(1)),'y',linspace(0,1/0.8,nBins(2)));
        PlotHVLines(0,'v','w--','linewidth',2);
        xlim([-1 1]*20); ylim([0 1]);
        title([conditionNames{condition} ' ' periodNames{k}]);
        ylabel('score'); xlabel('slope (m/s)');
        drawnow
    end
end
clims([0 20]);

SaveFig(fullfile('M:\home\raly\results\OML\histograms2d',[sessionID '-JointReplayScoreSlope-On-Off']))

%% Barplots
clf

for condition = 1:2
    if condition==1
        outputs = outputsOn;
        bad = cellfun(@isempty,outputs(:,1));
        bad = bad;% | ~(reOn);
    else
        outputs = outputsOff;
        bad = cellfun(@isempty,outputs(:,1));
        bad = bad;% | ~(reOff);
    end
    for k=1:2
        if k==1
            ok = ~bad & pre;
        else
            ok = ~bad & post;
        end
%         subplot(2,2,k+(condition-1)*2);
        scores = cell2mat(outputs(ok,1));
        pValues = cell2mat(outputs(ok,2));
        nShuffles = size(outputs{find(ok,1),5},1);
        shuffleID = repmat((1:nShuffles)',size(scores,1),1);
        shuffledScores = cell2mat(outputs(ok,5));
        duration = diff(r(ok,:),[],2);
        slopes = (cell2mat(outputs(ok,4))-cell2mat(outputs(ok,3)))/size(rEstimations,1)*4./duration; % in m/s
        shuffledSlopes = (cell2mat(outputs(ok,7))-cell2mat(outputs(ok,6)))/size(rEstimations,1)*4./repelem(duration,nShuffles);
        normalizedShuffledSlopes = (shuffledSlopes+100)./200; % from -100 to 100 m/s
        normalizedSlopes = (slopes+100)./200; % from -100 to 100 m/s
        these{condition,k} = [scores pValues slopes];
    end
end

g = Group(these{:});
gg = g; gg(ismember(gg(:,end),[1 3]),1) = (gg(ismember(gg(:,end),[1 3]),1) - nanmedian(gg(gg(:,end)==1,1)))./diff(quantile(gg(gg(:,end)==1,1),[0.25 0.75]));
gg(ismember(gg(:,end),[2 4]),1) = (gg(ismember(gg(:,end),[2 4]),1) - nanmedian(gg(gg(:,end)==2,1)))./diff(quantile(gg(gg(:,end)==2,1),[0.25 0.75]));

subplot(3,2,1);
ok = g(:,end)>2;
anovabar(gg(ok,1),g(ok,end)-2,'parametric',false)
title([strrep(sessionID,'_','-') ' all bursts (simplest, use this): ' num2str(round(nanmedian(gg(gg(:,end)==3,1))*1000)/1000) ', ' num2str(round(1000*nanmedian(gg(gg(:,end)==4,1)))/1000)]);
set(gca,'tickdir','out','xtick',1:2,'xticklabel',{'on','off'},'box','off');
ylabel('score relative to baseline');
% subplot(2,2,2);
% q = Accumulate(g(:,end),g(:,2)<0.05)./Accumulate(g(:,end),1); 
% q = Accumulate(g(:,end),g(:,3)>0)./Accumulate(g(:,end),1); 
% q = q(3:4)./q(1:2); bar(q-1);
% title(['all bursts: ' num2str(round(q(1)*1000)/1000) ', ' num2str(round(1000*q(2))/1000)]);
% set(gca,'tickdir','out','xtick',1:2,'xticklabel',{'on','off'},'box','off');
% set(gca,'yticklabel',get(gca,'ytick')+1);
% ylabel('proportion replay (post/baseline)');
% 

subplot(3,4,3);
ns = Accumulate(g(:,end),1);
q = Accumulate(g(:,end),g(:,2)<0.05);
z = zBinomialComparison(q(3),ns(3),q(1),ns(1));
z(2) = zBinomialComparison(q(4),ns(4),q(2),ns(2));
prediction = (q(3)./ns(3))/(q(1)./ns(1))*(q(2)/ns(2));
p = z2p(zBinomialComparison(q(4),ns(4),prediction));
q = q./ns; q = q(3:4)./q(1:2); 
bar(q-1);
text(1,(q(1)-1)/2,num2str(round(q(1)*1000)/1000),'FontSize',12,'HorizontalAlignment','center','Rotation',90,'Color','w')
text(2,(q(2)-1)/2,num2str(round(q(2)*1000)/1000),'FontSize',12,'HorizontalAlignment','center','Rotation',90,'Color','w')
title(['p_d=' num2str(p)]);
set(gca,'tickdir','out','xtick',1:2,'xticklabel',{['on,p=' num2str(z2p(z(1)))],['off,p=' num2str(z2p(z(2)))]},'box','off');
set(gca,'yticklabel',get(gca,'ytick')+1);
set(gca,'yticklabel',get(gca,'ytick')+1);
ylabel('proportion replay (post/baseline)');

subplot(3,4,4);
ns = Accumulate(g(:,end),1);
q = Accumulate(g(:,end),g(:,2)<0.05 & g(:,3)>0);
z = zBinomialComparison(q(3),ns(3),q(1),ns(1));
z(2) = zBinomialComparison(q(4),ns(4),q(2),ns(2));
prediction = (q(3)./ns(3))/(q(1)./ns(1))*(q(2)/ns(2));
p = z2p(zBinomialComparison(q(4),ns(4),prediction));
q = q./ns; q = q(3:4)./q(1:2); 
bar(q-1);
text(1,(q(1)-1)/2,num2str(round(q(1)*1000)/1000),'FontSize',12,'HorizontalAlignment','center','Rotation',90,'Color','w')
text(2,(q(2)-1)/2,num2str(round(q(2)*1000)/1000),'FontSize',12,'HorizontalAlignment','center','Rotation',90,'Color','w')
title(['p_d=' num2str(p)]);
set(gca,'tickdir','out','xtick',1:2,'xticklabel',{['on,p=' num2str(z2p(z(1)))],['off,p=' num2str(z2p(z(2)))]},'box','off');
set(gca,'yticklabel',get(gca,'ytick')+1);
ylabel('proportion forward replay (post/baseline)');

reOn = false(size(outputsOn,1),1);
reOff = false(size(outputsOn,1),1);
for i=1:size(outputsOn,1), try reOn(i,1) = outputsOn{i,1}>outputsOff{i,1}; end; end
for i=1:size(outputsOn,1), try reOff(i,1) = outputsOn{i,1}<outputsOff{i,1}; end; end
for condition = 1:2
    if condition==1
        outputs = outputsOn;
        bad = cellfun(@isempty,outputs(:,1));
        bad = bad | ~(reOn);
    else
        outputs = outputsOff;
        bad = cellfun(@isempty,outputs(:,1));
        bad = bad | ~(reOff);
    end
    for k=1:2
        if k==1
            ok = ~bad & pre;
        else
            ok = ~bad & post;
        end
%         subplot(2,2,k+(condition-1)*2);
        scores = cell2mat(outputs(ok,1));
        pValues = cell2mat(outputs(ok,2));
        nShuffles = size(outputs{find(ok,1),5},1);
        shuffleID = repmat((1:nShuffles)',size(scores,1),1);
        shuffledScores = cell2mat(outputs(ok,5));
        duration = diff(r(ok,:),[],2);
        slopes = (cell2mat(outputs(ok,4))-cell2mat(outputs(ok,3)))/size(rEstimations,1)*4./duration; % in m/s
        shuffledSlopes = (cell2mat(outputs(ok,7))-cell2mat(outputs(ok,6)))/size(rEstimations,1)*4./repelem(duration,nShuffles);
        normalizedShuffledSlopes = (shuffledSlopes+100)./200; % from -100 to 100 m/s
        normalizedSlopes = (slopes+100)./200; % from -100 to 100 m/s
        these{condition,k} = [scores pValues slopes];
    end
end

g = Group(these{:});
gg = g; gg(ismember(gg(:,end),[1 3]),1) = (gg(ismember(gg(:,end),[1 3]),1) - nanmedian(gg(gg(:,end)==1,1)))./diff(quantile(gg(gg(:,end)==1,1),[0.25 0.75]));
gg(ismember(gg(:,end),[2 4]),1) = (gg(ismember(gg(:,end),[2 4]),1) - nanmedian(gg(gg(:,end)==2,1)))./diff(quantile(gg(gg(:,end)==2,1),[0.25 0.75]));

subplot(3,2,3);
ok = g(:,end)>2;
anovabar(gg(ok,1),g(ok,end)-2,'parametric',false)
title(['better fit (on/off): ' num2str(round(nanmedian(gg(gg(:,end)==3,1))*1000)/1000) ', ' num2str(round(1000*nanmedian(gg(gg(:,end)==4,1)))/1000)]);
set(gca,'tickdir','out','xtick',1:2,'xticklabel',{'on','off'},'box','off');
ylabel('score relative to baseline');

subplot(3,4,7);
ns = Accumulate(g(:,end),1);
q = Accumulate(g(:,end),g(:,2)<0.05);
z = zBinomialComparison(q(3),ns(3),q(1),ns(1));
z(2) = zBinomialComparison(q(4),ns(4),q(2),ns(2));
prediction = (q(3)./ns(3))/(q(1)./ns(1))*(q(2)/ns(2));
p = z2p(zBinomialComparison(q(4),ns(4),prediction));
q = q./ns; q = q(3:4)./q(1:2); 
bar(q-1);
text(1,(q(1)-1)/2,num2str(round(q(1)*1000)/1000),'FontSize',12,'HorizontalAlignment','center','Rotation',90,'Color','w')
text(2,(q(2)-1)/2,num2str(round(q(2)*1000)/1000),'FontSize',12,'HorizontalAlignment','center','Rotation',90,'Color','w')
title(['p_d=' num2str(p)]);
set(gca,'tickdir','out','xtick',1:2,'xticklabel',{['on,p=' num2str(z2p(z(1)))],['off,p=' num2str(z2p(z(2)))]},'box','off');
set(gca,'yticklabel',get(gca,'ytick')+1);
ylabel('proportion replay (post/baseline)');


subplot(3,4,8);
ns = Accumulate(g(:,end),1);
q = Accumulate(g(:,end),g(:,2)<0.05 & g(:,3)>0);
z = zBinomialComparison(q(3),ns(3),q(1),ns(1));
z(2) = zBinomialComparison(q(4),ns(4),q(2),ns(2));
prediction = (q(3)./ns(3))/(q(1)./ns(1))*(q(2)/ns(2));
p = z2p(zBinomialComparison(q(4),ns(4),prediction));
q = q./ns; q = q(3:4)./q(1:2); 
bar(q-1);
text(1,(q(1)-1)/2,num2str(round(q(1)*1000)/1000),'FontSize',12,'HorizontalAlignment','center','Rotation',90,'Color','w')
text(2,(q(2)-1)/2,num2str(round(q(2)*1000)/1000),'FontSize',12,'HorizontalAlignment','center','Rotation',90,'Color','w')
title(['p_d=' num2str(p)]);
set(gca,'tickdir','out','xtick',1:2,'xticklabel',{['on,p=' num2str(z2p(z(1)))],['off,p=' num2str(z2p(z(2)))]},'box','off');
set(gca,'yticklabel',get(gca,'ytick')+1);
ylabel('proportion forward replay (post/baseline)');

% cla
% anovabar(g(:,3)>0,g(:,end));


outputs = outputsBoth;
reOn = false(size(outputs,1),1);
reOff = false(size(outputs,1),1);
for i=1:size(outputsOn,1), try reOn(i,1) = outputs{i,3}<=size(rEstimations,1)/2; end; end
for i=1:size(outputsOn,1), try reOff(i,1) = outputs{i,3}>size(rEstimations,1)/2; end; end
for condition = 1:2
    if condition==1
        bad = cellfun(@isempty,outputs(:,1));
        bad = bad | ~(reOn);
    else
        bad = cellfun(@isempty,outputs(:,1));
        bad = bad | ~(reOff);
    end
    for k=1:2
        if k==1
            ok = ~bad & pre;
        else
            ok = ~bad & post;
        end
%         subplot(2,2,k+(condition-1)*2);
        scores = cell2mat(outputs(ok,1));
        pValues = cell2mat(outputs(ok,2));
        nShuffles = size(outputs{find(ok,1),5},1);
        shuffleID = repmat((1:nShuffles)',size(scores,1),1);
        shuffledScores = cell2mat(outputs(ok,5));
        duration = diff(r(ok,:),[],2);
        slopes = (cell2mat(outputs(ok,4))-cell2mat(outputs(ok,3)))/size(rEstimations,1)*4./duration; % in m/s
        shuffledSlopes = (cell2mat(outputs(ok,7))-cell2mat(outputs(ok,6)))/size(rEstimations,1)*4./repelem(duration,nShuffles);
        normalizedShuffledSlopes = (shuffledSlopes+100)./200; % from -100 to 100 m/s
        normalizedSlopes = (slopes+100)./200; % from -100 to 100 m/s
        these{condition,k} = [scores pValues slopes];
    end
end

g = Group(these{:});
gg = g; gg(ismember(gg(:,end),[1 3]),1) = (gg(ismember(gg(:,end),[1 3]),1) - nanmedian(gg(gg(:,end)==1,1)))./diff(quantile(gg(gg(:,end)==1,1),[0.25 0.75]));
gg(ismember(gg(:,end),[2 4]),1) = (gg(ismember(gg(:,end),[2 4]),1) - nanmedian(gg(gg(:,end)==2,1)))./diff(quantile(gg(gg(:,end)==2,1),[0.25 0.75]));

subplot(3,2,5);
ok = g(:,end)>2;
anovabar(gg(ok,1),g(ok,end)-2,'parametric',false)
title(['decoded together (sorted posthoc): ' num2str(round(nanmedian(gg(gg(:,end)==3,1))*1000)/1000) ', ' num2str(round(1000*nanmedian(gg(gg(:,end)==4,1)))/1000)]);
set(gca,'tickdir','out','xtick',1:2,'xticklabel',{'on','off'},'box','off');
ylabel('score relative to baseline');

subplot(3,4,11);
ns = Accumulate(g(:,end),1);
q = Accumulate(g(:,end),g(:,2)<0.05);
z = zBinomialComparison(q(3),ns(3),q(1),ns(1));
z(2) = zBinomialComparison(q(4),ns(4),q(2),ns(2));
prediction = (q(3)./ns(3))/(q(1)./ns(1))*(q(2)/ns(2));
p = z2p(zBinomialComparison(q(4),ns(4),prediction));
q = q./ns; q = q(3:4)./q(1:2); 
bar(q-1);
text(1,(q(1)-1)/2,num2str(round(q(1)*1000)/1000),'FontSize',12,'HorizontalAlignment','center','Rotation',90,'Color','w')
text(2,(q(2)-1)/2,num2str(round(q(2)*1000)/1000),'FontSize',12,'HorizontalAlignment','center','Rotation',90,'Color','w')
title(['p_d=' num2str(p)]);
set(gca,'tickdir','out','xtick',1:2,'xticklabel',{['on,p=' num2str(z2p(z(1)))],['off,p=' num2str(z2p(z(2)))]},'box','off');
set(gca,'yticklabel',get(gca,'ytick')+1);
ylabel('proportion replay (post/baseline)');


subplot(3,4,12);
ns = Accumulate(g(:,end),1);
q = Accumulate(g(:,end),g(:,2)<0.05 & g(:,3)>0);
z = zBinomialComparison(q(3),ns(3),q(1),ns(1));
z(2) = zBinomialComparison(q(4),ns(4),q(2),ns(2));
prediction = (q(3)./ns(3))/(q(1)./ns(1))*(q(2)/ns(2));
p = z2p(zBinomialComparison(q(4),ns(4),prediction));
q = q./ns; q = q(3:4)./q(1:2); 
bar(q-1);
text(1,(q(1)-1)/2,num2str(round(q(1)*1000)/1000),'FontSize',12,'HorizontalAlignment','center','Rotation',90,'Color','w')
text(2,(q(2)-1)/2,num2str(round(q(2)*1000)/1000),'FontSize',12,'HorizontalAlignment','center','Rotation',90,'Color','w')
title(['p_d=' num2str(p)]);
set(gca,'tickdir','out','xtick',1:2,'xticklabel',{['on,p=' num2str(z2p(z(1)))],['off,p=' num2str(z2p(z(2)))]},'box','off');
set(gca,'yticklabel',get(gca,'ytick')+1);
ylabel('proportion forward replay (post/baseline)');


subplot(3,2,1); y = ylim; subplot(3,2,3); y = [min([min(ylim) y(1)]) max([max(ylim) y(2)])]; ylim(y); subplot(3,2,1); ylim(y);
subplot(3,4,3); y = ylim; 
for k=[4 7 8 11 12], subplot(3,4,k); y = [min([min(ylim) y(1)]) max([max(ylim) y(2)])]; end
for k=[3 4 7 8 11 12], subplot(3,4,k); ylim(y); set(gca,'yticklabel',get(gca,'ytick')+1); end

drawnow
SaveFig(fullfile('M:\home\raly\results\OML\replay',[sessionID '-Replay-barplots']))





















