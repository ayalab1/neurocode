function [errors,g,basepath,conditionType,taskType] = BatchCanCheeseReplay(basepath)

if ~exist('basepath','var') || isempty(basepath), basepath = pwd; end
basename = basenameFromBasepath(basepath);
[parentFolder,dayName] = fileparts(basepath);
[~,projectName] = fileparts(parentFolder);
sessionID = [projectName '_' dayName];

behavior = getStruct(basepath,'animal.behavior');
MergePoints = getStruct(basepath,'MergePoints');
sleep = ConsolidateIntervals(MergePoints.timestamps(cellfun(@(x) any(strfind(lower(x),'sleep')),MergePoints.foldernames),:),'epsilon',1);
try load(fullfile(basepath,'thetaCyclesTask.mat'),'cycles');
catch load(fullfile(basepath,'thetaCyclesTrack.mat'),'cycles');
end
try
    load(fullfile(basepath,[basename '.SleepState.states.mat']));
catch
    SleepScoreMaster(basepath,'noPrompts',true,'rejectchannels',[]);
    load(fullfile(basepath,[basename '.SleepState.states.mat']));
end
sws = SleepState.ints.NREMstate;
preSleep = SubtractIntervals(sleep(1,:), SubtractIntervals([0 Inf],sws));
postSleep = SubtractIntervals(sleep(2:end,:), SubtractIntervals([0 Inf],sws));

load(fullfile(basepath,[basename '_hmmLinearized.mat']));
l = behavior.position.linearized(:);
trials = behavior.trials;
run = behavior.run; nonrun = SubtractIntervals([0 Inf],run);
[spikes,regionID,regionNames,spikesCell,order] = GetAyaSpikes(pwd);
dd = sqrt((hmmLinearized.projection.x(:)-behavior.position.x(:)).^2 + (hmmLinearized.projection.y(:)-behavior.position.y(:)).^2);
t = behavior.timestamps(:);
figure; distanceThreshold = 5;
closeEnough = t(FindInterval(dd<distanceThreshold));
tooFar = t(FindInterval(dd>distanceThreshold));

pos = [t,l];
pos(:,2) = ZeroToOne(pos(:,2)); pos1 = pos;

clf
pieceSize = 0.5; step = 0.01;
x = [0 6]+15;

subplot(1,2,1);
intervals = SubtractIntervals(trials,ConsolidateIntervals(sortrows([nonrun;tooFar])));
bins = Bins(0,sum(diff(intervals,[],2)),pieceSize,step);
bins(:,1) = Unshift(bins(:,1),intervals);
bins(:,2) = Unshift(bins(:,2),intervals);
bins(diff(bins,[],2)>pieceSize*1.0001,:) = [];
[estimations,actual,errors] = ReconstructPosition(pos,spikes,bins,'training',intervals);
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

bad = IntervalsIntersect(cycles, SubtractIntervals([0 Inf],run));
bad = bad | interp1(behavior.timestamps(:),dd(:),mean(cycles,2))>10;
cycles(bad,:) = [];
[splitCycles] = SplitIntervals(cycles,'nPieces',6);
id = repmat((1:6)',length(cycles),1);
in = repelem(InIntervals(mean(cycles,2),trials) & InIntervals(mean(cycles,2),closeEnough),6,1);
subplot(1,2,2);
[estimations,actual,errors,average] = ReconstructPosition(pos1,spikes,splitCycles(in,:),'training',trials,'id',id(in));
errorArray = reshape(errors,size(errors,1),6,[]);
di = Shrink([nan;diff(actual)>0],6,1);
errorArray(:,:,di<0.5) = flipud(errorArray(:,:,di<0.5));
average = nanmean(errorArray,3);
imagesc(repmat(average,1,2))
title(['max distance ' num2str(distanceThreshold) 'cm'])

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

name = [sessionID '_OML_cheese_Cell.mat'];
[rEstimations] = ReconstructPosition(pos,spikes,rWindows,'training',ConsolidateIntervals(trials));
empty = (max(rEstimations)==min(rEstimations))';
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
outputsSaved = outputs; rEstimationsCellSaved = rEstimationsCell;

%%
for k=1:2
    bad = cellfun(@isempty,outputs(:,1));
    if k==1
        ok = ~bad & pre;
    else
        ok = ~bad & post;
    end
    seshID = ones(size(ok));
    %         subplot(2,2,k+(condition-1)*2);
    scores = cell2mat(outputs(ok,1));
    pValues = cell2mat(outputs(ok,2));
    nShuffles = size(outputs{find(ok,1),5},1);
    shuffleID = repmat((1:nShuffles)',size(scores,1),1);
    shuffledScores = cell2mat(outputs(ok,5));
    duration = durations(ok);
    slopes = (cell2mat(outputs(ok,4))-cell2mat(outputs(ok,3)))/size(average,1)*4./duration; % in m/s
    shuffledSlopes = (cell2mat(outputs(ok,7))-cell2mat(outputs(ok,6)))/size(average,1)*4./repelem(duration,nShuffles);
    normalizedShuffledSlopes = (shuffledSlopes+100)./200; % from -100 to 100 m/s
    normalizedSlopes = (slopes+100)./200; % from -100 to 100 m/s
    these{condition,k} = [scores pValues slopes seshID(ok)];
end

%%
subplot(1,3,1);
PlotColorMap(repmat(average,1,2))

%     these = saved{variety};
g0 = Group(these{:});
g = g0; % this will be normalized
% normalize replay
ok = ismember(g(:,end),[1 2]);
for i=1:max(seshID), ok = g(:,end)==1 & g(:,4)==i;
    normalizationFactor = [nanmedian(g0(ok,1)) diff(quantile(g0(ok,1),[0.25 0.75]))]; % subtract the median and divide by the quantile of the pre-sleep of the respective session
    g(ok & g(:,4)==i,1) = (g(ok & g(:,4)==i,1)-normalizationFactor(1))./normalizationFactor(2);
end

subplot(1,3,2);
ok = g0(:,end)>0;
anovabar(g(ok,1),g0(ok,end),'parametric',false)
title(['Replay: ' strrep(sessionID,'_','-') ' all bursts (simplest, use this): ' num2str(round(nanmedian(g(g(:,end)==3,1))*1000)/1000) ', ' num2str(round(1000*nanmedian(g(g(:,end)==4,1)))/1000)]);
set(gca,'tickdir','out','xtick',1:2,'xticklabel',{'on','off'},'box','off');
ylabel('score relative to baseline');

drawnow

subplot(1,6,5);
% normalise by the pre-sleep levels of significance for every session separately (to avoid inter-session influences from dominating the data)
sig = double(g(:,2)<0.05); sig(isnan(g(:,2))) = nan; ok = ismember(g(:,end),[1 2]);
for i=1:max(seshID), ok = g(:,end)==1 & g(:,4)==i; expected = nanmean(sig(ok)); ok = ismember(g(:,end),[1 2]) & g(:,4)==i;
    sig(g(:,4)==i & ok) = sig(g(:,4)==i & ok)./expected;
    meanSig(i,1) = expected; meanSig(i,1) = nanmean(sig(g(:,end)==2 & g(:,4)==i));
end
post = g(:,end)>0;
anovabar(sig(post)-1,g(post,end))
q = [nanmean(sig(g(:,end)==1)) nanmean(sig(g(:,end)==2))];
p = ranksum(sig(g(:,end)==1),sig(g(:,end)==2));
text(1,(q(1)-1)/2,num2str(round(q(1)*1000)/1000),'FontSize',12,'HorizontalAlignment','center','Rotation',90,'Color','w')
text(2,(q(2)-1)/2,num2str(round(q(2)*1000)/1000),'FontSize',12,'HorizontalAlignment','center','Rotation',90,'Color','w')
title(['p_d=' num2str(p)]);
p = [signrank(sig(g(:,end)==1)-1) signrank(sig(g(:,end)==2)-1)];
set(gca,'tickdir','out','xtick',1:2,'xticklabel',{['on,p=' num2str(p(1))],['off,p=' num2str(p(2))]},'box','off');
set(gca,'yticklabel',get(gca,'ytick')+1);
set(gca,'yticklabel',get(gca,'ytick')+1);
ylabel('proportion replay (post/baseline)');
drawnow

subplot(1,6,6);
ns = Accumulate(g0(:,end),g0(:,2)<0.05 & g0(:,3)~=0);
q = Accumulate(g0(:,end),g0(:,2)<0.05 & g0(:,3)>0);
z = zBinomialComparison(q(1),ns(1),0.5);
z(2) = zBinomialComparison(q(2),ns(2),0.5);
%     prediction = (q(3)./ns(3))/(q(1)./ns(1))*(q(2)/ns(2));
p = z2p(zBinomialComparison(q(2),ns(2),q(1),ns(1)));
q = q./ns; %q = q(3:4);%./q(1:2);
bar(q);
text(1,(q(1))/2,num2str(round(q(1)*1000)/1000),'FontSize',12,'HorizontalAlignment','center','Rotation',90,'Color','w')
text(2,(q(2))/2,num2str(round(q(2)*1000)/1000),'FontSize',12,'HorizontalAlignment','center','Rotation',90,'Color','w')
title(['p_d=' num2str(p)]);
set(gca,'tickdir','out','xtick',1:2,'xticklabel',{['on,p=' num2str(z2p(z(1)))],['off,p=' num2str(z2p(z(2)))]},'box','off');
set(gca,'yticklabel',get(gca,'ytick')+1);
ylabel('proportion forward replay (in post)');
drawnow









