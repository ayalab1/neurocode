%% Basic pipeline for replay
% These are the basic functions required to do replay analysis on linearized
% positions. This tutorial assumes you have followed the theta sequence
% tutorial, and that you have clean linearized positions.

%% Load spikes, positions, and theta cycles

basepath = 'Z:\Data\Can\OML22\day6'; % we will be using this session for the sake of this tutorial
basename = basenameFromBasepath(basepath);
spikeStructure = importSpikes('basepath',basepath,'CellType',"Pyramidal Cell",'brainRegion','CA1'); % load CA1 pyramidal cells
behavior = getStruct(basepath,'animal.behavior'); % load behavior structure

% We will be needing linearized positions:
positions = [behavior.timestamps(:) behavior.position.linearized(:)];
positions(any(isnan(positions),2),:) = [];
positions(:,2) = ZeroToOne(positions(:,2)); 
% Here the positions describe the position of the animal on a linear track.
% For the purposes of this tutorial, the direction of motion of the animal
% is ignored. Note that because of this, it would not make sense for us to 
% analyse forward / reverse slopes (because we will not be able to separate
% between reverse replay of the animal walking North and forward replay of
% the animal walking South). 
% For some projects, it may make more sense to concatenate the two directions, 
% so positions from 0 to 0.5 are when the animal is walking up the track and 
% from 0.5 to 1 - when the animal is going back. For others still, to score
% replay events for individual trajectories separately. In such cases, the 
% slope of the replay events would be meaningful, unlike the simple scoring 
% that this tutorial explores. Always keep in mind what makes sense for your 
% project. 

% We need running epochs, so make sure the behavior file contains them
if ~isfield(behavior,'run') % Here is some example code to detect running epochs:
    error('No "run" epochs in the behavior file.');
end

% We will need a spikes matrix:
try
    spikes = Group(spikeStructure.times{:},'sort',true);
catch 
    % if there are any problems executing "Group", this is how to produce the [timestamp id] matrix
    spikes = [];
    for i=1:length(spikeStructure.times)
        timestamps = spikeStructure.times{i};
        id = ones(size(timestamps))*i; % the id for that unit
        spikes = [spikes; timestamps id]; % add [timestamps id] for this unit to the matrix
    end
    spikes = sortrows(spikes); % sort spikes accoring to their time
end

%% Detect candidate intervals (bursts)

doRestrictToRipples = true; % this is a toggle that decides if you want your bursts to be restricted to detected ripples
if doRestrictToRipples
    try
        ripples = getStruct(basepath,'ripples');
    catch
        error(['No ripples found for session ' basepath]);
        % If you don't have ripples, you either need to detect them first
        % or set "doRestrictToRipples" to false and use all bursts
    end
end

bursts = FindBursts(spikes,'thresholds',[0 2],'smooth',0.01,'durations',[0 1]);
if doRestrictToRipples
    containsRipple = CountInIntervals(ripples.peaks(:),bursts(:,[1 3]))>0; % find burst that contains ripple in it
    bursts(~containsRipple,:) = []; % discard bursts that do not contain ripples
end
% Make sure none of the bursts are theta or gamma bursts associated with running:
inRunning = IntervalsIntersect(bursts(:,[1 3]),behavior.run); 
bursts(inRunning,:); 

% Now we have the final bursts, and this section will collect some information
% about them: the brain state (wake/sws) and the epoch (pre/task/post) that they took place

% Save the state for each burst (it makes sense to analyse awake replay and sleep replay separately):
load(fullfile(basepath,[basename '.SleepState.states.mat']),'SleepState');
awake = bsxfun(@plus,SleepState.ints.WAKEstate,[-1 0]);
sws = bsxfun(@plus,SleepState.ints.NREMstate,[-1 0]);
REM = bsxfun(@plus,SleepState.ints.REMstate,[-1 0]);
stateNames = {'WAKE','NREM','REM'}; % so state=1 means awake, state=2 means sws, state=3 means REM
state = InIntervals(bursts(:,2),awake)*1 + InIntervals(bursts(:,2),sws)*2 + InIntervals(bursts(:,2),REM)*3;

% Save the subsession in MergePoints for each burst:
MergePoints = getStruct(basepath,'MergePoints');
[~,epoch] = InIntervals(bursts(:,2),MergePoints.timestamps);

% Each of the burst intervals is a candidate replay event. To detect replay,
% we will split each burst in multiple 20ms bins and decode the position
% separately in each of those bins. The quality of the reconstruction will
% be tested for each candidate replay event later.
windowSize = 0.02; leeway = 0.001; 
% We will add 1ms to the start and the end of the burst. This is to make sure
% that the first spike is captured when we're reconstructing the position 
% within the interval. Sometimes, because of rounding errors, a spike at e.g.
% 5.3917 s may be considered outside of an interval starting at 5.3917 s. 
candidateIntervals = [bursts(:,1)-0.001 bursts(:,3)+0.001]; 
[rWindows,rID] = SplitIntervals(candidateIntervals,'pieceSize',windowSize);
% rWindows contains the 20ms bins for all the bursts, and rID contains 
% the ID of the burst that each "rWindows" bin belongs to. 

%% Compute Bayesian recostruction over the candidate events
nBins = 200; 
% As with theta cycles, we will use the "behavior.run" epoch to train the decoder
[estimations] = ReconstructPosition(positions,spikes,rWindows,'training',behavior.run,'nBins',nBins);
% While in theta cycles, we had the actual position of the animal to compare, 
% in the case of bursts, we will use the "estimations" probability matrix
empty = (max(estimations)==min(estimations))'; % keep in mind which bins are empty. Bursts with too many empty bins will not be scored

% Save each reconstructed event in a cell array. This makes it eacier to
% score many events in paralell later using parallel processing (parfor loop)
reconstructions = cell(size(bursts,1),1);
for i=1:size(bursts,1)
    reconstructions{i} = estimations(:,rID==i); 
end

%% Score candidate replay events

threshold = ceil(0.075*nBins); % As with theta sequences, use a threshold of 15 for the Davidson2009-style score
nShuffles = 2; % feel free to change to e.g. 500 to get p-values, but this will dramatically reduce the speed of this tutorial, so keep it 0 if you're just taking a look
[sequences,shuffled] = ReplayScore(reconstructions{1},'threshold',threshold,'nShuffles',nShuffles); % no need for shuffles as we don't require a p-value for each individual theta cycle
shuffleKind = 'column'; % In each shuffle, the function will apply a random shift in each of the columns and recompute the scores. 
% This is also referred to as "spatial shuffle" (because we are shufting the position coordinates in each column)
% Note that a temporal shuffle is also possible where we randomly shuffle the order of the columns instead.
for i=2:size(bursts,1) 
    if sum(rID==i & ~empty)>2 % minimum of 3 nonempty columns to score the event (otherwise the scores are meaningless)
        [sequence,shuffle] = ReplayScore(reconstructions{i},'threshold',threshold,'nShuffles',nShuffles,'shuffle',shuffleKind); 
    else
        [sequence,shuffle] = ReplayScore([],'threshold',threshold,'nShuffles',nShuffles,'shuffle',shuffleKind); % provide empty matrix for NaN scores
    end
    sequences = MergeStructures(sequences,sequence); % append the scores for each event
    shuffled = MergeStructures(shuffled,shuffle); % append the scores for each event
    % If you are impatient, you would be glad to see the progress of the scoring every 1000 events:
    if rem(i,1000)==0, disp([datestr(clock) ': ' num2str(i) '/' num2str(size(bursts,1)) ' done...']);  end 
end

%% Look at some example replay events

% Try to guess from "MergePoints" which is the task epoch. Feel free to change this for your project:
taskID = find(cellfun(@(x) ~contains(lower(x),'sleep') & ~contains(lower(x),'prob') & ~contains(lower(x),'baseline'),MergePoints.foldernames)); % the task epoch is the epoch that does not mention "sleep" or "prob" or "baseline"
% find the epoch most likely to correspond to pre-task sleep
sleepID = find(cellfun(@(x) any(strfind(lower(x),'sleep')),MergePoints.foldernames)); 
preID = sleepID(sleepID<taskID(1)); % any sleep epoch before the task is considered pre-task sleep
postID = sleepID(sleepID>taskID(1)); % any sleep epoch after the (first) task epoch is considered post-task sleep

pre = state==2 & ismember(epoch,preID); % all bursts within sws states in a "pre" epoch
task = state==1 & ismember(epoch,taskID); % all bursts within awake states in a "task" epoch
post = state==2 & ismember(epoch,postID); % all bursts within sws states in a "post" epoch

% Compute a single variable which tells us, for each event, if it is a sleep 
% event in pre-task sleep (1), an awake replay event during the task(2), or
% a sleep event in post-task sleep (3), or none of these (0);
group = pre*1 + task*2 + post*3; 

% an artibrary value that combines a couple of different measures: strong weighted correlation, good score, and low jump value
value = abs(sequences.weightedCorrelation) + sequences.score - sequences.jump.mean/nBins/2;

% Select examples in the top 20% of events (according to this arbitrary combined value we defined above). Also make sure event is not too short
examples = find(pre & (value > quantile(value(pre),0.8)) & diff(bursts(:,[1 3]),[],2)>0.08);
nExamples = 5;

clf
set(gcf,'position',[1 1 1920 964]);
for k=1:nExamples
    i = examples(k);
    subplot(3,nExamples,k);
    yPosition = linspace(0,1,size(reconstructions{i},1))'; x = ((1:size(reconstructions{i},2))-1)*windowSize + bursts(i,1);
    PlotColorMap(reconstructions{i},'y',yPosition,'x',x);
    hold on
    plot(x([1 end]),yPosition([sequences.lineStart(i) sequences.lineStop(i)]),'w--','linewidth',2);
    plot(x([1 end]),yPosition([sequences.lineStart(i) sequences.lineStop(i)])+threshold/nBins,'w');
    plot(x([1 end]),yPosition([sequences.lineStart(i) sequences.lineStop(i)])-threshold/nBins,'w');
    xlabel('Time (s)');
    clim([0 0.05])
    ylabel('Decoded position');
    set(gca,'box','off','TickDir','out','fontsize',12);
    title(['Event ' num2str(i) ' (pre); wcorr=' num2str(sequences.weightedCorrelation(i))]);
end

examples = find(task & (value > quantile(value(task),0.8)) & diff(bursts(:,[1 3]),[],2)>0.08);
for k=1:nExamples
    i = examples(k);
    subplot(3,nExamples,k+nExamples);
    yPosition = linspace(0,1,size(reconstructions{i},1))'; x = ((1:size(reconstructions{i},2))-1)*windowSize + bursts(i,1);
    PlotColorMap(reconstructions{i},'y',yPosition,'x',x);
    hold on
    plot(x([1 end]),yPosition([sequences.lineStart(i) sequences.lineStop(i)]),'w--','linewidth',2);
    plot(x([1 end]),yPosition([sequences.lineStart(i) sequences.lineStop(i)])+threshold/nBins,'w');
    plot(x([1 end]),yPosition([sequences.lineStart(i) sequences.lineStop(i)])-threshold/nBins,'w');
    xlabel('Time (s)');
    clim([0 0.05])
    ylabel('Decoded position');
    set(gca,'box','off','TickDir','out','fontsize',12);
    title(['Event ' num2str(i) ' (task); wcorr=' num2str(sequences.weightedCorrelation(i))]);
end

examples = find(post & (value > quantile(value(post),0.8)) & diff(bursts(:,[1 3]),[],2)>0.08);
for k=1:nExamples
    i = examples(k);
    subplot(3,nExamples,k+nExamples*2);
    yPosition = linspace(0,1,size(reconstructions{i},1))'; x = ((1:size(reconstructions{i},2))-1)*windowSize + bursts(i,1);
    PlotColorMap(reconstructions{i},'y',yPosition,'x',x);
    hold on
    plot(x([1 end]),yPosition([sequences.lineStart(i) sequences.lineStop(i)]),'w--','linewidth',2);
    plot(x([1 end]),yPosition([sequences.lineStart(i) sequences.lineStop(i)])+threshold/nBins,'w');
    plot(x([1 end]),yPosition([sequences.lineStart(i) sequences.lineStop(i)])-threshold/nBins,'w');
    xlabel('Time (s)');
    clim([0 0.05])
    ylabel('Decoded position');
    set(gca,'box','off','TickDir','out','fontsize',12);
    title(['Event ' num2str(i) ' (post); wcorr=' num2str(sequences.weightedCorrelation(i))]);
end

% The weighted correlation is plotted in the title of each panel
% The best bit line (used for computing the scores, slopes, and p-values) is superimposed.

%% Plot some replay statistics: compare the scores of pre-task sleep, task events, and post-task sleep

figure(1); clf;
set(gcf,'position',[600 600 900 300]);
% Compare the raw scores to the scores in shuffled data
values = [sequences.score shuffled.score];
labels = {'pre','task','post'};
for i=1:3,
    subplot(1,3,i);
    ok = group==i & ~any(isnan(values),2); % ignore other events (e.g. awake events in pre-task sleep)
    if size(values,2)<2
        handle = anovabox(values(ok,1:2),[],'alpha',[0 0.05],'paired',false);
        % We're setting the first alpha to 0 here because we don't want "anovabox"
        % to compare each of the groups to 0 (absolute values will by definition be
        % higher than zero).
        p = ranksum(values(group==i),values(group==i,2),'tail','right'); % compare if group "i" is higher than baseline
        PlotHVLines(nanmedian(values(group==i,2)),'h','k--','linewidth',2);
    else % if you have performed multiple shuffles, need to group the data differently to
        % perform the comparisons (rather than just having the shuffle be the 2nd column of "values")
        grouped = Group(values(ok,1),reshape(values(ok,2:end),[],1));
        handle = anovabox(grouped(:,1),grouped(:,2),'alpha',[0 0.05],'paired',false);
        p = ranksum(values(ok,1),reshape(values(ok,2:end),[],1),'tail','right'); % compare if group "i" is higher than baseline
        PlotHVLines(nanmedian(reshape(values(ok,2:end),[],1)),'h','k--','linewidth',2);
    end
    set(gca,'xtick',1:3,'xticklabel',{[labels{i} ', p=' num2str(p)],'shuffle'},'xticklabelrotation',45);
    ylabel('Score of linear fit');
    title(['Raw scores vs shuffles in ' labels{i}]);
end
EquateScales; % get the same xlim and ylim on all existing panels

figure(2); clf; 
set(gcf,'position',[600 200 900 300]);
% Compare quality of pre, task, and post replay, either normalizing by
% shuffled scores, or by the average score in baseline
subplot(1,2,1);
normalized = values;
for i=1:3, normalized(group==i,:) = (normalized(group==i,:) ./ nanmedian(reshape(values(group==i,2:end),[],1)) - 1)*100; end % measure as "% improvement" from shuffled scores in the respective condition
ok = group>0; % ignore other events (e.g. awake events in pre-task sleep)
anovabox(normalized(ok,1),group(ok),'alpha',[0 0.05]);
PlotHVLines(0,'h','k--','linewidth',2);
% p = signrank(sequences.quadrantScore); 
labels = {'pre','task','post'};
for k=1:3, p = ranksum(normalized(group==k,1),reshape(normalized(group==k,2:end),[],1),'tail','right'); % compare if group "k" is higher than baseline
    labels{k} = [labels{k} ', p=' num2str(p)]; 
end
set(gca,'xtick',1:3,'xticklabel',labels,'xticklabelrotation',45);
ylabel('% improvement from column shuffle');
title('Scores (relative to median score of respective shuffle)');

subplot(1,2,2);
normalized = (values./nanmedian(values(group==1))-1)*100; % measure as "% improvement" from scores in the pre-task sleep
ok = group>0; % ignore other events (e.g. awake events in pre-task sleep)
anovabox(normalized(ok),group(ok),'alpha',[0 0.05]);
PlotHVLines(0,'h','k--','linewidth',2);
% p = signrank(sequences.quadrantScore); 
labels = {'pre (baseline)','task','post'};
for k=2:3, p = ranksum(normalized(group==k),normalized(group==1),'tail','right'); % compare if group "k" is higher than baseline
    labels{k} = [labels{k} ', p=' num2str(p)]; 
end
set(gca,'xtick',1:3,'xticklabel',labels);
set(gca,'xtick',1:3,'xticklabel',labels,'xticklabelrotation',45);
ylabel('% improvement from baseline shuffle');
title('Scores (relative to median score of baseline)');

%% Perform the exact same tests, but using absolute value of the weighted correlation instead of scores

values = abs([sequences.weightedCorrelation shuffled.weightedCorrelation]);
figure(1); clf; 
set(gcf,'position',[600 600 900 300]);
labels = {'pre','task','post'};
for i=1:3,
    subplot(1,3,i);
    ok = group==i & ~any(isnan(values),2); % ignore other events (e.g. awake events in pre-task sleep)
    if size(values,2)<2
        handle = anovabox(values(ok,1:2),[],'alpha',[0 0.05],'paired',false);
        % We're setting the first alpha to 0 here because we don't want "anovabox"
        % to compare each of the groups to 0 (absolute values will by definition be
        % higher than zero).
        p = ranksum(values(group==i),values(group==i,2),'tail','right'); % compare if group "i" is higher than baseline
        PlotHVLines(nanmedian(values(group==i,2)),'h','k--','linewidth',2);
    else % if you have performed multiple shuffles, need to group the data differently to
        % perform the comparisons (rather than just having the shuffle be the 2nd column of "values")
        grouped = Group(values(ok,1),reshape(values(ok,2:end),[],1));
        handle = anovabox(grouped(:,1),grouped(:,2),'alpha',[0 0.05],'paired',false);
        p = ranksum(values(ok,1),reshape(values(ok,2:end),[],1),'tail','right'); % compare if group "i" is higher than baseline
        PlotHVLines(nanmedian(reshape(values(ok,2:end),[],1)),'h','k--','linewidth',2);
    end
    set(gca,'xtick',1:3,'xticklabel',{[labels{i} ', p=' num2str(p)],'shuffle'},'xticklabelrotation',45);
    ylabel('abs(Weighted correlation)');
    title('WCorr vs shuffles');
end

figure(2); clf; 
set(gcf,'position',[600 200 900 300]);
% Compare quality of pre, task, and post replay, either normalizing by
% shuffled values, or by the average score in baseline
subplot(1,2,1);
normalized = values;
for i=1:3, normalized(group==i,:) = (normalized(group==i,:) ./ nanmedian(reshape(values(group==i,2:end),[],1)) - 1)*100; end % measure as "% improvement" from shuffled scores in the respective condition
ok = group>0; % ignore other events (e.g. awake events in pre-task sleep)
anovabox(normalized(ok,1),group(ok),'alpha',[0 0.05]);
PlotHVLines(0,'h','k--','linewidth',2);
% p = signrank(sequences.quadrantScore); 
labels = {'pre','task','post'};
for k=1:3, p = ranksum(normalized(group==k,1),reshape(normalized(group==k,2:end),[],1),'tail','right'); % compare if group "k" is higher than baseline
    labels{k} = [labels{k} ', p=' num2str(p)]; 
end
set(gca,'xtick',1:3,'xticklabel',labels,'xticklabelrotation',45);
ylabel('% improvement from column shuffle');
title('Wcorrs (relative to median score of respective shuffle)');

subplot(1,2,2);
normalized = (values./nanmedian(values(group==1))-1)*100; % measure as "% improvement" from scores in the pre-task sleep
ok = group>0; % ignore other events (e.g. awake events in pre-task sleep)
anovabox(normalized(ok),group(ok),'alpha',[0 0.05]);
PlotHVLines(0,'h','k--','linewidth',2);
% p = signrank(sequences.quadrantScore); 
labels = {'pre (baseline)','task','post'};
for k=2:3, p = ranksum(normalized(group==k),normalized(group==1),'tail','right'); % compare if group "k" is higher than baseline
    labels{k} = [labels{k} ', p=' num2str(p)]; 
end
set(gca,'xtick',1:3,'xticklabel',labels,'xticklabelrotation',45);
ylabel('% improvement from baseline shuffle');
title('WCorr (relative to median score of baseline)');

% If we had used a larger number of shuffles (typical is n=500 shuffles), 
% we could be comparing % significant replay events using
% sequences.p, the p-values for each candidate event (not possible with
% n=1 shuffle), and sequences.zscores, which normalizes the score of
% each event with the mean score of all the shuffles of that event.


%% Compare the number of trajectory events as defined by Silva2015
% Silva2015 defined trajectory events as the events with weighted correlation 
% above 0.6 and jump distance below 0.4 (of the track):
trajectory = abs(sequences.weightedCorrelation)>0.6 & sequences.jump.max<0.4*nBins;
shTrajectory = abs(shuffled.weightedCorrelation)>0.6 & shuffled.jump.max<0.4*nBins;
nTrajectory = Accumulate(group(ok),trajectory(ok)); % How many trajectory events for each of the 3 groups
for i=1:size(shTrajectory,2) % for each shuffle
    nShTrajectory(:,i) = Accumulate(group(ok),shTrajectory(ok,i)); % How many trajectory events did the shuffles produce
end
nPossible = Accumulate(group(ok),~isnan(sequences.weightedCorrelation(ok))); % How many trajectory events are theoretically possible

figure(3)
set(gcf,'position',[600 300 900 600]);
clf 
subplot(1,2,1);
z = nan(3,1);
for k=1:3,
    z(k) = zBinomialComparison(nTrajectory(k),nPossible(k),mean(nShTrajectory(k,:)),nPossible(k));
    p = z2p(z(k)); % transform z-value into p-value
    labels{k} = [labels{k} ', p=' num2str(p)];
end

bar([nTrajectory mean(nShTrajectory,2)]);
set(gca,'xtick',1:3,'xticklabel',labels,'xticklabelrotation',45);
legend('data','shuffle');
legend('box','off','location','northwest');
set(gca,'fontsize',12,'box','off');
labels = {'pre (baseline)','task','post'};
ylabel('Number of trajectory events (wcorr>0.6, jump<0.4)');
title('Trajectory events');

subplot(1,2,2);
bar(z);
handle = PlotHVLines(p2z(0.05),'h','k--','linewidth',2);
legend(handle, 'p=0.05'); legend('box','off');
set(gca,'xtick',1:3,'xticklabel',labels,'xticklabelrotation',45);
ylabel('z-values, (s.d.-s above shuffle)');
set(gca,'fontsize',12,'box','off');
title('Effect size');

% This is a quick-and-dirty check for trajectory events: ideally you should
% perform multiple shuffles (Silva2015 performed 5000!) and actually count
% how many of them have more trajectory events than the observed data
% (whereas here we performed a binomial test based on a single shuffle).
% Note that results may be different with a temporal shuffle!

%% Save replay scores

% Before saving the replay sequence structure, I recommend you save some additional
% details which can be useful:
sequences.reconstructions = reconstructions;
sequences.intervals = bursts(:,[1 3]);
sequences.state = state; 
sequences.stateNames = {'WAKE','NREM','REM'};
sequences.epoch = epoch;
sequences.epochNames = MergePoints.foldernames';
sequences.events.timestamps = bursts(:,[1 3]);
sequences.events.peaks = bursts(:,2);
sequences.events.peakPower = bursts(:,4);
sequences.windowSize = windowSize;
sequences.params.windowSize = windowSize;
sequences.params.project = 'OML'; % feel free to label your own project
sequences.params.environment = 'LinearTrack';
sequences.training.intervals = behavior.run; % the training intervals that were used to decode the position
sequences.training.positions = positions; % the positions that were used to decode the position
sequences.params.windows = rWindows;
sequences.params.windowID = rID;

% Save sequences using the following names:
replay = sequences;
spatialShuffles = shuffled;
replayFilename = fullfile(basepath,[basename '.replay.mat']);

% The automatic pipeline script stops before the actual save to prevent
% the overwriting of any previously saved results with these limited 
% results (by default we only compute the results from a single shuffle)
return 
% Feel free to execute the following code in your own project
save(replayFilename,'replay','spatialShuffles');
% I like to perform a check if the file has saved properly (sometimes if variables are big it fails):
q = load(replayFilename,'replay','spatialShuffles');
if isempty(fieldnames(q)), % it didn't save properly, probably file size is too big
    save(replayFilename,'sequences','spatialShuffles','-v7.3'); % save using slower but more reliable protocol
end



