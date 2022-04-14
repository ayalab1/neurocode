blockSize = 5
randperm(numel(stimuli)*blockSize)
key = 2*ones(1:blockSize); key(1:end/2-1) = 1;
open key
key = 2*ones(1:blockSize,1);
key = 2*ones(blockSize,1); key(1:end/2-1) = 1;
key = 2*ones(blockSize,1);
key = 2*ones(blockSize,1); key(1:end/2) = 1;
key = 2*ones(blockSize*nStimuli,1); key(1:end/2) = 1;
stimuli = [0 90]; %stimulus angles: 8 stimuli; -1 is black
nStimuli = length(stimuli);
key = 2*ones(blockSize*nStimuli,1); key(1:end/2) = 1;
indices = randperm(numel(stimuli)*blockSize);
stimuli(randperm(numel(stimuli)))
key
theseAngles = stimuli(key(indices));
theseAngles = stimuli(key(indices))';
angles = [];
for j=1:nRepetitions/blockSize
key = 2*ones(blockSize*nStimuli,1); key(1:end/2) = 1;
indices = randperm(numel(stimuli)*blockSize);
theseAngles = stimuli(key(indices))';
angles = [angles; theseAngles(:)];
end
angles = repelem(angles,nFlashesPerStimulus+1); % add
nRepetitions/blockSize
blockSize
nRepetitions
nRepetitions = 30;
flashDuration = 1;
nFlashesPerStimulus = 2;
nProbeTrials = 5; % for each stimulus, this number of trials will be probe trials (presented sound, but no visual stimulus: this is to check if we get V1 representation of the missing visual cue associated with the presented sound!
firstPossibleProbeTrial = 15; % probe trials will be randomly distributed after the first e.g. 10 trials (it's good if training doesn't start with a probe trial)
data = nan(1+nRepetitions*nStimuli*(nFlashesPerStimulus+1),3);
% pre-set the pseudorandom stimulus order:
angles = [];
for j=1:nRepetitions/blockSize
key = 2*ones(blockSize*nStimuli,1); key(1:end/2) = 1;
indices = randperm(numel(stimuli)*blockSize);
theseAngles = stimuli(key(indices))';
angles = [angles; theseAngles(:)];
end
angles = repelem(angles,nFlashesPerStimulus+1); % add one extra line for the blank screen
ok = diff([-1;angles])~=0; angles(ok) = nan; % first line for each stimulus is actually the line for the blank screen
data(2:end,2) = angles;
task_audiovisual_pairing
%-- 2/2/2022 4:23 PM --%
v1
profile on
task_audiovisual_pairing
length(possibleLines)
possibleLines
stimuli(i)
data(2:end,2) = angles;
% Select random trials to be probe trials
if nProbeTrials>0
allPossibleLines = [1 strfind(data(:,2)'>=0,[0 1])+1 size(data,1)];
% write "-1" in 3rd column for probe trials:
ok = false;
while ~ok % start over if you got stuck in an infinite loop
data(:,3) = nan;
for i=1:length(stimuli)
if ~(stimuli(i)>-1) continue; end % there is no sound before the control stimulus anyway (or the blank screen (NaN))
possibleLines = strfind(data(:,2)'==stimuli(i),[0 1])+1; % the last line wasn't the stimulus (not a flash) but this line is
possibleLines(1:firstPossibleProbeTrial-1) = []; % first X trials are off limits and should not be probe trials
ok = false;
t = GetSecs;
while ~ok && GetSecs-t<2
indices = randperm(length(possibleLines),nProbeTrials);
lines = bsxfun(@plus,possibleLines(indices)',(0:nFlashesPerStimulus-1));
ok = true;
if any(diff(indices)==1), ok = false; end % avoid having 2 probe trials in a row (for the same stimulus)
if data(allPossibleLines(find(allPossibleLines==lines(1),1)-1),3)<0 || data(allPossibleLines(find(allPossibleLines==lines(1),1)+1),3)<0, ok = false; end % avoid having 2 probe trials in a row (even for different stimuli)
end
data(lines,3) = -1; % selected for a probe trial
end
end
end
i
stimuli(i)
possibleLines = strfind(data(:,2)'==stimuli(i),[0 1])+1;
possibleLines
clear all
task_audiovisual_pairing
length(possibleLines)
nProbeTrials
clear all
task_audiovisual_pairing
license
%-- 2/3/2022 2:23 PM --%
load('PeyracheShuffSpikes.mat')
raly
[assemblies,nMin] = SimSpikeAssemblies(spikes,0.02,3);
[td,triplets] = TolerantDefragmentSSA(origAssemblies,spikes,0.02,3,1,1,0.1);
[td,triplets] = TolerantDefragmentSSA(assemblies,spikes,0.02,3,1,1,0.1);
sum(td,2);
sum(assemblies,2);
[ntd,triplets] = DefragmentSSA(assemblies,spikes,0.02,3,1,1,0.1);
sum(ntd,2);
[td,triplets] = TolerantDefragmentSSA(assemblies,spikes,0.02,3,1,1,0.1);
combinations = cell2mat(cellCliques);
hits = cat(1,cellHits{:});
if tolerance>0
[~,u] = unique(combinations,'rows');
combinations = combinations(u,:);
hits = hits(u,:);
end
tolerance
combinations = cell2mat(cellCliques);
threshold
[~,u] = unique(combinations,'rows');
combinations = combinations(u,:);
hits = hits(u,:);
combinations = cell2mat(cellCliques);
hits = cat(1,cellHits{:});
[~,u] = unique(combinations,'rows');
combinations = combinations(u,:);
hits = hits(u,:);
n = floor(size(combinations,1)/1000)*1000;
pass = false(size(combinations));
threshold~=-Inf
verbose,tic;
loop = 00
n
n-1
n
display([datestr(clock) ': Last loop containing ' addComma(sum(cellfun(@(x) size(x,1),cellCliques((n+1):length(pass))))) ' cliques. toc=' addComma(round(toc))]);
n+1
length(pass)
qqqq = (n+1):length(pass);
open qqqq
length((n+1):size(combinations,1))
size(combinations,1)-n)
size(combinations,1)-n
display([datestr(clock) ': Last loop containing ' addComma(size(combinations,1)-n) ' cliques. toc=' addComma(round(toc))]);
j=(n+1)
parfor j=(n+1):size(combinations,1)% length(cellCliques)
thisCombination = combinations(j,:);
pass(j) = RateSSA2(thisCombination,spikes,windowSize,threshold,nMin);
if verbose && rem(j,100)==0, display([datestr(clock) ': cell passes ' num2str(j) ' computed.']); end
end
pass = false(size(combinations),1);
pass = false(size(combinations,1),1);
parfor j=(n+1):size(combinations,1)% length(cellCliques)
thisCombination = combinations(j,:);
pass(j) = RateSSA2(thisCombination,spikes,windowSize,threshold,nMin);
if verbose && rem(j,100)==0, display([datestr(clock) ': cell passes ' num2str(j) ' computed.']); end
end
sum(pass)
threshold
j=1
thisCombination = combinations(j,:);
sum(thisCombination)
find(sum(origAssemblies,2)>3,1)
find(origAssemblies(57,:))
find(sum(combinations(:,ans),2)==4)
j=23
thisCombination = combinations(j,:);
find(thisCombination)
RateSSA2(thisCombination,spikes,windowSize,threshold,nMin)
TRateSSA2(thisCombination,spikes,windowSize,threshold,nMin)
DRateSSA2(thisCombination,spikes,windowSize,threshold,nMin)
nonmemberMUA = spikes(~ok,1);
globalCount = sum(ExclusiveCountInIntervals(nonmemberMUA,activity))
globalCount = sum(ExclusiveCountInIntervals(spikes(:,1),activity)) - sum(ExclusiveCountInIntervals(s(:,1),activity))
count = sum(ExclusiveCountInIntervals(jSpikes,activity))
length(sum(~ok))
dbcont
[ntd,triplets] = DefragmentSSA(assemblies,spikes,0.02,3,1,1,0.1);
[td,triplets] = TolerantDefragmentSSA(assemblies,spikes,0.02,3,1,1,0.1);
dbcont
sum(td,2);
qqq = sum(td,2);
qq = sum(ntd,2);
Accumulate(td)
Accumulate(qqq)
Accumulate(qq)
[Accumulate(qq) Accumulate(qqq)]
edit mergeAssemblies.m
mergeAssemblies(td,ntd);
mergeAssemblies([td;ntd]);
isIncluded = bsxfun(@eq, commons,numNeurons);
PlotColorMap(isIncluded)
PlotColorMap(isIncluded+1)
close all
PlotColorMap(isIncluded+1)
open common
find(assemblies(3,:))
find(assemblies(4,:))
find(assemblies(2,:))
find(assemblies(15,:))
mergeAssemblies([td;ntd]');
isIncluded = bsxfun(@eq, commons,numNeurons);
mergeAssemblies([td;ntd]);
cats = [td;ntd];
cats = [td;ntd;groundTruthAssemblies'];
mergeAssemblies(cats');
assemblies = origAssemblies;
commons = double(assemblies)' * double(assemblies);
commons(1:length(commons)+1:end) = 0;
numNeurons = sum(assemblies);
isIncluded = bsxfun(@eq, commons,numNeurons);
PlotColorMap(isIncluded);
PlotColorMap(isIncluded+1);
plot(sum(isIncluded));
plot(sum(isIncluded(1:24,:)))
splitt = SSA_SplitAssemblies(td,3);
splitn = SSA_SplitAssemblies(ntd,3);
splitg = SSA_SplitAssemblies(groundTruthAssemblies',3);
cd D:
threshold = p2z(0.05)
threshold = p2z(0.01)
clear all
DefragmentSSA
4
num2str(max(sum(tolerant,2))
num2str(max(sum(tolerant,2)))
delete('temp-defragmented-03-Feb-2022-15h43-3.mat')
[fList,pList] = matlab.codetools.requiredFilesAndProducts('Tolerance_check.m');
cd D:
[fList,pList] = matlab.codetools.requiredFilesAndProducts('Tolerance_check.m');
open fList
fList = fList';
edit Portion
5
sum(assemblies,2);
tolerant = TolerantDefragmentSSA(assemblies,spikes,timescale,threshold,1,0,tolerance);
parfor j=(n+1):size(combinations,1)% length(cellCliques)
thisCombination = combinations(j,:);
pass(j) = RateSSA2(thisCombination,spikes,windowSize,threshold,nMin);
if verbose && rem(j,100)==0, display([datestr(clock) ': cell passes ' num2str(j) ' computed.']); end
end
n+1
size(combinations,1)
combinations = cell2mat(cellCliques);
j=1
possibleExtensions = [ j+find(sum(these(j+1:end,these(j,:)>0),2)==currentSize-1)];
find(these(1,:))
PlotColorMap(groungTruthAssemblies');
PlotColorMap(groundTruthAssemblies');
find(sum(groundTruthAssemblies([45 46 47 49],:))==4)
find(groundTruthAssemblies(10,:))
find(groundTruthAssemblies(10,:)>0)
groundTruthAssemblies(10,:);
find(groundTruthAssemblies(:,10))
[td,triplets] = TolerantDefragmentSSA(assemblies,spikes,0.02,3,1,1,0.1);
j=1
possibleExtensions = [ j+find(sum(these(j+1:end,these(j,:)>0),2)==currentSize-1)];
j
these(possibleExtensions,28)
possibleNeurons = sum(these([j;possibleExtensions],:),1)>=2
find(possibleNeurons)
possibleExtensions(any(these(possibleExtensions,~possibleNeurons),2)) = [];
possiblePartners = find(sum(these(:,possibleNeurons),2)==currentSize);
combinations = these(possibleExtensions,:); combinations(:,these(j,:)>0) = 1; combinations = unique(combinations,'rows');
find(combinations)(
find(combinations)
theseHits = cell(size(combinations,1),1); theseCliques = false(size(combinations,1),1);
i=1:size(combinations,1)
ok = combinations(i,:)>0;
comb = SSA_SplitAssemblies(combinations(i,:),3);
mean(~ismember(comb(:,ok),triplets(:,ok),'rows'))
PlotColorMap(comb(:,ok)+1);
subplot(1,2,1);
PlotColorMap(comb(:,ok)+1);
subplot(1,2,2);
PlotColorMap(triplets(:,ok)+1);
mean(~ismember(comb(:,ok),triplets(:,ok),'rows'))
~ismember(comb(:,ok),triplets(:,ok),'rows');
find(comb(end,ok))
triplets = SSA_SplitAssemblies(newAssemblies,3);
[ntd,triplets] = DefragmentSSA(assemblies,spikes,0.02,3,1,1,0.1);
sum(triplets(1,:))
2
dbcont
2
display(['Tolerant defragmentation found ' num2str(Portion(ismember(splitg,splitt,'rows'))*100) '% of the total triplets in the ground truth data']);
display(['False positive rate was ' num2str(1-Portion(ismember(splitt,splitg,'rows'))*100) '% of the found triplets']);
display(['False positive rate was ' num2str(100-Portion(ismember(splitt,splitg,'rows'))*100) '% of the found triplets']);
2
nchoosek(1:5,3);
tolerant0 = tolerant
tolerant0 = tolerant;
tolerant1 = tolerant;
tic; tolerant = TolerantDefragmentSSA(assemblies,spikes,timescale,threshold,1,0,tolerance); toc
tic; for i=1:100, tolerant = TolerantDefragmentSSA(assemblies,spikes,timescale,threshold,1,0,tolerance); end; toc
tic; for i=1:10, tolerant = TolerantDefragmentSSA(assemblies,spikes,timescale,threshold,1,0,tolerance); end; toc
gt = groundTruthAssemblies';
gtTriplets = SSA_SplitAssemblies(gt,3);
nMissed = round((1-sensitivity)*size(gtTriplets,1))
sensitivity = 0.9;
nMissed = round((1-sensitivity)*size(gtTriplets,1))
nTriplets = size(gtTriplets,1);
order = randperm(199,nMissed);
idxMissed = sort(randperm(199,nMissed));
defrag = TolerantDefragmentSSA(assemblies,[],timescale,-Inf,1,0,tolerance); toc
defrag = TolerantDefragmentSSA(detected,[],timescale,-Inf,1,0,tolerance); toc
tolerance
fp(i,1) = Portion(1-ismember(splitt,splitg,'rows'))
clear fp
i=1
fp(i,1) = Portion(1-ismember(splitt,splitg,'rows'))
clear tp
i
tp(i,1) = Portion(ismember(splitg,splitt,'rows'))
54
j
i
45
j
semplot(tolerances,tp,'k');
j
close all
semplot(tolerances,tp,'k');
semplot(tolerances,fp,'r');
clf
plot(tr,fp,'.');
plot(fp,tp,'.');
plot(fp,tp,'-.');
plot(fp(~isnan(fp)),tp(~isnan(fp)),'-.');
plot(fp(~isnan(fp)),tp(~isnan(fp)),'.-');
clf
semplot(tolerances,fp,'r');
semplot(tolerances,tp,'k');
open tp
semplot(tolerances,tp,'k');
clf
semplot(tolerances,tp,'k');
semplot(tolerances,fp,'r');
j=5;
j=6
i=6
tolerance = tolerances(j);
rng(i)
idxDetected = sort(randperm(199,nDetected));
detected = gtTriplets(idxDetected,:);
tolerant = TolerantDefragmentSSA(detected,[],timescale,-Inf,1,0,tolerance);
sum(tolerant,2);
open ans
find(~ismember(splitg,splitt,'rows'))
splitt = SSA_SplitAssemblies(tolerant,3);
find(~ismember(splitg,splitt,'rows'))
find(splitg(98,:))
find(splitg(99,:))
find(splitg(163,:))
find(splitg(169,:))
find(sum(gt(:,[4 17 24 32 36]),2)>4)
find(gt(13,:))
sum(detected(:,gt(13,:)),2);
sum(detected(:,gt(13,:)>0),2);
Accumulate(sum(detected(:,gt(13,:)>0),2)+1)
find(sum(detected(:,gt(13,:)>0),2)>1)
PlotColorMap(detected(find(sum(detected(:,gt(13,:)>0),2)>1),:));
clf
PlotColorMap(detected(find(sum(detected(:,gt(13,:)>0),2)>1),:));
hm =SSA_SplitAssemblies(gt(13,:),3);
size(hm)
Portion(ismember(hm,detected,'rows'));
Portion(ismember(hm,detected,'rows'))
tolerance
hmd = detected(find(sum(detected(:,gt(13,:)>0),2)>1),:);
open hmd
find(gt(13,:))
PlotHVLines(ans,'v');
size(hm)
open hmd
hmd = detected(ismember(detected,hm,'rows'),:);
tolerant = TolerantDefragmentSSA(hmd,[],timescale,-Inf,1,0,tolerance);
open tolerant
find(tolerant)
hmd = detected(find(sum(detected(:,gt(13,:)>0),2)>1),:);
tolerant = TolerantDefragmentSSA(hmd,[],timescale,-Inf,1,0,tolerance);
open tolerant
sum(tolerant,2)
sum(tolerant(:,[4 17 24 32 36]),2)
tolerant = TolerantDefragmentSSA(detected,[],timescale,-Inf,1,0,tolerance);
findmax(sum(tolerant(:,[4 17 24 32 36]),2))
max(sum(tolerant(:,[4 17 24 32 36]),2))
tolerant = TolerantDefragmentSSA(detected(find(sum(detected(:,gt(13,:)>0),2)>1),:),[],timescale,-Inf,1,0,tolerance);
max(sum(tolerant(:,[4 17 24 32 36]),2))
this = sortby(detected'
this = sortby(detected',gt(13,:))';
open this
PlotColorMap(this);
this = sortby(gt',-gt(13,:))';
PlotColorMap(this);
this = sortby(detected',-gt(13,:))';
PlotColorMap(this);
this = sortby(gt',-gt(13,:))';
PlotColorMap(this);
max(sum(this(:,1:5),2))
detected = gtTriplets(idxDetected,:);
this = sortby(gt',-gt(13,:))';
max(sum(this(:,1:5),2))
this = sortby(detected',-gt(13,:))';
max(sum(this(:,1:5),2))
sum(sum(this(:,1:5),2))==3)
sum(sum(this(:,1:5),2)==3)
PlotColorMap(this);
size(detected)
size(this)
clf
PlotColorMap(this+1);
this = sortrows(sortby(detected',-gt(13,:))');
PlotColorMap(this+1);
this = (sortby(detected',-gt(13,:))'); this = sortby(this,sum(this(:,1:5),2));
PlotColorMap(this+1);
this = (sortby(detected',-gt(13,:))'); this = sortby(this,this(:,1:5));
PlotColorMap(this+1);
this = (sortby(detected',-gt(13,:))'); this = sortby(this,[sum(this(:,1:5),2)*1000 + this(:,1)*100 + this(:,2)*10 + this(:,3)*1 + this(:,4)*0.1 + this(:,5)*0.001]);
PlotColorMap(this+1);
this = (sortby(detected',-gt(13,:))'); this = sortby(this,[sum((this(:,1:5),2)>1)*1000 + this(:,1)*100 + this(:,2)*10 + this(:,3)*1 + this(:,4)*0.1 + this(:,5)*0.001]);
this = (sortby(detected',-gt(13,:))'); this = sortby(this,[(sum(this(:,1:5),2))>1)*1000 + this(:,1)*100 + this(:,2)*10 + this(:,3)*1 + this(:,4)*0.1 + this(:,5)*0.001]);
this = (sortby(detected',-gt(13,:))'); this = sortby(this,[(sum(this(:,1:5),2)>1)*1000 + this(:,1)*100 + this(:,2)*10 + this(:,3)*1 + this(:,4)*0.1 + this(:,5)*0.001]);
PlotColorMap(this+1);
this = (sortby(detected',-gt(13,:))'); this = sortby(this,[((sum(this(:,1:5),2)>1)-0.5)*(this(:,1)*100 + this(:,2)*10 + this(:,3)*1 + this(:,4)*0.1 + this(:,5)*0.001)]);
this = (sortby(detected',-gt(13,:))'); this = sortby(this,[((sum(this(:,1:5),2)>1)-0.5).*(this(:,1)*100 + this(:,2)*10 + this(:,3)*1 + this(:,4)*0.1 + this(:,5)*0.001)]);
PlotColorMap(this+1);
this = (sortby(detected',-gt(13,:))'); this = sortby(this,[sum(this(:,1:5),2)*1000 + this(:,1)*100 + this(:,2)*10 + this(:,3)*1 + this(:,4)*0.1 + this(:,5)*0.001]);
PlotColorMap(this+1);
sum(sum(this(:,1:5),2)>1)
this = (sortby(detected',-gt(13,:))'); this = sortby(this,-[sum(this(:,1:5),2)*1000 + this(:,1)*100 + this(:,2)*10 + this(:,3)*1 + this(:,4)*0.1 + this(:,5)*0.001]);
PlotColorMap(this+1);
j
i
title(max(sum(tolerant(:,gt(13,:)>0),2)))
34
345
this = (sortby(detected',-gt(13,:))');
find(this(1,:))
this = Scramble(this);
find(this(1,:))
PlotColorMap(sortrows(this)+1);
this = Scramble(this);
PlotColorMap(sortrows(this)+1);
for i=1:100,
rng(i);
this = (sortby(detected',-gt(13,:))');
this = Scramble(this);
% this = this(1:179,:);
tolerant = TolerantDefragmentSSA(this,[],timescale,-Inf,1,0,tolerance);
% title(max(sum(tolerant(:,1:5),2)));
qq(i,1) = (max(sum(tolerant(:,1:5),2)));
end
this = (sortby(detected',gt(13,:))');
%     this = Scramble(this);
% this = this(1:179,:);
tolerant = TolerantDefragmentSSA(this,[],timescale,-Inf,1,0,tolerance);
(max(sum(tolerant(:,1:5),2)))
this = (sortby(detected',gt(13,:))');
this = sortby(this,-[sum(this(:,1:5),2)*1000 + this(:,1)*100 + this(:,2)*10 + this(:,3)*1 + this(:,4)*0.1 + this(:,5)*0.001]);
tolerant = TolerantDefragmentSSA(this,[],timescale,-Inf,1,0,tolerance);
(max(sum(tolerant(:,1:5),2)))
tolerant = TolerantDefragmentSSA(this,[],timescale,-Inf,1,0,tolerance);
qq(i,1) = (max(sum(tolerant(:,end-4:end),2)));
(max(sum(tolerant(:,end-4:end),2)))
5
45
find(seek)
keep = nan(100,6);
open keep
i
clf
anovabar(keep(:,1),keep(:,6)-keep(:,2));
that = sortby(keep,keep(:,1));
open that
PlotColorMap(that);
that = sortby(keep,keep(:,1)); that(:,1) = that(:,1)*5;
PlotColorMap(that);
that = sortby(keep,keep(:,1)); that(:,1) = that(:,1)*10;
PlotColorMap(that);
corr(that(:,1),that(:,2:end))
anovabar(that(:,1),that(:,2))
PlotXY(that(:,[2 1]),'.');
PlotXY(that(:,[2 1])+randn(size(that))/100,'.');
PlotXY(that(:,[2 1])+randn(size(that(:,[1 1])))/100,'.');
PlotXY(that(:,[2 1])+randn(size(that(:,[1 1])))/10,'.');
[~,orderUnits] = sort(-gt(13,:));
orderUnits = circshift(orderUnits(:),1)';
open keep2
i=1
[~,orderUnits] = sort(-gt(13,:));
orderUnits = circshift(orderUnits(:),i)';
detectedOrdered = detected(:,orderUnits);
seek = gt(13,orderUnits)>0;
find(seek)
orderedAssemblies = detectedOrdered;
tolerant = TolerantDefragmentSSA(orderedAssemblies,[],timescale,-Inf,1,0,tolerance);
(max(sum(tolerant(:,seek),2)))
tolerant = TolerantDefragmentSSA(orderedAssemblies(sum(orderedAssemblies(:,seek),2)>1,:),[],timescale,-Inf,1,0,tolerance);
(max(sum(tolerant(:,seek),2)))
tolerant = TolerantDefragmentSSA(orderedAssemblies(sum(orderedAssemblies(:,seek),2)>0,:),[],timescale,-Inf,1,0,tolerance);
(max(sum(tolerant(:,seek),2)))
tolerant = TolerantDefragmentSSA(orderedAssemblies(sum(orderedAssemblies(:,seek),2)>-1,:),[],timescale,-Inf,1,0,tolerance);
(max(sum(tolerant(:,seek),2)))
tolerant = TolerantDefragmentSSA(orderedAssemblies,[],timescale,-Inf,1,0,tolerance);
(max(sum(tolerant(:,seek),2)))
find(seek)
% the assembly is [2 3 4 5 6]
tolerant = TolerantDefragmentSSA(orderedAssemblies,[],timescale,-Inf,1,0,tolerance);
sum(these(:,2:6),2);
open ans
for j=1:size(these,1) % for assembly j
possibleExtensions = [ j+find(sum(these(j+1:end,these(j,:)>0),2)==currentSize-1)];
% remove neurons appearing too few times to possibly participate in a clique
possibleNeurons = sum(these([j;possibleExtensions],:),1)>=1;
possibleExtensions(any(these(possibleExtensions,~possibleNeurons),2)) = [];
possiblePartners = find(sum(these(:,possibleNeurons),2)==currentSize);
combinations = these(possibleExtensions,:); combinations(:,these(j,:)>0) = 1; combinations = unique(combinations,'rows');
theseHits = cell(size(combinations,1),1); theseCliques = false(size(combinations,1),1);
for i=1:size(combinations,1)
ok = combinations(i,:)>0;
comb = SSA_SplitAssemblies(combinations(i,:),3);
if mean(~ismember(comb(:,ok),triplets(:,ok),'rows'))<=tolerance
theseCliques(i) = true;
hitID = sum(these(possibleExtensions,ok),2)==currentSize;
theseHits{i} = possiblePartners(hitID);
end
end
cellCliques{j} = combinations(theseCliques,:); cellHits{j} = theseHits(theseCliques,:);
if verbose && (rem(j,500)==0)
display([datestr(clock) ': So far found ' addComma(sum(cellfun(@(x) size(x,1),cellCliques))) ' cliques after going through ' ...
addComma(j) '/' addComma(size(these,1)) ' groups. toc=' addComma(round(toc))]);
end
end
sum(these(:,2:6),2);
sums = sum(these(:,2:6),2);;
sums26 = sums; cliques26 = cellCliques;
save('temp26.mat','sums26','cliques26');
[~,orderUnits] = sort(-gt(13,:));
detectedOrdered = detected(:,orderUnits);
seek = gt(13,orderUnits)>0;
orderedAssemblies = detectedOrdered;
find(seek)
tolerant = TolerantDefragmentSSA(orderedAssemblies,[],timescale,-Inf,1,0,tolerance);
sums15 = sum(these(:,1:5),2);;
for j=1:size(these,1) % for assembly j
possibleExtensions = [ j+find(sum(these(j+1:end,these(j,:)>0),2)==currentSize-1)];
% remove neurons appearing too few times to possibly participate in a clique
possibleNeurons = sum(these([j;possibleExtensions],:),1)>=1;
possibleExtensions(any(these(possibleExtensions,~possibleNeurons),2)) = [];
possiblePartners = find(sum(these(:,possibleNeurons),2)==currentSize);
combinations = these(possibleExtensions,:); combinations(:,these(j,:)>0) = 1; combinations = unique(combinations,'rows');
theseHits = cell(size(combinations,1),1); theseCliques = false(size(combinations,1),1);
for i=1:size(combinations,1)
ok = combinations(i,:)>0;
comb = SSA_SplitAssemblies(combinations(i,:),3);
if mean(~ismember(comb(:,ok),triplets(:,ok),'rows'))<=tolerance
theseCliques(i) = true;
hitID = sum(these(possibleExtensions,ok),2)==currentSize;
theseHits{i} = possiblePartners(hitID);
end
end
cellCliques{j} = combinations(theseCliques,:); cellHits{j} = theseHits(theseCliques,:);
if verbose && (rem(j,500)==0)
display([datestr(clock) ': So far found ' addComma(sum(cellfun(@(x) size(x,1),cellCliques))) ' cliques after going through ' ...
addComma(j) '/' addComma(size(these,1)) ' groups. toc=' addComma(round(toc))]);
end
end
cliques15 = cellCliques;
save('temp15.mat','sums15','cliques15');
load('temp15.mat')
load('temp26.mat')
cliques = [cliques15 cliques26];
open cliques
size(cell2mat(cliques(:,1)))
size(cell2mat(cliques(:,2)))
size(unique(cell2mat(cliques(:,2)),'rows');
size(unique(cell2mat(cliques(:,2)),'rows'));
size(unique(cell2mat(cliques(:,2)),'rows'))
size(unique(cell2mat(cliques(:,1)),'rows'))
open sums26
sum(sums26==sums15)
sum(sums26~=sums15)
sums = [sums26 sums15];
open sums
tolerant = TolerantDefragmentSSA(orderedAssemblies,[],timescale,-Inf,1,0,tolerance);
these15 = these;
sums15 = sum(these(:,1:5),2);;
cliques15 = cellCliques;
load('temp15.mat')
save('temp15.mat','sums15','cliques15','these15');
[~,orderUnits] = sort(-gt(13,:));
i=1
orderUnits = circshift(orderUnits(:),i)';
% [~,orderUnits] = sort(1:50);
detectedOrdered = detected(:,orderUnits);
seek = gt(13,orderUnits)>0;
orderedAssemblies = detectedOrdered;
find(seek)
tolerant = TolerantDefragmentSSA(orderedAssemblies,[],timescale,-Inf,1,0,tolerance);
load('temp26.mat')
these26 = these;
save('temp26.mat','sums26','cliques26','these26');
load('temp15.mat')
load('temp26.mat')
size(these26)
size(these15)
find(these15(1,:))
find(these26(1,:))
find(these26(end,:))
find(these15(end,:))
tt = these15;
tt2 = these26(:,[2:end 1]);
sum(tt(1,:)==tt2(1,:))
sum(tt(end,:)==tt2(end,:))
PlotColorMap(tt-tt2)
sum(ismember(tt,tt2,'rows')
sum(ismember(tt,tt2,'rows'))
sum(~ismember(tt,tt2,'rows'))
[o1,o2,o3] = ismember(tt,tt2,'rows');
[o1,o2] = ismember(tt,tt2,'rows');
open o2
find(tt(164,:))
find(tt2(8,:))
find(tt2(164,:))
find(tt(8,:))
tt3 = tt2(o2,:);
find(tt3(8,:))
cliques = [cliques15 cliques26(order,:)];
cliques = [cliques15 cliques26(o2,:)];
ttt = [these15 these26(o2,:)];
sums = [sums15 sums26(o2,:)];
open sums
open cliques
open sums
find(cliques{178,1})
find(cliques{178,2})
find(cliques{178,2}(1,:))
find(cliques{178,1}(1,:))
find(cliques{178,2}(1,:))-1
find(cliques{178,2}(2,:))-1
pertinantCliques = cliques(sums(:,1)>1,:);
open pertinantCliques
find(pertinantCliques{7,1}(1,:))
find(pertinantCliques{7,1}(2,:))
find(pertinantCliques{7,1}(3,:))
temp = cell2mat(pertinantCliques);
size(cell2mat(pertinantCliques(:,1)))
size(cell2mat(pertinantCliques(:,2)))
temp1 = cell2mat(pertinantCliques(:,1));
temp2 = cell2mat(pertinantCliques(:,2));
open temp2
open temp1
temp2 = cell2mat(pertinantCliques(:,2)); temp2 = temp2(:,[2:end 1]);
sum(~ismember(temp1,temp2,'rows'))
sum(~ismember(temp2,temp1,'rows'))
find(~ismember(temp2,temp1,'rows'))
find(temp2(6,:))
open tt1
open tt
open cliques
open sums
these = these15;
possibleExtensions = [ j+find(sum(these(j+1:end,these(j,:)>0),2)==currentSize-1)];
% remove neurons appearing too few times to possibly participate in a clique
possibleNeurons = sum(these([j;possibleExtensions],:),1)>=1;
possibleExtensions(any(these(possibleExtensions,~possibleNeurons),2)) = [];
possiblePartners = find(sum(these(:,possibleNeurons),2)==currentSize);
combinations = these(possibleExtensions,:); combinations(:,these(j,:)>0) = 1; combinations = unique(combinations,'rows');
theseHits = cell(size(combinations,1),1); theseCliques = false(size(combinations,1),1);
for i=1:size(combinations,1)
ok = combinations(i,:)>0;
comb = SSA_SplitAssemblies(combinations(i,:),3);
if mean(~ismember(comb(:,ok),triplets(:,ok),'rows'))<=tolerance
theseCliques(i) = true;
hitID = sum(these(possibleExtensions,ok),2)==currentSize;
theseHits{i} = possiblePartners(hitID);
end
end
currentSize = 4;
possibleExtensions = [ j+find(sum(these(j+1:end,these(j,:)>0),2)==currentSize-1)];
% remove neurons appearing too few times to possibly participate in a clique
possibleNeurons = sum(these([j;possibleExtensions],:),1)>=1;
possibleExtensions(any(these(possibleExtensions,~possibleNeurons),2)) = [];
possiblePartners = find(sum(these(:,possibleNeurons),2)==currentSize);
combinations = these(possibleExtensions,:); combinations(:,these(j,:)>0) = 1; combinations = unique(combinations,'rows');
theseHits = cell(size(combinations,1),1); theseCliques = false(size(combinations,1),1);
for i=1:size(combinations,1)
ok = combinations(i,:)>0;
comb = SSA_SplitAssemblies(combinations(i,:),3);
if mean(~ismember(comb(:,ok),triplets(:,ok),'rows'))<=tolerance
theseCliques(i) = true;
hitID = sum(these(possibleExtensions,ok),2)==currentSize;
theseHits{i} = possiblePartners(hitID);
end
end
qq1 = combinations(theseCliques,:);
open qq1
j
these = these15;
j=1;
possibleExtensions = [ j+find(sum(these(j+1:end,these(j,:)>0),2)==currentSize-1)];
% remove neurons appearing too few times to possibly participate in a clique
possibleNeurons = sum(these([j;possibleExtensions],:),1)>=1;
possibleExtensions(any(these(possibleExtensions,~possibleNeurons),2)) = [];
possiblePartners = find(sum(these(:,possibleNeurons),2)==currentSize);
combinations = these(possibleExtensions,:); combinations(:,these(j,:)>0) = 1; combinations = unique(combinations,'rows');
theseHits = cell(size(combinations,1),1); theseCliques = false(size(combinations,1),1);
for i=1:size(combinations,1)
ok = combinations(i,:)>0;
comb = SSA_SplitAssemblies(combinations(i,:),3);
if mean(~ismember(comb(:,ok),triplets(:,ok),'rows'))<=tolerance
theseCliques(i) = true;
hitID = sum(these(possibleExtensions,ok),2)==currentSize;
theseHits{i} = possiblePartners(hitID);
end
end
cellCliques{j,1} = combinations(theseCliques,:);
possibleExtensions = [ j+find(sum(these(j+1:end,these(j,:)>0),2)==currentSize-1)];
max(sum(these(j+1:end,these(j,:)>0),2))
currentSize = 3;
possibleExtensions = [ j+find(sum(these(j+1:end,these(j,:)>0),2)==currentSize-1)];
% remove neurons appearing too few times to possibly participate in a clique
possibleNeurons = sum(these([j;possibleExtensions],:),1)>=1;
possibleExtensions(any(these(possibleExtensions,~possibleNeurons),2)) = [];
possiblePartners = find(sum(these(:,possibleNeurons),2)==currentSize);
combinations = these(possibleExtensions,:); combinations(:,these(j,:)>0) = 1; combinations = unique(combinations,'rows');
theseHits = cell(size(combinations,1),1); theseCliques = false(size(combinations,1),1);
for i=1:size(combinations,1)
ok = combinations(i,:)>0;
comb = SSA_SplitAssemblies(combinations(i,:),3);
if mean(~ismember(comb(:,ok),triplets(:,ok),'rows'))<=tolerance
theseCliques(i) = true;
hitID = sum(these(possibleExtensions,ok),2)==currentSize;
theseHits{i} = possiblePartners(hitID);
end
end
cellCliques{j} = combinations(theseCliques,:); cellHits{j} = theseHits(theseCliques,:);
triplets = SSA_SplitAssemblies(tt,3);
possibleExtensions = [ j+find(sum(these(j+1:end,these(j,:)>0),2)==currentSize-1)];
% remove neurons appearing too few times to possibly participate in a clique
possibleNeurons = sum(these([j;possibleExtensions],:),1)>=1;
possibleExtensions(any(these(possibleExtensions,~possibleNeurons),2)) = [];
possiblePartners = find(sum(these(:,possibleNeurons),2)==currentSize);
combinations = these(possibleExtensions,:); combinations(:,these(j,:)>0) = 1; combinations = unique(combinations,'rows');
theseHits = cell(size(combinations,1),1); theseCliques = false(size(combinations,1),1);
for i=1:size(combinations,1)
ok = combinations(i,:)>0;
comb = SSA_SplitAssemblies(combinations(i,:),3);
if mean(~ismember(comb(:,ok),triplets(:,ok),'rows'))<=tolerance
theseCliques(i) = true;
hitID = sum(these(possibleExtensions,ok),2)==currentSize;
theseHits{i} = possiblePartners(hitID);
end
end
cellCliques{j} = combinations(theseCliques,:); cellHits{j} = theseHits(theseCliques,:);
j=178l
j=178;
possibleExtensions = [ j+find(sum(these(j+1:end,these(j,:)>0),2)==currentSize-1)];
% remove neurons appearing too few times to possibly participate in a clique
possibleNeurons = sum(these([j;possibleExtensions],:),1)>=1;
possibleExtensions(any(these(possibleExtensions,~possibleNeurons),2)) = [];
possiblePartners = find(sum(these(:,possibleNeurons),2)==currentSize);
combinations = these(possibleExtensions,:); combinations(:,these(j,:)>0) = 1; combinations = unique(combinations,'rows');
theseHits = cell(size(combinations,1),1); theseCliques = false(size(combinations,1),1);
for i=1:size(combinations,1)
ok = combinations(i,:)>0;
comb = SSA_SplitAssemblies(combinations(i,:),3);
if mean(~ismember(comb(:,ok),triplets(:,ok),'rows'))<=tolerance
theseCliques(i) = true;
hitID = sum(these(possibleExtensions,ok),2)==currentSize;
theseHits{i} = possiblePartners(hitID);
end
end
cellCliques{j} = combinations(theseCliques,:); cellHits{j} = theseHits(theseCliques,:);
cellCliques = cellCliques';
these = these26;
possibleExtensions = [ j+find(sum(these(j+1:end,these(j,:)>0),2)==currentSize-1)];
% remove neurons appearing too few times to possibly participate in a clique
possibleNeurons = sum(these([j;possibleExtensions],:),1)>=1;
possibleExtensions(any(these(possibleExtensions,~possibleNeurons),2)) = [];
possiblePartners = find(sum(these(:,possibleNeurons),2)==currentSize);
combinations = these(possibleExtensions,:); combinations(:,these(j,:)>0) = 1; combinations = unique(combinations,'rows');
theseHits = cell(size(combinations,1),1); theseCliques = false(size(combinations,1),1);
for i=1:size(combinations,1)
ok = combinations(i,:)>0;
comb = SSA_SplitAssemblies(combinations(i,:),3);
if mean(~ismember(comb(:,ok),triplets(:,ok),'rows'))<=tolerance
theseCliques(i) = true;
hitID = sum(these(possibleExtensions,ok),2)==currentSize;
theseHits{i} = possiblePartners(hitID);
end
end
cellCliques{j} = combinations(theseCliques,:); cellHits{j} = theseHits(theseCliques,:);
these = tt3;
possibleExtensions = [ j+find(sum(these(j+1:end,these(j,:)>0),2)==currentSize-1)];
% remove neurons appearing too few times to possibly participate in a clique
possibleNeurons = sum(these([j;possibleExtensions],:),1)>=1;
possibleExtensions(any(these(possibleExtensions,~possibleNeurons),2)) = [];
possiblePartners = find(sum(these(:,possibleNeurons),2)==currentSize);
combinations = these(possibleExtensions,:); combinations(:,these(j,:)>0) = 1; combinations = unique(combinations,'rows');
theseHits = cell(size(combinations,1),1); theseCliques = false(size(combinations,1),1);
for i=1:size(combinations,1)
ok = combinations(i,:)>0;
comb = SSA_SplitAssemblies(combinations(i,:),3);
if mean(~ismember(comb(:,ok),triplets(:,ok),'rows'))<=tolerance
theseCliques(i) = true;
hitID = sum(these(possibleExtensions,ok),2)==currentSize;
theseHits{i} = possiblePartners(hitID);
end
end
cellCliques{j} = combinations(theseCliques,:); cellHits{j} = theseHits(theseCliques,:);
these = tt2;
possibleExtensions = [ j+find(sum(these(j+1:end,these(j,:)>0),2)==currentSize-1)];
% remove neurons appearing too few times to possibly participate in a clique
possibleNeurons = sum(these([j;possibleExtensions],:),1)>=1;
possibleExtensions(any(these(possibleExtensions,~possibleNeurons),2)) = [];
possiblePartners = find(sum(these(:,possibleNeurons),2)==currentSize);
combinations = these(possibleExtensions,:); combinations(:,these(j,:)>0) = 1; combinations = unique(combinations,'rows');
theseHits = cell(size(combinations,1),1); theseCliques = false(size(combinations,1),1);
for i=1:size(combinations,1)
ok = combinations(i,:)>0;
comb = SSA_SplitAssemblies(combinations(i,:),3);
if mean(~ismember(comb(:,ok),triplets(:,ok),'rows'))<=tolerance
theseCliques(i) = true;
hitID = sum(these(possibleExtensions,ok),2)==currentSize;
theseHits{i} = possiblePartners(hitID);
end
end
cellCliques{j} = combinations(theseCliques,:); cellHits{j} = theseHits(theseCliques,:);
2
open o2
j
these = tt2;
possibleExtensions = [ j+find(sum(these(j+1:end,these(j,:)>0),2)==currentSize-1)];
% remove neurons appearing too few times to possibly participate in a clique
possibleNeurons = sum(these([j;possibleExtensions],:),1)>=1;
possibleExtensions(any(these(possibleExtensions,~possibleNeurons),2)) = [];
possiblePartners = find(sum(these(:,possibleNeurons),2)==currentSize);
combinations = these(possibleExtensions,:); combinations(:,these(j,:)>0) = 1; combinations = unique(combinations,'rows');
theseHits = cell(size(combinations,1),1); theseCliques = false(size(combinations,1),1);
for i=1:size(combinations,1)
ok = combinations(i,:)>0;
comb = SSA_SplitAssemblies(combinations(i,:),3);
if mean(~ismember(comb(:,ok),triplets(:,ok),'rows'))<=tolerance
theseCliques(i) = true;
hitID = sum(these(possibleExtensions,ok),2)==currentSize;
theseHits{i} = possiblePartners(hitID);
end
end
cellCliques{j} = combinations(theseCliques,:); cellHits{j} = theseHits(theseCliques,:);
for j=1:size(these,1) % for assembly j
possibleExtensions = [ j+find(sum(these(j+1:end,these(j,:)>0),2)==currentSize-1)];
% remove neurons appearing too few times to possibly participate in a clique
possibleNeurons = sum(these([j;possibleExtensions],:),1)>=1;
possibleExtensions(any(these(possibleExtensions,~possibleNeurons),2)) = [];
possiblePartners = find(sum(these(:,possibleNeurons),2)==currentSize);
combinations = these(possibleExtensions,:); combinations(:,these(j,:)>0) = 1; combinations = unique(combinations,'rows');
theseHits = cell(size(combinations,1),1); theseCliques = false(size(combinations,1),1);
for i=1:size(combinations,1)
ok = combinations(i,:)>0;
comb = SSA_SplitAssemblies(combinations(i,:),3);
if mean(~ismember(comb(:,ok),triplets(:,ok),'rows'))<=tolerance
theseCliques(i) = true;
hitID = sum(these(possibleExtensions,ok),2)==currentSize;
theseHits{i} = possiblePartners(hitID);
end
end
cellCliques{j} = combinations(theseCliques,:); cellHits{j} = theseHits(theseCliques,:);
if verbose && (rem(j,500)==0)
display([datestr(clock) ': So far found ' addComma(sum(cellfun(@(x) size(x,1),cellCliques))) ' cliques after going through ' ...
addComma(j) '/' addComma(size(these,1)) ' groups. toc=' addComma(round(toc))]);
end
end
verbose = false;
for j=1:size(these,1) % for assembly j
possibleExtensions = [ j+find(sum(these(j+1:end,these(j,:)>0),2)==currentSize-1)];
% remove neurons appearing too few times to possibly participate in a clique
possibleNeurons = sum(these([j;possibleExtensions],:),1)>=1;
possibleExtensions(any(these(possibleExtensions,~possibleNeurons),2)) = [];
possiblePartners = find(sum(these(:,possibleNeurons),2)==currentSize);
combinations = these(possibleExtensions,:); combinations(:,these(j,:)>0) = 1; combinations = unique(combinations,'rows');
theseHits = cell(size(combinations,1),1); theseCliques = false(size(combinations,1),1);
for i=1:size(combinations,1)
ok = combinations(i,:)>0;
comb = SSA_SplitAssemblies(combinations(i,:),3);
if mean(~ismember(comb(:,ok),triplets(:,ok),'rows'))<=tolerance
theseCliques(i) = true;
hitID = sum(these(possibleExtensions,ok),2)==currentSize;
theseHits{i} = possiblePartners(hitID);
end
end
cellCliques{j} = combinations(theseCliques,:); cellHits{j} = theseHits(theseCliques,:);
if verbose && (rem(j,500)==0)
display([datestr(clock) ': So far found ' addComma(sum(cellfun(@(x) size(x,1),cellCliques))) ' cliques after going through ' ...
addComma(j) '/' addComma(size(these,1)) ' groups. toc=' addComma(round(toc))]);
end
end
these = these26;
for j=1:size(these,1) % for assembly j
possibleExtensions = [ j+find(sum(these(j+1:end,these(j,:)>0),2)==currentSize-1)];
% remove neurons appearing too few times to possibly participate in a clique
possibleNeurons = sum(these([j;possibleExtensions],:),1)>=1;
possibleExtensions(any(these(possibleExtensions,~possibleNeurons),2)) = [];
possiblePartners = find(sum(these(:,possibleNeurons),2)==currentSize);
combinations = these(possibleExtensions,:); combinations(:,these(j,:)>0) = 1; combinations = unique(combinations,'rows');
theseHits = cell(size(combinations,1),1); theseCliques = false(size(combinations,1),1);
for i=1:size(combinations,1)
ok = combinations(i,:)>0;
comb = SSA_SplitAssemblies(combinations(i,:),3);
if mean(~ismember(comb(:,ok),triplets(:,ok),'rows'))<=tolerance
theseCliques(i) = true;
hitID = sum(these(possibleExtensions,ok),2)==currentSize;
theseHits{i} = possiblePartners(hitID);
end
end
cellCliques{j,2} = combinations(theseCliques,:); cellHits{j} = theseHits(theseCliques,:);
if verbose && (rem(j,500)==0)
display([datestr(clock) ': So far found ' addComma(sum(cellfun(@(x) size(x,1),cellCliques))) ' cliques after going through ' ...
addComma(j) '/' addComma(size(these,1)) ' groups. toc=' addComma(round(toc))]);
end
end
triplets = SSA_SplitAssemblies(tt2,3);
triplets = SSA_SplitAssemblies(these,3);
for j=1:size(these,1) % for assembly j
possibleExtensions = [ j+find(sum(these(j+1:end,these(j,:)>0),2)==currentSize-1)];
% remove neurons appearing too few times to possibly participate in a clique
possibleNeurons = sum(these([j;possibleExtensions],:),1)>=1;
possibleExtensions(any(these(possibleExtensions,~possibleNeurons),2)) = [];
possiblePartners = find(sum(these(:,possibleNeurons),2)==currentSize);
combinations = these(possibleExtensions,:); combinations(:,these(j,:)>0) = 1; combinations = unique(combinations,'rows');
theseHits = cell(size(combinations,1),1); theseCliques = false(size(combinations,1),1);
for i=1:size(combinations,1)
ok = combinations(i,:)>0;
comb = SSA_SplitAssemblies(combinations(i,:),3);
if mean(~ismember(comb(:,ok),triplets(:,ok),'rows'))<=tolerance
theseCliques(i) = true;
hitID = sum(these(possibleExtensions,ok),2)==currentSize;
theseHits{i} = possiblePartners(hitID);
end
end
cellCliques{j,2} = combinations(theseCliques,:); cellHits{j} = theseHits(theseCliques,:);
if verbose && (rem(j,500)==0)
display([datestr(clock) ': So far found ' addComma(sum(cellfun(@(x) size(x,1),cellCliques))) ' cliques after going through ' ...
addComma(j) '/' addComma(size(these,1)) ' groups. toc=' addComma(round(toc))]);
end
end
tt = these15;
tt2 = these26(:,[2:end 1]);
tt3 = tt2(o2,:);
size(o2)
size(tt2)
size(tt1)
size(tt)
sum(ismember(tt,tt2,'rows'))
sum(~ismember(tt,tt2,'rows'))
sum(~ismember(tt2,tt,'rows'))
sum(ismember(tt2,tt,'rows'))
clear theseHits
open orderedAssemblies
find(seek)
these = these26; triplets = SSA_SplitAssemblies(orderedAssemblies,3);
for j=1:size(these,1) % for assembly j
possibleExtensions = [ j+find(sum(these(j+1:end,these(j,:)>0),2)==currentSize-1)];
% remove neurons appearing too few times to possibly participate in a clique
possibleNeurons = sum(these([j;possibleExtensions],:),1)>=1;
possibleExtensions(any(these(possibleExtensions,~possibleNeurons),2)) = [];
possiblePartners = find(sum(these(:,possibleNeurons),2)==currentSize);
combinations = these(possibleExtensions,:); combinations(:,these(j,:)>0) = 1; combinations = unique(combinations,'rows');
theseHits = cell(size(combinations,1),1); theseCliques = false(size(combinations,1),1);
for i=1:size(combinations,1)
ok = combinations(i,:)>0;
comb = SSA_SplitAssemblies(combinations(i,:),3);
if mean(~ismember(comb(:,ok),triplets(:,ok),'rows'))<=tolerance
theseCliques(i) = true;
hitID = sum(these(possibleExtensions,ok),2)==currentSize;
theseHits{i} = possiblePartners(hitID);
end
end
cellCliques{j,2} = combinations(theseCliques,:); cellHits{j} = theseHits(theseCliques,:);
if verbose && (rem(j,500)==0)
display([datestr(clock) ': So far found ' addComma(sum(cellfun(@(x) size(x,1),cellCliques))) ' cliques after going through ' ...
addComma(j) '/' addComma(size(these,1)) ' groups. toc=' addComma(round(toc))]);
end
end
find(cellCliques{178,2})
find(cellCliques{178,1})
find(cellCliques{177,1}(1,:))
find(cellCliques{177,1}(2,:))
size(triplets)
triplets26 = triplets;
[~,orderUnits] = sort(-gt(13,:));
seek26 = seek; ordered26 = orderedAssemblies;
[~,orderUnits] = sort(-gt(13,:));
detectedOrdered = detected(:,orderUnits);
seek = gt(13,orderUnits)>0;
orderedAssemblies = detectedOrdered;
find(seek)
seek15 = seek;
seek15 = seek; ordered15 = orderedAssemblies;
these = these15;
triplets = SSA_SplitAssemblies(ordered15);
for j=1:size(these,1) % for assembly j
possibleExtensions = [ j+find(sum(these(j+1:end,these(j,:)>0),2)==currentSize-1)];
% remove neurons appearing too few times to possibly participate in a clique
possibleNeurons = sum(these([j;possibleExtensions],:),1)>=1;
possibleExtensions(any(these(possibleExtensions,~possibleNeurons),2)) = [];
possiblePartners = find(sum(these(:,possibleNeurons),2)==currentSize);
combinations = these(possibleExtensions,:); combinations(:,these(j,:)>0) = 1; combinations = unique(combinations,'rows');
theseHits = cell(size(combinations,1),1); theseCliques = false(size(combinations,1),1);
for i=1:size(combinations,1)
ok = combinations(i,:)>0;
comb = SSA_SplitAssemblies(combinations(i,:),3);
if mean(~ismember(comb(:,ok),triplets(:,ok),'rows'))<=tolerance
theseCliques(i) = true;
hitID = sum(these(possibleExtensions,ok),2)==currentSize;
theseHits{i} = possiblePartners(hitID);
end
end
cellCliques{j,1} = combinations(theseCliques,:); cellHits{j} = theseHits(theseCliques,:);
if verbose && (rem(j,500)==0)
display([datestr(clock) ': So far found ' addComma(sum(cellfun(@(x) size(x,1),cellCliques))) ' cliques after going through ' ...
addComma(j) '/' addComma(size(these,1)) ' groups. toc=' addComma(round(toc))]);
end
end
sum(tt3(:)==tt(:))
sum(tt3(:)~=tt(:))
these = these15(:,[end 1:end-1]);
triplets15 = SSA_SplitAssemblies(ordered15);
triplets26 = SSA_SplitAssemblies(ordered26);
these = these15(:,[end 1:end-1]); triplets = triplets15(:,[end 1:end-1]);
for j=1:size(these,1) % for assembly j
possibleExtensions = [ j+find(sum(these(j+1:end,these(j,:)>0),2)==currentSize-1)];
% remove neurons appearing too few times to possibly participate in a clique
possibleNeurons = sum(these([j;possibleExtensions],:),1)>=1;
possibleExtensions(any(these(possibleExtensions,~possibleNeurons),2)) = [];
possiblePartners = find(sum(these(:,possibleNeurons),2)==currentSize);
combinations = these(possibleExtensions,:); combinations(:,these(j,:)>0) = 1; combinations = unique(combinations,'rows');
theseHits = cell(size(combinations,1),1); theseCliques = false(size(combinations,1),1);
for i=1:size(combinations,1)
ok = combinations(i,:)>0;
comb = SSA_SplitAssemblies(combinations(i,:),3);
if mean(~ismember(comb(:,ok),triplets(:,ok),'rows'))<=tolerance
theseCliques(i) = true;
hitID = sum(these(possibleExtensions,ok),2)==currentSize;
theseHits{i} = possiblePartners(hitID);
end
end
cellCliques{j,3} = combinations(theseCliques,:); cellHits{j} = theseHits(theseCliques,:);
if verbose && (rem(j,500)==0)
display([datestr(clock) ': So far found ' addComma(sum(cellfun(@(x) size(x,1),cellCliques))) ' cliques after going through ' ...
addComma(j) '/' addComma(size(these,1)) ' groups. toc=' addComma(round(toc))]);
end
end
these = these26; triplets = triplets26;
for j=1:size(these,1) % for assembly j
possibleExtensions = [ j+find(sum(these(j+1:end,these(j,:)>0),2)==currentSize-1)];
% remove neurons appearing too few times to possibly participate in a clique
possibleNeurons = sum(these([j;possibleExtensions],:),1)>=1;
possibleExtensions(any(these(possibleExtensions,~possibleNeurons),2)) = [];
possiblePartners = find(sum(these(:,possibleNeurons),2)==currentSize);
combinations = these(possibleExtensions,:); combinations(:,these(j,:)>0) = 1; combinations = unique(combinations,'rows');
theseHits = cell(size(combinations,1),1); theseCliques = false(size(combinations,1),1);
for i=1:size(combinations,1)
ok = combinations(i,:)>0;
comb = SSA_SplitAssemblies(combinations(i,:),3);
if mean(~ismember(comb(:,ok),triplets(:,ok),'rows'))<=tolerance
theseCliques(i) = true;
hitID = sum(these(possibleExtensions,ok),2)==currentSize;
theseHits{i} = possiblePartners(hitID);
end
end
cellCliques{j,2} = combinations(theseCliques,:); cellHits{j} = theseHits(theseCliques,:);
if verbose && (rem(j,500)==0)
display([datestr(clock) ': So far found ' addComma(sum(cellfun(@(x) size(x,1),cellCliques))) ' cliques after going through ' ...
addComma(j) '/' addComma(size(these,1)) ' groups. toc=' addComma(round(toc))]);
end
end
5
these = these26(o2,:); triplets = triplets26;
for j=1:size(these,1) % for assembly j
possibleExtensions = [ j+find(sum(these(j+1:end,these(j,:)>0),2)==currentSize-1)];
% remove neurons appearing too few times to possibly participate in a clique
possibleNeurons = sum(these([j;possibleExtensions],:),1)>=1;
possibleExtensions(any(these(possibleExtensions,~possibleNeurons),2)) = [];
possiblePartners = find(sum(these(:,possibleNeurons),2)==currentSize);
combinations = these(possibleExtensions,:); combinations(:,these(j,:)>0) = 1; combinations = unique(combinations,'rows');
theseHits = cell(size(combinations,1),1); theseCliques = false(size(combinations,1),1);
for i=1:size(combinations,1)
ok = combinations(i,:)>0;
comb = SSA_SplitAssemblies(combinations(i,:),3);
if mean(~ismember(comb(:,ok),triplets(:,ok),'rows'))<=tolerance
theseCliques(i) = true;
hitID = sum(these(possibleExtensions,ok),2)==currentSize;
theseHits{i} = possiblePartners(hitID);
end
end
cellCliques{j,2} = combinations(theseCliques,:); cellHits{j} = theseHits(theseCliques,:);
if verbose && (rem(j,500)==0)
display([datestr(clock) ': So far found ' addComma(sum(cellfun(@(x) size(x,1),cellCliques))) ' cliques after going through ' ...
addComma(j) '/' addComma(size(these,1)) ' groups. toc=' addComma(round(toc))]);
end
end
these = these26; triplets = triplets26;
for j=1:size(these,1) % for assembly j
possibleExtensions = [ j+find(sum(these(j+1:end,these(j,:)>0),2)==currentSize-1)];
% remove neurons appearing too few times to possibly participate in a clique
possibleNeurons = sum(these([j;possibleExtensions],:),1)>=1;
possibleExtensions(any(these(possibleExtensions,~possibleNeurons),2)) = [];
possiblePartners = find(sum(these(:,possibleNeurons),2)==currentSize);
combinations = these(possibleExtensions,:); combinations(:,these(j,:)>0) = 1; combinations = unique(combinations,'rows');
theseHits = cell(size(combinations,1),1); theseCliques = false(size(combinations,1),1);
for i=1:size(combinations,1)
ok = combinations(i,:)>0;
comb = SSA_SplitAssemblies(combinations(i,:),3);
if mean(~ismember(comb(:,ok),triplets(:,ok),'rows'))<=tolerance
theseCliques(i) = true;
hitID = sum(these(possibleExtensions,ok),2)==currentSize;
theseHits{i} = possiblePartners(hitID);
end
end
cellCliques{j,3} = combinations(theseCliques,:); cellHits{j} = theseHits(theseCliques,:);
if verbose && (rem(j,500)==0)
display([datestr(clock) ': So far found ' addComma(sum(cellfun(@(x) size(x,1),cellCliques))) ' cliques after going through ' ...
addComma(j) '/' addComma(size(these,1)) ' groups. toc=' addComma(round(toc))]);
end
end
cellCliques(:,1) = cellCliques(o2,3);
cellCliques(:,1) = cellCliques(o2,2);
these = these15; triplets = triplets15;
size(these)
cellCliques(:,4) = cellCliques(o2,2);
PlotColorMap(cellfun(@(x) size(x,1),cellCliques);
PlotColorMap(cellfun(@(x) size(x,1),cellCliques));
find(cellCliques{i,3}(1,:))
open o2
o2(102)
find(these26(172,:))
find(tt3(102,:))
find(tt(102,:))
these = these15; triplets = triplets15; k = 1;
j=102;
possibleExtensions = [ j+find(sum(these(j+1:end,these(j,:)>0),2)==currentSize-1)];
possibleNeurons = sum(these([j;possibleExtensions],:),1)>=1;
find(possibleNeurons)
possibleExtensions(any(these(possibleExtensions,~possibleNeurons),2)) = []
possiblePartners = find(sum(these(:,possibleNeurons),2)==currentSize)
combinations = these(possibleExtensions,:); combinations(:,these(j,:)>0) = 1; combinations = unique(combinations,'rows');
find(combinations(1,:))
find(combinations(2,:))
theseHits = cell(size(combinations,1),1); theseCliques = false(size(combinations,1),1);
for i=1:size(combinations,1)
ok = combinations(i,:)>0;
comb = SSA_SplitAssemblies(combinations(i,:),3);
if mean(~ismember(comb(:,ok),triplets(:,ok),'rows'))<=tolerance
theseCliques(i) = true;
hitID = sum(these(possibleExtensions,ok),2)==currentSize;
theseHits{i} = possiblePartners(hitID);
end
end
these = these26; triplets = triplets26; k = 2;
j = o2(102)
possibleExtensions = [ j+find(sum(these(j+1:end,these(j,:)>0),2)==currentSize-1)]
yy = cell2mat(cliques(:,1));
open yy
yy = unique(cell2mat(cliques(:,1)),'rows');
yy2 = unique(cell2mat(cliques(:,2)),'rows');
open yy2
sum(yy(:,1))
sum(yy2(:,1))
sum(yy2(:,2))
tolerant = TolerantDefragmentSSA(ordered26,[],timescale,-Inf,1,0,tolerance);
dbcont
load('temp26.mat')
set26 = setAside;
save('temp26.mat','sums26','cliques26','these26','set26');
tolerant = TolerantDefragmentSSA(ordered15,[],timescale,-Inf,1,0,tolerance);
load('temp-defragmented-04-Feb-2022-15h33-6.mat')
load('temp15.mat')
tolerant = TolerantDefragmentSSA(ordered15,[],timescale,-Inf,1,0,tolerance);
load('temp15.mat')
set15 = setAside;
save('temp15.mat','sums15','cliques15','these15','set15');
load('temp15.mat')
load('temp26.mat')
set = [set15 set26];
open set
tolerant = TolerantDefragmentSSA(ordered15,[],timescale,-Inf,1,0,tolerance);
open tolerant
max(sum(tolerant(:,seek26>0),2))
max(sum(tolerant(:,seek15>0),2))
tolerant = TolerantDefragmentSSA(ordered15,[],timescale,-Inf,1,0,tolerance);
max(sum(tolerant(:,seek15>0),2))
tolerant = TolerantDefragmentSSA(ordered15,[],timescale,-Inf,1,0,tolerance);
max(sum(tolerant(:,seek15>0),2))
tolerant = TolerantDefragmentSSA(ordered15,[],timescale,-Inf,1,0,tolerance);
possibleExtensions = [ j+find(sum(these(j+1:end,these(j,:)>0),2)==currentSize-1)]
possibleExtensions = find(sum(these(j+1:end,these(j,:)>0),2)==currentSize)
possibleExtensions = find(sum(these(j+1:end,these(j,:)>0),2)==currentSize-1)
possibleExtensions = find(sum(these(:,these(j,:)>0),2)==currentSize-1)
tolerant = TolerantDefragmentSSA(ordered15,[],timescale,-Inf,1,0,tolerance);
max(sum(tolerant(:,seek15>0),2))
tolerant = TolerantDefragmentSSA(ordered26,[],timescale,-Inf,1,0,tolerance);
max(sum(tolerant(:,seek26>0),2))
find(set{3,2}(1,:))
find(set{3,4}(1,:))
find(set{3,4}(1,[2:end 1]))
setset = set; for i=3:7, for k=1:2, setset{i,k} = set{i,k}(:,[2:end 1]); end; end
i
k
find(~ismember(setset{3,4},setset{3,2}))
find(~ismember(setset{3,4},setset{3,2},'rows'))
find(setset{3,4}(1,:))
find(setset{3,2}(1,:))
setset = set; for i=3:7, for k=3:4, setset{i,k} = set{i,k}(:,[2:end 1]); end; end
find(setset{3,2}(1,:))
find(setset{3,4}(1,:))
find(~ismember(setset{3,4},setset{3,2},'rows'))
find(~ismember(setset{3,3},setset{3,1},'rows'))
sum(~ismember(these15,set{3,2},'rows'))
q = these15*set{3,2};
q = these15*set{3,2}';
size(q)
sum(max(q,2)~=3)
sum(max(q,[],2)~=3)
tolerant = TolerantDefragmentSSA(ordered26,[],timescale,-Inf,1,0,tolerance);
load('temp15.mat', 'these15')
qq = cell2mat(hits(pass));
open qq
unique(qq);
load('temp26.mat', 'these26')
find(these26(15,:))
find(these26(14,:))
j=14;
these0 = these; these = these26;
possibleExtensions = find(sum(these(:,these(j,:)>0),2)==currentSize-1);
% remove neurons appearing too few times to possibly participate in a clique
possibleNeurons = sum(these([j;possibleExtensions],:),1)>=1;
possibleExtensions(any(these(possibleExtensions,~possibleNeurons),2)) = [];
possiblePartners = find(sum(these(:,possibleNeurons),2)==currentSize);
combinations = these(possibleExtensions,:); combinations(:,these(j,:)>0) = 1; combinations = unique(combinations,'rows');
theseHits = cell(size(combinations,1),1); theseCliques = false(size(combinations,1),1);
for i=1:size(combinations,1)
ok = combinations(i,:)>0;
comb = SSA_SplitAssemblies(combinations(i,:),3);
if mean(~ismember(comb(:,ok),triplets(:,ok),'rows'))<=tolerance
theseCliques(i) = true;
hitID = sum(these(possibleExtensions,ok),2)==currentSize;
theseHits{i} = possiblePartners(hitID);
end
end
cellCliques{j} = combinations(theseCliques,:); cellHits{j} = theseHits(theseCliques,:);
find(cellCliques{14}(1,:))
find(these(1,:))
find(combinations(1,:))
i=1
ok = combinations(i,:)>0;
comb = SSA_SplitAssemblies(combinations(i,:),3);
mean(~ismember(comb(:,ok),triplets(:,ok),'rows'))
tolerance
find(comb(1,ok))
find(comb(1,:))
find(these(15,:))
find(these(14,:))
find(comb(2,:))
find(comb(3,:))
find(comb(4,:))
i=4
i=1
~ismember(comb(:,ok),triplets(:,ok),'rows')
theseCliques(i)
hitID = sum(these(possibleExtensions,ok),2)==currentSize;
hitID = sum(these(possiblePartners,ok),2)==currentSize;
theseHits{i} = possiblePartners(hitID);
j
find(possibleNeurons)
possiblePartners(hitID)
345
tolerant = TolerantDefragmentSSA(ordered26,[],timescale,-Inf,1,0,tolerance);
max(sum(tolerant(:,seek26>0),2))
tolerant = TolerantDefragmentSSA(ordered26,[],timescale,-Inf,1,0,tolerance);
these0 = these;
open these
qq = cell2mat(hits(pass));
open qq
qq = unique(cell2mat(hits(pass)));
these0 = these;
tolerant = TolerantDefragmentSSA(ordered26,[],timescale,-Inf,1,0,tolerance);
these(cell2mat(hits(pass)),:) = [];
tolerant = TolerantDefragmentSSA(ordered26,[],timescale,-Inf,1,0,tolerance);
2
semplot(tolerances,tp,'k');
semplot(tolerances,fp,'r');
nTriplets
semplot(tolerances,tp,'k');
semplot(tolerances,fp,'r');
figure; semplot(tolerances,tp,'k');
semplot(tolerances,fp,'r');
nDetected
find(gt(1,:))
find(gt(2,:))
j=3; i=3
tolerance = tolerances(j);
rng(i)
idxDetected = sort(randperm(nTriplets,nDetected));
detected = gtTriplets(idxDetected,:);
% 3 17 missing
find(gtTriplets(3,:))
find(gtTriplets(17,:))
find(gt(1,:))
find(gt(2,:))
tolerant = TolerantDefragmentSSA(detected,[],timescale,-Inf,1,0,tolerance);
splitt = SSA_SplitAssemblies(tolerant,3);
find(~ismember(gtTriplets,splitt,'rows'))
sum(tolerant,2)
find(tolerant(5,:))
find(tolerant(1,:))
find(tolerant(2,:))
find(tolerant(4,:))
find(tolerant(6,:))
find(detected(:,33)>0)
SSA_SplitAssemblies([1 1 1 1 1],3)
SSA_SplitAssemblies([1 1 1 1 1],4)
SSA_SplitAssemblies([1 1 1 1],4)
SSA_SplitAssemblies([1 1 1 1],3)
SSA_SplitAssemblies([0 1 1 1 1; 1 0 1 1 1; 1 1 0 1 1; 1 1 1 0 1],3)
find(gtTriplets(3,:))
find(gtTriplets(17,:))
temp(detected,[],timescale,-Inf,1,0,tolerance);
currentSize
possibleExtensions = [ j+find(sum(these(j+1:end,these(j,:)>0),2)==currentSize-1)]
q = sum(these(j+1:end,these(j,:)>0),2);
open q
q = sum(these(:,these(j,:)>0),2);
possibleExtensions = find(sum(these(:,these(j,:)>0),2)>0;
possibleExtensions = find(sum(these(:,these(j,:)>0),2)>0);
possibleNeurons = sum(these([j;possibleExtensions],:),1)>=1;
sum(possibleNeurons)
temp(detected,[],timescale,-Inf,1,0,tolerance);
iteration
length(sizes2consider)
iteration
temp(detected,[],timescale,-Inf,1,0,tolerance);
currentSize
any(these(possibleExtensions,~possibleNeurons),2)
possiblePartners = find(sum(these(:,possibleNeurons),2)==currentSize);
possiblePartners = find(sum(these(:,possibleNeurons),2)>1);
combinations = SplitAssemblies(possibleNeurons,currentSize);
for j=1:size(these,1) % for assembly j
%             possibleExtensions = [ j+find(sum(these(j+1:end,these(j,:)>0),2)==currentSize-1)];
possibleExtensions = find(sum(these(:,these(j,:)>0),2)>0);
% remove neurons appearing too few times to possibly participate in a clique
possibleNeurons = sum(these([j;possibleExtensions],:),1)>=1;
if sum(possibleNeurons)<currentSize
continue;
end
possiblePartners = find(sum(these(:,possibleNeurons),2)>1);
combinations = SplitAssemblies(possibleNeurons,currentSize);
theseHits = cell(size(combinations,1),1); theseCliques = false(size(combinations,1),1);
for i=1:size(combinations,1)
ok = combinations(i,:)>0;
comb = SplitAssemblies(combinations(i,:),3);
if mean(~ismember(comb(:,ok),triplets(:,ok),'rows'))<=tolerance
theseCliques(i) = true;
hitID = sum(these(possiblePartners,ok),2)==currentSize;
theseHits{i} = possiblePartners(hitID);
end
end
cellCliques{j} = combinations(theseCliques,:); cellHits{j} = theseHits(theseCliques,:);
if verbose && (rem(j,500)==0)
display([datestr(clock) ': So far found ' addComma(sum(cellfun(@(x) size(x,1),cellCliques))) ' cliques after going through ' ...
addComma(j) '/' addComma(size(these,1)) ' groups. toc=' addComma(round(toc))]);
end
end
temp(detected,[],timescale,-Inf,1,0,tolerance);
cell2mat(cellCliques);
cellCliques(1)
cellCliques(2)
cellCliques(3)
cellCliques(end)
[cellCliques,cellHits] = deal(cell(size(these,1),1));
cellCliques(end)
cellCliques(1)
[cellCliques,cellHits] = deal(repmat({false},(size(these,1),1)));
[cellCliques,cellHits] = deal(repmat({false},size(these,1),1));
[cellCliques,cellHits] = deal(repmat({false(0,50)},size(these,1),1));
temp(detected,[],timescale,-Inf,1,0,tolerance);
currentSize
find(combinations(1,:))
find(combinations(2,:))
temp(detected,[],timescale,-Inf,1,0,tolerance);
if threshold~=-Inf
if verbose,tic; end
for loop = 000:1000:n-1 % do this in batches so that results are periodically saved in the cell, rather than kept by each worker until the end
if n==0, continue; end
if verbose
display([datestr(clock) ': Loop ' addComma(loop) ' out of ' addComma(n) '. toc=' addComma(round(toc))]);
end
parfor j=(1:1000)+loop % length(cellCliques)
thisCombination = combinations(j,:);
pass(j) = RateSSA2(thisCombination,spikes,windowSize,threshold,nMin);
if verbose && rem(j,100)==0, display([datestr(clock) ': ' num2str(j)]); end
end
end
if verbose
display([datestr(clock) ': Last loop containing ' addComma(size(combinations,1)-n) ' cliques. toc=' addComma(round(toc))]);
end
verbose = verbose; % make this variable accessible to the workers
parfor j=(n+1):size(combinations,1)% length(cellCliques)
thisCombination = combinations(j,:);
pass(j) = RateSSA2(thisCombination,spikes,windowSize,threshold,nMin);
if verbose && rem(j,100)==0, display([datestr(clock) ': cell passes ' num2str(j) ' computed.']); end
end
else % If the threshold is -Inf, don't check the first criterion
pass = true(size(combinations,1),1);
end
temp(detected,[],timescale,-Inf,1,0,tolerance);
if threshold~=-Inf
if verbose,tic; end
for loop = 000:1000:n-1 % do this in batches so that results are periodically saved in the cell, rather than kept by each worker until the end
if n==0, continue; end
if verbose
display([datestr(clock) ': Loop ' addComma(loop) ' out of ' addComma(n) '. toc=' addComma(round(toc))]);
end
parfor j=(1:1000)+loop % length(cellCliques)
thisCombination = combinations(j,:);
pass(j) = RateSSA2(thisCombination,spikes,windowSize,threshold,nMin);
if verbose && rem(j,100)==0, display([datestr(clock) ': ' num2str(j)]); end
end
end
if verbose
display([datestr(clock) ': Last loop containing ' addComma(size(combinations,1)-n) ' cliques. toc=' addComma(round(toc))]);
end
verbose = verbose; % make this variable accessible to the workers
parfor j=(n+1):size(combinations,1)% length(cellCliques)
thisCombination = combinations(j,:);
pass(j) = RateSSA2(thisCombination,spikes,windowSize,threshold,nMin);
if verbose && rem(j,100)==0, display([datestr(clock) ': cell passes ' num2str(j) ' computed.']); end
end
else % If the threshold is -Inf, don't check the first criterion
pass = true(size(combinations,1),1);
end
pass
newAssemblies = combinations(pass,:);
q = temp(detected,[],timescale,-Inf,1,0,tolerance);
iteration
currentSize
combinations(ismember(combinations,ignore,'rows'),:) = [];
ignore
size(ignore)
size(mergedCell{1})
size(mergedCell{iteration})
size(mergedCell{iteration-1})
cell2mat(mergedCell);
size(cell2mat(mergedCell))
ignorwe
ignore
islogical(ignore)
islogical(mergedCell{1})
if isempty(ignore), ignore = false(0,nUnits); end
combinations(ismember(combinations,ignore,'rows'),:) = [];
q = temp(detected,[],timescale,-Inf,1,0,tolerance);
open q
q = temp([detected; 1 2 3 4],[],timescale,-Inf,1,0,tolerance);
q = temp([detected; gt(end,:)],[],timescale,-Inf,1,0,tolerance);
size(q)
q = temp([detected; groundTruthAssemblies(:,3)'],[],timescale,-Inf,1,0,tolerance);
open q
tic; q = temp([detected; groundTruthAssemblies(:,3)'],[],timescale,-Inf,1,0,tolerance); end
tic; q = temp([detected; groundTruthAssemblies(:,3)'],[],timescale,-Inf,1,0,tolerance); toc
tic; q = TolerantDefragmentSSA([detected; groundTruthAssemblies(:,3)'],[],timescale,-Inf,1,0,tolerance); toc
sum(q,2)
tic; q = TolerantDefragmentSSA([detected; groundTruthAssemblies(:,3)'],[],timescale,-Inf,1,0,tolerance); toc
tic; q = temp([detected; groundTruthAssemblies(:,3)'],[],timescale,-Inf,1,0,tolerance); toc
profile on
q = temp([detected; groundTruthAssemblies(:,3)'],[],timescale,-Inf,1,0,tolerance);
profile viewer
q = temp([detected; groundTruthAssemblies(:,3)'],[],timescale,-Inf,1,0,tolerance);
mean(~ismember(comb(:,ok),triplets(:,ok),'rows'))
ok
comb = SplitAssemblies(combinations(i,ok),3);
mean(~ismember(comb,triplets(:,ok),'rows'))<=tolerance
mean(~ismember(comb,triplets(:,ok),'rows'))
sum(ok)
comb*triplets;
comb*triplets';
comb'*triplets;
comb'*triplets';
comb*triplets';
triplets*com;
triplets*comb;
triplets*comb';
triplets'*comb;
size(triplets)
comb*triplets(:,ok)';
max(comb*triplets(:,ok)',[],2);
sum(sum(triplets(:,~ok),2)==0)
sum(ismember(comb,triplets(:,ok),'rows'))
tic; for i=1:100, mean(~ismember(comb(:,ok),triplets(:,ok),'rows')); end; toc
tic; for i=1:100, mean(~ismember(comb,triplets(:,ok),'rows')); end; toc
tic; for i=1:1000, mean(~ismember(comb,triplets(:,ok),'rows')); end; toc
tic; for i=1:1000, sum(ismember(comb,triplets(:,ok),'rows')); end; toc
tic; for i=1:100000, sum(ismember(comb,triplets(:,ok),'rows')); end; toc
tic; for i=1:100000, sum(sum(triplets(:,~ok),2)==0); end; toc
(1- sum(sum(triplets(:,~ok),2)==0)/ size(comb,2))
(1- sum(sum(triplets(:,~ok),2)==0)/ size(comb,1))
sum(ok)
size(comb)
nchoosek(12,3)
tic; q = temp([detected; groundTruthAssemblies(:,3)'],[],timescale,-Inf,1,0,tolerance); toc
open q
tic; tolerant = TolerantDefragmentSSA(assemblies,spikes,timescale,-Inf,1,0,tolerance); toc
raly
clf
semplot(tolerances,fp,'r');
semplot(tolerances,tp,'k');
nanmean(tp)
p2z(0.01)
tolerant = TolerantDefragmentSSA(detected,[],timescale,p2z(0.01),1,0,tolerance);
tolerant = TolerantDefragmentSSA(detected,spikes,timescale,p2z(0.01),1,0,tolerance);
i
j
2
figure(1);
semplot(tolerances,tp,'k');
semplot(tolerances,fp,'r');
clf
semplot(tolerances,tp,'k');
semplot(tolerances,fp,'r');
figure; plot(fp(:),tp(:),'.-');
SSA_SplitAssemblies([1 1 1 1 1],4)
load('dataset_20m.mat')
clear all
load('dataset_20m.mat')
threshold = p2z(0.01);
timescale = 0.02;
assemblies = SimSpikeAssemblies(spikes,timescale,threshold);
tic; nontolerant = DefragmentSSA(assemblies,spikes,timescale,threshold,1,0); toc
load('dataset_5h.mat')
assemblies = SimSpikeAssemblies(spikes,timescale,threshold);
tic; nontolerant = DefragmentSSA(assemblies,spikes,timescale,threshold,1,0); toc
tolerance = 0;
tic; tolerant0 = TolerantDefragmentSSA(assemblies,spikes,timescale,threshold,1,0,tolerance); toc
cell2mat(cellCliques);
q = cellfun(@(x) size(x,2),cellCliques);
open q
max(spikes(:,2))
tic; tolerant0 = TolerantDefragmentSSA(assemblies,spikes,timescale,threshold,1,0,tolerance); toc
tolerance = 0.1;
tic; tolerant1 = TolerantDefragmentSSA(assemblies,spikes,timescale,threshold,1,0,tolerance); toc
profile on
tic; tolerant0 = TolerantDefragmentSSA(assemblies,spikes,timescale,threshold,1,1,0); toc
profile on
tic; tolerant0 = TolerantDefragmentSSA(assemblies,spikes,timescale,threshold,1,1,0); toc
tic; tolerant0 = TolerantDefragmentSSA(assemblies,spikes,timescale,threshold,1,1,0); toc
sums =  sum(these([j;possibleExtensions],:),1);
open sums
sums(possibleNeurons);
sum(these(1,:))
nchoosek(currentSize,3)
tolerance
SSA_SplitAssemblies([1 1 1 1],3)
SSA_SplitAssemblies([1 1 1 1 1],3)
nchoosek(5,3)-nchoosek(4,3)
nchoosek(4,3)/nchoosek(5,3)
nchoosek(currentSize-1,3)/nchoosek(currentSize,3)
minSize
minumumNumberOfTriplets = nchoosek(currentSize-1,minSize)/nchoosek(currentSize,minSize)
minumumNumberOfTriplets = nchoosek(currentSize-1,minSize)/nchoosek(currentSize,minSize) * nchoosek(currentSize,minSize)
nchoosek(currentSize-1,minSize)/nchoosek(currentSize,minSize)
minumumNumberOfTriplets = (1-nchoosek(currentSize-1,minSize)/nchoosek(currentSize,minSize)) * nchoosek(currentSize,minSize)
tolerance*nchoosek(currentSize,minSize)
tolerance = 0.1;
tolerance*nchoosek(currentSize,minSize)
(1-nchoosek(currentSize-1,minSize)/nchoosek(currentSize,minSize))
nchoosek(currentSize-1,minSize)/nchoosek(currentSize,minSize)
nchoosek(currentSize-1,minSize)
nchoosek(currentSize,minSize))
nchoosek(currentSize,minSize)
11480-10660
(nchoosek(currentSize,minSize) - nchoosek(currentSize-1,minSize))
minumumNumberOfTriplets = (nchoosek(currentSize,minSize) - nchoosek(currentSize-1,minSize)) ... % this is the number of fragments that may exclude a given neuron
- tolerance*nchoosek(currentSize,minSize); % this is the minumum number of fragments we need to detect to merge them
minumumNumberOfTriplets
currentSize = 5
minumumNumberOfTriplets = (nchoosek(currentSize,minSize) - nchoosek(currentSize-1,minSize)) ... % this is the number of fragments that may exclude a given neuron
- tolerance*nchoosek(currentSize,minSize); % this is the minumum number of fragments we need to detect to merge them
minumumNumberOfTriplets
tolerance
nchoosek(currentSize,minSize)
nchoosek(currentSize-1,minSize)
tolerance*nchoosek(currentSize,minSize)
(nchoosek(currentSize,minSize) - nchoosek(currentSize-1,minSize))
max([1 minumumNumberOfTriplets])
profile on
tic; tolerant0 = TolerantDefragmentSSA(assemblies,spikes,timescale,threshold,1,1,0); toc
profile on
tic; tolerant0 = TolerantDefragmentSSA(assemblies,spikes,timescale,threshold,1,1,0); toc
nchoosek(3,3)
profile on
tic; tolerant0 = TolerantDefragmentSSA(assemblies,spikes,timescale,threshold,1,0,0); toc
iteration
currentSize
assemblies = mergeAssemblies([assemblies origAssemblies]')';
assemblies = mergeAssemblies([assemblies; origAssemblies]')';
dbcont
profile viewer
profile on
tic; tolerant0 = TolerantDefragmentSSA(assemblies,spikes,timescale,threshold,1,0,0); toc
sum(ismember(nontolerant,tolerant0,'rows'))
sum(~ismember(nontolerant,tolerant0,'rows'))
find(~ismember(nontolerant,tolerant0,'rows'))
sum(nontolerant(~ismember(nontolerant,tolerant0,'rows'),:),2)
sum(~ismember(tolerant0,nontolerant,'rows'))
sum(tolerant0(~ismember(tolerant0,nontolerant,'rows'),:),2)
tolerance = 0.1;
tic; tolerant1 = TolerantDefragmentSSA(assemblies,spikes,timescale,threshold,1,0,tolerance); toc
sum(~ismember(tolerant0,tolerant1,'rows'),:),2)
sum(~ismember(tolerant0,tolerant1,'rows'))
tolerance
tic; tolerant3 = TolerantDefragmentSSA(assemblies,spikes,timescale,threshold,1,0,0.3); toc
profile viewer
sum(nontolerant(~ismember(nontolerant,tolerant0,'rows'),:),2)
find((~ismember(nontolerant,tolerant0,'rows'),:))
find((~ismember(nontolerant,tolerant0,'rows')))
find(nontolerant(86,:))
tic; tolerant3 = TolerantDefragmentSSA(assemblies,spikes,timescale,threshold,1,0,0); toc
Unfind([12 13 32 33 53],55)';
find(ans)
tic; tolerant3 = TolerantDefragmentSSA(assemblies,spikes,timescale,threshold,1,0,0); toc
max(sum(these(:,[12 13 32 33 53]),2))
find(sum(these(:,[12 13 32 33 53]),2)==3)
j=12;
possibleExtensions = find(sum(these(:,these(j,:)>0),2)>1);
% remove neurons appearing too few times to possibly participate in a clique
possibleNeurons = sum(these([j;possibleExtensions],:),1)>=max([1 minumumNumberOfTriplets]);
sum(possibleNeurons)
max([1 minumumNumberOfTriplets]
max([1 minumumNumberOfTriplets])
sum(these([j;possibleExtensions],:),1);
open ans
ans([12 13 32 33 53])
minumumNumberOfTriplets
minSize
find(sum(these(:,[12 13 32 33 53]),2)==3)
combinations = SplitAssemblies(Unfind([12 13 32 33 52],55)',currentSize);
currentSize
qq = SplitAssemblies(Unfind([12 13 32 33 52],55)',3);
Portion(ismember(qq,triplets,'rows')
Portion(ismember(qq,triplets,'rows'))
find(ismember(qq,triplets,'rows'))
find(qq(1,:))
find(qq(2,:))
find(qq(3,:))
missing = [Unfind([13 32 52],55)' Unfind([13 33 52],55)' Unfind([33 32 52],55)'];
open missing
missing = [Unfind([13 32 52],55)'; Unfind([13 33 52],55)'; Unfind([33 32 52],55)'];
s = assemblies*missing';
max(s)
cd D:
clear all
load('assemblies_final.mat')
load('triplets.mat')
tic; tolerant3 = TolerantDefragmentSSA(triplets(1:1000,:),spikes,timescale,threshold,1,0,0); toc
timescale = 0.02;
tic; tolerant3 = TolerantDefragmentSSA(triplets(1:1000,:),spikes,timescale,threshold,1,0,0); toc
tic; tolerant3 = TolerantDefragmentSSA(triplets(1:1000,:),spikes,timescale,threshold,1,1,0); toc
TolerantDefragmentSSA(triplets,spikes,timescale,threshold,1,1,0);
possibleExtensions = sum(these(:,these(j,:)>0),2)>1;
possibleNeurons = sum(these([j;possibleExtensions],:),1)>=max([1 minumumNumberOfTriplets]);
TolerantDefragmentSSA(triplets,spikes,timescale,threshold,1,1,0);
origAssemblies = assemblies;
origAssemblies = triplets;
nMin = 1;
verbose = 1;
tolerance
tolerance = 0.1;
TolerantDefragmentSSA
c = origAssemblies'*origAssemblies;
c  =c>0;
[ MC ] = ELSclique( c );
c(eye(size(c))==1) = 0;
[ MC ] = ELSclique( c );
% i started this around 12h
datestr(clock)
save('saved.mat','MC','origAssemblies','-v7.3');
spikes(end,10
spikes(end,10)
spikes(end,1)
spikes(end,1)/3600
clear all
clc
cd C:
v1
load('rez.mat')
cd D:
cd D:
34
raly
for i=1:size(rez.U), meanWaveform(:,:,i) = squeeze(rez.U(:,i,:)) * squeeze(rez.W(:,i,:))';end
meanWaveformGrouped = meanWaveform;
mWaveformRange = squeeze(max(meanWaveformGrouped,[],2) - min(meanWaveformGrouped,[],2));
PlotColorMap(mWaveformRange);
nClusters = rez.ops.Nfilt;
nClusters
for i = 1:nClusters, meanWaveform(:,:,i) = squeeze(rez.U(:,i,:)) * squeeze(rez.W(:,i,:))';end
meanWaveformGrouped = meanWaveform;
for i=1:nClusters
% Define the electrode group as the group with the widest range
thismax = mWaveformRange(:,i)==max(mWaveformRange(:,i)); % logical
groupID(i,1) = find(thismax,1); % index
end
i
mWaveformRange = squeeze(max(meanWaveformGrouped,[],2) - min(meanWaveformGrouped,[],2));
for i=1:nClusters
% Define the electrode group as the group with the widest range
thismax = mWaveformRange(:,i)==max(mWaveformRange(:,i)); % logical
groupID(i,1) = find(thismax,1); % index
end
[signal,noise] = deal(nan(nClusters,1));
for i=1:nClusters
thismax = mWaveformRange(:,i)==max(mWaveformRange(:,i)); % logical
% Get a matrix of differences between the ranges of the signal in each group.
% If the spike is a real spike, the range of the real electrode group will be A LOT higher than the rest.
% This is what we will use to detect noise clusters, where the ranges of electrode groups are more similar.
d = abs(bsxfun(@minus,mWaveformRange(:,i),mWaveformRange(:,i)'));
d(eye(size(d))==1) = nan;
% Difference between range of the electrode closest to the spike with the ranges of all the other electrodes
% When the spike is clearly visible in one elecrode group, this should be very large.
signal(i,1) = nanmean(reshape(d(thismax,~thismax),[],1));
% Difference between ranges of all the other electrodes (where no real spike should be visible)
% This is our estimate of differences to be expected due to noise, rather than a real spike
noise(i,1) = nanmean(reshape(d(~thismax,~thismax),[],1));
end
% signal to noise ratio (snr) below 4 is noisy (imperically defined)
snr = signal./noise;
noisy = snr<3.5;
hist(snr)
hist(snr,100)
Portion(snr<3.5)
findmin(snr)
w = meanWaveform;
PlotCOlorMap(w(:,:,124));
PlotColorMap(w(:,:,124));
PlotColorMap(w(:,:,125));
PlotColorMap(w(:,:,126));
PlotColorMap(w(:,:,127));
plot(w(42,:,127));
plot(sum(sum(w,3),2))
plot(nansum(nansum(w,3),2))
plot(nansum(nansum(w,3),1))
plot(nansum(nansum(w,3)~=0,1))
size(w)
size(rez.U)
size(w)
size(meanWaveform)
clear meanWaveform
for i = 1:nClusters, meanWaveform(:,:,i) = squeeze(rez.U(:,i,:)) * squeeze(rez.W(:,i,:))';end
size(meanWaveform)
q = squeeze(rez.U(:,i,:))
q = squeeze(rez.U(:,i,:));
size(q)
size(q*q')
size(rez.U)
PlotColorMap(squeeze(rez.U(:,1,:))
PlotColorMap(squeeze(rez.U(:,1,:)))
PlotColorMap(squeeze(rez.U(:,1,:))*squeeze(rez.U(:,1,:))')
size(rez.W)
size(meanWaveform)
nClusters
nChannels = size(rez.U,1)
nBins = size(rez.W,1)
meanWaveform = nan(nChannels, nBins, nClusters);
size(meanWaveform)
for i = 1:nClusters, meanWaveform(:,:,i) = squeeze(rez.U(:,i,:)) * squeeze(rez.W(:,i,:))';end
PlotColorMap(squeeze(rez.W(:,1,:)))
PlotColorMap((squeeze(rez.W(:,1,:))==0)+1)
empty = sum(sum(rez.W~=,3),2);
empty = sum(sum(rez.W~=0,3),2);
empty = any(any(rez.W,3),2);
empty = ~any(any(rez.W,3),2);
size(meanWaveform)
nClusters = rez.ops.Nfilt;
nChannels = size(rez.U,1);
% Ignore temporal bins which are always empty (seems to be a Kilosort bug, I haven't investigated)
empty = ~any(any(rez.W,3),2);
nBins = sum(~empty);
% Get the average waveform for each electrode group. We will be comparing these.
meanWaveform = nan(nChannels, nBins, nClusters);
for i = 1:nClusters, meanWaveform(:,:,i) = squeeze(rez.U(:,i,:)) * squeeze(rez.W(~empty,i,:))';end
% Define the waveform range (difference between peak value and trough value).
% The electrode closest to the spike will have the widest range.
mWaveformRange = squeeze(max(meanWaveform,[],2) - min(meanWaveform,[],2));
for i=1:nClusters
% Define the electrode group as the group with the widest range
thismax = mWaveformRange(:,i)==max(mWaveformRange(:,i)); % logical
groupID(i,1) = find(thismax,1); % index
end
i
for i=1:nClusters
% Define the electrode group as the group with the widest range
thismax = mWaveformRange(:,i)==max(mWaveformRange(:,i)); % logical
groupID(i,1) = find(thismax,1); % index
end
i
thismax = mWaveformRange(:,i)==max(mWaveformRange(:,i)); % logical
PlotCOlorMap(meanWaveform(:,:,6));
PlotColorMap(meanWaveform(:,:,6));
sum(rez.st3(:,2)==6)
hist(snr,100)
Portion(noisy)
s = Accumulate(rez.st3(:,2));
plot(snr,s);
plot(snr,s,'.');
plot(s,snr,'.');
plot(log(s),snr,'.');
PlotHVLines(3.5,'h');
log(20)
PlotHVLines(3,'v');
sum(s<20)
sum(s<20 & s>0)
log(50)
sum(s<50 & s>0)
PlotHVLines(4,'h');
plot(nanmean(meanWaveform,3))
plot(nanmean(nanmean(meanWaveform,3),1))
plot(Smooth(nanmean(nanmean(meanWaveform,3),1),1))
mm = squeeze(mean(meanWaveform,1));
d = mm - Smooth(mm,[1 0]);
PlotColorMap(d);
PlotColorMap(mm);
PlotColorMap(mm');
PlotColorMap(zscore(mm)');
semplot(zscore(mm)');
clf
PlotColorMap(zscore(mm)');
PlotColorMap(zscore(mm)' - zscore(Smooth(mm,[1 0]))');
d = mm - Smooth(mm,[1 0]);
size(d)
plot(sum(d,2))
plot(nansum(d,2))
plot(nansum(d,1))
plot(nanmean(d,1))
mm = nanzscore(squeeze(mean(meanWaveform,1)));
d = mm - Smooth(mm,[1 0]);
plot(nanmean(d,1))
plot(mm(1,:))
plot(mm(:,1))
hold all
plot(Smooth(mm(:,1)))
plot(Smooth(mm(:,1),1))
size(d)
plot(d(:,1))
clf
PlotColorMap(d)
plot(nanmean(abs(d),1))
hist(nanmean(abs(d),1))
hist(nanmean(abs(d),1),100)
plot(s,snr,'.');
plot(snr,s);
plot(snr,nanmean(abs(d)));
plot(snr,nanmean(abs(d)),'.');
plot(s,nanmean(abs(d)),'.');
crooked = nanmean(abs(d));
mmax(d(:))
mmax(abs(d(:))_
mmax(abs(d(:)))
mmax(nanmean(abs(d)))
subplot(2,2,4);
plot(abs(d(i,:)))
size(d)
size(meanWaveform)
dd = squeeze(meanWaveform(groupID,:,:));
dd = nan(size(d));
for i=1:nClusters
dd(:,i) = squeeze(meanWaveform(groupID(i),:,i));
end
dd = nan(size(d));
for i=1:nClusters
m = squeeze(meanWaveform(groupID(i),:,i));
dd(:,i) = m - Smooth(m,1);
end
size(m)
size(dd)
i
m = squeeze(meanWaveform(groupID(i),:,i));
size(d)
dd(:,i) = m - Smooth(m,1);
size(m)
crooked2 = nanmean(abs(dd));
clf
hist(crooked,100);
hist(crooked2,100);
plot(crooked,crooked2,'.');
clf
plot(crooked,crooked2,'.');
crooked2 = max(abs(dd));
crooked2 = nanmean(abs(dd));
crooked3 = max(abs(dd));
plot(crooked2,crooked3,'.');
clf
plot(crooked2,crooked3,'.');
hist(crooked2,100);
hist(crooked3,100);
plot(crooked2,crooked3,'.');
ranks = tiedrank([crooked crooked2 crooked3]);
mmax(ranks)
nanmean(ranks)
ranks = tiedrank([crooked crooked2 crooked3]);
ranks = tiedrank([crooked(:) crooked2(:) crooked3(:)]);
nanmean(ranks)
max(ranks)
j
size(m)
crooked4 = max(abs(ddd));
clf
hist(crooked4,100);
plot(crooked2,crooked4,'.');
plot(crooked3,crooked4,'.');
hold all
plot(xlim,xlim,'k--');
lfit(crooked3,crooked4)
plot(xlim,xlim*ans(1)+ans(2),'k--');
[h,ht] = hist([crooked(:) crooked2(:) crooked3(:) crooked4(:)],100);
plot(ht,h);
subplot(2,2,3);
[h,ht] = hist([crooked(:) crooked2(:) crooked3(:) crooked4(:)],100);
plot(ht,h);
allCr = [crooked(:) crooked2(:) crooked3(:) crooked4(:)];
k=1
maxD = min(addCr(:,k)) - addCr(i,k)
allCr = [crooked(:) crooked2(:) crooked3(:) crooked4(:)];
maxD = min(allCr(:,k)) - allCr(i,k)
maxD = [allCr(i,k)-min(allCr(:,k)) max(allCr(:,k)) - allCr(i,k)]
subplot(2,2,3);
allCr = [crooked(:) crooked2(:) crooked3(:) crooked4(:)];
for k=1:4
maxD = max([allCr(i,k)-min(allCr(:,k)) max(allCr(:,k)) - allCr(i,k)]);
[h(:,i),ht] = hist(allCr(:,k),linspace(allCr(i,k)-maxD,allCr(i,k)+maxD,100);
end
plot(linspace(-1,1,length(ht)),h);
title(sum(h(ht<0,:)) ./ sum(h))
min(allCr(:,k))-range(min(allCr(:,k))
[h(:,k),ht] = hist(allCr(:,k),linspace(min(allCr(:,k))-range(allCr(:,k)),max(allCr(:,k))+range(allCr(:,k)),100));
clf; plot(h(:,k))
FindClosest(ht,allCr(i,k))
h(:,k) = circshift(h(:,k),FindClosest(ht,allCr(i,k)));
hold all
plot(h(:,k))
[h(:,k),ht] = hist(allCr(:,k),linspace(min(allCr(:,k))-range(allCr(:,k)),max(allCr(:,k))+range(allCr(:,k)),100));
h(:,k) = circshift(h(:,k),50-FindClosest(ht,allCr(i,k)));
plot(h(:,k))
PlotHVLines(67,'h');
PlotHVLines(67,'v');
j
z = nanzscore(squeeze(mean(meanWaveform,1)));
sum(z(:)-mm(:))
nansum(z(:)-mm(:))
z = [zscore([crooked crooked2 crooked3 crooked4]);
z = zscore([crooked crooked2 crooked3 crooked4]);
mmax(z)
z = nanzscore([crooked crooked2 crooked3 crooked4]);
mmax(z)
[sum(h(ht<0,:)) ./ sum(h); z(i,:)]
i
size(z)
open z
[sum(h(ht<0,:)) ./ sum(h); z(i,:)]
sum(h(ht<0,:)) ./ sum(h)
title([sum(h(ht<0,:)) ./ sum(h); z(i,:)]')
j
figure; hist(zscore(allCr),100);
hist(nanzscore(allCr),100);
[h,ht] = Dist(100,nanzscore(allCr));
plot(ht,h);
mmax(nanzscore(allCr))
max(nanzscore(allCr))
for k=1:5, subplot(2,3,i); plot(ht,h(:,i));
end
for k=1:5, subplot(2,3,k); plot(ht,h(:,i)); end
for k=1:5, subplot(2,3,k); plot(ht,h(:,k)); end
for k=1:5, subplot(2,3,k); hsit(nanzscore(allCr(:,k)),100); end
for k=1:5, subplot(2,3,k); hist(nanzscore(allCr(:,k)),100); end
q = Accumulate(bad+1);
open q
size(q)
q(1,1,1,1,1)
q(2,2,2,2,2)
sum(q(:))
512-430-52
hm = sum(bad,2)>0 & sum(bad,2)<5;
sum(hm)
open hm
find(hm,1)
p2z(0.01)
bad(i,:)
45
k
i
m = nanzscore(squeeze(meanWaveform(groupID(i),:,i))');
dd(:,i) = m - Smooth(m,1);
%     m = [0; m(:); 0];
substitute = m;
for k=2:length(substitute)-1
substitute(k) = mean(m([k-1 k+1]));
end
figure; plot(m); hold all; plot(substitute);
m = [0; m(:); 0];
substitute = m;
for k=2:length(substitute)-1
substitute(k) = mean(m([k-1 k+1]));
end
substitute = substitute(2:end-1);
m = m(2:end-1);
figure; plot(m); hold all; plot(substitute);
j
open cd
figure; hist(crooked6(:),100);
hist(log(crooked6(:)),100);
hist((crooked6(:)),100);
hm = sum(bad,2)>0 & sum(bad,2)<6;
sum(hm)
bad = zs<p2z(0.01);
hm = sum(bad,2)>0 & sum(bad,2)<6;
sum(hm)
sum(bad(hm,:))
figure; plot(bad(hm,:))
plot(sortrows(bad(hm,:)))
PlotColorMap(1+sortrows(bad(hm,:)))
plot(sortrows(bad(hm,:)))
plot(sortrows(bad(hm,:)*1:6))
plot(sortrows(bsxfun(@times,bad(hm,:),1:6)))
plot(sortrows(bsxfun(@times,bad(hm,:),1:6),-6))
plot(sortrows(bsxfun(@times,bad(hm,:),1:6),[-6 -5]))
corr(bad(hm,:))
bad = zs>p2z(0.01);
hm = sum(bad,2)>0 & sum(bad,2)<6;
sum(bad(hm,:))
sum(hm)
qq(i,1) = input
qq(i,1) = input('ok? ','n')
qq(i,1) = input('ok? ','s')
1
open qq
qq(i,1) = num2str(input('ok? ','s'))
2
clear qq
qq(1,1) = num2str(input('ok? ','s'))
5
clear qq
qq(1,1) = str2double(input('ok? ','s'))
2
1
-1
1
-1
1
-1
1
0
-1
0
-1
1
-1
1
-1
1
-1
0
-1
+1
1
0
-1
1
-1
open notes
clf
PlotColorMap(sortrows([notes bad]);
PlotColorMap(sortrows([notes bad]));
PlotColorMap(sortrows([notes bad(hm)]));
PlotColorMap(sortrows([notes bad(hm,:)]));
indices = find(bad);
these = indices(notes>-1);
open these
nn = ones(size(bad(:,1)))*Inf; nn(hm,1) = notes;
PlotColorMap(sortrows([nn(these) bad(these,:)]));
PlotColorMap(sortrows([bad(these,:) -nn(these)]));
PlotColorMap(sortrows([bad(hm,:) nn(hm)]));
PlotColorMap(sortrows([bad(hm,:) nn(hm)],-7));
semplot(bad(nn==1,:));
clf
semplot(1:6,bad(nn==1,:));
semplot(1:6,double(bad(nn==1,:)));
semplot(1:6,double(bad(nn==0,:)),'r');
open nn
sum(nn==0)
sum(notes==0)
bad(nn==0,:)
semplot(1:6,double(bad(nn==1,:)),'b');
sum(nn==1)
bad(nn==1,:)
find(nn==1,2)
find(order==90)
find(sum(bad,2)==0,5)
i=1; figure;
clf
subplot(2,2,1);
PlotColorMap(meanWaveform(:,:,i));
ranks = tiedrank([crooked(:) crooked2(:) crooked3(:) crooked4(:) crooked5(:) crooked6(:)]);
title(ranks(i,:));
subplot(2,4,3);
plot(mm(:,i),'.-','linewidth',2,'markersize',20);
hold all
plot(Smooth(mm(:,i),1),'linewidth',2);
plot(Smooth(mm(:,i),1),'linewidth',2);
title('mean waveform');
subplot(2,4,4);
plot(meanWaveform(groupID(i),:,i),'.-','linewidth',2,'markersize',20);
hold all
plot(Smooth(meanWaveform(groupID(i),:,i),1),'linewidth',2);
m = (squeeze(meanWaveform(groupID(i),:,i))');
dd(:,i) = m - Smooth(m,1);
substitute = m;
for k=2:length(substitute)-1
substitute(k) = mean(m([k-1 k+1]));
end
% plot(m,'linewidth',2);
plot(substitute,'linewidth',2);
plot(substitute,'linewidth',2);
title(['biggest waveform: ' num2str(groupID(i))]);
subplot(2,2,3);
allCr = [crooked(:) crooked2(:) crooked3(:) crooked4(:) crooked5(:) crooked6(:)];
clear h
for k=1:size(allCr,2)
maxD = max([allCr(i,k)-min(allCr(:,k)) max(allCr(:,k)) - allCr(i,k)]);
[h(:,k),ht] = hist(allCr(:,k),linspace(min(allCr(:,k))-range(allCr(:,k)),max(allCr(:,k))+range(allCr(:,k)),100));
h(:,k) = circshift(h(:,k),length(ht)/2-FindClosest(ht,allCr(i,k)));
end
ht =linspace(-1,1,length(ht));
plot(ht,h,'linewidth',2); ylim([0 100]);
PlotHVLines(0,'v','k--');
legend('1','2','3','4','5','6');
z = nanzscore([crooked crooked2 crooked3 crooked4]);
title([sum(h(ht<0,:)) ./ sum(h)])
xlabel(num2str(bad(i,:)));
subplot(2,2,4);
plot(abs(d(:,i)));
hold all
plot(abs(dd(:,i)));
plot(abs(ddd(:,i)),'k');
ylabel(j);
title(zs(i,:));
drawnow
findmin(max(zs,2))
findmin(max(zs,[],2))
min(sum(zs,[],2))
min(sum(zs,2))
findmin(sum(zs,2),50)
findmin(sum(zs,2),5)
sum(bad(:,2))
Portion(bad(:,2))
sum(bad)
p2z(0.01)
p2z([0 1 2 3 4])
p2z([0 0.01 0.001 0.05 0.1 1])
p2z = @(x) sqrt(2) * erfcinv(x);;
p2z = @(p) sqrt(2) * erfcinv(p);;
p2z([0 0.01 0.001 0.05 0.1 1])
size(d)
size(m)
nClusters
nBins
size(d)
size(m)
groupID
Portion(noisy)(
Portion(noisy)
Portion(crooked)
Portion(detectedOnAllElectrodes)
Portion(tooFew)
p2z(0.01)
Portion(crooked)
Portion(noisy)
plot(zvalue,snr,'.');
clf
plot(zvalue,snr,'.');
PlotHVLines(p2z(0.01),'v');
PlotHVLines(3.5,'h');
scatter(zvalue,snr,15,nSpikes,'filled');
clim([0 20]);
open meanWaveform
for i=1:nClusters,thismax = mWaveformRange(:,i)==max(mWaveformRange(:,i)); % logical
groupID = find(thismax,1); % index
waveform = meanWaveform(groupID,:,i);
[qc{i,1},qc{i,2}] = FindLocalMinima(waveform(:));
end
for i=1:nClusters,thismax = mWaveformRange(:,i)==max(mWaveformRange(:,i)); % logical
groupID = find(thismax,1); % index
waveform = meanWaveform(groupID,:,i);
[qc{i,1}] = FindLocalMinima(waveform(:));
end
clear qc
for i=1:nClusters,thismax = mWaveformRange(:,i)==max(mWaveformRange(:,i)); % logical
groupID = find(thismax,1); % index
waveform = meanWaveform(groupID,:,i);
[qc{i,1}] = FindLocalMinima([(1:nBins)' waveform(:)]);
end
size(waveform)
(nBins)
fmat
isdmatrix([1 2; 3 4])
isdmatrix([1 2; 3 4],'@2')
isdmatrix([1 2; 3 4],'>1')
isdmatrix([1 2; 3 4],'>0')
for i=1:nClusters,thismax = mWaveformRange(:,i)==max(mWaveformRange(:,i)); % logical
groupID = find(thismax,1); % index
waveform = meanWaveform(groupID,:,i);
[qc{i,1}] = FindLocalMinima([(1:nBins)' waveform(:)]);
end
open qc
for i=1:nClusters,thismax = mWaveformRange(:,i)==max(mWaveformRange(:,i)); % logical
groupID = find(thismax,1); % index
waveform = meanWaveform(groupID,:,i);
nTroughs = length(FindLocalMinima([(1:nBins)' waveform(:)]))(;
end
for i=1:nClusters,thismax = mWaveformRange(:,i)==max(mWaveformRange(:,i)); % logical
groupID = find(thismax,1); % index
waveform = meanWaveform(groupID,:,i);
nTroughs(i,1) = length(FindLocalMinima([(1:nBins)' waveform(:)]));
end
open nTroughs
thismax = mWaveformRange(:,i)==max(mWaveformRange(:,i)); % logical
groupID = find(thismax,1); % index
waveform = meanWaveform(groupID,:,i);
FindLocalMinima([(1:nBins)' waveform(:)])
figure; plot(waveform);
PlotHVLines([5 11],'v');
for i=1:nClusters,thismax = mWaveformRange(:,i)==max(mWaveformRange(:,i)); % logical
groupID = find(thismax,1); % index
waveform = meanWaveform(groupID,:,i);
nTroughs(i,1) = length(FindLocalMinima([(1:nBins)' Smooth(waveform(:),1,'mode','c')]));
end
for i=1:nClusters,thismax = mWaveformRange(:,i)==max(mWaveformRange(:,i)); % logical
groupID = find(thismax,1); % index
waveform = meanWaveform(groupID,:,i);
nTroughs(i,1) = length(FindLocalMinima([(1:nBins)' Smooth(waveform(:),1,'type','c')]));
end
thismax = mWaveformRange(:,i)==max(mWaveformRange(:,i)); % logical
groupID = find(thismax,1); % index
waveform = meanWaveform(groupID,:,i);
(FindLocalMinima([(1:nBins)' Smooth(waveform(:),1,'type','c')]));
(FindLocalMinima([(1:nBins)' Smooth(waveform(:),1,'type','c')]))
figure; plot(waveform);
plot(Smooth(waveform(:),1,'type','c'))
plot(Smooth(waveform(:),1))
(FindLocalMinima([(1:nBins)' Smooth(waveform(:),1,'type','c')]))
(FindLocalMinima([(1:nBins)' Smooth(waveform(:),1)]))
for i=1:nClusters,thismax = mWaveformRange(:,i)==max(mWaveformRange(:,i)); % logical
groupID = find(thismax,1); % index
waveform = meanWaveform(groupID,:,i);
nTroughs(i,1) = length(FindLocalMinima([(1:nBins)' Smooth(waveform(:),1)]));
end
% Some waveforms look like a noisy oscillation: no single trough:
for i=1:nClusters
thismax = mWaveformRange(:,i)==max(mWaveformRange(:,i)); % logical
groupID = find(thismax,1); % index
waveform = zscore(meanWaveform(groupID,:,i));
troughs = (FindLocalMinima([(1:nBins)' Smooth(waveform(:),1)]));
maxTrough = min(waveform(troughs));
troughs(waveform(troughs)>maxTrough/2) = []; % ignore local minima that are negligible fluctuations compared to the real trough
nTroughs(i,1) = length(troughs);
end
multipleTroughs = nTroughs>1;
% We can also delete clusters with less than a certain number of spikes, e.g. 20
nSpikes = accumarray(rez.st3(:,2),1);
tooFew = nSpikes<20;
noisy = detectedOnAllElectrodes | crooked | multipleTroughs | tooFew;
Portion(noisy)
noisy = detectedOnAllElectrodes | crooked | multipleTroughs | tooFew;
size(noisy)
title(noisy(i));
noisy(i)
figure; hist(nSpikes,1000);
hist(log10(nSpikes),1000);
PlotHVLines(log10(20),'v');
Accumulate(nSpikes(nSpikes<20)+1);
t = rez.st3(rez.st3(:,2)==i);
open t
mmax(t)
mmax(t)/3600
PETH(t,t,'durations',[-1 1]*0.1);
figure; hist(t,100);
diff([19724.2346500000;19725.2439500000;19725.4308500000;19725.7462000000;19726.1443000000;19727.1825500000])
PETH(sort(t),sort(t),'durations',[-1 1]);
PETH(sort(t),sort(t),'durations',[-10 10]);
PETH(sort(t),sort(t),'durations',[-10 10],'show','on');
[h,ht] = PETH(sort(t),sort(t),'durations',[-10 10],'show','on');
open h
semplot(ht,h);
clf;
[h,ht] = PETH(sort(t),sort(t),'durations',[-1 1]*0.1,'show','on');
semplot(ht,h);
h(:,51) = 0;
clf
semplot(ht,h);
i
session
2
open rez
min(rez.st3(:,2))
max(rez.st3(:,2))
i-1
bad = [detectedOnAllElectrodes | crooked | multipleTroughs | tooFew]'
bad = [detectedOnAllElectrodes | crooked | multipleTroughs | tooFew];
sum(bad)
bad = [detectedOnAllElectrodes, crooked, multipleTroughs, tooFew];
sum(bad)
i
find(crooked,5)
i
t = rez.st3(rez.st3(:,2)==i);
t = rez.st3(rez.st3(:,2)==i)/20000;
[h,ht] = PETH(sort(t),sort(t),'durations',[-1 1]*0.1,'show','on'); h(:,51) = 0;
figure; bar(ht,h);
figure; bar(ht,sum(h));
diff(ht(1:2))
[h,ht] = PETH(sort(t),sort(t),'durations',[-1 1]*0.05,'show','on'); h(:,51) = 0;
bar(ht,sum(h));
findmax(zvalue,5)
figure; plot(log(nSpikes),zvalue,'.');
PlotHVLines(p2z(0.01),'h','k--');
zz = zvalue; zz(isnan(zz)) = 0;
findmax(zz,5)
groupID
k
PlotHVLines(p2z(0.01),'v','k--');
PlotHVLines(p2z(0.01),'v','r--');
PlotHVLines(p2z(0.001),'v','r--');
badSpikes = ismember(rez.st3(:,2),find(noisy));
sum(badSpikes)
sum(badSpikes)/1000
sum(~badSpikes)
nPeaksAndTroughs = nan(nClusters,1);
for i=1:nClusters
thismax = mWaveformRange(:,i)==max(mWaveformRange(:,i)); % logical
groupID = find(thismax,1); % index
waveform = zscore(meanWaveform(groupID,:,i));
troughs = (FindLocalMinima([(1:nBins)' Smooth(waveform(:),1)]));
maxTrough = min(waveform(troughs));
troughs(waveform(troughs)>maxTrough/2) = []; % ignore local minima that are negligible fluctuations compared to the real trough
nTroughs = length(troughs);
peaks = (FindLocalMaxima([(1:nBins)' Smooth(waveform(:),1)]));
maxPeak = max(waveform(peaks));
peaks(waveform(peaks)<maxPeak/2) = []; % ignore local minima that are negligible fluctuations compared to the real trough
nPeaks = length(peaks);
nPeaksAndTroughs(i,1) = min([nPeaks nTroughs]);
end
multipleTroughs = nTroughs>1;
sum(multipleTroughs)
sum(nPeaksAndTroughs>1)
multiplePeaksAndTroughs = nPeaksAndTroughs>1;
Portion(nPeaksAndTroughs)
Portion(nPeaksAndTroughs>1)
[maxValue,idx] = max(abs(waveform))
si = sign(waveform(idx))
load('chanMap.mat', 'ycoords')
load('chanMap.mat', 'xcoords')
load('chanMap.mat', 'ycoords')
figure; plot(ycoord);
plot(ycoords);
plot(xcoords,ycoords);
plot(ycoords);
waveform(idx-1) - waveform(idx)
ma = [waveform(idx-1)-waveform(idx) waveform(idx+1)-waveform(idx)]
waveform(idx-1)
si*[waveform(idx-1)-waveform(idx) waveform(idx+1)-waveform(idx)]
-si*[waveform(idx-1)-waveform(idx) waveform(idx+1)-waveform(idx)]
abs(waveform(idx-1)-waveform(idx))<abs(waveform(idx+1)-waveform(idx))
if abs(waveform(idx-1)-waveform(idx))<abs(waveform(idx+1)-waveform(idx))
idx1 = idx-1; else; idx1 = idx+1;
end
waveform([idx1 idx])
waveform(idx) - waveform(idx1)
d = nan(nClusters,1);
for i=1:nClusters
% Define the (closest) channel with the widest waveform range
thismax = mWaveformRange(:,i)==max(mWaveformRange(:,i)); % logical
groupID = find(thismax,1); % index
if isempty(groupID), continue; end
waveform = nanzscore(meanWaveform(groupID,:,i));
[maxValue,idx] = max(abs(waveform));
si = sign(waveform(idx));
% make sure the nearest bins are also the same sign (not a single bin fluctuation)
if abs(waveform(idx-1)-waveform(idx))<abs(waveform(idx+1)-waveform(idx))
idx1 = idx-1; else; idx1 = idx+1;
end
d(i,1) = si*(waveform(idx) - waveform(idx1));
end
i
groupID
waveform = nanzscore(meanWaveform(groupID,:,i));
[maxValue,idx] = max(abs(waveform));
si = sign(waveform(idx));
if abs(waveform(idx-1)-waveform(idx))<abs(waveform(idx+1)-waveform(idx))
idx1 = idx-1; else; idx1 = idx+1;
end
waveform(idx-1)
idx
figure; plot(waveform)
edit CleanRez
d = nan(nClusters,1);
for i=1:nClusters
% Define the (closest) channel with the widest waveform range
thismax = mWaveformRange(:,i)==max(mWaveformRange(:,i)); % logical
groupID = find(thismax,1); % index
if isempty(groupID), continue; end
waveform = nanzscore(meanWaveform(groupID,:,i));
[maxValue,idx] = max(abs(waveform));
if idx==1 || idx==length(waveform), continue; end
% make sure the nearest bins are also the same sign (not a single bin fluctuation)
if abs(waveform(idx-1)-waveform(idx))<abs(waveform(idx+1)-waveform(idx))
idx1 = idx-1; else; idx1 = idx+1;
end
d(i,1) = si*(waveform(idx) - waveform(idx1));
end
open d
figure; hist(d,100);
findmin(d)
for i=124
% Define the (closest) channel with the widest waveform range
thismax = mWaveformRange(:,i)==max(mWaveformRange(:,i)); % logical
groupID = find(thismax,1); % index
if isempty(groupID), continue; end
waveform = nanzscore(meanWaveform(groupID,:,i));
[maxValue,idx] = max(abs(waveform));
if idx==1 || idx==length(waveform), continue; end
% make sure the nearest bins are also the same sign (not a single bin fluctuation)
if abs(waveform(idx-1)-waveform(idx))<abs(waveform(idx+1)-waveform(idx))
idx1 = idx-1; else; idx1 = idx+1;
end
d(i,1) = si*(waveform(idx) - waveform(idx1));
end
plot(waveform)
waveform(idx)
si
d = nan(nClusters,1);
for i=1:nClusters
% Define the (closest) channel with the widest waveform range
thismax = mWaveformRange(:,i)==max(mWaveformRange(:,i)); % logical
groupID = find(thismax,1); % index
if isempty(groupID), continue; end
waveform = nanzscore(meanWaveform(groupID,:,i));
[maxValue,idx] = max(abs(waveform));
si = sign(waveform(idx));
if idx==1 || idx==length(waveform), continue; end
% make sure the nearest bins are also the same sign (not a single bin fluctuation)
if abs(waveform(idx-1)-waveform(idx))<abs(waveform(idx+1)-waveform(idx))
idx1 = idx-1; else; idx1 = idx+1;
end
d(i,1) = si*(waveform(idx) - waveform(idx1));
end
hist(log10(nSpikes),1000);
hist(d,100);
d = nan(nClusters,2);
for i=1:nClusters
% Define the (closest) channel with the widest waveform range
thismax = mWaveformRange(:,i)==max(mWaveformRange(:,i)); % logical
groupID = find(thismax,1); % index
if isempty(groupID), continue; end
waveform = nanzscore(meanWaveform(groupID,:,i));
[maxValue,idx] = max(abs(waveform));
si = sign(waveform(idx));
if idx==1 || idx==length(waveform), continue; end
% make sure the nearest bins are also the same sign (not a single bin fluctuation)
if abs(waveform(idx-1)-waveform(idx))<abs(waveform(idx+1)-waveform(idx))
idx1 = idx-1; else; idx1 = idx+1;
end
d(i,:) = [si*(waveform(idx) - waveform(idx1)) si*waveform(idx1)];
end
PlotXY(d,'.');
PlotHVLines(0,'h');
PlotXY(d(:,[1 3]),'.');
clf
PlotXY(d(:,[1 3]),'.');
PlotXY(d(:,[1 3])./d(:,[3 3]),'.');
PlotXY(d(:,[1 2])./d(:,[3 3]),'.');
PlotHVLines(0,'h');
clf
hist(d(:,1),100);
hist(d(:,2)./d(:,3),100);
FindClosest(d(:,2)./d(:,3),0.5)
ratio= d(:,2)./d(:,3);
ratio(363)
FindClosest(d(:,2)./d(:,3),0.4)
ratio(103)
FindClosest(d(:,2)./d(:,3),0.1)
ratio(13)
this = waveform; this = this(~isnan(this)); this = this(this>quantile(this,0.1) & this<quantile(this,0.9));
z = (waveform-nanmean(this))./(nanstd(this));
% Define the (closest) channel with the widest waveform range
thismax = mWaveformRange(:,i)==max(mWaveformRange(:,i)); % logical
groupID = find(thismax,1); % index
waveform = nanzscore(meanWaveform(groupID,:,i));
[maxValue,idx] = max(abs(waveform));
si = sign(waveform(idx));
if abs(waveform(idx-1)-waveform(idx))<abs(waveform(idx+1)-waveform(idx))
idx1 = idx-1; else; idx1 = idx+1;
end
this = waveform; this = this(~isnan(this)); this = this(this>quantile(this,0.1) & this<quantile(this,0.9));
z = (waveform-nanmean(this))./(nanstd(this));
clf
plto(waveform);
plot(waveform);hold all
plot(z);
open z
this = waveform;
this = this(~isnan(this));
this = waveform'; this = this(~isnan(this));
this = waveform'; this = this(~isnan(this)); this = this(this>=quantile(this,0.1) & this<=quantile(this,0.9));
z = (waveform-nanmean(this))./(nanstd(this));
plot(z);
this = meanWaveform(groupID,:,i)'; this = this(~isnan(this)); this = this(this>=quantile(this,0.1) & this<=quantile(this,0.9));
waveform = (meanWaveform(groupID,:,i)'-nanmean(this))./(nanstd(this));
hist(d(:,2)./d(:,3),100);
clf
hist(d(:,2)./d(:,3),100);
i
ratio= d(:,2)./d(:,3);
ratio(13)
d(13,:)
i=13
this = meanWaveform(groupID,:,i)'; this = this(~isnan(this)); this = this(this>=quantile(this,0.1) & this<=quantile(this,0.9));
waveform = (meanWaveform(groupID,:,i)'-nanmean(this))./(nanstd(this));
[maxValue,idx] = max(abs(waveform));
si = sign(waveform(idx));
if abs(waveform(idx-1)-waveform(idx))<abs(waveform(idx+1)-waveform(idx))
idx1 = idx-1; else; idx1 = idx+1;
end
plot(waveform)
plot(waveform,'.')
nanmean(this)
range(this)
nanstd(this)
q = meanWaveform(groupID,:,i)';
open q
i
groupID
plot(mWaveformRange(:,i));
groupID
thismax = mWaveformRange(:,i)==max(mWaveformRange(:,i)); % logical
groupID = find(thismax,1)
d(13,:)
plot(waveform,'.')
q = meanWaveform(groupID,:,i)';
groupID
PlotColorMap(meanWaveform(:,:,i));
i
q(1:3)
size(q)
q = meanWaveform(1,:,i)';
q(1:3)
q = meanWaveform(2,:,i)';
q(1:3)
PlotColorMap(meanWaveform(:,:,i));
PlotColorMap(meanWaveform(1:5,:,i));
PlotColorMap(meanWaveform(1:5,1:3,i));
meanWaveform(3,2,i)
meanWaveform(1,2,i)
clim
diff(clim)
q = meanWaveform(groupID,:,i)';
plot(q);
plot(zscore(q))
ok =this>=quantile(this,0.1) & this<=quantile(this,0.9);
hold all; plot(ok)
hold all; plot(ok,'.')
d(13,:)
clf
plot(waveform,'.')
this = zscore(meanWaveform(groupID,:,i))'; this = this(~isnan(this));
ok =this>=quantile(this,0.1) & this<=quantile(this,0.9);
hold all
plot(ok,'.-')
this(1)
this(end)
waveform = (meanWaveform(groupID,:,i))'; waveform = waveform-median(waveform);
hold all
plot(waveform,'.')
clf
plot(waveform,'.')
d(13,:)
clf
hist(d(:,2)./d(:,3),100);
hist(d(:,2),100);
diff(quantile(waveform,[0.25 0.75]))
hist(d(:,2),100);
d(13,:)
sum(isnan(d))
findmin(d)
clf
hist(d,100);
open d
hist(d,100);
singleBin = ~ (d>0);
Portion(singleBin)
k
singleBin(i)
d(i)
thismax = mWaveformRange(:,i)==max(mWaveformRange(:,i)); % logical
groupID = find(thismax,1); % index
groupID
waveform = (meanWaveform(groupID,:,i))'; waveform = (waveform-median(waveform))./diff(quantile(waveform,[0.25 0.75]));
[maxValue,idx] = max(abs(waveform));
si = sign(waveform(idx));
idx
sum(singleBin)
k
Portion(noisy)
that = [detectedOnAllElectrodes , singleBin , multiplePeaksAndTroughs , tooFew];
sum(that)
corr(that)
sum(noisy)
edit preprocessV1.m
ops.export.neurosuite
save('rez0.mat','rez');
45
cleanRez = CleanRez(rez);
whos rez
whos cleanRez
1219412382/1374158142
rezToPhy_KSW(rez);
45
rezToPhy_KSW(cleanRez);
2
savepath = pwd;
exists('savepath','file')
exists('savepath','var')
exist('savepath','var')
exist('savepath','file')
exist('savepat5h','file')
fid = fopen(fullfile(savepath,'cluster_group.tsv'),'w');
fwrite(fid, sprintf('cluster_id\t%s\r\n', 'group'));
for i=1:length(indices)
fwrite(fid, sprintf('%d\t%s\r\n', indices(i), 'noise'));
end
open indices
indices = find(noisy);
if exists('savepath','var') % if the user has already exported to phy and wants added labels to noise clusters
fid = fopen(fullfile(savepath,'cluster_group.tsv'),'w');
fwrite(fid, sprintf('cluster_id\t%s\r\n', 'group'));
indices = find(noisy);
for i=1:length(indices)
fwrite(fid, sprintf('%d\t%s\r\n', indices(i), 'noise'));
end
fclose(fid);
end
if exist('savepath','var') % if the user has already exported to phy and wants added labels to noise clusters
fid = fopen(fullfile(savepath,'cluster_group.tsv'),'w');
fwrite(fid, sprintf('cluster_id\t%s\r\n', 'group'));
indices = find(noisy);
for i=1:length(indices)
fwrite(fid, sprintf('%d\t%s\r\n', indices(i), 'noise'));
end
fclose(fid);
end
exist('savepath','folder')
exist('savepath','dir')
exist(fullfile(savepath,'cluster_group.tsv'))
rezToPhy_KSW(rez)
rezToPhy_KSW(rez);
[3 4 6 13]
fclose('all');
% + 1 !!!!
open Ymaze
blockSize = 4;
maxAllowedNumberOfBlocks = 1000;
nLeftPerBlock = 2;
[~, blockMatrix] = sort(rand(maxAllowedNumberOfBlocks, blockSize), 2); % assign random order (the order of random numbers)
key = 2*ones(1:blockSize); key(1:nLeftPerBlock) = 1; % a "key" vector that translates that number order (from 1 to 10) into "left" and "right" (1=left, 2=right).
blockMatrix = key(blockMatrix); % translate the block matrix using the "key" arduinoOutputs into these [1, 2] format
presentationOrder = reshape(blockMatrix', [], 1);
is2 = FindInterval(presentationOrder==2);
is21
is2
open is2
nn = is2(:,2)-is2(:,1)+1;
open nn
hist(nn,10);
[h,ht] = Dist(10,nn);
bar(ht,h);
h
maxAllowedNumberOfBlocks = 100000;
[~, blockMatrix] = sort(rand(maxAllowedNumberOfBlocks, blockSize), 2); % assign random order (the order of random numbers)
key = 2*ones(1:blockSize); key(1:nLeftPerBlock) = 1; % a "key" vector that translates that number order (from 1 to 10) into "left" and "right" (1=left, 2=right).
blockMatrix = key(blockMatrix); % translate the block matrix using the "key" arduinoOutputs into these [1, 2] format
presentationOrder = reshape(blockMatrix', [], 1); % the
is2 = FindInterval(presentationOrder==2);
nn = is2(:,2)-is2(:,1)+1;
hist(nn,10);
[h,ht] = Dist(10,nn);
h
plot(presentationOrder,'.');
plot(presentationOrder);
xlim([0 100]);
plot(presentationOrder,'.-');
xlim([0 100]);
ylim([-1 2]);
ylim([0 3]);
blockSize
blockSize = 10;
[~, blockMatrix] = sort(rand(maxAllowedNumberOfBlocks, blockSize), 2); % assign random order (the order of random numbers)
key = 2*ones(1:blockSize); key(1:nLeftPerBlock) = 1; % a "key" vector that translates that number order (from 1 to 10) into "left" and "right" (1=left, 2=right).
blockMatrix = key(blockMatrix); % translate the block matrix using the "key" arduinoOutputs into these [1, 2] format
presentationOrder = reshape(blockMatrix', [], 1); % the final order is a single vector (we don't care about blocks any more), just "left" vs "right" trials one after the other​
is2 = FindInterval(presentationOrder==2);
nn = is2(:,2)-is2(:,1)+1;
[h,ht] = Dist(10,nn);
figure;  bar(ht,h);
open h
savepath
~exist(fullfile(savepath,'cluster_group.tsv'))
fid = fopen(fullfile(savepath,'cluster_group.tsv'),'a'); % the file already exists. Add labels to the end
indices = find(noisy)-1; % phy notation starts from 0
for i=1:length(indices)
fwrite(fid, sprintf('%d\t%s\r\n', indices(i), 'noise'));
end
fclose(fid);
fclose('all');
rezToPhy_KSW(rez);
open Ymaze
fclose('all');
CleanRez(rez,'savepath',pwd);
open Ymaze
CleanRez(rez,'savepath',pwd);
savepath = pwd
exist('savepath','dir')
exist('savepath','var')
exist('savepath','var') && exist(savepath,'folder')
exist('savepath','var') && exist(savepath,'dir')
CleanRez(rez,'savepath',pwd);
rezToPhy_KSW(rez);
CleanRez(rez,'savepath',pwd);
clean = ans;
sum(clean.st3(:,2)==2)
sum(clean.st3(:,2)==3)
sum(rez.st3(:,2)==3)
raly
kilosortFolder = savepath
fullfile(kilosortFolder,'rez.mat')
load(fullfile(kilosortFolder,'rez.mat'),'rez');
which preprocessSession
function rez = CleanRez(rez,varargin)
% Remove noisy clusters from the rez file. This is an optional preprocessing
% step the user may choose to perform before exporting the Kilosort results to
% phy. Alternatively, one may provide a path (of the Kilosort folder), where
% the function would export the "noise" labels for the noisy clusters in
% phy format. Attention, this would remove any previous labelling.
%
%  USAGE
%
%    [cleanRez] = CleanRez(rez, <options>)
%
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'savepath'        folder where phy files have already been exported
%                       (by default, undefined as you may run CleanRez before
%                       exporting the clean rez to phy). If provided, cluster
%                       labels will be saved in a cluster_groups.tsv file for phy.
%    =========================================================================
%
%  OUTPUT
%
%    rez                cleaned rez file, without the noisy clusters
%
%
% Copyright (C) 2022 by Ralitsa Todorova
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.
% Parse parameter list
for i = 1:2:length(varargin)
if ~ischar(varargin{i})
error(['Parameter ' num2str(i+2) ' is not a property (type ''help <a href="matlab:help CleanRez">CleanRez</a>'' for details).']);
end
switch(lower(varargin{i}))
case 'savepath'
savepath = varargin{i+1};
if ~isfolder(savepath)
error('Incorrect value for property ''savepath'' (type ''help <a href="matlab:help CleanRez">CleanRez</a>'' for details).');
end
otherwise
error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help CleanRez">CleanRez</a>'' for details).']);
end
end
nClusters = rez.ops.Nfilt;
nChannels = size(rez.U,1);
% Ignore temporal bins which are always empty (seems to be a Kilosort bug, I haven't investigated)
empty = ~any(any(rez.W,3),2);
nBins = sum(~empty);
% Get the average waveform for each electrode group. We will be comparing these.
meanWaveform = nan(nChannels, nBins, nClusters);
for i = 1:nClusters, meanWaveform(:,:,i) = squeeze(rez.U(:,i,:)) * squeeze(rez.W(~empty,i,:))';end
% Define the waveform range (difference between peak value and trough value).
% The electrode closest to the spike will have the widest range.
mWaveformRange = squeeze(max(meanWaveform,[],2) - min(meanWaveform,[],2));
[signal,noise] = deal(nan(nClusters,1));
for i=1:nClusters
thismax = mWaveformRange(:,i)==max(mWaveformRange(:,i)); % logical
% Get a matrix of differences between the ranges of the signal in each group.
% If the spike is a real spike, the range of the real electrode group will be A LOT higher than the rest.
% This is what we will use to detect noise clusters, where the ranges of electrode groups are more similar.
d = abs(bsxfun(@minus,mWaveformRange(:,i),mWaveformRange(:,i)'));
d(eye(size(d))==1) = nan;
% Difference between range of the electrode closest to the spike with the ranges of all the other electrodes
% When the spike is clearly visible in one elecrode group, this should be very large.
signal(i,1) = nanmean(reshape(d(thismax,~thismax),[],1));
% Difference between ranges of all the other electrodes (where no real spike should be visible)
% This is our estimate of differences to be expected due to noise, rather than a real spike
noise(i,1) = nanmean(reshape(d(~thismax,~thismax),[],1));
end
% signal to noise ratio (snr) below 4 is noisy (imperically defined)
snr = signal./noise;
detectedOnAllElectrodes = snr<3.5;
% Our measure of segregation is the relative difference diff(a,b)/sum(a,b)
% between signal and noise as defined above
% measure = nanmean((signal-noise)./(signal+noise)); noisy = measure<0.6;
% Make sure the trough/peak is not a single deviation, but part of at least a couple of points deviating together
% (peak/trough on first or last bin not allowed)
d = nan(nClusters,1);
for i=1:nClusters
% Define the (closest) channel with the widest waveform range
thismax = mWaveformRange(:,i)==max(mWaveformRange(:,i)); % logical
groupID = find(thismax,1); % index
if isempty(groupID), continue; end
waveform = (meanWaveform(groupID,:,i))'; waveform = (waveform-median(waveform))./diff(quantile(waveform,[0.25 0.75]));
[maxValue,idx] = max(abs(waveform));
si = sign(waveform(idx));
if idx==1 || idx==length(waveform), continue; end
% make sure the nearest bins are also the same sign (not a single bin fluctuation)
if abs(waveform(idx-1)-waveform(idx))<abs(waveform(idx+1)-waveform(idx))
idx1 = idx-1; else; idx1 = idx+1;
end
d(i,1) = waveform(idx1)/waveform(idx);
end
singleBin = ~ (d>0);
% Some waveforms look like a noisy oscillation: no single trough:
nPeaksAndTroughs = nan(nClusters,1);
for i=1:nClusters
thismax = mWaveformRange(:,i)==max(mWaveformRange(:,i)); % logical
groupID = find(thismax,1); % index
waveform = nanzscore(meanWaveform(groupID,:,i));
if any(~isnan(waveform))
troughs = (FindLocalMinima([(1:nBins)' Smooth(waveform(:),1)]));
maxTrough = min(waveform(troughs));
troughs(waveform(troughs)>maxTrough/2) = []; % ignore local minima that are negligible fluctuations compared to the real trough
nTroughs = length(troughs);
peaks = (FindLocalMaxima([(1:nBins)' Smooth(waveform(:),1)]));
maxPeak = max(waveform(peaks));
peaks(waveform(peaks)<maxPeak/2) = []; % ignore local minima that are negligible fluctuations compared to the real trough
nPeaks = length(peaks);
nPeaksAndTroughs(i,1) = min([nPeaks nTroughs]);
end
end
multiplePeaksAndTroughs = nPeaksAndTroughs>1;
% We can also delete clusters with less than a certain number of spikes, e.g. 20
nSpikes = accumarray(rez.st3(:,2),1);
tooFew = nSpikes<=5;
noisy = detectedOnAllElectrodes | singleBin | multiplePeaksAndTroughs | tooFew;
% Option to grey them out instead of permanently remove
% fid = fopen(fullfile(clusteringpath,'cluster_group.tsv'),'w');
% fwrite(fid, sprintf('cluster_id\t%s\r\n', 'group'));
% for ii=1:length(cids)
%     if any(clu==cids(ii))
%         if any(goodIx==ii)
%             fwrite(fid, sprintf('%d\t%s\r\n', cids(ii), 'good'));
%         elseif any(muaIx==ii)
%             fwrite(fid, sprintf('%d\t%s\r\n', cids(ii), 'mua'));
%         elseif any(noiseIx==ii)
%             fwrite(fid, sprintf('%d\t%s\r\n', cids(ii), 'noise'));
%         end
%     end
% end
% delete noisy clusters
rez.cProj(ismember(rez.st3(:,2),find(noisy)),:) = [];
rez.st3(ismember(rez.st3(:,2),find(noisy)),:) = [];
if exist('savepath','var') && exist(savepath,'dir') % if the user has already exported to phy and wants added labels to noise clusters
display(['Exporting results to ' fullfile(savepath,'cluster_group.tsv')]);
if ~exist(fullfile(savepath,'cluster_group.tsv'))
fid = fopen(fullfile(savepath,'cluster_group.tsv'),'w');
fwrite(fid, sprintf('cluster_id\t%s\r\n', 'group'));
else
fid = fopen(fullfile(savepath,'cluster_group.tsv'),'a'); % the file already exists. Add labels to the end
end
indices = find(noisy)-1; % phy notation starts from 0
for i=1:length(indices)
fwrite(fid, sprintf('%d\t%s\r\n', indices(i), 'noise'));
end
fclose(fid);
end
which preprocessSession
cd(fileparts(ans))
which preprocessSession
cd C:\Users\Cornell\Documents\GitHub\neurocode\preProcessing\
cd M:
cd N:
concatenateDats(pwd,0,1);
sortFiles
~exist('deleteoriginaldatsbool','var')
sortFiles
names2sort = cellfun(@(X) str2num(X(end-5:end)),recordingnames,'UniformOutput',false);
names2sort = cell2mat(names2sort);
isempty(names2sort{1}) && ~isempty(recordingnames{1})
names2sort
recordingnames
names2sort = cellfun(@(X) str2num(X(end-10:end)),recordingnames,'UniformOutput',false);
names2sort = cell2mat(names2sort);
isempty(names2sort)
recordingnames{1}
edit Ymaze
[fList,pList] = matlab.codetools.requiredFilesAndProducts('CleanRez.m');
fList = fList';
open fList
cd C:\Users\Cornell\Documents\GitHub\neurocode\utilities\
edit nanzscore
edit FindLocalMaxima.m
edit FindLocalMinima.m
troughs = (FindLocalMinima([(1:nBins)' waveform(:)]))
troughs = (FindLocalMinima([(1:nBins)' Smooth(waveform(:),1)]))
troughs = strfind([nan;diff(Smooth(waveform(:),1))>0]',[0 1])'
qq = multiplePeaksAndTroughs;
% Some waveforms look like a noisy oscillation: no single peak/trough:
nPeaksAndTroughs = nan(nClusters,1);
for i=1:nClusters
thismax = mWaveformRange(:,i)==max(mWaveformRange(:,i)); % logical
groupID = find(thismax,1); % index
waveform = nanzscore(meanWaveform(groupID,:,i));
if any(~isnan(waveform))
troughs = strfind([nan;diff(Smooth(waveform(:),1))>0]',[0 1])';
maxTrough = min(waveform(troughs));
troughs(waveform(troughs)>maxTrough/2) = []; % ignore local minima that are negligible fluctuations compared to the real trough
nTroughs = length(troughs);
peaks = strfind([nan;diff(Smooth(waveform(:),1))>0]',[1 0])';
maxPeak = max(waveform(peaks));
peaks(waveform(peaks)<maxPeak/2) = []; % ignore local minima that are negligible fluctuations compared to the real trough
nPeaks = length(peaks);
nPeaksAndTroughs(i,1) = min([nPeaks nTroughs]);
end
end
multiplePeaksAndTroughs = nPeaksAndTroughs>1;
sum(qq==multiplePeaksAndTroughs)
which nanmean
cd M:\home\raly\Documents\code\general\
[fList,pList] = matlab.codetools.requiredFilesAndProducts('CleanRez.m');
fList = fList';
clear all
clc
raly
cd N:
load('day7.MergePoints.events.mat')
load('day11.MergePoints.events.mat')
load('day12.MergePoints.events.mat')
load('day12.session.mat')
load('day12.MergePoints.events.mat')
open OJR_20min.batch
raly
N:\Praveen\uLED\uLED1\day7\
X = BatchReturn('OJR_longTraining.batch');
cd(X{1})
session = sessionTemplate(pwd,'showGUI',true); %
cd(X{2})
session = sessionTemplate(pwd,'showGUI',true); %
load('day11.session.mat')
session1 = session;
cd(X{1})
load('day7.session.mat')
sessionu = session;
open sessionu
session = sessionTemplate(pwd,'showGUI',true); %
load('day7.MergePoints.events.mat')
q = MergePoints.foldernames{1};
open q
isnumber(q)
isdouble(q)
ismember(q,'a')
ismember(q,'pre')
ismember(q,'1234567890')
isdate = strfind(ismember(q,'1234567890'),[ones(1,6) 0 ones(1,6)])
isdate = strfind(ismember(q,'1234567890'),[ones(1,6) 0 ones(1,6)]) + (1:13);
isdate = @(x) strfind(ismember(x,'1234567890'),[ones(1,6) 0 ones(1,6)]) + (1:13);
isdate(q)
q(~isdate(q))
isdate = @(x) Unfind(strfind(ismember(x,'1234567890'),[ones(1,6) 0 ones(1,6)]) + (1:13),length(q))
isdate = @(x) Unfind(strfind(ismember(x,'1234567890'),[ones(1,6) 0 ones(1,6)]) + (1:13),length(q));
q(~isdate(q))
isdate = @(x) Unfind(strfind(ismember(x,'1234567890'),[ones(1,6) 0 ones(1,6)]) + (0:12),length(q));
q(~isdate(q))
isdate = @(x) Unfind(strfind(ismember(x,'1234567890'),[ones(1,6) 0 ones(1,6)]) + (-1:12),length(q));
isdate = @(x) accumarray(strfind(ismember(x,'1234567890'),[ones(1,6) 0 ones(1,6)]) + (-1:12),length(q));
q(~isdate(q))
~isdate(q)
isdate(q)
isdate = @(x) accumarray(strfind(ismember(x,'1234567890'),[ones(1,6) 0 ones(1,6)]) + (-1:12));
isdate(q)
isdate = @(x) accumarray((strfind(ismember(x,'1234567890'),[ones(1,6) 0 ones(1,6)]) + (-1:12))',1);
isdate(q)
open accumarray.m
isdate = @(x) accumarray((strfind(ismember(x,'1234567890'),[ones(1,6) 0 ones(1,6)]) + (-1:12))',1,[length(q) 1]));
isdate = @(x) accumarray((strfind(ismember(x,'1234567890'),[ones(1,6) 0 ones(1,6)]) + (-1:12))',1,[length(x) 1]);
isdate(q)
isdate = @(x) logical(accumarray((strfind(ismember(x,'1234567890'),[ones(1,6) 0 ones(1,6)]) + (-1:12))',1,[length(x) 1]));
q(~isdate(q))
qq = q(~isdate(q))
isdate(qq)
strfind(ismember(qq,'1234567890'),[ones(1,6) 0 ones(1,6)])
strfind(ismember(qq,'1234567890'),[ones(1,6) 0 ones(1,6)])+1
strfind(ismember(qq,'1234567890'),[ones(1,6) 0 ones(1,6)])+(1:2)
bsxfun(@plus,strfind(ismember(qq,'1234567890'),[ones(1,6) 0 ones(1,6)]),1:2)
accumarray([],1,[2 1])
existsdate = any(strfind(ismember(x,'1234567890'),[ones(1,6) 0 ones(1,6)]))
existsdate = any(strfind(ismember(q,'1234567890'),[ones(1,6) 0 ones(1,6)]))
existsdate = any(strfind(ismember(qq,'1234567890'),[ones(1,6) 0 ones(1,6)]))
session = sessionTemplate(pwd,'showGUI',true); %
q = day16_220206_092815_Sleep;
q = 'day16_220206_092815_Sleep';
strfind(q,'day')
strfind(q,'day',1)
strfind(q,'day')
min(strfind(q,'day'))
find(~ismember(q(min(strfind(q,'day')):end),'1234567890'),1)
find(~ismember(q(min(strfind(q,'day')+1):end),'1234567890'),1)
find(~ismember(q(min(strfind(q,'day')+2):end),'1234567890'),1)
find(~ismember(q(min(strfind(q,'day')+3):end),'1234567890'),1)
min(strfind(q,'day'):find(~ismember(q(min(strfind(q,'day')+3):end),'1234567890'),1)
min(strfind(q,'day'):find(~ismember(q(min(strfind(q,'day')+3):end),'1234567890'),1))
min(strfind(q,'day')+(0:find(~ismember(q(min(strfind(q,'day')+3):end),'1234567890'),1)))
min(strfind(q,'day'))+(0:find(~ismember(q(min(strfind(q,'day')+3):end),'1234567890'),1)))
min(strfind(q,'day'))+(0:find(~ismember(q(min(strfind(q,'day')+3):end),'1234567890'),1))
min(strfind(q,'day'))+(0:find(~ismember(q(min(strfind(q,'day')+3):end),'1234567890'),1)+1)
accumarray(min(strfind(q,'day'))+(0:find(~ismember(q(min(strfind(q,'day')+3):end),'1234567890'),1)+1)',1,[length(q) 1])
session = sessionTemplate(pwd,'showGUI',true); %
existday(session.epochs{i}.name)
isday(session.epochs{i}.name)
x = session.epochs{i}.name;
isday = @(x) logical(accumarray(min(strfind(x,'day'))+(0:find(~ismember(x(min(strfind(x,'day')+3):end),'1234567890'),1)+1)',1,[length(x) 1]));
session.epochs{i}.name(~isday(session.epochs{i}.name))
i
logical(accumarray(min(strfind(x,'day'))+(0:find(~ismember(x(min(strfind(x,'day')+3):end),'1234567890'),1)+1)',1,[length(x) 1]))
(~isday(session.epochs{i}.name)
(~isday(session.epochs{i}.name))
isday(session.epochs{i}.name)
x
session.epochs{i}.name
session.epochs{i}.name =  MergePoints.foldernames{i};
if existsdate(session.epochs{i}.name), session.epochs{i}.name =  session.epochs{i}.name(~isdate(session.epochs{i}.name)); end % remove the date from the "epochs" suggestion
if existday(session.epochs{i}.name), session.epochs{i}.name =  session.epochs{i}.name(~isday(session.epochs{i}.name)); end % remove the date from the "epochs" suggestion
session.epochs{i}.startTime =  MergePoints.timestamps(i,1);
session.epochs{i}.name
isday = @(x) logical(accumarray(min(strfind(x,'day'))+(0:find(~ismember(x(min(strfind(x,'day')+3):end),'1234567890'),1)+2)',1,[length(x) 1])); % find which characters are part of the dayX/dayXX expression
session.epochs{i}.name =  MergePoints.foldernames{i};
if existsdate(session.epochs{i}.name), session.epochs{i}.name =  session.epochs{i}.name(~isdate(session.epochs{i}.name)); end % remove the date from the "epochs" suggestion
if existday(session.epochs{i}.name), session.epochs{i}.name =  session.epochs{i}.name(~isday(session.epochs{i}.name)); end % remove the date from the "epochs" suggestion
session.epochs{i}.name
session = sessionTemplate(pwd,'showGUI',true); %
cd(X{2})
session = sessionTemplate(pwd,'showGUI',true); %
cd(X{3})
session = sessionTemplate(pwd,'showGUI',true); %
q = session;l
q = session;
load('day12.session.mat')
cd(X{4})
session = sessionTemplate(pwd,'showGUI',true); %
foldername = pwd
basename = pwd
clear a
clear all
X = BatchReturn('OJR_longTraining.batch');
basename = X{1};
X = BatchReturn('OJR_longTraining.batch');
basename = X{1};
raly
exist(fullfile(folder,[dayName '.ripples.events.mat']),'file')
folder=basename;
exist(fullfile(folder,[dayName '.ripples.events.mat']),'file')
swrCh = swrChannels('basepath',basepath);
basepath = basename
swrCh = swrChannels('basepath',basepath);
which today
today = datestr(datenum(clock))
v1
load('day4.ripples.events.mat', 'ripples')
today = datestr(floor(datenum(clock)))
swrCh = swrChannels('basepath',basepath);
%-- 2/15/2022 11:26 PM --%
X = BatchReturn('OJR_longTraining.batch');
basepath = X{1};
load(fullfile(folder,[dayName '.session.mat']),'session');
swrCh = swrChannels('basepath',basepath);
56
cd C;
cd C:
456
which scal2frq
cd raly
fmat
open CCG
which CCG
open CCG
v1
fmat
cd D:
fmat
swrCh = swrChannels('basepath',basepath);
%-- 2/16/2022 8:32 PM --%
X = BatchReturn('OJR_longTraining.batch');
basepath = X{1};
load(fullfile(folder,[dayName '.session.mat']),'session');
swrCh = swrChannels('basepath',basepath);
X = BatchReturn('OJR_longTraining.batch');
basepath = X{1};
load(fullfile(basepath,[dayName '.session.mat']),'session');
[parentFolder,dayName] = fileparts(folder);
[parentFolder,dayName] = fileparts(basepath);
load(fullfile(folder,[dayName '.session.mat']),'session');
load(fullfile(basepath,[dayName '.session.mat']),'session');
swrCh = swrChannels('basepath',basepath);
cd C:
cd(basepath);
swrCh = swrChannels('basepath',basepath);
fmat
cd C:
which FindRipples
swrCh = swrChannels('basepath',basepath);
swrCh.Ripple_Channel = 21; swrCh.Sharpwave_Channel = 30;  swrCh.Noise_Channel = 24;
rippleChannels = swrCh;
rippleChannels.method = 'Manual scoring on neuroscope by Raly Feb 16 2022';
rippleChannels
save('day7.channelinfo.ripples.mat','rippleChannels');
load('day7.channelinfo.ripples.mat')
ripples = DetectSWR([swrCh.ripple swrCh.sharpwave],'saveMat',true);
which swrChannels
exist(fullfile(basepath,[dayName '.lfp']),'file')
exist(fullfile(basepath,[dayName '.deltaWaves.events.mat']),'file')
exist(fullfile(basepath,[dayName '.ripples.events.mat']),'file')
ripples = DetectSWR([swrCh.Ripple_Channel swrCh.Sharpwave_Channel],'saveMat',true);
raly
display('Please review the figure and change the idx1 (ripples) and idx2 (calm non-ripple periods) groups.)
display('Please review the figure and change the idx1 (ripples) and idx2 (calm non-ripple periods) groups.')
display({'Please review the figure and change the idx1 (ripples)', 'and idx2 (calm non-ripple periods) groups.'})
display('Please review the figure and change the idx1 (ripples)', 'and idx2 (calm non-ripple periods) groups.')
disp([''Please review the figure and change the idx1 (ripples)', 'and idx2 (calm non-ripple periods) groupdisp(['])
disp(['Please review the figure and change the idx1 (ripples)', 'and idx2 (calm non-ripple periods) groupdisp(['])
disp({'Please review the figure and change the idx1 (ripples)', 'and idx2 (calm non-ripple periods) groupdisp(['})
disp({'Please review the figure and change the idx1 (ripples)' newline 'and idx2 (calm non-ripple periods) groupdisp(['})
disp(['Please review the figure and change the idx1 (ripples)' newline 'and idx2 (calm non-ripple periods) groupdisp(['])
disp(['Please review the figure and change the idx1 (ripples) ' newline 'and idx2 (calm non-ripple periods) groups.'])
disp(['Please review the figure and change the idx1 (ripples) ' newline 'and idx2 (calm non-ripple periods) groups.' newline ...
'Note that noisy periods need not participate in either group.'])
ripples = DetectSWR([swrCh.Ripple_Channel swrCh.Sharpwave_Channel],'saveMat',true,'forceDetect',true,'check',true);
ripples = DetectSWR([swrCh.Ripple_Channel swrCh.Sharpwave_Channel],'basepath',basepath,'saveMat',true,'forceDetect',true,'check',true);
find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim))
featureTs(find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim)))
featureTs(find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim)))/20000
ripPowerAll(find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim)))/20000
ripPowerAll(find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim)))
featureTs(find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim)))/20000
featureTs(find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim))+1)/20000
featureTs(find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim))+2)/20000
featureTs(find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim))+5)/20000
FindClosest(featureTs/20000,442.394)
featureTs(31253)
featureTs(31253)/20000
ripPowerAll(31253)
swDiffAll(31253)
hold all
hold all; i=31253; plot(ripPowerAll(i),swDiffAll(i),'yo');
hold all; i=31253; plot(ripPowerAll(i),swDiffAll(i),'ko');
hold all; i=31253; plot(swDiffAll(i),ripPowerAll(i),'ko');
clear swrCh
rippleChannels.Ripple_Channel = 21+1; rippleChannels.Sharpwave_Channel = 30+1;  rippleChannels.Noise_Channel = 24+1;
ripples = DetectSWR([swrCh.Ripple_Channel swrCh.Sharpwave_Channel],'basepath',basepath,'saveMat',true,'forceDetect',true,'check',true);
save('day7.channelinfo.ripples.mat','rippleChannels');
ripples = DetectSWR([rippleChannels.Ripple_Channel rippleChannels.Sharpwave_Channel],'basepath',basepath,'saveMat',true,'forceDetect',true,'check',true);
figure; DensityMap(swDiffAll(:),ripPowerAll(:));
ok = swDiffAll<4000 & ripPowerAll<200;
DensityMap(swDiffAll(ok),ripPowerAll(ok));
ok = swDiffAll<3000 & ripPowerAll<150;
DensityMap(swDiffAll(ok),ripPowerAll(ok));
DensityMap(swDiffAll(ok),ripPowerAll(ok),'smooth',0);
ok = swDiffAll<3000 & ripPowerAll<100;
DensityMap(swDiffAll(ok),ripPowerAll(ok),'smooth',0);
ok = swDiffAll<1500 & ripPowerAll<100;
DensityMap(swDiffAll(ok),ripPowerAll(ok),'smooth',0);
clf
plot(swDiffAll(idx2),ripPowerAll(idx2),'k.','markersize',1); hold on
plot(swDiffAll(idx1),ripPowerAll(idx1),'r.','markersize',1); legend('non-ripples','ripples');
featureTs(find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim)))/20000
featureTs(end)/20000
hist(featureTs);
clf
hist(featureTs);
hist(featureTs,100);
max(blocks)
max(blocks)/20000
max(blocks)/1250
featureTs(find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim)))/1250
featureTs(find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim)))/1250clf
clf
clf
plot(swDiffAll(idx2),ripPowerAll(idx2),'k.','markersize',1); hold on
plot(swDiffAll(idx1),ripPowerAll(idx1),'r.','markersize',1); legend('non-ripples','ripples');
featureTs(find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim)))/1250
ans - 15344
featureTs(find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim))+1)/1250
ans - 15344
find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim))
featureTs(71894+(-1:1))/1250
featureTs(71894+(-1:1))/1250 - 15343
featureTs(find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim)))/1250
hold all
find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim))
i = 728; plot(swDiffAll(i),ripPowerAll(i),'yo','linewidth',2);
find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim))
featureTs(find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim)))/1250
featureTs(find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim)))/1250-3976
i = 18623; plot(swDiffAll(i),ripPowerAll(i),'rx','linewidth',2);
find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim))
featureTs(find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim)))/1250
featureTs(find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim)))/1250-18250
i = 85548; plot(swDiffAll(i),ripPowerAll(i),'yo','linewidth',2);
old = [idx1 idx2];
clf
hist(ripPowerAll);
hist(ripPowerAll,1000);
hist(swDiffAll,1000);
hist(ripPowerAll,1000);
PlotHVLines(68.9,'v');
idx1 = ripPowerAll>68.9;
idx2 = ripPowerAll<68.9;
idx2 = ripPowerAll<=68.9;
clf
plot(swDiffAll(idx2),ripPowerAll(idx2),'k.','markersize',1); hold on
plot(swDiffAll(idx1),ripPowerAll(idx1),'r.','markersize',1); legend('non-ripples','ripples');
find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim))
featureTs(39564)/1250
featureTs(39564)/1250-8447
featureTs(58815)/1250
featureTs(58815)/1250-12557
featureTs(58815)/1250-12556
hold all
i = 58815; plot(swDiffAll(i),ripPowerAll(i),'yo','linewidth',2);
i = 39564; plot(swDiffAll(i),ripPowerAll(i),'rx','linewidth',2);
find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim))
featureTs(35105)/1250
featureTs(35105)/1250 - 7490
find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim))
featureTs(31942)/1250
featureTs(31942)/1250 - 6813
featureTs(31942 + (-1:1))/1250
featureTs(31942 + (-1:1))/1250 - 6813
ripPowerAll(31942 + (-1:1))
idx1(ripPowerAll>250) = 0;
idx1(swDiffAll>6000) = 0;
idx2(swDiffAll>2500) = 0;
clf
plot(swDiffAll(idx2),ripPowerAll(idx2),'k.','markersize',1); hold on
plot(swDiffAll(idx1),ripPowerAll(idx1),'r.','markersize',1); legend('non-ripples','ripples');
save('
raly
save('temp_indices_uLED1_day7.mat','idx1','idx2');
dbcont
r = ripples.timestamps(:,1:2);
load(fullfile(basepath,[dayName '.cell_metrics.cellinfo.mat']),'cell_metrics');
nNeurons = length(cell_metrics.spikes.times);
regions = cell_metrics.brainRegion;
regionNames = unique(regions)
regionCell = zeros(length(regions),1);
for i=1:length(regionNames)
regionCell(strcmp(regions,regionNames{i})) = i;
end
spikesCell = cell_metrics.spikes.times';
PFCindex = find(strcmp(regionNames,'PFC'));
pfc = Group(spikesCell{regionCell==PFCindex});
hpc = Group(spikesCell{regionCell==find(strcmp(regionNames,'CA1'))});
neurons = true;
cell_metrics = ProcessCellMetrics('session',session,'manualAdjustMonoSyn',false);
spikesCell = cell_metrics.spikes.times';
spikes = sortrows(Group(spikesCell{:}));
close all
load(fullfile(basepath,[dayName '.cell_metrics.cellinfo.mat']),'cell_metrics');
nNeurons = length(cell_metrics.spikes.times);
regions = cell_metrics.brainRegion;
regionNames = unique(regions)
regionCell = zeros(length(regions),1);
for i=1:length(regionNames)
regionCell(strcmp(regions,regionNames{i})) = i;
end
spikesCell = cell_metrics.spikes.times';
PFCindex = find(strcmp(regionNames,'PFC'));
pfc = Group(spikesCell{regionCell==PFCindex});
hpc = Group(spikesCell{regionCell==find(strcmp(regionNames,'CA1'))});
neurons = true;
spikes = sortrows(Group(spikesCell{:}));
cell_metrics.brainRegion
PETH(spikes(:,1),r(:,1));
fmat
PETH(spikes(:,1),r(:,1));
clear all
X = BatchReturn('OJR_longTraining.batch');
cd(X{1})
cd(X{2})
load('day12.ripples.events.mat')
load('day12.cell_metrics.cellinfo.mat')
spikesCell = cell_metrics.spikes.times';
spikes = sortrows(Group(spikesCell{:}));
load('day12.ripples.events.mat')
fmat
PETH(spikes(:,1),r(:,1));
r = ripples.timestamps(:,1:2);
PETH(spikes(:,1),r(:,1));
nNeurons = length(cell_metrics.spikes.times);
regions = cell_metrics.brainRegion;
regionNames = unique(regions)
regionCell = zeros(length(regions),1);
for i=1:length(regionNames)
regionCell(strcmp(regions,regionNames{i})) = i;
end
hpc = Group(spikesCell{regionCell==find(strcmp(regionNames,'CA1'))});
PETH(hpc(:,1),r(:,1));
fmat
PETH(hpc(:,1),r(:,1));
basepath = X{2};
basepath = X{3};
basepath
rippleChannels = swrChannels('basepath',basepath);
1
3/0.011
3/0.027
5/0.028
SaveCustomEvents('ripples.rip.evt');
SaveCustomEvents('ripples.rip.evt',r,{'ripple start','ripple stop'});
4784
4787
d = diff(ripples,[],2);
d = diff(r,[],2);
clf
plot(f);
plot(d)
load('day12.session.mat')
load('day12.MergePoints.events.mat')
pre = r(:,1)<MergePoints.timestamps(4,1);
post = r(:,1)>MergePoints.timestamps(6,1);
anovabar(d,post);
tic; lfpstruct = getLFP(114+1); display(['loaded! @' num2str(toc)]);
%     tic; lfpstruct = getLFP(104); HMC is 104, OR18 is 97
lfp = [lfpstruct.timestamps]; lfp(:,2) = lfpstruct.data;
clean = CleanLFP(lfp);
deltas0 = FindDeltaPeaks(clean);
figure; PlotXY(lfp);
which ConsolidateIntervals
fmay
fmat
which CleanLFP
clean = CleanLFP(lfp);
diff(noisyInterval)
threshold2
thresholds
plot(t,z);
plot(t,abs(z));
dur(artefactInterval
dur(artefactInterval)
PlotIntervals(artefactInterval,'color','r');
hold on
in = InIntervals(t,artefactInterval);
plot(t(~in),abs(z(~in)),'k');p
clf
size(d)
size(t)
plot(t,d);
plot(t,abs(d));
thresholds
[noisyInterval(:,1)-aroundArtefact2 noisyInterval(:,2)+aroundArtefact2]
diff(ans)
ConsolidateIntervals([1 2])
ConsolidateIntervals([noisyInterval(:,1)-aroundArtefact2 noisyInterval(:,2)+aroundArtefact2])
if ~isempty(noisyInterval),noisyInterval = ConsolidateIntervals([noisyInterval(:,1)-aroundArtefact2 noisyInterval(:,2)+aroundArtefact2]); bad = bad | InIntervals(t,noisyInterval);end
dbcont
clf
PlotXY(clean);
CleanLFP(lfp);
ConsolidateIntervals([1 2])
deltas0 = FindDeltaPeaks(clean);
[q,qt] = PETH(spikes(:,1),deltas(:,2));
[q,qt] = PETH(spikes(:,1),deltas0(:,2));
semplot(qt,q);
clf
PlotColorMap(sortby(q,deltas(:,5)),'x',qt);
PlotColorMap(sortby(q,deltas0(:,5)),'x',qt);
sortby([1 2]',[1 2]')
sum(deltas0(:,5)<2)
PlotColorMap(Shrink(sortby(q,deltas0(:,5)),5,1),'x',qt);
PlotColorMap(Shrink(sortby(deltas0(:,5),deltas0(:,5)),5,1),'x',qt);
PlotColorMap(Shrink(sortby(q,deltas0(:,5)),5,1),'x',qt);
PlotColorMap(Shrink(sortby(q,deltas0(:,5)),100,1),'x',qt);
hist(deltas(:,5),100);
hist(deltas0(:,5),100);
PlotColorMap(Shrink(sortby(q,deltas0(:,5)),100,1),'x',qt);
PlotColorMap(Shrink(sortby(q,deltas0(:,5)),10,1),'x',qt);
PlotColorMap(Shrink(sortby(q,deltas0(:,5)),1,1),'x',qt);
PlotColorMap(Smooth(Shrink(sortby(q,deltas0(:,5)),1,1),[5 0]),'x',qt);
PlotColorMap(Smooth(Shrink(sortby(q,-deltas0(:,5)),1,1),[5 0]),'x',qt);
findmax(deltas0(:,5),150);
[findmax(deltas0(:,5),150) deltas0(findmax(deltas0(:,5),150),5)];
PlotColorMap(Smooth(Shrink(sortby(q,-deltas0(:,5)),1,1),[2 0]),'x',qt);
ylim([0 200])
deltas = deltas0(deltas0(:,5)>2 & deltas(:,5)<4.5,:);
deltas = deltas0(deltas0(:,5)>2 & deltas0(:,5)<4.5,:);
PlotColorMap(Smooth(Shrink(sortby(q,deltas0(:,5)),1,1),[2 0]),'x',qt);
[findmin(deltas0(:,5),100) deltas0(findmin(deltas0(:,5),100),5)];
figure; hist(deltas(:,5))
hist(deltas0(:,5))
hist(deltas0(:,5),100)
hist(deltas0(:,5),1000)
PlotHVLines(1.27,'v');
PlotHVLines(1.5,'v');
sum(deltas0(:,5)<1.5)
PlotHVLines(1975,'h');
deltas = deltas0(deltas0(:,5)>1.5 & deltas0(:,5)<4.5,:);
FindClosest(deltas0(:,2),8257.46)
deltas0(3617 + (-2:10),:);
q = deltas0(3617 + (-2:10),:); q(:,1:3) = q(:,1:3)-8257;
open q
q = deltas0(3617 + (-2:10),:); q(:,1:3) = q(:,1:3)-8256;
q = deltas0(3617 + (-2:10),:); q(:,1:3) = q(:,1:3)-8262;
q = deltas0(3617 + (-2:10),:); q(:,1:3) = q(:,1:3)-8264;
FindClosest(deltas0(:,2),8265.6)
deltas0(3628,:)
deltas0(3628,4:6)
figure; PlotXY(deltas0(:,5:6),'.');
DensityMap(deltas0(:,5),deltas0(:,6));
DensityMap(deltas0(:,5),deltas0(:,6),'nBins',500);
DensityMap(deltas0(:,5),deltas0(:,6),'nBins',1000);
q = deltas0(3617 + (-2:10),:); q(:,1:3) = q(:,1:3)-8256;
q = deltas0(3617 + (-2:10),:); q(:,1:3) = q(:,1:3)-8257;
deltas0(:,7) = (deltas(:,5)-deltas(:,6))./(deltas(:,3)-deltas(:,2));
deltas0(:,7) = (deltas0(:,5)-deltas0(:,6))./(deltas0(:,3)-deltas0(:,2));
q = deltas0(3617 + (-2:10),:); q(:,1:3) = q(:,1:3)-8257;
figure; hist(deltas0(:,end),100);
[q,qt] = PETH(lfp,deltas0(:,2));
clf
semplot(q)
clf
PlotColorMap(sortby(q,deltas0(:,5)))
PlotColorMap(Shrink(sortby(q,deltas0(:,5)),100,1))
PlotColorMap(Shrink(sortby(q,deltas0(:,6)),100,1))
PlotColorMap(Shrink(sortby(q,deltas0(:,7)),100,1))
PlotColorMap(Shrink(sortby(q,deltas0(:,7)),1,1))
clim
PlotColorMap(Shrink(sortby(q,deltas0(:,7)),1,1)); clim(clim/2)
PlotColorMap(Shrink(sortby(q,deltas0(:,7)),1,1)); clim(clim/10)
PlotColorMap(Shrink(sortby(q,deltas0(:,7)),1,1)); clim(clim/5)
findmin(deltas0(:,7))
deltas0(2288,2)
floor(deltas0(2288,2))
addComma(deltas0(2288,2))
addComma(deltas0(2288,1:3))
addComma(deltas0(2288,1))
addComma(deltas0(2288,3))
clf
hist(deltas0(:,5),1000)
hist(deltas0(:,7),1000)
[findmin(deltas0(:,7),100) deltas0(findmin(deltas0(:,7),100),7)];
[findmin(deltas0(:,7),10000) deltas0(findmin(deltas0(:,7),10000),7)];
PlotHVLines(1200,'h');
PlotHVLines(1200,'h','k--');
deltas = deltas0(((deltas0(:,5)-deltas0(:,6))./(deltas0(:,3)-deltas0(:,2)))>20,:);
deltaWaves.timestamps = deltas(:,[1 3]); deltaWaves.peaks = deltas(:,2);  deltaWaves.peakNormedPower = deltas(:,5); deltaWaves.detectorName = ['channel 114(+1), CleanLFP, FindDeltaPeaks, (peak-trough)/dt > 20'];
save(fullfile(basepath,[dayName '.deltaWaves.events.mat']),'deltaWaves');
basepath
[parentFolder,dayName] = fileparts(basepath);
[~,projectName] = fileparts(parentFolder);
sessionID = [projectName '_' dayName];
cd(basepath);
save(fullfile(basepath,[dayName '.deltaWaves.events.mat']),'deltaWaves');
SaveCustomEvents('deltas.del.evt',deltas(:,1:3),{'deltas start','delta peak','deltas stop'});
close all
edit PETH
fmat
edit PETH
EPTH(spikes(:,1),deltas(:,2));
PETH(spikes(:,1),deltas(:,2));
which Bin
fmat
figure;  PETH(spikes(:,1),deltas(:,2));
which Bin
figure;  PETH(spikes(:,1),deltas(:,2));
figure; [q,qt] = PETH(spikes(:,1),deltas(:,2),'show','on');
open qt
qtBad = qt;
figure; [q,qt] = PETH(spikes(:,1),deltas(:,2),'show','on');
profile on
figure; [q,qt] = PETH(spikes(:,1),deltas(:,2),'show','on');
profile viewer
which Sync
which SyncHist
which PETH
figure; [q,qt] = PETH(spikes(:,1),deltas(:,2),'show','on');
which Accumulate
figure; [q,qt] = PETH(spikes(:,1),deltas(:,2),'show','on');
close all
which unique
which Sync
which SyncHist
which Bin
which Accumulate
figure; [q,qt] = PETH(spikes(:,1),deltas(:,2),'show','on');
open Smooth
which isstring
which isastring
today
close all
PETH(hpc(:,1),deltas(:,2));
PETH(hpc(:,1),r(:,1));
PETH(hpc(:,1),deltas(:,2));
cd(basepath)
interval = [0 1]+18044.1;
hist(Restrict(hpc(:,1),interval),100);
% try to run ripple detection with channels [9 10]+1
hist(Restrict(spikes(:,1),interval),100);
PETH(spikes(:,1),deltas(:,2));
PETH(spikes(:,1),r(:,1));
rippleChannels.Ripple_Channel = 9+1; rippleChannels.Sharpwave_Channel = 10+1;
rippleChannels.method = 'Manual scoring on neuroscope by Raly Feb 17 2022';
save('day12.channelinfo.ripples.mat','rippleChannels');
ripples = DetectSWR([rippleChannels.Ripple_Channel rippleChannels.Sharpwave_Channel],'basepath',basepath,'saveMat',true,'forceDetect',true,'check',true);
close all
plot(swDiffAll(idx2),ripPowerAll(idx2),'k.','markersize',1); hold on
plot(swDiffAll(idx1),ripPowerAll(idx1),'r.','markersize',1); legend('non-ripples','ripples');
disp(['Please review the figure and change the idx1 (ripples) ' newline 'and idx2 (calm non-ripple periods) groups.' newline ...
'Note that noisy periods need not participate in either group.']);
keyboard
num2str(t(find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim))))
rp =ripPowerAll; t = featureTs/1250; s = swDiffAll;
FindClosest(t,7978.32)
num2str(t(35499))
rp(35499)
hold all
i = 35499; plot(swDiffAll(i),ripPowerAll(i),'yo','linewidth',2);
figure; hist(s(:),1000);
hist(log(s(:)),1000);
sw = interp1(tl,lfpLow(:,2),t);
tl = (1:length(lfp))'/1250;
lfpLow     = firfilt( lfp, hSw2 );      % lowpass filter
eegLo      = firfilt( lfpLow, hSw1 );   % highpass filter
% lfpLow     = lfpLow -.lfpLo;            % difference of Gaussians
% clear .lfpLo                          %DL: not sure what this '.' was doing here...
lfpLow     = lfpLow - eegLo;            % di
sw = interp1(tl,lfpLow(:,2),t);
plot(sw,s,'.');
hold all
i = 35499; plot(sw(i),s(i),'yo','linewidth',2);
clf
plot(sw,s,'.','markersize',1);
hold all
i = 35499; plot(sw(i),s(i),'yo','linewidth',2);
clf
b = interp1(tl,lfpLow(:,1),t);
plot(b,sw,'.','markersize',1);
hold all
i = 35499; plot(b(i),sw(i),'yo','linewidth',2);
PlotHVLines(0,'h','k--');
PlotHVLines(0,'v','k--');
plot(xlim,xlim,'k');
find(InIntervals(b,xlim) & InIntervals(sw,ylim))
num2str(t(128485))
4/0.037
figure; plot(tl,lfpLow);
hm = FindLocalMinima([tl lfpLow(:,2))]);
hm = FindLocalMinima([tl lfpLow(:,2)]);
open hm
hm(:,2) = interp1(tl,lfpLow(:,2)-lfpLow(:,1));
hm(:,2) = interp1(tl,lfpLow(:,2)-lfpLow(:,1),hm);
PlotHVLines(Restrict(hm,xlim),'v');
hm(:,3) = interp1(tl,lfpLow(:,2),hm);
hm(:,3) = interp1(tl,lfpLow(:,2),hm(:,1));
clf
hist(hm(:,3),100)
hist(hm(:,3),1000)
hm0 = hm;
hm(hm(:,3)>0,:) = [];
hist(hm(:,2))
hist(hm(:,2),100)
hist(hm(:,2),1000)
hm(hm(:,2)>0,:) = [];
clf
PlotXY(hm(:,2:3),'.');
size(rp)
load('day12.cell_metrics.cellinfo.mat')
spikesCell = cell_metrics.spikes.times';
spikes = sortrows(Group(spikesCell{:}));
nNeurons = length(cell_metrics.spikes.times);
regions = cell_metrics.brainRegion;
regionNames = unique(regions)
regionCell = zeros(length(regions),1);
for i=1:length(regionNames)
regionCell(strcmp(regions,regionNames{i})) = i;
end
hpc = Group(spikesCell{regionCell==find(strcmp(regionNames,'CA1'))});
PETH(hpc(:,1),hm(:,1))
PETH(hpc(:,1),hm0(:,1))
PETH(hpc(:,1),hm(:,1))
load('day12.deltaWaves.events.mat')
deltas = repmat(deltaWaves.peaks,1,2);
PETH(hm(:,1),deltas(:,2));
rr = interp1(tl,ripPower0,hm(:,1));
plot(hm(:,2),rr,'.');
plot(-hm(:,2),rr,'.');
plot(-hm(:,2),rr,'.','markersize',1);
figure; DensityMap(-hm(:,2),rr);
DensityMap(-hm(:,2),rr,'nBins',500);
DensityMap(-hm(:,2),rr,'nBins',5000);
ok = hm(:,2)>-5000 & rr<200;
DensityMap(-hm(ok,2),rr(ok),'nBins',500);
DensityMap(-hm(ok,2),rr(ok),'nBins',500,'smooth',0);
dbquit
rippleChannels.Ripple_Channel = 4+1; rippleChannels.Sharpwave_Channel = 10+1;
ripples = DetectSWR([rippleChannels.Ripple_Channel rippleChannels.Sharpwave_Channel],'basepath',basepath,'saveMat',true,'forceDetect',true,'check',true);
% golden ripple: @27725.208
ok = swDiffAll<1500 & ripPowerAll<300;
figure; DensityMap(swDiffAll(ok),ripPowerAll(ok),'smooth',0);
min(swDiffAll)
ok = swDiffAll<15000 & ripPowerAll<300;
DensityMap(-hm(ok,2),rr(ok),'nBins',500,'smooth',0);
DensityMap(swDiffAll(ok),ripPowerAll(ok),'smooth',0);
DensityMap(swDiffAll(ok),ripPowerAll(ok),'smooth',0,'nBins',500);
ok = swDiffAll<5000 & ripPowerAll<300 & swDiffAll>0;
DensityMap(swDiffAll(ok),ripPowerAll(ok),'smooth',0,'nBins',500);
DensityMap(swDiffAll(ok),ripPowerAll(ok),'smooth',1,'nBins',500);
num2str(t(find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim))))
t = rez.st3(rez.st3(:,2)==i)/20000;
t = featureTs/1250;
num2str(t(find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim))))
open lfp
clf
tl = (1:length(lfp))'/1250;
plot(tl,lfp(:,1));
[clean,bad,badIntervals] = CleanLFP([tl,lfp(:,1)]);
[clean,bad,badIntervals] = CleanLFP([tl,double(lfp(:,1))]);
PlotHVLines(11607.8208,'v','k--');
PlotIntervals(Restrict(badIntervals,xlim));
this = [tl,double(lfp(:,1))];
that = Restrict(this,xlim);
[w,wt,wf] = WaveletSpectrogram(that,'range',[50 625]);
figure; PlotColorMap(w,'x',wt,'y',wf);
PlotHVLines(11607.8208,'v','k--');
Portion(InIntervals(t,badIntervals))
figure; anovabar(rp,InIntervals(t,badIntervals));
rp =ripPowerAll; t = featureTs/1250; s = swDiffAll;
figure; anovabar(rp,InIntervals(t,badIntervals));
Channels
Channels = [5 11 76+1];
lfp        = LoadBinary(lfp_file,'frequency',SR,'nChannels',nChan,...
'channels',Channels);
which CleanRez
raly
which DetectSWR
rippleChannels
ripples = DetectSWR([rippleChannels.Ripple_Channel rippleChannels.Sharpwave_Channel 76+1],'basepath',basepath,'saveMat',true,'forceDetect',true,'check',true);
lfpM       = lfp - int16(repmat(mean(lfp(:,1:2),2),[1 Nchan]));
plto(lfp);
plot(lfp);
clf
plot(lfp);
lfpM       = lfp(:,1) - int16(repmat(mean(lfp(:,1),2),[1 Nchan]));
rip        = firfilt( lfpM, hRip2 );    % highpass filter
eegLo      = firfilt( rip, hRip1 );     % lowpass filter
rip        = rip - eegLo;               % difference of gaussians
lfpM       = lfp(:,1) - int16(repmat(mean(lfp(:,1),2),[1 Nchan]));
lfpM       = lfp - int16(repmat(mean(lfp(:,1:2),2),[1 Nchan]));
rip        = firfilt( lfpM, hRip2 );    % highpass filter
eegLo      = firfilt( rip, hRip1 );     % lowpass filter
rip        = rip - eegLo;               % difference of gaussians
ripWindow  = pi / mean( ripBP );
powerWin   = makegausslpfir( 1 / ripWindow, SR, 6 );
rip        = abs(rip);
ripPower0  = firfilt( rip, powerWin );
clf
plot(ripPower0(:,1:2));
% subtract mean before filtering
lfpM       = lfp - int16(repmat(mean(lfp(:,1:2),2),[1 Nchan]));
%clear lfp
rip        = firfilt( lfpM, hRip2 );    % highpass filter
clear lfpM
eegLo      = firfilt( rip, hRip1 );     % lowpass filter
% rip        = rip -.lfpLo;               % difference of gaussians
% clear.lfpLo
rip        = rip - eegLo;               % difference of gaussians
clear eegLo
ripWindow  = pi / mean( ripBP );
powerWin   = makegausslpfir( 1 / ripWindow, SR, 6 );
rip        = abs(rip);
ripPower0  = firfilt( rip, powerWin );
plot(ripPower0);
plot(ripPower0(:,1)-ripPower0(:,end));
xlim(27725.208+[-1 1]);
plot(tl,ripPower0(:,1)-ripPower0(:,end));
tl = (1:length(lfp))'/1250;
plot(tl,ripPower0(:,1)-ripPower0(:,end));
xlim(27725.208+[-1 1]);
PlotHVLines(27725.208,'v');
nChan
Nchan
rip        = firfilt( lfp(:,3), hRip2 );    % highpass filter
eegLo      = firfilt( rip, hRip1 );     % lowpass filter
rip        = rip - eegLo;               % difference of gaussians
rip        = abs(rip);
noise  = firfilt( rip, powerWin );
hold off
interval = (27725.208+[-1 1]);
in = InIntervals(tl,interval);
plot(tl(in),ripPower0(in,1),'r');
hold all
plot(tl(in),noise(in,1),'k');
rip        = firfilt( lfp(:,1), hRip2 );    % highpass filter
eegLo      = firfilt( rip, hRip1 );     % lowpass filter
rip        = rip - eegLo;               % difference of gaussians
rip        = abs(rip);
signal  = firfilt( rip, powerWin );
plot(tl(in),signal(in,1),'k');
plot(tl(in),signal(in,1),'b');
% We can't use the difference because if we do, the values are much lower than the values of the (non-substracted) noise channel,
% and subtracting signal from the far-away noise channel doesn't make sense
rip        = firfilt( lfp(:,3), hRip2 );    % highpass filter
eegLo      = firfilt( rip, hRip1 );     % lowpass filter
rip        = rip - eegLo;               % difference of gaussians
rip        = abs(rip);
noise  = firfilt( rip, powerWin );
rip        = firfilt( lfp(:,1), hRip2 );    % highpass filter
eegLo      = firfilt( rip, hRip1 );     % lowpass filter
rip        = rip - eegLo;               % difference of gaussians
rip        = abs(rip);
signal  = firfilt( rip, powerWin );
ripPower0 = signal - noise;
plot(tl(in),ripPower0(in,1),'r');
interval = (23393.208+[-1 1]);
in = InIntervals(tl,interval);
plot(tl(in),ripPower0(in,1),'r');
clf
plot(tl(in),ripPower0(in,1),'r');
ripples = DetectSWR([rippleChannels.Ripple_Channel rippleChannels.Sharpwave_Channel 76+1],'basepath',basepath,'saveMat',true,'forceDetect',true,'check',true);
dbcont
ok = swDiffAll<4000 & ripPowerAll<300 & swDiffAll>0;
figure; DensityMap(swDiffAll(ok),ripPowerAll(ok),'smooth',1,'nBins',500);
rp =ripPowerAll; t = featureTs/1250; s = swDiffAll;
num2str(t(find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim))))
hold all
i = find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim)); plot(swDiffAll(i),ripPowerAll(i),'bx','linewidth',2);
num2str(t(find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim))))
interval = 19884.4152 + [0 1];
basepath
load('day12.deltaWaves.events.mat')
deltas = repmat(deltaWaves.peaks,1,2);
i = find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim)); plot(swDiffAll(i),ripPowerAll(i),'bx','linewidth',2);
tl = (1:length(lfp))'/1250;
sw = interp1(tl,lfpLow(:,2),t);
lfpLow     = firfilt( lfp, hSw2 );      % lowpass filter
eegLo      = firfilt( lfpLow, hSw1 );   % highpass filter
% lfpLow     = lfpLow -.lfpLo;            % difference of Gaussians
% clear .lfpLo                          %DL: not sure what this '.' was doing here...
lfpLow     = lfpLow - eegLo;            % difference o
sw = interp1(tl,lfpLow(:,2),t);
sw(i)
figure; scatter(s,rp,10,sw);
clim
clim([-1000 0]);
clim([-1 0]);
figure; hist(sw(s>2000),100;
figure; hist(sw(s>2000),100);
hist(sw(s>2000),1000);
bad = sw>0;
hold on
plot(swDiffAll(i),ripPowerAll(i),'kx','linewidth',2);
num2str(t(find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim))))
i = find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim)); plot(swDiffAll(i),ripPowerAll(i),'mo','linewidth',2);
plot(swDiffAll(i),ripPowerAll(i),'yo','linewidth',2);
idx2 = idx2 & ~bad;
idx1 = idx1 & ~bad;
clf
plot(swDiffAll(idx2),ripPowerAll(idx2),'k.','markersize',1); hold on
plot(swDiffAll(idx1),ripPowerAll(idx1),'r.','markersize',1); legend('non-ripples','ripples');
clf
plot(swDiffAll(idx2),ripPowerAll(idx2),'k.','markersize',1); hold on
plot(swDiffAll(idx1),ripPowerAll(idx1),'r.','markersize',1);
num2str(t(find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim))))
hold all
plot(swDiffAll(idx1),ripPowerAll(idx1),'go','markersize',1);
plot(swDiffAll(idx1),ripPowerAll(idx1),'go','markersize',20);
i
t(i+(0:3)) - floor(t(i))
basepath
mmax(t)
t = featureTs/1250;
FindClosest(t,11018.858)
i
floor(t(i))
i=51754
t(i+(-2:2)) - floor(t(i))
num2str(t(i+2))
interval = 11018.395 + [0 1];
in = InIntervals(tl,interval);
figure; plot(tl(in),lfpLow(in,:));
PlotHVLines(Restrict(t,xlim),'v','k--')
legend('1','2','3')
which DetectSWR
ripples = DetectSWR([rippleChannels.Ripple_Channel rippleChannels.Sharpwave_Channel 76+1],'basepath',basepath,'saveMat',true,'forceDetect',true,'check',true);
% More golden ripples around 11018.8
close all
raly
ripples = DetectSWR([rippleChannels.Ripple_Channel rippleChannels.Sharpwave_Channel 76+1],'basepath',basepath,'saveMat',true,'forceDetect',true,'check',true);
dbcont
sw = interp1(tl,lfpLow(:,2),t);
tl = (1:length(lfp))'/1250;
sw = interp1(tl,lfpLow(:,2),t);
t = featureTs/1250;
rp =ripPowerAll; t = featureTs/1250; s = swDiffAll;
sw = interp1(tl,lfpLow(:,2),t);
Portion(sw>0)
bad = sw>0;
idx1 = idx1 & ~bad;
idx2 = idx2 & ~bad;
clf
plot(swDiffAll(idx2),ripPowerAll(idx2),'k.','markersize',1); hold on
plot(swDiffAll(idx1),ripPowerAll(idx1),'r.','markersize',1);
FindClosest(t,11018.85)
i==49884;
num2str(t(i))
i=49884;
num2str(t(i))
hold all
plot(swDiffAll(i),ripPowerAll(i),'yo','linewidth',2);
num2str(t(i+1))
i=i+1;
plot(swDiffAll(i),ripPowerAll(i),'yo','linewidth',2);
num2str(t(i+1))
plot(swDiffAll(i),ripPowerAll(i),'yo','linewidth',2);
i=i+1
plot(swDiffAll(i),ripPowerAll(i),'yo','linewidth',2);
plot(swDiffAll(i),ripPowerAll(i),'go','linewidth',2);
num2str(t(find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim))))
i = find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim)); plot(swDiffAll(i),ripPowerAll(i),'co','linewidth',2);
num2str(t(find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim))))
figure; plot(lfpLow(:,2));
PlotHVLines(340.18,'v');
clf
figure; plot(tl,lfpLow(:,2));
PlotHVLines(340.18,'v');
[clean,bad,badIntervals] = CleanLFP([tl,double(lfp(:,2))]);
PlotIntervals(Restrict(badIntervals,xlim));
bad = sw>0 | InIntervals(t,badIntervals);
Portion(bad)
idx1 = idx1 & ~bad;
idx2 = idx1 & ~bad;
clf
plot(swDiffAll(idx2),ripPowerAll(idx2),'k.','markersize',1); hold on
plot(swDiffAll(idx1),ripPowerAll(idx1),'r.','markersize',1);
idx2 = ~idx1 & ~bad;
clf
plot(swDiffAll(idx2),ripPowerAll(idx2),'k.','markersize',1); hold on
plot(swDiffAll(idx1),ripPowerAll(idx1),'r.','markersize',1);
find(InIntervals(b,xlim) & InIntervals(sw,ylim))
num2str(t(find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim))))
hold on
plot(swDiffAll(i),ripPowerAll(i),'co','linewidth',2);
i = find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim)); plot(swDiffAll(i),ripPowerAll(i),'co','linewidth',2);
num2str(t(find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim))))
figure; plot(tl,swDiff0);
plot(tl,swDiff);
plot(tl,Smooth(swDiff,50));
PlotHVLines(1086.78,'v');
xlim(11018.8+[-1 1]);
PlotHVLines([-1 1]*1000,'h','k-');
plot(tl,Smooth(swDiff,1250*10));
clf
plot(tl,Smooth(abs(swDiff),1250*10));
PlotHVLines(1086.78,'v');
ylim(ylim);
xlim(11018.8+[-1 1]);
PlotHVLines(500,'h','k');
PlotHVLines(11018.8,'v','r');
PlotHVLines(300,'h','k');
smoothed = Smooth(abs(swDiff),1250*10);
noisyPeriods = tl(FindInterval(sm>300));
noisyPeriods = tl(FindInterval(smoothed>300));
dur(noisyPeriods)
noisyPeriods = Consolidate(noisyPeriods,'epsilon',1);
noisyPeriods = ConsolidateIntervals(noisyPeriods,'epsilon',1);
dur(noisyPeriods)
PlotIntervals(Restrict(noisyPeriods,xlim));
PlotIntervals(Restrict(bad,xlim));
PlotIntervals(Restrict(bad,xlim),'r');
PlotIntervals(Restrict(bad,xlim),'color','r'');
PlotIntervals(Restrict(bad,xlim),'color','r');
PlotIntervals(Restrict(badIntervals,xlim),'color','r');
bad = sw>0 | InIntervals(t,badIntervals) | InIntervals(t,noisyPeriods);
idx1 = idx1 & ~bad;
idx2 = idx2 & ~bad;
clf
plot(swDiffAll(idx2),ripPowerAll(idx2),'k.','markersize',1); hold on
plot(swDiffAll(idx1),ripPowerAll(idx1),'r.','markersize',1); legend('non-ripples','ripples');
num2str(t(find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim))))
hold all
i = find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim)); plot(swDiffAll(i),ripPowerAll(i),'yo','linewidth',2);
num2str(t(find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim))))
i = find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim)); plot(swDiffAll(i),ripPowerAll(i),'yo','linewidth',2);
num2str(t(find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim))))
i = find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim)); plot(swDiffAll(i),ripPowerAll(i),'yo','linewidth',2);
num2str(t(find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim))))
i = find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim)); plot(swDiffAll(i),ripPowerAll(i),'go','linewidth',2);
num2str(t(find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim))))
i = find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim)); plot(swDiffAll(i),ripPowerAll(i),'go','linewidth',2);
i+1
i=i+1
num2str(t(i))
plot(swDiffAll(i),ripPowerAll(i),'yo','linewidth',2);
num2str(t(find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim))))
i = find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim)); plot(swDiffAll(i),ripPowerAll(i),'bx','linewidth',2);
legend('off');
num2str(t(find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim))))
i = find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim)); plot(swDiffAll(i),ripPowerAll(i),'go','linewidth',2);
num2str(t(find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim))))
i = find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim)); plot(swDiffAll(i),ripPowerAll(i),'go','linewidth',2);
i=i-1
num2str(t(i))
plot(swDiffAll(i),ripPowerAll(i),'go','linewidth',2);
i=i-1
num2str(t(i))
plot(swDiffAll(i),ripPowerAll(i),'yo','linewidth',2);
num2str(t(i))
num2str(t(find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim))))
i = find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim)); plot(swDiffAll(i),ripPowerAll(i),'co','linewidth',2);
i=i-1
num2str(t(i))
plot(swDiffAll(i),ripPowerAll(i),'yo','linewidth',2);
num2str(t(i))
num2str(t(find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim))))
i = find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim)); plot(swDiffAll(i),ripPowerAll(i),'bx','linewidth',2);
num2str(t(find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim))))
i = find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim)); plot(swDiffAll(i),ripPowerAll(i),'bx','linewidth',2);
i=i+1
num2str(t(i))
plot(swDiffAll(i),ripPowerAll(i),'yo','linewidth',2);
num2str(t(find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim))))
i = find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim)); plot(swDiffAll(i),ripPowerAll(i),'bx','linewidth',2);
num2str(t(find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim))))
i = find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim)); plot(swDiffAll(i),ripPowerAll(i),'co','linewidth',2);
num2str(t(find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim))))
i = find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim)); plot(swDiffAll(i),ripPowerAll(i),'bx','linewidth',2);
num2str(t(find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim))))
i = find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim)); plot(swDiffAll(i),ripPowerAll(i),'co','linewidth',2);
num2str(t(find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim))))
i = find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim)); plot(swDiffAll(i),ripPowerAll(i),'bx','linewidth',2);
num2str(t(find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim))))
i = find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim)); plot(swDiffAll(i),ripPowerAll(i),'co','linewidth',2);
num2str(t(find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim))))
i = find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim)); plot(swDiffAll(i),ripPowerAll(i),'co','linewidth',2);
i = find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim)); plot(swDiffAll(i),ripPowerAll(i),'bx','linewidth',2);
i = find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim)); plot(swDiffAll(i),ripPowerAll(i),'wo','linewidth',2);
i = find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim)); plot(swDiffAll(i),ripPowerAll(i),'bx','linewidth',2);
num2str(t(find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim))))
i = find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim)); plot(swDiffAll(i),ripPowerAll(i),'go','linewidth',2);
i=i+1
num2str(t(i))
plot(swDiffAll(i),ripPowerAll(i),'yo','linewidth',2);
plot(swDiffAll(i),ripPowerAll(i),'ro','linewidth',2);
plot(swDiffAll(i),ripPowerAll(i),'yo','linewidth',2);
i=i+1
num2str(t(i))
plot(swDiffAll(i),ripPowerAll(i),'yo','linewidth',2);
i=i+1
num2str(t(i))
plot(swDiffAll(i),ripPowerAll(i),'yo','linewidth',2);
num2str(t(i))
i = find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim)); plot(swDiffAll(i),ripPowerAll(i),'yo','linewidth',2);
num2str(t(i))
i = find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim)); plot(swDiffAll(i),ripPowerAll(i),'yo','linewidth',2);
num2str(t(i))
i = find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim)); plot(swDiffAll(i),ripPowerAll(i),'go','linewidth',2);
num2str(t(i))
num2str(t(find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim))))
i = find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim)); plot(swDiffAll(i),ripPowerAll(i),'go','linewidth',2);
num2str(t(find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim))))
i = find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim)); plot(swDiffAll(i),ripPowerAll(i),'yo','linewidth',2);
i=i+1
num2str(t(i))
plot(swDiffAll(i),ripPowerAll(i),'yo','linewidth',2);
i=i+1
num2str(t(i))
plot(swDiffAll(i),ripPowerAll(i),'co','linewidth',2);
i=i+1
num2str(t(i))
plot(swDiffAll(i),ripPowerAll(i),'yo','linewidth',2);
num2str(t(i))
num2str(t(find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim))))
i = find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim)); plot(swDiffAll(i),ripPowerAll(i),'co','linewidth',2);
i=i-1
num2str(t(i))
i=i-1
num2str(t(i))
num2str(t(find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim))))
i = find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim)); plot(swDiffAll(i),ripPowerAll(i),'co','linewidth',2);
num2str(t(find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim))))
i = find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim)); plot(swDiffAll(i),ripPowerAll(i),'co','linewidth',2);
num2str(t(find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim))))
i = find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim)); plot(swDiffAll(i),ripPowerAll(i),'go','linewidth',2);
num2str(t(find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim))))
i = find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim)); plot(swDiffAll(i),ripPowerAll(i),'bx','linewidth',2);
num2str(t(find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim))))
i = find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim)); plot(swDiffAll(i),ripPowerAll(i),'bx','linewidth',2);
num2str(t(find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim))))
i = find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim)); plot(swDiffAll(i),ripPowerAll(i),'co','linewidth',2);
num2str(t(find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim))))
i = find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim)); plot(swDiffAll(i),ripPowerAll(i),'bx','linewidth',2);
num2str(t(find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim))))
i = find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim)); plot(swDiffAll(i),ripPowerAll(i),'gx','linewidth',2);
i = find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim)); plot(swDiffAll(i),ripPowerAll(i),'wx','linewidth',2);
i = find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim)); plot(swDiffAll(i),ripPowerAll(i),'go','linewidth',2);
num2str(t(find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim))))
i = find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim)); plot(swDiffAll(i),ripPowerAll(i),'co','linewidth',2);
num2str(t(find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim))))
i = find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim)); plot(swDiffAll(i),ripPowerAll(i),'co','linewidth',2);
num2str(t(find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim))))
i = find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim)); plot(swDiffAll(i),ripPowerAll(i),'go','linewidth',2);
save('temp_ripple_detection_OJR42_day12.mat','idx1','idx2','swDiffAll','ripPowerAll');
num2str(t(find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim))))
plot(swDiffAll(i),ripPowerAll(i),'yo','linewidth',2);
plot(swDiffAll(i),ripPowerAll(i),'go','linewidth',2);
i = find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim)); plot(swDiffAll(i),ripPowerAll(i),'yo','linewidth',2);
num2str(t(find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim))))
i = find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim)); plot(swDiffAll(i),ripPowerAll(i),'yo','linewidth',2);
num2str(t(find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim))))
i = find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim)); plot(swDiffAll(i),ripPowerAll(i),'yo','linewidth',2);
num2str(t(find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim))))
i = find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim)); plot(swDiffAll(i),ripPowerAll(i),'yo','linewidth',2);
num2str(t(find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim))))
i = find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim)); plot(swDiffAll(i),ripPowerAll(i),'go','linewidth',2);
num2str(t(find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim))))
i = find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim)); plot(swDiffAll(i),ripPowerAll(i),'yo','linewidth',2);
num2str(t(find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim))))
i = find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim)); plot(swDiffAll(i),ripPowerAll(i),'yo','linewidth',2);
num2str(t(find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim))))
i = find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim)); plot(swDiffAll(i),ripPowerAll(i),'yo','linewidth',2);
i = i+1
num2str(t(i))
plot(swDiffAll(i),ripPowerAll(i),'yo','linewidth',2);
num2str(t(i))
i = find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim)); plot(swDiffAll(i),ripPowerAll(i),'yo','linewidth',2);
num2str(t(i))
plot(swDiffAll(i),ripPowerAll(i),'go','linewidth',2);
num2str(t(i))
num2str(t(find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim))))
i = find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim)); plot(swDiffAll(i),ripPowerAll(i),'yo','linewidth',2);
num2str(t(find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim))))
i = find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim)); plot(swDiffAll(i),ripPowerAll(i),'go','linewidth',2);
num2str(t(find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim))))
i = find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim)); plot(swDiffAll(i),ripPowerAll(i),'co','linewidth',2);
num2str(t(find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim))))
i = find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim)); plot(swDiffAll(i),ripPowerAll(i),'co','linewidth',2);
num2str(t(find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim))))
i = find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim)); plot(swDiffAll(i),ripPowerAll(i),'yo','linewidth',2);
num2str(t(find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim))))
i = find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim)); plot(swDiffAll(i),ripPowerAll(i),'go','linewidth',2);
num2str(t(find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim))))
i = find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim)); plot(swDiffAll(i),ripPowerAll(i),'yo','linewidth',2);
i = i+1
num2str(t(i))
num2str(t(find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim))))
i = find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim)); plot(swDiffAll(i),ripPowerAll(i),'co','linewidth',2);
num2str(t(find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim))))
i = find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim)); plot(swDiffAll(i),ripPowerAll(i),'go','linewidth',2);
num2str(t(find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim))))
i = find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim)); plot(swDiffAll(i),ripPowerAll(i),'bx','linewidth',2);
num2str(t(find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim))))
i = find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim)); plot(swDiffAll(i),ripPowerAll(i),'go','linewidth',2);
num2str(t(find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim))))
i = find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim)); plot(swDiffAll(i),ripPowerAll(i),'bx','linewidth',2);
num2str(t(find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim))))
i = find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim)); plot(swDiffAll(i),ripPowerAll(i),'co','linewidth',2);
num2str(t(find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim))))
i = find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim)); plot(swDiffAll(i),ripPowerAll(i),'bx','linewidth',2);
num2str(t(find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim))))
i = find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim)); plot(swDiffAll(i),ripPowerAll(i),'co','linewidth',2);
num2str(t(find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim))))
i = find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim)); plot(swDiffAll(i),ripPowerAll(i),'yo','linewidth',2);
num2str(t(find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim))))
i = find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim)); plot(swDiffAll(i),ripPowerAll(i),'bx','linewidth',2);
num2str(t(find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim))))
i = find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim)); plot(swDiffAll(i),ripPowerAll(i),'bx','linewidth',2);
num2str(t(find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim))))
i = find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim)); plot(swDiffAll(i),ripPowerAll(i),'go','linewidth',2);
num2str(t(find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim))))
i = find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim)); plot(swDiffAll(i),ripPowerAll(i),'bx','linewidth',2);
num2str(t(find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim))))
i = find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim)); plot(swDiffAll(i),ripPowerAll(i),'go','linewidth',2);
q = UIInPolygon(swDiffAll,ripPowerAll);
idx1 = q & ~bad;
open q
idx2 = ~q & ~bad;
figure
clf
plot(swDiffAll(idx2),ripPowerAll(idx2),'k.','markersize',1); hold on
plot(swDiffAll(idx1),ripPowerAll(idx1),'r.','markersize',1); legend('non-ripples','ripples');
num2str(t(find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim))))
i = FindClosest(t,9273.1);
num2str(t(i))
plot(swDiffAll(i),ripPowerAll(i),'yo','linewidth',2);
hold on
plot(swDiffAll(i),ripPowerAll(i),'yo','linewidth',2);
num2str(t(i))
i = find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim)); plot(swDiffAll(i),ripPowerAll(i),'yo','linewidth',2);
num2str(t(i))
num2str(t(find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim))))
hold on; i = find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim)); plot(swDiffAll(i),ripPowerAll(i),'co','linewidth',2);
num2str(t(find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim))))
i = find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim)); plot(swDiffAll(i),ripPowerAll(i),'go','linewidth',2);
num2str(t(find(InIntervals(ripPowerAll,ylim) & InIntervals(swDiffAll,xlim))))
q = UIInPolygon(swDiffAll,ripPowerAll);
idx1 = q & ~bad;
idx2 = ~q & ~bad;
figure
clf
plot(swDiffAll(idx2),ripPowerAll(idx2),'k.','markersize',1); hold on
plot(swDiffAll(idx1),ripPowerAll(idx1),'r.','markersize',1); legend('non-ripples','ripples');
save('temp_ripple_detection_OJR42_day12.mat','idx1','idx2','swDiffAll','ripPowerAll');
dbcont
2
cd basepath
basepath
cd(basepath)
r = ripples.timestamps(:,1:2);
PETH(deltas(:,2),r(:,1));
figure; PETH(deltas(:,2),r(:,1));
PETH(deltas(:,2),r(:,1),'nBins',501);
PETH(hpc(:,1),r(:,1));
close all
pre = r(:,1)<MergePoints.timestamps(4,1);
post = r(:,1)>MergePoints.timestamps(6,1);
anovabar(diff(r,[],2),[pre post]);
anovabar(diff(r,[],2),[-pre + post]);
anovabar(diff(r,[],2)>0.1,[-pre + post]);
regionNames
pfc = Group(spikesCell{regionCell~=find(strcmp(regionNames,'CA1'))});
PETH(pfc(:,1),hpc(:,1))
PETH(pfc(:,1),r(:,1))
hold all
long = r(:,end)-r(:,1)>0.1;
PETH(pfc(:,1),r(log,1))
PETH(pfc(:,1),r(long,1))
PETH(pfc(:,1),r(~long,1))
clf
PETH(deltas(:,1),r(~long,1));
hold all
PETH(deltas(:,1),r(long,1));
open deltas
clf
PETH(deltas(:,2),r(~long,1));
PETH(deltas(:,2),r(~long,1),'nBins',501);
PETH(deltas(:,2),r(~long,1),'nBins',501); PlotHVLines(0,'v','k--');
hold all
PETH(deltas(:,2),r(long,1),'nBins',501);
PETH(deltas(:,2),r(long,1),'nBins',201);
clf
hist(deltas(:,2))
hist(deltas(:,2),100)
subplot(2,1,1);
hist(deltas(:,2),100)
subplot(2,1,2);
subplot(3,1,1);
hist(deltas(:,2),100)
subplot(3,1,2);
hist(r(~long,1),100)
subplot(3,1,3);
hist(r(long,1),100)
PlotHVLines(MergePoints.timestamps(:,1),'v','r');
load('day12.SleepState.states.mat')
PlotIntervals(SleepState.ints.NREMstate,'color','k');
sws = SleepState.ints.NREMstate;
figure; PETH(deltas(:,2),Restrict(r(~long,1),sws),'nBins',201);
hold all
PETH(deltas(:,2),Restrict(r(long,1),sws),'nBins',201);
clf
PETH(pfc(:,1),Restrict(r(~long,1),sws),'nBins',201);
hold all
PETH(pfc(:,1),Restrict(r(long,1),sws),'nBins',201);
r = [ripples.timestamps(:,1) ripples.peaks ripples.timestamps(:,2) ripples.peakNormedPower];
open ripples
r = [ripples.timestamps(:,1) ripples.peaks ripples.timestamps(:,2) ripples.ripMax ripples.SwMax];
r = [ripples.timestamps(:,1) ripples.peaks ripples.timestamps(:,2) ripples.RipMax ripples.SwMax];
open r
lfpstruct = getLFP(9+1);
lfp = [lfpstruct.timestamps]; lfp(:,2) = lfpstruct.data;
PETH(lfp,deltas(:,2));
SaveCustomEvents('ripples.rip.evt',r(:,1:3),{'ripple start','ripple peak','ripple stop'});
open lfp
PETH(lfp,r(:,2));
PETH(lfp,r(:,2),'durations',[-1 1]*0.5);
PETH(lfp,deltas(:,2),'durations',[-1 1]*0.5);
[q,qt] = PETH(lfp,deltas(:,2),'durations',[-1 1]*0.5);
PlotColorMap(sortby(q,deltas(:,5)),'x',qt);
clim
PlotColorMap(sortby(q,deltas(:,5)),'x',qt); clim(clim/2)
PlotColorMap(sortby(q,deltas(:,5)),'x',qt); clim(clim/5)
PlotColorMap(sortby(q,deltas(:,5)),'x',qt); clim(clim/10)
PlotColorMap(sortby(deltas(:,5),deltas(:,5)),'x',qt); clim(clim/10)
PlotColorMap(sortby(deltas(:,5),deltas(:,5)),'x',qt);
PlotColorMap(sortby(q,deltas(:,5)),'x',qt); clim(clim/10)
PlotColorMap(Shrink(sortby(q,deltas(:,5)),100,1),'x',qt); clim(clim/10)
figure; PETH(r(:,1),deltas(:,2),'durations',[-1 1]*0.5);
PETH(r(:,2),deltas(:,2),'durations',[-1 1]*0.5);
PETH(r(:,2),deltas(:,2),'durations',[-1 1]*0.2);
PETH(r(:,2),deltas(:,2),'durations',[-1 1]*0.5);
PETH(r(:,2),deltas(:,2),'durations',[-1 1]*0.5,'smooth',0);
PETH(r(long,2),deltas(:,2),'durations',[-1 1]*0.5,'smooth',0);
PETH(r(long,2),deltas(:,2),'durations',[-1 1],'smooth',0);
PETH(r(long,2),deltas(:,2),'durations',[-1 1],'smooth',1);
PETH(deltas(:,2),Restrict(r(long,2),sws),'durations',[-1 1],'smooth',1);
hold all
PETH(deltas(:,2),Restrict(r(~long,2),sws),'durations',[-1 1],'smooth',1);
followed = InIntervals(r(:,2),[deltas(:,2)-0.2 deltas(:,2)]);
Portion(followed)
clf
PETH(deltas(:,2),Restrict(r(:,2),sws),'durations',[-1 1],'smooth',1);
[q,qt] = PETH(deltas(:,2),Restrict(r(:,2),sws),'durations',[-1 1],'smooth',1);
semplot(qt,q(followed,:))
[q,qt] = PETH(deltas(:,2),r(:,2),'durations',[-1 1],'smooth',1);
in = InIntervals(r(:,2),sws);
semplot(qt,q(followed & in,:))
semplot(qt,q(~followed & in,:))
clf
anovavar(long(in),followed)
anovabar(long(in),followed(in))
anovabar(followed(in),long(in))
clf
[q,qt] = PETH(pfc(:,1),r(:,2),'durations',[-1 1],'smooth',1);
semplot(qt,q(followed & in,:))
semplot(qt,q(~followed & in,:))
semplot(qt,q(~followed & in & long,:))
clf
semplot(qt,q(~followed & in & long,:))
semplot(qt,q(~followed & in & ~long,:))
max(spikes(:,2))
clf
[q,qt] = PETH(hpc(:,1),r(:,2),'durations',[-1 1],'smooth',1);
semplot(qt,q(~followed & in & ~long,:))
semplot(qt,q(~followed & in & long,:))
semplot(qt,q(followed & in & long,:))
clf
semplot(qt,q(followed & in & long,:),'r')
semplot(qt,q(~followed & in & long,:),'b')
clf
semplot(qt,q(followed & in & ~long,:),'r')
semplot(qt,q(~followed & in & ~long,:),'b')
clf
semplot(qt,q(followed & in & long,:),'r')
semplot(qt,q(~followed & in & long,:),'b')
interval = [r(:,2)-0.05 r(:,2)+0.05];
empty = CountInIntervals(hpc(:,1),interval)==0;
sum(empty)
Portion(empty)
clf
semplot(qt,q(~followed & in & long & ~empty,:),'b')
semplot(qt,q(followed & in & long & ~empty,:),'r')
semplot(qt,q(followed & in & long & empty,:),'r')
clf
interval = [r(:,2)-0.1 r(:,2)+0.1];
empty = CountInIntervals(hpc(:,1),interval)==0;
semplot(qt,q(empty,:),'r')
semplot(qt,q(~empty,:),'r')
clf
find(empty,1)
r(1,2)
CountInIntervals(hpc(:,1),r(1,2)+[-1 1]*0.05)
CountInIntervals(hpc(:,1),r(1,2)+[-1 1]*0.1)
CountInIntervals(hpc(:,1),r(1,2)+[-1 1]*0.2)
plot(qt,q(1,:))
findmin(abs(qt))
clf
plot(q(:,51))
plot(sortby(q(:,51),empty))
find(empty & q(:,51)>0,1)
i=6;
plot(q(i,:))
plot(qt,q(i,:))
open r
diff(r(i,1:3))
CountInIntervals(hpc(:,1),r(i,2)+[-1 1]*0.2)
CountInIntervals(hpc(:,1),r(i,2)+[-1 1]*0.1)
[q,qt] = PETH(hpc(:,1),r(:,2),'durations',[-1 1],'smooth',0);
plot(qt,q(i,:))
open hpc
hpc = sortrows(Group(spikesCell{regionCell==find(strcmp(regionNames,'CA1'))}));
spikes = sortrows(Group(spikesCell{:}));
pfc = sortrows(Group(spikesCell{regionCell~=find(strcmp(regionNames,'CA1'))}));
empty = CountInIntervals(hpc(:,1),[r(:,2)-0.1 r(:,2)+0.1])==0;
Portion(empty)
empty = CountInIntervals(hpc(:,1),[r(:,2)-0.05 r(:,2)+0.05])==0;
Portion(empty)
clf
semplot(qt,q(~empty,:),'r')
semplot(qt,q(empty,:),'r')
clf
semplot(qt,q(in & long & ~empty,:),'r')
semplot(qt,q(in & ~long & ~empty,:),'b')
SaveCustomEvents('longRipples.ril.evt',r(long,1:3),{'ripple start','ripple peak','ripple stop'});
% golden ripple @ 18123.1 missed!
Ripple_Channel
rippleChannels.Ripple_Channel
rippleChannels.Sharpwave_Channel
ripples = DetectSWR([rippleChannels.Ripple_Channel rippleChannels.Sharpwave_Channel-1 76+1],'basepath',basepath,'saveMat',true,'forceDetect',true,'check',true);
which DetectSWR
raly
which DetectSWR
ripples = DetectSWR([rippleChannels.Ripple_Channel rippleChannels.Sharpwave_Channel-1 76+1],'basepath',basepath,'saveMat',true,'forceDetect',true,'check',true);
dbcont
[clean,bad,badIntervals] = CleanLFP([tl,double(lfp(:,2))]);
tl = (1:length(lfp))'/1250;
[clean,bad,badIntervals] = CleanLFP([tl,double(lfp(:,2))]);
rp =ripPowerAll; t = featureTs/1250; s = swDiffAll;
bad = InIntervals(t,badIntervals);
sum(bad)
hold on
plot(swDiffAll(bad),ripPowerAll(bad),'rx','markersize',1); legend('non-ripples','ripples');
idx1 = idx1 & ~bad;
idx2 = idx2 & ~bad;
clf
plot(swDiffAll(idx2),ripPowerAll(idx2),'k.','markersize',1); hold on
plot(swDiffAll(idx1),ripPowerAll(idx1),'r.','markersize',1);
[x,y,button] = ginput(1);
x
y
i = findmin(abs(x-swDiffAll)+abs(y-ripPowerAll))
[x,y,button] = ginput(1);
i = findmin(abs(x-swDiffAll)+abs(y-ripPowerAll));
t = featureTs/1250;
figure(3); clf; interval = t + [-1 1]*0.5; in = InIntervals(t,interval); plot(t(in),lfp(in,:));
sum(in)
interval
figure(2)
clf
plot(swDiffAll(idx2),ripPowerAll(idx2),'k.','markersize',1); hold on
plot(swDiffAll(idx1),ripPowerAll(idx1),'r.','markersize',1);
while true
figure(2);
[x,y,button] = ginput(1);
i = findmin(abs(x-swDiffAll)+abs(y-ripPowerAll));
t = featureTs/1250;
figure(3); clf; interval = t(i) + [-1 1]*0.5; in = InIntervals(t,interval); plot(t(in),lfp(in,:));
end
figure(3); clf; interval = t(i) + [-1 1]*0.5; in = InIntervals(t,interval); plot(t(in) - t(i),lfp(in,:)); xlabel(num2str(t(i)));
figure(3); clf; interval = t(i) + [-1 1]*0.5; in = InIntervals(tl,interval); plot(tl(in) - t(i),lfp(in,:)); xlabel(num2str(t(i)));
edit Ymaze
v1
which KbCheck
isnumber('1)
isnumber('1')
isdouble('1')
double('1')\
double('1')
double('0')
str2double('0')
str2double('1')
str2double('0.1325')
str2double(' ')
colors = Bright(1000,3);
colors = Bright(1000);
open colors
edit nothing
45
edit nothing
nothing
1
0.9
0.1
findmin(abs(score)-linspace(0,1,1000)'
findmin(abs(score)-linspace(0,1,1000)')
score
nothing
0.1
0
1
0.1
0.5
0.01
0.001
1
nothing
figure(2)
clf
plot(swDiffAll(idx2),ripPowerAll(idx2),'k.','markersize',1); hold on
plot(swDiffAll(idx1),ripPowerAll(idx1),'r.','markersize',1);
nothing
0
smoothed = Smooth(abs(swDiff),1250*10);
noisyPeriods = tl(FindInterval(smoothed>300));
figure; plot(smoothed))
figure; plot(smoothed)
PlotHVLines(300,'h');
noisyPeriods = ConsolidateIntervals(tl(FindInterval(smoothed>200)),'epsilon',1);
PlotIntervals(noisyPeriods)
clf
plot(tl,smoothed);
PlotIntervals(noisyPeriods)
bad = InIntervals(t,badIntervals) | InIntervals(t,noisyPeriods);
idx2 = idx2 & ~bad;
idx1 = idx1 & ~bad;
figure(2)
clf
plot(swDiffAll(idx2),ripPowerAll(idx2),'k.','markersize',1); hold on
plot(swDiffAll(idx1),ripPowerAll(idx1),'r.','markersize',1);
nothing
sw = interp1(tl,lfpLow(:,2),t);
bad = InIntervals(t,badIntervals) | InIntervals(t,noisyPeriods) | sw>0;
figure(2)
clf
plot(swDiffAll(idx2),ripPowerAll(idx2),'k.','markersize',1); hold on
plot(swDiffAll(idx1),ripPowerAll(idx1),'r.','markersize',1);
nothing
0
1
0.8
1
0
1
0.6
0.7
0
1
0
1
0
0.5
0.2
0.6
0.7
0.8
0.2
0.9
0.5
0.7
1
0
0.6
0.8
1
0.5
1
0.1
0.8
1
0.1
0.8
0.6
1
0.8
0.4
0.1
0
0.8
0.5
0.8
0.1
0.7
0.1
nothing
0.1
1
0.1
1
0.6
0.2
0.8
0.7
0.75
0.6
save('temp_ripple_detection_OJR42_day12.mat','idx1','idx2','swDiffAll','ripPowerAll','noisyPeriods','badPeriods','sw');
save('temp_ripple_detection_OJR42_day12_channel9.mat','idx1','idx2','swDiffAll','ripPowerAll','noisyPeriods','badIntervals','sw');
q = UIInPolygon(swDiffAll,ripPowerAll);
figure;
clf
plot(swDiffAll(idx2),ripPowerAll(idx2),'k.','markersize',1); hold on
plot(swDiffAll(idx1),ripPowerAll(idx1),'r.','markersize',1);
idx1 = q & ~bad;
idx2 = ~q & ~bad;
clf
plot(swDiffAll(idx2),ripPowerAll(idx2),'k.','markersize',1); hold on
plot(swDiffAll(idx1),ripPowerAll(idx1),'r.','markersize',1);
save('temp_ripple_detection_OJR42_day12_channel9.mat','idx1','idx2','swDiffAll','ripPowerAll','noisyPeriods','badIntervals','sw');
dbcont
close all
r0 = r;
r = [ripples.timestamps(:,1) ripples.peaks ripples.timestamps(:,2) ripples.RipMax ripples.SwMax];
PETH(deltas(:,2),Restrict(r(:,2),sws),'durations',[-1 1],'smooth',1);
Dist(1000,r(:,1),r0(:,1))
PETH(r(:,1),r0(:,1))
PETH(deltas(:,2),Restrict(r(:,2),sws),'durations',[-1 1],'smooth',1);
hold all
PETH(deltas(:,2),Restrict(r0(:,2),sws),'durations',[-1 1],'smooth',1);
clf
PETH(deltas(:,2),Restrict(r(:,2),sws),'durations',[-1 1],'smooth',1);
PETH(deltas(:,2),Restrict(r(:,1),sws),'durations',[-1 1],'smooth',1);
PETH(pfc(:,1),Restrict(r(:,1),sws),'durations',[-1 1],'smooth',1);
followed = InIntervals(r(:,2),[deltas(:,2)-0.2 deltas(:,2)]);
followed = InIntervals(r(:,1),[deltas(:,2)-0.2 deltas(:,2)]);
PETH(pfc(:,1),Restrict(r(~followed,1),sws),'durations',[-1 1],'smooth',1);
PETH(pfc(:,1),Restrict(r(followed,1),sws),'durations',[-1 1],'smooth',1);
hold all
PETH(pfc(:,1),Restrict(r(~followed,1),sws),'durations',[-1 1],'smooth',1);
long = r(:,end)-r(:,1)>0.1;
pre = r(:,1)<MergePoints.timestamps(4,1);
post = r(:,1)>MergePoints.timestamps(6,1);
clf
anovabar(long,[-pre +post]);
size(pre)
size(post)
size(long)
sum(pre)
sum(post)
sum(long)
long = r(:,3)-r(:,1)>0.1;
anovabar(long,[-pre +post]);
open r
Portion(long)
Portion(r0(:,3)-r0(:,1))
Portion(r0(:,3)-r0(:,1)>0.1)
Portion(r(:,3)-r(:,1)>0.1)
sum(r(:,3)-r(:,1)>0.1)
sum(r0(:,3)-r0(:,1)>0.1)
PETH(pfc(:,1),Restrict(r(:,1),sws)))
PETH(pfc(:,1),Restrict(r(:,1),sws))
PETH(pfc(:,1),Restrict(r0(:,1),sws))
PETH(pfc(:,1),Restrict(r(~followed,1),sws))
[q,qt] = PETH(pfc(:,1),Restrict(r(:,1),sws));
PlotColorMap(Shrink(sortby(q,diff(r(:,[1 3]),[],2)),100,1),'x',qt); clim(clim/10)
size(r)
size(q)
in = InIntervals(r(:,2),sws);
PlotColorMap(Shrink(sortby(q,diff(r(in,[1 3]),[],2)),100,1),'x',qt);
sum(in)
size(q)
in = InIntervals(r(:,1),sws);
PlotColorMap(Shrink(sortby(q,diff(r(in,[1 3]),[],2)),100,1),'x',qt);
PlotColorMap(Shrink(sortby(q,diff(r(in,[1 3]),[],2)),10,1),'x',qt);
PlotColorMap(Smooth(sortby(q,diff(r(in,[1 3]),[],2)),[10,1]),'x',qt);
PlotColorMap(Smooth(sortby(q,diff(r(in,[1 3]),[],2)),[50,1]),'x',qt);
PlotColorMap(Smooth(sortby(q,diff(r(in,[1 3]),[],2)),[20,1]),'x',qt);
[q,qt] = PETH(pfc(:,1),Restrict(r(:,2),sws));
PlotColorMap(Smooth(sortby(q,diff(r(in,[1 3]),[],2)),[20,1]),'x',qt);
[q,qt] = PETH(pfc(:,1),Restrict(r(:,2),sws),'durations',[-1 1]*0.5);
PlotColorMap(Smooth(sortby(q,diff(r(in,[1 3]),[],2)),[20,1]),'x',qt);
PlotColorMap(Smooth(sortby(q,diff(r(in,[1 3]),[],2)),[10,1]),'x',qt);
PlotColorMap(Shrink(sortby(q,diff(r(in,[1 3]),[],2)),10,1),'x',qt);
PlotColorMap(Shrink(sortby(q,diff(r(in,[1 3]),[],2)),20,1),'x',qt);
PlotColorMap(Shrink(sortby(q,diff(r(in,[1 3]),[],2)),50,1),'x',qt);
PlotColorMap(Shrink(sortby(q,diff(r(in,[1 3]),[],2)),50,1),'x',qt); clim(clim/2)
PlotColorMap(Shrink(sortby(q,diff(r(in,[1 3]),[],2)),50,1),'x',qt);
q(:) = zscore(q(:));
PlotColorMap(Shrink(sortby(q,diff(r(in,[1 3]),[],2)),50,1),'x',qt);
clim
clim(clim/2);
max(pfc(:,2))
for i=1:38
[q,qt] = PETH(pfc(pfc(:,2)==i,1),Restrict(r(:,2),sws),'durations',[-1 1]*0.5);
peths{i,1} = q;
end
open pethds
open peths
for i=1:max(pfc(:,2))
[q,qt] = PETH(pfc(pfc(:,2)==i,1),Restrict(r(:,2),sws),'durations',[-1 1]*0.5);
peths{i,1} = q;
q(:) = zscore(q(:));
peths{i,2} = q;
q = nanmean(q);
peths{i,3} = q;
end
q = cell2mat(peths(:,3));
PlotColorMap(q);
PlotColorMap(q,'x',qt);
restricted = Restrict(r(:,2),sws);
q = cell2mat(cellfun(@(x) nanmean(x(restricted>0.8,:)),'uniformoutput',0));
q = cell2mat(cellfun(@(x) nanmean(x(restricted>0.8,:)),peths(:,2),'uniformoutput',0));
PlotColorMap(q,'x',qt);
q = cell2mat(cellfun(@(x) nanmean(x(restricted(:,3)-restricted(:,1)>0.8,:)),peths(:,2),'uniformoutput',0));
restricted = Restrict([r(:,2) r(:,3)-r(:,1)],sws);
q = cell2mat(cellfun(@(x) nanmean(x(restricted(:,2)>0.8,:)),peths(:,2),'uniformoutput',0));
PlotColorMap(q,'x',qt);
open q
sum(restricted(:,2)>0.8)
q = cell2mat(cellfun(@(x) nanmean(x(restricted(:,2)>0.08,:)),peths(:,2),'uniformoutput',0));
PlotColorMap(q,'x',qt);
subplot(2,2,1);
q = cell2mat(cellfun(@(x) nanmean(x(restricted(:,2)>0.08,:)),peths(:,2),'uniformoutput',0));
PlotColorMap(q,'x',qt);
clim
subplot(2,2,2);
q = cell2mat(cellfun(@(x) nanmean(x(restricted(:,2)<0.08,:)),peths(:,2),'uniformoutput',0));
PlotColorMap(q,'x',qt);
clim
clims
clims([-1 1]*0.2)
clf
hsit(restricted(:,2))
hist(restricted(:,2))
hist(restricted(:,2),100)
subplot(2,2,1);
q = cell2mat(cellfun(@(x) nanmean(x(restricted(:,2)>=0.06,:)),peths(:,2),'uniformoutput',0));
PlotColorMap(q,'x',qt);
subplot(2,2,2);
q = cell2mat(cellfun(@(x) nanmean(x(restricted(:,2)<0.06,:)),peths(:,2),'uniformoutput',0));
PlotColorMap(q,'x',qt);
clims
clim
v1
Ymaze('nleftperblock',4,'screen',3,
raly
handle = figure(1);
savefig(handle,'test.pdf','-append');
help exportgraphics
exportgraphics(handle,'test.pdf','append',true);
xlabel('23');
exportgraphics('test.pdf','append',true);
exportgraphics(gcf,'test.pdf','append',true);
sum(r(:,3)-r(:,1)>0.1)
figure; hist(r(:,3)-r(:,1),100);
hist(r(:,3)-r(:,1),100);
hist(r(:,3)-r(:,1),200);
hist(r(:,3)-r(:,1),250);
hist(r(:,3)-r(:,1),80);
hist(r(:,3)-r(:,1),70);
hist(r(:,3)-r(:,1),\60);
hist(r(:,3)-r(:,1),60);
d = diff(r(:,[1 3]),[],2);
figure; Dist(1000,[r(:,3)-r(:,1) post+1],'grouped');
[h,ht] =  Dist(1000,[r(:,3)-r(:,1) post+1],'grouped');
plot(ht,Smooth(h,[2 0]));
plot(ht,Smooth(h,[5 0]));
figure; anovabar(d,[-pre post]);
d = diff(r(:,[1 3]),[],2);
size(d)
size(post)
size(pre)
figure; anovabar(d,[-pre +post]);
anovabar(d,[+post]);
raly
which Smooth
which CleanLFP
which ConsolidateIntervals
cd C:\Users\Cornell\Documents\GitHub\neurocode\utilities
edit Smooth
which Bright
which UIInPolygon
cd D
cd D:
which UIInPolygon
function in = UIInPolygon(X,Y)
%UIInPolygon - Find points in interactively defined polygon zone.
%
%  Mouse buttons:
%   left    =   add polygon point
%   right   =   remove last polygon point
%   middle  =   close polygon
%
%  USAGE
%
%    in = UIInPolygon(X,Y)
%
%    X,Y            polygon coordinates
%    in             logical vector indicating whether points at (X,Y) are inside
%
%  SEE
%
%    See also UISelect
%
% Copyright (C) 2009-2011 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.
[Xp,Yp,p] = UISelect;
in = inpolygon(X,Y,Xp,Yp);
hold on;
pts = plot(X(in),Y(in),'+r');
hold off;
pause(0.4);
delete(p);
delete(pts);
cd C:\Users\Cornell\Documents\GitHub\neurocode\utilities\
which UISelect
plot(100,4);
UIInPolygon(rand(10,1));
edit UIInPolygon.m
UIInPolygon(rand(10,1));
which inpolygon
UIInPolygon(rand(10,1),rand(10,1));
raly
X = BatchReturn('OJR_longTraining.batch');
cd N:\OJRproject\OJR42\day12
lfpstruct = getLFP(114+1);
lfp = [lfpstruct.timestamps]; lfp(:,2) = lfpstruct.data;
cd C:
lfpstructH = getLFP(59+1);
lfpH = [lfpstructH.timestamps]; lfpH(:,2) = lfpstructH.data;
speed = (0:0.5:lfp(end,1))'; speed(:,2) = 0;
[freezing,quietWake,SWS,REM] = behavioralStates(lfp,lfpH,speed,1);
addpath 'M:\home\raly\FMAToolbox\Analyses'
cd  M:\home\raly\
startup
which FilterLFP
[freezing,quietWake,SWS,REM] = behavioralStates(lfp,lfpH,speed,1);
PlotIntervals(SWS,'color','b');
PlotIntervals(REM,'color','g');
PlotIntervals(quietWake,'color',[1 1 1]*0.4);
PlotIntervals(freezing,'color','r');
[clean,bad,badIntervals] = CleanLFP(lfp);
PlotIntervals(badIntervals,'color','r');
dur(REM)
dur(SWS)
dur(REM)/60
dur(SWS)/60
dur(quietWake)/60
dur(freezing)/60
[freezing,quietWake,SWS,REM] = behavioralStates(lfp,lfpH,speed,1);
[freezing,quietWake,SWS,REM] = behavioralStates(lfp,lfpH,speed,1);
[freezing,quietWake,SWS,REM] = behavioralStates(lfp,lfpH,speed,1);
PlotIntervals(SWS,'color','b');
PlotIntervals(REM,'color','g');
PlotIntervals(quietWake,'color',[1 1 1]*0.4);
PlotIntervals(freezing,'color','r');
dur(freezing)
dur(freezing)/60
load('day12.EMGFromLFP.LFP.mat')
speed = [EMGFromLFP.timestamps EMGFromLFP.data];
figure; PlotXY(speed);
PlotXY(Smooth(speed,[0 10]));
PlotXY(Smooth(speed,[10 0]));
hold all
PlotXY(Smooth(speed,[15 0]));
smoothed = Smooth(speed,[15 0]);
[freezing,quietWake,SWS,REM] = behavioralStates(lfp,lfpH,speed,0.6);
dur(freezing)/60
PlotIntervals(SWS,'color','b');
PlotIntervals(REM,'color','g');
PlotIntervals(quietWake,'color',[1 1 1]*0.4);
PlotIntervals(freezing,'color','r');
dur(SWS)/60
dur(quietWake)/60
dur(REM)/60
PlotXY(Smooth(speed,[10 0]),'k');
[freezing,quietWake,SWS,REM] = behavioralStates(lfp,lfpH,speed,0.6);
figure; PlotXY(smoothedPower);
x = xlim;
xlim(x);
hold all
PlotXY(spindlePower);
PlotXY(smoothedPower,'k');
1.5/6
6/1.5
min(diff(freezing,[],2))
clf
hist(smoothedPower(:,2),100)
hist(smoothedPower(:,end),100)
mmax(smoothedPower(k==1))
mmax(smoothedPower(k==2))
PlotHVLines(336.21,'v');
PlotHVLines(335,'v');
PlotHVLines(330,'v');
PlotHVLines(325,'v');
clf
PlotXY(smoothedPower,'k');
plot(lfp(:,1),smoothedPower,'k');
plot(tLFP(:,1),smoothedPower,'k');
plot(hpcLFP(:,1),smoothedPower,'k');
plot(tSpindles,smoothedPower,'k');
PlotHVLines([336.21;325]);
PlotHVLines([336.21;325],'h');
PlotXY(Smooth(speed,[10 0]),'k');
PlotXY(Smooth(speed,[20 0]),'k');
PlotIntervals(SWS,'color','b');
PlotIntervals(REM,'color','g');
PlotIntervals(quietWake,'color',[1 1 1]*0.4);
PlotIntervals(freezing,'color','r');
legend('speed','SWS','REM','quietWake','freezing');
clf
handle(1) = PlotXY(Smooth(speed,[20 0]),'k');
handle(2) = PlotIntervals(SWS,'color','b');
handle(3) = PlotIntervals(REM,'color','g');
handle(4) = PlotIntervals(quietWake,'color',[1 1 1]*0.4);
handle(5) = PlotIntervals(freezing,'color','r');
legend('speed','SWS','REM','quietWake','freezing');
legend(handle,{'speed','SWS','REM','quietWake','freezing'});
handle(1) = PlotXY(Smooth(speed,[20 0]),'k','linewidth',2);
handle(1) = PlotXY(Smooth(speed,[20 0]),'k','linewidth',1.5);
clf
handle(1) = PlotXY(Smooth(speed,[20 0]),'k','linewidth',1.5);
handle(2) = PlotIntervals(SWS,'color','b');
handle(3) = PlotIntervals(REM,'color','g');
handle(4) = PlotIntervals(quietWake,'color',[1 1 1]*0.4);
handle(5) = PlotIntervals(freezing,'color','r');
names = {'speed',[num2str(dur(SWS)/3600) 'h SWS'],[num2str(dur(REM)/3600) 'h REM'],...
[num2str(dur(quietWake)/3600) 'h quietWake'],[num2str(dur(freezing)/3600) 'h freezing']};
legend(handle,names);
clf
handle(1) = PlotXY(Smooth(speed,[20 0]),'k','linewidth',1.5);
handle(2) = PlotIntervals(SWS,'color','b');
handle(3) = PlotIntervals(REM,'color','g');
handle(4) = PlotIntervals(quietWake,'color',[1 1 1]*0.4);
handle(5) = PlotIntervals(freezing,'color','r');
names = {'speed',[num2str(dur(SWS)/3600) 'h SWS'],[num2str(dur(REM)/60) 'min REM'],...
[num2str(dur(quietWake)/60) 'min quietWake'],[num2str(dur(freezing)/60) 'min freezing']};
legend(handle,names);
clf
handle(1) = PlotXY(Smooth(speed,[20 0]),'k','linewidth',1.5);
handle(2) = PlotIntervals(SWS,'color','b');
handle(3) = PlotIntervals(REM,'color','g');
handle(4) = PlotIntervals(quietWake,'color',[1 1 1]*0.4);
handle(5) = PlotIntervals(freezing,'color','r');
names = {'speed',[num2str(dur(SWS)/3600) 'h SWS'],[num2str(round(10*dur(REM)/60)/10) 'min REM'],...
[num2str(round(dur(quietWake)/60)) 'min quietWake'],[num2str(round(dur(freezing)/60)) 'min freezing']};
legend(handle,names);
clf
handle(1) = PlotXY(Smooth(speed,[20 0]),'k','linewidth',1.5);
handle(2) = PlotIntervals(SWS,'color','b');
handle(3) = PlotIntervals(REM,'color','g');
handle(4) = PlotIntervals(quietWake,'color',[1 1 1]*0.4);
handle(5) = PlotIntervals(freezing,'color','r');
names = {'speed',[num2str(round(10*dur(SWS)/3600)/10) 'h SWS'],[num2str(round(dur(REM)/60)) 'min REM'],...
[num2str(round(dur(quietWake)/60)) 'min quietWake'],[num2str(round(dur(freezing)/60)) 'min freezing']};
legend(handle,names);
legend(handle,names,'fontsize',15);
set(gca,'fontsize',15
set(gca,'fontsize',15)
xlabel('session time (s)');
set(gca,'ytick',[]);
raly
SaveFig('behavioralStates_OJR42_day12');
min(diff(freezing,[],2))
figure; hist(diff(freezing,[],2));
hist(diff(freezing,[],2),100);
%-- 2/22/2022 11:43 AM --%
p = inputParser;
addParameter(p,'basepath',pwd,@isfolder); % by default, current folder
addParameter(p,'fillMissingDatFiles',false,@islogical);
addParameter(p,'fillTypes',[],@iscellstr);
addParameter(p,'analogInputs',false,@islogical);
addParameter(p,'analogChannels',[],@isnumeric);
addParameter(p,'digitalInputs',false,@islogical);
addParameter(p,'digitalChannels',[],@isnumeric);
addParameter(p,'getAcceleration',false,@islogical);
addParameter(p,'cleanArtifacts',false,@islogical);
addParameter(p,'stateScore',false,@islogical);
addParameter(p,'spikeSort',true,@islogical);
addParameter(p,'cleanRez',false,@islogical);
addParameter(p,'getPos',true,@islogical);
addParameter(p,'removeNoise',false,@islogical); % raly: noise removal is bad, it removes periods 20ms after (because of the filter shifting) a peak in high gamma. See ayadata1\home\raly\Documents\notes\script_NoiseRemoval_bad.m for details.
addParameter(p,'runSummary',false,@islogical);
addParameter(p,'SSD_path','D:\KiloSort',@ischar)    % Path to SSD disk. Make it empty to disable SSD
% addParameter(p,'pullData',[],@isdir); To do...
parse(p,varargin{:});
basepath = p.Results.basepath;
fillMissingDatFiles = p.Results.fillMissingDatFiles;
fillTypes = p.Results.fillTypes;
analogInputs = p.Results.analogInputs;
analogChannels = p.Results.analogChannels;
digitalInputs = p.Results.digitalInputs;
digitalChannels = p.Results.digitalChannels;
getAcceleration = p.Results.getAcceleration;
cleanArtifacts = p.Results.cleanArtifacts;
stateScore = p.Results.stateScore;
spikeSort = p.Results.spikeSort;
cleanRez = p.Results.cleanRez;
getPos = p.Results.getPos;
removeNoise = p.Results.removeNoise;
runSummary = p.Results.runSummary;
SSD_path = p.Results.SSD_path;
if ~exist(basepath,'dir')
error('path provided does not exist')
end
cd(basepath)
p = inputParser;
addParameter(p,'basepath',pwd,@isfolder); % by default, current folder
addParameter(p,'fillMissingDatFiles',false,@islogical);
addParameter(p,'fillTypes',[],@iscellstr);
addParameter(p,'analogInputs',false,@islogical);
addParameter(p,'analogChannels',[],@isnumeric);
addParameter(p,'digitalInputs',false,@islogical);
addParameter(p,'digitalChannels',[],@isnumeric);
addParameter(p,'getAcceleration',false,@islogical);
addParameter(p,'cleanArtifacts',false,@islogical);
addParameter(p,'stateScore',false,@islogical);
addParameter(p,'spikeSort',true,@islogical);
addParameter(p,'cleanRez',false,@islogical);
addParameter(p,'getPos',true,@islogical);
addParameter(p,'removeNoise',false,@islogical); % raly: noise removal is bad, it removes periods 20ms after (because of the filter shifting) a peak in high gamma. See ayadata1\home\raly\Documents\notes\script_NoiseRemoval_bad.m for details.
addParameter(p,'runSummary',false,@islogical);
addParameter(p,'SSD_path','D:\KiloSort',@ischar)    % P
parse(p);
basepath = p.Results.basepath;
fillMissingDatFiles = p.Results.fillMissingDatFiles;
fillTypes = p.Results.fillTypes;
analogInputs = p.Results.analogInputs;
analogChannels = p.Results.analogChannels;
digitalInputs = p.Results.digitalInputs;
digitalChannels = p.Results.digitalChannels;
getAcceleration = p.Results.getAcceleration;
cleanArtifacts = p.Results.cleanArtifacts;
stateScore = p.Results.stateScore;
spikeSort = p.Results.spikeSort;
cleanRez = p.Results.cleanRez;
getPos = p.Results.getPos;
removeNoise = p.Results.removeNoise;
runSummary = p.Results.runSummary;
SSD_path = p.Results.SSD_path;
if ~exist(basepath,'dir')
error('path provided does not exist')
end
cd(basepath)
% Get session names
if strcmp(basepath(end),filesep)
basepath = basepath(1:end-1);
end
[~,basename] = fileparts(basepath);
% Get xml file in order
xmlFile = checkFile('fileType','.xml','searchSubdirs',true);
xmlFile = xmlFile(1);
if ~(strcmp(xmlFile.folder,basepath)&&strcmp(xmlFile.name(1:end-4),basename))
copyfile([xmlFile.folder,filesep,xmlFile.name],[basepath,filesep,basename,'.xml'])
end
end
xmlFile = checkFile('fileType','.xml','searchSubdirs',true);
filename
isempty(filename) && ~isempty(fileType)
file = dir([basepath,filesep,'*',fileType]);
file
isempty(file) && searchSubdirs
subFile = dir([basepath,filesep,'*',filesep,'*',fileType])
subFile.name
subFile(1).name = [];
bad = subFile.name
bad = subFile(1:end).name
bad = strcmp(subFile(1:end).name,'settings.xml');
bad = cellfun(@(x) strcmp(x,'settings.xml'),{subFile(1:end).name});
bad = cellfun(@(x) strcmp(x,'settings.xml'),{subFile(1:end).name})';
bad
subFile = dir([basepath,filesep,'*',filesep,'*',fileType]);
settings(bad,:) = [];
subFile(bad) = [];
subFile = dir([basepath,filesep,'*',filesep,'*',fileType]);
ignore = cellfun(@(x) strcmp(x,'settings.xml'),{subFile(1:end).name})'; % ignore "settings.xml" files
subFile(ignore) = [];
file = [file; subFile];
~isempty(file)
file
searchSuperdirs
mydir  = pwd;
idcs   = strfind(mydir,'\');
idcs
newdir = mydir(1:idcs(end)-1);
newdir
newdir = fileparts(basepath)
newdir = fileparts(basepath);
superFile = dir([newdir,filesep,'*',fileType]);
superFile
xmlFile = checkFile('fileType','.xml','searchSubdirs',true);
xmlFile = xmlFile(1);
xmlFile
~(strcmp(xmlFile.folder,basepath)&&strcmp(xmlFile.name(1:end-4),basename))
% Check info.rhd
% (assumes this will be the same across subsessions)
rhdFile = checkFile('fileType','.rhd','searchSubdirs',true);
rhdFile = rhdFile(1);
if ~(strcmp(rhdFile.folder,basepath)&&strcmp(rhdFile.name(1:end-4),basename))
copyfile([rhdFile.folder,filesep,rhdFile.name],[basepath,filesep,basename,'.rhd'])
end
session = sessionTemplate(basepath,'showGUI',false);
sessionInfo = LoadXml(fullfile(session.general.basePath,[session.general.name, '.xml']));
sessionInfo
sessionInfo = LoadXml(fullfile(session.general.basePath,[session.general.name, '.xml']));
isfield(sessionInfo,'SpkGrps')
session.extracellular
sessionInfo.AnatGrps
! ipython
evalc(['! ipython']); evalc(['! print(''hello'')'])
evalc(['! ipython']);
eval(['! ipython']);
eval(['! ipython']); eval(['! print(''hello'')'])
eval(['! ipython; print(''hello'')'])
eval(['! ipython | print(''hello'')'])
foo = sprintf('\n\n\n5')
a = eval('foo')
foo = sprintf('ipython \n print(''hello'')')
foo = sprintf('! ipython \n print(''hello'')')
eval(['! ipython | print(''hello'')'])
eval(foo)
foo = sprintf('! ipython \n print(''hello'')')
eval(foo)
foo = sprintf('! ipython \n print(''hello'') \n print(''hello2'')')
eval(foo)
eval('! ipython; print(''hello2'')
eval('! ipython; print(''hello2'')')
eval('! ipython ; print(''hello2'')')
eval('! ipython , print(''hello2'')')
eval('! ipython \n print(''hello2'')')
eval('! ipython print(''hello2'')')
eval('! ipython')
foo = sprintf('! ipython \n print(''hello'') \n print(''hello2'')')
eval(foo)
eval('! ipython')
eval('echo "hello2"')
pyrun('echo "hello2"')
pyrun('print(''hello2'')')
pyrun('version')
pyrun('print(''hello''), print(''hello2''')
pyrun('print(''hello'') , print(''hello2''')
pyrun('print(''hello'') ; print(''hello2''')
pyrun('print(''hello''); print(''hello2''')
pyrun('print(''hello'') \n print(''hello2''')
foo = 'print(''hello'') \n print(''hello2''';
pyrun(foo)
foo
foo = sprintf('print('hello') \n print('hello2'')')
foo = sprintf('print('hello') \n print('hello2'')
foo = sprintf('print('hello') \n print('hello2')')
foo = sprintf('print(''hello'') \n print(''hello2'')')
pyrun(foo)
foo = sprintf('print(''hello'')\nprint(''hello2'')')
pyrun(foo)
script = 'print(''hello'')\nprint(''hello2'')\nsubprocess.run(''conda activate C:\Users\Cornell\anaconda3\envs\DEEPLABCUT && "print(''hello2'')" && conda deactivate'', shell=True)';
foo
which pyrun
pyversion
cd C:
pathToCode = 'C:\Users\Cornell\Documents\GitHub\neurocode\future_script.py'
count(py.sys.path,pathToCode)==0
if count(py.sys.path,pathToCode)==0
insert(py.sys.path,int32(0),pathToCode)
end
pathToCode
% future_script = mySpeechRecognizer.py; fun_name = audioToText
pathToCode = fileparts(which('future_script.py'))
%
% start
pyversion
pathToCode = fileparts(which('future_script.py'))
if count(py.sys.path,pathToCode)==0
insert(py.sys.path,int32(0),pathToCode)
end
% start
pyversion
pathToCode = fileparts(which('future_script.py'))
if count(py.sys.path,pathToCode)==0
insert(py.sys.path,int32(0),pathToCode) % adds the code to the python path
end
pyOut = py.future_script.fun_name(basepath);
pathToCode = fileparts(which('test_py.py'))
edit test_py.py
sessionInfo = LoadXml(fullfile(session.general.basePath,[session.general.name, '.xml']));
session = sessionTemplate(basepath,'showGUI',false);
disp('Loading Neurosuite xml file metadata')
sessionInfo = LoadXml(fullfile(session.general.basePath,[session.general.name, '.xml']));
rxml.child
rxml.child.tag
rxml.child(3).tag
rxml.child(3)
rxml.child(3).attribs
xml
v1
sessionInfo0 = LoadXml('N:\V1test\V1Jean\day7\day7.xml');
dbcont
isfield(xml,'AnatGrps')
isfield(sessionInfo0,'AnatGrps')
session = sessionTemplate(basepath,'showGUI',false);
disp('Loading Neurosuite xml file metadata')
sessionInfo = LoadXml(fullfile(session.general.basePath,[session.general.name, '.xml']));
if isfield(sessionInfo,'SpkGrps')
session.extracellular.nSpikeGroups = length(sessionInfo.SpkGrps); % Number of spike groups
session.extracellular.spikeGroups.channels = {sessionInfo.SpkGrps.Channels}; % Spike groups
else
warning('No spike groups exist in the xml. Anatomical groups used instead')
session.extracellular.nSpikeGroups = size(sessionInfo.AnatGrps,2); % Number of spike groups
session.extracellular.spikeGroups.channels = {sessionInfo.AnatGrps.Channels}; % Spike groups
end
session = sessionTemplate(basepath,'showGUI',false);
session.extracellular
133/4
33.2/(1/20000)
9*60+47
session = sessionTemplate(basepath,'showGUI',false);
session
session.extracellular
session = sessionTemplate(basepath,'showGUI',false);
xml
rxml
rxml.child
rxml.child.tag
for i=1:length(rxml.child)
switch lower(rxml.child(i).tag)
case 'generalinfo'
xml.Date = rxml.child(i).child(1).value; % date of xml file creation?
case 'acquisitionsystem'
xml.nBits = str2num(rxml.child(i).child(1).value); % number of bits of the file
xml.nChannels = str2num(rxml.child(i).child(2).value);
xml.SampleRate = str2num(rxml.child(i).child(3).value);
xml.SampleTime = 1e6/xml.SampleRate; %to make backwards compatible
xml.VoltageRange = str2num(rxml.child(i).child(4).value);
xml.Amplification = str2num(rxml.child(i).child(5).value);
xml.Offset =  str2num(rxml.child(i).child(6).value);
case 'fieldpotentials'
xml.lfpSampleRate = str2num(rxml.child(i).child.value);
case 'anatomicaldescription'
tmp = rxml.child(i).child.child;
for grpI =1:length(tmp)
for chI=1:length(tmp(grpI).child)
xml.AnatGrps(grpI).Channels(chI) = str2num(tmp(grpI).child(chI).value);
xml.AnatGrps(grpI).Skip(chI) = str2num(tmp(grpI).child(chI).attribs.value);
end
end
case 'spikedetection'
if ~isempty(rxml.child(i).child)
tmp =rxml.child(i).child.child;
for grpI =1:length(tmp)
for chI=1:length(tmp(grpI).child(1).child)
xml.SpkGrps(grpI).Channels(chI) = str2num(tmp(grpI).child(1).child(chI).value);
end
if length(tmp(grpI).child)>1
xml.SpkGrps(grpI).nSamples = str2num(tmp(grpI).child(2).value);
xml.SpkGrps(grpI).PeakSample = str2num(tmp(grpI).child(3).value);
xml.SpkGrps(grpI).nFeatures = str2num(tmp(grpI).child(4).value);
end
%backwards compatibility
xml.nElecGps = length(tmp);
xml.ElecGp{grpI} = xml.SpkGrps(grpI).Channels;
end
else
xml.nElecGps = 0;
end
case 'programs'
tmp = rxml.child(i).child;
for i=1:length(tmp)
if strcmp(tmp(i).child(1).value,'process_mhipass')
for j=1:length(tmp(i).child(2).child )
if strcmp(tmp(i).child(2).child(j).child(1).value,'frequency')
xml.HiPassFreq = str2num(tmp(i).child(2).child(j).child(2).value);
break
end
end
end
end
end
end
xml
xml.AnatGrps.Channels
xml.SampleRate
xml.lfpSampleRate
xml.nChannels
session = sessionTemplate(basepath,'showGUI',false);
digitalInp = getDigitalIn('all','fs',session.extracellular.sr);
session = sessionTemplate(basepath,'showGUI',true);
digitalInp = getDigitalIn('all','fs',session.extracellular.sr);
filename
filename=dir('digitalIn.dat');
filename
session
session.epochs
sessionP = session;
load('day7.session.mat')
sessionV = session;
session= sessionP;
open preprocessV1
digitalInp = getDigitalIn('all','fs',sessionV.extracellular.sr);
isempty(filename)
filename=dir('digitalIn.dat')
filename = filename.name;
~isempty(dir('*.xml'))
sess = getSession;
sess
~isempty(dir('*DigitalIn.events.mat'))
isempty(dir('*DigitalIn.events.mat'))
isempty(filename)
filename=dir('digitalIn.dat');
filename = filename.name;
filename
digitalInp = getDigitalIn('all','fs',session.extracellular.sr);
isempty(dir('*DigitalIn.events.mat'))
isempty(filename)
filename=dir('digitalIn.dat');
filename
basepath
concatenateDats(basepath,0,1);
datpaths
otherdattypes
newpaths
datpaths
isempty(datpaths)
%% Get the XML
try
%Look for xml/sessionInfo in topfolder
%sessionInfo = bz_getSessionInfo(basepath,'noPrompts',true);
load([basename '.session.mat']); % Peter's sessionInfo
catch
%If none exists, look for xml in any of the subpaths
disp('No .xml or .sessionInfo in top folder, trying subfolders')
for ff = 1:length(recordingnames)
try
sessionInfo = LoadParameters(fullfile(basepath,recordingnames{ff}));
xmlfilename = fullfile(sessionInfo.session.path,[sessionInfo.session.name,'.xml']);
[SUCCESS,MESSAGE,MESSAGEID] = copyfile(...
xmlfilename,...
fullfile(basepath,[basename,'.xml']),'f');
display(['Copied xml from ',recordingnames{ff}])
break
catch
end
end
end
datpaths
otherdattypes
bad_otherdattypes
%% Handling inputs
% basic session name and and path
if ~exist('basepath','var')
basepath = cd;
end
basename = basenameFromBasepath(basepath);
if ~exist('deleteoriginaldatsbool','var')
deleteoriginaldatsbool = 0;
end
if ~exist('sortFiles','var')
sortFiles = 1;
end
otherdattypes = {'analogin';'digitalin';'auxiliary';'time';'supply'};
bad_otherdattypes = [];
for odidx = 1:length(otherdattypes)
%eval(['new' otherdattypes{odidx} 'path = fullfile(basepath,''' otherdattypes{odidx} '.dat'');'])
newpaths.(otherdattypes{odidx}) = fullfile(basepath,[otherdattypes{odidx}, '.dat']);
end
d = dir(basepath);
datpaths = {};
datsizes.amplifier = [];
recordingnames = {};
rcount = 0; %Count of good subfolders
d
i=6
a=6
d(a).isdir
exist(fullfile(basepath,d(a).name,[d(a).name,'.dat']),'file')
fullfile(basepath,d(a).name,[d(a).name,'.dat'])
fullfile(basepath,d(a).name,'amplifier_analogin_auxiliary_int16.dat')
concatenateDats(basepath,0,1);
a=6
d(a).isdir
if exist(fullfile(basepath,d(a).name,[d(a).name,'.dat']),'file')
ampfile = fullfile(basepath,d(a).name,[d(a).name,'.dat']);
elseif exist(fullfile(basepath,d(a).name,'amplifier_analogin_auxiliary_int16.dat'),'file')%Luke Sjulson's Modified code to record all 16bit signals in one file
ampfile = fullfile(basepath,d(a).name,'amplifier_analogin_auxiliary_int16.dat');
else
ampfile = fullfile(basepath,d(a).name,'amplifier.dat');
end
ampfile
exist(ampfile,'file')
otherdattypes = {'analogin';'digitalin';'auxiliary';'time';'supply'};
bad_otherdattypes = [];
for odidx = 1:length(otherdattypes)
%eval(['new' otherdattypes{odidx} 'path = fullfile(basepath,''' otherdattypes{odidx} '.dat'');'])
newpaths.(otherdattypes{odidx}) = fullfile(basepath,[otherdattypes{odidx}, '.dat']);
end
d = dir(basepath);
datpaths = {};
datsizes.amplifier = [];
recordingnames = {};
rcount = 0; %Count of good subfolders
for a = 1:length(d)
%look in each subfolder
if d(a).isdir
%Check for amplifier.dat or subfolderbaseName.dat
if exist(fullfile(basepath,d(a).name,[d(a).name,'.dat']),'file')
ampfile = fullfile(basepath,d(a).name,[d(a).name,'.dat']);
elseif exist(fullfile(basepath,d(a).name,'amplifier_analogin_auxiliary_int16.dat'),'file')%Luke Sjulson's Modified code to record all 16bit signals in one file
ampfile = fullfile(basepath,d(a).name,'amplifier_analogin_auxiliary_int16.dat');
else
ampfile = fullfile(basepath,d(a).name,'amplifier.dat');
end
if exist(ampfile,'file')
rcount = rcount+1;
datpaths.amplifier{rcount} = ampfile;
t = dir(ampfile);
datsizes.amplifier(rcount) = t.bytes;
recordingnames{rcount} = d(a).name;
end
for odidx = 1:length(otherdattypes)%loop through other .dat types found here
% eval([otherdattypes{odidx} 'datpaths{rcount} = fullfile(basepath,recordingnames{rcount},''' otherdattypes{odidx} '.dat'');'])
datpaths.(otherdattypes{odidx}){rcount} = fullfile(basepath,recordingnames{rcount},[otherdattypes{odidx} '.dat']);
%eval(['d2 = dir(' otherdattypes{odidx} 'datpaths.amplifier{rcount});'])
d2 = dir(datpaths.(otherdattypes{odidx}){rcount});
if isempty(d2)
bad_otherdattypes(odidx) = 1;
else
%eval([otherdattypes{odidx} 'datsizes(rcount) = d2(1).bytes;'])
datsizes.(otherdattypes{odidx})(rcount) = d2(1).bytes;
end
end
end
end
a
d(a).isdir
open d
d(a).isdir
a
d(a)
isdir(d(a))
isfolder(fileparts(d(a).folder,d(a).name))
isfolder(fullfile(d(a).folder,d(a).name))
fullfile(d(a).folder,d(a).name)
any(~ismember(d(a).name,'.'))
~ismember(d(a).name,'.')
d = dir(basepath);
datpaths = {};
datsizes.amplifier = [];
recordingnames = {};
rcount = 0; %Count of good subfolders
for a = 1:length(d)
%look in each subfolder
if any(~ismember(d(a).name,'.')) && d(a).isdir
%Check for amplifier.dat or subfolderbaseName.dat
if exist(fullfile(basepath,d(a).name,[d(a).name,'.dat']),'file')
ampfile = fullfile(basepath,d(a).name,[d(a).name,'.dat']);
elseif exist(fullfile(basepath,d(a).name,'amplifier_analogin_auxiliary_int16.dat'),'file')%Luke Sjulson's Modified code to record all 16bit signals in one file
ampfile = fullfile(basepath,d(a).name,'amplifier_analogin_auxiliary_int16.dat');
else
ampfile = fullfile(basepath,d(a).name,'amplifier.dat');
end
if exist(ampfile,'file')
rcount = rcount+1;
datpaths.amplifier{rcount} = ampfile;
t = dir(ampfile);
datsizes.amplifier(rcount) = t.bytes;
recordingnames{rcount} = d(a).name;
end
for odidx = 1:length(otherdattypes)%loop through other .dat types found here
% eval([otherdattypes{odidx} 'datpaths{rcount} = fullfile(basepath,recordingnames{rcount},''' otherdattypes{odidx} '.dat'');'])
datpaths.(otherdattypes{odidx}){rcount} = fullfile(basepath,recordingnames{rcount},[otherdattypes{odidx} '.dat']);
%eval(['d2 = dir(' otherdattypes{odidx} 'datpaths.amplifier{rcount});'])
d2 = dir(datpaths.(otherdattypes{odidx}){rcount});
if isempty(d2)
bad_otherdattypes(odidx) = 1;
else
%eval([otherdattypes{odidx} 'datsizes(rcount) = d2(1).bytes;'])
datsizes.(otherdattypes{odidx})(rcount) = d2(1).bytes;
end
end
end
end
a
any(~ismember(d(a).name,'.')) && d(a).isdir
if exist(fullfile(basepath,d(a).name,[d(a).name,'.dat']),'file')
ampfile = fullfile(basepath,d(a).name,[d(a).name,'.dat']);
elseif exist(fullfile(basepath,d(a).name,'amplifier_analogin_auxiliary_int16.dat'),'file')%Luke Sjulson's Modified code to record all 16bit signals in one file
ampfile = fullfile(basepath,d(a).name,'amplifier_analogin_auxiliary_int16.dat');
else
ampfile = fullfile(basepath,d(a).name,'amplifier.dat');
end
if exist(ampfile,'file')
rcount = rcount+1;
datpaths.amplifier{rcount} = ampfile;
t = dir(ampfile);
datsizes.amplifier(rcount) = t.bytes;
recordingnames{rcount} = d(a).name;
end
odidx
datpaths.(otherdattypes{odidx}){rcount} = fullfile(basepath,recordingnames{rcount},[otherdattypes{odidx} '.dat']);
datpaths
fullfile(basepath,recordingnames{rcount},[otherdattypes{odidx} '.dat']);
recordingnames
d(a).name
fullfile(basepath,d(a).name,[otherdattypes{odidx} '.dat'])
datpaths.(otherdattypes{odidx}){rcount} = fullfile(basepath,d(a).name,[otherdattypes{odidx} '.dat']);
odidx
fullfile(basepath,d(a).name,[otherdattypes{odidx} '.dat'])
datpaths.(otherdattypes{odidx}){odidx} = fullfile(basepath,d(a).name,[otherdattypes{odidx} '.dat']);
d2 = dir(datpaths.(otherdattypes{odidx}){rcount})
d2 = dir(datpaths.(otherdattypes{odidx}){odidx})
digitalinfile = fullfile(basepath,d(a).name,'digitalin.dat');
exist(digitalinfile,'file')
otherdattypes = {'analogin';'digitalin';'auxiliary';'time';'supply'};
bad_otherdattypes = [];
for odidx = 1:length(otherdattypes)
%eval(['new' otherdattypes{odidx} 'path = fullfile(basepath,''' otherdattypes{odidx} '.dat'');'])
newpaths.(otherdattypes{odidx}) = fullfile(basepath,[otherdattypes{odidx}, '.dat']);
end
d = dir(basepath);
datpaths = {};
datsizes.amplifier = [];
recordingnames = {};
rcount = 0; %Count of good subfolders
for a = 1:length(d)
%look in each subfolder
if any(~ismember(d(a).name,'.')) && d(a).isdir
%Check for amplifier.dat or subfolderbaseName.dat
if exist(fullfile(basepath,d(a).name,[d(a).name,'.dat']),'file')
ampfile = fullfile(basepath,d(a).name,[d(a).name,'.dat']);
elseif exist(fullfile(basepath,d(a).name,'amplifier_analogin_auxiliary_int16.dat'),'file')%Luke Sjulson's Modified code to record all 16bit signals in one file
ampfile = fullfile(basepath,d(a).name,'amplifier_analogin_auxiliary_int16.dat');
else
ampfile = fullfile(basepath,d(a).name,'amplifier.dat');
digitalinfile = fullfile(basepath,d(a).name,'digitalin.dat');
end
if exist(ampfile,'file')
rcount = rcount+1;
datpaths.amplifier{rcount} = ampfile;
t = dir(ampfile);
datsizes.amplifier(rcount) = t.bytes;
recordingnames{rcount} = d(a).name;
for odidx = 1:length(otherdattypes)%loop through other .dat types found here
% eval([otherdattypes{odidx} 'datpaths{rcount} = fullfile(basepath,recordingnames{rcount},''' otherdattypes{odidx} '.dat'');'])
datpaths.(otherdattypes{odidx}){rcount} = fullfile(basepath,recordingnames{rcount},[otherdattypes{odidx} '.dat']);
%eval(['d2 = dir(' otherdattypes{odidx} 'datpaths.amplifier{rcount});'])
d2 = dir(datpaths.(otherdattypes{odidx}){rcount});
if isempty(d2)
bad_otherdattypes(odidx) = 1;
else
%eval([otherdattypes{odidx} 'datsizes(rcount) = d2(1).bytes;'])
datsizes.(otherdattypes{odidx})(rcount) = d2(1).bytes;
end
end
elseif exist(digitalinfile,'file')
rcount = rcount+1;
for odidx = 1:length(otherdattypes)%loop through other .dat types found here
datpaths.(otherdattypes{odidx}){odidx} = fullfile(basepath,d(a).name,[otherdattypes{odidx} '.dat']);
d2 = dir(datpaths.(otherdattypes{odidx}){rcount});
if isempty(d2)
bad_otherdattypes(odidx) = 1;
else
%eval([otherdattypes{odidx} 'datsizes(rcount) = d2(1).bytes;'])
datsizes.(otherdattypes{odidx})(rcount) = d2(1).bytes;
end
end
end
end
end
a
exist(ampfile,'file')
exist(digitalinfile,'file')
rcount
odidx
datpaths.(otherdattypes{odidx}){odidx} = fullfile(basepath,d(a).name,[otherdattypes{odidx} '.dat']);
d2 = dir(datpaths.(otherdattypes{odidx}){rcount});
datpaths.(otherdattypes{odidx}){rcount}
otherdattypes{odidx}
fullfile(basepath,d(a).name,[otherdattypes{odidx} '.dat'])
otherdattypes = {'analogin';'digitalin';'auxiliary';'time';'supply'};
bad_otherdattypes = [];
for odidx = 1:length(otherdattypes)
%eval(['new' otherdattypes{odidx} 'path = fullfile(basepath,''' otherdattypes{odidx} '.dat'');'])
newpaths.(otherdattypes{odidx}) = fullfile(basepath,[otherdattypes{odidx}, '.dat']);
end
d = dir(basepath);
datpaths = {};
datsizes.amplifier = [];
recordingnames = {};
rcount = 0; %Count of good subfolders
for a = 1:length(d)
%look in each subfolder
if any(~ismember(d(a).name,'.')) && d(a).isdir
%Check for amplifier.dat or subfolderbaseName.dat
if exist(fullfile(basepath,d(a).name,[d(a).name,'.dat']),'file')
ampfile = fullfile(basepath,d(a).name,[d(a).name,'.dat']);
elseif exist(fullfile(basepath,d(a).name,'amplifier_analogin_auxiliary_int16.dat'),'file')%Luke Sjulson's Modified code to record all 16bit signals in one file
ampfile = fullfile(basepath,d(a).name,'amplifier_analogin_auxiliary_int16.dat');
else
ampfile = fullfile(basepath,d(a).name,'amplifier.dat');
digitalinfile = fullfile(basepath,d(a).name,'digitalin.dat');
end
if exist(ampfile,'file')
rcount = rcount+1;
datpaths.amplifier{rcount} = ampfile;
t = dir(ampfile);
datsizes.amplifier(rcount) = t.bytes;
recordingnames{rcount} = d(a).name;
for odidx = 1:length(otherdattypes)%loop through other .dat types found here
% eval([otherdattypes{odidx} 'datpaths{rcount} = fullfile(basepath,recordingnames{rcount},''' otherdattypes{odidx} '.dat'');'])
datpaths.(otherdattypes{odidx}){rcount} = fullfile(basepath,recordingnames{rcount},[otherdattypes{odidx} '.dat']);
%eval(['d2 = dir(' otherdattypes{odidx} 'datpaths.amplifier{rcount});'])
d2 = dir(datpaths.(otherdattypes{odidx}){rcount});
if isempty(d2)
bad_otherdattypes(odidx) = 1;
else
%eval([otherdattypes{odidx} 'datsizes(rcount) = d2(1).bytes;'])
datsizes.(otherdattypes{odidx})(rcount) = d2(1).bytes;
end
end
elseif exist(digitalinfile,'file')
rcount = rcount+1;
for odidx = 1:length(otherdattypes)%loop through other .dat types found here
datpaths.(otherdattypes{odidx}){rcount} = fullfile(basepath,d(a).name,[otherdattypes{odidx} '.dat']);
d2 = dir(datpaths.(otherdattypes{odidx}){rcount});
if isempty(d2)
bad_otherdattypes(odidx) = 1;
else
%eval([otherdattypes{odidx} 'datsizes(rcount) = d2(1).bytes;'])
datsizes.(otherdattypes{odidx})(rcount) = d2(1).bytes;
end
end
end
end
end
bad_otherdattypes
otherdattypes(find(bad_otherdattypes)) = [];%if there weren't analogin or digitalin in some recording
isempty(datpaths) || isempty(datpaths.amplifier)
otherdattypes = {'analogin';'digitalin';'auxiliary';'time';'supply'};
bad_otherdattypes = [];
for odidx = 1:length(otherdattypes)
%eval(['new' otherdattypes{odidx} 'path = fullfile(basepath,''' otherdattypes{odidx} '.dat'');'])
newpaths.(otherdattypes{odidx}) = fullfile(basepath,[otherdattypes{odidx}, '.dat']);
end
d = dir(basepath);
datpaths = {};
datsizes.amplifier = [];
recordingnames = {};
rcount = 0; %Count of good subfolders
for a = 1:length(d)
%look in each subfolder
if any(~ismember(d(a).name,'.')) && d(a).isdir
%Check for amplifier.dat or subfolderbaseName.dat
if exist(fullfile(basepath,d(a).name,[d(a).name,'.dat']),'file')
ampfile = fullfile(basepath,d(a).name,[d(a).name,'.dat']);
elseif exist(fullfile(basepath,d(a).name,'amplifier_analogin_auxiliary_int16.dat'),'file')%Luke Sjulson's Modified code to record all 16bit signals in one file
ampfile = fullfile(basepath,d(a).name,'amplifier_analogin_auxiliary_int16.dat');
else
ampfile = fullfile(basepath,d(a).name,'amplifier.dat');
digitalinfile = fullfile(basepath,d(a).name,'digitalin.dat');
end
if exist(ampfile,'file')
rcount = rcount+1;
datpaths.amplifier{rcount} = ampfile;
t = dir(ampfile);
datsizes.amplifier(rcount) = t.bytes;
recordingnames{rcount} = d(a).name;
for odidx = 1:length(otherdattypes)%loop through other .dat types found here
% eval([otherdattypes{odidx} 'datpaths{rcount} = fullfile(basepath,recordingnames{rcount},''' otherdattypes{odidx} '.dat'');'])
datpaths.(otherdattypes{odidx}){rcount} = fullfile(basepath,recordingnames{rcount},[otherdattypes{odidx} '.dat']);
%eval(['d2 = dir(' otherdattypes{odidx} 'datpaths.amplifier{rcount});'])
d2 = dir(datpaths.(otherdattypes{odidx}){rcount});
if isempty(d2)
bad_otherdattypes(odidx) = 1;
else
%eval([otherdattypes{odidx} 'datsizes(rcount) = d2(1).bytes;'])
datsizes.(otherdattypes{odidx})(rcount) = d2(1).bytes;
end
end
elseif exist(digitalinfile,'file')
rcount = rcount+1;
datpaths.amplifier{rcount} = [];
for odidx = 1:length(otherdattypes)%loop through other .dat types found here
datpaths.(otherdattypes{odidx}){rcount} = fullfile(basepath,d(a).name,[otherdattypes{odidx} '.dat']);
d2 = dir(datpaths.(otherdattypes{odidx}){rcount});
if isempty(d2)
bad_otherdattypes(odidx) = 1;
else
%eval([otherdattypes{odidx} 'datsizes(rcount) = d2(1).bytes;'])
datsizes.(otherdattypes{odidx})(rcount) = d2(1).bytes;
end
end
end
end
end
otherdattypes(find(bad_otherdattypes)) = [];%if there weren't analogin or digitalin in some recording
isempty(datpaths)
isempty(datpaths.amplifier)
datpaths.amplifier
if sortFiles
try
names2sort = cellfun(@(X) str2num(X(end-5:end)),recordingnames,'UniformOutput',false);
names2sort = cell2mat(names2sort);
%         if isempty(names2sort{1}) && ~isempty(recordingnames{1})
%             error('Last 6 digits were not numeric and therefore do not reflect the recording time.');
%         end
disp('Assuming the last 6 digits reflect recording time.')
%disp('Don''t like it? Write in some new options for sorting.')
catch
% names2sort = 1:length(recordingnames);
disp('Last 6 digits not numeric... sorting alphabetically')
end
[~,I] = sort(names2sort);
recordingnames = recordingnames(I);
datpaths.amplifier = datpaths.amplifier(I);
datsizes.amplifier = datsizes.amplifier(I);
for odidx = 1:length(otherdattypes)
datpaths.(otherdattypes{odidx}) = datpaths.(otherdattypes{odidx})(I);
datsizes.(otherdattypes{odidx}) = datsizes.(otherdattypes{odidx})(I);
end
%     % save txt with order of files to concatenate (moved to events.mat)
%     fid = fopen(fullfile(basepath,'concatORDER.txt'),'w');
%     for idx = 1:length(I)
%         fprintf(fid,[recordingnames{idx} '\n']);
%     end
%     fclose(fid);
end
newdatpath = fullfile(basepath,[basename,'.dat']);
isunix
ispc
length(datpaths.amplifier)>1
datpathsplus
datpaths.amplifier
datpathsplus = datpaths.amplifier;
cs = strjoin(datpathsplus);
catstring = ['! copy /b ', cs, ' ',newdatpath];
disp('Concatenating Amplifier Dats... be patient')
catstring
eval(catstring)%execute concatention
t = dir(newdatpath);
t.bytes
sum(datsizes.amplifier)
t.bytes ~= sum(datsizes.amplifier)
amplifier
ampfile
fclose(fopen(ampfile, 'w'))
5
exist(ampfile,'file')
45
concatenateDats(basepath,0,1);
dbcont
isempty(datpaths)
isempty(datpaths.amplifier)
sortFiles
names2sort = cellfun(@(X) str2num(X(end-5:end)),recordingnames,'UniformOutput',false);
concatenateDats(basepath,0,1);
ff
datpaths.time
for ff = 1:length(datpaths.time)
f = fopen(datpaths.time{ff},'r');
% Determine total number of samples in file
fileStart = ftell(f);
%Read the first time point
firsttimepoint = fread(f,1,'int32');
status = fseek(f,-4,'eof'); %int32 = 4 bytes
lasttimepoint = fread(f,1,'int32');
fileStop = ftell(f);
firstlasttimepoints(ff,:) = [firsttimepoint lasttimepoint];
numsamples(ff) = fileStop./4;
if ff==1
transitiontimes_samp = firstlasttimepoints(ff,:);
else
transitiontimes_samp(ff,:) = firstlasttimepoints(ff,:)+transitiontimes_samp(ff-1,2)+1;
end
end
disp(['Calculating merge times based on wideband samplingRate of ',num2str(session.extracellular.sr),'Hz.'])
transitiontimes_sec = transitiontimes_samp./session.extracellular.sr; %convert to seconds
eventsfilename = fullfile(basepath,[basename,'.MergePoints.events.mat']);
MergePoints.timestamps = transitiontimes_sec;
MergePoints.timestamps_samples = transitiontimes_samp;
MergePoints.firstlasttimpoints_samples = firstlasttimepoints;
MergePoints.foldernames = recordingnames;
MergePoints.filesmerged = datpaths;
MergePoints.filesizes = datsizes;
MergePoints.sizecheck = sizecheck;
MergePoints.detectorinfo.detectorname = 'bz_ConcatenateDats';
MergePoints.detectorinfo.detectiondate = datestr(now,'yyyy-mm-dd');
%Saving SleepStates
save(eventsfilename,'MergePoints');
MergePoints
MergePoints.timestamps
deleteoriginaldatsbool
clear all
basepath = pwdd;
basepath = pwd;
concatenateDats(basepath,0,1);
sortFiles
names2sort = cellfun(@(X) str2num(X(end-5:end)),recordingnames,'UniformOutput',false);
names2sort = cell2mat(names2sort);
disp('Assuming the last 6 digits reflect recording time.')
[~,I] = sort(names2sort)
recordingnames = recordingnames(I)
datpaths.amplifier
datpaths.amplifier = datpaths.amplifier(I);
datsizes.amplifier = datsizes.amplifier(I);
for odidx = 1:length(otherdattypes)
datpaths.(otherdattypes{odidx}) = datpaths.(otherdattypes{odidx})(I);
datsizes.(otherdattypes{odidx}) = datsizes.(otherdattypes{odidx})(I);
end
concatenateDats(basepath,0,1);
load('day2.MergePoints.events.mat')
session = sessionTemplate(basepath,'showGUI',true);
session
open session
subdit
subdir
dir
d = dir;
d = dir(basepath)
open d
names = d.name;
names = d(:).name;
names = {d(:).name};
names = {d(:).name}';
names = {d(:).name}'; names(cellfun(@(x),~any(~ismember(x),'.'),names)) = [];
cellfun(@(x),~any(~ismember(x),'.'),names)
names = {d(:).name}'; names(cellfun(@(x) ~any(~ismember(x),'.'),names)) = [];
cellfun(@(x) ~any(~ismember(x),'.'),names)
cellfun(@(x) ~any(~ismember(x,'.')),names)
names = {d(:).name}'; names(cellfun(@(x) ~any(~ismember(x,'.')),names)) = [];
names = {d(d.isdir).name}'; names(cellfun(@(x) ~any(~ismember(x,'.')),names)) = [];
names = {d(d.isdir).name}';
d.isdir
names = {d(cat(1,d.isdir)).name}';
names = {d(cat(1,d.isdir)).name}'; names(cellfun(@(x) ~any(~ismember(x,'.')),names)) = [];
cellfun(fun_existsdate,subdirs)
subdirs = {d(cat(1,d.isdir)).name}'; subdirs(cellfun(@(x) ~any(~ismember(x,'.')),subdirs)) = [];
cellfun(fun_existsdate,subdirs)
all(cellfun(fun_existsdate,subdirs))
cellfun(fun_isdate,subdirs)
cellfun(fun_isdate,subdirs,'UniformOutput',false)
getdate = @(x) str2double(x(fun_isdate(x)));
getDate(subdirs{1})
getdate(subdirs{1})
getdate = @(x) (x(fun_isdate(x)));
getdate(subdirs{1})
datenum(getdate(subdirs{1}))
datenum('220217_124850')
open concatenateDats.m
names2sort = cellfun(@(X) str2num(subdirs),recordingnames,'UniformOutput',false);
names2sort = cellfun(@(X) str2num(X(end-5:end),subdirs,'UniformOutput',false);
names2sort = cellfun(@(X) str2num(X(end-5:end),subdirs,'UniformOutput',false));
names2sort = cellfun(@(X) str2num(X(end-5:end)),subdirs,'UniformOutput',false));
names2sort = cellfun(@(X) str2num(X(end-5:end)),subdirs,'UniformOutput',false);
names2sort = cellfun(@(X) str2num(X(end-5:end)),subdirs,'UniformOutput',false);
names2sort = cell2mat(names2sort)
names2sort = cellfun(@(X) str2num(X(end-5:end)),subdirs,'UniformOutput',false);
names2sort = cell2mat(names2sort);
[~,I] = sort(names2sort);
subdirs = subdirs(I);
subdirs
session = sessionTemplate(basepath,'showGUI',true);
concatenateDats(basepath,0,1);
session = sessionTemplate(basepath,'showGUI',true);
digitalInp = getDigitalIn('all','fs',session.extracellular.sr);
dbcont
load('digitalIn.events.mat')
fullfile(basepath,[dayName '.xml'])
[~,dayName] = fileparts(basepath);
fullfile(basepath,[dayName '.xml'])
basename = pwd
'N:\Praveen\SocialBehavior\SM3\day1\day1.xml'
fullfile(basepath,[dayName '.xml'])
clear all
basepath = pwd
45
clear all
session = sessionTemplate(basepath,'showGUI',true);
save(fullfile(basepath,[basename, '.session.mat']),'session');
digitalInp = getDigitalIn('all','fs',session.extracellular.sr);
raly
load('digitalIn.events.mat')
data = digitalIn.its;
data = digitalIn.ints;
for i=1:length(data),data{i} = data{i}'; data{i}(:,end+1) = i; ned
for i=1:length(data),data{i} = data{i}'; data{i}(:,end+1) = i; end
data = digitalIn.ints';
for i=1:length(data),data{i} = data{i}'; data{i}(:,end+1) = i; end
intervals = sortrows(cell2mat(data));
cell2mat(@(x) sum(diff(x)),digitalIn.ints)
cellfun(@(x) sum(diff(x)),digitalIn.ints)
load(fullfile(basepath,[basename '.MergePoints.events.mat']),'MergePoints')
MergePoints.timestamps;
data = digitalIn.ints';
for i=1:length(data),data{i} = data{i}';end
subdata = data{i}(data{i}(:,1)>subsessions(i,1) & data{i}(:,1)<subsessions(i,2),:);
script_Praveen
subsessions = MergePoints.timestamps;
i=1
j=1
subdata = data{i}(data{i}(:,1)>subsessions(i,1) & data{i}(:,1)<subsessions(i,2),:);
PlotColorMap(clickerTotal)
PlotColorMap(clickerTotal(:,2:end))
bar(clickerTotal)
bar(clickerTotal(:,2:end))
session.epochs
session.epochs.names
session.epochs
bar(clickerTotal(:,2:end));
clear all
script_Praveen
clear all
script_Praveen
%-- 2/22/2022 6:31 PM --%
cd N:\Praveen\uLED\uLED1\day7
load('day7.ripples.events.mat')
r = [ripples.timestamps(:,1) ripples.peaks ripples.timestamps(:,2) ripples.RipMax ripples.SwMax];
X = BatchReturn('OJR_longTraining.batch');
rCell{1} = r;
load('day7.MergePoints.events.mat')
subsessions = MergePoints.timestamps;
pre = r(:,1)<MergePoints.timestamps(3,1);
subsessions = MergePoints.timestamps;
sleep = MergePoints.timestamps(cellfun(@(x) any(strcmp(lower(x),'sleep')),MergePoints.foldernames),:);
cellfun(@(x) any(strcmp(lower(x),'sleep')),MergePoints.foldernames)
sleep = MergePoints.timestamps(cellfun(@(x) any(strfind(lower(x),'sleep')),MergePoints.foldernames),:);
pre = InIntervals(sleep(1,:));
sleep = ConsolidateIntervals(MergePoints.timestamps(cellfun(@(x) any(strfind(lower(x),'sleep')),MergePoints.foldernames),:));
pre = InIntervals(r,sleep(1,:));
post = InIntervals(r,sleep(end,:));
anovabar(diff(r,[],2),[-pre post]);
anovabar(diff(r,[],2),[-pre+post]);
anovabar(diff(r(:,[1 3]),[],2),[-pre+post]);
anovabar(diff(r(:,[1 3]),[],2)>0.1,[-pre+post]);
sum(diff(r(:,[1 3]),[],2)>0.1)
id = (-pre+post)+2;
mmax(id)
[num2str(sum(diff(r(id==j,[1 3]),[],2)>0.1)) '/' num2str(id==j)]
sum(diff(r(id==j,[1 3]),[],2)>0.1)
tr
34
i=2
basepath = X{i};
[~,basename] = fileparts(basepath);
subplot(3,3,i);title(basepath);
load(fullfile(basepath,[basename '.MergePoints.events.mat']),'MergePoints');
load(fullfile(basepath,[basename '.ripples.events.mat']),'ripples');
r = [ripples.timestamps(:,1) ripples.peaks ripples.timestamps(:,2) ripples.RipMax ripples.SwMax];
sleep = ConsolidateIntervals(MergePoints.timestamps(cellfun(@(x) any(strfind(lower(x),'sleep')),MergePoints.foldernames),:));
MergePoints.foldernames
sleep
raly
cd(X{5})
cd(X{4})
load('day7.channelinfo.ripples.mat')
rippleChannels.Ripple_Channel = 32; rippleChannels.Sharpwave_Channel = 10;
rippleChannels.method = 'Manual scoring on neuroscope by Raly Feb 22 2022';
ripples = DetectSWR([rippleChannels.Ripple_Channel rippleChannels.Sharpwave_Channel 76+1],'basepath',basepath,'saveMat',true,'forceDetect',true,'check',true);
% First, remove points that occur in noisy lfp periods
tl = featureTs(:)/1250; % timestamps for lfp
[clean,bad,badIntervals] = CleanLFP([tl,double(lfp(:,2))]);
smoothed = Smooth(abs(swDiff),1250*10);
noisyPeriods = ConsolidateIntervals(tl(FindInterval(smoothed>200)),'epsilon',1);
sw = interp1(tl,lfpLow(:,2),t);
bad = InIntervals(t,badIntervals) | InIntervals(t,noisyPeriods) | sw>0;
idx1 = idx1 & ~bad; idx2 = idx2 & ~bad; % remove bad points from "idx1" and "idx2"
figure(1); % plot a clean figure without these points and chech indivudual ripples
clf
plot(swDiffAll(idx2),ripPowerAll(idx2),'k.','markersize',1); hold on
plot(swDiffAll(idx1),ripPowerAll(idx1),'r.','markersize',1);
while true % click "Ctrl+C" to exit when you feel you're done
figure(2);
colors = Bright(1000);
[x,y,button] = ginput(1);
i = findmin(abs(x-swDiffAll)+abs(y-ripPowerAll));
t = featureTs/1250;
tl = (1:length(lfp))'/1250;
figure(3); clf; interval = t(i) + [-1 1]*0.5; in = InIntervals(tl,interval); plot(tl(in) - t(i),lfp(in,:)); xlabel(num2str(t(i)));
legend('ripple channel','SW channel','noise channel');
%     x = input( prompt )
score = str2double(input('good(1)/bad(0)?','s'));
figure(2);
if ~isnan(score) && score>0.1
hold on
plot(swDiffAll(i),ripPowerAll(i),'o','linewidth',2,'color',colors(findmin(abs(score-linspace(0,1,1000)')),:));
else
plot(swDiffAll(i),ripPowerAll(i),'x','linewidth',2,'color',colors(findmin(abs(score-linspace(0,1,1000)')),:));
end
end
tl = (1:length(lfp))'/1250; % timestamps for lfp
[clean,bad,badIntervals] = CleanLFP([tl,double(lfp(:,2))]);
smoothed = Smooth(abs(swDiff),1250*10);
noisyPeriods = ConsolidateIntervals(tl(FindInterval(smoothed>200)),'epsilon',1);
sw = interp1(tl,lfpLow(:,2),t);
bad = InIntervals(t,badIntervals) | InIntervals(t,noisyPeriods) | sw>0;
idx1 = idx1 & ~bad; idx2 = idx2 & ~bad; % remove bad points from "idx1" and "idx2"
figure(1); % plot a clean figure without these points and chech indivudual ripples
clf
plot(swDiffAll(idx2),ripPowerAll(idx2),'k.','markersize',1); hold on
plot(swDiffAll(idx1),ripPowerAll(idx1),'r.','markersize',1);
while true % click "Ctrl+C" to exit when you feel you're done
figure(2);
colors = Bright(1000);
[x,y,button] = ginput(1);
i = findmin(abs(x-swDiffAll)+abs(y-ripPowerAll));
t = featureTs/1250;
tl = (1:length(lfp))'/1250;
figure(3); clf; interval = t(i) + [-1 1]*0.5; in = InIntervals(tl,interval); plot(tl(in) - t(i),lfp(in,:)); xlabel(num2str(t(i)));
legend('ripple channel','SW channel','noise channel');
%     x = input( prompt )
score = str2double(input('good(1)/bad(0)?','s'));
figure(2);
if ~isnan(score) && score>0.1
hold on
plot(swDiffAll(i),ripPowerAll(i),'o','linewidth',2,'color',colors(findmin(abs(score-linspace(0,1,1000)')),:));
else
plot(swDiffAll(i),ripPowerAll(i),'x','linewidth',2,'color',colors(findmin(abs(score-linspace(0,1,1000)')),:));
end
end
ripples = DetectSWR([rippleChannels.Ripple_Channel rippleChannels.Sharpwave_Channel 76+1],'basepath',basepath,'saveMat',true,'forceDetect',true,'check',true);
which ConsolidateIntervals
fmat
% First, remove points that occur in noisy lfp periods
tl = (1:length(lfp))'/1250; % timestamps for lfp
[clean,bad,badIntervals] = CleanLFP([tl,double(lfp(:,2))]);
smoothed = Smooth(abs(swDiff),1250*10);
noisyPeriods = ConsolidateIntervals(tl(FindInterval(smoothed>200)),'epsilon',1);
noisyPeriods = ConsolidateIntervals(tl(FindInterval(smoothed>200)),'epsilon',1);
tl(FindInterval(smoothed>200))
noisyPeriods = ConsolidateIntervals(tl(FindInterval(smoothed>200)),'epsilon',1);
sw = interp1(tl,lfpLow(:,2),t);
bad = InIntervals(t,badIntervals) | InIntervals(t,noisyPeriods) | sw>0;
idx1 = idx1 & ~bad; idx2 = idx2 & ~bad; % remove bad points from "idx1" and "idx2"
figure(1); % plot a clean figure without these points and chech indivudual ripples
clf
plot(swDiffAll(idx2),ripPowerAll(idx2),'k.','markersize',1); hold on
plot(swDiffAll(idx1),ripPowerAll(idx1),'r.','markersize',1);
while true % click "Ctrl+C" to exit when you feel you're done
figure(2);
colors = Bright(1000);
[x,y,button] = ginput(1);
i = findmin(abs(x-swDiffAll)+abs(y-ripPowerAll));
t = featureTs/1250;
tl = (1:length(lfp))'/1250;
figure(3); clf; interval = t(i) + [-1 1]*0.5; in = InIntervals(tl,interval); plot(tl(in) - t(i),lfp(in,:)); xlabel(num2str(t(i)));
legend('ripple channel','SW channel','noise channel');
%     x = input( prompt )
score = str2double(input('good(1)/bad(0)?','s'));
figure(2);
if ~isnan(score) && score>0.1
hold on
plot(swDiffAll(i),ripPowerAll(i),'o','linewidth',2,'color',colors(findmin(abs(score-linspace(0,1,1000)')),:));
else
plot(swDiffAll(i),ripPowerAll(i),'x','linewidth',2,'color',colors(findmin(abs(score-linspace(0,1,1000)')),:));
end
end
ripples = DetectSWR([rippleChannels.Ripple_Channel rippleChannels.Sharpwave_Channel 76+1],'basepath',basepath,'saveMat',true,'forceDetect',true,'check',true);
% First, remove points that occur in noisy lfp periods
tl = (1:length(lfp))'/1250; % timestamps for lfp
[clean,bad,badIntervals] = CleanLFP([tl,double(lfp(:,2))]);
smoothed = Smooth(abs(swDiff),1250*10);
noisyPeriods = ConsolidateIntervals(tl(FindInterval(smoothed>200)),'epsilon',1);
try sw = interp1(tl,lfpLow(:,2),t);
bad = InIntervals(t,badIntervals) | InIntervals(t,noisyPeriods) | sw>0;
end
idx1 = idx1 & ~bad; idx2 = idx2 & ~bad; % remove bad points from "idx1" and "idx2"
figure(1); % plot a clean figure without these points and chech indivudual ripples
clf
plot(swDiffAll(idx2),ripPowerAll(idx2),'k.','markersize',1); hold on
plot(swDiffAll(idx1),ripPowerAll(idx1),'r.','markersize',1);
while true % click "Ctrl+C" to exit when you feel you're done
figure(2);
colors = Bright(1000);
[x,y,button] = ginput(1);
i = findmin(abs(x-swDiffAll)+abs(y-ripPowerAll));
t = featureTs/1250;
tl = (1:length(lfp))'/1250;
figure(3); clf; interval = t(i) + [-1 1]*0.5; in = InIntervals(tl,interval); plot(tl(in) - t(i),lfp(in,:)); xlabel(num2str(t(i)));
legend('ripple channel','SW channel','noise channel');
%     x = input( prompt )
score = str2double(input('good(1)/bad(0)?','s'));
figure(2);
if ~isnan(score) && score>0.1
hold on
plot(swDiffAll(i),ripPowerAll(i),'o','linewidth',2,'color',colors(findmin(abs(score-linspace(0,1,1000)')),:));
else
plot(swDiffAll(i),ripPowerAll(i),'x','linewidth',2,'color',colors(findmin(abs(score-linspace(0,1,1000)')),:));
end
end
ripples = DetectSWR([rippleChannels.Ripple_Channel rippleChannels.Sharpwave_Channel 76+1],'basepath',basepath,'saveMat',true,'forceDetect',true,'check',true);
% First, remove points that occur in noisy lfp periods
tl = (1:length(lfp))'/1250; % timestamps for lfp
bad = false(size(tl));
try % optionally, remove periods of noisy LFP (large deflections)
[clean,bad,badIntervals] = CleanLFP([tl,double(lfp(:,2))]);
bad = bad | InIntervals(t,badIntervals);
end
try % optionally, remove periods when your channels are diverging abnormally for a long time (a channel went dead for some time)
smoothed = Smooth(abs(swDiff),1250*10);
noisyPeriods = ConsolidateIntervals(tl(FindInterval(smoothed>200)),'epsilon',1);
bad = bad | InIntervals(t,noisyPeriods);
end
idx1 = idx1 & ~bad; idx2 = idx2 & ~bad; % remove bad points from "idx1" and "idx2"
sum(InIntervals(t,badIntervals))
t = featureTs/1250; % timestamps for candidate ripple events
try % optionally, remove periods of noisy LFP (large deflections)
[clean,bad,badIntervals] = CleanLFP([tl,double(lfp(:,2))]);
bad = bad | InIntervals(t,badIntervals);
end
try % optionally, remove periods when your channels are diverging abnormally for a long time (a channel went dead for some time)
smoothed = Smooth(abs(swDiff),1250*10);
noisyPeriods = ConsolidateIntervals(tl(FindInterval(smoothed>200)),'epsilon',1);
bad = bad | InIntervals(t,noisyPeriods);
end
try % optionally, remove ripples in which the sharp wave channel was positive
% Only do this if you have a sharp-wave channel what's sufficiently deep for the sharp wave to be negative
sw = interp1(tl,lfpLow(:,2),t);
bad = bad | sw>0;
end
Portion(InIntervals(t,badIntervals))
Portion(InIntervals(t,noisyPeriods))
Portion(sw>0)
clf
plot(tl,smoothed);
PlotHVLines(6761,'v');
plot(tl,swDiff);
in = InIntervals(tl,12619 + [0 1]);
figure; plot(t(in),lfp(in,[1 2]));
plot(tl(in),lfp(in,[1 2]));
subplot(2,1,1);
plot(tl(in),lfp(in,:));
subplot(2,1,2);
plot(tl(in),swDiff(in,:));
hold all
plot(tl(in),lfp(in,2),'r','linewidth',2);
session
basepath
Channels
xlim
dbquit
basepath
rippleChannels3 =ripplesChannels;
rippleChannels3 =rippleChannels;
load('day3.channelinfo.ripples.mat')
rippleChannels3
rippleChannels
rippleChannels.method = 'Manual scoring on neuroscope by Raly Feb 23 2022';
rippleChannels.Ripple_Channel = 60+1; rippleChannels.Sharpwave_Channel = 34+1; rippleChannels.Noise_Channel = 76+1;
save('day10.channelinfo.ripples.mat','rippleChannels');
close all
ripples = DetectSWR([rippleChannels.Ripple_Channel rippleChannels.Sharpwave_Channel rippleChannels.Noise_Channel],'basepath',basepath,'saveMat',true,'forceDetect',true,'check',true);
edit preprocessV1
% First, remove points that occur in noisy lfp periods
tl = (1:length(lfp))'/1250; % timestamps for lfp
t = featureTs/1250; % timestamps for candidate ripple events
bad = false(size(tl));
try % optionally, remove periods of noisy LFP (large deflections)
[clean,bad,badIntervals] = CleanLFP([tl,double(lfp(:,2))]);
bad = bad | InIntervals(t,badIntervals);
end
try % optionally, remove periods when your channels are diverging abnormally for a long time (a channel went dead for some time)
smoothed = Smooth(abs(swDiff),1250*10);
noisyPeriods = ConsolidateIntervals(tl(FindInterval(smoothed>200)),'epsilon',1);
bad = bad | InIntervals(t,noisyPeriods);
end
try % optionally, remove ripples in which the sharp wave channel was positive
% Only do this if you have a sharp-wave channel what's sufficiently deep for the sharp wave to be negative
sw = interp1(tl,lfpLow(:,2),t);
bad = bad | sw>0;
end
idx1 = idx1 & ~bad; idx2 = idx2 & ~bad; % remove bad points from "idx1" and "idx2"
Portion(bad)
Portion(InIntervals(t,badIntervals))
Portion(InIntervals(t,noisyPeriods))
Portion(sw>0)
Portion(InIntervals(t,badIntervals))
Portion(InIntervals(t,noisyPeriods))
Portion(sw>0)
in = InIntervals(tl,12619 + [0 1]);
figure; plot(tl(in),lfp(in,:));
SideAxes(gca,'top',0.5);
plot(tl(in),sw(in));
plot(tl(in),lfpLow(in,2));
these = find(InIntervals(t,interval));
these = find(InIntervals(t,12619 + [0 1]));
open these
PlotHVLines(t(these),'v','k');
plot(tl(in),lfpLow(in,:));
raly
edit DetectSWR
Nchan
ripples = DetectSWR([rippleChannels.Ripple_Channel rippleChannels.Sharpwave_Channel rippleChannels.Noise_Channel],'basepath',basepath,'saveMat',true,'forceDetect',true,'check',true);
which DetectSWR
ripples = DetectSWR([rippleChannels.Ripple_Channel rippleChannels.Sharpwave_Channel rippleChannels.Noise_Channel],'basepath',basepath,'saveMat',true,'forceDetect',true,'check',true);
% First, remove points that occur in noisy lfp periods
tl = (1:length(lfp))'/1250; % timestamps for lfp
t = featureTs/1250; % timestamps for candidate ripple events
bad = false(size(tl));
try % optionally, remove periods of noisy LFP (large deflections)
[clean,bad,badIntervals] = CleanLFP([tl,double(lfp(:,2))]);
bad = bad | InIntervals(t,badIntervals);
end
try % optionally, remove periods when your channels are diverging abnormally for a long time (a channel went dead for some time)
smoothed = Smooth(abs(swDiff),1250*10);
noisyPeriods = ConsolidateIntervals(tl(FindInterval(smoothed>200)),'epsilon',1);
bad = bad | InIntervals(t,noisyPeriods);
end
try % optionally, remove ripples in which the sharp wave channel was positive
% Only do this if you have a sharp-wave channel what's sufficiently deep for the sharp wave to be negative
sw = interp1(tl,lfpLow(:,2),t);
bad = bad | sw>0;
end
idx1 = idx1 & ~bad; idx2 = idx2 & ~bad; % remove bad points fr
Portion(InIntervals(t,badIntervals))
Portion(InIntervals(t,noisyPeriods))
Portion(sw>0)
in = InIntervals(tl,12619 + [0 1]);
figure; plot(t(in),lfp(in,[1 2]));
plot(tl(in),lfp(in,[1 2]));
subplot(2,1,1);
plot(tl(in),lfp(in,:));
subplot(2,1,2);
plot(tl(in),swDiff(in,:));
subplot(2,1,1);
plot(tl(in),lfp(in,:));
subplot(2,1,2);
plot(tl(in),swDiff(in,:));
these = find(InIntervals(t,12619 + [0 1]));
PlotHVLines(t(these),'v','k');
figure; plot(tl,smoothed);
PlotHVLines(t(these),'v','k');
PlotXY(Restrict([tl sw],xlim));
PlotXY(Restrict([tl abs(swDiff)],xlim));
PlotXY(Restrict([tl Smooth(abs(swDiff),1250*10)],xlim));
interval =  12619 + [0 1]);
interval =  12619 + [0 1];
xlim(interval)
PlotHVLines(t(these),'v','y');
max(smoothed)
hm = tl(FindInterval(smoothed>500));
PlotIntervals(hm,'color','r');
% run a basic k-means clustering to partition SWR from ~SWR
[idx,Ctrs] = kmeans([swDiffAll(:) ripPowerAll(:)],2);
% we know there will be far fewer SWR than not, so let's assign idx1 to
% the cluster with the fewest members i.e. SWR
if sum(idx == 1) > sum(idx == 2)
idx1      = logical(idx == 1);
idx2      = logical(idx == 2);
% reverse labels
temp      = idx1;
idx1      = idx2;
idx2      = temp;
clear temp
else
idx1      = logical(idx == 1);
idx2      = logical(idx == 2);
end
bad = false(size(tl));
Portion(InIntervals(t,badIntervals))
Portion(InIntervals(t,noisyPeriods))
Portion(sw>0)
open badIntervals
figure; PlotXY([tl lfp(:,2)]);
PlotIntervals(badIntervals,'color','r');
find(tl>2220,1)
figure; PlotXY([tl double(lfp(:,2))]);
PlotIntervals(badIntervals,'color','r');
[clean,bad,badIntervals] = CleanLFP([tl,double(lfp(:,2))],'thresholds',[5 Inf]);
PlotIntervals(badIntervals,'color','y');
[clean,bad,badIntervals] = CleanLFP([tl,double(lfp(:,2))],'thresholds',[10 Inf]);
PlotIntervals(badIntervals,'color','k');
clf
figure; PlotXY([tl double(lfp(:,2))]);
PlotIntervals(badIntervals,'color','k');
[ii,w] = InIntervals(tl,badIntervals);
s = Accumulate(w(ii),abs(lfp(ii,2)),'mode','max');
s = Accumulate(w(ii),abs(double(lfp(ii,2))),'mode','max');
open s
findmin(s)
badIntervals(50,:)
xlim(badIntervals(50,:)+[-1 1]*2);
[1.5528    1.5532] - 15529
badIntervals(50,:) - 15529
badIntervals(50,:) - 15527
Channels
findmin(s,5)
addComma(badIntervals(35,:))
addComma(badIntervals(35,1))
addComma(badIntervals(35,2))
dur(badIntervals)
Portion(InIntervals(t,badIntervals))
Portion(InIntervals(t,noisyPeriods))
Portion(sw>0)
bad = false(size(tl));
bad = bad | InIntervals(t,badIntervals);
size(bad)
bad = false(size(t));
bad = bad | InIntervals(t,badIntervals);
Portion(bad)
idx1 = idx1 & ~bad; idx2 = idx2 & ~bad; % remove bad points from "idx1" and "idx2"
close all
figure(1); % plot a clean figure without these points and chech indivudual ripples
clf
plot(swDiffAll(idx2),ripPowerAll(idx2),'k.','markersize',1); hold on
plot(swDiffAll(idx1),ripPowerAll(idx1),'r.','markersize',1);
while true % click "Ctrl+C" to exit when you feel you're done
figure(2);
colors = Bright(1000);
[x,y,button] = ginput(1);
i = findmin(abs(x-swDiffAll)+abs(y-ripPowerAll));
t = featureTs/1250;
tl = (1:length(lfp))'/1250;
figure(3); clf; interval = t(i) + [-1 1]*0.5; in = InIntervals(tl,interval); plot(tl(in) - t(i),lfp(in,:)); xlabel(num2str(t(i)));
legend('ripple channel','SW channel','noise channel');
%     x = input( prompt )
score = str2double(input('good(1)/bad(0)?','s'));
figure(2);
if ~isnan(score) && score>0.1
hold on
plot(swDiffAll(i),ripPowerAll(i),'o','linewidth',2,'color',colors(findmin(abs(score-linspace(0,1,1000)')),:));
else
plot(swDiffAll(i),ripPowerAll(i),'x','linewidth',2,'color',colors(findmin(abs(score-linspace(0,1,1000)')),:));
end
end
figure(1); % plot a clean figure without these points and check indivudual ripples
clf
plot(swDiffAll(idx2),ripPowerAll(idx2),'k.','markersize',1); hold on
plot(swDiffAll(idx1),ripPowerAll(idx1),'r.','markersize',1);
while true % click "Ctrl+C" to exit when you feel you're done
figure(1);
colors = Bright(1000);
[x,y,button] = ginput(1);
i = findmin(abs(x-swDiffAll)+abs(y-ripPowerAll));
t = featureTs/1250;
tl = (1:length(lfp))'/1250;
figure(3); clf; interval = t(i) + [-1 1]*0.5; in = InIntervals(tl,interval); plot(tl(in) - t(i),lfp(in,:)); xlabel(num2str(t(i)));
legend('ripple channel','SW channel','noise channel');
%     x = input( prompt )
score = str2double(input('good(1)/bad(0)?','s'));
figure(2);
if ~isnan(score) && score>0.1
hold on
plot(swDiffAll(i),ripPowerAll(i),'o','linewidth',2,'color',colors(findmin(abs(score-linspace(0,1,1000)')),:));
else
plot(swDiffAll(i),ripPowerAll(i),'x','linewidth',2,'color',colors(findmin(abs(score-linspace(0,1,1000)')),:));
end
end
1
while true % click "Ctrl+C" to exit when you feel you're done
figure(1);
colors = Bright(1000);
[x,y,button] = ginput(1);
i = findmin(abs(x-swDiffAll)+abs(y-ripPowerAll));
t = featureTs/1250;
tl = (1:length(lfp))'/1250;
figure(3); clf; interval = t(i) + [-1 1]*0.5; in = InIntervals(tl,interval); plot(tl(in) - t(i),lfp(in,:)); xlabel(num2str(t(i)));
legend('ripple channel','SW channel','noise channel');
%     x = input( prompt )
score = str2double(input('good(1)/bad(0)?','s'));
figure(1);
if ~isnan(score) && score>0.1
hold on
plot(swDiffAll(i),ripPowerAll(i),'o','linewidth',2,'color',colors(findmin(abs(score-linspace(0,1,1000)')),:));
else
plot(swDiffAll(i),ripPowerAll(i),'x','linewidth',2,'color',colors(findmin(abs(score-linspace(0,1,1000)')),:));
end
end
0
0.05
1
0.8
0
0.05
1
0.5
1
0.4
0
1
0.1
1
0
1
0
while true % click "Ctrl+C" to exit when you feel you're done
figure(1);
colors = Bright(1000);
[x,y,button] = ginput(1);
i = findmin(abs(x-swDiffAll)+abs(y-ripPowerAll));
t = featureTs/1250;
tl = (1:length(lfp))'/1250;
figure(3); clf; interval = t(i) + [-1 1]*0.5; in = InIntervals(tl,interval); plot(tl(in) - t(i),lfp(in,:)); xlabel(num2str(t(i)));
legend('ripple channel','SW channel','noise channel');
%     x = input( prompt )
score = str2double(input('good(1)/bad(0)?','s'));
figure(1);
if ~isnan(score) && score>0.1
hold on
plot(swDiffAll(i),ripPowerAll(i),'o','linewidth',2,'color',colors(findmin(abs(score-linspace(0,1,1000)')),:));
else
plot(swDiffAll(i),ripPowerAll(i),'x','linewidth',2,'color',colors(findmin(abs(score-linspace(0,1,1000)')),:));
end
end
0
0.4
1
0
0.4
1
0.4
0
0.3
1
0.5
0.7
1
0.8
0
1
0.2
0.5
while true % click "Ctrl+C" to exit when you feel you're done
figure(1);
colors = Bright(1000);
[x,y,button] = ginput(1);
i = findmin(abs(x-swDiffAll)+abs(y-ripPowerAll));
t = featureTs/1250;
tl = (1:length(lfp))'/1250;
figure(3); clf; interval = t(i) + [-1 1]*0.5; in = InIntervals(tl,interval); plot(tl(in) - t(i),lfp(in,:)); xlabel(num2str(t(i)));
legend('ripple channel','SW channel','noise channel');
%     x = input( prompt )
score = str2double(input('good(1)/bad(0)?','s'));
figure(1);
if ~isnan(score) && score>0.1
hold on
plot(swDiffAll(i),ripPowerAll(i),'o','linewidth',2,'color',colors(findmin(abs(score-linspace(0,1,1000)')),:));
else
plot(swDiffAll(i),ripPowerAll(i),'x','linewidth',2,'color',colors(findmin(abs(score-linspace(0,1,1000)')),:));
end
end
1
0
1
0
1
0.1
1
0
0.9
0.8
0.95
0
0.2
0
0.2
1
0.5
1
0.3
0
0.5
0
1
0
i
ts(89993)
addComma(t(89993))
addComma(t(89994))
scores = nan(size(ripPowerAll,1),1);
while true % click "Ctrl+C" to exit when you feel you're done
figure(1);
colors = Bright(1000);
[x,y,button] = ginput(1);
i = findmin(abs(x-swDiffAll)+abs(y-ripPowerAll));
t = featureTs/1250;
tl = (1:length(lfp))'/1250;
figure(3); clf; interval = t(i) + [-1 1]*0.5; in = InIntervals(tl,interval); plot(tl(in) - t(i),lfp(in,:)); xlabel(num2str(t(i)));
legend('ripple channel','SW channel','noise channel');
%     x = input( prompt )
score = str2double(input('good(1)/bad(0)?','s'));
scores(i) = score; % save scores to optionally save
figure(1);
if ~isnan(score) && score>0.1
hold on
plot(swDiffAll(i),ripPowerAll(i),'o','linewidth',2,'color',colors(findmin(abs(score-linspace(0,1,1000)')),:));
else
plot(swDiffAll(i),ripPowerAll(i),'x','linewidth',2,'color',colors(findmin(abs(score-linspace(0,1,1000)')),:));
end
end
0.1
0.4
0
1
0
0.3
1
0.05
addComma(t(i))
addComma(t(i-1))
score
score = 0.8
score = 1;
score = 0.9;
scores(i) = score; % save scores to optionally save
figure(1);
if ~isnan(score) && score>0.1
hold on
plot(swDiffAll(i),ripPowerAll(i),'o','linewidth',2,'color',colors(findmin(abs(score-linspace(0,1,1000)')),:));
else
plot(swDiffAll(i),ripPowerAll(i),'x','linewidth',2,'color',colors(findmin(abs(score-linspace(0,1,1000)')),:));
end
plot(swDiffAll(i),ripPowerAll(i),'x','linewidth',2,'color','w');
plot(swDiffAll(i),ripPowerAll(i),'o','linewidth',2,'color',colors(findmin(abs(score-linspace(0,1,1000)')),:));
while true % click "Ctrl+C" to exit when you feel you're done
figure(1);
colors = Bright(1000);
[x,y,button] = ginput(1);
i = findmin(abs(x-swDiffAll)+abs(y-ripPowerAll));
t = featureTs/1250;
tl = (1:length(lfp))'/1250;
figure(3); clf; interval = t(i) + [-1 1]*0.5; in = InIntervals(tl,interval); plot(tl(in) - t(i),lfp(in,:)); xlabel(num2str(t(i)));
legend('ripple channel','SW channel','noise channel');
%     x = input( prompt )
score = str2double(input('good(1)/bad(0)?','s'));
scores(i) = score; % save scores to optionally save
figure(1);
if ~isnan(score) && score>0.1
hold on
plot(swDiffAll(i),ripPowerAll(i),'o','linewidth',2,'color',colors(findmin(abs(score-linspace(0,1,1000)')),:));
else
plot(swDiffAll(i),ripPowerAll(i),'x','linewidth',2,'color',colors(findmin(abs(score-linspace(0,1,1000)')),:));
end
end
1
0.6
1
0.6
1
addComma(t(i))
addComma(t(i+1))
score = str2double(input('good(1)/bad(0)?','s'));
scores(i) = score; % save scores to optionally save
figure(1);
if ~isnan(score) && score>0.1
hold on
plot(swDiffAll(i),ripPowerAll(i),'o','linewidth',2,'color',colors(findmin(abs(score-linspace(0,1,1000)')),:));
else
plot(swDiffAll(i),ripPowerAll(i),'x','linewidth',2,'color',colors(findmin(abs(score-linspace(0,1,1000)')),:));
end
0
while true % click "Ctrl+C" to exit when you feel you're done
figure(1);
colors = Bright(1000);
[x,y,button] = ginput(1);
i = findmin(abs(x-swDiffAll)+abs(y-ripPowerAll));
t = featureTs/1250;
tl = (1:length(lfp))'/1250;
figure(3); clf; interval = t(i) + [-1 1]*0.5; in = InIntervals(tl,interval); plot(tl(in) - t(i),lfp(in,:)); xlabel(num2str(t(i)));
legend('ripple channel','SW channel','noise channel');
%     x = input( prompt )
score = str2double(input('good(1)/bad(0)?','s'));
scores(i) = score; % save scores to optionally save
figure(1);
if ~isnan(score) && score>0.1
hold on
plot(swDiffAll(i),ripPowerAll(i),'o','linewidth',2,'color',colors(findmin(abs(score-linspace(0,1,1000)')),:));
else
plot(swDiffAll(i),ripPowerAll(i),'x','linewidth',2,'color',colors(findmin(abs(score-linspace(0,1,1000)')),:));
end
end
0.8
addComma(t(i+1))
addComma(t(i))
score = str2double(input('good(1)/bad(0)?','s'));
scores(i) = score; % save scores to optionally save
figure(1);
if ~isnan(score) && score>0.1
hold on
plot(swDiffAll(i),ripPowerAll(i),'o','linewidth',2,'color',colors(findmin(abs(score-linspace(0,1,1000)')),:));
else
plot(swDiffAll(i),ripPowerAll(i),'x','linewidth',2,'color',colors(findmin(abs(score-linspace(0,1,1000)')),:));
end
while true % click "Ctrl+C" to exit when you feel you're done
figure(1);
colors = Bright(1000);
[x,y,button] = ginput(1);
i = findmin(abs(x-swDiffAll)+abs(y-ripPowerAll));
t = featureTs/1250;
tl = (1:length(lfp))'/1250;
figure(3); clf; interval = t(i) + [-1 1]*0.5; in = InIntervals(tl,interval); plot(tl(in) - t(i),lfp(in,:)); xlabel(num2str(t(i)));
legend('ripple channel','SW channel','noise channel');
%     x = input( prompt )
score = str2double(input('good(1)/bad(0)?','s'));
scores(i) = score; % save scores to optionally save
figure(1);
if ~isnan(score) && score>0.1
hold on
plot(swDiffAll(i),ripPowerAll(i),'o','linewidth',2,'color',colors(findmin(abs(score-linspace(0,1,1000)')),:));
else
plot(swDiffAll(i),ripPowerAll(i),'x','linewidth',2,'color',colors(findmin(abs(score-linspace(0,1,1000)')),:));
end
end
0.95
0
0.8
1
0
0.5
1
0.5
0.4
1
0
1
0
1
0
1
addComma(t(i+`1))
addComma(t(i+1))
while true % click "Ctrl+C" to exit when you feel you're done
figure(1);
colors = Bright(1000);
[x,y,button] = ginput(1);
i = findmin(abs(x-swDiffAll)+abs(y-ripPowerAll));
t = featureTs/1250;
tl = (1:length(lfp))'/1250;
figure(3); clf; interval = t(i) + [-1 1]*0.5; in = InIntervals(tl,interval); plot(tl(in) - t(i),lfp(in,:)); xlabel(num2str(t(i)));
legend('ripple channel','SW channel','noise channel');
%     x = input( prompt )
score = str2double(input('good(1)/bad(0)?','s'));
scores(i) = score; % save scores to optionally save
figure(1);
if ~isnan(score) && score>0.1
hold on
plot(swDiffAll(i),ripPowerAll(i),'o','linewidth',2,'color',colors(findmin(abs(score-linspace(0,1,1000)')),:));
else
plot(swDiffAll(i),ripPowerAll(i),'x','linewidth',2,'color',colors(findmin(abs(score-linspace(0,1,1000)')),:));
end
end
0
1
0.4
i
addComma(t(i))
remember = (t(i))
interval = t(i)+[-0.5 0.5];
in = InIntervals(tl,interval);
figure; plot(tl(in),double(lfp(in,:)));
subplot(3,1,1);
plot(tl(in),double(lfp(in,:)));
subplot(3,1,2);
plot(tl(in),double(swDiffAll(in)));
plot(tl(in),double(swDiff(in)));
PlotHVLines(remember,'v','k--');
PlotHVLines([t(i) t(i+1)],'v','k--');
t(i+1)
addComma(t(i+1))
addComma(t(i))
PlotHVLines([t(i) t(i+1)],'v','k--');
subplot(3,1,3);
plot(tl(in),double(ripPower0(in)));
PlotHVLines([t(i) t(i+1)],'v','k--');
figure(2);
clf; interval = t(i) + [-1 1]*0.5; in = InIntervals(tl,interval);
subplot(3,1,1);
plot(tl(in) - t(i),lfp(in,:)); xlabel(num2str(t(i)));
PlotHVLines(Restrict(t,interval),'v','k--');
legend('ripple channel','SW channel','noise channel');
subplot(3,1,2);
plot(tl(in),double(swDiff(in)));
PlotHVLines(Restrict(t,interval),'v','k--');
subplot(3,1,3);
plot(tl(in),double(ripPower0(in)));
PlotHVLines([t(i) t(i+1)],'v','k--');
figure(2);
clf; interval = t(i) + [-1 1]*0.5; in = InIntervals(tl,interval);
subplot(3,1,1);
plot(tl(in) - t(i),lfp(in,:)); xlabel(num2str(t(i)));
PlotHVLines(Restrict(t,interval) - t(i),'v','k--');
legend('ripple channel','SW channel','noise channel');
subplot(3,1,2);
plot(tl(in)PlotHVLines(Restrict(t,interval) - t(i),'v','k--');,double(swDiff(in)));
PlotHVLines(Restrict(t,interval) - t(i),'v','k--');
subplot(3,1,3);
plot(tl(in)PlotHVLines(Restrict(t,interval) - t(i),'v','k--');,double(ripPower0(in)));
PlotHVLines(Restrict(t,interval) - t(i),'v','k--');
figure(2);
clf; interval = t(i) + [-1 1]*0.5; in = InIntervals(tl,interval);
subplot(3,1,1);
plot(tl(in) - t(i),lfp(in,:)); xlabel(num2str(t(i)));
PlotHVLines(Restrict(t,interval) - t(i),'v','k--');
legend('ripple channel','SW channel','noise channel');
subplot(3,1,2);
plot(tl(in) - t(i),double(swDiff(in)));
PlotHVLines(Restrict(t,interval) - t(i),'v','k--');
subplot(3,1,3);
plot(tl(in) - t(i),double(ripPower0(in)));
PlotHVLines(Restrict(t,interval) - t(i),'v','k--');
while true % click "Ctrl+C" to exit when you feel you're done
figure(1);
colors = Bright(1000);
[x,y,button] = ginput(1);
i = findmin(abs(x-swDiffAll)+abs(y-ripPowerAll));
t = featureTs/1250;
tl = (1:length(lfp))'/1250;
figure(2);
clf; interval = t(i) + [-1 1]*0.5; in = InIntervals(tl,interval);
subplot(3,1,1);
plot(tl(in) - t(i),lfp(in,:)); xlabel(num2str(t(i)));
PlotHVLines(Restrict(t,interval) - t(i),'v','k--');
legend('ripple channel','SW channel','noise channel');
subplot(3,1,2);
plot(tl(in) - t(i),double(swDiff(in)));
PlotHVLines(Restrict(t,interval) - t(i),'v','k--');
subplot(3,1,3);
plot(tl(in) - t(i),double(ripPower0(in)));
PlotHVLines(Restrict(t,interval) - t(i),'v','k--');
%     x = input( prompt )
score = str2double(input('good(1)/bad(0)?','s'));
scores(i) = score; % save scores to optionally save
figure(1);
if ~isnan(score) && score>0.1
hold on
plot(swDiffAll(i),ripPowerAll(i),'o','linewidth',2,'color',colors(findmin(abs(score-linspace(0,1,1000)')),:));
else
plot(swDiffAll(i),ripPowerAll(i),'x','linewidth',2,'color',colors(findmin(abs(score-linspace(0,1,1000)')),:));
end
end
0
1
0
ripPowerAll(i)
interp1(tl,ripPower0,t(i))
nancorr(ripPowerAll,interp1(tl,ripPower0,t))
figure; plot(ripPowerAll,interp1(tl,ripPower0,t))
figure; plot(ripPowerAll,interp1(tl,ripPower0,t),'.')
plot(ripPowerAll-interp1(tl,ripPower0,t),'.')
sum(~isnan(scores))
figure;  plot(swDiffAll,ripPowerAll,'k.','markersize',1); hold on
disp(['Please review the figure and change the idx1 (ripples) ' newline 'and idx2 (calm non-ripple periods) groups.' newline ...
'Note that noisy periods need not participate in either group.']);
these = find(~isnan(scores));
i
for j=1:length(these), i=these(j);
score = scores(i);
if ~isnan(score) && score>0.1
hold on
plot(swDiffAll(i),ripPowerAll(i),'o','linewidth',2,'color',colors(findmin(abs(score-linspace(0,1,1000)')),:));
else
plot(swDiffAll(i),ripPowerAll(i),'x','linewidth',2,'color',colors(findmin(abs(score-linspace(0,1,1000)')),:));
end
end
sum(ripPowerAll==interp1(tl,ripPower0,t))
Portion(ripPowerAll==interp1(tl,ripPower0,t))
Portion(ripPowerAll-interp1(tl,ripPower0,t) <0.01)
plot(ripPowerAll-interp1(tl,ripPower0,t),'.','markersize',1)
HalfWinSize
HalfWinSize/1250
[ii,w] = InIntervals(tl,[t-0.1 t+0.1]);
estPower = Accumulate(w(ii),ripPower0(ii),'mode','max');
plot(ripPowerAll-estPower,'.','markersize',1)
plot(ripPowerAll,estPower,'.','markersize',1)
plot(interp1(tl,ripPower0,t),estPower,'.','markersize',1)
size(idx)
size(ids)
ids = bsxfun(@plus,-HalfWinSize:HalfWinSize,featureTs);
open ids
ids = bsxfun(@plus,-HalfWinSize:HalfWinSize,featureTs); ids(ids<1)=1; ids(ids>length(tl)) = length(tl);
size(values)
values = ripPower0(ids);
size(values)
values = max(ripPower0(ids),[],2);
plot(values,estPower,'.','markersize',1)
plot(interp1(tl,ripPower0,t),values,'.','markersize',1)
plot(ripPowerAll,values,'.','markersize',1)
values(66452)
interval
in = InIntervals(tl,interval);
plot(tl(in),values(in));
plot(tl(in),ripPower0(in));
plot(tl(in)-t(i),ripPower0(in));
plot(tl(in)-t(66452),ripPower0(in));
PlotHVLines([-1 1]*HalfWinSize/1250);
WinSize
plot(interp1(tl,ripPower0,t),score);
plot(interp1(tl,ripPower0,t),scores);
clf
plot(interp1(tl,ripPower0,t),scores);
plot(interp1(tl,ripPower0,t),scores,'.');
scatter(ripPowerAll,scores,15,scores,'filled');
scatter(ripPowerAll,interp1(tl,ripPower0,t),15,scores,'filled');
hold all
plot(xlim,xlim,'k');
while true % click "Ctrl+C" to exit when you feel you're done
figure(1);
colors = Bright(1000);
[x,y,button] = ginput(1);
i = findmin(abs(x-swDiffAll)+abs(y-ripPowerAll));
t = featureTs/1250;
tl = (1:length(lfp))'/1250;
figure(2);
clf; interval = t(i) + [-1 1]*0.5; in = InIntervals(tl,interval);
subplot(3,1,1);
plot(tl(in) - t(i),lfp(in,:)); xlabel(num2str(t(i)));
PlotHVLines(Restrict(t,interval) - t(i),'v','k--');
legend('ripple channel','SW channel','noise channel');
subplot(3,1,2);
plot(tl(in) - t(i),double(swDiff(in)));
PlotHVLines(Restrict(t,interval) - t(i),'v','k--');
subplot(3,1,3);
plot(tl(in) - t(i),double(ripPower0(in)));
PlotHVLines(Restrict(t,interval) - t(i),'v','k--');
%     x = input( prompt )
score = str2double(input('good(1)/bad(0)?','s'));
scores(i) = score; % save scores to optionally save
figure(1);
if ~isnan(score) && score>0.1
hold on
plot(swDiffAll(i),ripPowerAll(i),'o','linewidth',2,'color',colors(findmin(abs(score-linspace(0,1,1000)')),:));
else
plot(swDiffAll(i),ripPowerAll(i),'x','linewidth',2,'color',colors(findmin(abs(score-linspace(0,1,1000)')),:));
end
end
1
0
1
i
figure; PETH([tl ripPower0],t);
PETH([tl ripPower0],t,'durations',[-1 1]*0.2,'nBins',501);
PlotHVLines([-1 1]*0.05,'v','k--');
PlotHVLines([-1 1]*0.05,'v','r--');
figure; hist(diff(t),100);
save('ripple_workspace.mat');
close all
ripples = DetectSWR([rippleChannels.Ripple_Channel rippleChannels.Sharpwave_Channel rippleChannels.Noise_Channel],'basepath',basepath,'saveMat',true,'forceDetect',true,'check',true);
% First, remove points that occur in noisy lfp periods
tl = (1:length(lfp))'/1250; % timestamps for lfp
t = featureTs/1250; % timestamps for candidate ripple events
bad = false(size(t));
[clean,bad,badIntervals] = CleanLFP([tl,double(lfp(:,2))],'thresholds',[10 Inf]);
bad = bad | InIntervals(t,badIntervals);
idx1 = idx1 & ~bad; idx2 = idx2 & ~bad; % remove bad poi
figure(1); % plot a clean figure without these points and check indivudual ripples
clf
plot(swDiffAll(idx2),ripPowerAll(idx2),'k.','markersize',1); hold on
plot(swDiffAll(idx1),ripPowerAll(idx1),'r.','markersize',1);
scores = nan(size(ripPowerAll,1),1);
while true % click "Ctrl+C" to exit when you feel you're done
figure(1);
colors = Bright(1000);
[x,y,button] = ginput(1);
i = findmin(abs(x-swDiffAll)+abs(y-ripPowerAll));
t = featureTs/1250;
tl = (1:length(lfp))'/1250;
figure(2);
clf; interval = t(i) + [-1 1]*0.5; in = InIntervals(tl,interval);
subplot(3,1,1);
plot(tl(in) - t(i),lfp(in,:)); xlabel(num2str(t(i)));
PlotHVLines(Restrict(t,interval) - t(i),'v','k--');
legend('ripple channel','SW channel','noise channel');
subplot(3,1,2);
plot(tl(in) - t(i),double(swDiff(in)));
PlotHVLines(Restrict(t,interval) - t(i),'v','k--');
subplot(3,1,3);
plot(tl(in) - t(i),double(ripPower0(in)));
PlotHVLines(Restrict(t,interval) - t(i),'v','k--');
%     x = input( prompt )
score = str2double(input('good(1)/bad(0)?','s'));
scores(i) = score; % save scores to optionally save
figure(1);
if ~isnan(score) && score>0.1
hold on
plot(swDiffAll(i),ripPowerAll(i),'o','linewidth',2,'color',colors(findmin(abs(score-linspace(0,1,1000)')),:));
else
plot(swDiffAll(i),ripPowerAll(i),'x','linewidth',2,'color',colors(findmin(abs(score-linspace(0,1,1000)')),:));
end
end
1
figure(1); % plot a clean figure without these points and check indivudual ripples
clf
plot(swDiffAll(idx1 | idx2),ripPowerAll(idx1 | idx2),'k.','markersize',1); hold on
scores = nan(size(ripPowerAll,1),1);
while true % click "Ctrl+C" to exit when you feel you're done
figure(1);
colors = Bright(1000);
[x,y,button] = ginput(1);
i = findmin(abs(x-swDiffAll)+abs(y-ripPowerAll));
t = featureTs/1250;
tl = (1:length(lfp))'/1250;
figure(2);
clf; interval = t(i) + [-1 1]*0.5; in = InIntervals(tl,interval);
subplot(3,1,1);
plot(tl(in) - t(i),lfp(in,:)); xlabel(num2str(t(i)));
PlotHVLines(Restrict(t,interval) - t(i),'v','k--');
legend('ripple channel','SW channel','noise channel');
subplot(3,1,2);
plot(tl(in) - t(i),double(swDiff(in)));
PlotHVLines(Restrict(t,interval) - t(i),'v','k--');
subplot(3,1,3);
plot(tl(in) - t(i),double(ripPower0(in)));
PlotHVLines(Restrict(t,interval) - t(i),'v','k--');
%     x = input( prompt )
score = str2double(input('good(1)/bad(0)?','s'));
scores(i) = score; % save scores to optionally save
figure(1);
if ~isnan(score) && score>0.1
hold on
plot(swDiffAll(i),ripPowerAll(i),'o','linewidth',2,'color',colors(findmin(abs(score-linspace(0,1,1000)')),:));
else
plot(swDiffAll(i),ripPowerAll(i),'x','linewidth',2,'color',colors(findmin(abs(score-linspace(0,1,1000)')),:));
end
end
1
0
1
0
1
0.7
0
1
0.1
1
0.6
1
0.5
1
0.8
1
0.9
1
0.8
1
0.3
0.6
0.8
0.2
0.5
1
while true % click "Ctrl+C" to exit when you feel you're done
figure(1);
colors = Bright(1000);
[x,y,button] = ginput(1);
i = findmin(abs(x-swDiffAll)+abs(y-ripPowerAll));
t = featureTs/1250;
tl = (1:length(lfp))'/1250;
figure(2);
clf; interval = t(i) + [-1 1]*0.5; in = InIntervals(tl,interval);
subplot(3,1,1);
plot(tl(in) - t(i),lfp(in,:)); xlabel(num2str(t(i)));
PlotHVLines(Restrict(t,interval) - t(i),'v','k--');
legend('ripple channel','SW channel','noise channel');
subplot(3,1,2);
plot(tl(in) - t(i),double(swDiff(in)));
PlotHVLines(Restrict(t,interval) - t(i),'v','k--');
subplot(3,1,3);
plot(tl(in) - t(i),double(ripPower0(in)));
PlotHVLines(Restrict(t,interval) - t(i),'v','k--');
%     x = input( prompt )
score = str2double(input('good(1)/bad(0)?','s'));
scores(i) = score; % save scores to optionally save
figure(1);
if ~isnan(score) && score>0.1
hold on
plot(swDiffAll(i),ripPowerAll(i),'o','linewidth',2,'color',colors(findmin(abs(score-linspace(0,1,1000)')),:));
elseif ~isnan(score)
plot(swDiffAll(i),ripPowerAll(i),'x','linewidth',2,'color',colors(findmin(abs(score-linspace(0,1,1000)')),:));
end
end
0.9
0.99
1
0.1
0
1
0
0.1
0.9
0
0.5
0.9
0.7
0.3
0.9
0.7
0.4
0.9
1
0.7
1
0.7
1
0.1
0.8
0
0.6
1
0.8
0.6
1
0.4
1
0.4
0.5
0.4
1
0.8
0
unselected = UIInPolygon(swDiffAll,ripPowerAll); % draw polygon to encircle the points you believe are ripples
unselected2 = UIInPolygon(swDiffAll,ripPowerAll); % draw polygon to encircle the points you believe are ripples
unselected = unselected | unselected2;
selected = ~unselected;
idx1 = selected & ~bad; % final ripples
idx2 = ~selected & ~bad; % non-ripples (yet free
size(bad)
size(selected)
size(t)
bad = false(size(t));
bad = bad | InIntervals(t,badIntervals);
size(bad)
idx1 = selected & ~bad; % final ripples
idx2 = ~selected & ~bad; % non-ripples (ye
clf
plot(swDiffAll(idx1 | idx2),ripPowerAll(idx1 | idx2),'k.','markersize',1); hold on
clf
plot(swDiffAll(idx2),ripPowerAll(idx2),'k.','markersize',1); hold on
plot(swDiffAll(idx1),ripPowerAll(idx1),'r.','markersize',1); legend('non-ripples','ripples');
save('detecting_ripples.mat','swDiffAll','ripPowerAll','idx1','idx2','bad','selected','scores');
ok = scores>0.1;
figure; DensityMap(swDiffAll(ok),ripPowerAll(ok));
ok = scores<=0.1;
DensityMap(swDiffAll(ok),ripPowerAll(ok));
map DensityMap(swDiffAll(ok),ripPowerAll(ok));
map =  DensityMap(swDiffAll(ok),ripPowerAll(ok));
open map
map0 =  DensityMap(swDiffAll(ok),ripPowerAll(ok));
map1 =  DensityMap(swDiffAll(ok)/max(swDiffAll),ripPowerAll(ok)/max(ripPowerAll));
ok = scores<=0.1; map0 =  DensityMap(swDiffAll(ok)/max(swDiffAll),ripPowerAll(ok)/max(ripPowerAll));
ok = scores>0.1; map1 =  DensityMap(swDiffAll(ok)/max(swDiffAll),ripPowerAll(ok)/max(ripPowerAll));
PlotColorMap(map1);
PlotColorMap(map0);
PlotColorMap(map1);
ok = scores>0.1; [map1,o1,o2] =  DensityMap(swDiffAll(ok)/max(swDiffAll),ripPowerAll(ok)/max(ripPowerAll));
open o1
open o2
ok = scores>0.1 & min([swDiffAll ripPowerAll],[],2)>0; [map1,o1,o2] =  DensityMap(swDiffAll(ok)/max(swDiffAll),ripPowerAll(ok)/max(ripPowerAll));
ok = scores<=0.1 & min([swDiffAll ripPowerAll],[],2)>0; [map0,o1,o2] =  DensityMap(swDiffAll(ok)/max(swDiffAll),ripPowerAll(ok)/max(ripPowerAll));
PlotColorMap(map1)
PlotColorMap(map0)
PlotColorMap(map1-map0)
clim
clim([-1 1]*0.1);
PlotColorMap(map1./mean(map1(:))-map0./mean(map0(:)))
hold all
plot(swDiffAll(ok)/max(swDiffAll)*50,ripPowerAll(ok)/max(ripPowerAll)*50,'ko')
ok = scores>0.1 & min([swDiffAll ripPowerAll],[],2)>0;
plot(swDiffAll(ok)/max(swDiffAll)*50,ripPowerAll(ok)/max(ripPowerAll)*50,'kx')
plot(swDiffAll(ok)/max(swDiffAll)*50,ripPowerAll(ok)/max(ripPowerAll)*50,'k*')
ok = scores<=0.1 & min([swDiffAll ripPowerAll],[],2)>0; [map0,o1,o2] =  DensityMap(swDiffAll(ok)/max(swDiffAll),ripPowerAll(ok)/max(ripPowerAll),'nBins',500);
ok = scores<=0.1 & min([swDiffAll ripPowerAll],[],2)>0; [map0,o1,o2] =  DensityMap(swDiffAll(ok)/max(swDiffAll),ripPowerAll(ok)/max(ripPowerAll),'nBins',500,'smooth',2);
ok = scores<=0.1 & min([swDiffAll ripPowerAll],[],2)>0; [map0,o1,o2] =  DensityMap(swDiffAll(ok)/max(swDiffAll),ripPowerAll(ok)/max(ripPowerAll),'nBins',500,'smooth',5);
ok = scores<=0.1 & min([swDiffAll ripPowerAll],[],2)>0; [map0,o1,o2] =  DensityMap(swDiffAll(ok)/max(swDiffAll),ripPowerAll(ok)/max(ripPowerAll),'nBins',500,'smooth',10);
ok = scores>0.1 & min([swDiffAll ripPowerAll],[],2)>0; [map1,o1,o2] =  DensityMap(swDiffAll(ok)/max(swDiffAll),ripPowerAll(ok)/max(ripPowerAll),'nBins',500,'smooth',10);
PlotColorMap(map1./mean(map1(:))-map0./mean(map0(:)))
PlotColorMap((map1./mean(map1(:))-map0./mean(map0(:))>0)+1)
hold all
plot(swDiffAll(ok)/max(swDiffAll)*50,ripPowerAll(ok)/max(ripPowerAll)*50,'k*')
clf
PlotColorMap((map1./mean(map1(:))-map0./mean(map0(:))>0)+1)
hold all
plot(swDiffAll(ok)/max(swDiffAll)*500,ripPowerAll(ok)/max(ripPowerAll)*500,'k*')
ok = scores<=0.1 & min([swDiffAll ripPowerAll],[],2)>0;
plot(swDiffAll(ok)/max(swDiffAll)*500,ripPowerAll(ok)/max(ripPowerAll)*500,'kx')
plot(swDiffAll(ok)/max(swDiffAll)*500,ripPowerAll(ok)/max(ripPowerAll)*500,'ko')
smooth
nBins
b = Bin([swDiffAll ripPowerAll],5);
dbcont
clf
close all
tic; lfpstruct = getLFP(69+1); display(['loaded! @' num2str(toc)]);
%     tic; lfpstruct = getLFP(104); HMC is 104, OR18 is 97
lfp = [lfpstruct.timestamps]; lfp(:,2) = lfpstruct.data;
r = [ripples.timestamps(:,1) ripples.peaks ripples.timestamps(:,2) ripples.RipMax ripples.SwMax];
PETH(lfp,r(:,1));
PETH(lfp,r(:,1));
FiguresOnMonitor2
lfp_pfc = lfp;
tic; lfpstruct = getLFP(60+1); display(['loaded! @' num2str(toc)]);
%     tic; lfpstruct = getLFP(104); HMC is 104, OR18 is 97
lfp = [lfpstruct.timestamps]; lfp(:,2) = lfpstruct.data;
lfp_ripple = lfp;
subplot(2,2,1);
PETH(lfp_pfc,r(:,1));
subplot(2,2,2);
PETH(lfp_ripple,r(:,1),'nBins',501);
subplot(2,2,1);
PETH(lfp_pfc,r(:,1),'nBins',501);
drawnow
tic; lfpstruct = getLFP(34+1); display(['loaded! @' num2str(toc)]);
lfp = [lfpstruct.timestamps]; lfp(:,2) = lfpstruct.data;
lfp_sw = lfp;
subplot(2,2,3);
PETH(lfp_sw,r(:,1),'nBins',501);
drawnow
load('day12.ripples.events.mat')
clear ripples
try load(fullfile(basepath,[basename '.MergePoints.events.mat']),'MergePoints');
postID = find(cellfun(@(x) ~isempty(strfind(x,'post')), MergePoints.foldernames));
sleep = ConsolidateIntervals(MergePoints.timestamps(cellfun(@(x) any(strfind(lower(x),'sleep')),MergePoints.foldernames),:));
catch
MergePoints = [];
display('No merge points saved. Keeping all ripples (impossible to restrict to post-task sleep');
end
MergePoints.timestamps(6,:)
diff(sleep,[],2)
diff(sleep,[],2)/3600
pre = InIntervals(r,sleep(1,:));
post = InIntervals(r,sleep(end,:));
open r
anovabar(diff(r(:,[1 3]),[],2)>0.1,[-pre+post]);
Accumulate((pre+post*2)+1,diff(r(:,[1 3]),[],2)>0.1)
Accumulate((pre+post*2)+1,1)
zBinomialComparison(26,2345,85,85)
zBinomialComparison(26,2345,85,6700)
Dist(1000,[diff(r(:,[1 3]),[],2),[-pre+post]+2],'grouped');
[h,ht] = Dist(1000,[diff(r(:,[1 3]),[],2),[-pre+post]+2],'grouped');
plot(ht,Smooth(h,[10 0]))
plot(ht,Smooth(cumsum(h),[0 0]))
clear all
load('day3.channelinfo.ripples.mat')
load('day10.ripples.events.mat')
SaveCustomEvents('ripples.rip.evt',r(:,1:3),{'ripple start','ripple peak','ripple stop'});
r = [ripples.timestamps(:,1) ripples.peaks ripples.timestamps(:,2) ripples.RipMax ripples.SwMax];
SaveCustomEvents('ripples.rip.evt',r(:,1:3),{'ripple start','ripple peak','ripple stop'});
9000/23000=
9000/23000
which GetLFP
raly
which getlfp
lfp = GetAyaLFP(122);
lfp_deep = lfp;
lfp_sup = GetAyaLFP(96);
lfp_deep = GetAyaLFP(122);
lfp_sup = GetAyaLFP(109);
d = [lfp(:,1) lfp_deep(:,2)-lfp_sup(:,2)];
d = [lfp_deep(:,1) lfp_deep(:,2)-lfp_sup(:,2)];
PlotXY(d);
PlotHVLines(4445.386,'v');
9/1.8926
hold all
SideAxes(gca,'top',0.5);
PlotXY(Restrict(lfp_deep,xlim));
hold all
PlotXY(Restrict(lfp_sup,xlim),'k');
clf
interval = 4442 + [0 5];
subplot(2,1,1);
PlotXY(Restrict(lfp_sup,interval),'k');
hold all
PlotXY(Restrict(lfp_deep,xlim));
subplot(2,1,2);
PlotXY(Restrict(d,interval));
10/
X = [4442.3961   4444.5836]   dX =
10/2.187
3/1.4
interval = 19106 + [0 5];
clf
subplot(2,1,1);
PlotXY(Restrict(lfp_sup,interval),'k');
hold all
PlotXY(Restrict(lfp_deep,xlim));
subplot(2,1,2);
PlotXY(Restrict(d,interval));
PlotHVLines(0,'h');
figure; PETH(d,r(:,1));
figure;PETH(lfp_deep,r(:,1));
peaks = FindLocalMaxima(d);
open peaks
PETH(peaks,r(:,1))
PETH(peaks,r(:,1),'durations',[-1 1]*2)
value = interp1(d(:,1),d(:,2),peaks);
clf
hist(value)
hist(value,100)
value0 = value; value = value(value>200);
Portion(value0>200)
PETH(peaks,r(:,1),'durations',[-1 1]*2)
PETH(peaks(value0>200),r(:,1),'durations',[-1 1]*2)
load('day10.MergePoints.events.mat')
sleep = ConsolidateIntervals(MergePoints.timestamps(cellfun(@(x) any(strfind(lower(x),'sleep')),MergePoints.foldernames),:));
pre = InIntervals(r,sleep(1,:));
post = InIntervals(r,sleep(end,:));
PETH(peaks(value0>200),r(pre,1),'durations',[-1 1]*2)
hold all
PETH(peaks(value0>200),r(post,1),'durations',[-1 1]*2)
peaks = FindLocalMaxima(lfp_d);
peaks = FindLocalMaxima(lfp_deep);
value = interp1(lfp_deep(:,1),lfp_deep(:,2),peaks);
figure; PETH(peaks(value0>500),r(pre,1),'durations',[-1 1]*2)
hold all
PETH(peaks(value0>500),r(post,1),'durations',[-1 1]*2)
figure; PETH(peaks(value>500),r(pre,1),'durations',[-1 1]*2)
hold all
PETH(peaks(value>500),r(post,1),'durations',[-1 1]*2)
peaks = FindLocalMaxima(lfp_d);
peaks = FindLocalMaxima(d);
deltas = peaks(values0>200);
deltas = peaks(value0>200);
clf
hist(deltas,1000)
hist(r(:,1),1000)
deltaWaves.peaks = deltas; deltaWaves.peakNormedPower = zscore(value0(value0>200));
deltaWaves.detectorName = ['peaks in difference signal b/n deep channel 122 and sup. channel 109, thresholded at >200 signal units'];
save(fullfile(basepath,[basename '.deltaWaves.events.mat']),'deltaWaves');
basepath = pwd
save(fullfile(basepath,[basename '.deltaWaves.events.mat']),'deltaWaves');
basename = pwd
basename = basenameFromBasepath(basepath);
save(fullfile(basepath,[basename '.deltaWaves.events.mat']),'deltaWaves');
SaveCustomEvents('deltas.del.evt',deltas,{'putative delta peak'});
interval
in = InIntervals(d,interval);
figure; plot(d(in,1), Smooth(d(in,2),10));
plot(d(in,1), Smooth(d(in,2),20));
20/1250
plot(d(in,1), Smooth(d(in,2),1250));
plot(d(in,1), Smooth(d(in,2),1250*0.1));
plot(d(in,1), Smooth(d(in,2),1250*0.01));
plot(d(in,1), Smooth(d(in,2),1250*0.02));
peaks = FindLocalMaxima([d(:,1) Smooth(d(:,2),0.02*1250)]);
value = interp1(d(:,1),d(:,2),peaks);
hist(value)
hist(value,100)
plot(d(in,1), Smooth(d(in,2),1250*0.02));
PETH(peaks(value>200),r(pre,1),'durations',[-1 1]*2)
hold all
PETH(peaks(value>200),r(post,1),'durations',[-1 1]*2)
deltas = peaks(value>200);
figure; plot(d(in,1), Smooth(lfp_deep(in,2),10));
peaks = FindLocalMaxima([d(:,1) Smooth(lfp_deep(:,2),0.02*1250)]);
value = interp1(lfp_deep(:,1),lfp_deep(:,2),peaks);
PETH(peaks(value>500),r(post,1),'durations',[-1 1]*2)
hold all
PETH(peaks(value>200),r(pre,1),'durations',[-1 1]*2)
clf
PETH(peaks(value>500),r(pre,1),'durations',[-1 1]*2)
hold all
PETH(peaks(value>500),r(post,1),'durations',[-1 1]*2)
PETH(deltas,r(pre,1),'durations',[-1 1]*2)
PETH(deltas,r(post,1),'durations',[-1 1]*2)
clf
PETH(deltas,peaks(value>500))
PETH(deltas,peaks(value>500),'durations',[-1 1])
PETH(deltas,peaks(value>500),'durations',[-1 1]*0.2)
deltaWaves.detectorName = ['peaks in difference signal (smooth 20ms) b/n deep channel 122 and sup. channel 109, thresholded at >200 signal units'];
deltaWaves.peaks = deltas; deltaWaves.peakNormedPower = zscore(value(value>200));
save(fullfile(basepath,[basename '.deltaWaves.events.mat']),'deltaWaves');
SaveCustomEvents('deltas.del.evt',deltas,{'putative delta peak'});
clf
hist(peaks(value>500),100);
figure
hist(deltas,100);
int = interp1(d(:,1),lfp_deep(:,2),deltas);
figure; hist(int,100);
hist(deltas(int>0),100);
hist(deltas(int<0),100);
PETH(deltas,r(:,1),'durations',[-1 1]*2)
PETH(deltas(int<0),r(:,1),'durations',[-1 1]*2)
PETH(deltas(int>0),r(:,1),'durations',[-1 1]*2)
deltas = deltas(int>0);
Portion(int>0)
deltaWaves.peaks = deltas; deltaWaves.peakNormedPower = zscore(value(value>200));
deltaWaves.detectorName = ['peaks in difference signal (smooth 20ms) b/n deep channel 122 and sup. channel 109, thresholded at >200 signal units as long as deep signal is positive'];
save(fullfile(basepath,[basename '.deltaWaves.events.mat']),'deltaWaves');
SaveCustomEvents('deltas.del.evt',deltas,{'putative delta peak'});
PETH(deltas,r(:,1),'durations',[-1 1]*2)
PETH(deltas,r(pre,1),'durations',[-1 1]*2)
hold all
PETH(deltas,r(post,1),'durations',[-1 1]*2)
clear all
edit preprocessSession.m
load('day3.channelinfo.ripples.mat')
rippleChannels.Ripple_Channel = 31+1; rippleChannels.Sharpwave_Channel = 9+1; rippleChannels.Noise_Channel = 71+1;
ripples = DetectSWR([rippleChannels.Ripple_Channel rippleChannels.Sharpwave_Channel rippleChannels.Noise_Channel],'basepath',basepath,'saveMat',true,'forceDetect',true,'check',true);
basepath = pwd;
ripples = DetectSWR([rippleChannels.Ripple_Channel rippleChannels.Sharpwave_Channel rippleChannels.Noise_Channel],'basepath',basepath,'saveMat',true,'forceDetect',true,'check',true);
tl = (1:length(lfp))'/1250; % timestamps for lfp
t = featureTs/1250; % timestamps for candidate ripple events
bad = false(size(t));
[clean,bad,badIntervals] = CleanLFP([tl,double(lfp(:,2))],'thresholds',[10 Inf]);
smoothed = Smooth(abs(swDiff),1250*10);
noisyPeriods = ConsolidateIntervals(tl(FindInterval(smoothed>500)),'epsilon',1);
sw = interp1(tl,lfpLow(:,2),t);
Portion(InIntervals(t,badIntervals))
Portion(InIntervals(t,noisyPeriods))
Portion(sw>0)
max(ripPowerAll(sw>0))
find(ripPowerAll>3774 & sw>0)
addComma(t(68959))
quantile(ripPowerAll(sw>0),0.8)
sum(ripPowerAll>93 & sw>0)
FindClosest(ripPowerAll.*(sw>0),100)
addComma(t(18454))
clf
plot(swDiffAll(sw>0),ripPowerAll(sw>0),'k.','markersize',1); hold on
hold all
plot(swDiffAll(sw<0),ripPowerAll(sw<0),'b.','markersize',1); hold on
plot(swDiffAll(sw>0),ripPowerAll(sw>0),'k.','markersize',1); hold on
plot(swDiffAll(sw>0),ripPowerAll(sw>0),'ko','markersize',1); hold on
plot(swDiffAll(sw>0),ripPowerAll(sw>0),'ko','markersize',10); hold on
FindClosest(ripPowerAll.*(sw>0),525.46)
addComma(t(103114))
sw(103114)
bad = false(size(t));
bad = bad | InIntervals(t,badIntervals);
bad = bad | InIntervals(t,noisyPeriods);
idx1 = idx1 & ~bad; idx2 = idx2 & ~bad; % remove bad points from "idx1" and "idx2"
figure(1); % plot a clean figure without these points and check indivudual ripples
clf
plot(swDiffAll(idx1 | idx2),ripPowerAll(idx1 | idx2),'k.','markersize',1); hold on
scores = nan(size(ripPowerAll,1),1);
while true % click "Ctrl+C" to exit when you feel you're done
figure(1);
colors = Bright(1000);
[x,y,button] = ginput(1);
i = findmin(abs(x-swDiffAll)+abs(y-ripPowerAll));
t = featureTs/1250;
tl = (1:length(lfp))'/1250;
figure(2);
clf; interval = t(i) + [-1 1]*0.5; in = InIntervals(tl,interval);
subplot(3,1,1);
plot(tl(in) - t(i),lfp(in,:)); xlabel(num2str(t(i)));
PlotHVLines(Restrict(t,interval) - t(i),'v','k--');
legend('ripple channel','SW channel','noise channel');
subplot(3,1,2);
plot(tl(in) - t(i),double(swDiff(in)));
PlotHVLines(Restrict(t,interval) - t(i),'v','k--');
subplot(3,1,3);
plot(tl(in) - t(i),double(ripPower0(in)));
PlotHVLines(Restrict(t,interval) - t(i),'v','k--');
%     x = input( prompt )
score = str2double(input('good(1)/bad(0)?','s'));
scores(i) = score; % save scores to optionally save
figure(1);
if ~isnan(score) && score>0.1
hold on
plot(swDiffAll(i),ripPowerAll(i),'o','linewidth',2,'color',colors(findmin(abs(score-linspace(0,1,1000)')),:));
elseif ~isnan(score)
plot(swDiffAll(i),ripPowerAll(i),'x','linewidth',2,'color',colors(findmin(abs(score-linspace(0,1,1000)')),:));
end
end
1
0
1
score(i) = nan;
plot(swDiffAll(i),ripPowerAll(i),'o','linewidth',2,'color','w');
score = [];
scores(i) = nan;
while true % click "Ctrl+C" to exit when you feel you're done
figure(1);
colors = Bright(1000);
[x,y,button] = ginput(1);
i = findmin(abs(x-swDiffAll)+abs(y-ripPowerAll));
t = featureTs/1250;
tl = (1:length(lfp))'/1250;
figure(2);
clf; interval = t(i) + [-1 1]*0.5; in = InIntervals(tl,interval);
subplot(3,1,1);
plot(tl(in) - t(i),lfp(in,:)); xlabel(num2str(t(i)));
PlotHVLines(Restrict(t,interval) - t(i),'v','k--');
legend('ripple channel','SW channel','noise channel');
subplot(3,1,2);
plot(tl(in) - t(i),double(swDiff(in)));
PlotHVLines(Restrict(t,interval) - t(i),'v','k--');
subplot(3,1,3);
plot(tl(in) - t(i),double(ripPower0(in)));
PlotHVLines(Restrict(t,interval) - t(i),'v','k--');
%     x = input( prompt )
score = str2double(input('good(1)/bad(0)?','s'));
scores(i) = score; % save scores to optionally save
figure(1);
if ~isnan(score) && score>0.1
hold on
plot(swDiffAll(i),ripPowerAll(i),'o','linewidth',2,'color',colors(findmin(abs(score-linspace(0,1,1000)')),:));
elseif ~isnan(score)
plot(swDiffAll(i),ripPowerAll(i),'x','linewidth',2,'color',colors(findmin(abs(score-linspace(0,1,1000)')),:));
end
end
1
0
1
ripples = DetectSWR([rippleChannels.Ripple_Channel rippleChannels.Sharpwave_Channel rippleChannels.Noise_Channel],'basepath',basepath,'saveMat',true,'forceDetect',true,'check',true);
commandwindow
ripples = DetectSWR([rippleChannels.Ripple_Channel rippleChannels.Sharpwave_Channel rippleChannels.Noise_Channel],'basepath',basepath,'saveMat',true,'forceDetect',true,'check',true);
% First, remove points that occur in noisy lfp periods
tl = (1:length(lfp))'/1250; % timestamps for lfp
t = featureTs/1250; % timestamps for candidate ripple events
bad = false(size(t));
[clean,bad,badIntervals] = CleanLFP([tl,double(lfp(:,2))],'thresholds',[10 Inf]);
bad = bad | InIntervals(t,badIntervals);
smoothed = Smooth(abs(swDiff),1250*10);
noisyPeriods = ConsolidateIntervals(tl(FindInterval(smoothed>500)),'epsilon',1);
bad = bad | InIntervals(t,noisyPeriods);
0
1
figure; plot(tl,lfp);
PlotIntervals(badIntervals,'color','k');
[clean,bad,badIntervals] = CleanLFP([tl,double(lfp(:,2))],'thresholds',[5 Inf]);
PlotIntervals(badIntervals,'color','k');
dur(badIntervals)
sum(bad)
Portion(bad)
bad = bad | InIntervals(t,badIntervals);
size(bad)
size(t)
bad = false(size(t));
bad = bad | InIntervals(t,badIntervals);
bad = bad | InIntervals(t,noisyPeriods);
Portion(bad)
1-Portion(idx1|idx2)
Portion(idx1|idx2)
idx1 = idx1 & ~bad; idx2 = idx2 & ~bad; % remove bad points from "idx1" and "idx2"
Portion(idx1|idx2)
bad(i)
interval = t(i)+[-1 1]*10;
in = InIntervals(tl,interval);
figure; plot(tl(in),lfp(in,:));
PlotIntervals(Restrict(badIntervals,interval))
[clean,~,badIntervals] = CleanLFP([tl,double(lfp(:,2))],'thresholds',[5 10]);
PlotIntervals(Restrict(badIntervals,interval))
[clean,~,badIntervals] = CleanLFP([tl,double(lfp(:,2))],'thresholds',[5 10]);
interval = 1142.6+[-1 1]*10;
in = InIntervlas(t,interval);
in = InIntervals(t,interval);
figure; plot(t(in),abs(z(in)));
figure; plot(t(in),abs(d(in)),'r');
plot(t,d);
PlotHVLines(interval);
dbcont
[clean2,~,badIntervals2] = CleanLFP([tl,double(lfp(:,3))],'thresholds',[5 10]);
PlotIntervals(Restrict(badIntervals2,interval))
dbcont
PlotIntervals(Restrict(badIntervals2,interval))
dur(badIntervals)
dur(badIntervals2)
dur(badIntervals2)/3600
[clean2,~,badIntervals2] = CleanLFP([tl,double(lfp(:,3))],'thresholds',[10 10]);
dur(badIntervals2)/3600
PlotIntervals(Restrict(badIntervals2,interval))
PlotIntervals(Restrict(badIntervals2,interval),'color','r')
[clean2,~,badIntervals2] = CleanLFP([tl,double(lfp(:,3))],'thresholds',[7 10]);
PlotIntervals(Restrict(badIntervals2,interval),'color','r')
[clean2,~,badIntervals2] = CleanLFP([tl,double(lfp(:,3))],'thresholds',[6 10]);
PlotIntervals(Restrict(badIntervals2,interval),'color','r')
Restrict(badIntervals2,interval)
interval
[clean2,~,badIntervals2] = CleanLFP([tl,double(lfp(:,3))],'thresholds',[5 10]);
PlotIntervals(Restrict(badIntervals2,interval),'color','r')
[clean2,~,badIntervals2] = CleanLFP([tl,double(lfp(:,3))],'thresholds',[5 10],'aroundartefact',[0.5 0.1]);
PlotIntervals(Restrict(badIntervals2,interval),'color','k')
clf
plot(tl(in),lfp(in,:));
PlotIntervals(Restrict(badIntervals2,interval),'color','k')
[clean2,~,badIntervals2] = CleanLFP([tl,double(lfp(:,3))],'thresholds',[5 10],'aroundartefact',[1 0.1]);
PlotIntervals(Restrict(badIntervals2,interval),'color','k')
[clean2,~,badIntervals2] = CleanLFP([tl,double(lfp(:,3))],'thresholds',[4 10],'aroundartefact',[1 0.1]);
PlotIntervals(Restrict(badIntervals2,interval),'color','k')
[clean2,~,badIntervals2] = CleanLFP([tl,double(lfp(:,3))],'thresholds',[4 10],'aroundartefact',[0.5 0.1]);
dur(badIntervals2)/3600
clf
Portion(InIntervals(t,badIntervals2))
FindClosest(ripPowerAll.*(InIntervals(t,badIntervals2)),600)
addComma(t(350461))
[clean2,~,badIntervals2] = CleanLFP([tl,double(lfp(:,3))],'thresholds',[4 10],'aroundartefact',[0.1 0.1]);
dur(badIntervals2)/3600
FindClosest(ripPowerAll.*(InIntervals(t,badIntervals2)),600)
addComma(t(197093))
ripPowerAll(197093)
FindClosest(ripPowerAll.*(InIntervals(t,badIntervals2)),500)
addComma(t(312832))
FindClosest(ripPowerAll.*(InIntervals(t,badIntervals2)),498)
FindClosest(ripPowerAll.*(InIntervals(t,badIntervals2)),518)
addComma(t(232288))
FindClosest(ripPowerAll.*(InIntervals(t,badIntervals2)),450)
addComma(t(237685))
FindClosest(ripPowerAll.*(InIntervals(t,badIntervals2)),480)
addComma(t(197089))
size(bad)
size(t)
Portion(bad)
bad = bad | InIntervals(t,badIntervals2);
Portion(bad)
idx1 = idx1 & ~bad; idx2 = idx2 & ~bad; % remove bad points from "idx1" and "idx2"
bad(i)
t(i)
FindClosest(badIntervals2(:,1),ans)
badIntervals2(328,:)
PlotInterlals(ans);
PlotIntervals(ans-t(i));
0
1
0
1
0
1
0
1
0
size(z)
size(matrix)
matrix = [swDiffAll ripPowerAll];
z = bsxfun(@rdivide,bsxfun(@minus,matrix,mean(matrix(~bad,:))),std(matrix(~bad,:)));
open z
size(distances)
distances = bsxfun(@minus,z,z(~isnan(scores),:));
distances = bsxfun(@minus,z,z(~isnan(scores),:)');
distances = bsxfun(@minus,z(:,1),z(~isnan(scores),1));
distances = bsxfun(@minus,z(:,1),z(~isnan(scores),1)');
distances = sqrt(bsxfun(@minus,z(:,1),z(~isnan(scores),1)').^2 + bsxfun(@minus,z(:,2),z(~isnan(scores),2)').^2);
figure; PlotColorMap(distances)
[~,nearestNeighbour] = min(distances,[],2);
estimated = scores(nearestNeighbour);
figure; distances = sqrt(bsxfun(@minus,z(:,1),z(~isnan(scores),1)').^2 + bsxfun(@minus,z(:,2),z(~isnan(scores),2)').^2);
[~,nearestNeighbour] = min(distances,[],2);
estimated = scores(nearestNeighbour);
scatter(matrix(:,1),matrix(:,2),5,estimated,'filled'); clim([0 1]);
estimated = scores(scored(nearestNeighbour));
scored = find(~isnan(scores));
estimated = scores(scored(nearestNeighbour));
scatter(matrix(:,1),matrix(:,2),5,estimated,'filled'); clim([0 1]);
scatter(matrix(:,1),matrix(:,2),5,estimated,'filled'); colormap(Bright); clim([0 1]);
scatter(matrix(:,1),matrix(:,2),5,estimated,'filled'); colormap(jet); clim([0 1]);
scatter(matrix(:,1),matrix(:,2),1,estimated,'filled'); colormap(Bright); clim([0 1]);
size(x)
size(y)
x
y
xlims = xlim; ylims = ylims; clf
% extrapolate the score of each ripple based on the score of the nearest scored ripple
scored = find(~isnan(scores));
distances = sqrt(bsxfun(@minus,z(:,1),z(scored,1)').^2 + bsxfun(@minus,z(:,2),z(~isnan(scores),2)').^2);
[~,nearestNeighbour] = min(distances,[],2);
estimated = scores(scored(nearestNeighbour));
scatter(matrix(:,1),matrix(:,2),1,estimated,'filled'); colormap(Bright); clim([0 1]);
xlim(xlims); ylim(ylims);
xlims = xlim; ylims = ylim; clf
% extrapolate the score of each ripple based on the score of the nearest scored ripple
scored = find(~isnan(scores));
distances = sqrt(bsxfun(@minus,z(:,1),z(scored,1)').^2 + bsxfun(@minus,z(:,2),z(~isnan(scores),2)').^2);
[~,nearestNeighbour] = min(distances,[],2);
estimated = scores(scored(nearestNeighbour));
scatter(matrix(:,1),matrix(:,2),1,estimated,'filled'); colormap(Bright); clim([0 1]);
xlim(xlims); ylim(ylims);
hold on; ok = scores>0.1;
scatter(swDiffAll(ok),ripPowerAll(i),10,scores(ok));
scatter(swDiffAll(ok),ripPowerAll(ok),10,scores(ok));
hold on; ok = scores>0.1;
scatter(swDiffAll(ok),ripPowerAll(ok),20,scores(ok));
hold on; ok = scored
scatter(swDiffAll(ok),ripPowerAll(ok),20,scores(ok));
scatter(swDiffAll(ok),ripPowerAll(ok),20,scores(ok)); scatter(swDiffAll(ok),ripPowerAll(ok),15,scores(ok));
0.8
figure(1);
0.5
xlim([-2000 4000]); ylim([-200 800]);
1
0
0.1
0.9
0.1
0.7
0.6
1
edit nothing
profile on
nothing
0.555
0
profile viewer
nothing
1
12
scores(i)
scores(i) = 1;
i = find(scores==0.555)
scores(i) = nan;
xlim([-2000 4000]); ylim([-200 500]);
xlim([-1000 4000]); ylim([-200 500]);
xlim([-1000 3000]); ylim([-200 500]);
nothing
0
1
figure(1);
xlims = xlim; ylims = ylim; clf
% extrapolate the score of each ripple based on the score of the nearest scored ripple
scored = find(~isnan(scores));
distances = sqrt(bsxfun(@minus,z(:,1),z(scored,1)').^2 + bsxfun(@minus,z(:,2),z(~isnan(scores),2)').^2);
[~,nearestNeighbour] = min(distances,[],2);
estimated = scores(scored(nearestNeighbour));
scatter(matrix(:,1),matrix(:,2),1,estimated,'filled'); colormap(Bright); clim([0 1]);
xlim(xlims); ylim(ylims);
hold on;
scatter(swDiffAll(scored),ripPowerAll(scored),20,scores(scored)); scatter(swDiffAll(scored),ripPowerAll(scored),15,scores(scored));
nothing
0
1
0
0.1
1
0
1
0.8
0
0.8
0
1
0.8
0
0.8
1
0
0.3
0.4
0.2
0.4
0
1
0
1
0.1
0.5
0.1
0.5
0
0.8
0.5
0
1
0.2
0
1
saved(i)(
scores(i)
scores(i) = 0;
nothing
0
0.4
0
1
0
1
0.5
0.1
0
0.05
0.5
0
0.5
0
1
0
1
0
1
save(fullfile(basepath,'detecting_ripples_temp.mat'),'swDiffAll','ripPowerAll','bad','scores');
scored = find(~isnan(scores));
distances = sqrt(bsxfun(@minus,z(:,1),z(scored,1)').^2 + bsxfun(@minus,z(:,2),z(~isnan(scores),2)').^2);
[~,nearestNeighbour] = min(distances,[],2);
estimated = scores(scored(nearestNeighbour));
selected = estimated>0.1;
idx1 = selected & ~bad; % final ripples
idx2 = ~selected & ~bad; % non-ripples (yet free from noise as well)
display('saving...');
save(fullfile(basepath,'DetectSWR_manual_scoring.mat'),'swDiffAll','ripPowerAll','idx1','idx2','bad','selected','scores');
figure;  clf
plot(swDiffAll(idx2),ripPowerAll(idx2),'k.','markersize',1); hold on
plot(swDiffAll(idx1),ripPowerAll(idx1),'r.','markersize',1); legend('non-ripples','ripples');
disp(['Please review the figure and change the idx1 (ripples) ' newline 'and idx2 (calm non-ripple periods) groups.' newline ...
'Note that noisy periods need not participate in either group.']);
edit SaveFig
gcf
saveas(figure(1),fullfile(basepath,'DetectSWR_manual_scoring.fig');
saveas(figure(1),fullfile(basepath,'DetectSWR_manual_scoring.fig'));
clabel('ripple score');
xlabel('Sharp wave depth'); ylabel('Ripple power');
set(gca,'fontsize',15);
clabel('Estimated ripple score');
title(['Manual scoring for ' basepath]);
'ripple channel: ' num2str(Channels(1))
title({['Manual scoring for ' basepath],['inputted channels: ' num2str(Channels)]);
title({['Manual scoring for ' basepath],['inputted channels: ' num2str(Channels)]});
title({['Manual scoring for ' basepath],['called with channels: ' num2str(Channels) ' (ripple, sharp wave, and (optionally) noise']});
title({['Manual scoring for ' basepath],['called with channels: ' num2str(Channels) ' (ripple, sharp wave, and (optionally) noise)']});
set(gca,'fontsize',15,'TickDir','out');
saveas(figure(1),fullfile(basepath,'DetectSWR_manual_scoring.fig'));
45
figure(1); % plot a clean figure without these points and check indivudual ripples
clf
plot(swDiffAll(idx1 | idx2),ripPowerAll(idx1 | idx2),'k.','markersize',1); hold on
scores = nan(size(ripPowerAll,1),1);
clf
plot(swDiffAll(idx2),ripPowerAll(idx2),'k.','markersize',1); hold on
plot(swDiffAll(idx1),ripPowerAll(idx1),'r.','markersize',1); legend('non-ripples','ripples');
dbcont
r = [ripples.timestamps(:,1) ripples.peaks ripples.timestamps(:,2) ripples.RipMax ripples.SwMax];
hist(r(:,1),1000)
6/0.057
11/3.2
load('day3.MergePoints.events.mat')
sleep = ConsolidateIntervals(MergePoints.timestamps(cellfun(@(x) any(strfind(lower(x),'sleep')),MergePoints.foldernames),:));
pre = InIntervals(r,sleep(1,:));
post = InIntervals(r,sleep(end,:));
anovabar(diff(r(:,[1 3]),[],2)>0.1,[-pre+post]);
anovabar(diff(r(:,[1 3]),[],2),[-pre+post]);
[h,ht] = Dist(1000,[diff(r(:,[1 3]),[],2),[-pre+post]+2],'grouped');
plot(ht,Smooth(cumsum(h),[0 0]))
Portion(diff(:,[1 3]),[],2)>0.1)
Portion(diff(r(:,[1 3]),[],2)>0.1)
plot(ht,Smooth(cumsum(h),[0 0]))
figure; anovabar(diff(r(:,[1 3]),[],2),[-pre+post]);
anovabar(diff(r(:,[1 3]),[],2)>0.8,[-pre+post]);
anovabar(diff(r(:,[1 3]),[],2)>0.08,[-pre+post]);
plot(ht,Smooth((h),[0 0]))
plot(log(ht),Smooth((h),[0 0]))
plot(log(ht),Smooth((h),[1 0]))
plot(log(ht),Smooth((h),[2 0]))
plot(log(ht),Smooth((h),[5 0]))
plot((ht),Smooth((h),[5 0])); set(gca,'XScale','log');
plot((ht),Smooth((h),[5 0])); set(gca,'XScale','log'); xlim([0.01 1])
plot((ht),Smooth((h),[5 0])); set(gca,'XScale','log'); xlim([0.001 1])
plot((ht),Smooth((h),[10 0])); set(gca,'XScale','log'); xlim([0.001 1])
PlotHVLines(0.08,'v');
plot((ht),Smooth((h),[0 0])); set(gca,'XScale','log'); xlim([0.001 1])
PlotHVLines(0.08,'v');
plot((ht),Smooth((h),[1 0])); set(gca,'XScale','log'); xlim([0.001 1])
clf
plot((ht),Smooth((h),[1 0])); set(gca,'XScale','log'); xlim([0.001 1])
PlotHVLines(0.08,'v');
[h,ht] = Dist(100,[diff(r(:,[1 3]),[],2),[-pre+post]+2],'grouped');
clf
plot((ht),Smooth((h),[1 0])); set(gca,'XScale','log'); xlim([0.001 1])
PlotHVLines(0.08,'v');
clf
plot((ht),Smooth((h),[0 0])); set(gca,'XScale','log'); xlim([0.001 1])
PlotHVLines(0.08,'v');
anovabar(diff(r(:,[1 3]),[],2)>0.07,[-pre+post]);
lfp = GetAyaLFP(31);
[s, wf, median, variability] = EventCenteredWavelet(lfp,[r(:,1)-0.5 r(:,1)+0.5]);
PlotColorMap(median);
open wf
hist(wf)
wf1 = wf(1):0.1:2f(end);
wf1 = wf(1):0.1:wf(end);
inter = interp1(wf(:),s,wf1(:));
inter = interp1(wf(:),median,wf1(:));
PlotColorMap(inter,'y',wf1);
clim
clim(clim/2)
wt = linspace(-1,1,size(s,2))*0.5;
PlotColorMap(inter,'y',wf1,'x',wt);
PlotColorMap(inter,'y',wf1,'x',wt); clim([0 20000]);
PlotColorMap(inter,'y',wf1,'x',wt); clim([0 200000]);
PlotColorMap(inter,'y',wf1,'x',wt); clim([0 2000000]);
median1 = median(s(:,:,pre),3);
median1 = nanmedian(s(:,:,pre),3);
median2 = nanmedian(s(:,:,~pre & ~post),3);
median3 = nanmedian(s(:,:,post),3);
subplot(2,2,1);
inter = interp1(wf(:),median1,wf1(:));
PlotColorMap(inter,'y',wf1,'x',wt); clim([0 2000000]);
subplot(2,2,2);
inter = interp1(wf(:),median2,wf1(:));
PlotColorMap(inter,'y',wf1,'x',wt); clim([0 2000000]);
subplot(2,2,3);
inter = interp1(wf(:),median3,wf1(:));
PlotColorMap(inter,'y',wf1,'x',wt); clim([0 2000000]);
clim
clim([0 1000000])
54
q = s(:,:,1);
open q
s(s==0) = nan;
45
clims(clim*0.8)
clims(clim*1.2)
raly
5
ripples = DetectSWR([rippleChannels.Ripple_Channel rippleChannels.Sharpwave_Channel rippleChannels.Noise_Channel],'basepath',basepath,'saveMat',true,'forceDetect',true,'check',true);
load('DetectSWR_manual_scoring.mat')
plot(swDiffAll(idx2),ripPowerAll(idx2),'k.','markersize',1); hold on
plot(swDiffAll(idx1),ripPowerAll(idx1),'r.','markersize',1); legend('non-ripples','ripples');
disp(['Please review the figure and change the idx1 (ripples) ' newline 'and idx2 (calm non-ripple periods) groups.' newline ...
'Note that noisy periods need not participate in either group.']);
dbcont
t = SWRs/1250; % timestamps for candidate ripple events
FindClosest(t,8371)
Restrict(t,[0 1]+8370.8)
Restrict(t,[0 1]+8370.8)-8371
interval = [0 1]+8370.8;
in = InIntervals(t,interval);
tl = (1:length(lfp))'/1250; % timestamps for lfp
in = InIntervals(tl,interval);
PlotXY(tl(in),lfp(in,1:2));
clf
plot(tl(in),lfp(in,1:2));
PlotHVLines(Restrict(t,[0 1]+8370.8));
PlotHVLines(Restrict(t,[0 1]+8370.8),'v--','k','linewidth',2);
PlotHVLines(Restrict(t,[0 1]+8370.8),'v','k--','linewidth',2);
clf
plot(tl(in),lfp(in,1:2));
PlotHVLines(Restrict(t,[0 1]+8370.8),'v','k--','linewidth',2);
interval = [0 1]+8370.8;
in = InIntervals(tl,interval);
clf
plot(tl(in),lfp(in,1:2));
PlotHVLines(Restrict(t,[0 0.5]+8370.8),'v','k--','linewidth',2);
interval = [0 0.5]+8370.8;
in = InIntervals(tl,interval);
clf
plot(tl(in),lfp(in,1:2));
PlotHVLines(Restrict(t,[0 0.5]+8370.8),'v','k--','linewidth',2);
interval = [0 0.5]+8371;
in = InIntervals(tl,interval);
clf
plot(tl(in),lfp(in,1:2));
PlotHVLines(Restrict(t,[0 0.5]+8370.8),'v','k--','linewidth',2);
interval = [-0.1 0.4]+8371;
in = InIntervals(tl,interval);
clf
plot(tl(in),lfp(in,1:2));
PlotHVLines(Restrict(t,[0 0.5]+8370.8),'v','k--','linewidth',2);
ylim
ylim(ylim*1.2);
\
ylim(ylim*1.2);
interval = [-0.1 0.4]+8371;
in = InIntervals(tl,interval);
clf
plot(tl(in),lfp(in,1:2)); ylim([-1 1]*5000);
PlotHVLines(Restrict(t,[0 0.5]+8370.8),'v','k--','linewidth',2);
FindClosest(t,8371)
find(InIntervals(t,interval)
find(InIntervals(t,interval))
for ii = 1:5636-1
if DEBUG
set(gcf,'Name',sprintf('Detection # %d out of %d',ii,Nswr))
end
% threshold feature distributions to localize about detection
% 1) detected sharp wave difference magnitude must be > thresSD locally
swDiff_chk   = swDiff(SWRs(ii)-bound:SWRs(ii)+bound);
MswDiff  = median(swDiff_chk);
SDswDiff = std(swDiff_chk);
if swDiff_chk(bound+1) < MswDiff + thresSDswD(2)*SDswDiff && useSPW
if DEBUG
disp('Sharp Wave Difference magnitude too small locally')
end
SWRindF(SWRid(ii)) = false;
continue
end
% 2) detected ripple power must be > thresSD locally
ripPow_chk    = ripPower0(SWRs(ii)-bound:SWRs(ii)+bound);
MripPow   = median(ripPow_chk);
SDripPow  = std(ripPow_chk);
if ripPowerAll(SWRid(ii)) < MripPow + thresSDrip(2)*SDripPow
if DEBUG
disp('Ripple Power too low locally')
end
SWRindF(SWRid(ii)) = false;
continue
end
% 3) detections can't be within a minimum inter-event interval thresh
if ii < Nswr
if SWR_diff(ii) < minIsi
if DEBUG
disp('IEI is too small')
end
SWRindF(SWRid(ii)) = false;
continue
end
end
% Estimate duration of detected Sharp-Wave event
startSWR = find(swDiff_chk(1:bound + 1) < MswDiff + thresSDswD(1)*SDswDiff,1,'last');
stopSWR  = find(swDiff_chk(bound + 1: end) < MswDiff + thresSDswD(1)*SDswDiff,1,'first');
durSWR   = (bound + 1 - startSWR) + stopSWR;
% Estimate duration of detected ripple power increase event
% find index of maximum ripple power increase around Sharp-Wave
[~,maxRpInd]  = max(ripPow_chk(bound+1-HalfWinSize:bound+1+HalfWinSize));
startRP       = find(ripPow_chk(1:bound - HalfWinSize + maxRpInd) ...
< MripPow + thresSDrip(1)*SDripPow,1,'last');
stopRP        = find(ripPow_chk(bound - HalfWinSize + maxRpInd: end) ...
< MripPow + thresSDrip(1)*SDripPow,1,'first');
durRP         = (bound - HalfWinSize + maxRpInd - startRP) + stopRP;
% 4) events associated with detections must exceed a minimum duration
if durSWR < minDurSWs && useSPW
if DEBUG
disp('Sharp-Wave Event Duration is too short.')
end
SWRindF(SWRid(ii)) = false;
continue
end
if durRP < minDurRPs
if DEBUG
disp('Ripple Event Duration is too short.')
end
SWRindF(SWRid(ii)) = false;
continue
end
% 5) sharp-wave must not exceed a maximum duration
if durSWR > maxDurSWs && useSPW
if DEBUG
disp('Sharp-Wave Event Duration is too long.')
end
SWRindF(SWRid(ii)) = false;
continue
end
% We have a valid detection
valid                     = valid + 1;
SWR_valid.Ts(valid,:)     = [SWRs(ii) (SWRs(ii)-(bound+1-startSWR)) (SWRs(ii)+stopSWR-1)];
swMax                     = swDiffAll(SWRid(ii));
ripMax                    = ripPowerAll(SWRid(ii));
SWR_valid.SwMax(valid,:)  = ...
[(swMax - MswDiff)/SDswDiff sum(swDiff_chk < swMax)/(2*bound + 1)];
SWR_valid.RipMax(valid,:) = ...
[(ripMax - MripPow)/SDripPow sum(ripPow_chk < ripMax)/(2*bound + 1)];
if DEBUG
h1 = subplot(2,3,[1 2]);
plot(timeline,swDiff_chk), axis tight, hold on
plot([timeline(startSWR) timeline(bound+stopSWR)],...
[swDiff_chk(startSWR) swDiff_chk(bound+stopSWR)],'og')
line([0 0],get(gca,'ylim'), ...
'linestyle','--','color','k')
line(get(gca,'xlim'),repmat(MswDiff,[1 2]), ...
'linestyle','--','color','k')
line(get(gca,'xlim'),repmat(MswDiff+thresSDswD(2)*SDswDiff,[1 2]),...
'linestyle','--','color','r')
line(get(gca,'xlim'),repmat(MswDiff-thresSDswD(2)*SDswDiff,[1 2]), ...
'linestyle','--','color','r')
xlabel('Samples'),
ylabel('A.U.'),
title(sprintf(['Sharp Wave Difference Magnitude\n',...
'Sharp Wave Duration:%5.1f'],durSWR*(1000/SR))), hold off
subplot(2,3,3), hist(swDiff_chk,100), ...
line(repmat(MswDiff+thresSDswD(2)*SDswDiff, [1 2]),get(gca,'ylim'), ...
'linestyle','--','color','k')
line(repmat(swDiff_chk(bound+1), [1 2]),get(gca,'ylim'), ...
'linestyle','--','color','r')
xlabel('A.U.'),
ylabel('Counts')
title(sprintf(['Histogram of Sharp Wave Difference Magnitude\n', ...
'%3.1f SD threshold (black) and detected value (red)'],thresSDswD(2)))
hold off
h2 = subplot(2,3,[4 5]);
plot(timeline,ripPow_chk), axis tight, hold on
plot([timeline(startSWR) timeline(bound+stopSWR)],...
[ripPow_chk(startSWR) ripPow_chk(bound+stopSWR)],'og')
line([0 0],get(gca,'ylim'), ...
'linestyle','--','color','k')
line(get(gca,'xlim'),repmat(MripPow,[1 2]), ...
'linestyle','--','color','k')
line(get(gca,'xlim'),repmat(MripPow+thresSDrip(2)*SDripPow,[1 2]), ...
'linestyle','--','color','r')
line(get(gca,'xlim'),repmat(MripPow-thresSDrip(2)*SDripPow,[1 2]), ...
'linestyle','--','color','r')
xlabel('Samples'),ylabel('A.U.'),title('Rip Pow'), hold off
subplot(2,3,6), hist(ripPow_chk,100), ...
line(repmat(MripPow+thresSDrip(2)*SDripPow, [1 2]),get(gca,'ylim'), ...
'linestyle','--','color','k')
line(repmat(ripPow_chk(bound+1), [1 2]),get(gca,'ylim'), ...
'linestyle','--','color','r')
xlabel('A.U.'),
ylabel('Counts')
title(sprintf(['Histogram of Ripple Power\n', ...
'%3.1f SD threshold (black) and detected value (red)'],thresSDrip(2)))
hold off
linkaxes([h1 h2],'x')
%         keyboard
pause(.1)
clf
end
end
ii
ii=ii+1
set(gcf,'Name',sprintf('Detection # %d out of %d',ii,Nswr))
figure
DEBUG = 1;
DEBUG = true
if DEBUG
set(gcf,'Name',sprintf('Detection # %d out of %d',ii,Nswr))
end
% threshold feature distributions to localize about detection
% 1) detected sharp wave difference magnitude must be > thresSD locally
swDiff_chk   = swDiff(SWRs(ii)-bound:SWRs(ii)+bound);
MswDiff  = median(swDiff_chk);
SDswDiff = std(swDiff_chk);
MswDiff
SDswDiff
figure; plot(swDiff_chk);
plot((-bound:bound)./SR,swDiff_chk);
4/0.61
in = InIntervals(tl,interval);
figure; plot(tl(in),swDiffp(in));
figure; plot(tl(in),swDiff0(in));
figure; plot(tl(in),swDiff(in));
plot((-bound:bound)./SR + t(ii),swDiff_chk);
xlim(interval);
MswDiff
PlotHVLines(MswDiff,'h');
PlotHVLines(MswDiff+[-1 1]*SDswDiff,'h','r');
swDiff_chk(bound+1) < MswDiff + thresSDswD(2)*SDswDiff && useSPW
swDiff_chk(bound+1)
bound
PlotHVLines(swDiff_chk(bound+1),'h','y');
MswDiff + thresSDswD(2)*SDswDiff
PlotHVLines(MswDiff+[-1 1]*thresSDswD(2)*SDswDiff,'h','r');
useSPW
% 2) detected ripple power must be > thresSD locally
ripPow_chk    = ripPower0(SWRs(ii)-bound:SWRs(ii)+bound);
MripPow   = median(ripPow_chk);
SDripPow  = std(ripPow_chk);
if ripPowerAll(SWRid(ii)) < MripPow + thresSDrip(2)*SDripPow
if DEBUG
disp('Ripple Power too low locally')
end
SWRindF(SWRid(ii)) = false;
continue
end
% 2) detected ripple power must be > thresSD locally
ripPow_chk    = ripPower0(SWRs(ii)-bound:SWRs(ii)+bound);
MripPow   = median(ripPow_chk);
SDripPow  = std(ripPow_chk);
if ripPowerAll(SWRid(ii)) < MripPow + thresSDrip(2)*SDripPow
if DEBUG
disp('Ripple Power too low locally')
end
SWRindF(SWRid(ii)) = false;
end
figure; plot((-bound:bound)./SR + t(ii),ripPow_chk);
PlotHVLines(MripPow + thresSDrip(2)*SDripPow,'h','r');
xlim(interval);
PlotHVLines(ripPowerAll(SWRid(ii)) ,'h',y');
PlotHVLines(ripPowerAll(SWRid(ii)) ,'h','y');
ii < Nswr
SWR_diff(ii)
minIsi
if SWR_diff(ii) < minIsi
if DEBUG
disp('IEI is too small')
end
SWRindF(SWRid(ii)) = false;
end
% Estimate duration of detected Sharp-Wave event
startSWR = find(swDiff_chk(1:bound + 1) < MswDiff + thresSDswD(1)*SDswDiff,1,'last');
stopSWR  = find(swDiff_chk(bound + 1: end) < MswDiff + thresSDswD(1)*SDswDiff,1,'first');
durSWR   = (bound + 1 - startSWR) + stopSWR;
startSWR
ttt = ((-bound:bound)./SR + t(ii));
size(ttt)
size(swDiff_chk)
PlotIntervals(ttt([startSWR stopSWR]));
stopSWR
(bound + 1 - startSWR)
durSWR
durSWR/1250
PlotIntervals(ttt([startSWR stopSWR+bound+1]));
PlotIntervals(ttt([startSWR stopSWR+bound+1]),'color','r');
% Estimate duration of detected ripple power increase event
% find index of maximum ripple power increase around Sharp-Wave
[~,maxRpInd]  = max(ripPow_chk(bound+1-HalfWinSize:bound+1+HalfWinSize));
startRP       = find(ripPow_chk(1:bound - HalfWinSize + maxRpInd) ...
< MripPow + thresSDrip(1)*SDripPow,1,'last');
stopRP        = find(ripPow_chk(bound - HalfWinSize + maxRpInd: end) ...
< MripPow + thresSDrip(1)*SDripPow,1,'first');
durRP         = (bound - HalfWinSize + maxRpInd - startRP) + stopRP;
PlotIntervals(ttt([startRP stopRP+bound+1]),'color','r');
PlotIntervals(ttt([startRP stopRP+bound+1]),'color','b');
durRP/1250
durSWR < minDurSWs
durRP < minDurRPs
maxDurSWs
maxDurSWs/1250
durSWR > maxDurSWs && useSPW
for ii = 1:Nswr
if DEBUG
set(gcf,'Name',sprintf('Detection # %d out of %d',ii,Nswr))
end
% threshold feature distributions to localize about detection
% 1) detected sharp wave difference magnitude must be > thresSD locally
swDiff_chk   = swDiff(SWRs(ii)-bound:SWRs(ii)+bound);
MswDiff  = median(swDiff_chk);
SDswDiff = std(swDiff_chk);
if swDiff_chk(bound+1) < MswDiff + thresSDswD(2)*SDswDiff && useSPW
if DEBUG
disp('Sharp Wave Difference magnitude too small locally')
end
SWRindF(SWRid(ii)) = false;
continue
end
% 2) detected ripple power must be > thresSD locally
ripPow_chk    = ripPower0(SWRs(ii)-bound:SWRs(ii)+bound);
MripPow   = median(ripPow_chk);
SDripPow  = std(ripPow_chk);
if ripPowerAll(SWRid(ii)) < MripPow + thresSDrip(2)*SDripPow
if DEBUG
disp('Ripple Power too low locally')
end
SWRindF(SWRid(ii)) = false;
continue
end
% 3) detections can't be within a minimum inter-event interval thresh
if ii < Nswr
if SWR_diff(ii) < minIsi
if DEBUG
disp('IEI is too small')
end
SWRindF(SWRid(ii)) = false;
continue
end
end
% Estimate duration of detected Sharp-Wave event
startSWR = find(swDiff_chk(1:bound + 1) < MswDiff + thresSDswD(1)*SDswDiff,1,'last');
stopSWR  = find(swDiff_chk(bound + 1: end) < MswDiff + thresSDswD(1)*SDswDiff,1,'first');
durSWR   = (bound + 1 - startSWR) + stopSWR;
% Estimate duration of detected ripple power increase event
% find index of maximum ripple power increase around Sharp-Wave
[~,maxRpInd]  = max(ripPow_chk(bound+1-HalfWinSize:bound+1+HalfWinSize));
startRP       = find(ripPow_chk(1:bound - HalfWinSize + maxRpInd) ...
< MripPow + thresSDrip(1)*SDripPow,1,'last');
stopRP        = find(ripPow_chk(bound - HalfWinSize + maxRpInd: end) ...
< MripPow + thresSDrip(1)*SDripPow,1,'first');
durRP         = (bound - HalfWinSize + maxRpInd - startRP) + stopRP;
% 4) events associated with detections must exceed a minimum duration
if durSWR < minDurSWs && useSPW
if DEBUG
disp('Sharp-Wave Event Duration is too short.')
end
SWRindF(SWRid(ii)) = false;
continue
end
if durRP < minDurRPs
if DEBUG
disp('Ripple Event Duration is too short.')
end
SWRindF(SWRid(ii)) = false;
continue
end
% 5) sharp-wave must not exceed a maximum duration
if durSWR > maxDurSWs && useSPW
if DEBUG
disp('Sharp-Wave Event Duration is too long.')
end
SWRindF(SWRid(ii)) = false;
continue
end
% We have a valid detection
valid                     = valid + 1;
SWR_valid.Ts(valid,:)     = [SWRs(ii) (SWRs(ii)-(bound+1-min(startSWR,startRP))) (SWRs(ii)+max(stopSWR,stopRP)-1)];
swMax                     = swDiffAll(SWRid(ii));
ripMax                    = ripPowerAll(SWRid(ii));
SWR_valid.SwMax(valid,:)  = ...
[(swMax - MswDiff)/SDswDiff sum(swDiff_chk < swMax)/(2*bound + 1)];
SWR_valid.RipMax(valid,:) = ...
[(ripMax - MripPow)/SDripPow sum(ripPow_chk < ripMax)/(2*bound + 1)];
if DEBUG
h1 = subplot(2,3,[1 2]);
plot(timeline,swDiff_chk), axis tight, hold on
plot([timeline(startSWR) timeline(bound+stopSWR)],...
[swDiff_chk(startSWR) swDiff_chk(bound+stopSWR)],'og')
line([0 0],get(gca,'ylim'), ...
'linestyle','--','color','k')
line(get(gca,'xlim'),repmat(MswDiff,[1 2]), ...
'linestyle','--','color','k')
line(get(gca,'xlim'),repmat(MswDiff+thresSDswD(2)*SDswDiff,[1 2]),...
'linestyle','--','color','r')
line(get(gca,'xlim'),repmat(MswDiff-thresSDswD(2)*SDswDiff,[1 2]), ...
'linestyle','--','color','r')
xlabel('Samples'),
ylabel('A.U.'),
title(sprintf(['Sharp Wave Difference Magnitude\n',...
'Sharp Wave Duration:%5.1f'],durSWR*(1000/SR))), hold off
subplot(2,3,3), hist(swDiff_chk,100), ...
line(repmat(MswDiff+thresSDswD(2)*SDswDiff, [1 2]),get(gca,'ylim'), ...
'linestyle','--','color','k')
line(repmat(swDiff_chk(bound+1), [1 2]),get(gca,'ylim'), ...
'linestyle','--','color','r')
xlabel('A.U.'),
ylabel('Counts')
title(sprintf(['Histogram of Sharp Wave Difference Magnitude\n', ...
'%3.1f SD threshold (black) and detected value (red)'],thresSDswD(2)))
hold off
h2 = subplot(2,3,[4 5]);
plot(timeline,ripPow_chk), axis tight, hold on
plot([timeline(startSWR) timeline(bound+stopSWR)],...
[ripPow_chk(startSWR) ripPow_chk(bound+stopSWR)],'og')
line([0 0],get(gca,'ylim'), ...
'linestyle','--','color','k')
line(get(gca,'xlim'),repmat(MripPow,[1 2]), ...
'linestyle','--','color','k')
line(get(gca,'xlim'),repmat(MripPow+thresSDrip(2)*SDripPow,[1 2]), ...
'linestyle','--','color','r')
line(get(gca,'xlim'),repmat(MripPow-thresSDrip(2)*SDripPow,[1 2]), ...
'linestyle','--','color','r')
xlabel('Samples'),ylabel('A.U.'),title('Rip Pow'), hold off
subplot(2,3,6), hist(ripPow_chk,100), ...
line(repmat(MripPow+thresSDrip(2)*SDripPow, [1 2]),get(gca,'ylim'), ...
'linestyle','--','color','k')
line(repmat(ripPow_chk(bound+1), [1 2]),get(gca,'ylim'), ...
'linestyle','--','color','r')
xlabel('A.U.'),
ylabel('Counts')
title(sprintf(['Histogram of Ripple Power\n', ...
'%3.1f SD threshold (black) and detected value (red)'],thresSDrip(2)))
hold off
linkaxes([h1 h2],'x')
%         keyboard
pause(.1)
clf
end
end
SWR_valid.Ts     = SWR_valid.Ts(1:valid,:);
SWR_valid.SwMax  = SWR_valid.SwMax(1:valid,:);
SWR_valid.RipMax = SWR_valid.RipMax(1:valid,:);
ii
DEBUG = false
for ii = 1:Nswr
if DEBUG
set(gcf,'Name',sprintf('Detection # %d out of %d',ii,Nswr))
end
% threshold feature distributions to localize about detection
% 1) detected sharp wave difference magnitude must be > thresSD locally
swDiff_chk   = swDiff(SWRs(ii)-bound:SWRs(ii)+bound);
MswDiff  = median(swDiff_chk);
SDswDiff = std(swDiff_chk);
if swDiff_chk(bound+1) < MswDiff + thresSDswD(2)*SDswDiff && useSPW
if DEBUG
disp('Sharp Wave Difference magnitude too small locally')
end
SWRindF(SWRid(ii)) = false;
continue
end
% 2) detected ripple power must be > thresSD locally
ripPow_chk    = ripPower0(SWRs(ii)-bound:SWRs(ii)+bound);
MripPow   = median(ripPow_chk);
SDripPow  = std(ripPow_chk);
if ripPowerAll(SWRid(ii)) < MripPow + thresSDrip(2)*SDripPow
if DEBUG
disp('Ripple Power too low locally')
end
SWRindF(SWRid(ii)) = false;
continue
end
% 3) detections can't be within a minimum inter-event interval thresh
if ii < Nswr
if SWR_diff(ii) < minIsi
if DEBUG
disp('IEI is too small')
end
SWRindF(SWRid(ii)) = false;
continue
end
end
% Estimate duration of detected Sharp-Wave event
startSWR = find(swDiff_chk(1:bound + 1) < MswDiff + thresSDswD(1)*SDswDiff,1,'last');
stopSWR  = find(swDiff_chk(bound + 1: end) < MswDiff + thresSDswD(1)*SDswDiff,1,'first');
durSWR   = (bound + 1 - startSWR) + stopSWR;
% Estimate duration of detected ripple power increase event
% find index of maximum ripple power increase around Sharp-Wave
[~,maxRpInd]  = max(ripPow_chk(bound+1-HalfWinSize:bound+1+HalfWinSize));
startRP       = find(ripPow_chk(1:bound - HalfWinSize + maxRpInd) ...
< MripPow + thresSDrip(1)*SDripPow,1,'last');
stopRP        = find(ripPow_chk(bound - HalfWinSize + maxRpInd: end) ...
< MripPow + thresSDrip(1)*SDripPow,1,'first');
durRP         = (bound - HalfWinSize + maxRpInd - startRP) + stopRP;
% 4) events associated with detections must exceed a minimum duration
if durSWR < minDurSWs && useSPW
if DEBUG
disp('Sharp-Wave Event Duration is too short.')
end
SWRindF(SWRid(ii)) = false;
continue
end
if durRP < minDurRPs
if DEBUG
disp('Ripple Event Duration is too short.')
end
SWRindF(SWRid(ii)) = false;
continue
end
% 5) sharp-wave must not exceed a maximum duration
if durSWR > maxDurSWs && useSPW
if DEBUG
disp('Sharp-Wave Event Duration is too long.')
end
SWRindF(SWRid(ii)) = false;
continue
end
% We have a valid detection
valid                     = valid + 1;
SWR_valid.Ts(valid,:)     = [SWRs(ii) (SWRs(ii)-(bound+1-min(startSWR,startRP))) (SWRs(ii)+max(stopSWR,stopRP)-1)];
swMax                     = swDiffAll(SWRid(ii));
ripMax                    = ripPowerAll(SWRid(ii));
SWR_valid.SwMax(valid,:)  = ...
[(swMax - MswDiff)/SDswDiff sum(swDiff_chk < swMax)/(2*bound + 1)];
SWR_valid.RipMax(valid,:) = ...
[(ripMax - MripPow)/SDripPow sum(ripPow_chk < ripMax)/(2*bound + 1)];
if DEBUG
h1 = subplot(2,3,[1 2]);
plot(timeline,swDiff_chk), axis tight, hold on
plot([timeline(startSWR) timeline(bound+stopSWR)],...
[swDiff_chk(startSWR) swDiff_chk(bound+stopSWR)],'og')
line([0 0],get(gca,'ylim'), ...
'linestyle','--','color','k')
line(get(gca,'xlim'),repmat(MswDiff,[1 2]), ...
'linestyle','--','color','k')
line(get(gca,'xlim'),repmat(MswDiff+thresSDswD(2)*SDswDiff,[1 2]),...
'linestyle','--','color','r')
line(get(gca,'xlim'),repmat(MswDiff-thresSDswD(2)*SDswDiff,[1 2]), ...
'linestyle','--','color','r')
xlabel('Samples'),
ylabel('A.U.'),
title(sprintf(['Sharp Wave Difference Magnitude\n',...
'Sharp Wave Duration:%5.1f'],durSWR*(1000/SR))), hold off
subplot(2,3,3), hist(swDiff_chk,100), ...
line(repmat(MswDiff+thresSDswD(2)*SDswDiff, [1 2]),get(gca,'ylim'), ...
'linestyle','--','color','k')
line(repmat(swDiff_chk(bound+1), [1 2]),get(gca,'ylim'), ...
'linestyle','--','color','r')
xlabel('A.U.'),
ylabel('Counts')
title(sprintf(['Histogram of Sharp Wave Difference Magnitude\n', ...
'%3.1f SD threshold (black) and detected value (red)'],thresSDswD(2)))
hold off
h2 = subplot(2,3,[4 5]);
plot(timeline,ripPow_chk), axis tight, hold on
plot([timeline(startSWR) timeline(bound+stopSWR)],...
[ripPow_chk(startSWR) ripPow_chk(bound+stopSWR)],'og')
line([0 0],get(gca,'ylim'), ...
'linestyle','--','color','k')
line(get(gca,'xlim'),repmat(MripPow,[1 2]), ...
'linestyle','--','color','k')
line(get(gca,'xlim'),repmat(MripPow+thresSDrip(2)*SDripPow,[1 2]), ...
'linestyle','--','color','r')
line(get(gca,'xlim'),repmat(MripPow-thresSDrip(2)*SDripPow,[1 2]), ...
'linestyle','--','color','r')
xlabel('Samples'),ylabel('A.U.'),title('Rip Pow'), hold off
subplot(2,3,6), hist(ripPow_chk,100), ...
line(repmat(MripPow+thresSDrip(2)*SDripPow, [1 2]),get(gca,'ylim'), ...
'linestyle','--','color','k')
line(repmat(ripPow_chk(bound+1), [1 2]),get(gca,'ylim'), ...
'linestyle','--','color','r')
xlabel('A.U.'),
ylabel('Counts')
title(sprintf(['Histogram of Ripple Power\n', ...
'%3.1f SD threshold (black) and detected value (red)'],thresSDrip(2)))
hold off
linkaxes([h1 h2],'x')
%         keyboard
pause(.1)
clf
end
end
56
SWR_valid.Ts     = SWR_valid.Ts(1:valid,:);
SWR_valid.SwMax  = SWR_valid.SwMax(1:valid,:);
SWR_valid.RipMax = SWR_valid.RipMax(1:valid,:);
FIGS
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% TIME MINIMUM TROUGH %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filter pyramidal lfp in ripple band and get the envelope
lfpfilt = bz_Filter(lfp(:,1),'passband',ripBP,'filter','fir1','nyquist',SR/2);
signEnv = 2*lfpfilt.*lfpfilt;
envelope = sgolayfilt(signEnv,4,101);
envelope = abs(sqrt(envelope));
% Compute timestap of the trough nearest to the envelope maximum
% That timestamp will be the middle timestamp of the ripple output
TsTrough = SWR_valid.Ts(:,1)*0;
parfor irip = 1:size(SWR_valid.Ts,1)
idxs = SWR_valid.Ts(irip,2):SWR_valid.Ts(irip,3);
[~, imax] = max(envelope(idxs));
imax = imax + idxs(1)-1;
idxsmax = round([imax-0.01*SR:imax+0.01*SR]);
[~, imin] = min(lfp(idxsmax));
imin = imin+idxsmax(1)-1;
TsTrough(irip) = imin;
%     plot((idxs-imin)/SR, lfp(idxs,1), 'k', 'linewidth', 2)
%     hold on
%     plot((idxs-imin)/SR, lfpfilt(idxs))
%     plot((idxs-imin)/SR, 1500*signEnv(idxs)/max(signEnv(idxs)))
%     plot(0/SR, lfp(imin), '*')
%     hold off
%     xlim([-0.05 0.05])
%     pause
end
% Keep ripples whose computed troughs are inside the start and end timestamps
idxsIn = find(TsTrough>SWR_valid.Ts(:,2) & TsTrough<SWR_valid.Ts(:,3));
% Compute instantaneous phase and amplitude
h = hilbert(lfpfilt);
phase = angle(h);
amplitude = abs(h);
unwrapped = unwrap(phase);
% Compute instantaneous frequency
frequency = bz_Diff(medfilt1(unwrapped,12),[0:length(lfpfilt)-1]'/SR,'smooth',0);
frequency = frequency/(2*pi);
%%%%%%%%%%%%%%%%%%%%%%%
%%% GENERATE OUTPUT %%%
%%%%%%%%%%%%%%%%%%%%%%%
% 1) swr detections
% - Start and end timestamps of each SWR
SWR.timestamps = SWR_valid.Ts(idxsIn,[2 3])/SR;
% Duration
SWR.duration = SWR.timestamps(:,2) - SWR.timestamps(:,1);
% amplitude at trough
SWR.amplitude = amplitude(TsTrough(idxsIn));
% Frequency at trough
SWR.frequency = frequency(TsTrough(idxsIn));
% - Peak power of ripple timestamp of each SWR
SWR.peaks = TsTrough(idxsIn)/SR; %SWR_valid.Ts(:,1)/SR;
% - Peak power values of each SWR
SWR.peakNormedPower = [];
% - Name of this function
SWR.detectorName = 'DetectSWR';
% - Standard deviation used to threshold ripple power
SWR.stdev = [thresSDrip(2) thresSDswD(2)];
SWR.SwMax  = SWR_valid.SwMax(idxsIn);
SWR.RipMax = SWR_valid.RipMax(idxsIn);
% Make sure all ripple epochs are unique.
% check start ts
[~,unique_idx,~] = unique(SWR.timestamps(:,1));
SWR.timestamps = SWR.timestamps(unique_idx,:);
SWR.duration = SWR.duration(unique_idx,:);
SWR.amplitude = SWR.amplitude(unique_idx,:);
SWR.frequency = SWR.frequency(unique_idx,:);
SWR.peaks = SWR.peaks(unique_idx,:);
SWR.SwMax = SWR.SwMax(unique_idx,:);
SWR.RipMax = SWR.RipMax(unique_idx,:);
SWR_valid.Ts = SWR_valid.Ts(unique_idx,:);
SWR_valid.SwMax = SWR_valid.SwMax(unique_idx,:);
SWR_valid.RipMax = SWR_valid.RipMax(unique_idx,:);
% check peak ts
[~,unique_idx,~] = unique(SWR.peaks);
SWR.timestamps = SWR.timestamps(unique_idx,:);
SWR.duration = SWR.duration(unique_idx,:);
SWR.amplitude = SWR.amplitude(unique_idx,:);
SWR.frequency = SWR.frequency(unique_idx,:);
SWR.peaks = SWR.peaks(unique_idx,:);
SWR.SwMax = SWR.SwMax(unique_idx,:);
SWR.RipMax = SWR.RipMax(unique_idx,:);
SWR_valid.Ts = SWR_valid.Ts(unique_idx,:);
SWR_valid.SwMax = SWR_valid.SwMax(unique_idx,:);
SWR_valid.RipMax = SWR_valid.RipMax(unique_idx,:);
size(idxsIn)
open idxsIn
7985/8098
% Compute instantaneous phase and amplitude
h = hilbert(lfpfilt);
phase = angle(h);
amplitude = abs(h);
unwrapped = unwrap(phase);
% Compute instantaneous frequency
frequency = bz_Diff(medfilt1(unwrapped,12),[0:length(lfpfilt)-1]'/SR,'smooth',0);
frequency = frequency/(2*pi);
%%%%%%%%%%%%%%%%%%%%%%%
%%% GENERATE OUTPUT %%%
%%%%%%%%%%%%%%%%%%%%%%%
% 1) swr detections
% - Start and end timestamps of each SWR
SWR.timestamps = SWR_valid.Ts(idxsIn,[2 3])/SR;
% Duration
SWR.duration = SWR.timestamps(:,2) - SWR.timestamps(:,1);
% amplitude at trough
SWR.amplitude = amplitude(TsTrough(idxsIn));
% Frequency at trough
SWR.frequency = frequency(TsTrough(idxsIn));
% - Peak power of ripple timestamp of each SWR
SWR.peaks = TsTrough(idxsIn)/SR; %SWR_valid.Ts(:,1)/SR;
% - Peak power values of each SWR
SWR.peakNormedPower = [];
% - Name of this function
SWR.detectorName = 'DetectSWR';
% - Standard deviation used to threshold ripple power
SWR.stdev = [thresSDrip(2) thresSDswD(2)];
SWR.SwMax  = SWR_valid.SwMax(idxsIn);
SWR.RipMax = SWR_valid.RipMax(idxsIn);
% Make sure all ripple epochs are unique.
% check start ts
[~,unique_idx,~] = unique(SWR.timestamps(:,1));
SWR.timestamps = SWR.timestamps(unique_idx,:);
SWR.duration = SWR.duration(unique_idx,:);
SWR.amplitude = SWR.amplitude(unique_idx,:);
SWR.frequency = SWR.frequency(unique_idx,:);
SWR.peaks = SWR.peaks(unique_idx,:);
SWR.SwMax = SWR.SwMax(unique_idx,:);
SWR.RipMax = SWR.RipMax(unique_idx,:);
SWR_valid.Ts = SWR_valid.Ts(unique_idx,:);
SWR_valid.SwMax = SWR_valid.SwMax(unique_idx,:);
SWR_valid.RipMax = SWR_valid.RipMax(unique_idx,:);
% check peak ts
[~,unique_idx,~] = unique(SWR.peaks);
SWR.timestamps = SWR.timestamps(unique_idx,:);
SWR.duration = SWR.duration(unique_idx,:);
SWR.amplitude = SWR.amplitude(unique_idx,:);
SWR.frequency = SWR.frequency(unique_idx,:);
SWR.peaks = SWR.peaks(unique_idx,:);
SWR.SwMax = SWR.SwMax(unique_idx,:);
SWR.RipMax = SWR.RipMax(unique_idx,:);
SWR_valid.Ts = SWR_valid.Ts(unique_idx,:);
SWR_valid.SwMax = SWR_valid.SwMax(unique_idx,:);
SWR_valid.RipMax = SWR_valid.RipMax(unique_idx,:);
size(idxsIn)
size(SWR_valid.Ts,1)
Nswr
SWR_valid.Ts     = zeros(size(SWRs,1),3);   % peak, start, stop in samples
SWR_valid.SwMax  = zeros(size(SWRs,1),2);   % z-score & percentile
SWR_valid.RipMax = zeros(size(SWRs,1),2);   % z-score & percentile
valid            = 0;
SWRindF          = SWRind;
for ii = 1:Nswr
if DEBUG
set(gcf,'Name',sprintf('Detection # %d out of %d',ii,Nswr))
end
% threshold feature distributions to localize about detection
% 1) detected sharp wave difference magnitude must be > thresSD locally
swDiff_chk   = swDiff(SWRs(ii)-bound:SWRs(ii)+bound);
MswDiff  = median(swDiff_chk);
SDswDiff = std(swDiff_chk);
if swDiff_chk(bound+1) < MswDiff + thresSDswD(2)*SDswDiff && useSPW
if DEBUG
disp('Sharp Wave Difference magnitude too small locally')
end
SWRindF(SWRid(ii)) = false;
continue
end
% 2) detected ripple power must be > thresSD locally
ripPow_chk    = ripPower0(SWRs(ii)-bound:SWRs(ii)+bound);
MripPow   = median(ripPow_chk);
SDripPow  = std(ripPow_chk);
if ripPowerAll(SWRid(ii)) < MripPow + thresSDrip(2)*SDripPow
if DEBUG
disp('Ripple Power too low locally')
end
SWRindF(SWRid(ii)) = false;
continue
end
% 3) detections can't be within a minimum inter-event interval thresh
if ii < Nswr
if SWR_diff(ii) < minIsi
if DEBUG
disp('IEI is too small')
end
SWRindF(SWRid(ii)) = false;
continue
end
end
% Estimate duration of detected Sharp-Wave event
startSWR = find(swDiff_chk(1:bound + 1) < MswDiff + thresSDswD(1)*SDswDiff,1,'last');
stopSWR  = find(swDiff_chk(bound + 1: end) < MswDiff + thresSDswD(1)*SDswDiff,1,'first');
durSWR   = (bound + 1 - startSWR) + stopSWR;
% Estimate duration of detected ripple power increase event
% find index of maximum ripple power increase around Sharp-Wave
[~,maxRpInd]  = max(ripPow_chk(bound+1-HalfWinSize:bound+1+HalfWinSize));
startRP       = find(ripPow_chk(1:bound - HalfWinSize + maxRpInd) ...
< MripPow + thresSDrip(1)*SDripPow,1,'last');
stopRP        = find(ripPow_chk(bound - HalfWinSize + maxRpInd: end) ...
< MripPow + thresSDrip(1)*SDripPow,1,'first');
durRP         = (bound - HalfWinSize + maxRpInd - startRP) + stopRP;
% 4) events associated with detections must exceed a minimum duration
if durSWR < minDurSWs && useSPW
if DEBUG
disp('Sharp-Wave Event Duration is too short.')
end
SWRindF(SWRid(ii)) = false;
continue
end
if durRP < minDurRPs
if DEBUG
disp('Ripple Event Duration is too short.')
end
SWRindF(SWRid(ii)) = false;
continue
end
% 5) sharp-wave must not exceed a maximum duration
if durSWR > maxDurSWs && useSPW
if DEBUG
disp('Sharp-Wave Event Duration is too long.')
end
SWRindF(SWRid(ii)) = false;
continue
end
% We have a valid detection
valid                     = valid + 1;
SWR_valid.Ts(valid,:)     = [SWRs(ii) (SWRs(ii)-(bound+1-min(startSWR,startRP))) (SWRs(ii)+max(stopSWR,stopRP)-1)];
swMax                     = swDiffAll(SWRid(ii));
ripMax                    = ripPowerAll(SWRid(ii));
SWR_valid.SwMax(valid,:)  = ...
[(swMax - MswDiff)/SDswDiff sum(swDiff_chk < swMax)/(2*bound + 1)];
SWR_valid.RipMax(valid,:) = ...
[(ripMax - MripPow)/SDripPow sum(ripPow_chk < ripMax)/(2*bound + 1)];
if DEBUG
h1 = subplot(2,3,[1 2]);
plot(timeline,swDiff_chk), axis tight, hold on
plot([timeline(startSWR) timeline(bound+stopSWR)],...
[swDiff_chk(startSWR) swDiff_chk(bound+stopSWR)],'og')
line([0 0],get(gca,'ylim'), ...
'linestyle','--','color','k')
line(get(gca,'xlim'),repmat(MswDiff,[1 2]), ...
'linestyle','--','color','k')
line(get(gca,'xlim'),repmat(MswDiff+thresSDswD(2)*SDswDiff,[1 2]),...
'linestyle','--','color','r')
line(get(gca,'xlim'),repmat(MswDiff-thresSDswD(2)*SDswDiff,[1 2]), ...
'linestyle','--','color','r')
xlabel('Samples'),
ylabel('A.U.'),
title(sprintf(['Sharp Wave Difference Magnitude\n',...
'Sharp Wave Duration:%5.1f'],durSWR*(1000/SR))), hold off
subplot(2,3,3), hist(swDiff_chk,100), ...
line(repmat(MswDiff+thresSDswD(2)*SDswDiff, [1 2]),get(gca,'ylim'), ...
'linestyle','--','color','k')
line(repmat(swDiff_chk(bound+1), [1 2]),get(gca,'ylim'), ...
'linestyle','--','color','r')
xlabel('A.U.'),
ylabel('Counts')
title(sprintf(['Histogram of Sharp Wave Difference Magnitude\n', ...
'%3.1f SD threshold (black) and detected value (red)'],thresSDswD(2)))
hold off
h2 = subplot(2,3,[4 5]);
plot(timeline,ripPow_chk), axis tight, hold on
plot([timeline(startSWR) timeline(bound+stopSWR)],...
[ripPow_chk(startSWR) ripPow_chk(bound+stopSWR)],'og')
line([0 0],get(gca,'ylim'), ...
'linestyle','--','color','k')
line(get(gca,'xlim'),repmat(MripPow,[1 2]), ...
'linestyle','--','color','k')
line(get(gca,'xlim'),repmat(MripPow+thresSDrip(2)*SDripPow,[1 2]), ...
'linestyle','--','color','r')
line(get(gca,'xlim'),repmat(MripPow-thresSDrip(2)*SDripPow,[1 2]), ...
'linestyle','--','color','r')
xlabel('Samples'),ylabel('A.U.'),title('Rip Pow'), hold off
subplot(2,3,6), hist(ripPow_chk,100), ...
line(repmat(MripPow+thresSDrip(2)*SDripPow, [1 2]),get(gca,'ylim'), ...
'linestyle','--','color','k')
line(repmat(ripPow_chk(bound+1), [1 2]),get(gca,'ylim'), ...
'linestyle','--','color','r')
xlabel('A.U.'),
ylabel('Counts')
title(sprintf(['Histogram of Ripple Power\n', ...
'%3.1f SD threshold (black) and detected value (red)'],thresSDrip(2)))
hold off
linkaxes([h1 h2],'x')
%         keyboard
pause(.1)
clf
end
end
SWR_valid.Ts     = SWR_valid.Ts(1:valid,:);
SWR_valid.SwMax  = SWR_valid.SwMax(1:valid,:);
SWR_valid.RipMax = SWR_valid.RipMax(1:valid,:);
if FIGS
% Evaluate fixed point Precision/Recall performance
figure('Color','w'),
% marginal histograms
[n1, ctr1] = hist(swDiffAll,1000);
[n2, ctr2] = hist(ripPowerAll,1000);
h1 = subplot(2,2,2);
plot(swDiffAll(idx1), ripPowerAll(idx1),'.r'), hold on
plot(swDiffAll(idx2), ripPowerAll(idx2),'.b'),
plot(swDiffAll(SWRindF),ripPowerAll(SWRindF),'og','markerfacecolor','none')
axis([min(swDiffAll(:)) max(swDiffAll(:)) min(ripPowerAll(:)) max(ripPowerAll(:))])
line(get(h1,'xlim'), repmat(thresL_Rip,[1 2]),'linestyle','--','color','k')
line(repmat(thresL_swD,[1 2]),get(h1,'ylim'),'linestyle','--','color','k')
title('Sharp Wave Difference vs. Ripple Power','fontsize',12)
xlabel('Sharp Wave Diff','fontsize',12)
ylabel('Ripple Power','fontsize',12)
h2 = subplot(2,2,4);
bar(ctr1,-n1,1),
axis([min(swDiffAll(:)) max(swDiffAll(:)) -max(n1)*1.1 0]),
axis('off'),
h3 = subplot(2,2,1);
barh(ctr2,-n2,1),
axis([-max(n2)*1.1 0 min(ripPowerAll(:)) max(ripPowerAll(:))])
axis('off'),
set(h1,'Position',[0.25 0.35 0.55 0.55],'fontsize',12);
set(h2,'Position',[.25 .1 .55 .15],'fontsize',12);
set(h3,'Position',[.05 .35 .15 .55],'fontsize',12);
string = sprintf(['K-means Clustering:\n\n',...
'    Center1: X  = %7.2f, Y = %7.2f\n', ...
'    Center2: X  = %7.2f, Y = %7.2f\n\n', ...
'Number of Detected SWR:%5d\n'],...
Ctrs(1,1),Ctrs(1,2),Ctrs(2,1),Ctrs(2,2), valid);
annotation('textbox',[0.81 0.70 .1 .1],'String',string, ...
'edgecolor','none','fontsize',12)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% TIME MINIMUM TROUGH %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filter pyramidal lfp in ripple band and get the envelope
lfpfilt = bz_Filter(lfp(:,1),'passband',ripBP,'filter','fir1','nyquist',SR/2);
signEnv = 2*lfpfilt.*lfpfilt;
envelope = sgolayfilt(signEnv,4,101);
envelope = abs(sqrt(envelope));
% Compute timestap of the trough nearest to the envelope maximum
% That timestamp will be the middle timestamp of the ripple output
TsTrough = SWR_valid.Ts(:,1)*0;
parfor irip = 1:size(SWR_valid.Ts,1)
idxs = SWR_valid.Ts(irip,2):SWR_valid.Ts(irip,3);
[~, imax] = max(envelope(idxs));
imax = imax + idxs(1)-1;
idxsmax = round([imax-0.01*SR:imax+0.01*SR]);
[~, imin] = min(lfp(idxsmax));
imin = imin+idxsmax(1)-1;
TsTrough(irip) = imin;
%     plot((idxs-imin)/SR, lfp(idxs,1), 'k', 'linewidth', 2)
%     hold on
%     plot((idxs-imin)/SR, lfpfilt(idxs))
%     plot((idxs-imin)/SR, 1500*signEnv(idxs)/max(signEnv(idxs)))
%     plot(0/SR, lfp(imin), '*')
%     hold off
%     xlim([-0.05 0.05])
%     pause
end
size(SWR_valid.Ts,1)
irip
size(SWR_valid.Ts,1)
irip=  ans
idxs = SWR_valid.Ts(irip,2):SWR_valid.Ts(irip,3);
idxs = SWR_valid.Ts(irip,2):SWR_valid.Ts(irip,3);
[~, imax] = max(envelope(idxs));
imax = imax + idxs(1)-1
idxsmax = round([imax-0.01*SR:imax+0.01*SR]);
0.01*SR
TsTrough2 = SWR_valid.Ts(:,1)*0;
parfor irip = 1:size(SWR_valid.Ts,1)
idxs = SWR_valid.Ts(irip,1):SWR_valid.Ts(irip,3);
[~, imax] = max(envelope(idxs));
imax = imax + idxs(1)-1;
idxsmax = round([imax-0.01*SR:imax+0.01*SR]);
[~, imin] = min(lfp(idxsmax));
imin = imin+idxsmax(1)-1;
TsTrough2(irip) = imin;
%     plot((idxs-imin)/SR, lfp(idxs,1), 'k', 'linewidth', 2)
%     hold on
%     plot((idxs-imin)/SR, lfpfilt(idxs))
%     plot((idxs-imin)/SR, 1500*signEnv(idxs)/max(signEnv(idxs)))
%     plot(0/SR, lfp(imin), '*')
%     hold off
%     xlim([-0.05 0.05])
%     pause
end
sum(TsTrough~=TsTrough2)
Portion(TsTrough>SWR_valid.Ts(:,2) & TsTrough<SWR_valid.Ts(:,3))
Portion(TsTrough2>SWR_valid.Ts(:,2) & TsTrough2<SWR_valid.Ts(:,3))
mmax(diff(SWR_valid.Ts(TsTrough2>SWR_valid.Ts(:,2) & TsTrough2<SWR_valid.Ts(:,3),[1 3]),[],2))
mmax(diff(SWR_valid.Ts(~(TsTrough2>SWR_valid.Ts(:,2) & TsTrough2<SWR_valid.Ts(:,3),[1 3])),[],2))
mmax(diff(SWR_valid.Ts(~(TsTrough2>SWR_valid.Ts(:,2) & TsTrough2<SWR_valid.Ts(:,3)),[1 3]),[],2))
mmax(diff(SWR_valid.Ts(~(TsTrough2>SWR_valid.Ts(:,2) & TsTrough2<SWR_valid.Ts(:,3)),[1 3]),[],2))/1250*1000
[SWRs(ii) (SWRs(ii)-(bound+1-min(startSWR,startRP))) (SWRs(ii)+max(stopSWR,stopRP)-1)]
[SWRs(ii) (SWRs(ii)-(bound+1-min(startSWR))) (SWRs(ii)+max(stopSWR)-1)];
[SWRs(ii) (SWRs(ii)-(bound+1-min(startSWR))) (SWRs(ii)+max(stopSWR)-1)]
[17842099    17842081    17842127] - 17842099
[ 17842099    17842081    17842153]-17842099
% Keep ripples whose computed troughs are inside the start and end timestamps
idxsIn = find(TsTrough>SWR_valid.Ts(:,2) & TsTrough<SWR_valid.Ts(:,3));
% Compute instantaneous phase and amplitude
h = hilbert(lfpfilt);
phase = angle(h);
amplitude = abs(h);
unwrapped = unwrap(phase);
% Compute instantaneous frequency
frequency = bz_Diff(medfilt1(unwrapped,12),[0:length(lfpfilt)-1]'/SR,'smooth',0);
frequency = frequency/(2*pi);
%%%%%%%%%%%%%%%%%%%%%%%
%%% GENERATE OUTPUT %%%
%%%%%%%%%%%%%%%%%%%%%%%
% 1) swr detections
% - Start and end timestamps of each SWR
SWR.timestamps = SWR_valid.Ts(idxsIn,[2 3])/SR;
% Duration
SWR.duration = SWR.timestamps(:,2) - SWR.timestamps(:,1);
% amplitude at trough
SWR.amplitude = amplitude(TsTrough(idxsIn));
% Frequency at trough
SWR.frequency = frequency(TsTrough(idxsIn));
% - Peak power of ripple timestamp of each SWR
SWR.peaks = TsTrough(idxsIn)/SR; %SWR_valid.Ts(:,1)/SR;
% - Peak power values of each SWR
SWR.peakNormedPower = [];
% - Name of this function
SWR.detectorName = 'DetectSWR';
% - Standard deviation used to threshold ripple power
SWR.stdev = [thresSDrip(2) thresSDswD(2)];
SWR.SwMax  = SWR_valid.SwMax(idxsIn);
SWR.RipMax = SWR_valid.RipMax(idxsIn);
% Make sure all ripple epochs are unique.
% check start ts
[~,unique_idx,~] = unique(SWR.timestamps(:,1));
SWR.timestamps = SWR.timestamps(unique_idx,:);
SWR.duration = SWR.duration(unique_idx,:);
SWR.amplitude = SWR.amplitude(unique_idx,:);
SWR.frequency = SWR.frequency(unique_idx,:);
SWR.peaks = SWR.peaks(unique_idx,:);
SWR.SwMax = SWR.SwMax(unique_idx,:);
SWR.RipMax = SWR.RipMax(unique_idx,:);
SWR_valid.Ts = SWR_valid.Ts(unique_idx,:);
SWR_valid.SwMax = SWR_valid.SwMax(unique_idx,:);
SWR_valid.RipMax = SWR_valid.RipMax(unique_idx,:);
% check peak ts
[~,unique_idx,~] = unique(SWR.peaks);
SWR.timestamps = SWR.timestamps(unique_idx,:);
SWR.duration = SWR.duration(unique_idx,:);
SWR.amplitude = SWR.amplitude(unique_idx,:);
SWR.frequency = SWR.frequency(unique_idx,:);
SWR.peaks = SWR.peaks(unique_idx,:);
SWR.SwMax = SWR.SwMax(unique_idx,:);
SWR.RipMax = SWR.RipMax(unique_idx,:);
SWR_valid.Ts = SWR_valid.Ts(unique_idx,:);
SWR_valid.SwMax = SWR_valid.SwMax(unique_idx,:);
SWR_valid.RipMax = SWR_valid.RipMax(unique_idx,:);
%Write event file: Generate Output and Write out Event file
if EVENTFILE
rippleFiles = dir('*.R*.evt');
if isempty(rippleFiles)
fileN = 1;
else
%set file index to next available value\
pat = '.R[0-9].';
fileN = 0;
for ii = 1:length(rippleFiles)
token  = regexp(rippleFiles(ii).name,pat);
val    = str2double(rippleFiles(ii).name(token+2:token+4));
fileN  = max([fileN val]);
end
fileN = fileN + 1;
end
fid = fopen(sprintf('%s%s%s.R%02d.evt',pathname,filesep,filename,fileN),'w');
% convert detections to milliseconds
SWR_final   = SWR_valid.Ts*(1000/SR);
fprintf(1,'Writing event file ...\n');
for ii = 1:size(SWR_final,1)
fprintf(fid,'%9.1f\tstart\n',SWR_final(ii,2));
fprintf(fid,'%9.1f\tpeak\n',SWR_final(ii,1));
fprintf(fid,'%9.1f\tstop\n',SWR_final(ii,3));
end
fclose(fid);
end
% 3) analysis params
params.filename     = filename;
params.Channels     = Channels;
params.SR           = SR;
params.nChan        = nChan;
params.swBP         = swBP;
params.ripBP        = ripBP;
params.per_thresswD = per_thresswD;
params.per_thresRip = per_thresRip;
params.WinSize      = WinSize;
params.Ns_chk       = Ns_chk;
params.thresSDswD   = thresSDswD;
params.thresSDrip   = thresSDrip;
params.minIsi       = minIsi;
params.minDurSW     = minDurSW;
params.maxDurSW     = maxDurSW;
params.minDurRP     = minDurRP;
params.EVENTFILE    = EVENTFILE;
params.TRAINING     = TRAINING;
params.DEBUG        = DEBUG;
SWR.detectorinfo.detectionparms = params;
SWR.detectorinfo.detectorname = 'bz_DetectSWR';
SWR.detectorinfo.detectiondate = datestr(datetime('now'),'dd-mmm-yyyy');
SWR.detectorinfo.detectionintervals = Epochs;
% write out log file
logAnalysisFile(mfilename('fullpath'),pathname);
ripples = SWR;
if saveMat
save([Filebase '.ripples.events.mat'],'ripples')
end
r = [ripples.timestamps(:,1) ripples.peaks ripples.timestamps(:,2) ripples.RipMax ripples.SwMax];
load('day3.MergePoints.events.mat')
pre = InIntervals(r,sleep(1,:));
sleep = ConsolidateIntervals(MergePoints.timestamps(cellfun(@(x) any(strfind(lower(x),'sleep')),MergePoints.foldernames),:));
pre = InIntervals(r,sleep(1,:));
post = InIntervals(r,sleep(end,:));
figure;  anovabar(diff(r(:,[1 3]),[],2)>0.07,[-pre+post]);
anovabar(diff(r(:,[1 3]),[],2)>0.1,[-pre+post]);
anovabar(diff(r(:,[1 3]),[],2),[-pre+post]);
anovabar(diff(r(:,[1 3]),[],2)>0.1,[-pre+post]);
r0 = r;
load('day3.ripples.events.mat')
r = [ripples.timestamps(:,1) ripples.peaks ripples.timestamps(:,2) ripples.RipMax ripples.SwMax];
figure; PETH(r(:,1),r0(:,1));
PETH(r(:,1),r0(:,1),'durations',[-1 1]*0.2);
PETH(r(:,1),r0(:,1),'durations',[-1 1]*0.1);
PETH(r(:,3),r0(:,3),'durations',[-1 1]*0.1);
PETH(r(:,2),r0(:,2),'durations',[-1 1]*0.1);
matched = FindClosest(r(:,2),r0(:,2));
size(matched)
size(r)
plot(diff(r(matched,[1 3]),[],2),diff(r0(:,[1 3]),[],2));
plot(diff(r(matched,[1 3]),[],2),diff(r0(:,[1 3]),[],2),'.');
mean(diff(r(matched,[1 3]),[],2)-diff(r0(:,[1 3]),[],2))
hist((diff(r(matched,[1 3]),[],2)-diff(r0(:,[1 3]),[],2),100);
hist(diff(r(matched,[1 3]),[],2)-diff(r0(:,[1 3]),[],2),100);
plot(diff(r(matched,[1 3]),[],2),diff(r0(:,[1 3]),[],2),'.');
PlotHVLines(0.1,'v','k--');
[h,ht] = Dist(100,[diff(r(:,[1 3]),[],2),[-pre+post]+2],'grouped');
[h,ht] = Dist(100,[diff(r0(:,[1 3]),[],2),[-pre+post]+2],'grouped');
figure;  plot(ht,Smooth((h),[0 0]))
figure;  plot(ht*1000,Smooth((h),[0 0]))
figure;  plot(ht*1000,Smooth((h),[0 0])); set(gca,'ytick',[]); xlabel('ripple duration (s)');
plot(ht*1000,Smooth((h),[0 0]),'linewidth',2); set(gca,'ytick',[]); xlabel('ripple duration (s)');
legend('pre sleep','task','post sleep');
SideAxes(gca,'right',0.5);
pre = InIntervals(r,sleep(1,:));
post = InIntervals(r,sleep(end,:));
[h,ht] = Dist(100,[diff(r(:,[1 3]),[],2),[-pre+post]+2],'grouped');
plot(ht*1000,Smooth((h),[0 0]),'linewidth',2); set(gca,'ytick',[]); xlabel('ripple duration (s)');
title('New start/stop');
title('Sharp wave start/stop (ignore ripple)');
figure; anovabar(diff(r(:,[1 3]),[],2),[-pre+post]);
anovabar(diff(r(:,[1 3]),[],2),[-pre+post]);
anovabar(diff(r(:,[1 3])>0.1,[],2),[-pre+post]);
anovabar(diff(r(:,[1 3]),[],2)>0.1,[-pre+post]);
ripples = SWR;
ripples = DetectSWR([rippleChannels.Ripple_Channel rippleChannels.Sharpwave_Channel rippleChannels.Noise_Channel],'basepath',basepath,'saveMat',true,'forceDetect',true,'check',true);
load('DetectSWR_manual_scoring.mat')
clf
plot(swDiffAll(idx2),ripPowerAll(idx2),'k.','markersize',1); hold on
plot(swDiffAll(idx1),ripPowerAll(idx1),'r.','markersize',1); legend('non-ripples','ripples');
disp(['Please review the figure and change the idx1 (ripples) ' newline 'and idx2 (calm non-ripple periods) groups.' newline ...
'Note that noisy periods need not participate in either group.']);
dbcont
hold all
plot(swDiffAll(SWRind),ripPowerAll(SWRind),'b.','markersize',1);
check
ripples = DetectSWR([rippleChannels.Ripple_Channel rippleChannels.Sharpwave_Channel rippleChannels.Noise_Channel],'basepath',basepath,'saveMat',true,'forceDetect',true,'check',true,'useSPW',false);
load('DetectSWR_manual_scoring.mat')
clf
plot(swDiffAll(idx2),ripPowerAll(idx2),'k.','markersize',1); hold on
plot(swDiffAll(idx1),ripPowerAll(idx1),'r.','markersize',1); legend('non-ripples','ripples');
disp(['Please review the figure and change the idx1 (ripples) ' newline 'and idx2 (calm non-ripple periods) groups.' newline ...
'Note that noisy periods need not participate in either group.']);
dbcont
r2 = r;
r = [ripples.timestamps(:,1) ripples.peaks ripples.timestamps(:,2) ripples.RipMax ripples.SwMax];
open r
open r2
open r0
pre = InIntervals(r,sleep(1,:));
post = InIntervals(r,sleep(end,:));
figure; anovabar(diff(r(:,[1 3]),[],2),[-pre+post]);
anovabar(diff(r(:,[1 3]),[],2)>0.1,[-pre+post]);
[h,ht] = Dist(100,[diff(r(:,[1 3]),[],2),[-pre+post]+2],'grouped');
plot(ht*1000,Smooth((h),[0 0]),'linewidth',2); set(gca,'ytick',[]); xlabel('ripple duration (s)');
clf
PETH(lfp,r(:,1));
hold all
PETH(lfp,r0(:,1));
PETH(lfp,r2(:,1));
matched = FindClosest(r(:,2),r0(:,2));
open matched
d = r(matched,2)-r0(:,2);
open d
size(r0)
size(r2)
size(r)
size(matched)
mmax(d)
Portion(d==0)
matched = FindClosest(r0(:,2),r(:,2));
d = r(matched,2)-r0(:,2);
d = r0(matched,2)-r(:,2);
Portion(d==0)
figure; anovabar(r(:,5),d==0);
anovabar(r(:,4),d==0);
anovabar(r(:,[4 5]),d==0);
anovabar(r(:,[3 4 5]),d==0);
anovabar(r(:,[4 5]),d==0);
Dist(100,[r(:,4) (d==0)+1],'grouped');
Dist(100,[r(:,5) (d==0)+1],'grouped');
PETH(lfp,r(d~=0,1));
clf
hist(d,100)
hist(abs(d),100)
hist(abs(d),1000)
new = r(abs(d)>0.1,:);
size(new)
SaveCustomEvents('temp.tem.evt',new(:,1:3),{'ripple start','ripple peak','ripple stop'});
anovabar(diff(r(:,[1 3]),[],2)>0.1,[-pre+post]);
anovabar(diff(r(d==0,[1 3]),[],2)>0.1,[-pre+post]);
g = [-pre+post];
anovabar(diff(r(d==0,[1 3]),[],2)>0.1,g(d==0));
anovabar(diff(r(d~=0,[1 3]),[],2)>0.1,g(d~=0));
Dist(100,[diff(r(:,[1 3]),[],2) (d==0)+1],'grouped');
anovabar(diff(r(d~=0,[1 3]),[],2)>0.1,g(d~=0));
ok = d~=0 & g~=0;
anovabar(diff(r(ok,[1 3]),[],2)>0.1,g(ok));
anovabar(diff(r(ok,[1 3]),[],2),g(ok));
ok = d~=0;
[h,ht] = Dist(100,[diff(r(ok,[1 3]),[],2),g(ok)+2],'grouped');
plot(ht*1000,Smooth((h),[0 0]),'linewidth',2); set(gca,'ytick',[]); xlabel('ripple duration (s)');
plot(ht*1000,Smooth((h),[1 0]),'linewidth',2); set(gca,'ytick',[]); xlabel('ripple duration (s)');
plot(ht*1000,Smooth(cumsum(h),[0 0]),'linewidth',2); set(gca,'ytick',[]); xlabel('ripple duration (s)');
ok = d==0;
[h,ht] = Dist(100,[diff(r(ok,[1 3]),[],2),g(ok)+2],'grouped');
plot(ht*1000,Smooth(cumsum(h),[0 0]),'linewidth',2); set(gca,'ytick',[]); xlabel('ripple duration (s)');
Dist(100,[r(:,4) (d==0)+1],'grouped');
Dist(100,[r(:,3) (d==0)+1],'grouped');
Dist(100,[r(:,5) (d==0)+1],'grouped');
Portion(r(:,5)>2)
Portion(r(ok,5)>2)
Portion(r(~ok,5)>2)
Dist(1000,[r(:,5) (d==0)+1],'grouped');
Portion(r(~ok,5)>2.5)
Portion(r(ok,5)>2.5)
DensityMap(r(:,4),r(:,5));
DensityMap(r(:,4),r(:,5),'smooth',0);
ok = g==-1; DensityMap(r(ok,4),r(ok,5),'smooth',0);
ok = g==0; DensityMap(r(ok,4),r(ok,5),'smooth',0);
ok = g==-1; DensityMap(r(ok,4),r(ok,5),'smooth',0);
ok = g==1; DensityMap(r(ok,4),r(ok,5),'smooth',0);
b = Bins(r(:,4),20);
b1 = Bins(r(:,4),20); b2 = Bins(r(:,5),20); b = [b1 b2];
b1 = Bin(r(:,4),20); b2 = Bin(r(:,5),20); b = [b1 b2];
for i=1:3,
qs{i,1} = Accumulate(b,1,'size',max(b));
end
open qs
for i=1:3,
qs{i,1} = Accumulate(b,1,'size',max(b));
subplot(2,2,i);
PlotColorMap(qs{i,1});
end
for i=1:3,
qs{i,1} = Accumulate(b(g==i,:),1,'size',max(b));
subplot(2,2,i);
PlotColorMap(qs{i,1});
end
open g
g = [-pre+post]+2;
for i=1:3,
qs{i,1} = Accumulate(b(g==i,:),1,'size',max(b));
subplot(2,2,i);
PlotColorMap(qs{i,1});
end
figure; ok = d==0;
anovabar(diff(r(ok,[1 3]),[],2)>0.1,g(ok));
ok = d~=0;
anovabar(diff(r(ok,[1 3]),[],2)>0.1,g(ok));
clear all
X = BatchReturn('OJR_longTraining.batch');
cd(X{1}0;
cd(X{1});
load('day7.channelinfo.ripples.mat')
load('day3.ripples.events.mat')
r = [ripples.timestamps(:,1) ripples.peaks ripples.timestamps(:,2) ripples.RipMax ripples.SwMax];
load('day3.MergePoints.events.mat')
sleep = ConsolidateIntervals(MergePoints.timestamps(cellfun(@(x) any(strfind(lower(x),'sleep')),MergePoints.foldernames),:));
pre = InIntervals(r,sleep(1,:));
post = InIntervals(r,sleep(end,:));
anovabar(r(:,5),[-pre post]);
sum(pre)
sum(post)
anovabar(r(1,5),[-pre post]);
anovabar(r(:,1),[-pre post]);
g = [-pre post]+2;
anovabar(r(:,1),g);
anovabar(diff(r(ok,[1 3]),[],2)>0.1,g(ok));
anovabar(diff(r(:,[1 3]),[],2)>0.1,[-pre+post]);
anovabar(r(:,1),[-pre+post]);
anovabar(r(:,5),[-pre+post]);
anovabar(r(:,4),[-pre+post]);
anovabar(r(:,5),[-pre+post]);
anovabar(r(:,5)>2.5,[-pre+post]);
anovabar(diff(r(:,[1 3]),[],2)>0.1,[-pre+post]);
anovabar(diff(r(:,[1 3]),[],2)>0.1,g);
g = [-pre+post]+2;
anovabar(diff(r(:,[1 3]),[],2)>0.1,g);
ok = r(:,5)<2.5; anovabar(diff(r(ok,[1 3]),[],2)>0.1,g(ok));
ok = r(:,5)>2.5; anovabar(diff(r(ok,[1 3]),[],2)>0.1,g(ok));
[b1,q] = Bin(r(:,5),nBins);
open q
limits = Accumulate(b1,values(:,1),'mode','min');
values = [r(:,5) d];
limits = Accumulate(b1,values(:,1),'mode','min');
limits = [Accumulate(b1,values(:,1),'mode','min') Accumulate(b1,values(:,1),'mode','max')];
b1 = Bin(r(:,4),20); b2 = Bin(r(:,5),20); b = [b1 b2];
limits = [min(x) max(x)]
limits = mmax(values(:,1));
limits = floor((-limits(1))/(limits(2)-limits(1))*nBins)+1;
limits = mmax(values(:,1));
(limits(2)-limits(1))*nBins)
(limits(2)-limits(1))*nBins
(-limits(1))/(limits(2)-limits(1))
(0)/(limits(2)-limits(1))
(limits(1))/(limits(2)-limits(1))
(limits(1)-limits(1))/(limits(2)-limits(1))
(limits(2)-limits(1))/(limits(2)-limits(1))
range(values(:,1))
range(values(:,1))/nBins;
range(values(:,1))/nBins
min(values(:,1)):range(values(:,1))/nBins:max(values(:,1));
limits = (min(values(:,1)):range(values(:,1))/nBins:max(values(:,1)))'; [limits = limits(1:end-1) limits(2:end)];
limits = (min(values(:,1)):range(values(:,1))/nBins:max(values(:,1)))'; limits = [limits(1:end-1) limits(2:end)];
midbins = mean(limits,2);
limits = (min(values(:,1)):range(values(:,1))/nBins:max(values(:,1)))'; limits = [limits(1:end-1) limits(2:end)];
midbins = mean(limits,2);
limits = (min(values(:,2)):range(values(:,2))/nBins:max(values(:,2)))'; limits = [limits(1:end-1) limits(2:end)];
midbins(:,2) = mean(limits,2);
clear all
load('day7.channelinfo.ripples.mat')
ripples = DetectSWR([rippleChannels.Ripple_Channel rippleChannels.Sharpwave_Channel rippleChannels.Noise_Channel],'basepath',basepath,'saveMat',true,'forceDetect',true,'check',true,'useSPW',false);
basepath = pwd;
ripples = DetectSWR([rippleChannels.Ripple_Channel rippleChannels.Sharpwave_Channel rippleChannels.Noise_Channel],'basepath',basepath,'saveMat',true,'forceDetect',true,'check',true,'useSPW',false);
% First, remove points that occur in noisy lfp periods
tl = (1:length(lfp))'/1250; % timestamps for lfp
t = featureTs/1250; % timestamps for candidate ripple events
bad = false(size(t));
try % optionally, remove periods of noisy LFP (large deflections)
[clean,bad,badIntervals] = CleanLFP([tl,double(lfp(:,2))],'thresholds',[10 Inf]);
bad = bad | InIntervals(t,badIntervals);
end
try % optionally, remove periods when your channels are diverging abnormally for a long time (a channel went dead for some time)
% if you need to use this, make sure your threshold if appropriate for your session. Plot the smoothed signal to see what seems
% reasonable to you!
smoothed = Smooth(abs(swDiff),1250*10);
noisyPeriods = ConsolidateIntervals(tl(FindInterval(smoothed>500)),'epsilon',1);
bad = bad | InIntervals(t,noisyPeriods);
end
try % optionally, remove ripples in which the sharp wave channel was positive
% Only do this if you have a sharp-wave channel what's sufficiently deep for the sharp wave to be negative
sw = interp1(tl,lfpLow(:,2),t);
bad = bad | sw>0;
end
idx1 = idx1 & ~bad; idx2 = idx2 & ~bad; % remove bad points from "idx1" and "idx2"
figure(1); % plot a clean figure without these points and check indivudual ripples
clf
plot(swDiffAll(idx1 | idx2),ripPowerAll(idx1 | idx2),'k.','markersize',1); hold on
scores = nan(size(ripPowerAll,1),1);
% figure(1); % plot a clean figure without these points and check indivudual ripples
% clf
% plot(swDiffAll(idx1 | idx2),ripPowerAll(idx1 | idx2),'k.','markersize',1); hold on
% scores = nan(size(ripPowerAll,1),1);
matrix = [swDiffAll ripPowerAll];
z = bsxfun(@rdivide,bsxfun(@minus,matrix,mean(matrix(~bad,:))),std(matrix(~bad,:)));
while true % click "Ctrl+C" to exit when you feel you're done
figure(1);
colors = Bright(1000);
[x,y,button] = ginput(1);
i = findmin(abs(x-swDiffAll)+abs(y-ripPowerAll));
t = featureTs/1250;
tl = (1:length(lfp))'/1250;
figure(2);
clf; interval = t(i) + [-1 1]*0.5; in = tl>interval(1) & tl<interval(end);
subplot(3,1,1);
plot(tl(in) - t(i),lfp(in,:)); title(num2str(t(i)));
PlotHVLines(Restrict(t,interval) - t(i),'v','k--');
legend('ripple channel','SW channel','noise channel');
subplot(3,1,2);
plot(tl(in) - t(i),double(swDiff(in)));
PlotHVLines(Restrict(t,interval) - t(i),'v','k--');
subplot(3,1,3);
plot(tl(in) - t(i),double(ripPower0(in)));
PlotHVLines(Restrict(t,interval) - t(i),'v','k--');
%     x = input( prompt )
commandwindow % select the command window
score = str2double(input('good(1)/bad(0)?','s'));
scores(i) = score; % save scores to optionally save
figure(1);
xlims = xlim; ylims = ylim; clf
% extrapolate the score of each ripple based on the score of the nearest scored ripple
scored = find(~isnan(scores));
distances = sqrt(bsxfun(@minus,z(:,1),z(scored,1)').^2 + bsxfun(@minus,z(:,2),z(~isnan(scores),2)').^2);
[~,nearestNeighbour] = min(distances,[],2);
estimated = scores(scored(nearestNeighbour));
scatter(matrix(:,1),matrix(:,2),1,estimated,'filled'); colormap(Bright); clim([0 1]);
xlim(xlims); ylim(ylims);
hold on;
scatter(swDiffAll(scored),ripPowerAll(scored),20,scores(scored)); scatter(swDiffAll(scored),ripPowerAll(scored),15,scores(scored));
clabel('Estimated ripple score');
xlabel('Sharp wave depth'); ylabel('Ripple power');
title({['Manual scoring for ' basepath],['called with channels: ' num2str(Channels) ' (ripple, sharp wave, and (optionally) noise)']});
if rem(i,10)==0
display('saving...');
save(fullfile(basepath,'DetectSWR_manual_scoring.mat'),'swDiffAll','ripPowerAll','bad','scores');
end
end
% First, remove points that occur in noisy lfp periods
tl = (1:length(lfp))'/1250; % timestamps for lfp
t = featureTs/1250; % timestamps for candidate ripple events
bad = false(size(t));
try % optionally, remove periods of noisy LFP (large deflections)
[clean,bad,badIntervals] = CleanLFP([tl,double(lfp(:,2))],'thresholds',[10 Inf]);
bad = bad | InIntervals(t,badIntervals);
end
try % optionally, remove periods when your channels are diverging abnormally for a long time (a channel went dead for some time)
% if you need to use this, make sure your threshold if appropriate for your session. Plot the smoothed signal to see what seems
% reasonable to you!
smoothed = Smooth(abs(swDiff),1250*10);
noisyPeriods = ConsolidateIntervals(tl(FindInterval(smoothed>500)),'epsilon',1);
bad = bad | InIntervals(t,noisyPeriods);
end
try % optionally, remove ripples in which the sharp wave channel was positive
% Only do this if you have a sharp-wave channel what's sufficiently deep for the sharp wave to be negative
sw = interp1(tl,lfpLow(:,2),t);
%                 bad = bad | sw>0;
end
idx1 = idx1 & ~bad; idx2 = idx2 & ~bad; % remove bad points from "idx1" and "idx2"
figure(1); % plot a clean figure without these points and check indivudual ripples
clf
plot(swDiffAll(idx1 | idx2),ripPowerAll(idx1 | idx2),'k.','markersize',1); hold on
scores = nan(size(ripPowerAll,1),1);
% figure(1); % plot a clean figure without these points and check indivudual ripples
% clf
% plot(swDiffAll(idx1 | idx2),ripPowerAll(idx1 | idx2),'k.','markersize',1); hold on
% scores = nan(size(ripPowerAll,1),1);
matrix = [swDiffAll ripPowerAll];
z = bsxfun(@rdivide,bsxfun(@minus,matrix,mean(matrix(~bad,:))),std(matrix(~bad,:)));
while true % click "Ctrl+C" to exit when you feel you're done
figure(1);
colors = Bright(1000);
[x,y,button] = ginput(1);
i = findmin(abs(x-swDiffAll)+abs(y-ripPowerAll));
t = featureTs/1250;
tl = (1:length(lfp))'/1250;
figure(2);
clf; interval = t(i) + [-1 1]*0.5; in = tl>interval(1) & tl<interval(end);
subplot(3,1,1);
plot(tl(in) - t(i),lfp(in,:)); title(num2str(t(i)));
PlotHVLines(Restrict(t,interval) - t(i),'v','k--');
legend('ripple channel','SW channel','noise channel');
subplot(3,1,2);
plot(tl(in) - t(i),double(swDiff(in)));
PlotHVLines(Restrict(t,interval) - t(i),'v','k--');
subplot(3,1,3);
plot(tl(in) - t(i),double(ripPower0(in)));
PlotHVLines(Restrict(t,interval) - t(i),'v','k--');
%     x = input( prompt )
commandwindow % select the command window
score = str2double(input('good(1)/bad(0)?','s'));
scores(i) = score; % save scores to optionally save
figure(1);
xlims = xlim; ylims = ylim; clf
% extrapolate the score of each ripple based on the score of the nearest scored ripple
scored = find(~isnan(scores));
distances = sqrt(bsxfun(@minus,z(:,1),z(scored,1)').^2 + bsxfun(@minus,z(:,2),z(~isnan(scores),2)').^2);
[~,nearestNeighbour] = min(distances,[],2);
estimated = scores(scored(nearestNeighbour));
scatter(matrix(:,1),matrix(:,2),1,estimated,'filled'); colormap(Bright); clim([0 1]);
xlim(xlims); ylim(ylims);
hold on;
scatter(swDiffAll(scored),ripPowerAll(scored),20,scores(scored)); scatter(swDiffAll(scored),ripPowerAll(scored),15,scores(scored));
clabel('Estimated ripple score');
xlabel('Sharp wave depth'); ylabel('Ripple power');
title({['Manual scoring for ' basepath],['called with channels: ' num2str(Channels) ' (ripple, sharp wave, and (optionally) noise)']});
if rem(i,10)==0
display('saving...');
save(fullfile(basepath,'DetectSWR_manual_scoring.mat'),'swDiffAll','ripPowerAll','bad','scores');
end
end
bad = false(size(t));
bad = bad | InIntervals(t,badIntervals);
smoothed = Smooth(abs(swDiff),1250*10);
noisyPeriods = ConsolidateIntervals(tl(FindInterval(smoothed>500)),'epsilon',1);
bad = bad | InIntervals(t,noisyPeriods);
Portion(bad)
Portion(InIntervals(t,noisyPeriods))
idx1 = idx1 & ~bad; idx2 = idx2 & ~bad; % remove bad points from "idx1" and "idx2"
figure(1); % plot a clean figure without these points and check indivudual ripples
clf
plot(swDiffAll(idx1 | idx2),ripPowerAll(idx1 | idx2),'k.','markersize',1); hold on
scores = nan(size(ripPowerAll,1),1);
% figure(1); % plot a clean figure without these points and check indivudual ripples
% clf
% plot(swDiffAll(idx1 | idx2),ripPowerAll(idx1 | idx2),'k.','markersize',1); hold on
% scores = nan(size(ripPowerAll,1),1);
matrix = [swDiffAll ripPowerAll];
z = bsxfun(@rdivide,bsxfun(@minus,matrix,mean(matrix(~bad,:))),std(matrix(~bad,:)));
edit nothing
nothing
1
0
1
0
1
0
figure; plot(t,double(lfp(:,2)));
plot(tl,double(lfp(:,2)));
PlotIntervals(badIntervals,'color','k');
PlotHVLines(8578.7,'v');
t(i)
addComma(t(i))
bad(i)
idx1(i)
idx2(i)
bad(i)
t(i)
xlim([-500 500]+t(i))
nothing
0
1
0.6
1
0.2
0
1
0
1
0
1
0
0.2
0
1
0.8
0.2
1
0
1
save(fullfile(basepath,'DetectSWR_manual_scoring.mat'),'swDiffAll','ripPowerAll','bad','scores');
dbquit
ripples = DetectSWR([rippleChannels.Ripple_Channel rippleChannels.Sharpwave_Channel],'basepath',basepath,'saveMat',true,'forceDetect',true,'check',true,'useSPW',false);
% First, remove points that occur in noisy lfp periods
tl = (1:length(lfp))'/1250; % timestamps for lfp
t = featureTs/1250; % timestamps for candidate ripple events
bad = false(size(t));
try % optionally, remove periods of noisy LFP (large deflections)
[clean,~,badIntervals] = CleanLFP([tl,double(lfp(:,2))],'thresholds',[10 Inf]);
bad = bad | InIntervals(t,badIntervals);
end
try % optionally, remove periods when your channels are diverging abnormally for a long time (a channel went dead for some time)
% if you need to use this, make sure your threshold if appropriate for your session. Plot the smoothed signal to see what seems
% reasonable to you!
smoothed = Smooth(abs(swDiff),1250*10);
noisyPeriods = ConsolidateIntervals(tl(FindInterval(smoothed>500)),'epsilon',1);
bad = bad | InIntervals(t,noisyPeriods);
end
try % optionally, remove ripples in which the sharp wave channel was positive
% Only do this if you have a sharp-wave channel what's sufficiently deep for the sharp wave to be negative
sw = interp1(tl,lfpLow(:,2),t);
%                 bad = bad | sw>0;
end
idx1 = idx1 & ~bad; idx2 = idx2 & ~bad; % remove bad points from "idx1" and "idx2"
figure(1); % plot a clean figure without these points and check indivudual ripples
clf
plot(swDiffAll(idx1 | idx2),ripPowerAll(idx1 | idx2),'k.','markersize',1); hold on
scores = nan(size(ripPowerAll,1),1);
% figure(1); % plot a clean figure without these points and check indivudual ripples
% clf
% plot(swDiffAll(idx1 | idx2),ripPowerAll(idx1 | idx2),'k.','markersize',1); hold on
% scores = nan(size(ripPowerAll,1),1);
matrix = [swDiffAll ripPowerAll];
z = bsxfun(@rdivide,bsxfun(@minus,matrix,mean(matrix(~bad,:))),std(matrix(~bad,:)));
Portion(bad)
nothing
0
nothing
0
Portion(bad)
x
y
i
t(i)
i = findmin(abs(x-swDiffAll.*(1-bad))+abs(y-ripPowerAll.*(1-bad)))
swDiffAll(i)
nothing
0
1
0
1
0
1
0
0.1
0
1
0
1
Browse('all');
nothing
1
0.1
0.8
0
0.8
1
0
1
1
1
0
[w,wt,wf] = WaveletSpectrogram(Restrict([tl double(lfp(:,1))],interval),'range',[50 625]);
figure; PlotColorMap(w,'x',wt,'y',wf);
nothing
1
0.9
1
0.2
0.1
0
1
scores(i) = nan;
nothing
0
1
0
1
0
1
0
1
0
peaks = FindLocalMaxima([tl(:,1) ripPower0]);
troughs = FindLocalMinima([tl(:,1) ripPower0]);
candidate = peak(FindClosest(peak,t));
candidate = peaks(FindClosest(peaks,t));
open candidate
PlotHVLines(candidate(i));
PlotHVLines(candidate(i)-t(i));
candidate = [troughs(FindClosest(candidate,t),'lower') candidate troughs(FindClosest(candidate,t),'higher')];
candidate = [troughs(FindClosest(troughs,candidate),'lower') candidate troughs(FindClosest(troughs,candidate),'higher')];
troughs(FindClosest(troughs,candidate),'lower');
candidate = [troughs(FindClosest(troughs,candidate,'lower')) candidate troughs(FindClosest(troughs,candidate,'higher'))];
PlotHVLines(candidate(i,:)-t(i));
dd = candidate(:,3)-candidate(:,1);
figure; plot(dd,score,'.');
clf
plot(dd,scores,'.');
dd(i)
nancorr(dd,scores)
min(dd(scores==1))
find(dd==ans & scores==1)
i=139998
t = featureTs/1250;
tl = (1:length(lfp))'/1250;
figure(2);
clf; interval = t(i) + [-1 1]*0.5; in = tl>interval(1) & tl<interval(end);
subplot(3,1,1);
plot(tl(in) - t(i),lfp(in,:)); title(num2str(t(i)));
PlotHVLines(Restrict(t,interval) - t(i),'v','k--');
legend('ripple channel','SW channel','noise channel');
subplot(3,1,2);
plot(tl(in) - t(i),double(swDiff(in)));
PlotHVLines(Restrict(t,interval) - t(i),'v','k--');
subplot(3,1,3);
plot(tl(in) - t(i),double(ripPower0(in)));
PlotHVLines(Restrict(t,interval) - t(i),'v','k--');
i = findmin(abs(x-swDiffAll.*(1-bad))+abs(y-ripPowerAll.*(1-bad)));
t = featureTs/1250;
tl = (1:length(lfp))'/1250;
figure(2);
clf; interval = t(i) + [-1 1]*0.5; in = tl>interval(1) & tl<interval(end);
subplot(3,1,1);
plot(tl(in) - t(i),lfp(in,:)); title(num2str(t(i)));
PlotHVLines(Restrict(t,interval) - t(i),'v','k--');
legend('ripple channel','SW channel','noise channel');
subplot(3,1,2);
plot(tl(in) - t(i),double(swDiff(in)));
PlotHVLines(Restrict(t,interval) - t(i),'v','k--');
subplot(3,1,3);
plot(tl(in) - t(i),double(ripPower0(in)));
PlotHVLines(Restrict(t,interval) - t(i),'v','k--');
starts = strfind(ripPower0'>40,[0 1]);
stops = strfind(ripPower0'>40,[1 0]);
starts = tl(strfind(ripPower0'>40,[0 1]));
stops = tl(strfind(ripPower0'>40,[1 0]));
candidate = peaks(FindClosest(peaks,t));
candidate = [starts(FindClosest(starts,candidate,'lower')) candidate stops(FindClosest(stops,candidate,'higher'))];
dd2 = candidate(:,3)-candidate(:,1);
nancorr(dd2,scores)
plot(dd2,scores,'.');
dd2(i)
min(dd2(scores==1))
find(dd2==ans & scores==1)
i=132139
t = featureTs/1250;
tl = (1:length(lfp))'/1250;
figure(2);
clf; interval = t(i) + [-1 1]*0.5; in = tl>interval(1) & tl<interval(end);
subplot(3,1,1);
plot(tl(in) - t(i),lfp(in,:)); title(num2str(t(i)));
PlotHVLines(Restrict(t,interval) - t(i),'v','k--');
legend('ripple channel','SW channel','noise channel');
subplot(3,1,2);
plot(tl(in) - t(i),double(swDiff(in)));
PlotHVLines(Restrict(t,interval) - t(i),'v','k--');
subplot(3,1,3);
plot(tl(in) - t(i),double(ripPower0(in)));
PlotHVLines(Restrict(t,interval) - t(i),'v','k--');
candidate(i,:)
candidate(i,:)-candidate(i,2)
nothing
1
0
1
0
1
0.2
1
nothing
1
nothing
1
0.8
1
nothing
1
0
1
while true % click "Ctrl+C" to exit when you feel you're done
figure(1);
colors = Bright(1000);
[x,y,button] = ginput(1);
i = findmin(abs(x-swDiffAll)+abs(y-ripPowerAll));
t = featureTs/1250;
tl = (1:length(lfp))'/1250;
figure(2);
clf; interval = t(i) + [-1 1]*0.5; in = tl>interval(1) & tl<interval(end);
subplot(3,1,1);
plot(tl(in) - t(i),lfp(in,:)); title(num2str(t(i)));
PlotHVLines(Restrict(t,interval) - t(i),'v','k--');
legend('ripple channel','SW channel','noise channel');
subplot(3,1,2);
plot(tl(in) - t(i),double(swDiff(in)));
PlotHVLines(Restrict(t,interval) - t(i),'v','k--');
subplot(3,1,3);
plot(tl(in) - t(i),double(ripPower0(in)));
PlotHVLines(Restrict(t,interval) - t(i),'v','k--');
%     x = input( prompt )
commandwindow % select the command window
score = str2double(input('good(1)/bad(0)?','s'));
scores(i) = score; % save scores to optionally save
figure(1);
xlims = xlim; ylims = ylim; clf
% extrapolate the score of each ripple based on the score of the nearest scored ripple
scored = find(~isnan(scores));
distances = sqrt(bsxfun(@minus,z(:,1),z(scored,1)').^2 + bsxfun(@minus,z(:,2),z(~isnan(scores),2)').^2);
[~,nearestNeighbour] = min(distances,[],2);
estimated = scores(scored(nearestNeighbour));
scatter(matrix(:,1),matrix(:,2),1,estimated,'filled'); colormap(Bright); clim([0 1]);
xlim(xlims); ylim(ylims);
hold on;
scatter(swDiffAll(scored),ripPowerAll(scored),20,scores(scored)); scatter(swDiffAll(scored),ripPowerAll(scored),15,scores(scored));
clabel('Estimated ripple score');
xlabel('Sharp wave depth'); ylabel('Ripple power');
title({['Manual scoring for ' basepath],['called with channels: ' num2str(Channels) ' (ripple, sharp wave, and (optionally) noise)']});
if rem(i,10)==0
display('saving...');
save(fullfile(basepath,'DetectSWR_manual_scoring.mat'),'swDiffAll','ripPowerAll','bad','scores');
end
end
1
nothing
0
1
0
1
0
1
0
1
0
1
0
0.7
1
0.9
1
0.8
1
31
save(fullfile(basepath,'DetectSWR_manual_scoring.mat'),'swDiffAll','ripPowerAll','bad','scores');
% (I tend to label an ugly ripple that is still clearly a ripple with a score of 0.2)
scored = find(~isnan(scores));
distances = sqrt(bsxfun(@minus,z(:,1),z(scored,1)').^2 + bsxfun(@minus,z(:,2),z(~isnan(scores),2)').^2);
[~,nearestNeighbour] = min(distances,[],2);
estimated = scores(scored(nearestNeighbour));
selected = estimated>0.1;
idx1 = selected & ~bad; % final ripples
idx2 = ~selected & ~bad; % non-ripples (yet free from noise as well)
display('saving...');
save(fullfile(basepath,'DetectSWR_manual_scoring.mat'),'swDiffAll','ripPowerAll','idx1','idx2','bad','selected','scores');
saveas(figure(1),fullfile(basepath,'DetectSWR_manual_scoring.fig'));
dbcont
r = [ripples.timestamps(:,1) ripples.peaks ripples.timestamps(:,2) ripples.RipMax ripples.SwMax];
load('day7.MergePoints.events.mat')
sleep = ConsolidateIntervals(MergePoints.timestamps(cellfun(@(x) any(strfind(lower(x),'sleep')),MergePoints.foldernames),:));
pre = InIntervals(r,sleep(1,:));
post = InIntervals(r,sleep(end,:));
g = [-pre+post]+2;
figure;  anovabar(diff(r(:,[1 3]),[],2)>0.1,g);
size(r)./r(end,1)
load('day7.cell_metrics.cellinfo.mat')
spikesCell = cell_metrics.spikes.times';
nNeurons = length(cell_metrics.spikes.times);
regions = cell_metrics.brainRegion;
regionNames = unique(regions)
regionCell = zeros(length(regions),1);
for i=1:length(regionNames)
regionCell(strcmp(regions,regionNames{i})) = i;
end
regionCell
spikesCell = cell_metrics.spikes.times';
spikes = sortrows(Group(spikesCell{:}));
[q,qt] = PETH(spikes(:,1),r(:,1));
semplot(qt,q);
clf
PlotColorMap(sortby(q,sum(q(:,qt>0 & qt<0.2)),'x',qt);
PlotColorMap(sortby(q,sum(q(:,qt>0 & qt<0.2,2)),'x',qt));
PlotColorMap(sortby(q,sum(q(:,qt>0 & qt<0.2,2))),'x',qt);
PlotColorMap(sortby(q,sum(q(:,qt>0 & qt<0.2),2)),'x',qt);
PlotColorMap(sortby(q,sum(q(:,qt>0 & qt<0.1),2)),'x',qt);
PlotColorMap(sortby(q,sum(q(:,qt>0 & qt<0.05),2)),'x',qt);
PlotColorMap(sortby(zscore(q,[],2),sum(q(:,qt>0 & qt<0.05),2)),'x',qt);
PlotColorMap(Shrink(sortby(zscore(q,[],2),sum(q(:,qt>0 & qt<0.05),2)),100,1),'x',qt);
hist(r(:,5),100)
ok = r(:<5)<2.5; anovabar(diff(r(ok,[1 3]),[],2)>0.1,g(ok));
ok = r(:,5)<2.5; anovabar(diff(r(ok,[1 3]),[],2)>0.1,g(ok));
ok = r(:,5)>2.5; anovabar(diff(r(ok,[1 3]),[],2)>0.1,g(ok));
ok = r(:,5)>2.5; anovabar(diff(r(ok,[1 3]),[],2),g(ok));
[h,ht] = Dist(100,[diff(r(ok,[1 3]),[],2),g(ok)+2],'grouped');
plot(ht*1000,Smooth(cumsum(h),[0 0]),'linewidth',2); set(gca,'ytick',[]); xlabel('ripple duration (s)');
clf
plot(ht*1000,Smooth(cumsum(h),[0 0]),'linewidth',2); set(gca,'ytick',[]); xlabel('ripple duration (s)');
open g
size(h)
mmax(g)
[h,ht] = Dist(100,[diff(r(ok,[1 3]),[],2),g(ok)],'grouped');
plot(ht*1000,Smooth(cumsum(h),[0 0]),'linewidth',2); set(gca,'ytick',[]); xlabel('ripple duration (s)');
plot(ht*1000,Smooth((h),[0 0]),'linewidth',2); set(gca,'ytick',[]); xlabel('ripple duration (s)');
ok = r(:,5)<2.5; anovabar(diff(r(ok,[1 3]),[],2),g(ok));
[h,ht] = Dist(100,[diff(r(ok,[1 3]),[],2),g(ok)],'grouped');
plot(ht*1000,Smooth((h),[0 0]),'linewidth',2); set(gca,'ytick',[]); xlabel('ripple duration (s)');
clear all
baepath = 'N:\OJRproject\OJR42\day11';
basepath = 'N:\OJRproject\OJR42\day11';
cd(basepath)
5/0.9+
5/0.9
load('Copy_of_day11.ripples.events.mat')
3/0.028
7/0.054
4/0.035
ripples = DetectSWR([43 21 114]+1,'basepath',basepath,'saveMat',true,'forceDetect',true,'check',true,'useSPW',false);
tl = (1:length(lfp))'/1250; % timestamps for lfp
t = featureTs/1250; % timestamps for candidate ripple events
bad = false(size(t));
[clean,~,badIntervals] = CleanLFP([tl,double(lfp(:,2))],'thresholds',[10 Inf]);
bad = bad | InIntervals(t,badIntervals);
smoothed = Smooth(abs(swDiff),1250*10);
noisyPeriods = ConsolidateIntervals(tl(FindInterval(smoothed>500)),'epsilon',1);
bad = bad | InIntervals(t,noisyPeriods);
Portion(bad)
sw = interp1(tl,lfpLow(:,2),t);
Portion(sw>0)
figure; plot(tl,lfp(:,2));
PlotIntervals(badIntervals,'color','k');
[clean,~,badIntervals] = CleanLFP([tl,double(lfp(:,2))],'thresholds',[8 Inf]);
PlotIntervals(badIntervals,'color','k');
[clean,~,badIntervals] = CleanLFP([tl,double(lfp(:,2))],'thresholds',[5 Inf]);
PlotIntervals(badIntervals,'color','k');
bad = bad | InIntervals(t,badIntervals);
Portion(InIntervals(t,noisyPeriods)_
Portion(InIntervals(t,noisyPeriods))
Portion(InIntervals(t,badIntervals))
PlotIntervals(noisyPeriods,'color','r');
clf
plot(tl,smoothed);
PlotIntervals(badIntervals,'color','k');
bad = false(size(t));
bad = bad | InIntervals(t,badIntervals);
idx1 = idx1 & ~bad; idx2 = idx2 & ~bad; % remove bad points from "idx1" and "idx2"
figure(1); % plot a clean figure without these points and check indivudual ripples
clf
plot(swDiffAll(idx1 | idx2),ripPowerAll(idx1 | idx2),'k.','markersize',1); hold on
scores = nan(size(ripPowerAll,1),1);
% figure(1); % plot a clean figure without these points and check indivudual ripples
% clf
% plot(swDiffAll(idx1 | idx2),ripPowerAll(idx1 | idx2),'k.','markersize',1); hold on
% scores = nan(size(ripPowerAll,1),1);
matrix = [swDiffAll ripPowerAll];
z = bsxfun(@rdivide,bsxfun(@minus,matrix,mean(matrix(~bad,:))),std(matrix(~bad,:)));
open nothing
nothing
clf
plot(swDiffAll(idx1 | idx2),ripPowerAll(idx1 | idx2),'k.','markersize',1); hold on
nothing
0
1
0.5
0
i
3/0.01901
[w,wt,wf] = WaveletSpectrogram(tl(in) - t(i),lfp(in,1),'range',[50 625]);
[w,wt,wf] = WaveletSpectrogram([tl(in) - t(i),lfp(in,1)],'range',[50 625]);
[w,wt,wf] = WaveletSpectrogram([tl(in) - t(i)+5,lfp(in,1)],'range',[50 625]);
[w,wt,wf] = WaveletSpectrogram([tl(in) - t(i)+5,double(lfp(in,1))],'range',[50 625]);
figure; PlotColorMap(w,'x',wt,'y',wf);
help FilterLFP
filtered = Filter0([tl(in) - t(i)+5,double(lfp(in,1))],[100 250]);
figure; PlotXY(filtered);
filtered = Filter0([tl(in) - t(i),double(lfp(in,1))],[100 250]);
PlotXY(filtered);
[~,a] = Phase(filtered);
hold all; PlotXY(a);
PlotXY([a(:,1), sqrt(a(:,2))]);
hold all
sqrt
PlotXY([a(:,1), (a(:,2))/2]);
hold all
PlotXY([a(:,1), (a(:,2))/2]);
powerWin
ripBP = [100 250];
hRip1      = makegausslpfir( ripBP( 1 ), SR, 6 );
hRip2      = makegausslpfir( ripBP( 2 ), SR, 6 );
ripPower00 = ripPower0;
% ripple channel
rip        = firfilt( lfp(:,1), hRip2 );    % highpass filter
eegLo      = firfilt( rip, hRip1 );     % lowpass filter
rip        = rip - eegLo;               % difference of gaussians
ripWindow  = pi / mean( ripBP );
powerWin   = makegausslpfir( 1 / ripWindow, SR, 6 );
rip        = abs(rip);
signal  = firfilt( rip, powerWin );
rip        = firfilt( lfp(:,3), hRip2 );    % highpass filter
eegLo      = firfilt( rip, hRip1 );     % lowpass filter
rip        = rip - eegLo;               % difference of gaussians
ripWindow  = pi / mean( ripBP );
powerWin   = makegausslpfir( 1 / ripWindow, SR, 6 );
rip        = abs(rip);
noise  = firfilt( rip, powerWin );
ripPower0 = signal-noise;
figure(2);
clf; interval = t(i) + [-1 1]*0.5; in = tl>interval(1) & tl<interval(end);
subplot(3,1,1);
plot(tl(in) - t(i),lfp(in,:)); title(num2str(t(i)));
PlotHVLines(Restrict(t,interval) - t(i),'v','k--');
legend('ripple channel','SW channel','noise channel');
subplot(3,1,2);
plot(tl(in) - t(i),double(swDiff(in)));
PlotHVLines(Restrict(t,interval) - t(i),'v','k--');
subplot(3,1,3);
plot(tl(in) - t(i),double(ripPower0(in)));
PlotHVLines(Restrict(t,interval) - t(i),'v','k--');
clims(clim/2);
PlotHVLines([100 250],'h','w--');
clims(clim*0.5);
clims(clim*1.5);
filteredHigh = Filter0([tl(in) - t(i),double(lfp(in,1))],[300 600]);
[~,ah] = Phase(filteredHigh);
PlotXY(a);
clf
PlotXY(a);
hold all
PlotXY(ah);
detrended = Detrend(lfp);
detrended = Detrend([t double(lfp(:,1))]);
detrended = Detrend([tl double(lfp(:,1))]);
filteredHigh = Filter0([tl(in) - t(i),double(lfp(in,1))],[300 600]);
filtered = Filter0([tl(in) - t(i),double(lfp(in,1))],[100 250]);
[~,ah] = Phase(filteredHigh);
[~,a] = Phase(filtered);
clf
PlotXY(a);
hold all
PlotXY(ah);
filteredHigh = Filter0(detrended(in,:),[300 600]);
filtered = Filter0(detrended(in,:),[100 250]);
[~,ah] = Phase(filteredHigh);
[~,a] = Phase(filtered);
clf
PlotXY(a);
hold all
PlotXY(ah);
zs = zscore([a(:,2) ah(:,2)]);
plot(zs(in,:))
size(zs)
clf
plot(zs)
filteredHigh = Filter0(detrended,[300 600]);
filtered = Filter0(detrended,[100 250]);
[~,ah] = Phase(filteredHigh);
[~,a] = Phase(filtered);
zs = zscore([a(:,2) ah(:,2)]);
plot(zs(in,:))
5
HalfWinSize
indices = bsxfun(@plus,featureTs,-HalfWinSize:HalfWinSize);
open indices
wrongRipplePowerAll = ripplePowerAll;
wrongRipPowerAll = ripPowerAll;
signal = z(:,2)-z(:,1);
signal = zs(:,2)-zs(:,1);
in = InIntervals(tl,interval);
plot([a(in,2) ah(in,2)]);
signal = a(:,2)-ah(:,1);
zs = zscore(ah(:,2));
noisy = max(zs(indices),[],2);
open noisy
clf
hist(noisy,100);
Portion(noisy>2)
noisy = max(zs(indices),[],2)>2;
noisy = max(zs(indices),[],2)>p2z(0.05);
Portion(noisy>2)
Portion(noisy)
ripPower0 = max(signal(indices),[],2);
anovabar(ripPower0,noisy)
Portion(bad)
bad = bad | InIntervals(t,badIntervals) | noisy;
i
clf; interval = t(i) + [-1 1]*0.5; in = tl>interval(1) & tl<interval(end);
subplot(3,1,1);
plot(tl(in) - t(i),lfp(in,:)); title(num2str(t(i)));
PlotHVLines(Restrict(t,interval) - t(i),'v','k--');
legend('ripple channel','SW channel','noise channel');
subplot(3,1,2);
plot(tl(in) - t(i),double(swDiff(in)));
PlotHVLines(Restrict(t,interval) - t(i),'v','k--');
subplot(3,1,3);
plot(tl(in) - t(i),double(ripPower0(in)));
PlotHVLines(Restrict(t,interval) - t(i),'v','k--');
clf; interval = t(i) + [-1 1]*0.5; in = tl>interval(1) & tl<interval(end);
subplot(3,1,1);
plot(tl(in) - t(i),lfp(in,:)); title(num2str(t(i)));
PlotHVLines(Restrict(t,interval) - t(i),'v','k--');
legend('ripple channel','SW channel','noise channel');
subplot(3,1,2);
plot(tl(in) - t(i),double(swDiff(in)));
PlotHVLines(Restrict(t,interval) - t(i),'v','k--');
subplot(3,1,3);
plot(tl(in) - t(i),double(ripPower0(in)));
size(ripPower0)
size(tl)
ripPowerAll = max(signal(indices),[],2);
ripPower0 = signal-noise;
figure(2);
clf; interval = t(i) + [-1 1]*0.5; in = tl>interval(1) & tl<interval(end);
subplot(3,1,1);
plot(tl(in) - t(i),lfp(in,:)); title(num2str(t(i)));
PlotHVLines(Restrict(t,interval) - t(i),'v','k--');
legend('ripple channel','SW channel','noise channel');
subplot(3,1,2);
plot(tl(in) - t(i),double(swDiff(in)));
PlotHVLines(Restrict(t,interval) - t(i),'v','k--');
subplot(3,1,3);
plot(tl(in) - t(i),double(ripPower0(in)));
PlotHVLines(Restrict(t,interval) - t(i),'v','k--');
nanmean(ripPower0)
nanmean(signal)
signal = a(:,2)-ah(:,2);
nanmean(signal)
ripPower0 = signal-noise;
ripPowerAll = max(signal(indices),[],2);
figure(2);
clf; interval = t(i) + [-1 1]*0.5; in = tl>interval(1) & tl<interval(end);
subplot(3,1,1);
plot(tl(in) - t(i),lfp(in,:)); title(num2str(t(i)));
PlotHVLines(Restrict(t,interval) - t(i),'v','k--');
legend('ripple channel','SW channel','noise channel');
subplot(3,1,2);
plot(tl(in) - t(i),double(swDiff(in)));
PlotHVLines(Restrict(t,interval) - t(i),'v','k--');
subplot(3,1,3);
plot(tl(in) - t(i),double(ripPower0(in)));
PlotHVLines(Restrict(t,interval) - t(i),'v','k--');
figure(1); % plot a clean figure without these points and check indivudual ripples
clf
plot(swDiffAll(idx1 | idx2),ripPowerAll(idx1 | idx2),'k.','markersize',1); hold on
scores = nan(size(ripPowerAll,1),1);
% figure(1); % plot a clean figure without these points and check indivudual ripples
% clf
% plot(swDiffAll(idx1 | idx2),ripPowerAll(idx1 | idx2),'k.','markersize',1); hold on
% scores = nan(size(ripPowerAll,1),1);
matrix = [swDiffAll ripPowerAll];
z = bsxfun(@rdivide,bsxfun(@minus,matrix,mean(matrix(~bad,:))),std(matrix(~bad,:)));
nothing
1
0
1
0
1
0
0.9
0
0.2
swDiff = -swDiff;
swDiffAll= -swDiffAll;
ripples = DetectSWR([43 21 114]+1,'basepath',basepath,'saveMat',true,'forceDetect',true,'check',true,'useSPW',false);
tl = (1:length(lfp))'/1250; % timestamps for lfp
detrended = Detrend([tl double(lfp(:,1))]);
filtered = Filter0(detrended,[100 250]);
filteredHigh = Filter0(detrended,[300 600]);
[~,a] = Phase(filtered);
[~,ah] = Phase(filteredHigh);
signal = a(:,2)-ah(:,2);
ripPower0 = signal;
indices = bsxfun(@plus,featureTs,-HalfWinSize:HalfWinSize);
ripPowerAll = max(signal(indices),[],2);
zs = zscore(ah(:,2));
noisy = max(zs(indices),[],2)>p2z(0.05);
t = featureTs/1250; % timestamps for candidate ripple events
bad = false(size(t));
[clean,~,badIntervals] = CleanLFP([tl,double(lfp(:,2))],'thresholds',[10 Inf]);
bad = bad | InIntervals(t,badIntervals);
Portion(bad)
Portion(noisy)
bad = bad | noisy;
idx1 = idx1 & ~bad; idx2 = idx2 & ~bad; % remove bad points from "idx1"
Portion(idx1|idx2)
figure(1); % plot a clean figure without these points and check indivudual ripples
clf
plot(swDiffAll(idx1 | idx2),ripPowerAll(idx1 | idx2),'k.','markersize',1); hold on
scores = nan(size(ripPowerAll,1),1);
% figure(1); % plot a clean figure without these points and check indivudual ripples
% clf
% plot(swDiffAll(idx1 | idx2),ripPowerAll(idx1 | idx2),'k.','markersize',1); hold on
% scores = nan(size(ripPowerAll,1),1);
matrix = [swDiffAll ripPowerAll];
z = bsxfun(@rdivide,bsxfun(@minus,matrix,mean(matrix(~bad,:))),std(matrix(~bad,:)));
nothing
0
1
0.8
0.1
1
0.2
save(fullfile(basepath,'DetectSWR_manual_scoring.mat'),'swDiffAll','ripPowerAll','bad','scores');
save(fullfile(basepath,'DetectSWR_manual_scoring.mat'),'t','swDiffAll','ripPowerAll','bad','scores');
ripples = FindRipples(basepath,43+1);
ripples = FindRipples(basepath,43+1,'noise',114+1);
ripples = FindRipples(basepath,43+1,'noise',114+1,'show','on','passband',[110 250]);
open BatchFindRipples.m
open preprocessSession.m
SleepScoreMaster(basepath,'noPrompts',true,'ignoretime',pulses.intsPeriods); % try to sleep score
SleepScoreMaster(basepath,'noPrompts',true);
thetaEpochs(basepath);
load('day11.EMGFromLFP.LFP.mat')
load('day11.SleepState.states.mat')
SleepScoreMaster(basepath,'noPrompts',true,'overwrite',true,'SWChannels',1+[114 98 91 110],'ThetaChannels',1+[48 20 23 4],'rejectChannels',1+[49 54 52 55 62 37 63 38 60 40 61 41 32 36 33 34 29 35 7 12 13 14 15 8 11 84 94 69 95 70 92 72 93 73 64 68 65 66 78 106 67 71 121 125]);
load('day11.EMGFromLFP.LFP.mat')
load('day11.SleepState.states.mat')
load('day11.SleepStateEpisodes.states.mat')
ripples = BatchFindRipples(pwd,43,[],SleepState.ints.WAKEstate);
ripples = BatchFindRipples(pwd,[21 43 17 47],[],SleepState.ints.WAKEstate);
45
r = ripples;
SaveCustomEvents('ripples.rip.evt',r(:,1:3),{'ripple start','ripple peak','ripple stop'});
lfp = GetAyaLFP(114);
display(['loaded! @' num2str(toc)]);
clean = CleanLFP(lfp);
deltas0 = FindDeltaPeaks(clean);
close all
PETH(r(:,2),deltas(:,2));
PETH(r(:,2),deltas0(:,2));
[h,ht] = PETH(r(:,2),deltas0(:,2));
PlotColorMap(Shrink(h,100,1));
PlotColorMap(Shrink(sortby(h,deltas0(:,4)),100,1));
PlotColorMap(Shrink(sortby(h,deltas0(:,5)),100,1));
PlotColorMap(Shrink(sortby(h,deltas0(:,5)),1000,1));
PlotColorMap(Shrink(sortby(h,deltas0(:,4)),1000,1));
PlotColorMap(Shrink(sortby(h,deltas0(:,4)),1000,4));
PlotColorMap(Shrink(sortby(h,deltas0(:,4)),100,4));
PlotColorMap(Shrink(sortby(h,deltas0(:,5)),100,4));
PlotColorMap(Shrink(sortby(h,deltas0(:,4)-deltas0(:,5)),100,4));
deltas0(:,7) = (deltas0(:,5)-deltas0(:,6))./(deltas0(:,3)-deltas0(:,2));
PlotColorMap(Shrink(sortby(h,deltas0(:,6)),100,4));
PlotColorMap(Shrink(sortby(h,deltas0(:,7)),100,4));
PlotColorMap(Shrink(sortby(h,deltas0(:,7)),1,4));
PlotColorMap(Shrink(sortby(h,deltas0(:,7)),1000,4));
PlotColorMap(Shrink(sortby(deltas0(:,7),deltas0(:,7)),1000,1));
PlotColorMap(Shrink(sortby(h,deltas0(:,7)),1000,4));
PlotColorMap(Shrink(sortby(h,deltas0(:,7)),500,4));
figure; hist(deltas0(:,7),100);
semplot(ht,h(deltas(:,7)<20,:))
semplot(ht,h(deltas0(:,7)<20,:))
b = Bin(deltas0(:,7));
b = Bin(deltas0(:,7),5);
semplot(ht,h(b==1,:))
clf
semplot(ht,h(b==2,:))
quantile(deltas0(:,7),0.2)
b = Bin(deltas0(:,7),10);
quantile(deltas0(:,7),0.1)
quantile(deltas0(:,7),0.05)
b = Bin(deltas0(:,7),20);
clf
semplot(ht,h(b==1,:))
semplot(ht,h(b==2,:))
clf
semplot(ht,h(b==2,:))
semplot(ht,h(b==3,:))
semplot(ht,h(b==4,:))
semplot(ht,h(b==20,:))
semplot(ht,h(b>10,:))
clf
semplot(ht,h(b==2,:))
semplot(ht,h(b==3,:))
semplot(ht,h(b==4,:))
semplot(ht,h(b>10,:))
b = Bin(tiedrank(deltas0(:,7)),20);
clf
semplot(ht,h(b==1,:))
semplot(ht,h(b==20,:))
semplot(ht,h(b==2,:))
clf
semplot(ht,h(b==2,:))
semplot(ht,h(b==1,:))
semplot(ht,h(b==3,:))
clf
semplot(ht,h(b==3,:))
semplot(ht,h(b==20,:))
max(deltas0(b==3,7))
max(deltas0(b==2,7))
max(deltas0(b==4,7))
semplot(ht,h(b==3,:))
semplot(ht,h(b==20,:))
clf
max(deltas0(b==1,7))
semplot(ht,h(b==2 & deltas0(:,7)<20,:))
semplot(ht,h(b==2 & deltas0(:,7)>20,:))
clf
PlotColorMap(Shrink(sortby(h,deltas0(:,7)),500,4));
PlotColorMap(Shrink(sortby(h,deltas0(:,5)-deltas0(:,6)),500,4));
PlotColorMap(Shrink(sortby(h,deltas0(:,7)),500,4));
d = @(x,y) (x-y)./(x+y);
clf
semplot(ht,h);
response = d(mean(h(:,InIntervals(ht(:),[-0.2 -0.1; 0.1 0.2])),2),mean(h(:,InIntervals(ht(:),[-0.05 0.05])),2),
response = d(mean(h(:,InIntervals(ht(:),[-0.2 -0.1; 0.1 0.2])),2),mean(h(:,InIntervals(ht(:),[-0.05 0.05])),2));
PlotColorMap(Shrink(sortby(h,response),500,4));
semplot(ht,h);
clf
semplot(ht,h);
response = d(mean(h(:,InIntervals(ht(:),[-0.2 -0.1; 0.1 0.4])),2),mean(h(:,InIntervals(ht(:),[-0.05 0.05])),2));
PlotColorMap(Shrink(sortby(h,response),500,4));
PlotColorMap(Shrink(sortby(h,response),500,1));
z = zscore(h,[],2);
response = d(mean(z(:,InIntervals(ht(:),[-0.2 -0.1; 0.1 0.4])),2),mean(z(:,InIntervals(ht(:),[-0.05 0.05])),2));
PlotColorMap(Shrink(sortby(h,response),500,1));
response = d(mean(h(:,InIntervals(ht(:),[-0.2 -0.1; 0.1 0.4])),2),mean(h(:,InIntervals(ht(:),[-0.05 0.05])),2));
PlotColorMap(Shrink(sortby(h,response),500,1));
plot(response,deltas(:,7))
clf
plot(response,deltas0(:,7))
plot(response,deltas0(:,7),'.')
nancorr(response,deltas0(:,7))
[c,p] = nancorr(response,deltas0(:,7))
[c,p] = nancorr(response,deltas0(:,5)-deltas0(:,6))
clf
hist(deltas0(:,5)-deltas0(:,6),100);
quantile(deltas0(:,7),0.1)
quantile(deltas0(:,5)-deltas0(:,6),0.1)
p2z(0.01)
p2z(0.001)
hm = deltas0(:,5)-deltas0(:,6);
semplot(ht,h(hm<3.3,:));
semplot(ht,h(hm<3.5,:));
semplot(ht,h(hm<4,:));
clf
semplot(ht,h(hm<3.5,:));
clf
semplot(ht,h(hm<3.3,:));
semplot(ht,h(hm>3.3,:));
clf
semplot(ht,h(hm>3 * hm<3.3,:));
semplot(ht,h(hm>3 & hm<3.3,:));
clf
semplot(ht,h(hm>3 & hm<3.3,:));
Portion(hm>3)
deltas = deltas0(deltas0(:,5)-deltas0(:,6)>3,:); % these thresholds should be manually refined for each session
deltas = deltas0(deltas0(:,5)-deltas0(:,6)>3,:); % these thresholds should be manually refined for each session
deltaWaves.timestamps = deltas(:,[1 3]); deltaWaves.peaks = deltas(:,2);  deltaWaves.peakNormedPower = deltas(:,5); deltaWaves.detectorName = ['channel 114(+1), CleanLFP, FindDeltaPeaks, peak-trough>3'];
save(fullfile(basepath,[basename '.deltaWaves.events.mat']),'deltaWaves');
basename = basenameFromBasepath(basepath);
save(fullfile(basepath,[basename '.deltaWaves.events.mat']),'deltaWaves');
[h,ht] = PETH(deltas(:,2),r(:,1));
clf
PlotColorMap(sortby(zscore(h,[],2),r(:,4)),'x',ht);
PlotColorMap(Shrink(sortby(zscore(h,[],2),r(:,4)),100,1),'x',ht);
PlotColorMap(Shrink(sortby(zscore(h,[],2),r(:,4)),500,1),'x',ht);
PlotColorMap(Shrink(sortby(zscore(h,[],2),r(:,4)),1000,1),'x',ht);
dd = r(:,3)-r(:,1);
PlotColorMap(Shrink(sortby(zscore(h,[],2),dd),1000,1),'x',ht);
PlotColorMap(Shrink(sortby(dd,dd),1000,1),'x',ht);
PlotColorMap(Shrink(sortby(zscore(h,[],2),dd),1000,1),'x',ht);
hist(dd)
hist(dd,100)
sleep = ConsolidateIntervals(MergePoints.timestamps(cellfun(@(x) any(strfind(lower(x),'sleep')),MergePoints.foldernames),:));
load('day11.MergePoints.events.mat')
sleep = ConsolidateIntervals(MergePoints.timestamps(cellfun(@(x) any(strfind(lower(x),'sleep')),MergePoints.foldernames),:));
pre = InIntervals(r,sleep(1,:));
post = InIntervals(r,sleep(end,:));
anovabar(dd,[-pre+post]);
anovabar(dd>0.1,[-pre+post]);
r(findmax(dd))
addComma(r(findmax(dd)))
r(findmax(dd,5))
max(dd_
max(dd)_
max(dd)
sum(d==0.2)
sum(dd==0.2)
sum(dd>0.19)
sum(dd>0.199)
PlotColorMap(Shrink(sortby(zscore(h,[],2),dd),1000,1),'x',ht);
PlotColorMap(Shrink(sortby(dd,dd),1000,1),'x',ht);
PlotColorMap(Shrink(sortby(zscore(h,[],2),dd),1000,1),'x',ht);
PlotColorMap(Shrink(sortby(zscore(h,[],2),dd),1000,2),'x',ht);
PlotColorMap(Shrink(sortby(zscore(h,[],2),r(:,4)),1000,2),'x',ht);
PlotColorMap(Shrink(sortby(zscore(h,[],2),r(:,4)),2000,2),'x',ht);
plot(r(:,4),dd)
plot(r(:,4),dd,'.'')
plot(r(:,4),dd,'.')
nancorr(r(:,4),dd)
PlotColorMap(Shrink(sortby(zscore(h,[],2),r(:,4)),1000,2),'x',ht);
PlotColorMap(Shrink(sortby(zscore(h,[],2),dd),1000,2),'x',ht);
figure; anovabar(dd>0.1,[-pre+post]);
anovabar(dd,[-pre+post]);
load('day12.cell_metrics.cellinfo.mat')
session
nowRipples = ripples;
load('day11.ripples.events.mat')
oldRipples = ripples;
clf
PETH(deltas,nowRipples(:,2));
PETH(deltas(:,2),nowRipples(:,2));
PETH(deltas(:,2),oldRipples.peaks);
hold all
PETH(deltas(:,2),nowRipples(:,2));
open oldRipples
datestr(clock)
datestr(floor(datenum(clock)))
ripplesAny = BatchFindRipples(pwd,[21 43 17 47],[],SleepState.ints.WAKEstate);
figure; hist(d,100);
max(d)
sum(d>maxRippleDuration)
sum((d<minRippleDuration))
minRippleDuration
d(d>maxRippleDuration)
figure; hist(d,100);
hist(d(d<maxRippleDuration),100);
sum(toMerge)
r = [starts stops];
hist(diff(r,[],2),100);
hist(Restrict(diff(r,[],2),[0 1]),100);
% merge close ones
iri = r(2:end,1) - r(1:end-1,2);
newSize =  r(2:end,2) - r(1:end-1,1);
toMerge = iri<minInterRippleInterval & newSize<maxRippleDuration;
rippleStart = strfind([0 toMerge'],[0 1])';
rippleEnd = rippleStart+1; % next ripple
r0 = r;
% Adjust end (incorporate second ripple into first)
r(rippleStart,2) = r(rippleEnd,2);
% Remove now-redundant second ripple
r(rippleEnd,:) = [];
iri = r(2:end,1) - r(1:end-1,2);
newSize =  r(2:end,2) - r(1:end-1,1);
toMerge = iri<minInterRippleInterval & newSize<maxRippleDuration;
hist(Restrict(diff(r,[],2),[0 1]),100);
hist(Restrict(diff(r0,[],2),[0 1]),100);
hist(Restrict(diff(r,[],2),[0 1]),100);
minInterRippleInterval
minInterRippleInterval =0.01;
r = [starts stops];
% merge close ones
iri = r(2:end,1) - r(1:end-1,2);
newSize =  r(2:end,2) - r(1:end-1,1);
toMerge = iri<minInterRippleInterval & newSize<maxRippleDuration;
while sum(toMerge)>0,
% Get the index of the first ripple in a pair to be merged
rippleStart = strfind([0 toMerge'],[0 1])';
rippleEnd = rippleStart+1; % next ripple
% Adjust end (incorporate second ripple into first)
r(rippleStart,2) = r(rippleEnd,2);
% Remove now-redundant second ripple
r(rippleEnd,:) = [];
iri = r(2:end,1) - r(1:end-1,2);
newSize =  r(2:end,2) - r(1:end-1,1);
toMerge = iri<minInterRippleInterval & newSize<maxRippleDuration;
end
hist(Restrict(diff(r,[],2),[0 1]),100);
clf
hist(Restrict(diff(r,[],2),[0 1]),100);
sum(toMerge)
basepath
session
%SaveCustomEvents('temp.tem.evt',r(:,1:2),{'ripple start','ripple peak','ripple stop'});
toMerge = iri<0.3 & newSize<maxRippleDuration;
toMerge = iri<0.03 & newSize<maxRippleDuration;
Portion(toMerge)
SaveCustomEvents('ripples.rip.evt',r(toMerge,1:2),{'ripple start','ripple stop'});
SaveCustomEvents('temp.tem.evt',r(toMerge,1:2),{'ripple start','ripple peak','ripple stop'});
SaveCustomEvents('temp.tem.evt',r(toMerge,1:2),{'ripple start','ripple stop'});
SaveCustomEvents('temp.te2.evt',r(find(toMerge)+1,1:2),{'start','stop'});
while sum(toMerge)>0,
% Get the index of the first ripple in a pair to be merged
rippleStart = strfind([0 toMerge'],[0 1])';
rippleEnd = rippleStart+1; % next ripple
% Adjust end (incorporate second ripple into first)
r(rippleStart,2) = r(rippleEnd,2);
% Remove now-redundant second ripple
r(rippleEnd,:) = [];
iri = r(2:end,1) - r(1:end-1,2);
newSize =  r(2:end,2) - r(1:end-1,1);
toMerge = iri<minInterRippleInterval & newSize<maxRippleDuration;
end
% peak threshold check
[in,w] = InIntervals(t,r);
[peakV,peak] = Accumulate(w(in),signal(in),'mode','max');
ok = peakV>peakThreshold;
tIn = t(in);
peak = tIn(peak(ok));
r = [r(ok,1) peak r(ok,end) peakV(ok)];
% duration check
d = r(:,3) - r(:,1);
r(d<minRippleDuration | d>maxRippleDuration,:) = [];
open r
SaveCustomEvents('ripples.rip.evt',r(:,1:3),{'ripple start','ripple peak','ripple stop'});
channels
Channels
BatchFindRipples(pwd,[21 43 17 47],[],SleepState.ints.WAKEstate);
basepath
BatchFindRipples(basepath,[21 43 17 47],[],SleepState.ints.WAKEstate);
ripples = helper_function(t,signal,threshold,minInterRippleInterval,minRippleDuration,maxRippleDuration,peakThreshold);
SaveCustomEvents('ripples.rip.evt',r(:,1:3),{'ripple start','ripple peak','ripple stop'});
SaveCustomEvents('ripples.rip.evt',ripples(:,1:3),{'ripple start','ripple peak','ripple stop'});
ripples_uncut = helper_function(t,signal,threshold,0.01,minRippleDuration,maxRippleDuration,peakThreshold);
r = ripples_uncut(:,[1 end]);
SaveCustomEvents('temp.tem.evt',r(:,1:2),{'ripple start','ripple stop'});
open r
ripples_uncut = helper_function(t,signal,threshold,0.01,0,1,peakThreshold);
r = ripples_uncut(:,[1 end]);
SaveCustomEvents('temp.tem.evt',r(:,1:2),{'ripple start','ripple stop'});
ripples
open ripples
FindClosest(ripples(:,2)2040.28)
FindClosest(ripples(:,2)*2040.28)
FindClosest(ripples(:,2),2040.28)
ripples(1649,2)-2040.28
FindClosest(r(:,1),2040.28)
r(1907,1)-2040.28
r(1907,2)-2040.28
addComma(r(1907,1))
r = ripples_uncut(:,[1 3]);
SaveCustomEvents('temp.tem.evt',r(:,1:2),{'ripple start','ripple stop'});
r(1907,1)
addComma(r(1907,1))
addComma(r(1907,2))
addComma(r(1907+1,1))
addComma(r(1907-1,1))
addComma(r(1907-1,2))
edit temp
ripples_uncut = helper_function(t,signal,threshold,0.01,0,1,0);
r = ripples_uncut(:,[1 3]);
FindClosest(r(:,1),2040.28)
addComma(r(8828,1))
addComma(ripples(1649,1))
addComma(ripples(1649,2))
addComma(r(8828-1,1))
SaveCustomEvents('temp.tem.evt',r(:,1:2),{'ripple start','ripple stop'});
size(ok)
ok = InIntervals(t,2040+[0 1]);
sum(ok)
plot(t(ok),signal(ok));
PlotHVLines(threshold,'h','r--');
PlotHVLines(Restrict(r,xlim),'v','k-');
PlotHVLines(Restrict(ripples(:,[1 3]),xlim),'v','r--');
xlim(2040+[0 1]);
interval =  2040+[0 1];
clf;
ok = InIntervals(t,interval);
plot(t(ok),signal(ok));
PlotHVLines(threshold,'h','r--');
PlotHVLines(Restrict(r,interval),'v','k-');
PlotHVLines(Restrict(ripples(:,[1 3]),interval),'v','r--');
xlim(interval);
interval =  2042+[0 1];
clf;
ok = InIntervals(t,interval);
plot(t(ok),signal(ok));
PlotHVLines(threshold,'h','r--');
PlotHVLines(Restrict(r,interval),'v','k-');
PlotHVLines(Restrict(ripples(:,[1 3]),interval),'v','r--');
xlim(interval);
hold all
plot(t(ok),noisy(ok),'y');
plot(t(ok),me(ok),'r');
clf;
ok = InIntervals(t,interval);
plot(t(ok),signal(ok));
PlotHVLines(threshold,'h','r--');
PlotHVLines(Restrict(r,interval),'v','k-');
PlotHVLines(Restrict(ripples(:,[1 3]),interval),'v','r--');
xlim(interval);
dd = me-noisy;
plot(t(ok),dd(ok));
dd = nanzscore(me-noisy);
clf;
ok = InIntervals(t,interval);
plot(t(ok),signal(ok));
PlotHVLines(threshold,'h','r--');
PlotHVLines(Restrict(r,interval),'v','k-');
PlotHVLines(Restrict(ripples(:,[1 3]),interval),'v','r--');
xlim(interval);
hold all
plot(t(ok),dd(ok),'y');
PlotHVLines(0,'h','k');
signal0 = signal;
% Find ripples
signal = me-noisy; % remove noise from envelope
signal(me>noisy) = nanzscore(signal(me>noisy));
% signal(nanzscore(noisy)>threshold) = 0;
signal(~(me>noisy)) = 0;
plot(t(ok),signal(ok),'r');
signal = me-noisy; % remove noise from envelope
signal(me>noisy) = nanzscore(signal(me>noisy));
plot(t(ok),signal(ok),'r');
clf;
ok = InIntervals(t,interval);
plot(t(ok),signal(ok));
PlotHVLines(threshold,'h','r--');
PlotHVLines(Restrict(r,interval),'v','k-');
PlotHVLines(Restrict(ripples(:,[1 3]),interval),'v','r--');
xlim(interval);
hold all
plot(t(ok),dd(ok),'y');
PlotHVLines(0,'h','k');
threshold
peakThreshold
interval =  2049+[0 1];
clf;
ok = InIntervals(t,interval);
plot(t(ok),signal(ok));
PlotHVLines(threshold,'h','r--');
PlotHVLines(Restrict(r,interval),'v','k-');
PlotHVLines(Restrict(ripples(:,[1 3]),interval),'v','r--');
xlim(interval);
signal(~(me>noisy)) = 0;
clf;
ok = InIntervals(t,interval);
plot(t(ok),signal(ok));
PlotHVLines(threshold,'h','r--');
PlotHVLines(Restrict(r,interval),'v','k-');
PlotHVLines(Restrict(ripples(:,[1 3]),interval),'v','r--');
xlim(interval);
PlotHVLines(2,'h','y--');
ripples1 = ripples;
ripples = helper_function(t,signal,2,minInterRippleInterval,minRippleDuration,maxRippleDuration,peakThreshold);
ripples = helper_function(t,signal,2,0.02,minRippleDuration,maxRippleDuration,peakThreshold);
PlotHVLines(Restrict(ripples(:,[1 3]),interval),'v','b--','linewidth',2);
interval
xlim(interval);
figure; hist(diff(ripples(:,[1 3]),[],2));
hist(diff(ripples(:,[1 3]),[],2),100);
load('day11.MergePoints.events.mat')
sleep = ConsolidateIntervals(MergePoints.timestamps(cellfun(@(x) any(strfind(lower(x),'sleep')),MergePoints.foldernames),:));
pre = InIntervals(ripples,sleep(1,:));
post = InIntervals(ripples,sleep(1,:));
anovabar(diff(ripples(:,[1 3]),[],2)>0.1,[pre+post]);
anovabar(diff(ripples(:,[1 3]),[],2)>0.1,[-pre+post]);
post = InIntervals(ripples,sleep(2,:));
anovabar(diff(ripples(:,[1 3]),[],2)>0.1,[-pre+post]);
clf
anovabar(diff(ripples(:,[1 3]),[],2)>0.1,[-pre+post]);
load('day11.deltaWaves.events.mat')
deltas = deltaWaves.peaks;
rd = InIntervals(ripples(:,2),[delta-0.2 delta]);
rd = InIntervals(ripples(:,2),[deltas-0.2 deltas]);
PETH(ripples(:,2),deltas)
dr = InIntervals(ripples(:,2),[deltas+0.05 deltas+0.35]);
anovabar(diff(ripples(:,[1 3]),[],2)>0.1,[-pre+post]);
anovabar(rd,[-pre+post]);
anovabar(dr,[-pre+post]);
PETH(ripples(:,2),deltas)
PETH(ripples(:,2),deltas,'nBins',501)
rd = InIntervals(ripples(:,2),[deltas-0.16 deltas-0.04]);
dr = InIntervals(ripples(:,2),[deltas+0.1 deltas+0.4]);
anovabar(dr,[-pre+post]);
anovabar(rd,[-pre+post]);
PETH(ripples(:,2),deltas,'nBins',501)
PETH(ripples(:,2),Restrict(deltas,sleep(1,:)),'nBins',501)
hold all
PETH(ripples(:,2),Restrict(deltas,sleep(2,:)),'nBins',501)
sleep
PETH(ripples(:,2),Restrict(deltas,sleep(3,:)),'nBins',501)
post = InIntervals(ripples,sleep(3,:));
clf
anovabar(rd,[-pre+post]);
anovabar(dr,[-pre+post]);
anovabar(diff(ripples(:,[1 3]),[],2)>0.1,[-pre+post]);
corr(diff(ripples(:,[1 3])),dr)
Portion(dr)
size(dr)
size(diff(ripples(:,[1 3])))
corr(diff(ripples(:,[1 3]),[],2),dr)
anovabar(diff(ripples(:,[1 3]),[],2),dr)
anovabar(diff(ripples(:,[1 3]),[],2)>0.1,dr)
anovabar(diff(ripples(:,[1 3]),[],2),dr)
anovabar(diff(ripples(:,[1 3]),[],2),rd)
anovabar(diff(ripples(:,[1 3]),[],2)>0.1,rd)
[h,ht] = PETH(deltas(:,2),r(:,1));
[h,ht] = PETH(deltas(:,2),ripples(:,1));
[h,ht] = PETH(deltas,ripples(:,1));
dd = diff(ripples(:,[1 3]),[],2);
PlotColorMap(Shrink(sortby(zscore(h,[],2),dd),1000,2),'x',ht);
ripples = helper_function(t,signal,2,0.025,minRippleDuration,maxRippleDuration,peakThreshold);
dd = diff(ripples(:,[1 3]),[],2);
[h,ht] = PETH(deltas,ripples(:,1));
figure; PlotColorMap(Shrink(sortby(zscore(h,[],2),dd),1000,2),'x',ht);
dr = InIntervals(ripples(:,2),[deltas+0.1 deltas+0.4]);
corr(diff(ripples(:,[1 3]),[],2),dr)
size(ripples)
size(ripples1)
ripples = helper_function(t,signal,1.5,0.025,minRippleDuration,maxRippleDuration,peakThreshold);
[h,ht] = PETH(deltas,ripples(:,1));
dd = diff(ripples(:,[1 3]),[],2);
figure; PlotColorMap(Shrink(sortby(zscore(h,[],2),dd),1000,2),'x',ht);
PlotColorMap(Shrink(sortby(zscore(h,[],2),dd),10000,2),'x',ht);
PlotColorMap(Shrink(sortby(zscore(h,[],2),dd),1000,2),'x',ht);
dr = InIntervals(ripples(:,2),[deltas+0.1 deltas+0.4]);
rd = InIntervals(ripples(:,2),[deltas-0.16 deltas-0.04]);
anovabar(diff(ripples(:,[1 3]),[],2)>0.1,rd)
anovabar(diff(ripples(:,[1 3]),[],2),rd)
anovabar(diff(ripples(:,[1 3]),[],2),dr)
anovabar(diff(ripples(:,[1 3]),[],2),[-pre post])
pre = InIntervals(ripples,sleep(1,:));
post = InIntervals(ripples,sleep(3,:));
anovabar(diff(ripples(:,[1 3]),[],2),[-pre+post])
PlotHVLines(Restrict(ripples(:,[1 3]),interval),'v','k-.','linewidth',2);
clf;
interval =  17800+[0 1];
clf;
ok = InIntervals(t,interval);
plot(t(ok),signal(ok));
PlotHVLines(threshold,'h','r--');
PlotHVLines(Restrict(r,interval),'v','k-');
PlotHVLines(Restrict(ripples(:,[1 3]),interval),'v','r--');
xlim(interval);
PlotHVLines(Restrict(ripples(:,[1 3]),interval),'v','b-.','linewidth',2);
PlotHVLines(1.5,'h','y--');
plot(t(ok),Smooth(signal(ok),[1250*0.1]));
plot(t(ok),Smooth(signal(ok),[1250*0.01]));
plot(t(ok),Smooth(signal(ok),[1250*0.01])*2);
plot(t(ok),Smooth(signal(ok),[1250*0.01])*2,'linewidth',2);
interval =  20489+[0 1];
clf;
ok = InIntervals(t,interval);
plot(t(ok),signal(ok));
PlotHVLines(threshold,'h','r--');
PlotHVLines(Restrict(r,interval),'v','k-');
PlotHVLines(Restrict(ripples(:,[1 3]),interval),'v','r--');
xlim(interval);
PlotHVLines(Restrict(ripples(:,[1 3]),interval),'v','b-.','linewidth',2);
PlotHVLines(1.5,'h','y--');
plot(t(ok),Smooth(signal(ok),[1250*0.01])*2,'linewidth',2);
dt = mode(diff(dt))
dt = mode(diff(t))
0.01/dt
smooth
dt = mode(diff(t)); smooth = 0.01/dt; % 10ms smooth
% Find ripples
signal = me-noisy; % remove noise from envelope
signal(me>noisy) = nanzscore(Smooth(signal(me>noisy),smooth));
% signal(nanzscore(noisy)>threshold) = 0;
signal(~(me>noisy)) = 0;
plot(t(ok),signal(ok));
figure; plot(t(ok),signal(ok));
plot(t(ok),signal(ok),'linewidth',2);
dt = mode(diff(t)); smooth = 0.005/dt; % 10ms smooth
% Find ripples
signal = me-noisy; % remove noise from envelope
signal(me>noisy) = nanzscore(Smooth(signal(me>noisy),smooth));
% signal(nanzscore(noisy)>threshold) = 0;
signal(~(me>noisy)) = 0;
plot(t(ok),signal(ok),'linewidth',2);
clf
pre = InIntervals(ripples,sleep(1,:));
post = InIntervals(ripples,sleep(3,:));
anovabar(diff(ripples(:,[1 3]),[],2),[-pre+post])
rd = InIntervals(ripples(:,2),[deltas-0.16 deltas-0.04]);
anovabar(rd,[-pre+post])
anovabar(dr,[-pre+post])
r(end,1)
21118/r(end,1)
ripples = helper_function(t,signal,threshold,minInterRippleInterval,minRippleDuration,maxRippleDuration,peakThreshold);
%%Save ripples
threshold
threshold = 1.5
minInterRippleInterval
minInterRippleInterval = 0.02;
maxRippleDuration = 0.5;
ripples = helper_function(t,signal,threshold,minInterRippleInterval,minRippleDuration,maxRippleDuration,peakThreshold);
hist(diff(ripples(:,[1 3]),[],2),100);
anovabar(diff(ripples(:,[1 3]),[],2),[-pre+post])
pre = InIntervals(ripples,sleep(1,:));
post = InIntervals(ripples,sleep(3,:));
anovabar(diff(ripples(:,[1 3]),[],2),[-pre+post])
anovabar(diff(ripples(:,[1 3]),[],2)>0.1,[-pre+post])
hist(diff(ripples(:,[1 3]),[],2),100);
maxRippleDuration
minInterRippleInterval
threshold
minRippleDuration
hist(log(diff(ripples(:,[1 3]),[],2)),100);
hist(log(diff(ripples(:,[1 3]),[],2)),1000);
clf
hist(log(diff(ripples(:,[1 3]),[],2)),1000);
hist(log10(diff(ripples(:,[1 3]),[],2)),1000);
this = log10(diff(ripples(:,[1 3]),[],2));
hist(this,1000);
hist(this,100);
hist(this,1000);
hist(this,100);
hist(this,50);
% here = max(strfind(session,'/'));
% try cd([session(1:here) '/events/']); catch cd(session(1:here));end
cd(session);
threshold = 1.5;
peakThreshold = 3;
if str2double(session(end-15:end-13)) == 272; peakThreshold = 4; end
if str2double(session(end-15:end-13)) == 291; peakThreshold = 5; end %noisy
minInterRippleInterval = 0.020;
minRippleDuration = 0.020;
maxRippleDuration = 0.5;
aroundArtefact = 0.05; epsilon = 0.15; %for consolidatin
ripples = helper_function(t,signal,threshold,minInterRippleInterval,minRippleDuration,maxRippleDuration,peakThreshold);
45
length(ripples)/ripples(end,1)
hist(log10(diff(ripples(:,[1 3]),[],2)),1000);
hist(log10(diff(ripples(:,[1 3]),[],2)),10);
hist(log10(diff(ripples(:,[1 3]),[],2)),100);
hist((diff(ripples(:,[1 3]),[],2)),100);
PETH(deltas(:,2),nowRipples(:,2));
PETH(deltas,ripples(:,2));
PETH(deltas,ripples(:,2),'nBins',501);
[h,ht] = PETH(deltas,ripples(:,1));
clf
dd = diff(ripples(:,[1 3]),[],2);
PlotColorMap(Shrink(sortby(zscore(h,[],2),dd),1000,2),'x',ht);
PlotColorMap(Shrink(sortby(zscore(h,[],2),dd),100,2),'x',ht);
PlotColorMap(Shrink(sortby(zscore(h,[],2),dd),500,1),'x',ht);
SaveCustomEvents('temp.tem.evt',ripples(:,[1 3]),{'ripple start','ripple stop'});
bad = ~InIntervals(ripples1(:,2),ripples(:,[1 3]));
Portion(bad)
SaveCustomEvents('temp.te2.evt',ripples1(bad,[1 3]),{'ripple start','ripple stop'});
5/0.037
interval =  20489+[0 1];
clf;
ok = InIntervals(t,interval);
plot(t(ok),signal(ok));
PlotHVLines(threshold,'h','r--');
PlotHVLines(Restrict(r,interval),'v','k-');
PlotHVLines(Restrict(ripples(:,[1 3]),interval),'v','r--');
xlim(interval);
PlotHVLines(Restrict(ripples(:,[1 3]),interval),'v','b-.','linewidth',2);
PlotHVLines(1.5,'h','y--');
% Find ripples
signal = me-noisy; % remove noise from envelope
signal(me>noisy) = nanzscore(Smooth(signal(me>noisy),smooth));
% signal(nanzscore(noisy)>threshold) = 0;
signal(~(me>noisy)) = 0;
ripples = helper_function(t,signal,threshold,minInterRippleInterval,minRippleDuration,maxRippleDuration,peakThreshold);
5
interval =  20489+[0 1];
clf;
ok = InIntervals(t,interval);
plot(t(ok),signal(ok));
PlotHVLines(threshold,'h','r--');
PlotHVLines(Restrict(r,interval),'v','k-');
PlotHVLines(Restrict(ripples(:,[1 3]),interval),'v','r--');
xlim(interval);
PlotHVLines(Restrict(ripples(:,[1 3]),interval),'v','b-.','linewidth',2);
PlotHVLines(1.5,'h','y--');
PlotHVLines(3,'h','y-');
signal = me-noisy; % remove noise from envelope
signal(me>noisy) = nanzscore(signal(me>noisy));
% signal(nanzscore(noisy)>threshold) = 0;
signal(~(me>noisy)) = 0;
ripples = helper_function(t,signal,threshold,minInterRippleInterval,minRippleDuration,maxRippleDuration,peakThreshold);
interval =  20489+[0 1];
clf;
ok = InIntervals(t,interval);
plot(t(ok),signal(ok));
PlotHVLines(threshold,'h','r--');
PlotHVLines(Restrict(r,interval),'v','k-');
PlotHVLines(Restrict(ripples(:,[1 3]),interval),'v','r--');
xlim(interval);
PlotHVLines(Restrict(ripples(:,[1 3]),interval),'v','b-.','linewidth',2);
PlotHVLines(1.5,'h','y--');
PlotHVLines(3,'h','y-');
bad = ~InIntervals(ripples1(:,2),ripples(:,[1 3]));
Portion(bad)
SaveCustomEvents('temp.te2.evt',ripples1(bad,[1 3]),{'ripple start','ripple stop'});
SaveCustomEvents('temp.tem.evt',ripples(:,[1 3]),{'ripple start','ripple stop'});
SaveCustomEvents('ripples.rip.evt',ripples(:,1:3),{'ripple start','ripple peak','ripple stop'});
params.basepath = basepath;
params.channels = channels;
params.noisyChannel = noisyChannel;
params.excludeInterval = excludeInterval;
params.threshold = threshold;
params.peakThreshold = peakThreshold;
params.minInterRippleInterval = minInterRippleInterval;
params.minRippleDuration = minRippleDuration;
params.maxRippleDuration = maxRippleDuration;
params.aroundArtefact = aroundArtefact;
params.epsilon = epsilon;
params.excludeInterval = excludeInterval;
structure.timestamps = ripples(:,[1 3]);
structure.duration = diff(ripples(:,[1 3]),[],2);
structure.amplitude = ripples(:,4);
structure.peaks = ripples(:,2);
structure.detectorName = 'BatchFindRipples';
structure.detectorinfo.detectionparms = params;
structure.detectorinfo.detectiondate = datestr(floor(datenum(clock))); % inputs
structure.detectorinfo.detectionintervals = SubtractIntervals([t(1) t(end,1)],excludeInterval);
session
basepath = session;
params.basepath = basepath;
params.channels = channels;
params.noisyChannel = noisyChannel;
params.excludeInterval = excludeInterval;
params.threshold = threshold;
params.peakThreshold = peakThreshold;
params.minInterRippleInterval = minInterRippleInterval;
params.minRippleDuration = minRippleDuration;
params.maxRippleDuration = maxRippleDuration;
params.aroundArtefact = aroundArtefact;
params.epsilon = epsilon;
params.excludeInterval = excludeInterval;
structure.timestamps = ripples(:,[1 3]);
structure.duration = diff(ripples(:,[1 3]),[],2);
structure.amplitude = ripples(:,4);
structure.peaks = ripples(:,2);
structure.detectorName = 'BatchFindRipples';
structure.detectorinfo.detectionparms = params;
structure.detectorinfo.detectiondate = datestr(floor(datenum(clock))); % inputs
structure.detectorinfo.detectionintervals = SubtractIntervals([t(1) t(end,1)],excludeInterval);
rrr = ripples;
ripples = structure;
save('day11.ripples.events.mat','ripples');
ripples = rrr;
dbcont
ripples = ans;
[fList,pList] = matlab.codetools.requiredFilesAndProducts('BatchFindRipples.m');
fList = fList';
open fList
load('day11.ripples.events.mat', 'ripples')
load('day12.ripples.events.mat')
load('day12.ripples.events.mat', 'ripples')
clear all
load('day12.ripples.events.mat')
X = BatchReturn('OJR_longTraining.batch');
cd(X{1})
cd(X{2})
cd(X{3})
cd(X{4})
cd(X{5})
cd(X{6})
clear all
load('day10.channelinfo.ripples.mat')
save('day11.channelinfo.ripples.mat','ripplesChannels');
save('day11.channelinfo.ripples.mat','rippleChannels');
load('detecting_ripples.mat')
plot(swDiffAll,ripPowerAll,'k');
plot(swDiffAll,ripPowerAll,'k.');
hold all
plot(swDiffAll(idx1),ripPowerAll(idx1),'r.');
ripples = DetectSWR([rippleChannels.Ripple_Channel rippleChannels.Sharpwave_Channel rippleChannels.Noise_Channel],'basepath',basepath,'saveMat',true,'forceDetect',true,'check',true,'useSPW',false);
basepath = pwd
ripples = DetectSWR([rippleChannels.Ripple_Channel rippleChannels.Sharpwave_Channel rippleChannels.Noise_Channel],'basepath',basepath,'saveMat',true,'forceDetect',true,'check',true,'useSPW',false);
% First, remove points that occur in noisy lfp periods
tl = (1:length(lfp))'/1250; % timestamps for lfp
t = featureTs/1250; % timestamps for candidate ripple events
bad = false(size(t));
try % optionally, remove periods of noisy LFP (large deflections)
[clean,~,badIntervals] = CleanLFP([tl,double(lfp(:,2))],'thresholds',[10 Inf]);
bad = bad | InIntervals(t,badIntervals);
end
try % optionally, remove periods when your channels are diverging abnormally for a long time (a channel went dead for some time)
% if you need to use this, make sure your threshold if appropriate for your session. Plot the smoothed signal to see what seems
% reasonable to you!
smoothed = Smooth(abs(swDiff),1250*10);
noisyPeriods = ConsolidateIntervals(tl(FindInterval(smoothed>500)),'epsilon',1);
bad = bad | InIntervals(t,noisyPeriods);
end
try % optionally, remove ripples in which the sharp wave channel was positive
% Only do this if you have a sharp-wave channel what's sufficiently deep for the sharp wave to be negative
sw = interp1(tl,lfpLow(:,2),t);
%                 bad = bad | sw>0;
end
idx1 = idx1 & ~bad; idx2 = idx2 & ~bad; % remove bad points from "idx1" and "idx2"
figure(1); % plot a clean figure without these points and check indivudual ripples
clf
plot(swDiffAll(idx1 | idx2),ripPowerAll(idx1 | idx2),'k.','markersize',1); hold on
scores = nan(size(ripPowerAll,1),1);
% figure(1); % plot a clean figure without these points and check indivudual ripples
% clf
% plot(swDiffAll(idx1 | idx2),ripPowerAll(idx1 | idx2),'k.','markersize',1); hold on
% scores = nan(size(ripPowerAll,1),1);
matrix = [swDiffAll ripPowerAll];
z = bsxfun(@rdivide,bsxfun(@minus,matrix,mean(matrix(~bad,:))),std(matrix(~bad,:)));
while true % click "Ctrl+C" to exit when you feel you're done
figure(1);
colors = Bright(1000);
[x,y,button] = ginput(1);
i = findmin(abs(x-swDiffAll)+abs(y-ripPowerAll));
t = featureTs/1250;
tl = (1:length(lfp))'/1250;
figure(2);
clf; interval = t(i) + [-1 1]*0.5; in = tl>interval(1) & tl<interval(end);
subplot(3,1,1);
plot(tl(in) - t(i),lfp(in,:)); title(num2str(t(i)));
PlotHVLines(Restrict(t,interval) - t(i),'v','k--');
legend('ripple channel','SW channel','noise channel');
subplot(3,1,2);
plot(tl(in) - t(i),double(swDiff(in)));
PlotHVLines(Restrict(t,interval) - t(i),'v','k--');
subplot(3,1,3);
plot(tl(in) - t(i),double(ripPower0(in)));
PlotHVLines(Restrict(t,interval) - t(i),'v','k--');
%     x = input( prompt )
commandwindow % select the command window
score = str2double(input('good(1)/bad(0)?','s'));
scores(i) = score; % save scores to optionally save
figure(1);
xlims = xlim; ylims = ylim; clf
% extrapolate the score of each ripple based on the score of the nearest scored ripple
scored = find(~isnan(scores));
distances = sqrt(bsxfun(@minus,z(:,1),z(scored,1)').^2 + bsxfun(@minus,z(:,2),z(~isnan(scores),2)').^2);
[~,nearestNeighbour] = min(distances,[],2);
estimated = scores(scored(nearestNeighbour));
scatter(matrix(:,1),matrix(:,2),1,estimated,'filled'); colormap(Bright); clim([0 1]);
xlim(xlims); ylim(ylims);
hold on;
scatter(swDiffAll(scored),ripPowerAll(scored),20,scores(scored)); scatter(swDiffAll(scored),ripPowerAll(scored),15,scores(scored));
clabel('Estimated ripple score');
xlabel('Sharp wave depth'); ylabel('Ripple power');
title({['Manual scoring for ' basepath],['called with channels: ' num2str(Channels) ' (ripple, sharp wave, and (optionally) noise)']});
if rem(i,10)==0
display('saving...');
save(fullfile(basepath,'DetectSWR_manual_scoring.mat'),'t','swDiffAll','ripPowerAll','bad','scores');
end
end
channels
Channels
rippleChannels.Ripple_Channel = 41+1; rippleChannels.Sharpwave_Channel = 24+1; rippleChannels.Noise_Channel = 71+1;
save('day11.channelinfo.ripples.mat','rippleChannels');
ripples = DetectSWR([rippleChannels.Ripple_Channel rippleChannels.Sharpwave_Channel rippleChannels.Noise_Channel],'basepath',basepath,'saveMat',true,'forceDetect',true,'check',true,'useSPW',false);
close all
ripples = DetectSWR([rippleChannels.Ripple_Channel rippleChannels.Sharpwave_Channel rippleChannels.Noise_Channel],'basepath',basepath,'saveMat',true,'forceDetect',true,'check',true,'useSPW',false);
title('11');
% First, remove points that occur in noisy lfp periods
tl = (1:length(lfp))'/1250; % timestamps for lfp
t = featureTs/1250; % timestamps for candidate ripple events
bad = false(size(t));
try % optionally, remove periods of noisy LFP (large deflections)
[clean,~,badIntervals] = CleanLFP([tl,double(lfp(:,2))],'thresholds',[6 Inf]);
bad = bad | InIntervals(t,badIntervals);
end
try % optionally, remove periods when your channels are diverging abnormally for a long time (a channel went dead for some time)
% if you need to use this, make sure your threshold if appropriate for your session. Plot the smoothed signal to see what seems
% reasonable to you!
smoothed = Smooth(abs(swDiff),1250*10);
noisyPeriods = ConsolidateIntervals(tl(FindInterval(smoothed>500)),'epsilon',1);
bad = bad | InIntervals(t,noisyPeriods);
end
try % optionally, remove ripples in which the sharp wave channel was positive
% Only do this if you have a sharp-wave channel what's sufficiently deep for the sharp wave to be negative
sw = interp1(tl,lfpLow(:,2),t);
%                 bad = bad | sw>0;
end
idx1 = idx1 & ~bad; idx2 = idx2 & ~bad; % remove bad points from "idx1" and "idx2"
Portion(bad)
figure(1); % plot a clean figure without these points and check indivudual ripples
clf
plot(swDiffAll(idx1 | idx2),ripPowerAll(idx1 | idx2),'k.','markersize',1); hold on
scores = nan(size(ripPowerAll,1),1);
% figure(1); % plot a clean figure without these points and check indivudual ripples
% clf
% plot(swDiffAll(idx1 | idx2),ripPowerAll(idx1 | idx2),'k.','markersize',1); hold on
% scores = nan(size(ripPowerAll,1),1);
matrix = [swDiffAll ripPowerAll];
z = bsxfun(@rdivide,bsxfun(@minus,matrix,mean(matrix(~bad,:))),std(matrix(~bad,:)));
figure(1); % plot a clean figure without these points and check indivudual ripples
clf
plot(swDiffAll(idx1 | idx2),ripPowerAll(idx1 | idx2),'k.','markersize',1); hold on
scores = nan(size(ripPowerAll,1),1);
% figure(1); % plot a clean figure without these points and check indivudual ripples
% clf
% plot(swDiffAll(idx1 | idx2),ripPowerAll(idx1 | idx2),'k.','markersize',1); hold on
% scores = nan(size(ripPowerAll,1),1);
matrix = [swDiffAll ripPowerAll];
z = bsxfun(@rdivide,bsxfun(@minus,matrix,mean(matrix(~bad,:))),std(matrix(~bad,:)));
while true % click "Ctrl+C" to exit when you feel you're done
figure(1);
colors = Bright(1000);
[x,y,button] = ginput(1);
i = findmin(abs(x-swDiffAll)+abs(y-ripPowerAll));
t = featureTs/1250;
tl = (1:length(lfp))'/1250;
figure(2);
clf; interval = t(i) + [-1 1]*0.5; in = tl>interval(1) & tl<interval(end);
subplot(3,1,1);
plot(tl(in) - t(i),lfp(in,:)); title(num2str(t(i)));
PlotHVLines(Restrict(t,interval) - t(i),'v','k--');
legend('ripple channel','SW channel','noise channel');
subplot(3,1,2);
plot(tl(in) - t(i),double(swDiff(in)));
PlotHVLines(Restrict(t,interval) - t(i),'v','k--');
subplot(3,1,3);
plot(tl(in) - t(i),double(ripPower0(in)));
PlotHVLines(Restrict(t,interval) - t(i),'v','k--');
%     x = input( prompt )
commandwindow % select the command window
score = str2double(input('good(1)/bad(0)?','s'));
scores(i) = score; % save scores to optionally save
figure(1);
xlims = xlim; ylims = ylim; clf
% extrapolate the score of each ripple based on the score of the nearest scored ripple
scored = find(~isnan(scores));
distances = sqrt(bsxfun(@minus,z(:,1),z(scored,1)').^2 + bsxfun(@minus,z(:,2),z(~isnan(scores),2)').^2);
[~,nearestNeighbour] = min(distances,[],2);
estimated = scores(scored(nearestNeighbour));
scatter(matrix(:,1),matrix(:,2),1,estimated,'filled'); colormap(Bright); clim([0 1]);
xlim(xlims); ylim(ylims);
hold on;
scatter(swDiffAll(scored),ripPowerAll(scored),20,scores(scored)); scatter(swDiffAll(scored),ripPowerAll(scored),15,scores(scored));
clabel('Estimated ripple score');
xlabel('Sharp wave depth'); ylabel('Ripple power');
title({['Manual scoring for ' basepath],['called with channels: ' num2str(Channels) ' (ripple, sharp wave, and (optionally) noise)']});
if rem(i,10)==0
display('saving...');
save(fullfile(basepath,'DetectSWR_manual_scoring.mat'),'t','swDiffAll','ripPowerAll','bad','scores');
end
end
1
0
1
0
t(i)
figure; PlotXY([tl,double(lfp(:,2))]);
PlotHVLines(t(i),'v','k--');
PlotIntervals(badIntervals,'color','k');
bad(i)
idx1(i)
idx2(i)
idx1 = idx1 & ~bad; idx2 = idx2 & ~bad; % remove bad points from "idx1" and "idx2"
figure(1); % plot a clean figure without these points and check indivudual ripples
clf
plot(swDiffAll(idx1 | idx2),ripPowerAll(idx1 | idx2),'k.','markersize',1); hold on
scores = nan(size(ripPowerAll,1),1);
edit nothing
nothing
1
0
1
0
1
0.5
0
0.1
0.5
0
figure(1); % plot a clean figure without these points and check indivudual ripples
clf
plot(swDiffAll(idx1 | idx2),ripPowerAll(idx1 | idx2),'k.','markersize',1); hold on
scores = nan(size(ripPowerAll,1),1);
% either draw a polyglon:
selected = UIInPolygon(swDiffAll,ripPowerAll); % draw polygon to encircle the points you believe are ripples
% or use the nearest neighbours with a threshold you've used
% (I tend to label an ugly ripple that is still clearly a ripple with a score of 0.2)
scored = find(~isnan(scores));
distances = sqrt(bsxfun(@minus,z(:,1),z(scored,1)').^2 + bsxfun(@minus,z(:,2),z(~isnan(scores),2)').^2);
[~,nearestNeighbour] = min(distances,[],2);
estimated = scores(scored(nearestNeighbour));
selected = estimated>0.1;
idx1 = selected & ~bad; % final ripples
idx2 = ~selected & ~bad; % non-ripples (yet free from noise as well)
display('saving...');
save(fullfile(basepath,'DetectSWR_manual_scoring.mat'),'t','swDiffAll','ripPowerAll','idx1','idx2','bad','selected','scores');
saveas(figure(1),fullfile(basepath,'DetectSWR_manual_scoring.fig'));
dbcont
load('DetectSWR_manual_scoring.mat')
ripples = DetectSWR([rippleChannels.Ripple_Channel rippleChannels.Sharpwave_Channel rippleChannels.Noise_Channel],'basepath',basepath,'saveMat',true,'forceDetect',true,'check',true,'useSPW',false);
% First, remove points that occur in noisy lfp periods
tl = (1:length(lfp))'/1250; % timestamps for lfp
t = featureTs/1250; % timestamps for candidate ripple events
bad = false(size(t));
try % optionally, remove periods of noisy LFP (large deflections)
[clean,~,badIntervals] = CleanLFP([tl,double(lfp(:,2))],'thresholds',[6 Inf]);
bad = bad | InIntervals(t,badIntervals);
end
try % optionally, remove periods when your channels are diverging abnormally for a long time (a channel went dead for some time)
% if you need to use this, make sure your threshold if appropriate for your session. Plot the smoothed signal to see what seems
% reasonable to you!
smoothed = Smooth(abs(swDiff),1250*10);
noisyPeriods = ConsolidateIntervals(tl(FindInterval(smoothed>500)),'epsilon',1);
bad = bad | InIntervals(t,noisyPeriods);
end
try % optionally, remove ripples in which the sharp wave channel was positive
% Only do this if you have a sharp-wave channel what's sufficiently deep for the sharp wave to be negative
sw = interp1(tl,lfpLow(:,2),t);
%                 bad = bad | sw>0;
end
idx1 = idx1 & ~bad; idx2 = idx2 & ~bad; % remove bad points from "idx1" and "idx2"
figure(1); % plot a clean figure without these points and check indivudual ripples
clf
plot(swDiffAll(idx1 | idx2),ripPowerAll(idx1 | idx2),'k.','markersize',1); hold on
scores = nan(size(ripPowerAll,1),1);
% figure(1); % plot a clean figure without these points and check indivudual ripples
% clf
% plot(swDiffAll(idx1 | idx2),ripPowerAll(idx1 | idx2),'k.','markersize',1); hold on
% scores = nan(size(ripPowerAll,1),1);
matrix = [swDiffAll ripPowerAll];
z = bsxfun(@rdivide,bsxfun(@minus,matrix,mean(matrix(~bad,:))),std(matrix(~bad,:)));
Portion(bad)
figure(1); % plot a clean figure without these points and check indivudual ripples
clf
plot(swDiffAll(idx1 | idx2),ripPowerAll(idx1 | idx2),'k.','markersize',1); hold on
scores = nan(size(ripPowerAll,1),1);
% figure(1); % plot a clean figure without these points and check indivudual ripples
% clf
% plot(swDiffAll(idx1 | idx2),ripPowerAll(idx1 | idx2),'k.','markersize',1); hold on
% scores = nan(size(ripPowerAll,1),1);
matrix = [swDiffAll ripPowerAll];
z = bsxfun(@rdivide,bsxfun(@minus,matrix,mean(matrix(~bad,:))),std(matrix(~bad,:)));
nothing
0
1
0.7
0.4
0
1
0.7
0
selected = UIInPolygon(swDiffAll,ripPowerAll); % draw polygon to encircle the points you believe are ripples
% (I tend to label an ugly ripple that is still clearly a ripple with a score of 0.2)
scored = find(~isnan(scores));
distances = sqrt(bsxfun(@minus,z(:,1),z(scored,1)').^2 + bsxfun(@minus,z(:,2),z(~isnan(scores),2)').^2);
[~,nearestNeighbour] = min(distances,[],2);
estimated = scores(scored(nearestNeighbour));
selected = estimated>0.1;
idx1 = selected & ~bad; % final ripples
idx2 = ~selected & ~bad; % non-ripples (yet free from noise as well)
display('saving...');
save(fullfile(basepath,'DetectSWR_manual_scoring.mat'),'t','swDiffAll','ripPowerAll','idx1','idx2','bad','selected','scores');
saveas(figure(1),fullfile(basepath,'DetectSWR_manual_scoring.fig'));
dbcont
close all
clear all
X = BatchReturn('OJR_longTraining.batch');
cd(X{1})
cd(X{2})
cd(X{1})
cd(X{2})
load('day7.chanCoords.channelInfo.mat')
load('day7.cell_metrics.cellinfo.mat')
channel_mapping
cd(X{2})
channel_mapping
load('day12.cell_metrics.cellinfo.mat')
channel_mapping
cd(X{2})
channel_mapping
cd(X{3})
cd(X{4})
load('day3.cell_metrics.cellinfo.mat')
channel_mapping
cd(X{5})
channel_mapping
cd(X{6})
close all
channel_mapping
cd(X{6})
channel_mapping
load('day12.ripples.events.mat')
load('day12.MergePoints.events.mat')
sleep = ConsolidateIntervals(MergePoints.timestamps(cellfun(@(x) any(strfind(lower(x),'sleep')),MergePoints.foldernames),:));
pre = InIntervals(ripples,sleep(1,:));
load('day12.ripples.events.mat'); rippleStructure = ripples; ripples = [ripples.timestamps(:,1) ripples.peaks ripples.timestamps(:,2) ripples.RipMax ripples.SwMax];
pre = InIntervals(ripples,sleep(1,:));
post = InIntervals(ripples,sleep(3,:));
post = InIntervals(ripples,sleep(2,:));
dd = diff(ripples(:,[1 3]),[],2);
anovabar(dd,[-pre+post]);
anovabar(dd>0.1,[-pre+post]);
load('day11.ripples.events.mat')
load('day11.MergePoints.events.mat')
rippleStructure = ripples; ripples = [ripples.timestamps(:,1) ripples.peaks ripples.timestamps(:,2) ripples.RipMax ripples.SwMax];
sleep = ConsolidateIntervals(MergePoints.timestamps(cellfun(@(x) any(strfind(lower(x),'sleep')),MergePoints.foldernames),:));
pre = InIntervals(ripples,sleep(1,:));
post = InIntervals(ripples,sleep(2,:));
dd = diff(ripples(:,[1 3]),[],2);
anovabar(dd,[-pre+post]);
2
anovabar(dd>0.1,[-pre+post]);
cd(X{6})
load('day10.ripples.events.mat')
load('day10.MergePoints.events.mat')
rippleStructure = ripples; ripples = [ripples.timestamps(:,1) ripples.peaks ripples.timestamps(:,2) ripples.RipMax ripples.SwMax];
sleep = ConsolidateIntervals(MergePoints.timestamps(cellfun(@(x) any(strfind(lower(x),'sleep')),MergePoints.foldernames),:));
pre = InIntervals(ripples,sleep(1,:));
post = InIntervals(ripples,sleep(2,:));
dd = diff(ripples(:,[1 3]),[],2);
anovabar(dd,[-pre+post]);
anovabar(dd>0.1,[-pre+post]);
cd(X{5})
load('day7.ripples.events.mat')
load('day7.MergePoints.events.mat')
rippleStructure = ripples; ripples = [ripples.timestamps(:,1) ripples.peaks ripples.timestamps(:,2) ripples.RipMax ripples.SwMax];
sleep = ConsolidateIntervals(MergePoints.timestamps(cellfun(@(x) any(strfind(lower(x),'sleep')),MergePoints.foldernames),:));
pre = InIntervals(ripples,sleep(1,:));
post = InIntervals(ripples,sleep(2,:));
dd = diff(ripples(:,[1 3]),[],2);
anovabar(dd,[-pre+post]);
anovabar(dd>0.1,[-pre+post]);
anovabar(dd,[-pre+post]);
raly
open BatchFindRipples
batch = StartBatch(@BatchLoadHPCPFCData,'OJR_longTraining.batch');
54
batch = StartBatch(@BatchLoadHPCPFCData,'OJR_longTraining.batch');
close all
close all hidden
Hide
Hide('off');
Hide('none');
i
for i=1:3
for j=1:3
if i==j
%             PlotColorMap(Smooth(qs{i,1},smooth));
else
subplot(3,3,(i-1)*3+j);
PlotColorMap(Smooth(qs{i,1}./sum(qs{i}(:)) - qs{j,1}./sum(qs{j}(:)),smooth),'x',midbins(:,1),'y',midbins(:,2));
title([names{i} '-' names{j}]);
xlabel('SW depth'); ylabel('ripple duration (ms)');
end
end
end
clims([-1 1]*0.001);
r = [ripples.timestamps(:,1) ripples.peaks ripples.timestamps(:,2) ripples.RipMax ripples.SwMax];
sleep = ConsolidateIntervals(MergePoints.timestamps(cellfun(@(x) any(strfind(lower(x),'sleep')),MergePoints.foldernames),:));
d = diff(r(:,[1 3]),[],2);
nBins = 50; smooth = 3;
values = [r(:,5) d*1000];
[b1] = Bin(values(:,1),nBins); b2 = Bin(values(:,2),nBins); b = [b1 b2];
limits = (min(values(:,1)):range(values(:,1))/nBins:max(values(:,1)))'; limits = [limits(1:end-1) limits(2:end)];
midbins = mean(limits,2);
limits = (min(values(:,2)):range(values(:,2))/nBins:max(values(:,2)))'; limits = [limits(1:end-1) limits(2:end)];
midbins(:,2) = mean(limits,2);
pre = InIntervals(ripples,sleep(1,:));
post = InIntervals(ripples,sleep(2,:));
g = [-pre+post]+2;
r = [ripples.timestamps(:,1) ripples.peaks ripples.timestamps(:,2) ripples.RipMax ripples.SwMax];
sleep = ConsolidateIntervals(MergePoints.timestamps(cellfun(@(x) any(strfind(lower(x),'sleep')),MergePoints.foldernames),:));
d = diff(r(:,[1 3]),[],2);
nBins = 50; smooth = 3;
values = [r(:,5) d*1000];
[b1] = Bin(values(:,1),nBins); b2 = Bin(values(:,2),nBins); b = [b1 b2];
limits = (min(values(:,1)):range(values(:,1))/nBins:max(values(:,1)))'; limits = [limits(1:end-1) limits(2:end)];
limits = (min(values(:,2)):range(values(:,2))/nBins:max(values(:,2)))'; limits = [limits(1:end-1) limits(2:end)];
midbins(:,2) = mean(limits,2);
pre = InIntervals(ripples,sleep(1,:));
post = InIntervals(ripples,sleep(2,:));
limits = (min(values(:,2)):range(values(:,2))/nBins:max(values(:,2)))'; limits = [limits(1:end-1) limits(2:end)];
midbins(:,2) = mean(limits,2);
pre = InIntervals(ripples,sleep(1,:));
r = [ripples.timestamps(:,1) ripples.peaks ripples.timestamps(:,2) ripples.RipMax ripples.SwMax];
sleep = ConsolidateIntervals(MergePoints.timestamps(cellfun(@(x) any(strfind(lower(x),'sleep')),MergePoints.foldernames),:));
d = diff(r(:,[1 3]),[],2);
nBins = 50; smooth = 3;
values = [r(:,5) d*1000];
[b1] = Bin(values(:,1),nBins); b2 = Bin(values(:,2),nBins); b = [b1 b2];
limits = (min(values(:,1)):range(values(:,1))/nBins:max(values(:,1)))'; limits = [limits(1:end-1) limits(2:end)];
midbins = mean(limits,2);
limits = (min(values(:,2)):range(values(:,2))/nBins:max(values(:,2)))'; limits = [limits(1:end-1) limits(2:end)];
midbins(:,2) = mean(limits,2);
pre = InIntervals(r,sleep(1,:));
post = InIntervals(r,sleep(2,:));
g = [-pre+post]+2;
for i=1:3,
qs{i,1} = Accumulate(b(g==i,:),1,'size',max(b));
end
names = {'pre','task','post'};
clf
for i=1:3
for j=1:3
if i==j
%             PlotColorMap(Smooth(qs{i,1},smooth));
else
subplot(3,3,(i-1)*3+j);
PlotColorMap(Smooth(qs{i,1}./sum(qs{i}(:)) - qs{j,1}./sum(qs{j}(:)),smooth),'x',midbins(:,1),'y',midbins(:,2));
title([names{i} '-' names{j}]);
xlabel('SW depth'); ylabel('ripple duration (ms)');
end
end
end
clims([-1 1]*0.001);
% ylabel(clim);
for i=1:3, subplot(3,3,(i-1)*3+i); PlotColorMap(Smooth(qs{i}./sum(qs{i}(:)),smooth),'x',midbins(:,1),'y',midbins(:,2)); if i==1, c = clim;  else, clim(c); end; title([names{i}]); xlabel('SW depth'); ylabel('ripple duration (ms)'); end
; end
K>> dbcomnt
K>> dbconty\
K>> dbcont
dbcont
batch = StartBatch(@BatchLoadHPCPFCData,'OJR_longTraining.batch');
dbcont
batch = StartBatch(@BatchLoadHPCPFCData,'OJR_longTraining.batch');
dbcont
basepath
[parentFolder,basename] = fileparts(basepath);
[~,projectName] = fileparts(parentFolder);
sessionID = [projectName '_' basename];
cd(basepath);
load(fullfile(basepath,[basename '.session.mat']),'session');
try load(fullfile(basepath,[basename '.MergePoints.events.mat']),'MergePoints');
postID = find(cellfun(@(x) ~isempty(strfind(x,'post')), MergePoints.foldernames));
sleep = ConsolidateIntervals(MergePoints.timestamps(cellfun(@(x) any(strfind(lower(x),'sleep')),MergePoints.foldernames),:));
catch
MergePoints = [];
display('No merge points saved. Keeping all ripples (impossible to restrict to post-task sleep');
end
if ~exist(fullfile(basepath,[basename '.lfp']),'file')
'LFP does not exist!'
LFPfromDat(pwd,'outFs',1250,'useGPU',false);
end
% Detect delta waves
if exist(fullfile(basepath,[basename '.deltaWaves.events.mat']),'file')
load(fullfile(basepath,[basename '.deltaWaves.events.mat']),'deltaWaves');
deltas = repmat(deltaWaves.peaks,1,2);
else
'deltas do not exist!'
deltas = [];
%     tic;
% %     lfpstruct = getLFP(114+1);
% %     %     tic; lfpstruct = getLFP(104); HMC is 104, OR18 is 97
% %     lfp = [lfpstruct.timestamps]; lfp(:,2) = lfpstruct.data;
%
%     lfp = GetAyaLFP(114);
%     display(['loaded! @' num2str(toc)]);
%     clean = CleanLFP(lfp);
%     deltas0 = FindDeltaPeaks(clean);
% %     deltas = deltas0(deltas0(:,5)>1.5 & deltas0(:,5)<4.5,:); % these thresholds should be manually refined for each session
%     deltas = deltas0(deltas0(:,5)-deltas0(:,6)>3,:); % these thresholds should be manually refined for each session
%     deltaWaves.timestamps = deltas(:,[1 3]); deltaWaves.peaks = deltas(:,2);  deltaWaves.peakNormedPower = deltas(:,5); deltaWaves.detectorName = ['channel 114(+1), CleanLFP, FindDeltaPeaks, peak-trough>3'];
%     save(fullfile(basepath,[basename '.deltaWaves.events.mat']),'deltaWaves');
end
% Detect ripples
if exist(fullfile(basepath,[basename '.ripples.events.mat']),'file')
load(fullfile(basepath,[basename '.ripples.events.mat']),'ripples');
else
load(fullfile(basepath,[basename '.session.mat']),'session');
rippleChannels = swrChannels('basepath',basepath); % or better yet, pick them on neuroscope (don't forget to add 1! rippleChannels.Ripple_Channel = neuroscope_Channel + 1)
ripples = DetectSWR([rippleChannels.Ripple_Channel rippleChannels.Sharpwave_Channel],'basepath',basepath,'saveMat',true,'forceDetect',true,'check',true);
end
r = ripples.timestamps(:,1:2);
% try
load(fullfile(basepath,[basename '.cell_metrics.cellinfo.mat']),'cell_metrics');
nNeurons = length(cell_metrics.spikes.times);
regions = cell_metrics.brainRegion;
regionNames = unique(regions)
regionCell = zeros(length(regions),1);
for i=1:length(regionNames)
regionCell(strcmp(regions,regionNames{i})) = i;
end
spikesCell = cell_metrics.spikes.times';
PFCindex = find(strcmp(regionNames,'ILA'));
HPCindex = find(strcmp(regionNames,'CA1'));
pfc = []; if ~isempty(PFCindex), pfc = Group(spikesCell{regionCell==PFCindex}); end
hpc = []; if ~isempty(PFCindex), hpc = Group(spikesCell{regionCell==HPCindex}); end
neurons = true;
% catch
%     neurons = false; pfc = zeros(0,2); hpc = zeros(0,2);
% end
spikes = sortrows(Group(spikesCell{:}));
Group({})
Group
batch = StartBatch(@BatchLoadHPCPFCData,'OJR_longTraining.batch');
clear all
batch = StartBatch(@BatchLoadHPCPFCData,'OJR_longTraining.batch');
X = get(batch,'UserData');
open X
sleep = X{i,5};
ripples = X{i,4};
pre = InIntervals(ripples,sleep(1,:));
post = InIntervals(ripples,sleep(2,:));
subplot(2,3,i);
d = diff(ripples(:,[1 3]),[],2);
46
clf
anovabar(dd(:,1),dd(:,2))
m = cellfun(@mean,ds);
plot(m)
plot(m','.')
plot(m','.-')
clear m
plot(m','.-')
plot(m(:,[1 3])','.-')
friedman(m(:,[1 3]))
diff(m(:,[1 3]),[],2)
signrank(m(:,1),m(:,3))
signrank(m(:,1),m(:,2))
anovabar(m);
friedman(m);
stats = friedman(m);
[~,stats] = friedman(m);
[~,~,stats] = friedman(m);
multcompare(stats,1,2)
multcompare(stats);
close all
anovabar(m);
plot(m(:,[1 3])','.-')
plot(zm(:,[1 3])','.-')
plot(m(:,[1 3])','.-')
plot(zm(:,[1 3])','.-')
anovabar(m);
dd = cell2mat([cellfun(@zscore,ds(:,1),'uniformoutput',0) ds(:,2)]);
anovabar(m);
anovabar(dd);
anovabar(dd(:,1),dd(:,2))
Dist(1000,dd,'grouped');
Dist(1000,[dd(:,1) dd(:,2)+2],'grouped');
[h,ht] = Dist(1000,[dd(:,1) dd(:,2)+2],'grouped');
plot(ht,Smooth(h,[2 0]))
plot(ht,Smooth(h,[5 0]))
plot(ht,cumsum(h));
plot(m(:,[1 3])','.-')
plot(m(:,[1 2 3])','.-')
plot(m(:,[1 3])','.-')
plot(zm(:,[1 3])','.-')
kkeyboard
cdN:\OJRproject\OJR43\day11
cd N:\OJRproject\OJR43\day11
load('day11.ripples.events.mat')
load('day11.MergePoints.events.mat')
sleep = X{i,5};
sleep = ConsolidateIntervals(MergePoints.timestamps(cellfun(@(x) any(strfind(lower(x),'sleep')),MergePoints.foldernames),:));
pre = InIntervals(ripples,sleep(1,:));
structure.detectorinfo.detectionintervals = SubtractIntervals([t(1) t(end,1)],excludeInterval);
rippleStructure = ripples; ripples = [ripples.timestamps(:,1) ripples.peaks ripples.timestamps(:,2) ripples.RipMax ripples.SwMax];
pre = InIntervals(ripples,sleep(1,:));
post = InIntervals(ripples,sleep(2,:));
d =diff(ripples(:,[1 3]),[],2);
anovabar(d,-pre+post);
anovabar(d>0.1,-pre+post);
load('day12.ripples.events.mat')
load('day12.MergePoints.events.mat')
sleep = ConsolidateIntervals(MergePoints.timestamps(cellfun(@(x) any(strfind(lower(x),'sleep')),MergePoints.foldernames),:));
structure.detectorinfo.detectionintervals = SubtractIntervals([t(1) t(end,1)],excludeInterval);
rippleStructure = ripples; ripples = [ripples.timestamps(:,1) ripples.peaks ripples.timestamps(:,2) ripples.RipMax ripples.SwMax];
pre = InIntervals(ripples,sleep(1,:));
post = InIntervals(ripples,sleep(2,:));
d =diff(ripples(:,[1 3]),[],2);
anovabar(d,-pre+post);
anovabar(d>0.1,-pre+post);
i=2;
hist(X{i,2}(:,1));
i=2; hist(X{i,2}(:,1));
clf
i=2; hist(X{i,2}(:,1),100);
hist(X{i,2}(:,1),100);
hist(X{i,4}(:,1),100);
cd 'N:\OJRproject\OJR42\day11'
load('day11.SleepState.states.mat', 'SleepState')
PlotIntervals(SleepState.ints.NREMstate,'color','b');
clf
PlotColorMap(h');
hh = h(2,:)./h(1,:);
hh = [h(2,:)./h(1,:); h(3,:)];
plot(ht,hh);
plot(ht,zscore(hh));
plot(ht,zscore(hh,[],2));
plot(ht,nanzscore(hh,[],2));
h0 = h;
h = Smooth(h,[0 2]);
plot(ht,h);
hh = [h(2,:)./h(1,:); h(3,:)];
plot(ht,nanzscore(hh,[],2));
PlotIntervals(SleepState.ints.NREMstate,'color','b');
plot(ht,nanzscore(hh,[],2),'linewidth',2);
figure; plot(ht,nanzscore(h,[],2),'linewidth',2);
legend('all ripples','long ripples','delta');
PlotHVLines(sleep,'v','k--');
PlotIntervals(SleepState.ints.NREMstate,'color',[0.8 0.8 1]);
legend('all ripples','long ripples','delta');
[h,ht] = PETH(deltas,ripples(:,1));
close
PlotColorMap(Shrink(sortby(zscore(h,[],2),d),500,1),'x',ht);
PlotColorMap(Shrink(sortby(zscore(h,[],2),d),10,1),'x',ht);
[h,ht] = PETH(deltas(:,1),ripples(:,1));
PlotColorMap(Shrink(sortby(zscore(h,[],2),d),10,1),'x',ht);
PlotColorMap(Shrink(sortby(zscore(h,[],2),d),500,1),'x',ht);
PlotColorMap(Shrink(sortby(zscore(h,[],2),d),1000,1),'x',ht);
PlotColorMap(Smooth(Shrink(sortby(zscore(h,[],2),d),1000,1),2),'x',ht);
PlotColorMap(Smooth(Shrink(sortby(zscore(h,[],2),d),1000,1),[1 2]),'x',ht);
PlotColorMap(Smooth(Shrink(sortby(zscore(h,[],2),d),500,1),[1 2]),'x',ht);
i
cellfun(@(x) isempty(strfind(x,'train')),MergePoints.foldernames);
MergePoints = X{i,7};
cellfun(@(x) isempty(strfind(x,'train')),MergePoints.foldernames);
training = ConsolidateIntervals(MergePoints.timestamps(cellfun(@(x) ~isempty(strfind(x,'train')),MergePoints.foldernames),:));
65
PETH([re(:,1) nanmean(re(:,2:end),2)],ripples(:,1));
clf
PETH([re(:,1) nanmean(re(:,2:end),2)],ripples(:,1));
training
dur(training)
[templates,correlations,weights] = ActivityTemplatesICA(Restrict(pfc,training,'shift','on'));
i=4;
sleep = X{i,5};
ripples = X{i,4};
pre = InIntervals(ripples,sleep(1,:));
post = InIntervals(ripples,sleep(2,:));
d = diff(ripples,[],2);
MergePoints = X{i,7};
training = ConsolidateIntervals(MergePoints.timestamps(cellfun(@(x) ~isempty(strfind(x,'train')),MergePoints.foldernames),:));
ds{i,1} = d;
ds{i,2} = -pre+post;
%     subplot(2,3,i);
%     anovabar(d,[-pre+post]);
m(i,1) = nanmean(d(pre));
m(i,2) = nanmean(d(~pre & ~post));
m(i,3) = nanmean(d(post));
d = zscore(d);
zm(i,1) = nanmean(d(pre));
zm(i,2) = nanmean(d(~pre & ~post));
zm(i,3) = nanmean(d(post));
pfc = X{i,1};
hpc = X{i,2};
[templates,correlations,weights] = ActivityTemplatesICA(Restrict(pfc,training,'shift','on'));
re = ReactivationStrength(pfc,templates,'step',0.01);
clf
for k=1:size(re,2)-1
[q,qt] = PETH(re(:,[1 1+k]),ripples(:,1));
qq{k,1} = nanmean(q(pre,:));
qq{k,2} = nanmean(q(post,:));
end
kk = 0;
56
45
q = cell2mat(qq(:,1));
PlotColorMap(cell2mat(qq));
PlotColorMap(zscore(cell2mat(qq),[],2));
z = zscore(cell2mat(qq),[],2);
z1 = z(:,1:100);
z2 = z(:,101:end);
clf
semplot(qt,zq);
semplot(qt,z1);
semplot(qt,z2,'r');
for i=1:nSessions
sleep = X{i,5};
ripples = X{i,4};
pre = InIntervals(ripples,sleep(1,:));
post = InIntervals(ripples,sleep(2,:));
d = diff(ripples,[],2);
MergePoints = X{i,7};
training = ConsolidateIntervals(MergePoints.timestamps(cellfun(@(x) ~isempty(strfind(x,'train')),MergePoints.foldernames),:));
ds{i,1} = d;
ds{i,2} = -pre+post;
%     subplot(2,3,i);
%     anovabar(d,[-pre+post]);
m(i,1) = nanmean(d(pre));
m(i,2) = nanmean(d(~pre & ~post));
m(i,3) = nanmean(d(post));
d = zscore(d);
zm(i,1) = nanmean(d(pre));
zm(i,2) = nanmean(d(~pre & ~post));
zm(i,3) = nanmean(d(post));
pfc = X{i,1};
hpc = X{i,2};
[templates,correlations,weights] = ActivityTemplatesICA(Restrict(pfc,training,'shift','on'));
re = ReactivationStrength(pfc,templates,'step',0.01);
[i size(re)]
end
for i=1:nSessions
try
sleep = X{i,5};
ripples = X{i,4};
pre = InIntervals(ripples,sleep(1,:));
post = InIntervals(ripples,sleep(2,:));
d = diff(ripples,[],2);
MergePoints = X{i,7};
training = ConsolidateIntervals(MergePoints.timestamps(cellfun(@(x) ~isempty(strfind(x,'train')),MergePoints.foldernames),:));
ds{i,1} = d;
ds{i,2} = -pre+post;
%     subplot(2,3,i);
%     anovabar(d,[-pre+post]);
m(i,1) = nanmean(d(pre));
m(i,2) = nanmean(d(~pre & ~post));
m(i,3) = nanmean(d(post));
d = zscore(d);
zm(i,1) = nanmean(d(pre));
zm(i,2) = nanmean(d(~pre & ~post));
zm(i,3) = nanmean(d(post));
pfc = X{i,1};
hpc = X{i,2};
[templates,correlations,weights] = ActivityTemplatesICA(Restrict(pfc,training,'shift','on'));
re = ReactivationStrength(pfc,templates,'step',0.01);
[i size(re)]
end
end
kk=0;
int1 = rand(10,1)*100; int1 = [int1 int1+1];
PlotIntervals(int1);
clf
PlotIntervals(int1);
int2 = rand(10,1)*100; int2 = [int2 int2+1];
PlotIntervals(int2);
PlotIntervals(int2,'color','r');
PlotIntervals(int1);
overlap = SubtractIntervals(int1,SubtractIntervals([0 Inf],int2));
PlotColorMap(overlap,'color','b');
PlotIntervals(overlap,'color','b');
intersection = @(int1,int2) SubtractIntervals(int1,SubtractIntervals([0 Inf],int2))
clf
z = zscore(cell2mat(qq),[],2);
z1 = z(:,1:100);
z2 = z(:,101:end);
semplot(qt,z1);
semplot(qt,z2,'r');
figure; subplot(1,2,1); PlotColorMap(z1,'z',qt);
subplot(1,2,1); PlotColorMap(z1,'x',qt);
subplot(1,2,2); PlotColorMap(z2,'x',qt);
clims
clim
clim([-5 5]);
climx([-5 5]);
clims([-5 5]);
clims([-1 1]*7);
subplot(1,1,1); PlotColorMap(z2-z1,'x',qt);
open X
sleep
diff(sleep,[],2)
diff(sleep,[],2)/3600
kk=0;
clear qq
z = zscore(cell2mat(qq),[],2);
z1 = z(:,1:100);
z2 = z(:,101:end);
subplot(1,2,2); PlotColorMap(z2,'x',qt);
subplot(1,2,1); PlotColorMap(z1,'x',qt);
semplot(qt,z2,'r');
semplot(qt,z1,'k');
subplot(2,2,1);
semplot(qt,z1,'k');
subplot(2,2,2);
semplot(qt,z2,'r');
z = (cell2mat(qq));
z1 = z(:,1:100);
z2 = z(:,101:end);
subplot(2,2,3);
semplot(qt,z1,'k');
subplot(2,2,4);
semplot(qt,z2,'r');
subplot(1,2,2); PlotColorMap(z2,'x',qt);
subplot(1,2,1); PlotColorMap(z1,'x',qt);
figure; PETH(deltas(:,1),ripples(:,1));
i=6
cd {'N:\OJRproject\OJR43\day10'}
cd N:\OJRproject\OJR43\day10
load('day10.ripples.events.mat')
r = [ripples.timestamps(:,1) ripples.peaks ripples.timestamps(:,2) ripples.RipMax ripples.SwMax];
open r
anovabar(r(:,5),pre);
close
figure; anovabar(r(:,5),pre);
anovabar(r(:,5),post);
anovabar(r(:,4),post);
anovabar(r(:,5),post);
z = (cell2mat(qq),[],2); z = zBaseline(z,101:200,2);
z = (cell2mat(qq)); z = zBaseline(z,101:200,2);
z1 = z(:,1:100);
z2 = z(:,101:end);
cla
semplot(qt,z2,'r');
semplot(qt,z1,'k');
z = (cell2mat(qq)); z = zBaseline(z,Unfind(101:200,200),2);
z1 = z(:,1:100);
z2 = z(:,101:end);
cla
semplot(qt,z2,'r');
semplot(qt,z1,'k');
5
2
z = (cell2mat(qq)); z = zBaseline(z,Unfind(101:200,200),2);
open qq
raly
z = (cell2mat(qq(:,1:2))); z = zBaseline(z,Unfind(101:200,200),2);
subplot(2,2,3); cla;
semplot(qt,z1,'k');
subplot(2,2,4); cla;
semplot(qt,z2,'r');
for i=1:nSessions
sleep = X{i,5};
ripples = X{i,4};
sleep = sleep(:,1); sleep(:,2) = sleep(:,1)+3600;
pre = InIntervals(ripples,sleep(1,:));
post = InIntervals(ripples,sleep(2,:));
d = diff(ripples,[],2);
MergePoints = X{i,7};
%     try
load(fullfile(basepath,[basename '.SleepState.states.mat']));
sws = SleepState.ints.NREMstate;
end
i
preSleep = SubtractIntervals(sleep(1,:), SubtractIntervals([0 Inf],sws));
postSleep = SubtractIntervals(sleep(2,:), SubtractIntervals([0 Inf],sws));
Unshift(3600,postSleep)
sleep = X{i,5};
postSleep = SubtractIntervals(sleep(2:end,:), SubtractIntervals([0 Inf],sws));
Unshift(3600,postSleep)
preSleep(preSleep(:,1)>Unshift(3600,preSleep),:) = [];
dur(preSleep)
Unshift(3600,preSleep)
preSleep = SubtractIntervals(sleep(1,:), SubtractIntervals([0 Inf],sws));
Unshift(3600,preSleep)
dur(SubtractIntervals(sleep(1,:), SubtractIntervals([0 Inf],sws)))
preSleep(preSleep(:,1)>Unshift(3600,preSleep),:) = [];
dur(preSleep)
dur(postSleep)
preSleep = SubtractIntervals(sleep(1,:), SubtractIntervals([0 Inf],sws));
postSleep = SubtractIntervals(sleep(2:end,:), SubtractIntervals([0 Inf],sws));
% Restrict to 1h of sws only:
limit = Unshift(3600,preSleep); preSleep(preSleep(:,1)>limit,:) = []; if preSleep(end,2)>limit, preSleep(end,2) = limit; end
limit = Unshift(3600,postSleep); postSleep(postSleep(:,1)>limit,:) = []; if postSleep(end,2)>limit, postSleep(end,2) = limit; end
dur(preSleep)
dur(postSleep)
45
size(pre)
size(post)
size(r)
size(d)
z = (cell2mat(qq(:,1:2))); z = zBaseline(z,Unfind(101:200,200),2);
open qq
z = (cell2mat(qq(:,1:2))); z = zBaseline(z,Unfind(101:200,200),2);
subplot(1,2,2); PlotColorMap(z2,'x',qt);
subplot(1,2,1); PlotColorMap(z1,'x',qt);
subplot(2,2,2); PlotColorMap(z2,'x',qt);
subplot(2,2,1); PlotColorMap(z1,'x',qt);
subplot(2,2,4); semplot(qt,z2,'r');
subplot(2,2,3); semplot(qt,z1,'k');
45
length(qt)
2
i=1
i
sleep = X{i,5};
ripples = X{i,4};
d = diff(ripples,[],2);
MergePoints = X{i,7};
load(fullfile(basepath,[basename '.SleepState.states.mat']));
sws = SleepState.ints.NREMstate;
preSleep = SubtractIntervals(sleep(1,:), SubtractIntervals([0 Inf],sws));
postSleep = SubtractIntervals(sleep(2:end,:), SubtractIntervals([0 Inf],sws));
% Restrict to 1h of sws only:
limit = Unshift(3600,preSleep); preSleep(preSleep(:,1)>limit,:) = []; if preSleep(end,2)>limit, preSleep(end,2) = limit; end
limit = Unshift(3600,postSleep); postSleep(postSleep(:,1)>limit,:) = []; if postSleep(end,2)>limit, postSleep(end,2) = limit; end
training = ConsolidateIntervals(MergePoints.timestamps(cellfun(@(x) ~isempty(strfind(x,'train')),MergePoints.foldernames),:));
pre = InIntervals(ripples,preSleep);
post = InIntervals(ripples,postSleep);
ds{i,1} = d;
ds{i,2} = -pre+post;
%     subplot(2,3,i);
%     anovabar(d,[-pre+post]);
m(i,1) = nanmean(d(pre));
m(i,2) = nanmean(d(~pre & ~post));
m(i,3) = nanmean(d(post));
d = zscore(d);
zm(i,1) = nanmean(d(pre));
zm(i,2) = nanmean(d(~pre & ~post));
zm(i,3) = nanmean(d(post));
pfc = X{i,1};
hpc = X{i,2};
pfc = hpc;
isempty(pfc)
open X
cd N:\OJRproject\OJR42\day11
i=3;
pfc = X{i,1};
hpc = X{i,2};
max(pfc(:,2))
max(hpc(:,2))
38+3
cd N:\Praveen\uLED\uLED1\day7
batch2 = StartBatch(@BatchLoadHPCPFCData,'OJR_longTraining.batch');
q = get(batch,'UserData');
open q
q = get(batch2,'UserData');
open q
load('day7.cell_metrics.cellinfo.mat')
batch2 = StartBatch(@BatchLoadHPCPFCData,'OJR_longTraining.batch');
X = get(batch2,'UserData');
open DX
open X
4
load('digitalIn.events.mat')
figure; PlotIntervals(sws);
PlotHVLines(sleep,'v','k--','linewidth',2);
clickers = Group(digitalIn.timestampsOn{:});
open clickers
RasterPlot(clickers);
PlotHVLines(sleep,'v','k--','linewidth',2);
clickers = cellfun(@(x) x', digitalIn.ints);
clickers = cellfun(@(x) x', digitalIn.ints, 'UniformOutput', false);
clickers = cellfun(@(x) x', digitalIn.ints, 'UniformOutput', false)';
for i=1:5, ylim([-1 1]*0.5+i); PlotIntervals(clickers{i}); end
ylim([0 6]);
i
for i=2:5, ylim([-1 1]*0.5+i); PlotIntervals(clickers{i}); end
ylim([0 6]);
clf
for i=2:5, ylim([-1 1]*0.5+i); PlotIntervals(clickers{i}); end
ylim([0 6]);
for i=2:5, ylim([-1 1]*0.5+i); PlotIntervals(clickers{i},'color','r'); end
ylim([0 6]);
PlotHVLines(sleep,'v','k--','linewidth',2);
training = (MergePoints.timestamps(cellfun(@(x) ~isempty(strfind(x,'train')),MergePoints.foldernames),:));
PlotHVLines(training,'v','r--','linewidth',2);
PlotHVLines(training,'v','y--','linewidth',2);
338/60
294/60
sec2min(294)
460+22-41
sec2min(441)
dur(training(end,:))
sec2min(dur(training(end,:)))
RasterPlot(Group(digitalIn.timestampsOn{:}));
batch = StartBatch(@BatchLoadHPCPFCData,'OJR_longTraining.batch');
X2 = get(batch,'UserData');
open X2
digitalInp = getDigitalIn('all','fs',session.extracellular.sr);
load('day7.session.mat')
digitalInp = getDigitalIn('all','fs',session.extracellular.sr);
batch = StartBatch(@BatchLoadHPCPFCData,'OJR_longTraining.batch');
i=1; figure;  [pfc,hpc,deltas,ripples,sleep,session,MergePoints,clickers,basepath] = X{i,:}; for i=2:5, ylim([-1 1]*0.5+i); PlotIntervals(clickers{i},'color','r'); end
X = get(batch,'UserData');
i=1; figure;  [pfc,hpc,deltas,ripples,sleep,session,MergePoints,clickers,basepath] = X{i,:}; for i=2:5, ylim([-1 1]*0.5+i); PlotIntervals(clickers{i},'color','r'); end
training = ConsolidateIntervals(MergePoints.timestamps(cellfun(@(x) ~isempty(strfind(x,'train')),MergePoints.foldernames),:));
PlotHVLines(training,'v','y--','linewidth',2);
PlotHVLines(training,'v','y--','linewidth',2); ylim([0 6]);
PlotHVLines(training,'v','y--','linewidth',2); ylim([0 12]);
PlotHVLines(training,'v','y--','linewidth',2); ylim([0 16]);
i=1; figure;  [pfc,hpc,deltas,ripples,sleep,session,MergePoints,clickers,basepath] = X{i,:}; for i=2:length(clickers), ylim([-1 1]*0.5+i); PlotIntervals(clickers{i},'color','r'); end
i
ylim([0 6]);
ylim([0 16]);
ylim([0 16]);
training
diff(training,[],2)
diff(training,[],2)/60
MergePoints.foldernames'
cd basepath
cd(basepath)
dur(clickers{2})
dur(clickers{3})
close all
exploration = ConsolidateIntervals(clickers{2};clickers{3});
exploration = ConsolidateIntervals([clickers{2};clickers{3}]);
exploration = ConsolidateIntervals(sortrows([clickers{2};clickers{3}]));
exploration = sortrows([clickers{2};clickers{3}]); exploration(any(isnan(exploration),2),:) = []; exploration = ConsolidateIntervals(exploration);
dur(exploration)
i
N:\Praveen\uLED\uLED1\day7\
cd N:\Praveen\uLED\uLED1\day7\
2
23
clim
clims([-5 5])
clims([-1 1]*3)
clims([-1 1]*5)
clim
clims([-1 1]*5)
size(z)
close all
figure('name','exploration,100ms,100ms');
z = (cell2mat(qq(:,1:2)));
z = zscore(z,[],2);
% z = zBaseline(z,Unfind((length(qt)+1):(length(qt)*2),(length(qt)*2)),2);
% z([17 23],:) = [];
z1 = z(:,1:length(qt)); z2 = z(:,(length(qt)+1):end);
smooth = [0 2];
clf
subplot(2,3,1); PlotColorMap(Smooth(z1,smooth),'x',qt); PlotHVLines(0,'v','k--','linewidth',2);
subplot(2,3,2); PlotColorMap(Smooth(z2,smooth),'x',qt); PlotHVLines(0,'v','k--','linewidth',2);
subplot(2,3,3); PlotColorMap(Smooth(z2-z1,smooth),'x',qt); PlotHVLines(0,'v','k--','linewidth',2);
clims([-1 1]*5);
subplot(2,3,4); semplot(qt,z1,'k',smooth(end));  semplot(qt,z2,'r',smooth(end)); y = ylim; PlotHVLines(0,'v','k--','linewidth',2); PlotHVLines(0,'h','k--','linewidth',2); ylim(y);
subplot(2,3,5); semplot(qt,z2,'r',smooth(end)); ylim(y); PlotHVLines(0,'v','k--','linewidth',2); PlotHVLines(0,'h','k--','linewidth',2);ylim(y);
subplot(2,3,6); semplot(qt,z2-z1,'m',smooth(end)); PlotHVLines(0,'v','k--','linewidth',2); PlotHVLines(0,'h','k--','linewidth',2);
raly
lines = cumsum(Accumulate(cell2mat(qq(:,end))))+0.5;
lines
lines = cumsum(Accumulate(cell2mat(qq(:,end))));
lines
lines = Accumulate(cell2mat(qq(:,end)));
lines = Accumulate(cell2mat(qq(:,end)))
nComponents = Accumulate(cell2mat(qq(:,end)))
open X
[parentFolder,basename] = fileparts(basepath)
[parentFolder,basename] = fileparts(basepath);
[~,projectName] = fileparts(parentFolder);
sessionID = [projectName '_' basename];
sessionID
for i=1:nSessions
i
[pfc,hpc,deltas,ripples,sleep,session,MergePoints,clickers,basepath] = X{i,:};
d = diff(ripples,[],2);
[parentFolder,basename] = fileparts(basepath);
[~,projectName] = fileparts(parentFolder);
projectNames{i,1} = projectName;
end
for i=1:nSessions
%     i
[pfc,hpc,deltas,ripples,sleep,session,MergePoints,clickers,basepath] = X{i,:};
d = diff(ripples,[],2);
[parentFolder,basename] = fileparts(basepath);
[~,projectName] = fileparts(parentFolder);
sessionID = [projectName '_' basename];
sessionIDs{i,1} = sessionID;
end
nComponents = Accumulate(cell2mat(qq(:,end)));
sessions = sessionIDs{nComponents>0}
sessions = sessionIDs(nComponents>0)
sessions = sessionIDs(nComponents>0)
lines = cumsum(nComponents(nComponents>0));
lines
lines = cumsum(nComponents(nComponents>0))+0.5;
PlotHVLines(lines,'h','w--');
lines = cumsum(nComponents(nComponents>0))+0.5;
set(gca,'ytick',mean([[0;lines(1:end-1)] lines],2));
set(gca,'ytick',mean([[0;lines(1:end-1)] lines],2),'yticklabel',sessions);
addSlash = @(x) strfind(x,'_');
addSlash(sessions{1})
sessions{1}
x = OJR42_day12;
x = 'OJR42_day12';
strfind(x,'_')
x(strfind(x,'_'))
x(1:strfind(x,'_'))
x(sort([(1:length(x)) strfind(x,'_')]))
script_OJR_analyse_long_training
title(x(sort([(1:length(x)) strfind(x,'_')])))
addSlash = @(x) x(sort([(1:length(x)) strfind(x,'_')]));
addSlash = @(x) x(sort([(1:length(x)) strfind(x,'_')]));
sessions = cellfun(addSlash,sessionIDs(nComponents>0));
45
close all
45
close all
figure
clf
matrices = {z1,z2,z2-z1};
titles = {[name ' in pre-task sleep'],[name ' in pre-task sleep'],'difference'};
colors = {'k','r','m'};
j=1
subplot(2,3,j); PlotColorMap(Smooth(matrices{j},smooth),'x',qt); PlotHVLines(0,'v','k--','linewidth',2);
PlotHVLines(lines,'h','w--');
set(gca,'ytick',mean([[0;lines(1:end-1)] lines],2),'yticklabel',sessions);
for j=1:3
subplot(2,3,j); PlotColorMap(Smooth(matrices{j},smooth),'x',qt); PlotHVLines(0,'v','k--','linewidth',2);
PlotHVLines(lines,'h','w--');
set(gca,'ytick',mean([[0;lines(1:end-1)] lines],2),'yticklabel',sessions);
title(titles{j});
subplot(2,3,3+j); semplot(qt,matrices{j},colors{j},smooth(end));
if j==1, semplot(qt,matrices{j+1},colors{j+1},smooth(end)); y=ylim; end
ylim(y); PlotHVLines(0,'v','k--','linewidth',2); PlotHVLines(0,'h','k--','linewidth',2); ylim(y);
ylabel('mean +/ s.e.m.');
end
sort([(1:length(x)) strfind(x,'_')])
strrep(x,'_','\_')
name
strrep(name,' ','-')
strrep(name,' ','_')
SaveFig(fullfile('M:\home\raly\results\PFC\reactivation',strrep([name ' around ripples'],' ','_')))
meanResponse = nanmean(q(:,[find(InIntervals(qt,[-0.1 0.2])) find(InIntervals(qt,[-0.1 0.2]))+length(qt)]),2);
find(meanResponse==0)
find(meanResponse<0)
2
open sessions
i=3
[pfc,hpc,deltas,ripples,sleep,session,MergePoints,clickers,basepath] = X{i,:};
d = diff(ripples,[],2);
[parentFolder,basename] = fileparts(basepath);
[~,projectName] = fileparts(parentFolder);
sessionID = [projectName '_' basename];
sessionIDs{i,1} = sessionID;
sws = [0 Inf];
try
basename = basenameFromBasepath(basepath);
load(fullfile(basepath,[basename '.SleepState.states.mat']));
sws = SleepState.ints.NREMstate;
catch
'No sleep scoring'
end
2
[pfc,hpc,deltas,ripples,sleep,session,MergePoints,clickers,basepath] = X{i,:};
d = diff(ripples,[],2);
[parentFolder,basename] = fileparts(basepath);
[~,projectName] = fileparts(parentFolder);
sessionID = [projectName '_' basename];
sessionIDs{i,1} = sessionID;
sws = [0 Inf];
basename = basenameFromBasepath(basepath);
load(fullfile(basepath,[basename '.SleepState.states.mat']));
sws = SleepState.ints.NREMstate;
whos SleepState
124106233/1024
124106233/1024/1024
preSleep = SubtractIntervals(sleep(1,:), SubtractIntervals([0 Inf],sws));
postSleep = SubtractIntervals(sleep(2:end,:), SubtractIntervals([0 Inf],sws));
% Restrict to 1h of sws only:
limit = Unshift(3600,preSleep); preSleep(preSleep(:,1)>limit,:) = []; if preSleep(end,2)>limit, preSleep(end,2) = limit; end
limit = Unshift(3600,postSleep); postSleep(postSleep(:,1)>limit,:) = []; if postSleep(end,2)>limit, postSleep(end,2) = limit; end
training = ConsolidateIntervals(MergePoints.timestamps(cellfun(@(x) ~isempty(strfind(x,'train')),MergePoints.foldernames),:));
exploration = sortrows([clickers{2};clickers{3}]); exploration(any(isnan(exploration),2),:) = []; exploration = ConsolidateIntervals(exploration);
exploration = Restrict(exploration,training);
pre = InIntervals(ripples,preSleep);
post = InIntervals(ripples,postSleep);
ds{i,1} = d;
ds{i,2} = -pre+post;
%     subplot(2,3,i);
%     anovabar(d,[-pre+post]);
m(i,1) = nanmean(d(pre));
m(i,2) = nanmean(d(~pre & ~post));
m(i,3) = nanmean(d(post));
d = zscore(d);
zm(i,1) = nanmean(d(pre));
zm(i,2) = nanmean(d(~pre & ~post));
zm(i,3) = nanmean(d(post));
pfc = X{i,1};
hpc = X{i,2};
structure
structure = 2
isempty(pfc)
restriction=1
if restriction==1, intervals = training; else, intervals = exploration; end
[templates,correlations,weights] = ActivityTemplatesICA(Restrict(pfc,intervals,'shift','on'),'binsize',0.1);
re = ReactivationStrength(pfc,templates,'step',0.01,'binSize',0.1);
open post
open postSleep
in = InIntervals(re(:,1),postSleep(24,:));
plot(re(in,1),re(in,2:end));
xlim(9319+[0 5]);
xlim(9319+[0 5]+5);
xlim(9319+[0 5]);
xlim(9319+[0 5]+10);
xlim(9319+[0 5]+15);
xlim(9319+[0 5]+20);
xlim(9319+[0 5]+25);
xlim(9319+[0 5]+30);
xlim(9319+[0 5]+35);
figure; PETH(re(:,[1 2]),ripples(:,2));
PETH(re(:,[1 2]),ripples(:,2),'durations',[-1 1]*5);
PETH(re(:,[1 3]),ripples(:,2),'durations',[-1 1]*5);
PETH(hpc(:,1),ripples(:,2),'durations',[-1 1]*5);
PETH(pfc(:,1),ripples(:,2),'durations',[-1 1]*5);
open X
i=3; PETH(X{i,1}(:,1),X{i,4}(:,2),'durations',[-1 1]*5);
i=4; PETH(X{i,1}(:,1),X{i,4}(:,2),'durations',[-1 1]*5);
i=5; PETH(X{i,1}(:,1),X{i,4}(:,2),'durations',[-1 1]*5);
i=5; PETH(X{i,1}(:,1),X{i,4}(:,2),'durations',[-1 1]*5,'nBins',501);
i=6; PETH(X{i,1}(:,1),X{i,4}(:,2),'durations',[-1 1]*5,'nBins',501);
i=4; PETH(X{i,1}(:,1),X{i,4}(:,2),'durations',[-1 1]*5,'nBins',501);
i=3; PETH(X{i,1}(:,1),X{i,4}(:,2),'durations',[-1 1]*5,'nBins',501);
i=4; PETH(X{i,1}(:,1),X{i,4}(:,2),'durations',[-1 1]*10,'nBins',501);
i=3; PETH(X{i,1}(:,1),X{i,4}(:,2),'durations',[-1 1]*10,'nBins',501);
i=5; PETH(X{i,1}(:,1),X{i,4}(:,2),'durations',[-1 1]*10,'nBins',501);
i=6; PETH(X{i,1}(:,1),X{i,4}(:,2),'durations',[-1 1]*10,'nBins',501);
i=6; PETH(X{i,2}(:,1),X{i,4}(:,2),'durations',[-1 1]*10,'nBins',501);
i=6; PETH(X{i,2}(:,1),X{i,4}(:,2),'durations',[-1 1]*10,'nBins',101);
i=3; PETH(X{i,2}(:,1),X{i,4}(:,2),'durations',[-1 1]*10,'nBins',101);
i=6; PETH(X{i,2}(:,1),X{i,4}(:,2),'durations',[-1 1]*10,'nBins',101);
130/70
i=3; PETH(X{i,2}(:,1),X{i,4}(:,2),'durations',[-1 1]*10,'nBins',101);
13/8
i=4; PETH(X{i,2}(:,1),X{i,4}(:,2),'durations',[-1 1]*10,'nBins',101);
130/90
xlim(9319+[0 5]+40);
xlim(9319+[0 5]+45);
xlim(9319+[0 5]+50);
figure; PETH(re(:,[1 2]),deltas(:,2));
PETH(re(:,[1 2]),deltas(:,2),'durations',[-1 1]*2,'nBins',501);
PETH(re(:,[1 3]),deltas(:,2),'durations',[-1 1]*2,'nBins',501);
PETH(re(:,[1 4]),deltas(:,2),'durations',[-1 1]*2,'nBins',501);
PETH(re(:,[1 5]),deltas(:,2),'durations',[-1 1]*2,'nBins',501);
PETH(re(:,[1 6]),deltas(:,2),'durations',[-1 1]*2,'nBins',501);
PETH(re(:,[1 5]),deltas(:,2),'durations',[-1 1]*2,'nBins',501);
PETH(re(:,[1 5]),ripples(:,2),'durations',[-1 1]*2,'nBins',501);
PETH(ripples(:,2),deltas(:,2),'durations',[-1 1]*2,'nBins',501);
PETH(ripples(:,2),deltas(:,2),'durations',[-1 1]*1,'nBins',501);
i=3; PETH(X{i,3}(:,2),X{i,2}(:,2),'durations',[-1 1]*1,'nBins',501);
i=3; PETH(X{i,3}(:,1),X{i,2}(:,1),'durations',[-1 1]*1,'nBins',501);
i=3; PETH(X{i,3}(:,1),X{i,4}(:,1),'durations',[-1 1]*1,'nBins',501);
i=4; PETH(X{i,3}(:,1),X{i,4}(:,1),'durations',[-1 1]*1,'nBins',501);
i=2; PETH(X{i,3}(:,1),X{i,4}(:,1),'durations',[-1 1]*1,'nBins',501);
i=3; PETH(X{i,3}(:,1),X{i,4}(:,1),'durations',[-1 1]*1,'nBins',501);
i=2; PETH(X{i,3}(:,1),X{i,4}(:,1),'durations',[-1 1]*1,'nBins',501);
i=6; PETH(X{i,3}(:,1),X{i,4}(:,1),'durations',[-1 1]*1,'nBins',501);
hold all
i=2; PETH(X{i,3}(:,1),X{i,4}(:,1),'durations',[-1 1]*1,'nBins',501);
i=3; PETH(X{i,3}(:,1),X{i,4}(:,1),'durations',[-1 1]*1,'nBins',501);
xlim(9374+[0 5]);
xlim(9374+[0 5]+5);
session
figure; PETH(re(:,[1 5],deltas(:,2)));
figure; PETH(re(:,[1 5]),deltas(:,2));
PETH(re(:,[1 5]),X{i,3}(:,2));
raly
save('temp.mat','re');
figure; PETH(ripples(:,1),deltas(:,2));
clear all
batch = StartBatch(@BatchLoadHPCPFCData,'OJR_longTraining.batch');
X = get(batch,'UserData');
clear all
script_OJR_analyse_long_training
i
2
45
5
set(gcf,'position',[1 1 1920 1024])
67
set(gcf,'position',[1 1 1920 1024])
%%
name = [name ' ICs'];
dd = cell2mat([cellfun(@zscore,ds(:,1),'uniformoutput',0) ds(:,2)]);
nComponents = Accumulate(cell2mat(qq(:,end)));
addSlash = @(x) strrep(x,'_','\_');
sessions = cellfun(addSlash,sessionIDs(nComponents>0),'UniformOutput',0);
lines = cumsum(nComponents(nComponents>0))+0.5;
q = (cell2mat(qq(:,1:2)));
meanResponse = nanmean(q(:,[find(InIntervals(qt,[-0.1 0.2])) find(InIntervals(qt,[-0.1 0.2]))+length(qt)]),2);
%         flip = meanResponse<0; q(flip,:) = -q(flip,:);
% z = zBaseline(z,Unfind((length(qt)+1):(length(qt)*2),(length(qt)*2)),2);
z = zscore(q,[],2);
z1 = z(:,1:length(qt)); z2 = z(:,(length(qt)+1):end);
smooth = [0 2];
% figure('name',name);
clf
set(gcf,'position',[1 1 1920 1024])
matrices = {z1,z2,z2-z1};
titles = {[name ' in pre-task sleep'],[name ' in post-task sleep'],'(post-pre) difference'};
colors = {'k','r','m'};
for j=1:3
subplot(2,3,j); PlotColorMap(Smooth(matrices{j},smooth),'x',qt); PlotHVLines(0,'v','k--','linewidth',2);
PlotHVLines(lines,'h','w--','linewidth',2);
set(gca,'ytick',mean([[0;lines(1:end-1)] lines],2),'yticklabel',sessions);
title(titles{j});
subplot(2,3,3+j); semplot(qt,matrices{j},colors{j},smooth(end));
if j==1, semplot(qt,matrices{j+1},colors{j+1},smooth(end)); y=ylim; end
ylim(y); PlotHVLines(0,'v','k--','linewidth',2); PlotHVLines(0,'h','k--','linewidth',2); ylim(y);
ylabel('mean +/ s.e.m. (z-units)');
end
name = [name ' ICs'];
dd = cell2mat([cellfun(@zscore,ds(:,1),'uniformoutput',0) ds(:,2)]);
nComponents = Accumulate(cell2mat(qq(:,end)));
addSlash = @(x) strrep(x,'_','\_');
sessions = cellfun(addSlash,sessionIDs(nComponents>0),'UniformOutput',0);
lines = cumsum(nComponents(nComponents>0))+0.5;
q = (cell2mat(qq(:,1:2)));
meanResponse = nanmean(q(:,[find(InIntervals(qt,[-0.1 0.2])) find(InIntervals(qt,[-0.1 0.2]))+length(qt)]),2);
flip = meanResponse<0; q(flip,:) = -q(flip,:);
% z = zBaseline(z,Unfind((length(qt)+1):(length(qt)*2),(length(qt)*2)),2);
z = zscore(q,[],2);
z1 = z(:,1:length(qt)); z2 = z(:,(length(qt)+1):end);
smooth = [0 2];
% figure('name',name);
clf
set(gcf,'position',[1 1 1920 1024])
matrices = {z1,z2,z2-z1};
titles = {[name ' in pre-task sleep'],[name ' in post-task sleep'],'(post-pre) difference'};
colors = {'k','r','m'};
for j=1:3
subplot(2,3,j); PlotColorMap(Smooth(matrices{j},smooth),'x',qt); PlotHVLines(0,'v','k--','linewidth',2);
PlotHVLines(lines,'h','w--','linewidth',2);
set(gca,'ytick',mean([[0;lines(1:end-1)] lines],2),'yticklabel',sessions);
title(titles{j});
subplot(2,3,3+j); semplot(qt,matrices{j},colors{j},smooth(end));
if j==1, semplot(qt,matrices{j+1},colors{j+1},smooth(end)); y=ylim; end
ylim(y); PlotHVLines(0,'v','k--','linewidth',2); PlotHVLines(0,'h','k--','linewidth',2); ylim(y);
ylabel('mean +/ s.e.m. (z-units)');
end
name = [name ' ICs'];
dd = cell2mat([cellfun(@zscore,ds(:,1),'uniformoutput',0) ds(:,2)]);
nComponents = Accumulate(cell2mat(qq(:,end)));
addSlash = @(x) strrep(x,'_','\_');
sessions = cellfun(addSlash,sessionIDs(nComponents>0),'UniformOutput',0);
lines = cumsum(nComponents(nComponents>0))+0.5;
q = (cell2mat(qq(:,1:2)));
meanResponse = nanmean(q(:,[find(InIntervals(qt,[-0.1 0.2])) find(InIntervals(qt,[-0.1 0.2]))+length(qt)]),2);
%         flip = meanResponse<0; q(flip,:) = -q(flip,:);
% z = zBaseline(z,Unfind((length(qt)+1):(length(qt)*2),(length(qt)*2)),2);
z = zscore(q,[],2);
z1 = z(:,1:length(qt)); z2 = z(:,(length(qt)+1):end);
smooth = [0 2];
% figure('name',name);
clf
set(gcf,'position',[1 1 1920 1024])
matrices = {z1,z2,z2-z1};
titles = {[name ' in pre-task sleep'],[name ' in post-task sleep'],'(post-pre) difference'};
colors = {'k','r','m'};
for j=1:3
subplot(2,3,j); PlotColorMap(Smooth(matrices{j},smooth),'x',qt); PlotHVLines(0,'v','k--','linewidth',2);
PlotHVLines(lines,'h','w--','linewidth',2);
set(gca,'ytick',mean([[0;lines(1:end-1)] lines],2),'yticklabel',sessions);
title(titles{j});
subplot(2,3,3+j); semplot(qt,matrices{j},colors{j},smooth(end));
if j==1, semplot(qt,matrices{j+1},colors{j+1},smooth(end)); y=ylim; end
ylim(y); PlotHVLines(0,'v','k--','linewidth',2); PlotHVLines(0,'h','k--','linewidth',2); ylim(y);
ylabel('mean +/ s.e.m. (z-units)');
end
flip = meanResponse<0;
open flip
find(flip)
ylim([15 25]);
ylim([18 21]);
ylim([19 21]);
clim
figure; plot(q(20,:));
i=3
open X
[pfc,hpc,deltas,ripples,sleep,session,MergePoints,clickers,basepath] = X{i,:};
d = diff(ripples,[],2);
[parentFolder,basename] = fileparts(basepath);
[~,projectName] = fileparts(parentFolder);
sessionID = [projectName '_' basename];
sessionIDs{i,1} = sessionID;
sws = [0 Inf];
try
basename = basenameFromBasepath(basepath);
load(fullfile(basepath,[basename '.SleepState.states.mat']));
sws = SleepState.ints.NREMstate;
catch
'No sleep scoring'
end
preSleep = SubtractIntervals(sleep(1,:), SubtractIntervals([0 Inf],sws));
postSleep = SubtractIntervals(sleep(2:end,:), SubtractIntervals([0 Inf],sws));
% Restrict to 1h of sws only:
limit = Unshift(3600,preSleep); preSleep(preSleep(:,1)>limit,:) = []; if preSleep(end,2)>limit, preSleep(end,2) = limit; end
limit = Unshift(3600,postSleep); postSleep(postSleep(:,1)>limit,:) = []; if postSleep(end,2)>limit, postSleep(end,2) = limit; end
training = ConsolidateIntervals(MergePoints.timestamps(cellfun(@(x) ~isempty(strfind(x,'train')),MergePoints.foldernames),:));
exploration = sortrows([clickers{2};clickers{3}]); exploration(any(isnan(exploration),2),:) = []; exploration = ConsolidateIntervals(exploration);
exploration = Restrict(exploration,training);
%     exploration = SubtractIntervals(training,exploration);
pre = InIntervals(ripples,preSleep);
post = InIntervals(ripples,postSleep);
ds{i,1} = d;
ds{i,2} = -pre+post;
%     subplot(2,3,i);
%     anovabar(d,[-pre+post]);
m(i,1) = nanmean(d(pre));
m(i,2) = nanmean(d(~pre & ~post));
m(i,3) = nanmean(d(post));
d = zscore(d);
zm(i,1) = nanmean(d(pre));
zm(i,2) = nanmean(d(~pre & ~post));
zm(i,3) = nanmean(d(post));
pfc = X{i,1};
hpc = X{i,2};
if structure==1,
pfc = hpc;
end
if isempty(pfc), continue; end
if restriction==1, intervals = training; else, intervals = exploration; end
[templates,correlations,weights] = ActivityTemplatesICA(Restrict(pfc,intervals,'shift','on'),'binsize',0.1);
re = ReactivationStrength(pfc,templates,'step',0.01,'binSize',0.1);
% Restrict to 1h of sws only:
limit = Unshift(3600,preSleep); preSleep(preSleep(:,1)>limit,:) = []; if preSleep(end,2)>limit, preSleep(end,2) = limit; end
limit = Unshift(3600,postSleep); postSleep(postSleep(:,1)>limit,:) = []; if postSleep(end,2)>limit, postSleep(end,2) = limit; end
training = ConsolidateIntervals(MergePoints.timestamps(cellfun(@(x) ~isempty(strfind(x,'train')),MergePoints.foldernames),:));
exploration = sortrows([clickers{2};clickers{3}]); exploration(any(isnan(exploration),2),:) = []; exploration = ConsolidateIntervals(exploration);
exploration = Restrict(exploration,training);
%     exploration = SubtractIntervals(training,exploration);
pre = InIntervals(ripples,preSleep);
post = InIntervals(ripples,postSleep);
ds{i,1} = d;
ds{i,2} = -pre+post;
%     subplot(2,3,i);
%     anovabar(d,[-pre+post]);
m(i,1) = nanmean(d(pre));
m(i,2) = nanmean(d(~pre & ~post));
m(i,3) = nanmean(d(post));
d = zscore(d);
zm(i,1) = nanmean(d(pre));
zm(i,2) = nanmean(d(~pre & ~post));
zm(i,3) = nanmean(d(post));
pfc = X{i,1};
hpc = X{i,2};
if structure==1,
pfc = hpc;
end
if restriction==1, intervals = training; else, intervals = exploration; end
[templates,correlations,weights] = ActivityTemplatesICA(Restrict(pfc,intervals,'shift','on'),'binsize',0.1);
re = ReactivationStrength(pfc,templates,'step',0.01,'binSize',0.1);
[h,ht] = PETH(re(:,[1 8]),ripples(:,1));
figure; PlotColorMap(h);
semplot(ht,h);
clf
semplot(ht,h(pre,:);
semplot(ht,h(pre,:));
find(pre,10)
ripples(50,1)
preSleep = SubtractIntervals(sleep(1,:), SubtractIntervals([0 Inf],sws));
dur(preSleep)
open MergePoints
clear all
close all
batch = StartBatch(@BatchLoadHPCPFCData,'OJR_longTraining.batch');
X = get(batch,'UserData');
2
5
clim
clims([-1 1]*2);
name = [name ' ICs'];
dd = cell2mat([cellfun(@zscore,ds(:,1),'uniformoutput',0) ds(:,2)]);
nComponents = Accumulate(cell2mat(qq(:,end)));
addSlash = @(x) strrep(x,'_','\_');
sessions = cellfun(addSlash,sessionIDs(nComponents>0),'UniformOutput',0);
lines = cumsum(nComponents(nComponents>0))+0.5;
q = (cell2mat(qq(:,1:2)));
%         z = zBaseline(q,Unfind((length(qt)+1):(length(qt)*2),(length(qt)*2)),2);
z = zscore(q,[],2);
meanResponse = nanmean(z(:,[find(InIntervals(qt,[-0.1 0.2])) find(InIntervals(qt,[-0.1 0.1]))+length(qt)]),2);
flip = meanResponse<0; z(flip,:) = -z(flip,:);
z1 = z(:,1:length(qt)); z2 = z(:,(length(qt)+1):end);
smooth = [0 2];
% figure('name',name);
clf
set(gcf,'position',[1 1 1920 1024])
matrices = {z1,z2,z2-z1};
titles = {[name ' in pre-task sleep'],[name ' in post-task sleep'],'(post-pre) difference'};
colors = {'k','r','m'};
for j=1:3
subplot(2,3,j); PlotColorMap(Smooth(matrices{j},smooth),'x',qt); PlotHVLines(0,'v','k--','linewidth',2);
PlotHVLines(lines,'h','w--','linewidth',2);
set(gca,'ytick',mean([[0;lines(1:end-1)] lines],2),'yticklabel',sessions);
title(titles{j});
subplot(2,3,3+j); semplot(qt,matrices{j},colors{j},smooth(end));
if j==1, semplot(qt,matrices{j+1},colors{j+1},smooth(end)); y=ylim; end
ylim(y); PlotHVLines(0,'v','k--','linewidth',2); PlotHVLines(0,'h','k--','linewidth',2); ylim(y);
ylabel('mean +/ s.e.m. (z-units)');
end
meanResponse = nanmean(z(:,[find(InIntervals(qt,[-0.1 0.1])) find(InIntervals(qt,[-0.1 0.1]))+length(qt)]),2);
name = [name ' ICs'];
dd = cell2mat([cellfun(@zscore,ds(:,1),'uniformoutput',0) ds(:,2)]);
nComponents = Accumulate(cell2mat(qq(:,end)));
addSlash = @(x) strrep(x,'_','\_');
sessions = cellfun(addSlash,sessionIDs(nComponents>0),'UniformOutput',0);
lines = cumsum(nComponents(nComponents>0))+0.5;
q = (cell2mat(qq(:,1:2)));
%         z = zBaseline(q,Unfind((length(qt)+1):(length(qt)*2),(length(qt)*2)),2);
z = zscore(q,[],2);
meanResponse = nanmean(z(:,[find(InIntervals(qt,[-0.1 0.1])) find(InIntervals(qt,[-0.1 0.1]))+length(qt)]),2);
flip = meanResponse<0; z(flip,:) = -z(flip,:);
z1 = z(:,1:length(qt)); z2 = z(:,(length(qt)+1):end);
smooth = [0 2];
% figure('name',name);
clf
set(gcf,'position',[1 1 1920 1024])
matrices = {z1,z2,z2-z1};
titles = {[name ' in pre-task sleep'],[name ' in post-task sleep'],'(post-pre) difference'};
colors = {'k','r','m'};
for j=1:3
subplot(2,3,j); PlotColorMap(Smooth(matrices{j},smooth),'x',qt); PlotHVLines(0,'v','k--','linewidth',2);
PlotHVLines(lines,'h','w--','linewidth',2);
set(gca,'ytick',mean([[0;lines(1:end-1)] lines],2),'yticklabel',sessions);
title(titles{j});
subplot(2,3,3+j); semplot(qt,matrices{j},colors{j},smooth(end));
if j==1, semplot(qt,matrices{j+1},colors{j+1},smooth(end)); y=ylim; end
ylim(y); PlotHVLines(0,'v','k--','linewidth',2); PlotHVLines(0,'h','k--','linewidth',2); ylim(y);
ylabel('mean +/ s.e.m. (z-units)');
end
z = zBaseline(q,Unfind((length(qt)+1):(length(qt)*2),(length(qt)*2)),2);
%         z = zscore(q,[],2);
meanResponse = nanmean(z(:,[find(InIntervals(qt,[-0.1 0.1])) find(InIntervals(qt,[-0.1 0.1]))+length(qt)]),2);
flip = meanResponse<0; z(flip,:) = -z(flip,:);
z1 = z(:,1:length(qt)); z2 = z(:,(length(qt)+1):end);
smooth = [0 2];
% figure('name',name);
clf
set(gcf,'position',[1 1 1920 1024])
matrices = {z1,z2,z2-z1};
titles = {[name ' in pre-task sleep'],[name ' in post-task sleep'],'(post-pre) difference'};
colors = {'k','r','m'};
for j=1:3
subplot(2,3,j); PlotColorMap(Smooth(matrices{j},smooth),'x',qt); PlotHVLines(0,'v','k--','linewidth',2);
PlotHVLines(lines,'h','w--','linewidth',2);
set(gca,'ytick',mean([[0;lines(1:end-1)] lines],2),'yticklabel',sessions);
title(titles{j});
subplot(2,3,3+j); semplot(qt,matrices{j},colors{j},smooth(end));
if j==1, semplot(qt,matrices{j+1},colors{j+1},smooth(end)); y=ylim; end
ylim(y); PlotHVLines(0,'v','k--','linewidth',2); PlotHVLines(0,'h','k--','linewidth',2); ylim(y);
ylabel('mean +/ s.e.m. (z-units)');
end
z = zBaseline(q,Unfind((length(qt)+1):(length(qt)*2),(length(qt)*2)),2);
z = zscore(q,[],2);
meanResponse = nanmean(z(:,[find(InIntervals(qt,[-0.1 0.1])) find(InIntervals(qt,[-0.1 0.1]))+length(qt)]),2);
flip = meanResponse<0; z(flip,:) = -z(flip,:);
z1 = z(:,1:length(qt)); z2 = z(:,(length(qt)+1):end);
smooth = [0 2];
% figure('name',name);
clf
set(gcf,'position',[1 1 1920 1024])
matrices = {z1,z2,z2-z1};
titles = {[name ' in pre-task sleep'],[name ' in post-task sleep'],'(post-pre) difference'};
colors = {'k','r','m'};
for j=1:3
subplot(2,3,j); PlotColorMap(Smooth(matrices{j},smooth),'x',qt); PlotHVLines(0,'v','k--','linewidth',2);
PlotHVLines(lines,'h','w--','linewidth',2);
set(gca,'ytick',mean([[0;lines(1:end-1)] lines],2),'yticklabel',sessions);
title(titles{j});
subplot(2,3,3+j); semplot(qt,matrices{j},colors{j},smooth(end));
if j==1, semplot(qt,matrices{j+1},colors{j+1},smooth(end)); y=ylim; end
ylim(y); PlotHVLines(0,'v','k--','linewidth',2); PlotHVLines(0,'h','k--','linewidth',2); ylim(y);
ylabel('mean +/ s.e.m. (z-units)');
end
clims
clim
SaveFig(fullfile('M:\home\raly\results\PFC\reactivation',strrep(['Extended_training ' name ' around ripples'],' ','_')))
54
open qq
eopn qqs
indices = cell2mat(qq(:,end));
open qqs
this = cellfun(@(x) tiedrank(x)/length(x),ds(:,1),'UniformOutput',0);
this = [ds,cellfun(@(x) tiedrank(x)/length(x),ds(:,1),'UniformOutput',0)];
open this
that = cell2mat(dd(indices,1))
open dd
dd = cell2mat(this(indices,:));
q = cell2mat(qqs);
PlotColorMap(q);
PlotColorMap(Shrink(q,100,1))
PlotColorMap(Shrink(q,1000,1))
qqz = cellfun(@(x) (x-(nanmean(x(:))))./nanstd(x(:)),qqs);
qqz = cellfun(@(x) (x-(nanmean(x(:))))./nanstd(x(:)),qqs,'UniformOutput',false);
q = cell2mat(qqz);
PlotColorMap(Shrink(q,1000,1))
PlotColorMap(Shrink(sortby(q,dd(:,2)),1000,1))
PlotColorMap(Shrink(sortby(q,dd(:,2)),10000,1))
PlotColorMap(Shrink(sortby(q,dd(:,2)+rand(size(dd(:,1)))*0.001),10000,1))
Portion(dd(:,2))
Portion(dd(:,2)==1)
Portion(dd(:,2)==-1)
PlotColorMap(Shrink(sortby(q,dd(:,2)+rand(size(dd(:,1)))*0.001),50000,1))
PlotColorMap(Shrink(sortby(q,dd(:,2)+rand(size(dd(:,1)))*0.001),1000,1))
PlotColorMap(Shrink(sortby(q,dd(:,2)+rand(size(dd(:,1)))*0.001),2000,1))
anovabar(dd(:,3),dd(:,2))
PlotColorMap(Shrink(sortby(q,dd(:,2)+rand(size(dd(:,1)))*0.001),2000,1))
ok = dd(:,2)~=0;
PlotColorMap(Shrink(sortby(q,dd(ok,2)+rand(size(dd(ok,1)))*0.001),2000,1))
PlotColorMap(Shrink(sortby(dd,dd(ok,2)+rand(size(dd(ok,1)))*0.001),2000,1))
PlotColorMap(Shrink(sortby(q,dd(ok,2)+rand(size(dd(ok,1)))*0.001),2000,1))
PlotColorMap(Shrink(sortby(q,dd(ok,2)+rand(size(dd(ok,1)))*0.001),2000,1),'x',ht)
PlotColorMap(Shrink(sortby(q,dd(ok,2)+rand(size(dd(ok,1)))*0.001),2000,1),'x',qt)
semplot(q(dd(:,2)==1,:));
semplot(q(dd(:,2)==-1,:));
clf
PlotColorMap(Shrink(sortby(q,dd(ok,2)+rand(size(dd(ok,1)))*0.001),2000,1),'x',qt)
PlotColorMap(Shrink(sortby(q,dd(ok,2)*100+rand(size(dd(ok,1)))*0.001),2000,1),'x',qt)
PlotColorMap(Shrink(sortby(q,dd(ok,2)*100+rand(size(dd(ok,1)))*0.001),10000,1),'x',qt)
PlotColorMap(Shrink(sortby(q,dd(ok,2)<1),10000,1),'x',qt)
size(ok)
mmax(dd(ok,2))
mmax(dd(~ok,2))
semplot(q(dd(ok,2)<0,:))
semplot(q(dd(ok,2)>0,:))
semplot(q(dd(:,2)==-1,:));
clf
PlotColorMap(Shrink(sortby(q(ok,:),dd(ok,2)<1),10000,1),'x',qt)
PlotColorMap(Shrink(sortby(q(ok,:),dd(ok,2)),10000,1),'x',qt)
PlotColorMap(Shrink(sortby(q(ok,:),dd(ok,2)),100,1),'x',qt)
PlotColorMap(Shrink(sortby(q(ok,:),dd(ok,2)),1000,1),'x',qt)
PlotColorMap(Shrink(sortby(q(ok,:),dd(ok,2)),800,1),'x',qt)
PlotColorMap(Smooth(Shrink(sortby(q(ok,:),dd(ok,2)),800,1),1),'x',qt)
PlotColorMap(Smooth(Shrink(sortby(q(ok,:),dd(ok,1)),800,1),1),'x',qt)
PlotColorMap(Smooth(Shrink(sortby(q(ok,:),dd(ok,3)),800,1),1),'x',qt)
PlotColorMap(Smooth(Shrink(sortby(q(ok,:),dd(ok,3)),800,1),0),'x',qt)
PlotColorMap(Smooth(Shrink(sortby(q(ok,:),dd(ok,3)),8000,1),0),'x',qt)
PlotColorMap(Smooth(Shrink(sortby(q(ok,:),dd(ok,1)),8000,1),0),'x',qt)
PlotColorMap(Smooth(Shrink(sortby(q(ok,:),dd(ok,1)),8000,5),0),'x',qt)
PlotColorMap(Smooth(Shrink(sortby(q(ok,:),dd(ok,1)),8000,1),[0 5]),'x',qt)
PlotColorMap(Smooth(Shrink(sortby(q(ok,:),dd(ok,1)),8000,1),[0 10]),'x',qt)
PlotColorMap(Smooth(Shrink(sortby(q(ok,:),dd(ok,3)),8000,1),[0 10]),'x',qt)
PlotColorMap(Smooth(Shrink(sortby(q(ok,:),dd(ok,3)),4000,1),[0 10]),'x',qt)
PlotColorMap(Smooth(Shrink(sortby(q(ok,:),dd(ok,3)),4000,1),[2 10]),'x',qt)
PlotColorMap(Smooth(Shrink(sortby(q(ok,:),dd(ok,3)),4000,1),[1 20]),'x',qt)
PlotColorMap(Smooth(Shrink(sortby(dd(ok,:),dd(ok,3)),4000,1),[1 20]),'x',qt)
clim
PlotColorMap(Smooth(Shrink(sortby(dd(ok,:),dd(ok,3)),4000,1),[1 20]),'x',qt)
PlotColorMap(Smooth(Shrink(sortby(q(ok,:),dd(ok,3)),4000,1),[1 20]),'x',qt)
open qq
open qqz
id = repelem(indices,cellfun(@(x) size(x,1),qqz));
size(id)
ok = dd(:,2)~=0 & id==3;
PlotColorMap(Smooth(Shrink(sortby(q(ok,:),dd(ok,3)),4000,1),[1 20]),'x',qt)
PlotColorMap(Smooth(Shrink(sortby(q(ok,:),dd(ok,3)),4000,1),[0 20]),'x',qt)
PlotColorMap(Smooth(Shrink(sortby(q(ok,:),dd(ok,3)),4000,1),[0 5]),'x',qt)
ok = dd(:,2)~=0 & id==4;
PlotColorMap(Smooth(Shrink(sortby(q(ok,:),dd(ok,3)),4000,1),[0 5]),'x',qt)
PlotColorMap(Smooth(Shrink(sortby(q(ok,:),dd(ok,3)),sum(ok)/3,1),[0 5]),'x',qt)
PlotColorMap(Smooth(Shrink(sortby(q(ok,:),dd(ok,3)),ceil(sum(ok)/3),1),[0 5]),'x',qt)
PlotColorMap(Smooth(Shrink(sortby(q(ok,:),dd(ok,3)),ceil(sum(ok)/10),1),[0 5]),'x',qt)
ok = id==4;PlotColorMap(Smooth(Shrink(sortby(q(ok,:),dd(ok,3)),ceil(sum(ok)/10),1),[0 5]),'x',qt)
ok = id==4 & dd(:,2)~=0;PlotColorMap(Smooth(Shrink(sortby(q(ok,:),dd(ok,3)),ceil(sum(ok)/10),1),[0 5]),'x',qt)
ok = id==5 & dd(:,2)~=0;PlotColorMap(Smooth(Shrink(sortby(q(ok,:),dd(ok,3)),ceil(sum(ok)/10),1),[0 5]),'x',qt)
ok = id==6 & dd(:,2)~=0;PlotColorMap(Smooth(Shrink(sortby(q(ok,:),dd(ok,3)),ceil(sum(ok)/10),1),[0 5]),'x',qt)
ok = id==6;
anovabar(dd(ok,1),dd(ok,2));
boxplot(dd(ok,1),dd(ok,2));
anovabar(dd(ok,1),dd(ok,2));
anovabar(dd(ok,1),dd(ok,2)); ylim(quantile(dd(ok,1),[0.1 0.9]))
ok = id==5;
anovabar(dd(ok,1),dd(ok,2)); ylim(quantile(dd(ok,1),[0.1 0.9]))
open X
i=3; PETH(X{i,3}(:,1),X{i,4}(:,1),'durations',[-1 1]*1,'nBins',501);
i=5; PETH(X{i,3}(:,1),X{i,4}(:,1),'durations',[-1 1]*1,'nBins',501);
i=6; PETH(X{i,3}(:,1),X{i,4}(:,1),'durations',[-1 1]*1,'nBins',501);
i=5; PETH(X{i,3}(:,1),X{i,4}(:,1),'durations',[-1 1]*1,'nBins',501);
i=6; PETH(X{i,3}(:,1),X{i,4}(:,1),'durations',[-1 1]*1,'nBins',501);
i=3; PETH(X{i,3}(:,1),X{i,4}(:,1),'durations',[-1 1]*1,'nBins',501);
i=6; PETH(X{i,3}(:,1),X{i,4}(:,1),'durations',[-1 1]*1,'nBins',501);
cd(X{i,end})
load('day10.animal.behavior.mat', 'behavior')
figure; plot(behavior.position.x,behavior.position.y);
load('day10.ripples.events.mat', 'ripples')
clf
i=6; hist(diff(X{i,4},[],2),100);
i=3; hist(diff(X{i,4},[],2),100);
i=6; hist(diff(X{i,4},[],2),100);
g = Group(diff(X{3,4},[],2),diff(X{6,4},[],2));
anovabar(g(:,1),g(:,2))
PlotColorMap(Smooth(Shrink(sortby(q(ok,:),dd(ok,3)),ceil(sum(ok)/10),1),[0 5]),'x',qt)
ok = id==6 & dd(:,2)~=0;PlotColorMap(Smooth(Shrink(sortby(q(ok,:),dd(ok,3)),ceil(sum(ok)/10),1),[0 5]),'x',qt)
ok = id==5 & dd(:,2)~=0;PlotColorMap(Smooth(Shrink(sortby(q(ok,:),dd(ok,3)),ceil(sum(ok)/10),1),[0 5]),'x',qt)
ok = id==4 & dd(:,2)~=0;PlotColorMap(Smooth(Shrink(sortby(q(ok,:),dd(ok,3)),ceil(sum(ok)/10),1),[0 5]),'x',qt)
ok = id==3 & dd(:,2)~=0;PlotColorMap(Smooth(Shrink(sortby(q(ok,:),dd(ok,3)),ceil(sum(ok)/10),1),[0 5]),'x',qt)
ok = id==1 & dd(:,2)~=0;PlotColorMap(Smooth(Shrink(sortby(q(ok,:),dd(ok,3)),ceil(sum(ok)/10),1),[0 5]),'x',qt)
cd N:\OJRproject\OJR43\day7
unique(id)
for i=1:4
subplot(2,2,i)
ok = id==i+2 & dd(:,2)~=0;PlotColorMap(Smooth(Shrink(sortby(q(ok,:),dd(ok,3)),ceil(sum(ok)/10),1),[0 5]),'x',qt)
end
clim
clims
clim
clims([-1 1]*0.1);
clims([-1 1]*0.05);
clims([-1 1]*0.08);
clims([-1 1]*0.09);
clims([-1 1]*0.07);
[pfc,hpc,deltas,ripples,sleep,session,MergePoints,clickers,basepath] = X{i,:};
batch = StartBatch(@BatchLoadHPCPFCData,'OJR_longTraining.batch');
X = get(batch,'UserData');
raly
open ds
this = [ds,cellfun(@(x) tiedrank(x)/length(x),ds(:,1),'UniformOutput',0)];
dd = cell2mat(this(indices,:));
size(qqz)
q = cell2mat(qqz);
size(q)
open dd
dd = cell2mat(this(indices,:));
this = [dsi,cellfun(@(x) tiedrank(x)/length(x),ds(:,1),'UniformOutput',0)];
i
i=6
[pfc,hpc,deltas,ripples,sleep,session,MergePoints,clickers,basepath] = X{i,:};
d = diff(ripples(:,[1 3]),[],2);
[parentFolder,basename] = fileparts(basepath);
[~,projectName] = fileparts(parentFolder);
sessionID = [projectName '_' basename];
sessionIDs{i,1} = sessionID;
sws = [0 Inf];
basename = basenameFromBasepath(basepath);
load(fullfile(basepath,[basename '.SleepState.states.mat']));
sws = SleepState.ints.NREMstate;
preSleep = SubtractIntervals(sleep(1,:), SubtractIntervals([0 Inf],sws));
postSleep = SubtractIntervals(sleep(2:end,:), SubtractIntervals([0 Inf],sws));
% Restrict to 1h of sws only:
limit = Unshift(3600,preSleep); preSleep(preSleep(:,1)>limit,:) = []; if preSleep(end,2)>limit, preSleep(end,2) = limit; end
limit = Unshift(3600,postSleep); postSleep(postSleep(:,1)>limit,:) = []; if postSleep(end,2)>limit, postSleep(end,2) = limit; end
training = ConsolidateIntervals(MergePoints.timestamps(cellfun(@(x) ~isempty(strfind(x,'train')),MergePoints.foldernames),:));
exploration = sortrows([clickers{2};clickers{3}]); exploration(any(isnan(exploration),2),:) = []; exploration = ConsolidateIntervals(exploration);
exploration = Restrict(exploration,training);
%     exploration = SubtractIntervals(training,exploration);
pre = InIntervals(ripples,preSleep);
post = InIntervals(ripples,postSleep);
ds{i,1} = d;
ds{i,2} = -pre+post;
ds{i,3} = ripples(:,4:5);
%     subplot(2,3,i);
this = [dsi,cellfun(@(x) tiedrank(x)/length(x),ds(:,1),'UniformOutput',0)];
this = [ds,cellfun(@(x) tiedrank(x)/length(x),ds(:,1),'UniformOutput',0)];
dd = cell2mat(this(indices,:));
open dd
size(q)
q = cell2mat(qqz);
size(q)
qqz = cellfun(@(x) (x-(nanmean(x(:))))./nanstd(x(:)),qqs,'UniformOutput',false);
open re
isempty(pfc)
nSessions
open X
batch = StartBatch(@BatchLoadHPCPFCData,'OJR_longTraining.batch');
X = get(batch,'UserData');
X = get(batch,'UserData');
batch = StartBatch(@BatchLoadHPCPFCData,'OJR_longTraining.batch');
r = [ripples.timestamps(:,1) ripples.peaks ripples.timestamps(:,2) ripples.amplitude]; r(:,5) = nan;
dbcont
X = get(batch,'UserData');
435
qqz = cellfun(@(x) (x-(nanmean(x(:))))./nanstd(x(:)),qqs,'UniformOutput',false);
indices = cell2mat(qq(:,end));
this = [ds,cellfun(@(x) tiedrank(x)/length(x),ds(:,1),'UniformOutput',0)];
dd = cell2mat(this(indices,:));
q = cell2mat(qqz);
size(q)
size(dd)
for i=1:4
subplot(2,2,i)
ok = id==i+2 & dd(:,2)~=0;PlotColorMap(Smooth(Shrink(sortby(q(ok,:),dd(ok,3)),ceil(sum(ok)/10),1),[0 5]),'x',qt)
end
id = repelem(indices,cellfun(@(x) size(x,1),qqz));
for i=1:4
subplot(2,2,i)
ok = id==i+2 & dd(:,2)~=0;PlotColorMap(Smooth(Shrink(sortby(q(ok,:),dd(ok,3)),ceil(sum(ok)/10),1),[0 5]),'x',qt)
end
clims([-1 1]*0.07);
open dd
figure; ok = id>0 & dd(:,2)~=0;PlotColorMap(Smooth(Shrink(sortby(q(ok,:),dd(ok,3)),ceil(sum(ok)/10),1),[0 5]),'x',qt)
ok = id>0 & dd(:,2)~=0;PlotColorMap(Smooth(Shrink(sortby(q(ok,:),dd(ok,4)),ceil(sum(ok)/10),1),[0 5]),'x',qt)
ok = id>0 & dd(:,2)~=0;PlotColorMap(Smooth(Shrink(sortby(q(ok,:),dd(ok,5)),ceil(sum(ok)/10),1),[0 5]),'x',qt)
for i=1:4
subplot(2,2,i)
ok = id==i+2 & dd(:,2)~=0;PlotColorMap(Smooth(Shrink(sortby(q(ok,:),dd(ok,end)),ceil(sum(ok)/10),1),[0 5]),'x',qt)
end
for i=1:4
subplot(2,2,i)
ok = id==i+2 & dd(:,2)~=0;PlotColorMap(Smooth(Shrink(sortby(q(ok,:),dd(ok,end)),ceil(sum(ok)/10),1),[0 5]),'x',qt)
title(addSlash(sessionIDs{i+i}));
end
re = nanmean(q(:,qt>0 & qt<0.5);
re = nanmean(q(:,qt>0 & qt<0.5));
open re
re = nanmean(q(:,qt>0 & qt<0.5),2);
figure; plot(dd(:,3),re);
plot(dd(:,3),re,'.');
plot(dd(:,end),re,'.');
nancorr(dd(:,end),re)
[c,p] = nancorr(dd(:,end),re)
[c,p] = nancorr(dd(:,end),re,'type','spearman')
plot(dd(:,end),tiedrank(re),'.');
DensityMap(dd(:,end),tiedrank(re))
DensityMap(dd(:,end),tiedrank(re),'smooth',0)
DensityMap(dd(:,end),tiedrank(re),'smooth',0,'nBins',500)
DensityMap(dd(:,end),tiedrank(re),'smooth',0,'nBins',[50 500])
DensityMap(tiedrank(dd(:,end)+rand(size(dd(:,1)))*0.0001),tiedrank(re),'smooth',0,'nBins',[50 500])
DensityMap(tiedrank(dd(:,end)+rand(size(dd(:,1)))*0.0001),tiedrank(re),'smooth',0,'nBins',[500 500])
DensityMap(tiedrank(dd(:,end)+rand(size(dd(:,1)))*0.0001),tiedrank(re+rand(size(re(:,1)))*0.0001),'smooth',0,'nBins',[500 500])
DensityMap(tiedrank(dd(:,end)+rand(size(dd(:,1)))*0.0001),tiedrank(re+rand(size(re(:,1)))*0.0001),'smooth',5,'nBins',[500 500])
nancorr(tiedrank(dd(:,end)+rand(size(dd(:,1)))*0.0001),tiedrank(re+rand(size(re(:,1)))*0.0001))
[c,p] = nancorr(tiedrank(dd(:,end)+rand(size(dd(:,1)))*0.0001),tiedrank(re+rand(size(re(:,1)))*0.0001))
DensityMap(tiedrank(dd(:,end)+rand(size(dd(:,1)))*0.0001),tiedrank(re+rand(size(re(:,1)))*0.0001),'smooth',10,'nBins',[500 500])
DensityMap(tiedrank(dd(:,end)+rand(size(dd(:,1)))*0.0001),tiedrank(re+rand(size(re(:,1)))*0.0001),'smooth',1,'nBins',[500 500])
DensityMap(tiedrank(dd(:,end)+rand(size(dd(:,1)))*0.0001),tiedrank(re+rand(size(re(:,1)))*0.0001),'smooth',1,'nBins',[50 50])
DensityMap(tiedrank(dd(:,end)+rand(size(dd(:,1)))*0.0001),tiedrank(re+rand(size(re(:,1)))*0.0001),'smooth',0,'nBins',[50 50])
DensityMap(tiedrank(dd(:,end)+rand(size(dd(:,1)))*0.0001),tiedrank(re+rand(size(re(:,1)))*0.0001),'smooth',0,'nBins',[1 1]*20)
re = nanmean(q(:,qt>0 & qt<0.2),2);
DensityMap(tiedrank(dd(:,end)+rand(size(dd(:,1)))*0.0001),tiedrank(re+rand(size(re(:,1)))*0.0001),'smooth',0,'nBins',[1 1]*20)
[c,p] = nancorr(tiedrank(dd(:,end)+rand(size(dd(:,1)))*0.0001),tiedrank(re+rand(size(re(:,1)))*0.0001))
for i=1:4
subplot(2,2,i)
ok = id==i+2 & dd(:,2)~=0;PlotColorMap(Smooth(Shrink(sortby(q(ok,:),dd(ok,end)),ceil(sum(ok)/10),1),[0 5]),'x',qt)
title(addSlash(sessionIDs{i+i}));
end
clf
set(gcf,'position',[1 1 1920 1024])
matrices = {z1,z2,z2-z1};
titles = {[name ' in pre-task sleep'],[name ' in post-task sleep'],'(post-pre) difference'};
colors = {'k','r','m'};
for j=1:3
subplot(2,3,j); PlotColorMap(Smooth(matrices{j},smooth),'x',qt); PlotHVLines(0,'v','k--','linewidth',2);
PlotHVLines(lines,'h','w--','linewidth',2);
set(gca,'ytick',mean([[0;lines(1:end-1)] lines],2),'yticklabel',sessions);
title(titles{j});
subplot(2,3,3+j); semplot(qt,matrices{j},colors{j},smooth(end));
if j==1, semplot(qt,matrices{j+1},colors{j+1},smooth(end)); y=ylim; end
ylim(y); PlotHVLines(0,'v','k--','linewidth',2); PlotHVLines(0,'h','k--','linewidth',2); ylim(y);
ylabel('mean +/ s.e.m. (z-units)');
end
clims
45
5+
45
open nothing
5
25
name = [name ' ICs'];
dd = cell2mat([cellfun(@zscore,ds(:,1),'uniformoutput',0) ds(:,2)]);
nComponents = Accumulate(cell2mat(qq(:,end)));
addSlash = @(x) strrep(x,'_','\_');
sessions = cellfun(addSlash,sessionIDs(nComponents>0),'UniformOutput',0);
lines = cumsum(nComponents(nComponents>0))+0.5;
q = (cell2mat(qq(:,1:2)));
z = zBaseline(q,Unfind((length(qt)+1):(length(qt)*2),(length(qt)*2)),2);
%         z = zscore(q,[],2);
meanResponse = nanmean(z(:,[find(InIntervals(qt,[-0.1 0.1])) find(InIntervals(qt,[-0.1 0.1]))+length(qt)]),2);
flip = meanResponse<0; z(flip,:) = -z(flip,:);
z1 = z(:,1:length(qt)); z2 = z(:,(length(qt)+1):end);
smooth = [0 2];
% figure('name',name);
clf
set(gcf,'position',[1 1 1920 1024])
matrices = {z1,z2,z2-z1};
titles = {[name ' in pre-task sleep'],[name ' in post-task sleep'],'(post-pre) difference'};
colors = {'k','r','m'};
for j=1:3
subplot(2,3,j); PlotColorMap(Smooth(matrices{j},smooth),'x',qt); PlotHVLines(0,'v','k--','linewidth',2);
PlotHVLines(lines,'h','w--','linewidth',2);
set(gca,'ytick',mean([[0;lines(1:end-1)] lines],2),'yticklabel',sessions);
title(titles{j});
subplot(2,3,3+j); semplot(qt,matrices{j},colors{j},smooth(end));
if j==1, semplot(qt,matrices{j+1},colors{j+1},smooth(end)); y=ylim; end
ylim(y); PlotHVLines(0,'v','k--','linewidth',2); PlotHVLines(0,'h','k--','linewidth',2); ylim(y);
ylabel('mean +/ s.e.m. (z-units)');
end
clims
drawnow
raly
4
clear all
cd N:\OJRproject\OJR33\day7
load('digitalIn.events.mat')
load('day7.pulses.events.mat')
pulses.timestamps(pulses.eventGroupID==2,:);
stimHPC = pulses.timestamps(pulses.eventGroupID==2,:);
stimPFC = pulses.timestamps(pulses.eventGroupID==1,:);
[h,ht] = mPETH(stimPFC(:,1),stimHPC(:,1));
PETH(stimPFC(:,1),stimHPC(:,1));
PETH(stimPFC(:,1),stimHPC(:,1),'durations',[-1 1]*2);
PETH(stimPFC(:,1),stimHPC(:,1),'durations',[-1 1]*5);
open mPTH
open mPETH
which fastPETH
exist('Sync')
exist('Sync',2)
exist('Sync','fun')
exist('Sync','var')
exist('Sync','class')
mPETH(stimPFC(:,1),stimHPC(:,1),'durations',[-1 1]*5);
mPETH(stimPFC(:,1),stimHPC(:,1),'durations',[-1 1]*5,'nBins',20);
PETH(stimPFC(:,1),stimHPC(:,1),'durations',[-1 1]*5,'nBins',20);
PETH(stimPFC(:,1),stimHPC(:,1),'durations',[-1 1]*5,'nBins',10);
findmax(mPETH(stimPFC(:,1),stimHPC(:,1),'durations',[-1 1]*5,'nBins',10))
findmax(mPETH(stimPFC(:,1),stimHPC(:,1),'durations',[-1 1]*5,'nBins',10))>5
if findmax(mPETH(stimPFC(:,1),stimHPC(:,1),'durations',[-1 1]*5,'nBins',10))<5 % the PFC stimulation should follow the HPC one
warning('HPC stimulation appears to follow PFC stimulation... Check if your stim channels are correct')
end
raly
batch = StartBatch(@BatchLoadHPCPFCStimData,'OJR_longTraining.batch');
batch = StartBatch(@BatchLoadHPCPFCStimData,'OJR_shortTraining_optoRipples_delayedPFC.batch');
X = get(batch,'UserData');
85
5
clf
set(gcf,'position',[1 1 1920 1024])
matrices = {z1,z2,z2-z1};
titles = {[name ' in pre-task sleep'],[name ' in post-task sleep'],'(post-pre) difference'};
colors = {'k','r','m'};
for j=1:3
subplot(2,3,j); PlotColorMap(Smooth(matrices{j},smooth),'x',qt); PlotHVLines(0,'v','k--','linewidth',2);
PlotHVLines(lines,'h','w--','linewidth',2);
set(gca,'ytick',mean([[0;lines(1:end-1)] lines],2),'yticklabel',sessions);
title(titles{j});
subplot(2,3,3+j); semplot(qt,matrices{j},colors{j},smooth(end));
if j==1, semplot(qt,matrices{j+1},colors{j+1},smooth(end)); y=ylim; end
ylim(y); PlotHVLines(0,'v','k--','linewidth',2); PlotHVLines(0,'h','k--','linewidth',2); ylim(y);
ylabel('mean +/ s.e.m. (z-units)');
end
clims
close all
Hide('none');
clf
set(gcf,'position',[1 1 1920 1024])
matrices = {z1,z2,z2-z1};
titles = {[name ' in pre-task sleep'],[name ' in post-task sleep'],'(post-pre) difference'};
colors = {'k','r','m'};
for j=1:3
subplot(2,3,j); PlotColorMap(Smooth(matrices{j},smooth),'x',qt); PlotHVLines(0,'v','k--','linewidth',2);
PlotHVLines(lines,'h','w--','linewidth',2);
set(gca,'ytick',mean([[0;lines(1:end-1)] lines],2),'yticklabel',sessions);
title(titles{j});
subplot(2,3,3+j); semplot(qt,matrices{j},colors{j},smooth(end));
if j==1, semplot(qt,matrices{j+1},colors{j+1},smooth(end)); y=ylim; end
ylim(y); PlotHVLines(0,'v','k--','linewidth',2); PlotHVLines(0,'h','k--','linewidth',2); ylim(y);
ylabel('mean +/ s.e.m. (z-units)');
end
clims
5
i=2;
[pfc,hpc,deltas,ripples,~,sleep,session,MergePoints,clickers,basepath] = X{i,:};
open pfc
open X
i=3;
[pfc,hpc,deltas,ripples,~,sleep,session,MergePoints,clickers,basepath] = X{i,:};
figure; PETH(pfc(:,1),ripples(:,1));
PETH(hpc(:,1),ripples(:,1));
figure; PETH(pfc(:,1),ripples(:,1));
PETH(pfc(:,1),ripples(:,1),'durations',[-1 1]*0.2);
PETH(hpc(:,1),ripples(:,1),'durations',[-1 1]*0.2);
d = diff(ripples(:,[1 2]),[],2);
[parentFolder,basename] = fileparts(basepath);
[~,projectName] = fileparts(parentFolder);
sessionID = [projectName '_' basename];
sessionIDs{i,1} = sessionID;
sws = [0 Inf];
try
basename = basenameFromBasepath(basepath);
load(fullfile(basepath,[basename '.SleepState.states.mat']));
sws = SleepState.ints.NREMstate;
catch
['No sleep scoring for session ' num2str(i) ': ' basepath]
end
preSleep = SubtractIntervals(sleep(1,:), SubtractIntervals([0 Inf],sws));
postSleep = SubtractIntervals(sleep(2:end,:), SubtractIntervals([0 Inf],sws));
limit = Unshift(3600,preSleep); preSleep(preSleep(:,1)>limit,:) = []; if preSleep(end,2)>limit, preSleep(end,2) = limit; end
limit = Unshift(3600,postSleep); postSleep(postSleep(:,1)>limit,:) = []; if postSleep(end,2)>limit, postSleep(end,2) = limit; end
training = ConsolidateIntervals(MergePoints.timestamps(cellfun(@(x) ~isempty(strfind(x,'train')),MergePoints.foldernames),:));
exploration = sortrows([clickers{2};clickers{3}]); exploration(any(isnan(exploration),2),:) = []; exploration = ConsolidateIntervals(exploration);
exploration = Restrict(exploration,training);
pre = InIntervals(ripples,preSleep);
post = InIntervals(ripples,postSleep);
sum(pre)
sum(post)
size(ripples)
[templates,correlations,weights] = ActivityTemplatesICA(Restrict(pfc,intervals,'shift','on'),'binsize',0.1);
re = ReactivationStrength(pfc,templates,'step',0.01,'binSize',0.1);
plot(pfc(:,1));
pfc=  sortrows(pfc);\
pfc=  sortrows(pfc);
[templates,correlations,weights] = ActivityTemplatesICA(Restrict(pfc,intervals,'shift','on'),'binsize',0.1);
re = ReactivationStrength(pfc,templates,'step',0.01,'binSize',0.1);
open re
PETH(re(:,[1 2]),ripples(:,1));
clear all
batch = StartBatch(@BatchLoadHPCPFCStimData,'OJR_shortTraining_optoRipples_delayedPFC.batch');
X = get(batch,'UserData');
567
z = nanzscore(q,[],2);
meanResponse = nanmean(z(:,[find(InIntervals(qt,[-0.1 0.1])) find(InIntervals(qt,[-0.1 0.1]))+length(qt)]),2);
flip = meanResponse<0; z(flip,:) = -z(flip,:);
z1 = z(:,1:length(qt)); z2 = z(:,(length(qt)+1):end);
smooth = [0 2];
% figure('name',name);
clf
set(gcf,'position',[1 1 1920 1024])
matrices = {z1,z2,z2-z1};
titles = {[name ' in pre-task sleep'],[name ' in post-task sleep'],'(post-pre) difference'};
colors = {'k','r','m'};
for j=1:3
subplot(2,3,j); PlotColorMap(Smooth(matrices{j},smooth),'x',qt); PlotHVLines(0,'v','k--','linewidth',2);
PlotHVLines(lines,'h','w--','linewidth',2);
set(gca,'ytick',mean([[0;lines(1:end-1)] lines],2),'yticklabel',sessions);
title(titles{j});
subplot(2,3,3+j); semplot(qt,matrices{j},colors{j},smooth(end));
if j==1, semplot(qt,matrices{j+1},colors{j+1},smooth(end)); y=ylim; end
ylim(y); PlotHVLines(0,'v','k--','linewidth',2); PlotHVLines(0,'h','k--','linewidth',2); ylim(y);
ylabel('mean +/ s.e.m. (z-units)');
end
clims
drawnow
for j=1:3
subplot(2,3,j); PlotColorMap(Smooth(matrices{j},smooth),~isnan(matrices{j}),'x',qt); PlotHVLines(0,'v','k--','linewidth',2);
PlotHVLines(lines,'h','w--','linewidth',2);
set(gca,'ytick',mean([[0;lines(1:end-1)] lines],2),'yticklabel',sessions);
title(titles{j});
subplot(2,3,3+j); semplot(qt,matrices{j},colors{j},smooth(end));
if j==1, semplot(qt,matrices{j+1},colors{j+1},smooth(end)); y=ylim; end
ylim(y); PlotHVLines(0,'v','k--','linewidth',2); PlotHVLines(0,'h','k--','linewidth',2); ylim(y);
ylabel('mean +/ s.e.m. (z-units)');
end
clims
clf
set(gcf,'position',[1 1 1920 1024])
matrices = {z1,z2,z2-z1};
titles = {[name ' in pre-task sleep'],[name ' in post-task sleep'],'(post-pre) difference'};
colors = {'k','r','m'};
for j=1:3
subplot(2,3,j); PlotColorMap(Smooth(matrices{j},smooth),~isnan(matrices{j}),'x',qt); PlotHVLines(0,'v','k--','linewidth',2);
PlotHVLines(lines,'h','w--','linewidth',2);
set(gca,'ytick',mean([[0;lines(1:end-1)] lines],2),'yticklabel',sessions);
title(titles{j});
subplot(2,3,3+j); semplot(qt,matrices{j},colors{j},smooth(end));
if j==1, semplot(qt,matrices{j+1},colors{j+1},smooth(end)); y=ylim; end
ylim(y); PlotHVLines(0,'v','k--','linewidth',2); PlotHVLines(0,'h','k--','linewidth',2); ylim(y);
ylabel('mean +/ s.e.m. (z-units)');
end
clims
i=1
[pfc,hpc,deltas,ripples,~,sleep,session,MergePoints,clickers,basepath] = X{i,:};
hist(ripples(:,1))
hist(ripples(:,1),1000)
cd(basepath)
load('day7.pulses.events.mat')
q = pulses.timestamps(pulses.eventGroupID==2);
q1 = pulses.timestamps(pulses.eventGroupID==1);
hist(q1,100);
hist(q2,100);
PETH(pfc(:,1),q1)
[h,ht] = PETH(pfc(:,1),q1)
[h,ht] = PETH(pfc(:,1),q1);
PlotColorMap(h);
open q
1294
figure; hist(q,1000);
q(1294)
PlotHVLines(ans,'v','w--','linewidth',2);
clf
set(gcf,'position',[1 1 1920 1024])
matrices = {z1,z2,z2-z1};
titles = {[name ' in pre-task sleep'],[name ' in post-task sleep'],'(post-pre) difference'};
colors = {'k','r','m'};
for j=1:3
subplot(2,3,j); PlotColorMap(Smooth(matrices{j},smooth),~isnan(matrices{j}),'x',qt); PlotHVLines(0,'v','k--','linewidth',2);
PlotHVLines(lines,'h','w--','linewidth',2);
set(gca,'ytick',mean([[0;lines(1:end-1)] lines],2),'yticklabel',sessions);
title(titles{j});
subplot(2,3,3+j); semplot(qt,matrices{j},colors{j},smooth(end));
if j==1, semplot(qt,matrices{j+1},colors{j+1},smooth(end)); y=ylim; end
ylim(y); PlotHVLines(0,'v','k--','linewidth',2); PlotHVLines(0,'h','k--','linewidth',2); ylim(y);
ylabel('mean +/ s.e.m. (z-units)');
end
clims
PlotColorMap(h);
q(2976)
postSleep
open sleep
q(2976)-sleep(2,1)
q(1294)-sleep(2,1)
clear all
batch = StartBatch(@BatchLoadHPCPFCStimData,'OJR_shortTraining_optoRipples_delayedPFC.batch');
cd N:\OJRproject\OJR34\day7
load('day7.session.mat')
[pulses] = getAnalogPulses('samplingRate',session.extracellular.sr);
localPaths
filetarget
dbquit
[pulses] = getAnalogPulses('samplingRate',session.extracellular.sr);
saveas(gca,['pulses\analogPulsesDetection.png']);
filetarget = split(pwd,filesep); filetarget = filetarget{end};
if ~isempty(pul) % if no pulses, not save anything...
pulses.timestamps = stackCell(pul);
pulses.amplitude = stackCell(val);
pulses.duration = stackCell(dur);
pulses.eventGroupID = stackCell(eventGroupID);
pulses.analogChannel = stackCell(eventChannel);
intsPeriods = [];
for ii = 1:length(stimPer)
intsPeriods = [intsPeriods; stimPer{ii}];
end
pulses.intsPeriods = intsPeriods;
% sorting output
[~, idx] = sort(pulses.timestamps(:,1));
pulses.timestamps = pulses.timestamps(idx,:);
pulses.amplitude = pulses.amplitude(idx,:);
pulses.duration = pulses.duration(idx,:);
pulses.eventGroupID = pulses.eventGroupID(idx,:);
pulses.analogChannel = pulses.analogChannel(idx,:);
[~, idx] = sort(pulses.intsPeriods(:,1));
pulses.intsPeriods = pulses.intsPeriods(idx,:);
% remove pulses that are too short (noise)
if minDur
ind = find(pulses.duration > minDur);
pulses.timestamps = pulses.timestamps(ind,:);
pulses.amplitude = pulses.amplitude(ind,:);
pulses.duration = pulses.duration(ind,:);
pulses.eventGroupID = pulses.eventGroupID(ind,:);
pulses.analogChannel = pulses.analogChannel(ind,:);
% need to deal with intsPeriods
end
disp('Saving locally...');
save([filetarget '.pulses.events.mat'],'pulses');
else
pulses = [];
end
cd(prevPath);
load('day7.pulses.events.mat')
open pulses
pulses.timestamps(34)
addComma(pulses.timestamps(34))
addComma(pulses.timestamps(50))
pulses.eventGroupID(50))
pulses.eventGroupID(50)
pulses.eventGroupID(34)
pulses.eventGroupID(35)
pulses.eventGroupID(38)
lfp = GetAyaLFP(50);
PETH(lfp,pulses.timestamps);
PETH(lfp,pulses.timestamps(:,1));
PETH(lfp,pulses.timestamps(:,1),'durations',[-1 1]*5);
[pulses] = getAnalogPulses('samplingRate',session.extracellular.sr);
[pulses] = getAnalogPulses('samplingRate',session.extracellular.sr,'manualThr',true);
dbcont
PETH(lfp,pulses.timestamps(pulses.eventGroupID==1,1),'durations',[-1 1]*5);
PETH(lfp,pulses.timestamps(pulses.eventGroupID==2,1),'durations',[-1 1]*5);
batch = StartBatch(@BatchLoadHPCPFCStimData,'OJR_shortTraining_optoRipples_delayedPFC.batch');
r = [ripples.timestamps(:,1) ripples.peaks ripples.timestamps(:,2) ripples.amplitude]; r(:,5) = nan;
dbcont
r = [ripples.timestamps(:,1) ripples.peaks ripples.timestamps(:,2) ripples.amplitude]; r(:,5) = nan;
dbcont
r = [ripples.timestamps(:,1) ripples.peaks ripples.timestamps(:,2) ripples.amplitude]; r(:,5) = nan;
[ripples.timestamps(:,1) ripples.peaks ripples.timestamps(:,2) ripples.amplitude];
ripples.timestamps = zeros(0,2); ripples.peaks = []; ripples.amplitude = [];
r = [ripples.timestamps(:,1) ripples.peaks ripples.timestamps(:,2) ripples.amplitude]; r(:,5) = nan;
dbcont
cd N:\OJRproject\OJR43\day10
cd N:\OJRproject\OJR43\day7
clear all
45
figure; PlotXY(lfp);
PlotIntervals(badIntervals,'color','r');
load('day7.ripples.events.mat')
[clean,~,badIntervals] = CleanLFP(lfp,'thresholds',[6 7]);
PlotIntervals(badIntervals,'color','r');
[clean,~,badIntervals] = CleanLFP(lfp,'thresholds',[5 7]);
PlotIntervals(badIntervals,'color','r');
[clean,~,badIntervals] = CleanLFP(lfp,'thresholds',[5 7]);
figure; plot(t,abs(d));
% maybe 1.5
PlotHVLines(1.5,'h','r');
PlotHVLines(5269,'v','r');
PlotHVLines(5271.3,'v','k');
figure; plot(t,abs(z));
PlotHVLines(3,'h','k');
PlotHVLines(1.5,'h','k');
PlotHVLines(1,'h','k');
PlotHVLines(0.5,'h','k');
load('day7.EMGFromLFP.LFP.mat')
figure; plot(EMGFromLFP.timestamps,EMGFromLFP.data);
load('day7.SleepState.states.mat')
[clean,~,badIntervals] = CleanLFP(lfp,'thresholds',[3 1.5]);
dur(badIntervals)/3600
deltas0 = FindDeltaPeaks(clean);
%     deltas = deltas0(deltas0(:,5)>1.5 & deltas0(:,5)<4.5,:); % these thresholds should be manually refined for each session
deltas = deltas0(deltas0(:,5)-deltas0(:,6)>3,:); % these thresholds should be manually refined for each session
deltaWaves.timestamps = deltas(:,[1 3]); deltaWaves.peaks = deltas(:,2);  deltaWaves.peakNormedPower = deltas(:,5); deltaWaves.detectorName = ['channel ' num2str(channel) '(+1), CleanLFP, FindDeltaPeaks, peak-trough>3'];
% save(fullfile(basepath,[basename '.deltaWaves.events.mat']),'deltaWaves');
r = [ripples.timestamps(:,1) ripples.peaks ripples.timestamps(:,2) ripples.amplitude]; r(:,5) = nan;
PETH(r(:,[2]),deltas(:,2));
clf
PETH(r(:,[2]),deltas(:,2));
[parentFolder,basename] = fileparts(basepath);
[~,projectName] = fileparts(parentFolder);
load(fullfile(basepath,[basename '.cell_metrics.cellinfo.mat']),'cell_metrics');
nNeurons = length(cell_metrics.spikes.times);
regions = cell_metrics.brainRegion;
regionNames = unique(regions);
regionCell = zeros(length(regions),1);
for i=1:length(regionNames)
regionCell(strcmp(regions,regionNames{i})) = i;
end
spikesCell = cell_metrics.spikes.times';
PFCindex = [find(strcmp(regionNames,'ILA')) find(strcmp(regionNames,'PFC'))];
HPCindex = find(strcmp(regionNames,'CA1'));
pfc = []; if ~isempty(PFCindex), pfc = sortrows(Group(spikesCell{regionCell==PFCindex})); end
pwd
basepath = pwd
basename = basenameFromBasepath(basepath);
load(fullfile(basepath,[basename '.cell_metrics.cellinfo.mat']),'cell_metrics');
nNeurons = length(cell_metrics.spikes.times);
regions = cell_metrics.brainRegion;
regionNames = unique(regions);
regionCell = zeros(length(regions),1);
for i=1:length(regionNames)
regionCell(strcmp(regions,regionNames{i})) = i;
end
spikesCell = cell_metrics.spikes.times';
PFCindex = [find(strcmp(regionNames,'ILA')) find(strcmp(regionNames,'PFC'))];
HPCindex = find(strcmp(regionNames,'CA1'));
pfc = []; if ~isempty(PFCindex), pfc = sortrows(Group(spikesCell{regionCell==PFCindex})); end
[h,ht] = PETH(pfc(:,1),deltas0(:,2));
PlotColorMap(Smooth(Shrink(sortby(h,deltas(:,5)-deltas(:,6),100,1),[0 5]),'x',qt)
PlotColorMap(Smooth(Shrink(sortby(h,deltas(:,5)-deltas(:,6),100,1),[0 5]),'x',ht))
PlotColorMap(Smooth(Shrink(sortby(h,deltas(:,5)-deltas(:,6)),100,1),[0 5]),'x',ht)
PlotColorMap(Smooth(Shrink(sortby(h,deltas0(:,5)-deltas0(:,6)),100,1),[0 5]),'x',ht)
size(h)
size(deltas0)
PETH(pfc(:,1),deltas0(:,2));
s = deltas0(:,5)-deltas0(:,6);
clf
hist(s,100);
PETH(pfc(:,1),deltas0(s<3,2));
PETH(pfc(:,1),deltas0(s<5,2));
hist(s,100);
PETH(pfc(:,1),deltas0(s<5,2));
PETH(pfc(:,1),deltas0(s<4,2));
PETH(pfc(:,1),deltas0(s<3.5,2));
PlotColorMap(Smooth(Shrink(sortby(h,s),1,1),[0 5]),'x',ht)
PlotColorMap(Smooth(Shrink(sortby(h,s),10,1),[0 5]),'x',ht)
PlotColorMap(Smooth(Shrink(sortby(h,s),100,1),[0 5]),'x',ht)
figure; hist(deltas(:,2));
DensityMap(deltas0(:,2),s));
DensityMap(deltas0(:,2),s);
DensityMap(deltas0(:,2),s,'nBins',100);
DensityMap(deltas0(:,2),s,'nBins',100,'smooth',0);
[qq,qt] = DensityMap(deltas0(:,2),s,'nBins',100,'smooth',0);
open qt
[qq,qx,qy] = DensityMap(deltas0(:,2),s,'nBins',100,'smooth',0);
PlotColorMap(qq,'x',qx,'y',qy);
PlotColorMap(zscore(qq,[],2),'x',qx,'y',qy);
PlotColorMap(Smooth(zscore(qq,[],2),[0 5]),'x',qx,'y',qy);
[qq,qx,qy] = DensityMap(deltas0(:,2),tiedrank(s),'nBins',100,'smooth',0);
[qq,qx,qy] = DensityMap(deltas0(:,2),tiedrank(s),'nBins',100,'smooth',[0 5]);
[qq,qx,qy] = DensityMap(deltas0(:,2),tiedrank(s),'nBins',100,'smooth',[5 0]);
quantile(s,4000/14000)
b = Bin(s,100);
for i=1:100, hh(i,:) = nanmean(h(b==i,:)); end
PlotColorMap(hh);
b = Bin(tiedrank(s),100);
for i=1:100, hh(i,:) = nanmean(h(b==i,:)); end
PlotColorMap(hh);
semplot(hh(1:40,:))
clf
PlotColorMap(hh);
PlotColorMap(zscore(hh,[],2));
PlotColorMap(Smooth(zscore(hh,[],2),5))
PlotColorMap(Smooth(zscore(hh,[],2),[2 2))
PlotColorMap(Smooth(zscore(hh,[],2),[2 2]))
load('day7.SleepState.states.mat')
bad = InIntervals(deltas0(:,2),SleepState.ints.WAKEstate);
Portion(bad)
semplot(ht,h(bad,:))
clf
Dist(100,[s 2-bad],'grouped');
clf
PlotColorMap(Smooth(Shrink(sortby(h,s),100,1),[0 5]),'x',ht)
PlotColorMap(Smooth(Shrink(sortby(h(~bad,:),s(~bad)),100,1),[0 5]),'x',ht)
[h,ht] = PETH(pfc(:,1),deltas0(:,2),'durations',[-1 1]*5);
PlotColorMap(Smooth(Shrink(sortby(h(~bad,:),s(~bad)),100,1),[0 5]),'x',ht)
z = zscore(h,[],2);
PlotColorMap(Smooth(Shrink(sortby(z(~bad,:),s(~bad)),100,1),[0 5]),'x',ht)
semplot(ht,z);
semplot(ht,z(bad,:));
clf
semplot(ht,z(bad,:));
semplot(ht,z(~bad,:));
semplot(ht,z(~bad,:),'r');
clf
anovabar(s,bad)
anovabar(s>3,bad)
anovabar(s>4,bad)
sum(s>4)
sum(s>4 & ~bad)
deltas = deltas0(deltas0(:,5)-deltas0(:,6)>4,:); % these thresholds should be manually refined for each session
clf
[h,ht] = PETH(pfc(:,1),deltas(:,2),'durations',[-1 1]*5);
ss = deltas(:,5)-deltas(:,6);
deltas = deltas0(deltas0(:,5)-deltas0(:,6)>4 & ~bad,:); % these thresholds should be manually refined for each session
ss = deltas(:,5)-deltas(:,6);
deltas = deltas0(deltas0(:,5)-deltas0(:,6)>4 & ~bad,:); % these thresholds should be manually refined for each session
deltas = deltas0(deltas0(:,5)-deltas0(:,6)>4,:); % these thresholds should be manually refined for each session
ss = deltas(:,5)-deltas(:,6);
PlotColorMap(Smooth(Shrink(sortby(h,ss),100,1),[0 5]),'x',ht)
PlotColorMap(Smooth(Shrink(sortby(h,ss),100,1),[0 1]),'x',ht)
deltas = deltas0(deltas0(:,5)-deltas0(:,6)>4 & ~bad,:); % these thresholds should be manually refined for each session
ss = deltas(:,5)-deltas(:,6);
[h,ht] = PETH(pfc(:,1),deltas(:,2),'durations',[-1 1]*5);
PlotColorMap(Smooth(Shrink(sortby(h,ss),100,1),[0 1]),'x',ht)
PlotColorMap(Smooth(Shrink(sortby(h,deltas(:,5)),100,1),[0 1]),'x',ht)
PlotColorMap(Smooth(Shrink(sortby(h,deltas(:,6)),100,1),[0 1]),'x',ht)
PlotColorMap(Smooth(Shrink(sortby(h,deltas(:,4)),100,1),[0 1]),'x',ht)
PlotColorMap(Smooth(Shrink(sortby(h,deltas(:,5)),100,1),[0 1]),'x',ht)
PlotColorMap(Smooth(Shrink(sortby(h,deltas(:,5)-deltas(:,4)),100,1),[0 1]),'x',ht)
PlotColorMap(Smooth(Shrink(sortby(h,deltas(:,5)-deltas(:,6)),100,1),[0 1]),'x',ht)
clim
clim([4 6])
clim([4.5 6])
clim([5 6])
clim([5 7])
hist(deltas(:,5)-deltas(:,6))
hist(deltas(:,5)-deltas(:,6),100)
figure; DensityMap(deltas(:,2),deltas(:,5)-deltas(:,6),'nBins',100,'smooth',[5 0]);
list = findmin(deltas(:,5)-deltas(:,6),10)
i=1; addComma(deltas(list(i),2))
i=i+1; addComma(deltas(list(i),2))
list(i)
addComma(deltas(2738,2))
addComma(deltas(2738+1,2))
addComma(deltas(2738+2,2))
size(ss)
ss(2738:(2738+2))
ss = deltas(:,5)-mean([deltas(:,6) deltas(:,4)],2);
figure; hist(ss,100)
ss(2738:(2738+2))
quantile(s,127/144)
list = findmin(ss,10)
i=0;
i=i+1; addComma(deltas(list(i),2))
addComma(deltas(list(i)+1,2))
addComma(deltas(list(i)-1,2))
FindClosest(deltas0(:,2),19049.394)
deltas0(11985,:)
deltas0(11985,4:end)
2.4100+3.5631
clf
hist(s,100)
quantile(s,43,92)
quantile(s,43/92)
list = findmin(ss,10)
i=0
i=i+1; addComma(deltas(list(i),2))
addComma(deltas(list(i)+1,2))
FindClosest(deltas0(:,2),3500.375)
deltas0(2262:end)
deltas0(2262,4:end)
addComma(deltas(list(i)+1,2))
deltas0(FindClosest(deltas0(:,2),3511.239),4:6)
deltas0(FindClosest(deltas0(:,2),3514.337),4:6)
deltas0(FindClosest(deltas0(:,2),3514.337)+1,4:6)
deltas0(FindClosest(deltas0(:,2),3514.337)+1,2)-3514
deltas0(FindClosest(deltas0(:,2),3524.517)+1,4:6)
deltas0(FindClosest(deltas0(:,2),3524.517),4:6)
deltas0(FindClosest(deltas0(:,2),3524.517),2)-
deltas0(FindClosest(deltas0(:,2),3524.517),2)-3524
figure; PlotXY(clean); xlim([0 3]+3523);
deltas00 = FindDeltaPeaks(clean);
figure; PlotXY(z);
xlim([0 3]+3523);
this = Restrict(lfp,[0 3]+3523);
[w,wt,wf] = WaveletSpectrogram(this,'range',[0.05 50]);
figure; PlotColorMap(w,'x',wt,'y','wf');
PlotColorMap(w,'x',wt,'y',wf);
figure; PlotXY(Restrict(filtered,[0 3]+3523));
subplot(3,1,1);
PlotXY(Restrict(lfp,[0 3]+3523));
subplot(3,1,2);
PlotXY(Restrict(filtered,[0 3]+3523));
subplot(3,1,3);
PlotXY(Restrict(z,[0 3]+3523));
% List positions (in # samples) of successive up,down,up crossings in an Nx3 matrix
n = length(up);
where = [up(1:n-1) down(1:n-1) up(2:n)];
PltoHVLines(Restrict(up,xlim),'v','r');
PlotHVLines(Restrict(up,xlim),'v','r');
PlotHVLines(Restrict(t(up),xlim),'v','r');
t = lfp(:,1);
PlotHVLines(Restrict(t(up),xlim),'v','r');
PlotHVLines(Restrict(t(down),xlim),'v','b');
% List positions but also z-scored amplitudes in an Nx6 matrix (positions then amplitudes)
z = filtered;
z(:,2) = zscore(filtered(:,2));
subplot(3,1,3); cla
PlotXY(Restrict(z,[0 3]+3523));
deltas = z(where,:);
deltas = reshape(deltas,size(where,1),6);
open deltas
Restrict(deltas,xlim)
PlotIntervals(Restrict(deltas(:,1:3),xlim));
PlotIntervals(Restrict(deltas(:,[1 3]),xlim));
PlotHVLines(Restrict(deltas(:,[1 3]),xlim),'v','k--');
duration = deltas(:,3) - deltas(:,1);
in = find(InIntervals(deltas(:,2),xlim)
in = find(InIntervals(deltas(:,2),xlim))
in = find(InIntervals(deltas(:,2),[3524.4 3524.6])))
in = find(InIntervals(deltas(:,2),[3524.4 3524.6]))
duration(13596)
minDuration
maxDuration
% Discard waves that are too long or too short
duration = deltas(:,3) - deltas(:,1);
deltas(duration<minDuration/1000|duration>maxDuration/1000,:) = [];
in = find(InIntervals(deltas(:,2),[3524.4 3524.6]))
% Threshold z-scored peak and trough amplitudes
base = deltas(:,4);
peak = deltas(:,5);
trough = deltas(:,6);
case1 = peak > highPeak & trough <= -lowTrough;
case2 = peak >= lowPeak & trough < -highTrough;
peakabovebase = peak > 0.15 + base;
peakabovebase(13067)
case1(case1)
case1(13067)
case2(13067)
deltas(13067,4:end)
highPeak
lowTrough
lowPeak
highTrough
in = find(InIntervals(deltas(:,2),[3524.4 3524.6]))
in = find(InIntervals(deltas(:,2),[3525.4 3525.6]))
in = find(InIntervals(deltas(:,2),[3525.4 3526]))
addComma(deltas(13071,2))
deltas(13071,4:6)
lowPeak
highTrough
case1(13071)
case2(13071)
case2(13072)
deltas(13072,4:6)
subplot(3,1,1);
PlotXY(Restrict(lfp,[0 3]+3523.5));
subplot(3,1,2);
PlotXY(Restrict(filtered,[0 3]+3523.5));
subplot(3,1,3);
PlotXY(Restrict(z,[0 3]+3523.5));
addComma(deltas(13072,2))
PlotHVLines(Restrict(deltas(:,[1 3]),xlim),'v','k--');
PlotHVLines(Restrict(deltas(:,2),xlim),'v','y','linewidth',2);
slow = FilterLFP(lfp,'passband',[0 2],'order',8);
PlotXY(Restrict(slow,[0 3]+3523.5),'r');
slow = FilterLFP(lfp,'passband',[0 2]);
PlotXY(Restrict(slow,[0 3]+3523.5),'r');
cla
PlotXY(Restrict([lfp(:,1) lfp(:,2)-slow(:,2)],[0 3]+3523.5),'r');
filtered = FilterLFP(lfp,'passband',[0 9],'order',8);
cla
PlotXY(Restrict([filtered],[0 3]+3523.5),'r');
interval = p0 3]+3523.5;
interval = [0 3]+3523.5;
interval = [0 3]+3532;
subplot(3,1,1);
PlotXY(Restrict(lfp,interval));
subplot(3,1,2);
PlotXY(Restrict(filtered,interval));
subplot(3,1,3);
PlotXY(Restrict(z,interval));
xlim(interval);
[parentFolder,basename] = fileparts(basepath);
[~,projectName] = fileparts(parentFolder);
basepath = pwd
basename = basenameFromBasepath(basepath);
load(fullfile(basepath,[basename '.cell_metrics.cellinfo.mat']),'cell_metrics');
nNeurons = length(cell_metrics.spikes.times);
regions = cell_metrics.brainRegion;
regionNames = unique(regions);
regionCell = zeros(length(regions),1);
for i=1:length(regionNames)
regionCell(strcmp(regions,regionNames{i})) = i;
end
spikesCell = cell_metrics.spikes.times';
PFCindex = [find(strcmp(regionNames,'ILA')) find(strcmp(regionNames,'PFC'))];
HPCindex = find(strcmp(regionNames,'CA1'));
pfc = []; if ~isempty(PFCindex), pfc = sortrows(Group(spikesCell{regionCell==PFCindex})); end
RasterPlot(Restrict(pfc,interval));
RasterPlot(Restrict([pfc(:,1) pfc(:,2)+1],interval));
RasterPlot(Restrict([pfc(:,1) pfc(:,2)+3],interval));
RasterPlot(Restrict([pfc(:,1) pfc(:,2)+2],interval));
fgamma = FilterLFP(lfp,'passband',[20 100]);
PlotXY(Restrict(fgamma,interval));
[~,gamma]= Phase(fgamma);
open gamma
PlotXY(Restrict(gamma,interval));
cla
PlotXY(Restrict(gamma,interval));
xlim(interval);
interval = [0 3]+3523.5;
PlotXY(Restrict(gamma,interval));
xlim(interval);
PlotXY(Smooth(Restrict(gamma,interval),[5 0]));
cla
PlotXY(Smooth(Restrict(gamma,interval),[5 0]));
PlotXY(Smooth(Restrict(gamma,interval),[0.1*1250 0]));
cla
PlotXY(Smooth(Restrict(gamma,interval),[0.05*1250 0]));
PlotXY(Smooth(Restrict([gamma(:,1) -gamma(:,2)],interval),[0.01*1250 0]));
cla
PlotXY(Smooth(Restrict([gamma(:,1) -gamma(:,2)],interval),[0.01*1250 0]));
hold off
PlotXY(Smooth(Restrict([gamma(:,1) -gamma(:,2)],interval),[0.02*1250 0]));
hist(Restrict(pfc(:,1),interval),linspace(interval(1),interval(2),100));
xlim(interval);
hist(Restrict(pfc(:,1),interval),linspace(interval(1),interval(2),500));
xlim(interval);
PlotHVLines(-100,'h','k');
PlotHVLines(-200,'h','k');
subplot(4,1,4);
[h,ht] = hist(Restrict(pfc(:,1),interval),linspace(interval(1),interval(2),100));
plot(ht,Smooth(h,1));
[h,ht] = hist(Restrict(pfc(:,1),interval),linspace(interval(1),interval(2),1000));
plot(ht,Smooth(h,1));
subplot(4,1,4);
[h,ht] = hist(Restrict(pfc(:,1),interval),linspace(interval(1),interval(2),1000));
bar(ht,Smooth(h,1));
FindClosest(deltas(:,2),3527.663)
deltas(13078,4:6)
1.9975+3.2223
FindClosest(deltas(:,2),3540.89)
deltas(13128,4:6)
2.9690   +2.8751
deltas(FindClosest(deltas(:,2),3542.192),4:6)
deltas(FindClosest(deltas(:,2),3562.853),4:6)
deltas(FindClosest(deltas(:,2),3568.899),4:6)
2.5280   -2.1090
2.5280   +2.1090
deltas(FindClosest(deltas(:,2),3574.706),4:6)
subplot(4,1,2);
PlotXY(Restrict(z,interval));
ylim(ylim);
z0 = z;
% List positions but also z-scored amplitudes in an Nx6 matrix (positions then amplitudes)
z = filtered;
z(:,2) = zscore(filtered(:,2));
subplot(4,1,2); hold all
PlotXY(Restrict(z,interval));
deltas0 = FindDeltaPeaks(clean);
n = length(up);
where = [up(1:n-1) down(1:n-1) up(2:n)];
% List positions but also z-scored amplitudes in an Nx6 matrix (positions then amplitudes)
z = filtered;
z(:,2) = zscore(filtered(:,2));
deltas = z(where,:);
deltas = reshape(deltas,size(where,1),6);
deltas(FindClosest(deltas(:,2),3574.706),4:6)
2.4747   +1.4180
deltas0 = FindDeltaPeaks(clean);
45
clim
clim([0.6 1.6])
clim([0.6 1.7])
clim([0.7 1.7])
clim([0.8 1.7])
clim([1 1.7])
quantile(ss,18/73)
deltaWaves.timestamps = deltas(:,[1 3]); deltaWaves.peaks = deltas(:,2);  deltaWaves.peakNormedPower = deltas(:,5); deltaWaves.detectorName = ['channel 114(+1), CleanLFP, FindDeltaPeaks, peak-trough>3.5 & ~inAwake'];
save(fullfile(basepath,[basename '.deltaWaves.events.mat']),'deltaWaves');
r = [ripples.timestamps(:,1) ripples.peaks ripples.timestamps(:,2) ripples.RipMax ripples.SwMax];
load(fullfile(basepath,[basename '.MergePoints.events.mat']),'MergePoints');
sleep = ConsolidateIntervals(MergePoints.timestamps(cellfun(@(x) any(strfind(lower(x),'sleep')),MergePoints.foldernames),:))
clf
preSleep = SubtractIntervals(sleep(1,:), SubtractIntervals([0 Inf],sws));
sws = SleepState.ints.NREMstate;
preSleep = SubtractIntervals(sleep(1,:), SubtractIntervals([0 Inf],sws));
postSleep
postSleep = SubtractIntervals(sleep(2:end,:), SubtractIntervals([0 Inf],sws));
pre = InIntervals(ripples,preSleep);
pre = InIntervals(r,preSleep);
post = InIntervals(r,postSleep);
semplot(deltas(:,2),r(:,2));
PETH(deltas(:,2),r(:,2));
[h,ht] = PETH(deltas(:,2),r(:,2));
semplot(ht,h(pre,:));
semplot(ht,h(post,:),'r'
);
semplot(ht,h(post,:),'r')
clf
dpre = InIntervals(deltas(:,2),preSleep);
dpost = InIntervals(deltas(:,2),postSleep);
[h,ht] = PETH(r(:,2),deltas(:,2));
semplot(ht,h(dpost,:),'r')
semplot(ht,h(dpre,:),'k')
fol = InIntervals(r(:,2),[deltas(:,2) deltas(:,2)+0.2]);
prec = InIntervals(r(:,2),[deltas(:,2)-0.2 deltas(:,2)]);
anovabar(prec,[-pre+post]);
clf
anovabar(prec,[-pre+post]);
ok = pre|post;
anovabar(prec(ok),[-pre(ok)+post(ok)]);
anovabar(fol(ok),[-pre(ok)+post(ok)]);
anovabar(prec(ok),[-pre(ok)+post(ok)]);
dprec = InIntervals(deltas(:,2),[r(:,2)-0.2 r(:,2)]);
dfol = InIntervals(deltas(:,2),[r(:,2) r(:,2)+0.2]);
anovabar(dprec(ok),[-dpre(ok)+dpost(ok)]);
ok = dpre|dpost;
anovabar(dprec(ok),[-dpre(ok)+dpost(ok)]);
anovabar(dfol(ok),[-dpre(ok)+dpost(ok)]);
anovabar([dprec(ok) dfol(ok)],[-dpre(ok)+dpost(ok)]);
anovabar([dprec(ok) dfol(ok) dfol(ok)],[-dpre(ok)+dpost(ok)]);
close all
\clear all
clear all
cd N:\OJRproject\OJR43\day10
cd N:\OJRproject\OJR43\day3
basepath = pwd
basename = basenameFromBasepath(basepath);
45
basepath = pwd
basename = basenameFromBasepath(basepath);
load(fullfile(basepath,[basename '.cell_metrics.cellinfo.mat']),'cell_metrics');
nNeurons = length(cell_metrics.spikes.times);
regions = cell_metrics.brainRegion;
regionNames = unique(regions);
regionCell = zeros(length(regions),1);
for i=1:length(regionNames)
regionCell(strcmp(regions,regionNames{i})) = i;
end
spikesCell = cell_metrics.spikes.times';
PFCindex = [find(strcmp(regionNames,'ILA')) find(strcmp(regionNames,'PFC'))];
HPCindex = find(strcmp(regionNames,'CA1'));
pfc = []; if ~isempty(PFCindex), pfc = sortrows(Group(spikesCell{regionCell==PFCindex})); end
load('day3.SleepState.states.mat')
figure; PlotXY(lfp);
PlotIntervals(badIntervals,'color','r');
figure; PlotXY(clean);
PlotIntervals(badIntervals,'color','r');
[clean,~,badIntervals] = CleanLFP(lfp,'thresholds',[5 7]);
PlotIntervals(badIntervals,'color','r');
[clean,~,badIntervals] = CleanLFP(lfp,'thresholds',[4 7]);
PlotIntervals(badIntervals,'color','r');
[clean,~,badIntervals] = CleanLFP(lfp,'thresholds',[3 7]);
PlotIntervals(badIntervals,'color','r');
[clean,~,badIntervals] = CleanLFP(lfp,'thresholds',[2.5 7]);
PlotIntervals(badIntervals,'color','r');
[clean,~,badIntervals] = CleanLFP(lfp,'thresholds',[2 7]);
PlotIntervals(badIntervals,'color','r');
clf
PlotXY(clean);
[clean,~,badIntervals] = CleanLFP(lfp,'thresholds',[2 7]);
figure; plot(t,z);
% Substitute noisy signal with interpolated signal as if artefact did not exist
values(bad) = interp1(t(~bad),values(~bad),t(bad,1));
badIntervals = ConsolidateIntervals(sortrows([artefactInterval; noisyInterval]));
ok = ~isnan(values);
clean = [t(ok) values(ok)];
clf
PlotXY(clean);
plot(t,z);
xlim([-1 1]*10+1483.38);
[clean,~,badIntervals] = CleanLFP(lfp,'thresholds',[1 7]);
dbcont
PlotXY(clean);
deltas0 = FindDeltaPeaks(clean);
%     deltas = deltas0(deltas0(:,5)>1.5 & deltas0(:,5)<4.5,:); % these thresholds should be manually refined for each session
deltas = deltas0(deltas0(:,5)-deltas0(:,6)>3,:); % these thresholds should be manually refined for each session
deltaWaves.timestamps = deltas(:,[1 3]); deltaWaves.peaks = deltas(:,2);  deltaWaves.peakNormedPower = deltas(:,5); deltaWaves.detectorName = ['channel ' num2str(channel) '(+1), CleanLFP, FindDeltaPeaks, peak-trough>3'];
% save(fullfile(basepath,[basename '.deltaWaves.events.mat']),'deltaWaves');
dur(badIntervals)
figure; PlotXY(lfp);
plot(t,z);
[clean,~,badIntervals] = CleanLFP(lfp,'thresholds',[3 7]);
figure; PlotXY(lfp);
PlotIntervals(badIntervals,'color','r');
deltas0 = FindDeltaPeaks(clean);
%     deltas = deltas0(deltas0(:,5)>1.5 & deltas0(:,5)<4.5,:); % these thresholds should be manually refined for each session
deltas = deltas0(deltas0(:,5)-deltas0(:,6)>3,:); % these thresholds should be manually refined for each session
deltaWaves.timestamps = deltas(:,[1 3]); deltaWaves.peaks = deltas(:,2);  deltaWaves.peakNormedPower = deltas(:,5); deltaWaves.detectorName = ['channel ' num2str(channel) '(+1), CleanLFP, FindDeltaPeaks, peak-trough>3'];
% save(fullfile(basepath,[basename '.deltaWaves.events.mat']),'deltaWaves');
PlotXY(clean)
PETH(pfc(:,1),deltas(:,2))
FindClosest(deltas(:,2),9483.3)
addComma(deltas(1644,2))
FindClosest(deltas(:,2),9482.2)
addComma(deltas(1642,2))
deltas(1642,4:6)
size(h)
size(deltas)
plot(h(1642,:))
plot(ht,h(1642,:))
basepath = pwd
basename = basenameFromBasepath(basepath);
load(fullfile(basepath,[basename '.cell_metrics.cellinfo.mat']),'cell_metrics');
nNeurons = length(cell_metrics.spikes.times);
regions = cell_metrics.brainRegion;
regionNames = unique(regions);
regionCell = zeros(length(regions),1);
for i=1:length(regionNames)
regionCell(strcmp(regions,regionNames{i})) = i;
end
spikesCell = cell_metrics.spikes.times';
PFCindex = [find(strcmp(regionNames,'ILA')) find(strcmp(regionNames,'PFC'))];
HPCindex = find(strcmp(regionNames,'CA1'));
pfc = []; if ~isempty(PFCindex), pfc = sortrows(Group(spikesCell{regionCell==PFCindex})); end
load('day3.ripples.events.mat')
r = [ripples.timestamps(:,1) ripples.peaks ripples.timestamps(:,2) ripples.RipMax ripples.SwMax];
PETH(r(:,1),deltas(:,2))
PETH(pfc(:,1),deltas(:,2))
SaveCustomEvents('temp.tem.evt',deltas(:,2),'putative delta peak');
FindClosest(deltas(:,2),9483.94)
addComma(deltas(1646,2))
deltas(1646,4:6)
1.8456 + 2.92
plot(ht,h(1644,:))
plot(ht,h(1646,:))
PlotColorMap(h);
PlotColorMap(h(~bad,:));
size(h)
size(bad)
PlotColorMap(Shrink(h(~bad,:),100,1));
PlotColorMap(Shrink(h,100,1));
PlotColorMap(Shrink(h,500,1));
PlotColorMap(Shrink(h,500,1),'x',ht);
FindClosest(deltas(:,2),[9491.585; 9494.394])
addComma(deltas(1653,2))
addComma(deltas(1654,2))
plot(ht,h(1653,:))
9491.585-2.6
9491.585-2.98
addComma(9491.585-2.98)
addComma(9491.585-2.6)
PlotColorMap(h(1650 + (0:5),:));
semplot(ht,h(1650 + (0:5),:))
clf
addComma(deltas(1650)
addComma(deltas(1650,2))
addComma(deltas(1649,2))
semplot(ht,h(1649 + [0 1 4],:))
clf
PlotColorMap(h(1649 + [0 1 4],:));
PlotColorMap(h(1649 + 0:5,:));
PlotColorMap(h(1649 + (0:5),:));
lfps = GetAyaLFP(96);
lfpd = GetAyaLFP(109);
lfp105 = lfp;
lfp(:,2) = lfpd(:,2)-lfps(:,2);
PETH(lfp,deltas(:,2));
[h1,ht] = PETH(lfp,deltas(:,2));
[h,ht] = PETH(pfc(:,1),deltas(:,2));
[histogram2d,h0,difference] = JointPETH(h1,h,1);
PlotColorMap(histogram2d);
size(h1)
size(h)
semplot(ht,h)
semplot(ht,h1)
mmax(lfp(:,2))
[histogram2d,h0,difference] = JointPETH(h1+7570,h,1);
PlotColorMap(histogram2d);
deltas00 = FindDeltaPeaks(clean);
PlotXY(lfp);
clean0 = clean;
clean = Restrict(lfp,SubtractIntervals([0 Inf],badIntervals));
PlotXY(clean)
deltas00 = deltas0;
deltas0 = FindDeltaPeaks(clean);
PETH(pfc(:,1),deltas0(:,2))
PETH(pfc(:,1),deltas(:,2))
PETH(pfc(:,1),deltas0(:,2))
PETH(pfc(:,1),deltas0(:,2),'nBins',501)
PETH(pfc(:,1),deltas0(:,1),'nBins',501)
PETH(pfc(:,1),deltas0(:,2),'nBins',501)
PETH(pfc(:,1),deltas0(:,2)-0.02,'nBins',501)
PETH(pfc(:,1),deltas0(:,2)-0.03,'nBins',501)
PETH(pfc(:,1),deltas0(:,2)-0.05,'nBins',501)
PETH(pfc(:,1),deltas0(:,2)-0.06,'nBins',501)
PETH(r(:,1),deltas0(:,2)-0.06,'nBins',501)
hist(s,100);
PlotHVLines(sum(s>3)/100),'h','w');
PlotHVLines(sum(s>3)/100,'h','w');
sum(s>3)/100
SaveCustomEvents('temp.tem.evt',deltas(:,2),'putative delta peak');
SaveCustomEvents('temp.tem.evt',deltas(:,1:3),{'putative delta start','putative delta peak','putative delta stop'});
filtered = FilterLFP(lfpd,'passband',[0 9],'order',8);
z = filtered;
z(:,2) = [diff(filtered(:,2));0];
z = FilterLFP(z,'passband',[0 6],'order',8);
z(:,2) = zscore(z(:,2));
% Find positions (in # samples) of zero crossings
[up,down] = ZeroCrossings(z);
down = find(down);
up = find(up);
if down(1) < up(1), down(1) = []; end
% List positions (in # samples) of successive up,down,up crossings in an Nx3 matrix
n = length(up);
where = [up(1:n-1) down(1:n-1) up(2:n)];
open where
ok = InIntervals(deltas(:,2),t(where(:,[1 3])));
t = lfp(:,1);
ok = InIntervals(deltas(:,2),t(where(:,[1 3])));
Portion(ok)
FindClosest(deltas(:,2)9489.58)
FindClosest(deltas(:,2)*9489.58)
FindClosest(deltas(:,2),9489.58)
addComma(deltas(2631,2))
interval = [-1 1]+9489 ;
PlotXY(Restrict(filtered,interval));
PlotHVLines(deltas(2631,1:3));
PlotXY(Restrict(filtered,interval+[-1 1]*30));
PlotHVLines(Restrict(deltas(:,1),xlim),'y');
PlotHVLines(Restrict(deltas(:,3),xlim),'b');
PlotHVLines(Restrict(deltas(:,2),xlim),'r--');
xlim(9489.25 + [0 5]);
PlotHVLines(Restrict(deltas(:,1),xlim),'y');
PlotHVLines(Restrict(deltas(:,3),xlim),'b');
PlotHVLines(Restrict(deltas(:,2),xlim),'r--');
find(InIntervals(deltas(:,2),xlim))
figure; PlotColorMap(h(find(InIntervals(deltas(:,2),xlim)),:),'x',ht);
size(h)
size(deltas)
xlim
find(InIntervals(deltas(:,2),9489.25 + [0 5]))
figure; PlotColorMap(h(find(InIntervals(deltas(:,2),9489.25 + [0 5])),:),'x',ht);
PlotColorMap(h(1649 + (0:5),:));
PlotColorMap(h(1649 + (0:5),:),'x',ht);
figure; PlotColorMap(h(find(InIntervals(deltas(:,2),9489.25 + [0 5])),:),'x',ht);
PlotHVLines(Restrict(deltas(:,1),xlim),'y');
PlotHVLines(Restrict(deltas(:,3),xlim),'b--');
PlotHVLines(Restrict(deltas(:,1),xlim),'g');
PlotHVLines(Restrict(deltas(:,3),xlim),'b--');
figure; PlotXY(Restrict(filtered,interval+[-1 1]*30));
xlim(9489.25 + [0 5]);
figure; PlotXY(FilterLFP(Restrict(lfp,interval+[-1 1]*30),'passband',[0 9],'order',8));
xlim(9489.25 + [0 5]);
hold all
PlotXY(Restrict(filtered,interval+[-1 1]*30));
figure; PlotXY(FilterLFP(Restrict(lfpd,interval+[-1 1]*30),'passband',[0 9],'order',8));
xlim(9489.25 + [0 5]);
PlotXY(FilterLFP(Restrict(lfpd,interval+[-1 1]*30),'passband',[0 9],'order',8));
phase = Phase(filtered, deltas(:,2));
clf
PlotColorMap(Smooth(Shrink(sortby(h,deltas(:,5)-deltas(:,6)),100,1),[0 1]),'x',ht)
PlotColorMap(Smooth(Shrink(sortby(h,phase(:,1)),100,1),[0 1]),'x',ht)
PlotColorMap(Smooth(Shrink(sortby(h,phase(:,2)),100,1),[0 1]),'x',ht)
size(phase)
size(deltas)
PlotColorMap(Smooth(Shrink(sortby(h,phase(:,2)),1000,1),[0 1]),'x',ht)
PlotColorMap(Smooth(Shrink(sortby(h,phase(:,2)),500,1),[0 1]),'x',ht)
[hl,ht] = PETH(lfpd(:,1),deltas(:,2));
PlotColorMap(Smooth(Shrink(sortby(hl,phase(:,2)),500,1),[0 1]),'x',ht)
PlotColorMap(Smooth(Shrink(sortby(hl,phase(:,2)),100,1),[0 1]),'x',ht)
size(hl)
plot(hl(1,:))
plot(hl(3,:))
[hl,ht] = PETH(lfpd,deltas(:,2));
PlotColorMap(Smooth(Shrink(sortby(hl,phase(:,2)),100,1),[0 1]),'x',ht)
PlotColorMap(Smooth(Shrink(sortby(hl,wrap(phase(:,2))),100,1),[0 1]),'x',ht)
hist(phase(:,2))
hist(phase(:,2),100)
hist([phase(:,2);phase(:,2)+2*pi],100)
hist([phase(:,2);phase(:,2)+2*pi]/pi,100)
subplot(1,2,1);
[q,qt ] = hist([phase(:,2);phase(:,2)+2*pi]/pi,100)
[q,qt ] = hist([phase(:,2);phase(:,2)+2*pi]/pi,100);
bar(qt,q);
bar(q,qt);
barh)qt,q);
barh(qt,q);
barh(qt,q); set(gca,'xdir','reverse');
[q,qt ] = hist([phase(:,2)]/pi,100);
barh(qt,q); set(gca,'xdir','reverse');
[q,qt ] = hist([phase(:,2)]/pi,80);
barh(qt,q); set(gca,'xdir','reverse');
[q,qt ] = hist([phase(:,2)]/pi,50);
barh(qt,q); set(gca,'xdir','reverse');
subplot(1,2,2);
b = Bin(phase(:,2),50);
for k=1:50, qq(i,:) = nanmean(h(b==i,:)); end
PlotColorMap(qq);
for k=1:50, qq(k,:) = nanmean(h(b==k,:)); end
PlotColorMap(qq);
for k=1:50, qq(k,:) = nanmean(hl(b==k,:)); end
clear qq
for k=1:50, qq(k,:) = nanmean(hl(b==k,:)); end
PlotColorMap(qq);
PlotColorMap(qq,'x',qt);
PlotColorMap(qq,'x',ht);
PlotHVLines(0,'v','k--');
mmax(phase(:,2))
inphase = wrap(phase(:,2));
mmax(inphase(:,2))
mmax(inphase)
inphase = abs(wrap(phase(:,2)))<(pi/2);
Portion(inphase);
Portion(inphase)
[q,qt ] = hist([phase(inphase,2)]/pi,qt);
hold all
barh(qt,q); set(gca,'xdir','reverse');
figure; semplot(ht,h(inphase,:));
size(ht)
size(h)
[h,ht] = PETH(pfc(:,1),deltas(:,2));
semplot(ht,h(inphase,:));
semplot(ht,h(~inphase,:));
clf
[h,ht] = PETH(r(:,1),deltas(:,2));
semplot(ht,h(inphase,:));
semplot(ht,h(inphase,:),'r');
semplot(ht,h(~inphase,:));
find(InIntervals(deltas(:,2),9489.25 + [0 5]))
[find(InIntervals(deltas(:,2),9489.25 + [0 5])) inphase(ans)]
phase = Phase(filtered, deltas(:,2)+0.06);
semplot(ht,h(~inphase,:));
inphase = abs(wrap(phase(:,2)))<(pi/2);
semplot(ht,h(~inphase,:));
semplot(ht,h(inphase,:),'r'
semplot(ht,h(inphase,:),'r')
clf
semplot(ht,h(inphase,:),'r')
semplot(ht,h(~inphase,:));
clf
[h,ht] = PETH(r(:,1),deltas(:,2));
PlotColorMap(Smooth(Shrink(sortby(hl,wrap(phase(:,2))),100,1),[0 1]),'x',ht)
PlotColorMap(Smooth(Shrink(sortby(h,wrap(phase(:,2))),100,1),[0 1]),'x',ht)
quantile(wrap(phase(:,2)),1/2)
subplot(1,3,1);
PlotColorMap(Smooth(Shrink(sortby(h,wrap(phase(:,2))),100,1),[0 1]),'x',ht)
subplot(1,3,2);
PlotColorMap(Smooth(Shrink(sortby(hl,wrap(phase(:,2))),100,1),[0 1]),'x',ht)
[h,ht] = PETH(pfc(:,1),deltas(:,2));
subplot(1,3,3);
PlotColorMap(Smooth(Shrink(sortby(h,wrap(phase(:,2))),100,1),[0 1]),'x',ht)
PlotHVLines(0,'v','k--');
PlotHVLines(-0.1,'v','k--');
find(InIntervals(deltas(:,2),9509.5 + [0 5]))
inphase(ans)
phase(find(InIntervals(deltas(:,2),9509.5 + [0 5])),2)
xlim(9509.5 + [0 5])
figure; PlotXY(Restrict(filtered,9509.5 + [0 5]));
hold all
hold off
PlotXY(Phase(Restrict(filtered,9509.5 + [0 5])))
PlotHVLines(Restrict(deltas(:,2),xlim),'r--');
PlotHVLines([1 2]*pi,'h','k--');
PlotHVLines([0.5:0.5:2]*pi,'h','k--');
temp = 0.01:0.01:10;
figure; plot(temp);
plot(temp,sin(temp));
hold all
PlotXY(Phase([temp,sin(temp)]));
Phase([temp,sin(temp)]);
Phase([temp(:),sin(temp(:))]);
PlotXY(Phase([temp(:),sin(temp(:))]));
which z2
clf
this = Restrict(filtered,9509.5 + [0 5])));
this = Restrict(filtered,9509.5 + [0 5]));
this = Restrict(filtered,9509.5 + [0 5]);
PlotXY(this);
PlotXY(z2(this));
hold all
PlotXY(z2(Phase(this)));
PlotHVLines(Restrict(deltas(:,2),xlim),'k--','linewidth',2);
PlotHVLines(Restrict(deltas(:,2),xlim)-0.06,'r--','linewidth',2);
PlotColorMap(Smooth(Shrink(sortby(inphase,wrap(phase(:,2))),100,1),[0 1]),'x',ht)
xlim(9524.5 + [0 5])
PlotXY(Restrict(filtered,interval+[0 1]*120));
PlotHVLines(Restrict(deltas(:,2),xlim),'r--','linewidth',2);
find(InIntervals(deltas(:,2),xlim))
inphase(ans)
deltas = deltas0(s>3 & ~bad,:); % these thresholds should be manually refined for each session
PlotHVLines(Restrict(deltas(:,2),xlim),'k--','linewidth',2);
SaveCustomEvents('temp.te2.evt',deltas(:,1:3),{'putative delta start','putative delta peak','putative delta stop'});
deltas = deltas0(s>3 & ~bad,:); % these thresholds should be manually refined for each session
SaveCustomEvents('temp.te2.evt',deltas(:,1:3),{'putative delta start','putative delta peak','putative delta stop'});
phase = Phase(filtered, deltas(:,2)+0.06);
inphase = abs(wrap(phase(:,2)))<(pi/2);
figure; semplot(ht,hl(inphase,:));
semplot(ht,hl(~inphase,:));
clf
semplot(ht,hr(inphase,:));
[hr,ht] = PETH(r(:,1),deltas(:,2));
clf
semplot(ht,hr(inphase,:));
semplot(ht,hr(~inphase,:));
semplot(ht,hr(inphase,:),'r');
clf
phase = Phase(filtered, deltas(:,2));
inphase = abs(wrap(phase(:,2)))<(pi/2);
semplot(ht,hr(inphase,:),'r');
semplot(ht,hr(~inphase,:),'k');
deltaWaves.timestamps = deltas(:,[1 3]); deltaWaves.peaks = deltas(:,2);  deltaWaves.peakNormedPower = deltas(:,5); deltaWaves.detectorName = ['using difference between channels 96 and 109 (+1), peak-trough>3 & ~inAwake & corresponds to peak in 109 signal'];
save(fullfile(basepath,[basename '.deltaWaves.events.mat']),'deltaWaves');
deltas = deltas(inphase,:); % these thresholds should be manually refined for each session
deltaWaves.timestamps = deltas(:,[1 3]); deltaWaves.peaks = deltas(:,2);  deltaWaves.peakNormedPower = deltas(:,5); deltaWaves.detectorName = ['using difference between channels 96 and 109 (+1), peak-trough>3 & ~inAwake & corresponds to peak in 109 signal'];
save(fullfile(basepath,[basename '.deltaWaves.events.mat']),'deltaWaves');
close all
PETH(r(:,1),deltas(:,2))
PETH(pfc(:,1),deltas(:,2))
hist(Restrict(pfc(:,1),9529+[0 5]),100);
[q,qt] = hist(pfc(:,1),t(1):0.05:t(end));
plot(qt,q);
xlim(9529 + [0 5])
bar(qt,q);
xlim(9529 + [0 5])
bins = [0;qt(:)]; bins = [bins(1:end-1) bins(2:end)];
[~,b] = InIntervals(lfpd(:,1),bins);
s = Accumulate(b,lfpd(b>0,2),'mode','mean');
s = Accumulate(b(b>0),lfpd(b>0,2),'mode','mean');
size(s)
size(q)
[x,xt] = xcorr(s(:),q(:),1000);
plot(xt,x);
bar(qt,q);
xlim(9529 + [0 5])
nancorr(s(:),q(:))
phases = Phase(filtered,qt);
open pahses
open phases
b = Bin(phases(:,2),50);
for k=1:50, qqq(k,:) = nanmean(q(b==k,:)); end
size(b)
size(q)
clear qqq
for k=1:50, qqq(k) = nanmean(q(b==k)); end
figure; plot(qqq);
plot([qqq qqq]
plot([qqq qqq])
mmax(phases(:,2))
figure; plot(qqq);
b = Bin(wrap(phases(:,2)),50);
for k=1:50, qqq(k) = nanmean(q(b==k)); end
plot([qqq qqq])
mmax(phases(b==20,2))
inverted = lfpd; inverted(:,2) = -lfpd(:,2);
clean = Restrict(inverted,SubtractIntervals([0 Inf],badIntervals));
deltas0 = FindDeltaPeaks(clean);
45
subplot(1,3,3);
PlotColorMap(Smooth(Shrink(sortby(h,wrap(phase(:,2))),100,1),[0 1]),'x',ht)
PlotHVLines(0,'v','k--'); PlotHVLines(-0.1,'v','k--');
figure; PETH(pfc(:,1),deltas(:,2))
figure(2); clf
subplot(1,3,1);
PlotColorMap(Smooth(Shrink(sortby(hr,ss),100,1),[0 1]),'x',ht)
PlotHVLines(0,'v','k--'); PlotHVLines(-0.1,'v','k--');
subplot(1,3,2);
PlotColorMap(Smooth(Shrink(sortby(hl,ss),100,1),[0 1]),'x',ht)
PlotHVLines(0,'v','k--'); PlotHVLines(-0.1,'v','k--');
subplot(1,3,3);
PlotColorMap(Smooth(Shrink(sortby(h,ss),100,1),[0 1]),'x',ht)
PlotHVLines(0,'v','k--'); PlotHVLines(-0.1,'v','k--');
ss = deltas(:,5)-deltas(:,6);
figure(2); clf
subplot(1,3,1);
PlotColorMap(Smooth(Shrink(sortby(hr,ss),100,1),[0 1]),'x',ht)
PlotHVLines(0,'v','k--'); PlotHVLines(-0.1,'v','k--');
subplot(1,3,2);
PlotColorMap(Smooth(Shrink(sortby(hl,ss),100,1),[0 1]),'x',ht)
PlotHVLines(0,'v','k--'); PlotHVLines(-0.1,'v','k--');
subplot(1,3,3);
PlotColorMap(Smooth(Shrink(sortby(h,ss),100,1),[0 1]),'x',ht)
PlotHVLines(0,'v','k--'); PlotHVLines(-0.1,'v','k--');
5
PETH(pfc(:,1),deltas(:,2))
xlim(9537 + [0 5])
PlotHVLines(Restrict(deltas(:,2),xlim),'k--','linewidth',2);
PlotHVLines(Restrict(deltas(:,1),xlim),'g');
PlotHVLines(Restrict(deltas(:,3),xlim),'b--');
clean = Restrict(inverted,SubtractIntervals([0 Inf],badIntervals));
deltas0 = FindDeltaPeaks(clean);
s = deltas0(:,5)-deltas0(:,6);
bad = InIntervals(deltas0(:,2),SleepState.ints.WAKEstate);
deltas = deltas0(s>3 & ~bad,:); % these thresholds should be manually refined for each session
1
PlotHVLines(Restrict(deltas(:,3),xlim),'b--');
s = deltas0(:,5)-deltas0(:,6);
bad = InIntervals(deltas0(:,2),SleepState.ints.WAKEstate);
deltas = deltas0(s>3 & ~bad,:); % these thresholds should be manually refined for each session
xlim(9517 + [0 5])
PlotHVLines(Restrict(deltas(:,3),xlim),'b--');
PlotHVLines(Restrict(deltas(:,1),xlim),'g');
PlotHVLines(Restrict(deltas(:,2),xlim),'k--','linewidth',2);
difference = [t lfpd(:,2)-lfps(:,2)];
clean = Restrict(difference,SubtractIntervals([0 Inf],badIntervals));
deltas0 = FindDeltaPeaks(clean);
load('day3.SleepState.states.mat')
bad = InIntervals(deltas0(:,2),SleepState.ints.WAKEstate);
s = deltas0(:,5)-deltas0(:,6); % score delta wave strength
deltas = deltas0(s>3 & ~bad,:);
deltas(:,1:3) = deltas(:,1:3)-0.06; % make it correspond to pfc unit silence
% only allow putative delta waves occurring in the correct phase
filtered = FilterLFP(lfpd,'passband',[0 9],'order',8);
phase = Phase(filtered, deltas(:,2));
inphase = abs(wrap(phase(:,2)))>(pi/2); % "peak" should be closer to signal peak than to trough
deltas = deltas(inphase,:); % these thresholds should be manually refined for each session
deltaWaves.timestamps = deltas(:,[1 3]); deltaWaves.peaks = deltas(:,2);  deltaWaves.peakNormedPower = deltas(:,5); deltaWaves.detectorName = ['using difference between channels 96 and 109 (+1), peak-trough>3 & ~inAwake & corresponds to trough (inversion!) in 109 signal'];
figure; PETH(pfc(:,1),deltas(:,2))
PETH(r(:,1),deltas(:,2))
save(fullfile(basepath,[basename '.deltaWaves.events.mat']),'deltaWaves');
sleep = ConsolidateIntervals(MergePoints.timestamps(cellfun(@(x) any(strfind(lower(x),'sleep')),MergePoints.foldernames),:))
load('day3.MergePoints.events.mat')
sleep = ConsolidateIntervals(MergePoints.timestamps(cellfun(@(x) any(strfind(lower(x),'sleep')),MergePoints.foldernames),:))
presentationOrder = reshape(blockMatrix', [], 1); % the final order is a single vector (we don't care about blocks any more), just "left" vs "right" trials one after the other​
preSleep = SubtractIntervals(sleep(1,:), SubtractIntervals([0 Inf],sws));
sws = SleepState.ints.NREMstate;
preSleep = SubtractIntervals(sleep(1,:), SubtractIntervals([0 Inf],sws));
postSleep = SubtractIntervals(sleep(2:end,:), SubtractIntervals([0 Inf],sws));
pre = InIntervals(r,preSleep);
post = InIntervals(r,postSleep);
[hr,ht] = PETH(deltas(:,2),r(:,2));
[hr,ht] = PETH(deltas(:,2),r(:,2)); ht = fliplr(ht);
semplot(ht,hr);
semplot(ht,hr(pre,:));
clf
semplot(ht,hr(pre,:));
semplot(ht,hr(post,:),'r');
ht = fliplr(ht);
clf
semplot(ht,hr(pre,:));
semplot(ht,hr(post,:),'r');
save(fullfile(basepath,[basename '.deltaWaves.events.mat']),'deltaWaves');
cd N:\OJRproject\OJR42\day12
SaveCustomEvents('deltas.del.evt',deltas(:,1:3),{'deltas start','delta peak','deltas stop'});
cd N:\OJRproject\OJR42\day12
cd N:\OJRproject\OJR42\day11
clear all
close all
channel = 79;
lfp = GetAyaLFP(channel);
display(['loaded! @' num2str(toc)]);
figure; PlotXY(lfp);
PlotIntervals(badIntervals,'color','r');
45
[clean,~,badIntervals] = CleanLFP(lfp,'thresholds',[3 7]);
PlotIntervals(badIntervals,'color','r');
clf
PlotXY(clean)
[clean,~,badIntervals] = CleanLFP(lfp,'thresholds',[5 10]);
figure; PlotXY(lfp);
PlotIntervals(badIntervals,'color','r');
figure; PlotXY(clean)
[clean,~,badIntervals] = CleanLFP(lfp,'thresholds',[5 10]);
figure; plot(t,z);
PlotHVLines(3,'h','k--')
PlotHVLines(4,'h','k--')
xlim([-10 10]+4827);
[clean,~,badIntervals] = CleanLFP(lfp,'thresholds',[3 10]);
[clean,~,badIntervals] = CleanLFP(lfp,'thresholds',[4 10]);
PlotXY(clean);
deltas0 = FindDeltaPeaks(clean);
%     deltas = deltas0(deltas0(:,5)>1.5 & deltas0(:,5)<4.5,:); % these thresholds should be manually refined for each session
deltas = deltas0(deltas0(:,5)-deltas0(:,6)>3,:); % these thresholds should be manually refined for each session
deltaWaves.timestamps = deltas(:,[1 3]); deltaWaves.peaks = deltas(:,2);  deltaWaves.peakNormedPower = deltas(:,5); deltaWaves.detectorName = ['channel ' num2str(channel) '(+1), CleanLFP, FindDeltaPeaks, peak-trough>3'];
% save(fullfile(basepath,[basename '.deltaWaves.events.mat']),'deltaWaves');
basepath = pwd
basename = basenameFromBasepath(basepath);
load(fullfile(basepath,[basename '.cell_metrics.cellinfo.mat']),'cell_metrics');
nNeurons = length(cell_metrics.spikes.times);
regions = cell_metrics.brainRegion;
regionNames = unique(regions);
regionCell = zeros(length(regions),1);
for i=1:length(regionNames)
regionCell(strcmp(regions,regionNames{i})) = i;
end
spikesCell = cell_metrics.spikes.times';
PFCindex = [find(strcmp(regionNames,'ILA')) find(strcmp(regionNames,'PFC'))];
HPCindex = find(strcmp(regionNames,'CA1'));
pfc = []; if ~isempty(PFCindex), pfc = sortrows(Group(spikesCell{regionCell==PFCindex})); end
load('day11.ripples.events.mat')
r = [ripples.timestamps(:,1) ripples.peaks ripples.timestamps(:,2) ripples.RipMax ripples.SwMax];
r = [ripples.timestamps(:,1) ripples.peaks ripples.timestamps(:,2)];
load('day11.SleepState.states.mat')
figurr; hist(s,100);
figure; hist(s,100);
Portion(s>3)
deltaWaves.timestamps = deltas(:,[1 3]); deltaWaves.peaks = deltas(:,2);  deltaWaves.peakNormedPower = deltas(:,5); deltaWaves.detectorName = ['channel ' num2str(channel) '(+1), CleanLFP, FindDeltaPeaks, peak-trough>3'];
deltaWaves.detectorName
save(fullfile(basepath,[basename '.deltaWaves.events.mat']),'deltaWaves');
load('day11.MergePoints.events.mat')
sleep = ConsolidateIntervals(MergePoints.timestamps(cellfun(@(x) any(strfind(lower(x),'sleep')),MergePoints.foldernames),:))
sws = SleepState.ints.NREMstate;
preSleep = SubtractIntervals(sleep(1,:), SubtractIntervals([0 Inf],sws));
postSleep = SubtractIntervals(sleep(2:end,:), SubtractIntervals([0 Inf],sws));
[h,ht] = PETH(deltas(:,2),r(:,2));
semplot(ht,hr(post,:),'r');
[hr,ht] = PETH(deltas(:,2),r(:,2));
semplot(ht,hr(post,:),'r');
post = InIntervals(r,postSleep);
pre = InIntervals(r,preSleep);
semplot(ht,hr(post,:),'r');
semplot(ht,hr(pre,:));
clear all
cd N:\OJRproject\OJR33\day7
cd N:\OJRproject\OJR33\day9
cd N:\OJRproject\OJR34\day2
cd N:\OJRproject\OJR34\day7
open OJR_shortTraining.batch
cd N:\OJRproject\OJR33\day10
cd N:\OJRproject\OJR33\day11
cd N:\OJRproject\OJR34\day1
cd N:\OJRproject\OJR34\day8
cd N:\OJRproject\OJR34\day9
batchO = StartBatch(@BatchLoadHPCPFCStimData,'OJR_shortTraining_optoRipples_delayedPFC.batch');
Xo = get(batchO,'UserData');
batchS = StartBatch(@BatchLoadHPCPFCStimData,'OJR_shortTraining.batch');
Xs = get(batchS,'UserData');
batchL = StartBatch(@BatchLoadHPCPFCStimData,'OJR_longTraining.batch');
batchS = StartBatch(@BatchLoadHPCPFCData,'OJR_shortTraining.batch');
batchL = StartBatch(@BatchLoadHPCPFCData,'OJR_longTraining.batch');
Xs = get(batchS,'UserData');
Xl = get(batchL,'UserData');
XX = {Xs,Xl,Xo};
354
k
i
size(r)
size(deltas)
X{i,end}
k=2; PlotColorMap([mm0{k}  mm1{k}]);
k=2; PlotColorMap(zscore([mm0{k}  mm1{k}],[],2));
subplot(1,3,1);
PlotColorMap(mm0{k});
subplot(1,3,2);
PlotColorMap(mm1{k});
subplot(1,3,3);
PlotColorMap(mm1{k}-mm0{k});
subplot(1,3,1);
PlotColorMap(mm0{k},'x',ht);
subplot(1,3,2);
PlotColorMap(mm1{k},'x',ht);
subplot(1,3,3);
PlotColorMap(mm1{k}-mm0{k},'x',ht);
clim
clim([-1 1]*min(abs(clim)))
56
2
45
open mm1
values = nanmean(mm1{k}(:,InIntervals(ht,[0 0.2])),2);
open values
values = [nanmean(mm0{k}(:,InIntervals(ht,[0 0.2])),2) nanmean(mm1{k}(:,InIntervals(ht,[0 0.2])),2)];
figure; plot(values)
values(1,:) = nan;
plot(values)
anovabar(values)
k=1;
values_s = [nanmean(mm0{k}(:,InIntervals(ht,[0 0.2])),2) nanmean(mm1{k}(:,InIntervals(ht,[0 0.2])),2)];
g = Group(diff(values_s,[],2),diff(values,[],2));
anovabar(g(:,1),g(:,2));
figure; k=2; PlotColorMap(nanzscore([mm0{k} mm1{k}],[],2))
k=2; PlotColorMap(nanzscore([mm0{k} mm1{k}],[],2))
k=2; PlotColorMap(zscore([mm0{k} mm1{k}],[],2))
k=1; PlotColorMap(zscore([mm0{k} mm1{k}],[],2))
k=3; PlotColorMap(zscore([mm0{k} mm1{k}],[],2))
k=3; PlotColorMap(nanzscore([mm0{k} mm1{k}],[],2))
k=2; PlotColorMap(nanzscore([mm0{k} mm1{k}],[],2))
k=1; PlotColorMap(nanzscore([mm0{k} mm1{k}],[],2))
k=2; PlotColorMap(nanzscore([mm0{k} mm1{k}],[],2))
open Xl
k=3; PlotColorMap(nanzscore([mm0{k} mm1{k}],[],2))
open Xo
load('day10.channelinfo.ripples.mat')
PlotColorMap(mm1{k}-mm0{k},'x',ht);
k=1
PlotColorMap(mm1{k}-mm0{k},'x',ht);
k=2; PlotColorMap(nanzscore([mm0{k} mm1{k}],[],2))
k=1; PlotColorMap(nanzscore([mm0{k} mm1{k}],[],2))
k=3; PlotColorMap(nanzscore([mm0{k} mm1{k}],[],2))
clim
clim([-1 1]*min(abs(clim)))
ns(i,1) = sum(InIntervals(deltas(:,2),[r(:,2)+0.05 r(:,2)+0.25]))
clear ns
465
open ns
q = ns';
open q
q = Group(ns{:});
q = Group(ns{:}); q = [q(:,2)./q(:,1) q(:,end)];
anovabar(q(:,1),q(:,2));
clf
anovabar(q(:,1),q(:,2));
kruskalbar(q(:,1),q(:,2));
kruskalbar(q(:,1)-1,q(:,2));
anovabar(q(:,1)-1,q(:,2));
q = Group(ns{:}); q = [q(:,2)./q(:,1) q(:,end)];
q = Group(ns{:});
open q
45
q = Group(ns{:});
anovabar(q(:,1)-1,q(:,2));
clf
anovabar(q(:,1:2))
anovabar(q(:,1:2),q(:,end))
anovabar(diff(q(:,1:2),[],2),q(:,end))
q = Group(ns{:}); q = [q(:,2)./q(:,1) q(:,end)];
q = Group(ns{:}); q = [q(:,2)./q(:,1) q(:,end)]; q(q==Inf) = nan;
anovabar(q(:,1:2),q(:,end))
anovabar(q(:,1),q(:,end))
anovabar(q(:,1)-1,q(:,end))
plot(q(:,end),q(:,1)-1)
plot(q(:,end),q(:,1)-1,'.')
plot(q(:,end),q(:,1)-1,'.'); xlim([0 4]);
signrank([[1.07633660250121;1.12170358974211;1.53428010748644;1.16719606641281;0.910823515012342]])
qq = Group(ns{:}); q = [qq(:,2)./qq(:,1) qq(:,end)]; q(q==Inf) = nan;
open qq
open Xo
[pfc,hpc,deltas,r,sleep,session,MergePoints,clickers,basepath] = BatchLoadHPCPFCData('N:\OJRproject\OJR34\day2')
[pfc,hpc,deltas,r,sleep,session,MergePoints,clickers,basepath] = BatchLoadHPCPFCData('N:\OJRproject\OJR34\day2');
Xo(3,:) = [pfc,hpc,deltas,r,sleep,session,MergePoints,clickers,basepath];
Xo(3,:) = {pfc,hpc,deltas,r,sleep,session,MergePoints,clickers,basepath};
k
X = XX{k};
open X
i=3
deltas = X{i,3};
r = X{i,4};
sleep = X{i,5};
%         preSleep=sleep(1,:);
%         postSleep = sleep(2:end,:);
sws = [0 Inf];
basepath = X{i,end};
basename = basenameFromBasepath(basepath);
load(fullfile(basepath,[basename '.SleepState.states.mat']));
sws = SleepState.ints.NREMstate;
preSleep = SubtractIntervals(sleep(1,:), SubtractIntervals([0 Inf],sws));
postSleep = SubtractIntervals(sleep(2:end,:), SubtractIntervals([0 Inf],sws));
pre = InIntervals(r,preSleep);
post = InIntervals(r,postSleep);
sum(InIntervals(deltas(:,2),[r(pre,2)+0.05 r(pre,2)+0.25]))./dur(preSleep)
sum(InIntervals(deltas(:,2),[r(post,2)+0.05 r(post,2)+0.25]))./dur(postSleep)
ns{k}(i,1) = sum(InIntervals(deltas(:,2),[r(pre,2)+0.05 r(pre,2)+0.25]))./dur(preSleep);
ns{k}(i,2) = sum(InIntervals(deltas(:,2),[r(post,2)+0.05 r(post,2)+0.25]))./dur(postSleep);
qq = Group(ns{:}); q = [qq(:,2)./qq(:,1) qq(:,end)]; q(q==Inf) = nan;
open qq
batchS = StartBatch(@BatchLoadHPCPFCData,'OJR_shortTraining.batch');
Xs = get(batchS,'UserData'); % short
XX = {Xs,Xl,Xo};
qq = Group(ns{:}); q = [qq(:,2)./qq(:,1) qq(:,end)]; q(q==Inf) = nan;
open q
plot(q(:,end),q(:,1)-1,'.'); xlim([0 4]);
hold all
anovabar(q(:,1)-1,q(:,end))
45
qq = Group(ns{:}); q = [qq(:,2)./qq(:,1) qq(:,end)]; q(q==Inf) = nan;
anovabar(q(:,1)-1,q(:,end))
hold all
plot(q(:,end),q(:,1)-1,'.'); xlim([0 4]);
plot(q(:,end),q(:,1)-1,'.','markersize',20); xlim([0 4]);
kruskalbar(q(:,1)-1,q(:,end))
clf
kruskalbar(q(:,1)-1,q(:,end))
open q
open Xs
kruskalbar(q(:,1),q(:,end))
0.0814/0.0525
0.1014/0.0691
cd {'N:\OJRproject\OJR34\day9'}
cd('N:\OJRproject\OJR34\day9')
k=1
X = XX{k};
i=5
deltas = X{i,3};
r = X{i,4};
sleep = X{i,5};
%         preSleep=sleep(1,:);
%         postSleep = sleep(2:end,:);
sws = [0 Inf];
basepath = X{i,end};
try
basename = basenameFromBasepath(basepath);
load(fullfile(basepath,[basename '.SleepState.states.mat']));
sws = SleepState.ints.NREMstate;
catch
['No sleep scoring for session ' num2str(i) ': ' basepath]
%             continue
end
preSleep = SubtractIntervals(sleep(1,:), SubtractIntervals([0 Inf],sws));
postSleep = SubtractIntervals(sleep(2:end,:), SubtractIntervals([0 Inf],sws));
pre = InIntervals(r,preSleep);
post = InIntervals(r,postSleep);
sum(InIntervals(deltas(:,2),[r(pre,2)+0.05 r(pre,2)+0.25]))
sum(InIntervals(deltas(:,2),[r(pre,2)+0.05 r(pre,2)+0.25]))./dur(preSleep)
sum(InIntervals(deltas(:,2),[r(post,2)+0.05 r(post,2)+0.25]))
sum(InIntervals(deltas(:,2),[r(post,2)+0.05 r(post,2)+0.25]))./dur(postSleep)
dur(preSleep)
dur(postSleep)
figure; hist(r(:,2),10000);
PlotIntervals(preSleep,'color','b');
PlotIntervals(postSleep,'color','r');
figure; hist(deltas(:,2),10000);
PlotIntervals(preSleep,'color','b');
PlotIntervals(postSleep,'color','r');
clf
PETH(r(:,2),deltas(:,2));
PETH(r(:,2),Restrict(deltas(:,2),preSleep));
hold all
PETH(r(:,2),Restrict(deltas(:,2),postSleep));
clf
hist(r(:,2),100);
hist(deltas(:,2),100);
hist(r(:,2),100);
PlotIntervals(preSleep,'color','b');
load('day8.SleepState.states.mat')
sws = SleepState.ints.NREMstate;
figure; PlotIntervals(sws,'color','b');
load('day8.deltaWaves.events.mat')
load('day8.ripples.events.mat')
r = [ripples.timestamps(:,1) ripples.peaks ripples.timestamps(:,2)];
deltas = X{i,3};
i=4
deltas = X{i,3};
clf
PETH(r(:,2),Restrict(deltas(:,2),postSleep));
load('day8.MergePoints.events.mat')
sleep = X{i,5};
postSleep = SubtractIntervals(sleep(2:end,:), SubtractIntervals([0 Inf],sws));
preSleep = SubtractIntervals(sleep(1,:), SubtractIntervals([0 Inf],sws));
PETH(r(:,2),Restrict(deltas(:,2),postSleep));
hold all
PETH(r(:,2),Restrict(deltas(:,2),preSleep));
load('objScore.mat')
% good ripple effect. dis index = 29, preferece = -6
open q
cd {'N:\OJRproject\OJR34\day1'}
load('objScore.mat')
i=2
deltas = X{i,3};
r = X{i,4};
sleep = X{i,5};
%         preSleep=sleep(1,:);
%         postSleep = sleep(2:end,:);
sws = [0 Inf];
basepath = X{i,end};
try
basename = basenameFromBasepath(basepath);
load(fullfile(basepath,[basename '.SleepState.states.mat']));
sws = SleepState.ints.NREMstate;
catch
['No sleep scoring for session ' num2str(i) ': ' basepath]
%             continue
end
preSleep = SubtractIntervals(sleep(1,:), SubtractIntervals([0 Inf],sws));
postSleep = SubtractIntervals(sleep(2:end,:), SubtractIntervals([0 Inf],sws));
pre = InIntervals(r,preSleep);
post = InIntervals(r,postSleep);
sum(InIntervals(deltas(:,2),[r(pre,2)+0.05 r(pre,2)+0.25]))
sum(InIntervals(deltas(:,2),[r(pre,2)+0.05 r(pre,2)+0.25]))./dur(preSleep)
sum(InIntervals(deltas(:,2),[r(post,2)+0.05 r(post,2)+0.25]))
sum(InIntervals(deltas(:,2),[r(post,2)+0.05 r(post,2)+0.25]))./dur(postSleep)
clf
hist(r(:,2),100);
PlotIntervals(preSleep,'color','b');
sleep
diff(sleep,[],2)/3600
PlotIntervals(postSleep,'color','r');
size(r)
dur(preSleep)
preSleep = SubtractIntervals([0 2000], SubtractIntervals([0 Inf],sws));
sum(InIntervals(deltas(:,2),[r(post,2)+0.05 r(post,2)+0.25]))./dur(postSleep)
sum(InIntervals(deltas(:,2),[r(pre,2)+0.05 r(pre,2)+0.25]))./dur(preSleep)
pre = InIntervals(r,preSleep);
sum(InIntervals(deltas(:,2),[r(pre,2)+0.05 r(pre,2)+0.25]))./dur(preSleep)
preSleep = SubtractIntervals(sleep(1,:), SubtractIntervals([0 Inf],sws));
pre = InIntervals(r,preSleep);
sum(InIntervals(deltas(:,2),[r(pre,2)+0.05 r(pre,2)+0.25]))./dur(preSleep)
preSleep = SubtractIntervals([4000 9000], SubtractIntervals([0 Inf],sws));
pre = InIntervals(r,preSleep);
sum(InIntervals(deltas(:,2),[r(pre,2)+0.05 r(pre,2)+0.25]))./dur(preSleep)
bins = Bins(0,MergePoints.timestamps(end,end),0.1);
count = CountInIntervals(Restrict(deltas(:,2),[r(:,2)+0.05 r(:,2)+0.25]));
count = CountInIntervals(Restrict(deltas(:,2),[r(:,2)+0.05 r(:,2)+0.25]),bins);
plot(count)
clf
t = mean(bins,2);
plot(t,count);
plot(t,Smooth(count,100))
plot(t,Smooth(count,1000))
hm = [t count];
PlotXY(hm(:,1),Smooth(hm(:,2),200));
in = InIntervals(t,sws);
PlotXY(hm(in,1),Smooth(hm(in,2),200));
PlotXY(hm(in,1),Smooth(hm(in,2),1000));
PlotIntervals(preSleep,'color','b');
preSleep = SubtractIntervals(sleep(1,:), SubtractIntervals([0 Inf],sws));
PlotIntervals(preSleep,'color','b');
clf
PlotIntervals(preSleep,'color','b');
PlotXY(hm(in,1),Smooth(hm(in,2),1000));
clf
PlotXY(hm(in,1),Smooth(hm(in,2),1000));
PlotXY(hm(:,1),Smooth(hm(:,2),1000));
PlotIntervals(preSleep,'color','b');
PlotIntervals(postSleep,'color','r');
PlotXY(hm(:,1),Smooth(hm(:,2),100));
clf
PlotXY(hm(:,1),Smooth(hm(:,2),100));
PlotIntervals(postSleep,'color','r');
clf
PlotXY(Restrict(hm,sws,'shift','on'))
PlotXY(Smooth(Restrict(hm,sws,'shift','on'),[5 0]))
PlotXY(Smooth(Restrict(hm,sws,'shift','on'),[100 0]))
PlotXY(Smooth(Restrict(hm,sws,'shift','on'),[1000 0]))
dur(preSleep)
PlotHVLines(5036,'v','k--')
preSleep = SubtractIntervals(sleep(1,:), SubtractIntervals([0 Inf],sws));
preSleep = SubtractIntervals(preSleep,[Unshift(3600,preSleep) Inf]);
dur(preSleep)
pre = InIntervals(r,preSleep);
sum(InIntervals(deltas(:,2),[r(pre,2)+0.05 r(pre,2)+0.25]))./dur(preSleep)
sum(InIntervals(deltas(:,2),[r(post,2)+0.05 r(post,2)+0.25]))./dur(postSleep)
0.0517/   0.0406
dur(postSleep)
preSleep = SubtractIntervals(preSleep,[Unshift(dur(postSleep),preSleep) Inf]);
sum(InIntervals(deltas(:,2),[r(post,2)+0.05 r(post,2)+0.25]))./dur(postSleep)
sum(InIntervals(deltas(:,2),[r(pre,2)+0.05 r(pre,2)+0.25]))./dur(preSleep)
pre = InIntervals(r,preSleep);
post = InIntervals(r,postSleep);
sum(InIntervals(deltas(:,2),[r(pre,2)+0.05 r(pre,2)+0.25]))./dur(preSleep)
preSleep = SubtractIntervals(sleep(1,:), SubtractIntervals([0 Inf],sws));
preSleep = SubtractIntervals(preSleep,[Unshift(dur(postSleep),preSleep) Inf]);
dur(preSleep)
preSleep = SubtractIntervals(sleep(1,:), SubtractIntervals([0 Inf],sws));
preSleep = SubtractIntervals(preSleep,[Unshift(3600,preSleep) Inf]);
pre = InIntervals(r,preSleep);
sum(InIntervals(deltas(:,2),[r(pre,2)+0.05 r(pre,2)+0.25]))./dur(preSleep)
Unshift(938,sws)
figure; PlotXY(hm);
PlotIntervals(preSleep,'color','b');
load('day11.ripples.events.mat')
load('day11.deltaWaves.events.mat')
deltas = [deltaWaves.timestamps(:,1) deltaWaves.peaks deltaWaves.timestamps(:,2)];
SaveCustomEvents('deltas.del.evt',deltas(:,1:3),{'deltas start','delta peak','deltas stop'});
figure; PETH(r(:,2),deltas(:,2))
load('day11.channelinfo.ripples.mat')
load('day10.channelinfo.ripples.mat')
load('day7.channelinfo.ripples.mat')
rippleChannels.Ripple_Channel = 36+1; rippleChannels.Sharpwave_Channel = 9+1; rippleChannels.Noise_Channel = [];
save('day7.channelinfo.ripples.mat','rippleChannels');
ripples = DetectSWR([rippleChannels.Ripple_Channel rippleChannels.Sharpwave_Channel],'basepath',basepath,'saveMat',true,'forceDetect',true,'check',true);
% First, remove points that occur in noisy lfp periods
tl = (1:length(lfp))'/1250; % timestamps for lfp
t = featureTs/1250; % timestamps for candidate ripple events
bad = false(size(t));
[clean,~,badIntervals] = CleanLFP([tl,double(lfp(:,2))],'thresholds',[6 Inf]);
figure; plot(t,abs(d));
% ok ripple @ 10826.95
plot(t,abs(z))
PlotHVLines(8,'h','k--')
[clean,~,badIntervals] = CleanLFP([tl,double(lfp(:,2))],'thresholds',[7 1]);
bad = bad | InIntervals(t,badIntervals);
dbcont
bad = bad | InIntervals(t,badIntervals);
Portion(bad)
idx1 = idx1 & ~bad; idx2 = idx2 & ~bad; % remove bad points from "idx1" and "idx2"
figure(1); % plot a clean figure without these points and check indivudual ripples
clf
plot(swDiffAll(idx1 | idx2),ripPowerAll(idx1 | idx2),'k.','markersize',1); hold on
matrix = [swDiffAll ripPowerAll];
z = bsxfun(@rdivide,bsxfun(@minus,matrix,mean(matrix(~bad,:))),std(matrix(~bad,:)));
scores = nan(size(ripPowerAll,1),1);
while true % click "Ctrl+C" to exit when you feel you're done
figure(1);
colors = Bright(1000);
[x,y,button] = ginput(1);
[~,i] = min(abs(x-swDiffAll.*(1-bad))+abs(y-ripPowerAll.*(1-bad)));
t = featureTs/1250;
tl = (1:length(lfp))'/1250;
figure(2);
clf; interval = t(i) + [-1 1]*0.5; in = tl>interval(1) & tl<interval(end);
subplot(3,1,1);
plot(tl(in) - t(i),lfp(in,:)); title(num2str(t(i)));
hold on; plot(bsxfun(@times,ones(1,2),(Restrict(t,interval) - t(i))),ylim,'k--');
legend('ripple channel','SW channel','noise channel');
subplot(3,1,2);
plot(tl(in) - t(i),double(swDiff(in)));
hold on; plot(bsxfun(@times,ones(1,2),(Restrict(t,interval) - t(i))),ylim,'k--');
subplot(3,1,3);
plot(tl(in) - t(i),double(ripPower0(in)));
hold on; plot(bsxfun(@times,ones(1,2),(Restrict(t,interval) - t(i))),ylim,'k--');
%     x = input( prompt )
commandwindow % select the command window
score = str2double(input('good(1)/bad(0)?','s'));
scores(i) = score; % save scores to optionally save
figure(1);
xlims = xlim; ylims = ylim; clf
% extrapolate the score of each ripple based on the score of the nearest scored ripple
scored = find(~isnan(scores));
distances = sqrt(bsxfun(@minus,z(:,1),z(scored,1)').^2 + bsxfun(@minus,z(:,2),z(~isnan(scores),2)').^2);
[~,nearestNeighbour] = min(distances,[],2);
estimated = scores(scored(nearestNeighbour));
selected = estimated>0.1;
idx1 = selected & ~bad; % final ripples
idx2 = ~selected & ~bad; % non-ripples (yet free from noise as well)
scatter(matrix(~bad,1),matrix(~bad,2),1,estimated(~bad),'filled'); colormap(Bright); clim([0 1]);
xlim(xlims); ylim(ylims);
hold on;
scatter(swDiffAll(scored),ripPowerAll(scored),20,scores(scored)); scatter(swDiffAll(scored),ripPowerAll(scored),15,scores(scored));
clabel('Estimated ripple score');
xlabel('Sharp wave depth'); ylabel('Ripple power');
title({['Manual scoring for ' basepath],['called with channels: ' num2str(Channels) ' (ripple, sharp wave, and (optionally) noise)']});
if rem(i,10)==0
display('saving...');
save(fullfile(basepath,'DetectSWR_manual_scoring.mat'),'t','swDiffAll','ripPowerAll','bad','scores');
end
end
0
1
0
1
0
1
0
% Save your progress!
display('saving...');
save(fullfile(basepath,'DetectSWR_manual_scoring.mat'),'t','swDiffAll','ripPowerAll','idx1','idx2','bad','selected','scores');
saveas(figure(1),fullfile(basepath,'DetectSWR_manual_scoring.fig'));
sum(idx1)
dbcont
batchS = StartBatch(@BatchLoadHPCPFCData,'OJR_shortTraining.batch');
Xs = get(batchS,'UserData'); % short
XX = {Xs,Xl,Xo};
XX = {Xs,Xl,Xo};
q = Group(ns{:}); q = [q(:,2)./q(:,1) q(:,end)]; q(q==Inf) = nan;
qq = Group(ns{:}); q = [qq(:,2)./qq(:,1) qq(:,end)]; q(q==Inf) = nan;
open qq
open q
kruskalbar(q(:,1),q(:,end))
clf
kruskalbar(q(:,1),q(:,end))
kruskalbar(q(:,1)-1,q(:,end))
anovabar(q(:,1)-1,q(:,end))
q(11,1) = nan;
anovabar(q(:,1)-1,q(:,end))
kruskalbar(q(:,1)-1,q(:,end))
anovabar(q(:,1)-1,q(:,end))
open Xl
open Xs
raly
clear all
close all
34
open q
open Xs
anovabar(q(:,1)-1,q(:,end))
kruskalbar(q(:,1)-1,q(:,end))
open mm0
open Xpo
open Xo
batchO = StartBatch(@BatchLoadHPCPFCData,'OJR_shortTraining_optoRipples_delayedPFC.batch');
Xo = get(batchO,'UserData'); % opto
34
batchO = StartBatch(@BatchLoadHPCPFCData,'OJR_shortTraining_optoRipples_delayedPFC.batch');
Xo = get(batchO,'UserData'); % opto
open q
45.
45
batchO = StartBatch(@BatchLoadHPCPFCData,'OJR_shortTraining_optoRipples_delayedPFC.batch');
XX = {Xs,Xl,Xo};
batchO = StartBatch(@BatchLoadHPCPFCData,'OJR_shortTraining_optoRipples_delayedPFC.batch');
Xo = get(batchO,'UserData'); % opto
XX = {Xs,Xl,Xo};
clf
anovabar(q(:,1)-1,q(:,end))
kruskalbar(q(:,1)-1,q(:,end))
anovabar(q(:,1)-1,q(:,end))
45
clims([-1 1]*min(abs(clim)))
PlotXY(q);
PlotXY(q,'.');
plot(q(:,end),q(:,1)-1,'.','markersize',20); xlim([0 4]);
for i=1:size(q,1)
this = q; this(i,1) = nan;
qqq(i,1) = nanmean(this(q(:,2)==q(i,2)));
qqq(i,2) = q(i,2);
end
open qqq
plot(q(:,end),qqq(:,1)-1,'.','markersize',20); xlim([0 4]);
anovabar(qqq(:,1),q(:,end))
anovabar(qqq(:,1)-1,q(:,end))
plot(q(:,end),q(:,1)-1,'.','markersize',20); xlim([0 4]);
plot(q(:,end),q(:,1),'.','markersize',20); xlim([0 4]);
plot(q(:,end),q(:,1)*100,'.','markersize',20); xlim([0 4]);
set(gca,'xtick',1:3,'xticklabel',{'short','long','stim'});
ylabel('Joint ripple-delta occurrence relative to pre-sleep (%)');
set(gca,'box','off','fontsize',15);
PlotHVLines(100,'h','k--')
ylim([50 200]);
raly
SaveFig('M:\home\raly\results\PFC\delta\JointRippleDeltaOccurrence_short_long_opto_conditions')
clf
plot(q(:,end),q(:,1)*100,'.','markersize',20); xlim([0 4]);
set(gca,'xtick',1:3,'xticklabel',{'short','long','stim'});
ylabel('Joint ripple-delta occurrence relative to pre-sleep (%)');
set(gca,'box','off','fontsize',15);
PlotHVLines(100,'h','k--')
ylim([50 200]);
SaveFig('M:\home\raly\results\PFC\delta\JointRippleDeltaOccurrence_short_long_opto_conditions')
set(gca,'xtick',1:3,'xticklabel',{'short training','long training','opto prolongation'});
SaveFig('M:\home\raly\results\PFC\delta\JointRippleDeltaOccurrence_short_long_opto_conditions')
open XX
SaveFig('M:\home\raly\results\PFC\delta\JointRippleDeltaOccurrence_short_long_opto_conditions')
SaveFig('M:\home\raly\results\PFC\delta\DeltaOccurrence_around_ripples_short_long_opto_conditions')
open 'M:\home\raly\results\PFC\delta\JointRippleDeltaOccurrence_short_long_opto_conditions.fig'
5
2
clear all
X = Xl;
5
X = Xs;
54
5
X = Xo;
5
clim
clims([-1 1]*min(abs(clim)))
56
clim
clims([-1 1]*min(abs(clim)))
clims([-5 5]);
z1 = zscore(z1,[],2); z2 = zscore(z2,[],2);
% figure('name',name);
clf
set(gcf,'position',[1 1 1920 1024])
matrices = {z1,z2,z2-z1};
titles = {[name ' in pre-task sleep'],[name ' in post-task sleep'],'(post-pre) difference'};
colors = {'k','r','m'};
for j=1:3
subplot(2,3,j); PlotColorMap(Smooth(matrices{j},smooth),~isnan(matrices{j}),'x',qt); PlotHVLines(0,'v','k--','linewidth',2);
PlotHVLines(lines,'h','w--','linewidth',2);
set(gca,'ytick',mean([[0;lines(1:end-1)] lines],2),'yticklabel',sessions);
title(titles{j});
subplot(2,3,3+j); semplot(qt,matrices{j},colors{j},smooth(end));
if j==1, semplot(qt,matrices{j+1},colors{j+1},smooth(end)); y=ylim; end
ylim(y); PlotHVLines(0,'v','k--','linewidth',2); PlotHVLines(0,'h','k--','linewidth',2); ylim(y);
ylabel('mean +/ s.e.m. (z-units)');
end
clims
drawnow
name = [name ' spikes'];
dd = cell2mat([cellfun(@zscore,ds(:,1),'uniformoutput',0) ds(:,2)]);
nComponents = Accumulate(cell2mat(qq(:,end)));
addSlash = @(x) strrep(x,'_','\_');
sessions = cellfun(addSlash,sessionIDs(nComponents>0),'UniformOutput',0);
lines = cumsum(nComponents(nComponents>0))+0.5;
q = (cell2mat(qq(1:2,1:2)));
z = zBaseline(q,Unfind((length(qt)+1):(length(qt)*2),(length(qt)*2)),2);
%         z = nanzscore(q,[],2);
meanResponse = nanmean(z(:,[find(InIntervals(qt,[-0.1 0.1])) find(InIntervals(qt,[-0.1 0.1]))+length(qt)]),2);
flip = meanResponse<0; z(flip,:) = -z(flip,:);
z1 = z(:,1:length(qt)); z2 = z(:,(length(qt)+1):end);
smooth = [0 2];
z1 = zscore(z1,[],2); z2 = zscore(z2,[],2);
% figure('name',name);
clf
set(gcf,'position',[1 1 1920 1024])
matrices = {z1,z2,z2-z1};
titles = {[name ' in pre-task sleep'],[name ' in post-task sleep'],'(post-pre) difference'};
colors = {'k','r','m'};
for j=1:3
subplot(2,3,j); PlotColorMap(Smooth(matrices{j},smooth),~isnan(matrices{j}),'x',qt); PlotHVLines(0,'v','k--','linewidth',2);
PlotHVLines(lines,'h','w--','linewidth',2);
set(gca,'ytick',mean([[0;lines(1:end-1)] lines],2),'yticklabel',sessions);
title(titles{j});
subplot(2,3,3+j); semplot(qt,matrices{j},colors{j},smooth(end));
if j==1, semplot(qt,matrices{j+1},colors{j+1},smooth(end)); y=ylim; end
ylim(y); PlotHVLines(0,'v','k--','linewidth',2); PlotHVLines(0,'h','k--','linewidth',2); ylim(y);
ylabel('mean +/ s.e.m. (z-units)');
end
clims
drawnow
name = [name ' spikes'];
dd = cell2mat([cellfun(@zscore,ds(:,1),'uniformoutput',0) ds(:,2)]);
nComponents = Accumulate(cell2mat(qq(:,end)));
addSlash = @(x) strrep(x,'_','\_');
sessions = cellfun(addSlash,sessionIDs(nComponents>0),'UniformOutput',0);
lines = cumsum(nComponents(nComponents>0))+0.5;
q = (cell2mat(qq(1:2,1:2)));
z = zBaseline(q,Unfind((length(qt)+1):(length(qt)*2),(length(qt)*2)),2);
%         z = nanzscore(q,[],2);
meanResponse = nanmean(z(:,[find(InIntervals(qt,[-0.1 0.1])) find(InIntervals(qt,[-0.1 0.1]))+length(qt)]),2);
flip = meanResponse<0; z(flip,:) = -z(flip,:);
z1 = z(:,1:length(qt)); z2 = z(:,(length(qt)+1):end);
smooth = [0 2];
% figure('name',name);
clf
set(gcf,'position',[1 1 1920 1024])
matrices = {z1,z2,z2-z1};
titles = {[name ' in pre-task sleep'],[name ' in post-task sleep'],'(post-pre) difference'};
colors = {'k','r','m'};
for j=1:3
subplot(2,3,j); PlotColorMap(Smooth(matrices{j},smooth),~isnan(matrices{j}),'x',qt); PlotHVLines(0,'v','k--','linewidth',2);
PlotHVLines(lines,'h','w--','linewidth',2);
set(gca,'ytick',mean([[0;lines(1:end-1)] lines],2),'yticklabel',sessions);
title(titles{j});
subplot(2,3,3+j); semplot(qt,matrices{j},colors{j},smooth(end));
if j==1, semplot(qt,matrices{j+1},colors{j+1},smooth(end)); y=ylim; end
ylim(y); PlotHVLines(0,'v','k--','linewidth',2); PlotHVLines(0,'h','k--','linewidth',2); ylim(y);
ylabel('mean +/ s.e.m. (z-units)');
end
clims
drawnow
name = [name ' spikes'];
dd = cell2mat([cellfun(@zscore,ds(:,1),'uniformoutput',0) ds(:,2)]);
nComponents = Accumulate(cell2mat(qq(:,end)));
addSlash = @(x) strrep(x,'_','\_');
sessions = cellfun(addSlash,sessionIDs(nComponents>0),'UniformOutput',0);
lines = cumsum(nComponents(nComponents>0))+0.5;
q = (cell2mat(qq(151:end,1:2)));
z = zBaseline(q,Unfind((length(qt)+1):(length(qt)*2),(length(qt)*2)),2);
%         z = nanzscore(q,[],2);
meanResponse = nanmean(z(:,[find(InIntervals(qt,[-0.1 0.1])) find(InIntervals(qt,[-0.1 0.1]))+length(qt)]),2);
flip = meanResponse<0; z(flip,:) = -z(flip,:);
z1 = z(:,1:length(qt)); z2 = z(:,(length(qt)+1):end);
smooth = [0 2];
% figure('name',name);
clf
set(gcf,'position',[1 1 1920 1024])
matrices = {z1,z2,z2-z1};
titles = {[name ' in pre-task sleep'],[name ' in post-task sleep'],'(post-pre) difference'};
colors = {'k','r','m'};
for j=1:3
subplot(2,3,j); PlotColorMap(Smooth(matrices{j},smooth),~isnan(matrices{j}),'x',qt); PlotHVLines(0,'v','k--','linewidth',2);
PlotHVLines(lines,'h','w--','linewidth',2);
set(gca,'ytick',mean([[0;lines(1:end-1)] lines],2),'yticklabel',sessions);
title(titles{j});
subplot(2,3,3+j); semplot(qt,matrices{j},colors{j},smooth(end));
if j==1, semplot(qt,matrices{j+1},colors{j+1},smooth(end)); y=ylim; end
ylim(y); PlotHVLines(0,'v','k--','linewidth',2); PlotHVLines(0,'h','k--','linewidth',2); ylim(y);
ylabel('mean +/ s.e.m. (z-units)');
end
clims
drawnow
z1 = zscore(z1,[],2); z2 = zscore(z2,[],2);
% figure('name',name);
clf
set(gcf,'position',[1 1 1920 1024])
matrices = {z1,z2,z2-z1};
titles = {[name ' in pre-task sleep'],[name ' in post-task sleep'],'(post-pre) difference'};
colors = {'k','r','m'};
for j=1:3
subplot(2,3,j); PlotColorMap(Smooth(matrices{j},smooth),~isnan(matrices{j}),'x',qt); PlotHVLines(0,'v','k--','linewidth',2);
PlotHVLines(lines,'h','w--','linewidth',2);
set(gca,'ytick',mean([[0;lines(1:end-1)] lines],2),'yticklabel',sessions);
title(titles{j});
subplot(2,3,3+j); semplot(qt,matrices{j},colors{j},smooth(end));
if j==1, semplot(qt,matrices{j+1},colors{j+1},smooth(end)); y=ylim; end
ylim(y); PlotHVLines(0,'v','k--','linewidth',2); PlotHVLines(0,'h','k--','linewidth',2); ylim(y);
ylabel('mean +/ s.e.m. (z-units)');
end
clims
name = [name ' spikes'];
dd = cell2mat([cellfun(@zscore,ds(:,1),'uniformoutput',0) ds(:,2)]);
nComponents = Accumulate(cell2mat(qq(:,end)));
addSlash = @(x) strrep(x,'_','\_');
sessions = cellfun(addSlash,sessionIDs(nComponents>0),'UniformOutput',0);
lines = cumsum(nComponents(nComponents>0))+0.5;
q = (cell2mat(qq(151:end,1:2)));
%         z = zBaseline(q,Unfind((length(qt)+1):(length(qt)*2),(length(qt)*2)),2);
z = nanzscore(q,[],2);
meanResponse = nanmean(z(:,[find(InIntervals(qt,[-0.1 0.1])) find(InIntervals(qt,[-0.1 0.1]))+length(qt)]),2);
flip = meanResponse<0; z(flip,:) = -z(flip,:);
z1 = z(:,1:length(qt)); z2 = z(:,(length(qt)+1):end);
smooth = [0 2];
% figure('name',name);
clf
set(gcf,'position',[1 1 1920 1024])
matrices = {z1,z2,z2-z1};
titles = {[name ' in pre-task sleep'],[name ' in post-task sleep'],'(post-pre) difference'};
colors = {'k','r','m'};
for j=1:3
subplot(2,3,j); PlotColorMap(Smooth(matrices{j},smooth),~isnan(matrices{j}),'x',qt); PlotHVLines(0,'v','k--','linewidth',2);
PlotHVLines(lines,'h','w--','linewidth',2);
set(gca,'ytick',mean([[0;lines(1:end-1)] lines],2),'yticklabel',sessions);
title(titles{j});
subplot(2,3,3+j); semplot(qt,matrices{j},colors{j},smooth(end));
if j==1, semplot(qt,matrices{j+1},colors{j+1},smooth(end)); y=ylim; end
ylim(y); PlotHVLines(0,'v','k--','linewidth',2); PlotHVLines(0,'h','k--','linewidth',2); ylim(y);
ylabel('mean +/ s.e.m. (z-units)');
end
clims
drawnow
this = z2-z1;
clf
PlotColorMap(this);
PlotColorMap(this*this);
PlotColorMap(this*this');
PlotColorMap(this);
PlotColorMap(this'*this);
PlotColorMap(this'*this,'x',ht,'y',ht);
PlotColorMap(this'*this,'x',qt,'y',qt);
clikm
clim
open this
axis square
PlotColorMap(corr(this),'x',qt,'y',qt);
clim
PlotColorMap(corr(this));
PlotColorMap(corr(this),'x',qt,'y',qt);
clim
PlotColorMap(this);
[~,m] = max(this,[],2);
PlotColorMap(sortby(this,m));
PlotColorMap(sortby(zscore(this,[],2),m));
clim
clim([-1 1]*2);
PlotColorMap(sortby(this,m));
clim([-1 1]*2);
PlotColorMap(sortby(this,m));
PlotColorMap(sortby(this,sum(this,2)));
PlotColorMap(sortby(this,sum(this,2)),'x',qt);
clim
clim([-1 1]*2);
clim([-1 1]*5);
PlotColorMap(Smooth(sortby(this,sum(this,2)),[0 2])'x',qt);
PlotColorMap(Smooth(sortby(this,sum(this,2)),[0 2]),'x',qt);
clim
clim([-1 1]*5);
clim([-1 1]*4);
clim([-1 1]*3);
m = mean(this(:,qt<-0.5),2);
m = [mean(this(:,qt<-0.5),2) mean(this(:,InIntervals(qt,[0 0.1])),2)];
corr(m)
figure; PlotXY(m,'.');
hold all
plot(xlim,xlim,'k');
PlotHVLines(0,'h','k--','linewidth',2);
PlotHVLines(0,'k','k--','linewidth',2);
PlotHVLines(0,'v','k--','linewidth',2);
PlotColorMap(Smooth(sortby(this,sum(this,2)),[0 2]),'x',qt);
PlotColorMap(Smooth(sortby(this,sum(this,2)),[1 2]),'x',qt);
PlotColorMap(Smooth(sortby(this,sum(this,2)),[2 2]),'x',qt);
m = [mean(this(:,qt<-0.5),2) mean(this(:,InIntervals(qt,[0 0.1])),2)]; m(:,2) = m(:,2)-m(:,1);
figure; PlotXY(m,'.');
nancorr(m(:,1),m(:,2))
PlotHVLines(0,'v','k--','linewidth',2);
PlotHVLines(0,'k','k--','linewidth',2);
PlotHVLines(0,'h','k--','linewidth',2);
z1 = zscore(z1,[],2); z2 = zscore(z2,[],2);
% figure('name',name);
clf
set(gcf,'position',[1 1 1920 1024])
matrices = {z1,z2,z2-z1};
titles = {[name ' in pre-task sleep'],[name ' in post-task sleep'],'(post-pre) difference'};
colors = {'k','r','m'};
for j=1:3
subplot(2,3,j); PlotColorMap(Smooth(matrices{j},smooth),~isnan(matrices{j}),'x',qt); PlotHVLines(0,'v','k--','linewidth',2);
PlotHVLines(lines,'h','w--','linewidth',2);
set(gca,'ytick',mean([[0;lines(1:end-1)] lines],2),'yticklabel',sessions);
title(titles{j});
subplot(2,3,3+j); semplot(qt,matrices{j},colors{j},smooth(end));
if j==1, semplot(qt,matrices{j+1},colors{j+1},smooth(end)); y=ylim; end
ylim(y); PlotHVLines(0,'v','k--','linewidth',2); PlotHVLines(0,'h','k--','linewidth',2); ylim(y);
ylabel('mean +/ s.e.m. (z-units)');
end
clims
2
isempty(pfc)
max(pfc(:,2))
name = [name ' spikes'];
dd = cell2mat([cellfun(@zscore,ds(:,1),'uniformoutput',0) ds(:,2)]);
nComponents = Accumulate(cell2mat(qq(:,end)));
addSlash = @(x) strrep(x,'_','\_');
sessions = cellfun(addSlash,sessionIDs(nComponents>0),'UniformOutput',0);
lines = cumsum(nComponents(nComponents>0))+0.5;
q = (cell2mat(qq(:,1:2)));
%         z = zBaseline(q,Unfind((length(qt)+1):(length(qt)*2),(length(qt)*2)),2);
z = nanzscore(q,[],2);
meanResponse = nanmean(z(:,[find(InIntervals(qt,[-0.1 0.1])) find(InIntervals(qt,[-0.1 0.1]))+length(qt)]),2);
flip = meanResponse<0; z(flip,:) = -z(flip,:);
z1 = z(:,1:length(qt)); z2 = z(:,(length(qt)+1):end);
smooth = [0 2];
z1 = zscore(z1,[],2); z2 = zscore(z2,[],2);
% figure('name',name);
clf
set(gcf,'position',[1 1 1920 1024])
matrices = {z1,z2,z2-z1};
titles = {[name ' in pre-task sleep'],[name ' in post-task sleep'],'(post-pre) difference'};
colors = {'k','r','m'};
for j=1:3
subplot(2,3,j); PlotColorMap(Smooth(matrices{j},smooth),~isnan(matrices{j}),'x',qt); PlotHVLines(0,'v','k--','linewidth',2);
PlotHVLines(lines,'h','w--','linewidth',2);
set(gca,'ytick',mean([[0;lines(1:end-1)] lines],2),'yticklabel',sessions);
title(titles{j});
subplot(2,3,3+j); semplot(qt,matrices{j},colors{j},smooth(end));
if j==1, semplot(qt,matrices{j+1},colors{j+1},smooth(end)); y=ylim; end
ylim(y); PlotHVLines(0,'v','k--','linewidth',2); PlotHVLines(0,'h','k--','linewidth',2); ylim(y);
ylabel('mean +/ s.e.m. (z-units)');
end
clims
drawnow
345
raly
34
open X
Xo = get(batchO,'UserData'); % opto
Xs = get(batchS,'UserData'); % short
Xl = get(batchL,'UserData'); % long
XX = {Xs,Xl,Xo};
2
5
useClims(1,:) = clim;
useClims(2,:) = clim;
useClims(3,:) = clim;
useClims
for j=1:3
subplot(2,3,j); PlotColorMap(Smooth(matrices{j},smooth),~isnan(matrices{j}),'x',qt); PlotHVLines(0,'v','k--','linewidth',2);
PlotHVLines(lines,'h','w--','linewidth',2);
set(gca,'ytick',mean([[0;lines(1:end-1)] lines],2),'yticklabel',sessions);
usedClims(j,:) = clim;
title(titles{j});
subplot(2,3,3+j); semplot(qt,matrices{j},colors{j},smooth(end));
if j==1, semplot(qt,matrices{j+1},colors{j+1},smooth(end)); y=ylim; end
ylim(y); PlotHVLines(0,'v','k--','linewidth',2); PlotHVLines(0,'h','k--','linewidth',2); ylim(y);
ylabel('mean +/ s.e.m. (z-units)');
end
j
subplot(2,3,j); PlotColorMap(Smooth(matrices{j},smooth),~isnan(matrices{j}),'x',qt); PlotHVLines(0,'v','k--','linewidth',2);
clf
set(gcf,'position',[1 1 1920 1024])
matrices = {z1,z2,z2-z1};
titles = {[name ' in pre-task sleep'],[name ' in post-task sleep'],'(post-pre) difference'};
colors = {'k','r','m'};
for j=1:3
subplot(2,3,j); PlotColorMap(Smooth(matrices{j},smooth),~isnan(matrices{j}),'x',qt); PlotHVLines(0,'v','k--','linewidth',2);
PlotHVLines(lines,'h','w--','linewidth',2);
set(gca,'ytick',mean([[0;lines(1:end-1)] lines],2),'yticklabel',sessions);
usedClims(j,:) = clim;
title(titles{j});
subplot(2,3,3+j); semplot(qt,matrices{j},colors{j},smooth(end));
if j==1, semplot(qt,matrices{j+1},colors{j+1},smooth(end)); y=ylim; end
ylim(y); PlotHVLines(0,'v','k--','linewidth',2); PlotHVLines(0,'h','k--','linewidth',2); ylim(y);
ylabel('mean +/ s.e.m. (z-units)');
end
close all
46
56
nComponents = Accumulate(cell2mat(qq(:,end)));
addSlash = @(x) strrep(x,'_','\_');
sessions = cellfun(addSlash,sessionIDs(nComponents>0),'UniformOutput',0);
lines = cumsum(nComponents(nComponents>0))+0.5;
q = (cell2mat(qq(:,1:2)));
%         z = zBaseline(q,Unfind((length(qt)+1):(length(qt)*2),(length(qt)*2)),2);
z = nanzscore(q,[],2);
meanResponse = nanmean(z(:,[find(InIntervals(qt,[-0.1 0.1])) find(InIntervals(qt,[-0.1 0.1]))+length(qt)]),2);
flip = meanResponse<0; z(flip,:) = -z(flip,:);
z1 = z(:,1:length(qt)); z2 = z(:,(length(qt)+1):end);
smooth = [0 2];
z1 = zscore(z1,[],2); z2 = zscore(z2,[],2);
% figure('name',name);
clf
set(gcf,'position',[1 1 1920 1024])
matrices = {z1,z2,z2-z1};
titles = {[name ' in pre-task sleep'],[name ' in post-task sleep'],'(post-pre) difference'};
colors = {'k','r','m'};
for j=1:3
subplot(2,3,j); PlotColorMap(Smooth(matrices{j},smooth),~isnan(matrices{j}),'x',qt); PlotHVLines(0,'v','k--','linewidth',2);
PlotHVLines(lines,'h','w--','linewidth',2);
set(gca,'ytick',mean([[0;lines(1:end-1)] lines],2),'yticklabel',sessions);
usedClims(j,:) = clim;
title(titles{j});
subplot(2,3,3+j); semplot(qt,matrices{j},colors{j},smooth(end));
if j==1, semplot(qt,matrices{j+1},colors{j+1},smooth(end)); y=ylim; end
ylim(y); PlotHVLines(0,'v','k--','linewidth',2); PlotHVLines(0,'h','k--','linewidth',2); ylim(y);
ylabel('mean +/ s.e.m. (z-units)');
end
drawnow
usedClims
clims([-1 1]*min(abs(usedClims(:))))
clims([-1 1]*min(min(abs(usedClims(1:2,:)))))
clims([0 1]*min(min(abs(usedClims(1:2,:)))))
clims([max(usedClims(:,1)) min(usedClims(:,2))])
z = nanzscore(q,[],2);
meanResponse = nanmean(z(:,[find(InIntervals(qt,[-0.1 0.1])) find(InIntervals(qt,[-0.1 0.1]))+length(qt)]),2);
flip = meanResponse<0; z(flip,:) = -z(flip,:);
z1 = z(:,1:length(qt)); z2 = z(:,(length(qt)+1):end);
smooth = [0 2];
%         z1 = zscore(z1,[],2); z2 = zscore(z2,[],2);
% figure('name',name);
clf
set(gcf,'position',[1 1 1920 1024])
matrices = {z1,z2,z2-z1};
titles = {[name ' in pre-task sleep'],[name ' in post-task sleep'],'(post-pre) difference'};
colors = {'k','r','m'};
for j=1:3
subplot(2,3,j); PlotColorMap(Smooth(matrices{j},smooth),~isnan(matrices{j}),'x',qt); PlotHVLines(0,'v','k--','linewidth',2);
PlotHVLines(lines,'h','w--','linewidth',2);
set(gca,'ytick',mean([[0;lines(1:end-1)] lines],2),'yticklabel',sessions);
usedClims(j,:) = clim;
title(titles{j});
subplot(2,3,3+j); semplot(qt,matrices{j},colors{j},smooth(end));
if j==1, semplot(qt,matrices{j+1},colors{j+1},smooth(end)); y=ylim; end
ylim(y); PlotHVLines(0,'v','k--','linewidth',2); PlotHVLines(0,'h','k--','linewidth',2); ylim(y);
ylabel('mean +/ s.e.m. (z-units)');
end
clims([max(usedClims(:,1)) min(usedClims(:,2))])
drawnow
z = zBaseline(q,Unfind((length(qt)+1):(length(qt)*2),(length(qt)*2)),2);
%         z = nanzscore(q,[],2);
meanResponse = nanmean(z(:,[find(InIntervals(qt,[-0.1 0.1])) find(InIntervals(qt,[-0.1 0.1]))+length(qt)]),2);
flip = meanResponse<0; z(flip,:) = -z(flip,:);
z1 = z(:,1:length(qt)); z2 = z(:,(length(qt)+1):end);
smooth = [0 2];
%         z1 = zscore(z1,[],2); z2 = zscore(z2,[],2);
% figure('name',name);
clf
set(gcf,'position',[1 1 1920 1024])
matrices = {z1,z2,z2-z1};
titles = {[name ' in pre-task sleep'],[name ' in post-task sleep'],'(post-pre) difference'};
colors = {'k','r','m'};
for j=1:3
subplot(2,3,j); PlotColorMap(Smooth(matrices{j},smooth),~isnan(matrices{j}),'x',qt); PlotHVLines(0,'v','k--','linewidth',2);
PlotHVLines(lines,'h','w--','linewidth',2);
set(gca,'ytick',mean([[0;lines(1:end-1)] lines],2),'yticklabel',sessions);
usedClims(j,:) = clim;
title(titles{j});
subplot(2,3,3+j); semplot(qt,matrices{j},colors{j},smooth(end));
if j==1, semplot(qt,matrices{j+1},colors{j+1},smooth(end)); y=ylim; end
ylim(y); PlotHVLines(0,'v','k--','linewidth',2); PlotHVLines(0,'h','k--','linewidth',2); ylim(y);
ylabel('mean +/ s.e.m. (z-units)');
end
clims([max(usedClims(:,1)) min(usedClims(:,2))])
drawnow
q = (cell2mat(qq(17:end,1:2)));
z = zBaseline(q,Unfind((length(qt)+1):(length(qt)*2),(length(qt)*2)),2);
%         z = nanzscore(q,[],2);
meanResponse = nanmean(z(:,[find(InIntervals(qt,[-0.1 0.1])) find(InIntervals(qt,[-0.1 0.1]))+length(qt)]),2);
flip = meanResponse<0; z(flip,:) = -z(flip,:);
z1 = z(:,1:length(qt)); z2 = z(:,(length(qt)+1):end);
smooth = [0 2];
%         z1 = zscore(z1,[],2); z2 = zscore(z2,[],2);
% figure('name',name);
clf
set(gcf,'position',[1 1 1920 1024])
matrices = {z1,z2,z2-z1};
titles = {[name ' in pre-task sleep'],[name ' in post-task sleep'],'(post-pre) difference'};
colors = {'k','r','m'};
for j=1:3
subplot(2,3,j); PlotColorMap(Smooth(matrices{j},smooth),~isnan(matrices{j}),'x',qt); PlotHVLines(0,'v','k--','linewidth',2);
PlotHVLines(lines,'h','w--','linewidth',2);
set(gca,'ytick',mean([[0;lines(1:end-1)] lines],2),'yticklabel',sessions);
usedClims(j,:) = clim;
title(titles{j});
subplot(2,3,3+j); semplot(qt,matrices{j},colors{j},smooth(end));
if j==1, semplot(qt,matrices{j+1},colors{j+1},smooth(end)); y=ylim; end
ylim(y); PlotHVLines(0,'v','k--','linewidth',2); PlotHVLines(0,'h','k--','linewidth',2); ylim(y);
ylabel('mean +/ s.e.m. (z-units)');
end
clims([max(usedClims(:,1)) min(usedClims(:,2))])
drawnow
z1 = zscore(z1,[],2); z2 = zscore(z2,[],2);
% figure('name',name);
clf
set(gcf,'position',[1 1 1920 1024])
matrices = {z1,z2,z2-z1};
titles = {[name ' in pre-task sleep'],[name ' in post-task sleep'],'(post-pre) difference'};
colors = {'k','r','m'};
for j=1:3
subplot(2,3,j); PlotColorMap(Smooth(matrices{j},smooth),~isnan(matrices{j}),'x',qt); PlotHVLines(0,'v','k--','linewidth',2);
PlotHVLines(lines,'h','w--','linewidth',2);
set(gca,'ytick',mean([[0;lines(1:end-1)] lines],2),'yticklabel',sessions);
usedClims(j,:) = clim;
title(titles{j});
subplot(2,3,3+j); semplot(qt,matrices{j},colors{j},smooth(end));
if j==1, semplot(qt,matrices{j+1},colors{j+1},smooth(end)); y=ylim; end
ylim(y); PlotHVLines(0,'v','k--','linewidth',2); PlotHVLines(0,'h','k--','linewidth',2); ylim(y);
ylabel('mean +/ s.e.m. (z-units)');
end
clims([max(usedClims(:,1)) min(usedClims(:,2))])
drawnow
z = zBaseline(q,Unfind((length(qt)+1):(length(qt)*2),(length(qt)*2)),2);
%         z = nanzscore(q,[],2);
meanResponse = nanmean(z(:,[find(InIntervals(qt,[-0.1 0.1])) find(InIntervals(qt,[-0.1 0.1]))+length(qt)]),2);
flip = meanResponse<0; z(flip,:) = -z(flip,:);
z1 = z(:,1:length(qt)); z2 = z(:,(length(qt)+1):end);
smooth = [0 2];
z1 = zscore(z1,[],2); z2 = zscore(z2,[],2);
% figure('name',name);
clf
set(gcf,'position',[1 1 1920 1024])
matrices = {z1,z2,z2-z1};
titles = {[name ' in pre-task sleep'],[name ' in post-task sleep'],'(post-pre) difference'};
colors = {'k','r','m'};
for j=1:3
subplot(2,3,j); PlotColorMap(Smooth(matrices{j},smooth),~isnan(matrices{j}),'x',qt); PlotHVLines(0,'v','k--','linewidth',2);
PlotHVLines(lines,'h','w--','linewidth',2);
set(gca,'ytick',mean([[0;lines(1:end-1)] lines],2),'yticklabel',sessions);
usedClims(j,:) = clim;
title(titles{j});
subplot(2,3,3+j); semplot(qt,matrices{j},colors{j},smooth(end));
if j==1, semplot(qt,matrices{j+1},colors{j+1},smooth(end)); y=ylim; end
ylim(y); PlotHVLines(0,'v','k--','linewidth',2); PlotHVLines(0,'h','k--','linewidth',2); ylim(y);
ylabel('mean +/ s.e.m. (z-units)');
end
clims([max(usedClims(:,1)) min(usedClims(:,2))])
drawnow
z = zBaseline(q,Unfind((length(qt)+1):(length(qt)*2),(length(qt)*2)),2);
%         z = nanzscore(q,[],2);
meanResponse = nanmean(z(:,[find(InIntervals(qt,[-0.1 0.1])) find(InIntervals(qt,[-0.1 0.1]))+length(qt)]),2);
flip = meanResponse<0; z(flip,:) = -z(flip,:);
z1 = z(:,1:length(qt)); z2 = z(:,(length(qt)+1):end);
smooth = [0 2];
%         z1 = zscore(z1,[],2); z2 = zscore(z2,[],2);
% figure('name',name);
clf
set(gcf,'position',[1 1 1920 1024])
matrices = {z1,z2,z2-z1};
titles = {[name ' in pre-task sleep'],[name ' in post-task sleep'],'(post-pre) difference'};
colors = {'k','r','m'};
for j=1:3
subplot(2,3,j); PlotColorMap(Smooth(matrices{j},smooth),~isnan(matrices{j}),'x',qt); PlotHVLines(0,'v','k--','linewidth',2);
PlotHVLines(lines,'h','w--','linewidth',2);
set(gca,'ytick',mean([[0;lines(1:end-1)] lines],2),'yticklabel',sessions);
usedClims(j,:) = clim;
title(titles{j});
subplot(2,3,3+j); semplot(qt,matrices{j},colors{j},smooth(end));
if j==1, semplot(qt,matrices{j+1},colors{j+1},smooth(end)); y=ylim; end
ylim(y); PlotHVLines(0,'v','k--','linewidth',2); PlotHVLines(0,'h','k--','linewidth',2); ylim(y);
ylabel('mean +/ s.e.m. (z-units)');
end
clims([max(usedClims(:,1)) min(usedClims(:,2))])
drawnow
z1 = zscore(z1,[],2); z2 = zscore(z2,[],2);
% figure('name',name);
clf
set(gcf,'position',[1 1 1920 1024])
matrices = {z1,z2,z2-z1};
titles = {[name ' in pre-task sleep'],[name ' in post-task sleep'],'(post-pre) difference'};
colors = {'k','r','m'};
for j=1:3
subplot(2,3,j); PlotColorMap(Smooth(matrices{j},smooth),~isnan(matrices{j}),'x',qt); PlotHVLines(0,'v','k--','linewidth',2);
PlotHVLines(lines,'h','w--','linewidth',2);
set(gca,'ytick',mean([[0;lines(1:end-1)] lines],2),'yticklabel',sessions);
usedClims(j,:) = clim;
title(titles{j});
subplot(2,3,3+j); semplot(qt,matrices{j},colors{j},smooth(end));
if j==1, semplot(qt,matrices{j+1},colors{j+1},smooth(end)); y=ylim; end
ylim(y); PlotHVLines(0,'v','k--','linewidth',2); PlotHVLines(0,'h','k--','linewidth',2); ylim(y);
ylabel('mean +/ s.e.m. (z-units)');
end
clims([max(usedClims(:,1)) min(usedClims(:,2))])
drawnow
clim
clims([-3 3]);
nComponents = Accumulate(cell2mat(qq(:,end)));
addSlash = @(x) strrep(x,'_','\_');
sessions = cellfun(addSlash,sessionIDs(nComponents>0),'UniformOutput',0);
lines = cumsum(nComponents(nComponents>0))+0.5;
q = (cell2mat(qq(:,1:2)));
meanResponse = nanmean(q(:,[find(InIntervals(qt,[-0.1 0.1])) find(InIntervals(qt,[-0.1 0.1]))+length(qt)]),2);
%                 z = zBaseline(q,Unfind((length(qt)+1):(length(qt)*2),(length(qt)*2)),2);
%         z = nanzscore(q,[],2);
z=q;
flip = meanResponse<0; z(flip,:) = -z(flip,:);
z1 = z(:,1:length(qt)); z2 = z(:,(length(qt)+1):end);
smooth = [0 2];
%         z1 = zscore(z1,[],2); z2 = zscore(z2,[],2);
% figure('name',name);
clf
set(gcf,'position',[1 1 1920 1024])
matrices = {z1,z2,z2-z1};
titles = {[name ' in pre-task sleep'],[name ' in post-task sleep'],'(post-pre) difference'};
colors = {'k','r','m'};
for j=1:3
subplot(2,3,j); PlotColorMap(Smooth(matrices{j},smooth),~isnan(matrices{j}),'x',qt); PlotHVLines(0,'v','k--','linewidth',2);
PlotHVLines(lines,'h','w--','linewidth',2);
set(gca,'ytick',mean([[0;lines(1:end-1)] lines],2),'yticklabel',sessions);
usedClims(j,:) = clim;
title(titles{j});
subplot(2,3,3+j); semplot(qt,matrices{j},colors{j},smooth(end));
if j==1, semplot(qt,matrices{j+1},colors{j+1},smooth(end)); y=ylim; end
ylim(y); PlotHVLines(0,'v','k--','linewidth',2); PlotHVLines(0,'h','k--','linewidth',2); ylim(y);
ylabel('mean +/ s.e.m. (z-units)');
end
clims([max(usedClims(:,1)) min(usedClims(:,2))])
drawnow
clf
set(gcf,'position',[1 1 1920 1024])
matrices = {z1,z2,z2-z1};
titles = {[name ' in pre-task sleep'],[name ' in post-task sleep'],'(post-pre) difference'};
colors = {'k','r','m'};
for j=1:3
subplot(2,3,j); PlotColorMap(Smooth(matrices{j},smooth),~isnan(matrices{j}),'x',qt); PlotHVLines(0,'v','k--','linewidth',2);
PlotHVLines(lines,'h','w--','linewidth',2);
set(gca,'ytick',mean([[0;lines(1:end-1)] lines],2),'yticklabel',sessions);
usedClims(j,:) = clim;
title(titles{j});
subplot(2,3,3+j); semplot(qt,matrices{j},colors{j},smooth(end));
if j==1, semplot(qt,matrices{j+1},colors{j+1},smooth(end)); y=ylim; end
PlotHVLines(0,'v','k--','linewidth',2); PlotHVLines(0,'h','k--','linewidth',2); ylim(y);
ylabel('mean +/ s.e.m. (z-units)');
end
clims([max(usedClims(:,1)) min(usedClims(:,2))])
clf
set(gcf,'position',[1 1 1920 1024])
matrices = {z1,z2,z2-z1};
titles = {[name ' in pre-task sleep'],[name ' in post-task sleep'],'(post-pre) difference'};
colors = {'k','r','m'};
for j=1:3
subplot(2,3,j); PlotColorMap(Smooth(matrices{j},smooth),~isnan(matrices{j}),'x',qt); PlotHVLines(0,'v','k--','linewidth',2);
PlotHVLines(lines,'h','w--','linewidth',2);
set(gca,'ytick',mean([[0;lines(1:end-1)] lines],2),'yticklabel',sessions);
usedClims(j,:) = clim;
title(titles{j});
subplot(2,3,3+j); semplot(qt,matrices{j},colors{j},smooth(end));
if j==1, semplot(qt,matrices{j+1},colors{j+1},smooth(end)); y=ylim; end
%             ylim(y); PlotHVLines(0,'v','k--','linewidth',2); PlotHVLines(0,'h','k--','linewidth',2); ylim(y);
ylabel('mean +/ s.e.m. (z-units)');
end
clims([max(usedClims(:,1)) min(usedClims(:,2))])
z = zBaseline(q,Unfind((length(qt)+1):(length(qt)*2),(length(qt)*2)),2);
z1 = z(:,1:length(qt)); z2 = z(:,(length(qt)+1):end);
smooth = [0 2];
%         z1 = zscore(z1,[],2); z2 = zscore(z2,[],2);
% figure('name',name);
clf
set(gcf,'position',[1 1 1920 1024])
matrices = {z1,z2,z2-z1};
titles = {[name ' in pre-task sleep'],[name ' in post-task sleep'],'(post-pre) difference'};
colors = {'k','r','m'};
for j=1:3
subplot(2,3,j); PlotColorMap(Smooth(matrices{j},smooth),~isnan(matrices{j}),'x',qt); PlotHVLines(0,'v','k--','linewidth',2);
PlotHVLines(lines,'h','w--','linewidth',2);
set(gca,'ytick',mean([[0;lines(1:end-1)] lines],2),'yticklabel',sessions);
usedClims(j,:) = clim;
title(titles{j});
subplot(2,3,3+j); semplot(qt,matrices{j},colors{j},smooth(end));
if j==1, semplot(qt,matrices{j+1},colors{j+1},smooth(end)); y=ylim; end
%             ylim(y); PlotHVLines(0,'v','k--','linewidth',2); PlotHVLines(0,'h','k--','linewidth',2); ylim(y);
ylabel('mean +/ s.e.m. (z-units)');
end
clims([max(usedClims(:,1)) min(usedClims(:,2))])
drawnow
z = nanzscore(q,[],2);
z1 = z(:,1:length(qt)); z2 = z(:,(length(qt)+1):end);
smooth = [0 2];
%         z1 = zscore(z1,[],2); z2 = zscore(z2,[],2);
% figure('name',name);
clf
set(gcf,'position',[1 1 1920 1024])
matrices = {z1,z2,z2-z1};
titles = {[name ' in pre-task sleep'],[name ' in post-task sleep'],'(post-pre) difference'};
colors = {'k','r','m'};
for j=1:3
subplot(2,3,j); PlotColorMap(Smooth(matrices{j},smooth),~isnan(matrices{j}),'x',qt); PlotHVLines(0,'v','k--','linewidth',2);
PlotHVLines(lines,'h','w--','linewidth',2);
set(gca,'ytick',mean([[0;lines(1:end-1)] lines],2),'yticklabel',sessions);
usedClims(j,:) = clim;
title(titles{j});
subplot(2,3,3+j); semplot(qt,matrices{j},colors{j},smooth(end));
if j==1, semplot(qt,matrices{j+1},colors{j+1},smooth(end)); y=ylim; end
%             ylim(y); PlotHVLines(0,'v','k--','linewidth',2); PlotHVLines(0,'h','k--','linewidth',2); ylim(y);
ylabel('mean +/ s.e.m. (z-units)');
end
clims([max(usedClims(:,1)) min(usedClims(:,2))])
nComponents = Accumulate(cell2mat(qq(:,end)));
addSlash = @(x) strrep(x,'_','\_');
sessions = cellfun(addSlash,sessionIDs(nComponents>0),'UniformOutput',0);
lines = cumsum(nComponents(nComponents>0))+0.5;
q = (cell2mat(qq(:,1:2)));
meanResponse = nanmean(q(:,[find(InIntervals(qt,[-0.1 0.1])) find(InIntervals(qt,[-0.1 0.1]))+length(qt)]),2);
%                 z = zBaseline(q,Unfind((length(qt)+1):(length(qt)*2),(length(qt)*2)),2);
%         z = nanzscore(q,[],2);
%         z=q;
%         flip = meanResponse<0; z(flip,:) = -z(flip,:);
z1 = z(:,1:length(qt)); z2 = z(:,(length(qt)+1):end);
smooth = [0 2];
%         z1 = zscore(z1,[],2); z2 = zscore(z2,[],2);
% figure('name',name);
clf
set(gcf,'position',[1 1 1920 1024])
matrices = {z1,z2,z2-z1};
titles = {[name ' in pre-task sleep'],[name ' in post-task sleep'],'(post-pre) difference'};
colors = {'k','r','m'};
for j=1:3
subplot(2,3,j); PlotColorMap(Smooth(matrices{j},smooth),~isnan(matrices{j}),'x',qt); PlotHVLines(0,'v','k--','linewidth',2);
PlotHVLines(lines,'h','w--','linewidth',2);
set(gca,'ytick',mean([[0;lines(1:end-1)] lines],2),'yticklabel',sessions);
usedClims(j,:) = clim;
title(titles{j});
subplot(2,3,3+j); semplot(qt,matrices{j},colors{j},smooth(end));
if j==1, semplot(qt,matrices{j+1},colors{j+1},smooth(end)); y=ylim; end
%             ylim(y); PlotHVLines(0,'v','k--','linewidth',2); PlotHVLines(0,'h','k--','linewidth',2); ylim(y);
ylabel('mean +/ s.e.m. (z-units)');
end
clims([max(usedClims(:,1)) min(usedClims(:,2))])
drawnow
nComponents = Accumulate(cell2mat(qq(:,end)));
addSlash = @(x) strrep(x,'_','\_');
sessions = cellfun(addSlash,sessionIDs(nComponents>0),'UniformOutput',0);
lines = cumsum(nComponents(nComponents>0))+0.5;
q = (cell2mat(qq(:,1:2)));
meanResponse = nanmean(q(:,[find(InIntervals(qt,[-0.1 0.1])) find(InIntervals(qt,[-0.1 0.1]))+length(qt)]),2);
%                 z = zBaseline(q,Unfind((length(qt)+1):(length(qt)*2),(length(qt)*2)),2);
z = nanzscore(q,[],2);
%         flip = meanResponse<0; z(flip,:) = -z(flip,:);
z1 = z(:,1:length(qt)); z2 = z(:,(length(qt)+1):end);
smooth = [0 2];
%         z1 = zscore(z1,[],2); z2 = zscore(z2,[],2);
% figure('name',name);
clf
set(gcf,'position',[1 1 1920 1024])
matrices = {z1,z2,z2-z1};
titles = {[name ' in pre-task sleep'],[name ' in post-task sleep'],'(post-pre) difference'};
colors = {'k','r','m'};
for j=1:3
subplot(2,3,j); PlotColorMap(Smooth(matrices{j},smooth),~isnan(matrices{j}),'x',qt); PlotHVLines(0,'v','k--','linewidth',2);
PlotHVLines(lines,'h','w--','linewidth',2);
set(gca,'ytick',mean([[0;lines(1:end-1)] lines],2),'yticklabel',sessions);
usedClims(j,:) = clim;
title(titles{j});
subplot(2,3,3+j); semplot(qt,matrices{j},colors{j},smooth(end));
if j==1, semplot(qt,matrices{j+1},colors{j+1},smooth(end)); y=ylim; end
%             ylim(y); PlotHVLines(0,'v','k--','linewidth',2); PlotHVLines(0,'h','k--','linewidth',2); ylim(y);
ylabel('mean +/ s.e.m. (z-units)');
end
clims([max(usedClims(:,1)) min(usedClims(:,2))])
drawnow
nComponents = Accumulate(cell2mat(qq(:,end)));
addSlash = @(x) strrep(x,'_','\_');
sessions = cellfun(addSlash,sessionIDs(nComponents>0),'UniformOutput',0);
lines = cumsum(nComponents(nComponents>0))+0.5;
q = (cell2mat(qq(:,1:2)));
meanResponse = nanmean(q(:,[find(InIntervals(qt,[-0.1 0.1])) find(InIntervals(qt,[-0.1 0.1]))+length(qt)]),2);
%                 z = zBaseline(q,Unfind((length(qt)+1):(length(qt)*2),(length(qt)*2)),2);
z = nanzscore(q,[],2);
%         flip = meanResponse<0; z(flip,:) = -z(flip,:);
z1 = z(:,1:length(qt)); z2 = z(:,(length(qt)+1):end);
smooth = [0 2];
z1 = zscore(z1,[],2); z2 = zscore(z2,[],2);
% figure('name',name);
clf
set(gcf,'position',[1 1 1920 1024])
matrices = {z1,z2,z2-z1};
titles = {[name ' in pre-task sleep'],[name ' in post-task sleep'],'(post-pre) difference'};
colors = {'k','r','m'};
for j=1:3
subplot(2,3,j); PlotColorMap(Smooth(matrices{j},smooth),~isnan(matrices{j}),'x',qt); PlotHVLines(0,'v','k--','linewidth',2);
PlotHVLines(lines,'h','w--','linewidth',2);
set(gca,'ytick',mean([[0;lines(1:end-1)] lines],2),'yticklabel',sessions);
usedClims(j,:) = clim;
title(titles{j});
subplot(2,3,3+j); semplot(qt,matrices{j},colors{j},smooth(end));
if j==1, semplot(qt,matrices{j+1},colors{j+1},smooth(end)); y=ylim; end
%             ylim(y); PlotHVLines(0,'v','k--','linewidth',2); PlotHVLines(0,'h','k--','linewidth',2); ylim(y);
ylabel('mean +/ s.e.m. (z-units)');
end
clims([max(usedClims(:,1)) min(usedClims(:,2))])
drawnow
nComponents = Accumulate(cell2mat(qq(:,end)));
addSlash = @(x) strrep(x,'_','\_');
sessions = cellfun(addSlash,sessionIDs(nComponents>0),'UniformOutput',0);
lines = cumsum(nComponents(nComponents>0))+0.5;
q = (cell2mat(qq(:,1:2)));
meanResponse = nanmean(q(:,[find(InIntervals(qt,[-0.1 0.1])) find(InIntervals(qt,[-0.1 0.1]))+length(qt)]),2);
%                 z = zBaseline(q,Unfind((length(qt)+1):(length(qt)*2),(length(qt)*2)),2);
z = nanzscore(q,[],2);
flip = meanResponse<0; z(flip,:) = -z(flip,:);
z1 = z(:,1:length(qt)); z2 = z(:,(length(qt)+1):end);
smooth = [0 2];
z1 = zscore(z1,[],2); z2 = zscore(z2,[],2);
% figure('name',name);
clf
set(gcf,'position',[1 1 1920 1024])
matrices = {z1,z2,z2-z1};
titles = {[name ' in pre-task sleep'],[name ' in post-task sleep'],'(post-pre) difference'};
colors = {'k','r','m'};
for j=1:3
subplot(2,3,j); PlotColorMap(Smooth(matrices{j},smooth),~isnan(matrices{j}),'x',qt); PlotHVLines(0,'v','k--','linewidth',2);
PlotHVLines(lines,'h','w--','linewidth',2);
set(gca,'ytick',mean([[0;lines(1:end-1)] lines],2),'yticklabel',sessions);
usedClims(j,:) = clim;
title(titles{j});
subplot(2,3,3+j); semplot(qt,matrices{j},colors{j},smooth(end));
if j==1, semplot(qt,matrices{j+1},colors{j+1},smooth(end)); y=ylim; end
%             ylim(y); PlotHVLines(0,'v','k--','linewidth',2); PlotHVLines(0,'h','k--','linewidth',2); ylim(y);
ylabel('mean +/ s.e.m. (z-units)');
end
clims([max(usedClims(:,1)) min(usedClims(:,2))])
drawnow
nComponents = Accumulate(cell2mat(qq(:,end)));
addSlash = @(x) strrep(x,'_','\_');
sessions = cellfun(addSlash,sessionIDs(nComponents>0),'UniformOutput',0);
lines = cumsum(nComponents(nComponents>0))+0.5;
q = (cell2mat(qq(:,1:2)));
meanResponse = nanmean(q(:,[find(InIntervals(qt,[-0.1 0.1])) find(InIntervals(qt,[-0.1 0.1]))+length(qt)]),2);
%                 z = zBaseline(q,Unfind((length(qt)+1):(length(qt)*2),(length(qt)*2)),2);
z = nanzscore(q,[],2);
flip = meanResponse<0; z(flip,:) = -z(flip,:);
z1 = z(:,1:length(qt)); z2 = z(:,(length(qt)+1):end);
smooth = [0 2];
%         z1 = zscore(z1,[],2); z2 = zscore(z2,[],2);
% figure('name',name);
clf
set(gcf,'position',[1 1 1920 1024])
matrices = {z1,z2,z2-z1};
titles = {[name ' in pre-task sleep'],[name ' in post-task sleep'],'(post-pre) difference'};
colors = {'k','r','m'};
for j=1:3
subplot(2,3,j); PlotColorMap(Smooth(matrices{j},smooth),~isnan(matrices{j}),'x',qt); PlotHVLines(0,'v','k--','linewidth',2);
PlotHVLines(lines,'h','w--','linewidth',2);
set(gca,'ytick',mean([[0;lines(1:end-1)] lines],2),'yticklabel',sessions);
usedClims(j,:) = clim;
title(titles{j});
subplot(2,3,3+j); semplot(qt,matrices{j},colors{j},smooth(end));
if j==1, semplot(qt,matrices{j+1},colors{j+1},smooth(end)); y=ylim; end
%             ylim(y); PlotHVLines(0,'v','k--','linewidth',2); PlotHVLines(0,'h','k--','linewidth',2); ylim(y);
ylabel('mean +/ s.e.m. (z-units)');
end
clims([max(usedClims(:,1)) min(usedClims(:,2))])
drawnow
sum(flip)
find(flip)
figure; PlotColorMap(q);
plot(q(2,:))
plot(Smooth(q(2,:),1))
plot(Smooth(q(7,:),1))
i=1
X = XX{condition};
qq = {}; sessionIDs = {}; kk = 0; nSessions = size(X,1);
[pfc,hpc,deltas,ripples,sleep,session,MergePoints,clickers,sws,basepath] = X{i,:};
[parentFolder,basename] = fileparts(basepath);
[~,projectName] = fileparts(parentFolder);
sessionID = [projectName '_' basename];
sessionIDs{i,1} = sessionID;
preSleep = SubtractIntervals(sleep(1,:), SubtractIntervals([0 Inf],sws));
postSleep = SubtractIntervals(sleep(2:end,:), SubtractIntervals([0 Inf],sws));
pre = InIntervals(ripples,preSleep);
post = InIntervals(ripples,postSleep);
%         % Restrict to 1h of sws only:
%         limit = Unshift(3600,preSleep); preSleep(preSleep(:,1)>limit,:) = []; if preSleep(end,2)>limit, preSleep(end,2) = limit; end
%         limit = Unshift(3600,postSleep); postSleep(postSleep(:,1)>limit,:) = []; if postSleep(end,2)>limit, postSleep(end,2) = limit; end
training = ConsolidateIntervals(MergePoints.timestamps(cellfun(@(x) ~isempty(strfind(x,'train')),MergePoints.foldernames),:)); % training epochs are named with "train" somewhere in the folder name
exploration = sortrows([clickers{2};clickers{3}]); exploration(any(isnan(exploration),2),:) = []; exploration = ConsolidateIntervals(exploration);
exploration = Restrict(exploration,training);
if structure==1, spikes = hpc; else, spikes = pfc; end
kind
if kind==2, intervals = training; else, intervals = exploration; end
rng(0);
[templates,correlations,weights] = ActivityTemplatesICA(Restrict(spikes,intervals,'shift','on'),'binsize',0.1);
re = ReactivationStrength(spikes,templates,'step',0.01,'binSize',0.1);
for j=1:size(re,2)-1
[q,qt] = PETH(re(:,[1 1+j]),ripples(:,1),'durations',[-1 1]*2,'nBins',501);
kk = kk+1;
qq{kk,1} = nanmean(q(pre,:),1);
qq{kk,2} = nanmean(q(post,:),1);
qq{kk,3} = i;
%                             qqs{kk,1} = q;
end
q = (cell2mat(qq(:,1:2)));
meanResponse = nanmean(q(:,[find(InIntervals(qt,[-0.1 0.1])) find(InIntervals(qt,[-0.1 0.1]))+length(qt)]),2);
flip = meanResponse<0; q(flip,:) = -q(flip,:);
%                 z = zBaseline(q,Unfind((length(qt)+1):(length(qt)*2),(length(qt)*2)),2);
z = nanzscore(q,[],2);
z1 = z(:,1:length(qt)); z2 = z(:,(length(qt)+1):end);
smooth = [0 2];
z1 = zscore(z1,[],2); z2 = zscore(z2,[],2);
% figure('name',name);
clf
set(gcf,'position',[1 1 1920 1024])
matrices = {z1,z2,z2-z1};
titles = {[name ' in pre-task sleep'],[name ' in post-task sleep'],'(post-pre) difference'};
colors = {'k','r','m'};
for j=1:3
subplot(2,3,j); PlotColorMap(Smooth(matrices{j},smooth),~isnan(matrices{j}),'x',qt); PlotHVLines(0,'v','k--','linewidth',2);
PlotHVLines(lines,'h','w--','linewidth',2);
set(gca,'ytick',mean([[0;lines(1:end-1)] lines],2),'yticklabel',sessions);
usedClims(j,:) = clim;
title(titles{j});
subplot(2,3,3+j); semplot(qt,matrices{j},colors{j},smooth(end));
if j==1, semplot(qt,matrices{j+1},colors{j+1},smooth(end)); y=ylim; end
ylim(y); PlotHVLines(0,'v','k--','linewidth',2); PlotHVLines(0,'h','k--','linewidth',2); ylim(y);
ylabel('mean +/ s.e.m. (z-units)');
end
clims([max(usedClims(:,1)) min(usedClims(:,2))])
max(spikes(:,2))
bar(weights(:,2))
for j=1:max(spikes(:,2))
[q,qt] = PETH(spikes(spikes(:,2)==j),ripples(:,1),'durations',[-1 1]*2,'nBins',501);
this(j,:) = nanmean(q);
this0(j,:) = nanmean(q(pre,:));
this1(j,:) = nanmean(q(post,:));
end
open this
PlotColorMap(zscore(this,[],2));
PlotColorMap(sortby(zscore(this,[],2),weights(:,2)))
subplot(1,2,1) ;PlotColorMap(sortby(zscore(this,[],2),weights(:,2)))
subplot(1,2,2) ; barh(weights(:,2));
subplot(1,2,1) ;PlotColorMap(sortby(zscore(this,[],2),weights(:,2)))
subplot(1,2,2) ; barh(sortby(weights(:,2),weights(:,2)));
subplot(1,2,1) ;PlotColorMap(Smooth(sortby(zscore(this,[],2),weights(:,2)),2));
subplot(1,2,2) ; barh(sortby(weights(:,2),weights(:,2)));
subplot(1,2,1) ;PlotColorMap(Smooth(sortby(zscore(this,[],2),weights(:,2)),5));
subplot(1,2,1) ;PlotColorMap(Smooth(sortby(zscore(this,[],2),weights(:,2)),[2 5]));
subplot(1,2,1) ;PlotColorMap(Shrink(sortby(zscore(this,[],2),weights(:,2)),10,1))
subplot(1,2,1) ;PlotColorMap(Shrink(sortby(zscore(this,[],2),weights(:,2)),15,1))
subplot(1,2,1) ;PlotColorMap(Shrink(sortby(zscore(this,[],2),weights(:,2)),15,5),'x',qt)
subplot(1,2,1) ;PlotColorMap(Shrink(sortby(zscore(this,[],2),weights(:,2)),20,5),'x',qt)
subplot(1,2,1) ;PlotColorMap(Shrink(sortby(zscore(this,[],2),weights(:,2)),20,5),'x',qt); xlim([-1 1]*0.5);
figure; PETH(re(:,[1 3]),ripples(:,2))
PETH(re(:,[1 3]),ripples(:,2),'nbins',501)
PETH(re(:,[1 3]),ripples(:,2),'nbins',501,'durations',[-1 1]*2)
PETH(re(:,[1 3]),ripples(:,2),'nbins',501,'durations',[-1 1]*5)
PETH(re(:,[1 3]),ripples(:,2),'nbins',1001,'durations',[-1 1]*5)
re0 = ReactivationStrength(spikes,templates,'step',0.01,'binSize',0.1);
re1 = ReactivationStrength(spikes,-templates,'step',0.01,'binSize',0.1);
PETH(re0(:,[1 3]),ripples(:,2),'nbins',1001,'durations',[-1 1]*5)
hold all
PETH(re1(:,[1 3]),ripples(:,2),'nbins',1001,'durations',[-1 1]*5)
open templates
tt = weights(:,2); tt = tt*tt';
open tt
PlotColorMap(tt);
tt = weights(:,2); tt = tt*tt'; tt(eye(size(tt))==1) = 0;
re0 = ReactivationStrength(spikes,tt,'step',0.01,'binSize',0.1);
tt = -weights(:,2); tt = tt*tt'; tt(eye(size(tt))==1) = 0;
re1 = ReactivationStrength(spikes,tt,'step',0.01,'binSize',0.1);
PETH(re0,ripples(:,2),'nbins',1001,'durations',[-1 1]*5)
clf
PETH(re0,ripples(:,2),'nbins',1001,'durations',[-1 1]*5)
hold all
PETH(re1,ripples(:,2),'nbins',1001,'durations',[-1 1]*5)
[h,stats] = CompareDistributions(matrices{1},matrices{2})
for j=1:3
subplot(2,3,j); PlotColorMap(Smooth(matrices{j},smooth),~isnan(matrices{j}),'x',qt); PlotHVLines(0,'v','k--','linewidth',2);
PlotHVLines(lines,'h','w--','linewidth',2);
set(gca,'ytick',mean([[0;lines(1:end-1)] lines],2),'yticklabel',sessions);
usedClims(j,:) = clim;
title(titles{j});
subplot(2,3,3+j); semplot(qt,matrices{j},colors{j},smooth(end));
if j==1, semplot(qt,matrices{j+1},colors{j+1},smooth(end)); y=ylim; end
if j==2, [h,stats] = CompareDistributions(matrices{1},matrices{2}); end
ylim(y); PlotHVLines(0,'v','k--','linewidth',2); PlotHVLines(0,'h','k--','linewidth',2); ylim(y);
ylabel('mean +/ s.e.m. (z-units)');
end
clims([max(usedClims(:,1)) min(usedClims(:,2))])
sum(stats.above)
size(matrices{1})
size(ht)
size(qt)
ht(FindInterval(stats.below))
qt(FindInterval(stats.below))
if h>0, try plot(ht(FindInterval(stats.below)),ones(1,2)*min(ylim),'r','linewidth',5); end
try plot(ht(FindInterval(stats.above)),ones(1,2)*max(ylim),'r','linewidth',5); end
end
h>0
plot(ht(FindInterval(stats.below)),ones(1,2)*min(ylim),'r','linewidth',5);
try plot(qt(FindInterval(stats.below)),ones(1,2)*min(ylim),'r','linewidth',5); end
try plot(qt(FindInterval(stats.above)),ones(1,2)*max(ylim),'r','linewidth',5); end
h
[h,stats] = CompareDistributions(matrices{2},matrices{1});
if h>0, try plot(qt(FindInterval(stats.below)),ones(1,2)*min(ylim),'k','linewidth',5); end
try plot(qt(FindInterval(stats.above)),ones(1,2)*max(ylim),'k','linewidth',5); end
end
[h,stats] = CompareDistributions(matrices{2},matrices{1});
if h>0, try plot(qt(FindInterval(stats.below)),ones(1,2)*min(ylim),'k','linewidth',5); end
try plot(qt(FindInterval(stats.above)),ones(1,2)*max(ylim),'k','linewidth',5); end
h
h
[h,stats] = CompareDistributions(matrices{2},matrices{1});
h
q = matrices{3};
open q
shuffled = Scramble(matrices{3});
sum(shuffled(1,:))
sum(shuffled(2,:))
sum(shuffled(:,1))
sum(shuffled(:,2))
sum(q(:,1))
shuffled = Scramble(matrices{3}')';
shuffled = Shuffle(matrices{3});
sum(shuffled(:,1))
semplot(qt,shuffled)
cla
semplot(qt,shuffled)
figure; semplot(qt,shuffled)
rng(0); [h,stats] = CompareDistributions(matrices{3},Scramble(matrices{3}
rng(0); [h,stats] = CompareDistributions(matrices{3},Scramble(matrices{3}')');
h
5
h
qt(FindInterval(stats.above)
qt(FindInterval(stats.above))
plot(qt(FindInterval(stats.above)),ones(1,2)*max(ylim),'k','linewidth',5);
plot(qt(FindInterval(stats.above)),ones(size(FindInterval(stats.above)))*max(ylim),'k','linewidth',5);
clf
plot(qt(FindInterval(stats.above)),ones(size(FindInterval(stats.above)))*max(ylim),'k','linewidth',5);
plot(qt(FindInterval(stats.above))',ones(size(FindInterval(stats.above)))*max(ylim),'k','linewidth',5);
matrices = {z1,z2,z2-z1};
titles = {[name ' in pre-task sleep'],[name ' in post-task sleep'],'(post-pre) difference'};
colors = {'k','r','m'};
for j=1:3
subplot(2,3,j); PlotColorMap(Smooth(matrices{j},smooth),~isnan(matrices{j}),'x',qt); PlotHVLines(0,'v','k--','linewidth',2);
PlotHVLines(lines,'h','w--','linewidth',2);
set(gca,'ytick',mean([[0;lines(1:end-1)] lines],2),'yticklabel',sessions);
usedClims(j,:) = clim;
title(titles{j});
subplot(2,3,3+j); semplot(qt,matrices{j},colors{j},smooth(end));
if j==1, semplot(qt,matrices{j+1},colors{j+1},smooth(end)); y=ylim; end
ylim(y); PlotHVLines(0,'v','k--','linewidth',2); PlotHVLines(0,'h','k--','linewidth',2); ylim(y);
if j==3, rng(0); [h,stats] = CompareDistributions(matrices{3},Scramble(matrices{3}')');
if h>0, try plot(qt(FindInterval(stats.below))',ones(size(FindInterval(stats.below)))*min(ylim),'k','linewidth',5); end
try plot(qt(FindInterval(stats.above))',ones(size(FindInterval(stats.above)))*max(ylim),'k','linewidth',5); end
end
end
ylabel('mean +/ s.e.m. (z-units)');
end
clims([max(usedClims(:,1)) min(usedClims(:,2))])
drawnow
if h>0,
intervals = (FindInterval(stats.below));
if ~isempty(intervals), intervals(:,1) = intervals(:,1)-1; intervals(:,2) = intervals(:,2)+1;
plot(qt(intervals)',ones(size(intervals))'*min(ylim),'k','linewidth',5);
end
intervals = (FindInterval(stats.above));
if ~isempty(intervals), intervals(:,1) = intervals(:,1)-1; intervals(:,2) = intervals(:,2)+1;
plot(qt(intervals)',ones(size(intervals))'*max(ylim),'k','linewidth',5);
end
end
close all
strrep(strrep(name,':',''),' ','_')
condition
35
close all
cd M:\home\raly\results\PFC\reactivation
435
response = [mean(matrices{1}(:,InIntervals(qt,[0 0.5])),2) mean(matrices{2}(:,InIntervals(qt,[0 0.5])),2)];
clf
set(gcf,'position',[1 1 1920 1024])
matrices = {z1,z2,z2-z1};
titles = {[name ' in pre-task sleep'],[name ' in post-task sleep'],'(post-pre) difference'};
colors = {'k','r','m'};
for j=1:3
subplot(2,3,j); PlotColorMap(Smooth(matrices{j},smooth),~isnan(matrices{j}),'x',qt); PlotHVLines(0,'v','k--','linewidth',2);
PlotHVLines(lines,'h','w--','linewidth',2);
set(gca,'ytick',mean([[0;lines(1:end-1)] lines],2),'yticklabel',sessions);
usedClims(j,:) = clim;
title(titles{j});
subplot(2,3,3+j);
if j==2
response = [mean(matrices{1}(:,InIntervals(qt,[0 0.5])),2) mean(matrices{2}(:,InIntervals(qt,[0 0.5])),2)];
kruskalbar(response);
else
if j==1, semplot(qt,matrices{j+1},colors{j+1},smooth(end)); y=ylim; end
ylim(y); PlotHVLines(0,'v','k--','linewidth',2); PlotHVLines(0,'h','k--','linewidth',2); ylim(y);
if j==3, rng(0); [h,stats] = CompareDistributions(matrices{3},Scramble(matrices{3}')');
if h>0,
intervals = (FindInterval(stats.below));
if ~isempty(intervals), intervals(:,1) = intervals(:,1)-1; intervals(:,2) = intervals(:,2)+1;
plot(qt(intervals)',ones(size(intervals))'*min(ylim),'k','linewidth',5);
end
intervals = (FindInterval(stats.above));
if ~isempty(intervals), intervals(:,1) = intervals(:,1)-1; intervals(:,2) = intervals(:,2)+1;
plot(qt(intervals)',ones(size(intervals))'*max(ylim),'k','linewidth',5);
end
end
end
ylabel('mean +/ s.e.m. (z-units)');
end
end
clims([max(usedClims(:,1)) min(usedClims(:,2))])
drawnow
anovabar(response);
title(['p=' num2str(signrank(response(:,2)-response(:,1)))])
clf
set(gcf,'position',[1 1 1920 1024])
matrices = {z1,z2,z2-z1};
titles = {[name ' in pre-task sleep'],[name ' in post-task sleep'],'(post-pre) difference'};
colors = {'k','r','m'};
for j=1:3
subplot(2,3,j); PlotColorMap(Smooth(matrices{j},smooth),~isnan(matrices{j}),'x',qt); PlotHVLines(0,'v','k--','linewidth',2);
PlotHVLines(lines,'h','w--','linewidth',2);
set(gca,'ytick',mean([[0;lines(1:end-1)] lines],2),'yticklabel',sessions);
usedClims(j,:) = clim;
title(titles{j});
subplot(2,3,3+j);
if j==2
response = [mean(matrices{1}(:,InIntervals(qt,[0 0.5])),2) mean(matrices{2}(:,InIntervals(qt,[0 0.5])),2)];
kruskalbar(response);
ylabel('Median response 0-500ms (zz-units)');
title(['p=' num2str(signrank(response(:,2)-response(:,1)))]);
set(gca,'xticklabel',{'pre','post'});
else
if j==1, handle = semplot(qt,matrices{j+1},colors{j+1},smooth(end)); y=ylim; end
ylim(y); PlotHVLines(0,'v','k--','linewidth',2); PlotHVLines(0,'h','k--','linewidth',2); ylim(y);
%                 if j==3, rng(0); [h,stats] = CompareDistributions(matrices{3},Scramble(matrices{3}')');
%                     if h>0,
%                         intervals = (FindInterval(stats.below));
%                         if ~isempty(intervals), intervals(:,1) = intervals(:,1)-1; intervals(:,2) = intervals(:,2)+1;
%                             plot(qt(intervals)',ones(size(intervals))'*min(ylim),'k','linewidth',5);
%                         end
%                         intervals = (FindInterval(stats.above));
%                         if ~isempty(intervals), intervals(:,1) = intervals(:,1)-1; intervals(:,2) = intervals(:,2)+1;
%                             plot(qt(intervals)',ones(size(intervals))'*max(ylim),'k','linewidth',5);
%                         end
%                     end
%                 end
ylabel('mean +/ s.e.m. (z-units)');
end
end
clims([max(usedClims(:,1)) min(usedClims(:,2))])
handle = semplot(qt,matrices{j+1},colors{j+1},smooth(end));
open handle
handle
clf
set(gcf,'position',[1 1 1920 1024])
matrices = {z1,z2,z2-z1};
titles = {[name ' in pre-task sleep'],[name ' in post-task sleep'],'(post-pre) difference'};
colors = {'k','r','m'};
for j=1:3
subplot(2,3,j); PlotColorMap(Smooth(matrices{j},smooth),~isnan(matrices{j}),'x',qt); PlotHVLines(0,'v','k--','linewidth',2);
PlotHVLines(lines,'h','w--','linewidth',2);
set(gca,'ytick',mean([[0;lines(1:end-1)] lines],2),'yticklabel',sessions);
usedClims(j,:) = clim;
title(titles{j});
subplot(2,3,3+j);
if j==2
response = [mean(matrices{1}(:,InIntervals(qt,[0 0.5])),2) mean(matrices{2}(:,InIntervals(qt,[0 0.5])),2)];
kruskalbar(response);
ylabel('Median response 0-500ms (zz-units)');
title(['p=' num2str(signrank(response(:,2)-response(:,1)))]);
set(gca,'xticklabel',{'pre','post'});
else
handle = semplot(qt,matrices{j+1},colors{j+1},smooth(end)); y=ylim;
if j==1, handle2 = semplot(qt,matrices{j+1},colors{j+1},smooth(end)); y=ylim;
legend(handle,'pre',handle2,'post');
end
ylim(y); PlotHVLines(0,'v','k--','linewidth',2); PlotHVLines(0,'h','k--','linewidth',2); ylim(y);
%                 if j==3, rng(0); [h,stats] = CompareDistributions(matrices{3},Scramble(matrices{3}')');
%                     if h>0,
%                         intervals = (FindInterval(stats.below));
%                         if ~isempty(intervals), intervals(:,1) = intervals(:,1)-1; intervals(:,2) = intervals(:,2)+1;
%                             plot(qt(intervals)',ones(size(intervals))'*min(ylim),'k','linewidth',5);
%                         end
%                         intervals = (FindInterval(stats.above));
%                         if ~isempty(intervals), intervals(:,1) = intervals(:,1)-1; intervals(:,2) = intervals(:,2)+1;
%                             plot(qt(intervals)',ones(size(intervals))'*max(ylim),'k','linewidth',5);
%                         end
%                     end
%                 end
ylabel('mean +/ s.e.m. (z-units)');
end
end
clims([max(usedClims(:,1)) min(usedClims(:,2))])
handle
legend(handle,'pre');
legend(handle,'pre'); legend(handle2,'post');
handle = semplot(qt,matrices{j+1},colors{j+1},smooth(end)); y=ylim;
if j==1, handle2 = semplot(qt,matrices{j+1},colors{j+1},smooth(end)); y=ylim;
legend(handle,'pre'); legend(handle2,'post');
handle = semplot(qt,matrices{j},colors{j},smooth(end)); y=ylim;
handle = semplot(qt,matrices{j},colors{j},smooth(end)); y=ylim;
legend(handle,'pre'); legend(handle2,'post');
legend(handle,'pre',handle2,'post');
legend([handle handle2],{'pre','post'});
matrices = {z1,z2,z2-z1};
titles = {[name ' in pre-task sleep'],[name ' in post-task sleep'],'(post-pre) difference'};
colors = {'k','r','m'};
for j=1:3
subplot(2,3,j); PlotColorMap(Smooth(matrices{j},smooth),~isnan(matrices{j}),'x',qt); PlotHVLines(0,'v','k--','linewidth',2);
PlotHVLines(lines,'h','w--','linewidth',2);
set(gca,'ytick',mean([[0;lines(1:end-1)] lines],2),'yticklabel',sessions);
usedClims(j,:) = clim;
title(titles{j});
subplot(2,3,3+j);
if j==2
response = [mean(matrices{1}(:,InIntervals(qt,[0 0.5])),2) mean(matrices{2}(:,InIntervals(qt,[0 0.5])),2)];
kruskalbar(response);
ylabel('Median response 0-500ms (zz-units)');
title(['p=' num2str(signrank(response(:,2)-response(:,1)))]);
set(gca,'xticklabel',{'pre','post'});
else
handle = semplot(qt,matrices{j},colors{j},smooth(end)); y=ylim;
if j==1, handle2 = semplot(qt,matrices{j+1},colors{j+1},smooth(end)); y=ylim;
legend([handle handle2],{'pre','post'});
end
ylim(y); PlotHVLines(0,'v','k--','linewidth',2); PlotHVLines(0,'h','k--','linewidth',2); ylim(y);
%                 if j==3, rng(0); [h,stats] = CompareDistributions(matrices{3},Scramble(matrices{3}')');
%                     if h>0,
%                         intervals = (FindInterval(stats.below));
%                         if ~isempty(intervals), intervals(:,1) = intervals(:,1)-1; intervals(:,2) = intervals(:,2)+1;
%                             plot(qt(intervals)',ones(size(intervals))'*min(ylim),'k','linewidth',5);
%                         end
%                         intervals = (FindInterval(stats.above));
%                         if ~isempty(intervals), intervals(:,1) = intervals(:,1)-1; intervals(:,2) = intervals(:,2)+1;
%                             plot(qt(intervals)',ones(size(intervals))'*max(ylim),'k','linewidth',5);
%                         end
%                     end
%                 end
ylabel('mean +/ s.e.m. (z-units)');
end
end
clims([max(usedClims(:,1)) min(usedClims(:,2))])
drawnow
close all
2
Hide('none');
close all
profile on;
script_OJR_PFC_reactivation
profile viewer
13
close all
clim
clim([-3 3])
qq = {}; sessionIDs = {}; kk = 0; nSessions = size(X,1);
condition
X = XX{condition};
condition
i=1
[pfc,hpc,deltas,ripples,sleep,session,MergePoints,clickers,sws,basepath] = X{i,:};
[parentFolder,basename] = fileparts(basepath);
[~,projectName] = fileparts(parentFolder);
sessionID = [projectName '_' basename];
sessionIDs{i,1} = sessionID;
preSleep = SubtractIntervals(sleep(1,:), SubtractIntervals([0 Inf],sws));
postSleep = SubtractIntervals(sleep(2:end,:), SubtractIntervals([0 Inf],sws));
pre = InIntervals(ripples,preSleep);
post = InIntervals(ripples,postSleep);
%         % Restrict to 1h of sws only:
limit = Unshift(3600,preSleep); preSleep(preSleep(:,1)>limit,:) = []; if preSleep(end,2)>limit, preSleep(end,2) = limit; end
limit = Unshift(3600,postSleep); postSleep(postSleep(:,1)>limit,:) = []; if postSleep(end,2)>limit, postSleep(end,2) = limit; end
training = ConsolidateIntervals(MergePoints.timestamps(cellfun(@(x) ~isempty(strfind(x,'train')),MergePoints.foldernames),:)); % training epochs are named with "train" somewhere in the folder name
exploration = sortrows([clickers{2};clickers{3}]); exploration(any(isnan(exploration),2),:) = []; exploration = ConsolidateIntervals(exploration);
exploration = Restrict(exploration,training);
if structure==1, spikes = hpc; else, spikes = pfc; end
for j=1:max(spikes(:,2)) % use spiking activity
[q,qt] = PETH(spikes(spikes(:,2)==j),ripples(:,1),'durations',[-1 1]*2,'nBins',501);
kk = kk+1;
qq{kk,1} = nanmean(q(pre,:),1);
qq{kk,2} = nanmean(q(post,:),1);
qq{kk,3} = i;
%                             qqs{kk,1} = q;
end
i=2
[pfc,hpc,deltas,ripples,sleep,session,MergePoints,clickers,sws,basepath] = X{i,:};
[parentFolder,basename] = fileparts(basepath);
[~,projectName] = fileparts(parentFolder);
sessionID = [projectName '_' basename];
sessionIDs{i,1} = sessionID;
preSleep = SubtractIntervals(sleep(1,:), SubtractIntervals([0 Inf],sws));
postSleep = SubtractIntervals(sleep(2:end,:), SubtractIntervals([0 Inf],sws));
pre = InIntervals(ripples,preSleep);
post = InIntervals(ripples,postSleep);
%         % Restrict to 1h of sws only:
limit = Unshift(3600,preSleep); preSleep(preSleep(:,1)>limit,:) = []; if preSleep(end,2)>limit, preSleep(end,2) = limit; end
limit = Unshift(3600,postSleep); postSleep(postSleep(:,1)>limit,:) = []; if postSleep(end,2)>limit, postSleep(end,2) = limit; end
training = ConsolidateIntervals(MergePoints.timestamps(cellfun(@(x) ~isempty(strfind(x,'train')),MergePoints.foldernames),:)); % training epochs are named with "train" somewhere in the folder name
exploration = sortrows([clickers{2};clickers{3}]); exploration(any(isnan(exploration),2),:) = []; exploration = ConsolidateIntervals(exploration);
exploration = Restrict(exploration,training);
if structure==1, spikes = hpc; else, spikes = pfc; end
for j=1:max(spikes(:,2)) % use spiking activity
[q,qt] = PETH(spikes(spikes(:,2)==j),ripples(:,1),'durations',[-1 1]*2,'nBins',501);
kk = kk+1;
qq{kk,1} = nanmean(q(pre,:),1);
qq{kk,2} = nanmean(q(post,:),1);
qq{kk,3} = i;
%                             qqs{kk,1} = q;
end
figure; PETH(spikes(:,1),ripples(:,1),'durations',[-1 1]*2,'nBins',501);
i=3
[pfc,hpc,deltas,ripples,sleep,session,MergePoints,clickers,sws,basepath] = X{i,:};
[parentFolder,basename] = fileparts(basepath);
[~,projectName] = fileparts(parentFolder);
sessionID = [projectName '_' basename];
sessionIDs{i,1} = sessionID;
preSleep = SubtractIntervals(sleep(1,:), SubtractIntervals([0 Inf],sws));
postSleep = SubtractIntervals(sleep(2:end,:), SubtractIntervals([0 Inf],sws));
pre = InIntervals(ripples,preSleep);
post = InIntervals(ripples,postSleep);
%         % Restrict to 1h of sws only:
limit = Unshift(3600,preSleep); preSleep(preSleep(:,1)>limit,:) = []; if preSleep(end,2)>limit, preSleep(end,2) = limit; end
limit = Unshift(3600,postSleep); postSleep(postSleep(:,1)>limit,:) = []; if postSleep(end,2)>limit, postSleep(end,2) = limit; end
training = ConsolidateIntervals(MergePoints.timestamps(cellfun(@(x) ~isempty(strfind(x,'train')),MergePoints.foldernames),:)); % training epochs are named with "train" somewhere in the folder name
exploration = sortrows([clickers{2};clickers{3}]); exploration(any(isnan(exploration),2),:) = []; exploration = ConsolidateIntervals(exploration);
exploration = Restrict(exploration,training);
figure; PETH(spikes(:,1),ripples(:,1),'durations',[-1 1]*2,'nBins',501);
if structure==1, spikes = hpc; else, spikes = pfc; end
figure; PETH(spikes(:,1),ripples(:,1),'durations',[-1 1]*2,'nBins',501);
[h,ht] = PETH(spikes(:,1),ripples(:,1),'durations',[-1 1]*2,'nBins',501);
PlotColorMap(h,'x',ht)
PlotColorMap(Shrink(h,100,1),'x',ht)
PlotColorMap(Shrink(h,100,3),'x',ht)
PlotColorMap(Shrink(h,500,3),'x',ht)
PlotColorMap(Shrink(h(pre|post,:),500,3),'x',ht)
ripples = Restrict(ripples,sws);
pre = InIntervals(ripples,preSleep);
post = InIntervals(ripples,postSleep);
[h,ht] = PETH(spikes(:,1),ripples(:,1),'durations',[-1 1]*2,'nBins',501);
PlotColorMap(Shrink(h,500,3),'x',ht)
[hd,ht] = PETH(deltas(:,2),ripples(:,1),'durations',[-1 1]*2,'nBins',501);
PlotColorMap(Shrink(hd,500,3),'x',ht)
followed = InIntervals(ripples(:,2),deltas(:,2),[0 0.2]);
followed = InIntervals(ripples(:,2),bsxfun(@plus,deltas(:,2),[0 0.2]));
PlotColorMap(Shrink(h,500,3),'x',ht)
PlotColorMap(Shrink(sortby(h,followed),500,3),'x',ht)
followed = InIntervals(ripples(:,2),bsxfun(@plus,deltas(:,2),-[0 0.2]));
PlotColorMap(Shrink(sortby(h,followed),500,3),'x',ht)
Portion(followed)
followed = InIntervals(ripples(:,2),bsxfun(@plus,deltas(:,2),[-0.2 0]));
Portion(followed)
PlotColorMap(Shrink(sortby(h,followed),500,3),'x',ht)
PETH(spikes(:,1),ripples(:,2),'nbins',1001,'durations',[-1 1]*5)
PETH(spikes(:,1),ripples(:,2),'nbins',1001,'durations',[-1 1]*2)
PETH(spikes(:,1),ripples(:,1),'nbins',1001,'durations',[-1 1]*2)
PETH(spikes(:,1),ripples(:,1),'nbins',1001,'durations',[-1 1]*1)
PETH(spikes(:,1),ripples(:,1),'nbins',1001,'durations',[-1 1]*0.5)
PETH(spikes(:,1),ripples(,1),'nbins',1001,'durations',[-1 1]*0.5);
[h,ht] = PETH(spikes(:,1),ripples(,1),'nbins',1001,'durations',[-1 1]*0.5);
[h,ht] = PETH(spikes(:,1),ripples(:,1),'nbins',1001,'durations',[-1 1]*0.5);
semplot(ht,h(fol,:))
semplot(ht,h(followed,:))
semplot(ht,h(~followed,:))
clf
semplot(ht,hd)
[hd,hdt] = PETH(deltas(:,2),ripples(:,1),'durations',[-1 1]*2,'nBins',501);
semplot(hdt,hd)
[hd,hdt] = PETH(deltas(:,2),ripples(:,1),'durations',[-1 1]*0.5,'nBins',501);
semplot(hdt,hd)
clf
semplot(hdt,hd)
semplot(hdt,hd,'r',1)
semplot(ht,h,'b',1)
semplot(ht,zscore(h,[],2),'b',1)
clf
semplot(hdt,zscore(hd,[],2),'r',1)
semplot(ht,zscore(h,[],2),'b',1)
semplot(ht,zscore(h,[],2),'b',5)
clf
semplot(ht,zscore(h,[],2),'b',3)
semplot(ht,zscore(h(followed,:),[],2),'b',3)
semplot(ht,zscore(h(~followed,:),[],2),'r',3)
for j=1:max(spikes(:,2)) % use spiking activity
[q,qt] = PETH(spikes(spikes(:,2)==j),ripples(:,1),'durations',[-1 1]*2,'nBins',501);
kk = kk+1;
qq{kk,1} = nanmean(q(pre,:),1);
qq{kk,2} = nanmean(q(post,:),1);
qq{kk,3} = i;
qqs{kk,1} = q;
end
kk=0; qq = {];
kk=0; qq = {};
clear qqs
for j=1:max(spikes(:,2)) % use spiking activity
[q,qt] = PETH(spikes(spikes(:,2)==j),ripples(:,1),'durations',[-1 1]*2,'nBins',501);
kk = kk+1;
qq{kk,1} = nanmean(q(pre,:),1);
qq{kk,2} = nanmean(q(post,:),1);
qq{kk,3} = i;
qqs{kk,1} = q;
end
clf
q = (cell2mat(qq(:,1:2)));
PlotColorMap(q);
PlotColorMap(zscore(q,[],2));
semplot(zscore(q,[],2));
z = zscore(q,[],2);
z1 = z(:,1:length(qt)); z2 = z(:,(length(qt)+1):end);
clf
semplot(qt,z1,'k');
semplot(qt,z2,'r');
qqA = qq;
kk=0;
for j=1:max(spikes(:,2)) % use spiking activity
[q,qt] = PETH(spikes(spikes(:,2)==j),ripples(followed,1),'durations',[-1 1]*2,'nBins',501);
kk = kk+1;
qq{kk,1} = nanmean(q(pre,:),1);
qq{kk,2} = nanmean(q(post,:),1);
qq{kk,3} = i;
%                                                     qqs{kk,1} = q;
end
size(followed)
kk=0;
for j=1:max(spikes(:,2)) % use spiking activity
[q,qt] = PETH(spikes(spikes(:,2)==j),ripples(:,1),'durations',[-1 1]*2,'nBins',501);
kk = kk+1;
qq{kk,1} = nanmean(q(pre & followed,:),1);
qq{kk,2} = nanmean(q(post & followed,:),1);
qq{kk,3} = i;
%                                                     qqs{kk,1} = q;
end
qqF = qq;
z = zscore(q,[],2);
subplot(2,2,1);
semplot(qt,z1,'k');
semplot(qt,z2,'r');
subplot(2,2,2);
q = (cell2mat(qq(:,1:2)));
z = zscore(q,[],2);
z1 = z(:,1:length(qt)); z2 = z(:,(length(qt)+1):end);
semplot(qt,z1,'k');
semplot(qt,z2,'r');
for j=1:max(spikes(:,2)) % use spiking activity
[q,qt] = PETH(spikes(spikes(:,2)==j),ripples(:,1),'durations',[-1 1]*2,'nBins',501);
kk = kk+1;
qq{kk,1} = nanmean(q(pre & ~followed,:),1);
qq{kk,2} = nanmean(q(post & ~followed,:),1);
qq{kk,3} = i;
%                                                     qqs{kk,1} = q;
end
qqNF = qq;
subplot(2,2,3);
q = (cell2mat(qq(:,1:2)));
z = zscore(q,[],2);
z1 = z(:,1:length(qt)); z2 = z(:,(length(qt)+1):end);
semplot(qt,z1,'k');
semplot(qt,z2,'r');
open q
open qqF
clear qq
kk=0;
for j=1:max(spikes(:,2)) % use spiking activity
[q,qt] = PETH(spikes(spikes(:,2)==j),ripples(:,1),'durations',[-1 1]*2,'nBins',501);
kk = kk+1;
qq{kk,1} = nanmean(q(pre & ~followed,:),1);
qq{kk,2} = nanmean(q(post & ~followed,:),1);
qq{kk,3} = i;
%                                                     qqs{kk,1} = q;
end
qqNF = qq;
cla
q = (cell2mat(qq(:,1:2)));
z = zscore(q,[],2);
z1 = z(:,1:length(qt)); z2 = z(:,(length(qt)+1):end);
semplot(qt,z1,'k');
semplot(qt,z2,'r');
45
close all
52
j
deltaCondition
clfd
clf
EquateScales
SaveFig(fullfile('M:\home\raly\results\PFC\reactivation\semplots_delta_no_delta');
SaveFig(fullfile('M:\home\raly\results\PFC\reactivation\semplots_delta_no_delta'));
SaveFig(fullfile('M:\home\raly\results\PFC\reactivation\PFC_IC_reactivation_semplots_delta_no_delta'));
EquateScales
for i=1:9, subplot(3,3,i); PlotHVLines(0,'h','k--','linewidth',2);; end
for i=1:9, subplot(3,3,i); PlotHVLines(0,'v','k--','linewidth',2);; end
SaveFig(fullfile('M:\home\raly\results\PFC\reactivation\PFC_spikes_semplots_delta_no_delta'));
subplot(3,4,(deltaCondition-1)*3+condition);
clf
SaveFig(fullfile('M:\home\raly\results\PFC\reactivation\PFC_spikes_semplots_delta_no_delta'));
figure; hist(diff(ripples(:,2)),100);
hist(Restrict(diff(ripples(:,2)),[0 5]),100);
hist(Restrict(diff(ripples(:,2)),[0 1]),100);
[q,qt] = PETH(spikes(:,1),ripples(:,1),'durations',[-1 1]*2,'nBins',501);
figure; PlotColorMap(Shrink(q,100,1));
PlotColorMap(Shrink(q,100,1),'x',ht);
[pfc,hpc,deltas,ripples,sleep,session,MergePoints,clickers,sws,basepath] = X{i,:};
ripples = Restrict(ripples,sws);
[q,qt] = PETH(spikes(:,1),ripples(:,1),'durations',[-1 1]*2,'nBins',501);
PlotColorMap(Shrink(q,100,1),'x',ht);
Dist(100,Restrict(spikes,sws,'shift','on'),'grouped');
[h,ht] = Dist(100,Restrict(spikes,sws,'shift','on'),'grouped');
Dist(0:10:spikes(end,1),Restrict(spikes,sws,'shift','on'),'grouped');
PlotColorMap(h);
spikes(1,1)
max(spikes(:,1))
this = Restrict(spikes,sws,'shift','on');
Dist(0:10:this(end,1),this,'grouped');
[h,ht] = Dist(0:10:this(end,1),this,'grouped');
PlotColorMap(h);
PlotColorMap(h');
PlotColorMap(h'./sum(h,2));
PlotColorMap(h'./sum(h',2));
PlotColorMap(h'./sum(h',1));
[h,ht] = Dist(0:10:spikes(end,1),spikes,'grouped');
PlotColorMap(h'./sum(h',1));
sum(spikes(:,2)==34)
hist(spikes(spikes(:,2)==34))
hist(spikes(spikes(:,2)==34),100)
PlotIntervals(preSleep,'color','b');
PlotIntervals(postSleep,'color','b');
basepath
max(spikes(:,2))
ylim
PlotHVLines(241.5,'h','k--');
PlotHVLines(241+34 + .5,'h','k--');
basepath
edit CleanRez
SaveFig(fullfile('M:\home\raly\results\PFC\reactivation\HPC_spikes_semplots_delta_no_delta'));
2
SaveFig(fullfile('M:\home\raly\results\PFC\reactivation\PFC_IC_reactivation_semplots_delta_no_delta'));
[h,ht] = Dist(0:10:spikes(end,1),spikes,'grouped');
clf
PlotColorMap(double(h==0));
PlotColorMap(double(h==0)');
[h,ht] = Dist(0:100:spikes(end,1),spikes,'grouped');
PlotColorMap(double(h==0)');
[h,ht] = Dist(0:600:spikes(end,1),spikes,'grouped');
PlotColorMap(double(h==0)');
find(any(h==0,2))
find(any(h==0))
okUnits = ~any(h==0);
badUnits = any(h==0)';
newIDs = cumsum(~badUnits);
newIDs = cumsum(~badUnits); newIDs(badUnits) = 0;
script_OJR_PFC_reactivation
PlotColorMap(double(h==0)');
find(badIDs)
find(badUnits);
find(badUnits)
figure; hist(spikes(spikes(:,2)==10,1),1000);
hist(spikes(spikes(:,2)==10,1),1000);
hist(spikes(spikes(:,2)==13,1),1000);
hist(spikes(spikes(:,2)==10,1),1000);
PlotIntervals(postSleep,'color','b');
PlotIntervals(preSleep,'color','b');
dbcont
figure; hist(spikes(spikes(:,2)==10,1),1000);
find(badUnits)
clf
figure; hist(spikes(spikes(:,2)==7,1),1000);
clf
[q,qt] = PETH(spikes(spikes(:,2)==7,1),ripples(:,2));
PlotColorMap(q);
PlotColorMap(Shrink(q,100,1))
PlotColorMap(Shrink(q,100,1),'x',qt)
PlotColorMap(Shrink(q,500,1),'x',qt)
PlotColorMap(double(h==0)');
sum(h(:)==0)
size(h)
nUnits
max(spikes(:,2))
sum(h==0,2);
open ans
sum(h==0,1)';
badUnits = sum(h==0,1)'>1; % any empty 10-minute bins are signs of cells that should be excluded. Allow a single exception bin
sum(badUnits)
dbcont
sum(badUnits)
clf
PlotColorMap(double(h==0)');
find(badUnits)
hist(spikes(spikes(:,2)==3,1),1000);
hist(spikes(spikes(:,2)==9,1),1000);
hist(spikes(spikes(:,2)==1,1),1000);
dbcont
find(badUnits)
hist(spikes(spikes(:,2)==16,1),1000);
dbcont
find(badUnits)
hist(spikes(spikes(:,2)==6,1),1000);
hist(spikes(spikes(:,2)==12,1),1000);
hist(spikes(spikes(:,2)==13,1),1000);
hist(spikes(spikes(:,2)==24,1),1000);
hist(spikes(spikes(:,2)==26,1),1000);
basepath
dbcont
close all
nBad = sum(badUnits);
nBad)
nBad
ceil((nBad+1)/2)
basepath
nBad = sum(badUnits);
indices = find(bad);
figure;
subplot(ceil((nBad+1)/2),2,1);
PlotColorMap(double(h==0)');
for bad = 1:nBad
subplot(ceil((nBad+1)/2),2,bad+1);
hist(linspace(0,spikes(end,1),1000), spikes(spikes(:,2)==indices(bad),1));
end
title(basepath);
nBad = sum(badUnits);
indices = find(badUnits);
figure;
subplot(ceil((nBad+1)/2),2,1);
PlotColorMap(double(h==0)');
for bad = 1:nBad
subplot(ceil((nBad+1)/2),2,bad+1);
hist(linspace(0,spikes(end,1),1000), spikes(spikes(:,2)==indices(bad),1));
end
title(basepath);
basepath
nBad = sum(badUnits);
indices = find(badUnits);
figure;
subplot(ceil((nBad+1)/2),2,1);
PlotColorMap(double(h==0)');
for bad = 1:nBad
subplot(ceil((nBad+1)/2),2,bad+1);
hist(spikes(spikes(:,2)==indices(bad),1),linspace(0,spikes(end,1),1000));
end
title(basepath)
xlim([0 spikes(end,1)])
for bad = 1:nBad
subplot(ceil((nBad+1)/2),2,bad+1);
hist(spikes(spikes(:,2)==indices(bad),1),linspace(0,spikes(end,1),1000));
xlim([0 spikes(end,1)])
end
indices
close all
script_OJR_PFC_reactivation
subplot(ceil((nBad+1)/2),2,1);
title('Empty bins');
PlotColorMap(double(h==0)','x',ht);
xlabel('time (s)'); ylabel('Cell ID (HPC -> PFC)');
subplot(ceil((nBad+1)/2),2,1);
PlotColorMap(double(h==0)','x',ht);
title('Empty bins');
xlabel('time (s)'); ylabel('Cell ID (HPC -> PFC)');
title(['Empty bins of ' basepath]);
for bad = 1:nBad
subplot(ceil((nBad+1)/2),2,bad+1);
hist(spikes(spikes(:,2)==indices(bad),1),linspace(0,spikes(end,1),1000));
xlim([0 spikes(end,1)])
title(['Unit ' num2str(indices(bad)) ' from first plot']);
end
basepath
nBad = sum(badUnits);
indices = find(badUnits);
figure;
subplot(ceil((nBad+1)/2),2,1);
PlotColorMap(double(h==0)','x',ht);
title(['Empty bins of ' basepath]);
xlabel('time (s)'); ylabel('Cell ID (HPC -> PFC)');
for bad = 1:nBad
subplot(ceil((nBad+1)/2),2,bad+1);
hist(spikes(spikes(:,2)==indices(bad),1),linspace(0,spikes(end,1),1000));
xlim([0 spikes(end,1)])
title(['Unit ' num2str(indices(bad)) ' from first plot']);
end
spikes = [hpc; pfc(:,1) pfc(:,2)+max(hpc(:,2))];
[h,ht] = Dist(0:600:spikes(end,1),spikes,'grouped'); % bin spikes in 10-minute bins
badUnits = sum(h==0,1)'>1; % any em
basepath
nBad = sum(badUnits);
indices = find(badUnits);
figure;
subplot(ceil((nBad+1)/2),2,1);
PlotColorMap(double(h==0)','x',ht);
title(['Empty bins of ' basepath]);
xlabel('time (s)'); ylabel('Cell ID (HPC -> PFC)');
for bad = 1:nBad
subplot(ceil((nBad+1)/2),2,bad+1);
hist(spikes(spikes(:,2)==indices(bad),1),linspace(0,spikes(end,1),1000));
xlim([0 spikes(end,1)])
title(['Unit ' num2str(indices(bad)) ' from first plot']);
end
drawnow
spikes = [hpc; pfc(:,1) pfc(:,2)+max(hpc(:,2))];
1
close all
script_OJR_PFC_reactivation
close all
raly
deltaCondition = 2;
structure = 2
kind = 1
condition = 2;
X = XX{condition};
condition
qq = {}; sessionIDs = {}; kk = 0; nSessions = size(X,1);
cd M:\home\raly\results\PFC\reactivation
i=1
[pfc,hpc,deltas,ripples,sleep,session,MergePoints,clickers,sws,basepath] = X{i,:};
if isempty(pfc),pfc = zeros(0,2); end;  if isempty(hpc),hpc = zeros(0,2); end;
[parentFolder,basename] = fileparts(basepath);
[~,projectName] = fileparts(parentFolder);
sessionID = [projectName '_' basename];
sessionIDs{i,1} = sessionID;
preSleep = SubtractIntervals(sleep(1,:), SubtractIntervals([0 Inf],sws));
postSleep = SubtractIntervals(sleep(2:end,:), SubtractIntervals([0 Inf],sws));
if deltaCondition==1
preSleep = SubtractIntervals(sleep(1,:), sws);
postSleep = SubtractIntervals(sleep(2:end,:), sws);
end
try
if deltaCondition>2
followed = InIntervals(ripples(:,2),bsxfun(@plus,deltas(:,2),[-0.2 0]));
if deltaCondition==3
ripples(~followed,:) = [];
else
ripples(followed,:) = [];
end
end
catch
continue;
end
pre = InIntervals(ripples,preSleep);
post = InIntervals(ripples,postSleep);
% Restrict to 1h of sws only:
limit = Unshift(3600,preSleep); preSleep(preSleep(:,1)>limit,:) = []; if preSleep(end,2)>limit, preSleep(end,2) = limit; end
limit = Unshift(3600,postSleep); postSleep(postSleep(:,1)>limit,:) = []; if postSleep(end,2)>limit, postSleep(end,2) = limit; end
training = ConsolidateIntervals(MergePoints.timestamps(cellfun(@(x) ~isempty(strfind(x,'train')),MergePoints.foldernames),:)); % training epochs are named with "train" somewhere in the folder name
exploration = sortrows([clickers{2};clickers{3}]); exploration(any(isnan(exploration),2),:) = []; exploration = ConsolidateIntervals(exploration);
exploration = Restrict(exploration,training);
if structure==1, spikes = hpc; else, spikes = pfc; end
spikes = [hpc; pfc(:,1) pfc(:,2)+max(hpc(:,2))];
if isempty(spikes), continue; end
% remove lost units
[h,ht] = Dist(0:600:spikes(end,1),spikes,'grouped'); % bin spikes in 10-minute bins
badUnits = sum(h==0,1)'>1; % any empty 10-minute bins are signs of cells that should be excluded. Allow a single exception bin
newIDs = cumsum(~badUnits); newIDs(badUnits) = 0;
spikes(:,2) = newIDs(spikes(:,2));
spikes(spikes(:,2)==0,:) = []; % remove bad units
if kind==1
for j=1:max(spikes(:,2)) % use spiking activity
[q,qt] = PETH(spikes(spikes(:,2)==j),ripples(:,1),'durations',[-1 1]*2,'nBins',501);
kk = kk+1;
qq{kk,1} = nanmean(q(pre,:),1);
qq{kk,2} = nanmean(q(post,:),1);
qq{kk,3} = i;
%                                                     qqs{kk,1} = q;
end
elseif kind>1 % detect components and their reactivation
if kind==2, intervals = training; else, intervals = exploration; end
rng(0);
[templates,correlations,weights] = ActivityTemplatesICA(Restrict(spikes,intervals,'shift','on'),'binsize',0.1);
re = ReactivationStrength(spikes,templates,'step',0.01,'binSize',0.1);
for j=1:size(re,2)-1
[q,qt] = PETH(re(:,[1 1+j]),ripples(:,1),'durations',[-1 1]*2,'nBins',501);
kk = kk+1;
qq{kk,1} = nanmean(q(pre,:),1);
qq{kk,2} = nanmean(q(post,:),1);
qq{kk,3} = i;
%                             qqs{kk,1} = q;
end
end
if i==nSessions, qq(kk+1:end,:) = []; end
[pfc,hpc,deltas,ripples,sleep,session,MergePoints,clickers,sws,basepath] = X{i,:};
if isempty(pfc),pfc = zeros(0,2); end;  if isempty(hpc),hpc = zeros(0,2); end;
[parentFolder,basename] = fileparts(basepath);
[~,projectName] = fileparts(parentFolder);
sessionID = [projectName '_' basename];
sessionIDs{i,1} = sessionID;
preSleep = SubtractIntervals(sleep(1,:), SubtractIntervals([0 Inf],sws));
postSleep = SubtractIntervals(sleep(2:end,:), SubtractIntervals([0 Inf],sws));
if deltaCondition==1
preSleep = SubtractIntervals(sleep(1,:), sws);
postSleep = SubtractIntervals(sleep(2:end,:), sws);
end
if deltaCondition>2
followed = InIntervals(ripples(:,2),bsxfun(@plus,deltas(:,2),[-0.2 0]));
if deltaCondition==3
ripples(~followed,:) = [];
else
ripples(followed,:) = [];
end
end
pre = InIntervals(ripples,preSleep);
post = InIntervals(ripples,postSleep);
% Restrict to 1h of sws only:
limit = Unshift(3600,preSleep); preSleep(preSleep(:,1)>limit,:) = []; if preSleep(end,2)>limit, preSleep(end,2) = limit; end
limit = Unshift(3600,postSleep); postSleep(postSleep(:,1)>limit,:) = []; if postSleep(end,2)>limit, postSleep(end,2) = limit; end
training = ConsolidateIntervals(MergePoints.timestamps(cellfun(@(x) ~isempty(strfind(x,'train')),MergePoints.foldernames),:)); % training epochs are named with "train" somewhere in the folder name
exploration = sortrows([clickers{2};clickers{3}]); exploration(any(isnan(exploration),2),:) = []; exploration = ConsolidateIntervals(exploration);
exploration = Restrict(exploration,training);
if structure==1, spikes = hpc; else, spikes = pfc; end
[h,ht] = Dist(0:600:spikes(end,1),spikes,'grouped'); % bin spikes in 10-minute bins
badUnits = sum(h==0,1)'>1; % any empty 10-minute bins are signs of cells that should be excluded. Allow a single exception bin
newIDs = cumsum(~badUnits); newIDs(badUnits) = 0;
spikes(:,2) = newIDs(spikes(:,2));
spikes(spikes(:,2)==0,:) = []; % remove bad units
if structure==1, spikes = hpc; else, spikes = pfc; end
[h,ht] = Dist(0:600:spikes(end,1),spikes,'grouped'); % bin spikes in 10-minute bins
i=2
[pfc,hpc,deltas,ripples,sleep,session,MergePoints,clickers,sws,basepath] = X{i,:};
if isempty(pfc),pfc = zeros(0,2); end;  if isempty(hpc),hpc = zeros(0,2); end;
[parentFolder,basename] = fileparts(basepath);
[~,projectName] = fileparts(parentFolder);
sessionID = [projectName '_' basename];
sessionIDs{i,1} = sessionID;
preSleep = SubtractIntervals(sleep(1,:), SubtractIntervals([0 Inf],sws));
postSleep = SubtractIntervals(sleep(2:end,:), SubtractIntervals([0 Inf],sws));
if deltaCondition==1
preSleep = SubtractIntervals(sleep(1,:), sws);
postSleep = SubtractIntervals(sleep(2:end,:), sws);
end
try
if deltaCondition>2
followed = InIntervals(ripples(:,2),bsxfun(@plus,deltas(:,2),[-0.2 0]));
if deltaCondition==3
ripples(~followed,:) = [];
else
ripples(followed,:) = [];
end
end
catch
continue;
end
pre = InIntervals(ripples,preSleep);
post = InIntervals(ripples,postSleep);
% Restrict to 1h of sws only:
limit = Unshift(3600,preSleep); preSleep(preSleep(:,1)>limit,:) = []; if preSleep(end,2)>limit, preSleep(end,2) = limit; end
limit = Unshift(3600,postSleep); postSleep(postSleep(:,1)>limit,:) = []; if postSleep(end,2)>limit, postSleep(end,2) = limit; end
training = ConsolidateIntervals(MergePoints.timestamps(cellfun(@(x) ~isempty(strfind(x,'train')),MergePoints.foldernames),:)); % training epochs are named with "train" somewhere in the folder name
exploration = sortrows([clickers{2};clickers{3}]); exploration(any(isnan(exploration),2),:) = []; exploration = ConsolidateIntervals(exploration);
exploration = Restrict(exploration,training);
if structure==1, spikes = hpc; else, spikes = pfc; end
basepath
i=3
[pfc,hpc,deltas,ripples,sleep,session,MergePoints,clickers,sws,basepath] = X{i,:};
if isempty(pfc),pfc = zeros(0,2); end;  if isempty(hpc),hpc = zeros(0,2); end;
[parentFolder,basename] = fileparts(basepath);
[~,projectName] = fileparts(parentFolder);
sessionID = [projectName '_' basename];
sessionIDs{i,1} = sessionID;
preSleep = SubtractIntervals(sleep(1,:), SubtractIntervals([0 Inf],sws));
postSleep = SubtractIntervals(sleep(2:end,:), SubtractIntervals([0 Inf],sws));
if deltaCondition==1
preSleep = SubtractIntervals(sleep(1,:), sws);
postSleep = SubtractIntervals(sleep(2:end,:), sws);
end
if deltaCondition>2
followed = InIntervals(ripples(:,2),bsxfun(@plus,deltas(:,2),[-0.2 0]));
if deltaCondition==3
ripples(~followed,:) = [];
else
ripples(followed,:) = [];
end
end
pre = InIntervals(ripples,preSleep);
post = InIntervals(ripples,postSleep);
% Restrict to 1h of sws only:
limit = Unshift(3600,preSleep); preSleep(preSleep(:,1)>limit,:) = []; if preSleep(end,2)>limit, preSleep(end,2) = limit; end
limit = Unshift(3600,postSleep); postSleep(postSleep(:,1)>limit,:) = []; if postSleep(end,2)>limit, postSleep(end,2) = limit; end
training = ConsolidateIntervals(MergePoints.timestamps(cellfun(@(x) ~isempty(strfind(x,'train')),MergePoints.foldernames),:)); % training epochs are named with "train" somewhere in the folder name
exploration = sortrows([clickers{2};clickers{3}]); exploration(any(isnan(exploration),2),:) = []; exploration = ConsolidateIntervals(exploration);
exploration = Restrict(exploration,training);
if structure==1, spikes = hpc; else, spikes = pfc; end
[h,ht] = Dist(0:600:spikes(end,1),spikes,'grouped'); % bin spikes in 10-minute bins
badUnits = sum(h==0,1)'>1; % any empty 10-minute bins are signs of cells that should be excluded. Allow a single exception bin
newIDs = cumsum(~badUnits); newIDs(badUnits) = 0;
spikes(:,2) = newIDs(spikes(:,2));
spikes(spikes(:,2)==0,:) = []; % remove bad units
if kind==1
for j=1:max(spikes(:,2)) % use spiking activity
[q,qt] = PETH(spikes(spikes(:,2)==j),ripples(:,1),'durations',[-1 1]*2,'nBins',501);
kk = kk+1;
qq{kk,1} = nanmean(q(pre,:),1);
qq{kk,2} = nanmean(q(post,:),1);
qq{kk,3} = i;
%                                                     qqs{kk,1} = q;
end
elseif kind>1 % detect components and their reactivation
if kind==2, intervals = training; else, intervals = exploration; end
rng(0);
[templates,correlations,weights] = ActivityTemplatesICA(Restrict(spikes,intervals,'shift','on'),'binsize',0.1);
re = ReactivationStrength(spikes,templates,'step',0.01,'binSize',0.1);
for j=1:size(re,2)-1
[q,qt] = PETH(re(:,[1 1+j]),ripples(:,1),'durations',[-1 1]*2,'nBins',501);
kk = kk+1;
qq{kk,1} = nanmean(q(pre,:),1);
qq{kk,2} = nanmean(q(post,:),1);
qq{kk,3} = i;
%                             qqs{kk,1} = q;
end
end
if i==nSessions, qq(kk+1:end,:) = []; end
ripples = Restrict(ripples,sws);
PETH(spikes(:,1),ripples(:,1));
s = spikes(:,1);
r = ripples(:,1);
PETH(s,r);
PETH(s,r,'durations',[-1 1]*0.5);
PETH(s,r,'durations',[-1 1]*0.2);
[h,ht] = PlotColorMap(double(h==0)','x',ht);
[h,ht] = PETH(s,r,'durations',[-1 1]*0.2);
PlotColorMap(Shrink(h,500,1),'x',ht)
d = diff([0;r(:,1)]);
open d
PlotColorMap(Shrink(sortby(h,d),500,1),'x',ht)
PlotColorMap(Shrink(sortby(h,d),1000,1),'x',ht)
quantile(h,1-1/6)
quantile(d,1-1/6)
quantile(d,0.5)
semplot(ht,h);
clf
semplot(ht,h(d>0.5,:));
semplot(ht,h(d<0.5,:),'r');
clf
[hd,ht] = PETH(deltas(:,2),r,'durations',[-1 1]*0.2);
semplot(ht,h(d>0.5,:));
semplot(ht,h(d<0.5,:),'r');
clf
semplot(ht,hd(d>0.5,:));
semplot(ht,hd(d<0.5,:),'r');
clf
PlotColorMap(Shrink(sortby(hd,d),1000,1),'x',ht)
PlotColorMap(Smooth(Shrink(sortby(hd,d),1000,1),[0 2]),'x',ht)
PlotColorMap(Smooth(Shrink(sortby(hd,d),100,1),[0 2]),'x',ht)
PlotColorMap(Smooth(Shrink(sortby(hd,d),500,1),[0 2]),'x',ht)
[hd,ht] = PETH(deltas(:,2),r,'durations',[-1 1]*0.5);
PlotColorMap(Smooth(Shrink(sortby(hd,d),500,1),[0 2]),'x',ht)
PlotColorMap(Smooth(Shrink(sortby(hd,d),2000,1),[0 2]),'x',ht)
PlotColorMap(Smooth(Shrink(sortby(hd,d),2000,1),[0 2]),'x',ht); colormap(parula);
colormap(inferno)
colormap(hot)
colormap(inferno)
PlotColorMap(Smooth(Shrink(sortby(h,d),2000,1),[0 2]),'x',ht);
subplot(2,1,1); PlotColorMap(Smooth(Shrink(sortby(h,d),2000,1),[0 2]),'x',ht);
subplot(1,1,1); PlotColorMap(Smooth(Shrink(sortby(h,d),2000,1),[0 2]),'x',ht);
subplot(1,23,1); PlotColorMap(Smooth(Shrink(sortby(h,d),2000,1),[0 2]),'x',ht);
subplot(1,2,1); PlotColorMap(Smooth(Shrink(sortby(h,d),2000,1),[0 2]),'x',ht);
subplot(1,2,2); PlotColorMap(Smooth(Shrink(sortby(hd,d),2000,1),[0 2]),'x',ht);
colormap(inferno)
subplot(2,1,1); PlotColorMap(Smooth(Shrink(sortby(h,d),2000,1),[0 2]),'x',ht);
subplot(2,1,2); PlotColorMap(Smooth(Shrink(sortby(hd,d),2000,1),[0 2]),'x',ht);
followed = InIntervals(ripples(:,2),bsxfun(@plus,deltas(:,2),[-0.2 0]));
subplot(2,1,1); PlotColorMap(Smooth(Shrink(sortby(h,d + followed*10000),2000,1),[0 2]),'x',ht);
subplot(2,1,2); PlotColorMap(Smooth(Shrink(sortby(hd,d + followed*10000),2000,1),[0 2]),'x',ht);
subplot(2,1,1); PlotColorMap(Smooth(Shrink(sortby(h,d + followed*10000),1000,1),[0 2]),'x',ht);
subplot(2,1,2); PlotColorMap(Smooth(Shrink(sortby(hd,d + followed*10000),1000,1),[0 2]),'x',ht);
clim
clim([[0.2 0.4])
clim([[0.2 0.4]])
clim([[0.2 0.3]])
clim([[0.25 0.3]])
clim([[0.25 0.4]])
clim([[0.25 0.35]])
clim([[0.25 0.33]])
clim([[0.25 0.32]])
size(h)
ok = ~followed;
subplot(2,2,3); PlotColorMap(Smooth(Shrink(sortby(h(ok,:),d(ok)),2000,1),[0 2]),'x',ht);
clim
basepath
condition
size(ripples)
size(pre)
pre = InIntervals(ripples,preSleep);
post = InIntervals(ripples,postSleep);
clf
semplot(ht,h(pre & ~followed,:),'b');
semplot(ht,h(pre & ~followed,:),'b',5);
clf
semplot(ht,h(pre & ~followed,:),'b',2);
semplot(ht,h(post & ~followed,:),'r',2);
clear this
for i=1:max(spikes(:,2))
s = spikes(spikes(:,2)==this,1);
[h,ht] = PETH(s,r,'durations',[-1 1]*0.2);
for i=1:max(spikes(:,2))
s = spikes(spikes(:,2)==this,1);
[h,ht] = PETH(s,r,'durations',[-1 1]*0.5,'nBins',201);
this0(i,:) = nanmean(h(pre &~followed,:));
this1(i,:) = nanmean(h(post &~followed,:));
end
open this
raly
kkeyboard
for i=1:max(spikes(:,2))
s = spikes(spikes(:,2)==i,1);
[h,ht] = PETH(s,r,'durations',[-1 1]*0.5,'nBins',201);
this0(i,:) = nanmean(h(pre &~followed,:));
this1(i,:) = nanmean(h(post &~followed,:));
end
clf
PlotColorMap(this0)
PlotColorMap(this1)
PlotColorMap(this1-this0)
semplot(ht,this1-this0)
semplot(ht,zscore(this1,[],2)-zscore(this0,[],2);
semplot(ht,zscore(this1,[],2)-zscore(this0,[],2));
clf
z = [this0 this1]; z=zscore(z,[],2);
qt = ht;
z1 = z(:,1:length(qt)); z2 = z(:,(length(qt)+1):end);
semplot(ht,z1);
semplot(ht,z2,'r');
clf
semplot(ht,z2-z1,'m');
semplot(ht,z2-z1,'m',5);
clf
semplot(ht,z2-z1,'m',2);
semplot(ht,zscore(this1,[],2)-zscore(this0,[],2),'k');
semplot(ht,zscore(this1,[],2)-zscore(this0,[],2),'k',2);
clf
semplot(ht,z2-z1,'m',2);
semplot(ht,zscore(this1,[],2)-zscore(this0,[],2),'k',2);
response = [mean(this0(:,InIntervals(qt,[0 0.5])),2) mean(this0(:,InIntervals(qt,[0 0.5])),2)];
clf
anovabar(response)
plot(response)
response = [mean(this0(:,InIntervals(qt,[0 0.5])),2) mean(this1(:,InIntervals(qt,[0 0.5])),2)];
plot(response)
hist(diff(response,[],2))
hist(diff(response,[],2),100)
hist(diff(response,[],2),20)
signrank(diff(response,[],2))
ratio = diff(response,[],2)./sum(response,2);
hist(ratio)
hist(ratio,20)
signrank(ratio)
hist(ratio,50)
max(spikes(:,2))
find(badUnits)
[h,ht] = Dist(0:600:spikes(end,1),spikes,'grouped'); % bin spikes in 10-minute bins
badUnits = sum(h==0,1)'>1; % any empty 10-min
find(badUnits)
matrices = {zscore(this0,[],2),zscore(this1,[],2),zscore(this1,[],2)-zscore(this0,[],2)};
for i=1:3, subplot(1,3,i); PlotColorMap(matrices{i},'x',ht);
end
for i=1:3, subplot(1,3,i); PlotColorMap(matrices{i},'x',qt); end
clim
clims([-5 5]);
clims([-1 1]*2);
for i=1:3, subplot(1,3,i); PlotColorMap(Smooth(matrices{i},[0 2]),'x',qt); end
for i=1:3, subplot(1,3,i); PlotColorMap(Smooth(sortby(matrices{i},ratio),[0 2]),'x',qt); end
for i=1:3, subplot(1,3,i); PlotColorMap(Shrink(Smooth(sortby(matrices{i},ratio),[0 2]),5,1),'x',qt); end
response = [mean(this0(:,InIntervals(qt,[0 0.2])),2) mean(this1(:,InIntervals(qt,[0 0.2])),2)];
ratio = diff(response,[],2)./sum(response,2);
signrank(ratio)
signrank(diff(response,[],2))
for i=1:3, subplot(1,3,i); PlotColorMap(Shrink(Smooth(sortby(matrices{i},ratio),[0 2]),5,1),'x',qt); end
for i=1:3, subplot(1,3,i); PlotColorMap(Shrink(Smooth(sortby(matrices{i},ratio),[0 2]),5,1),'x',qt); PlotHVLines([0 0.2],'v','k--'); end
clim
clims([-1 1])
figure; hist(ratio,100)
cd M:\home\raly\results\PFC\reactivation
2
SaveFig(fullfile('M:\home\raly\results\PFC\reactivation\PFC_spikes_post_colormaps_delta_no_delta'));
clim
clims
clims([-3 3])
clims([-1 1])
clims([-1 1]*2)
SaveFig(fullfile('M:\home\raly\results\PFC\reactivation\PFC_spikes_post_colormaps_delta_no_delta'));
45
figure(3);
subplot(3,4,(condition-1)*4+deltaCondition);
matrices = {z1,z2,z2-z1};
response = [mean(matrices{1}(:,InIntervals(qt,[0 0.2])),2) mean(matrices{2}(:,InIntervals(qt,[0 0.2])),2)];
ratio = diff(response,[],2)./sum(response,2);
kruskalbar(response);
ylabel('Median response 0-200ms (z-units)');
title(['p=' num2str(signrank(response(:,2)-response(:,1)))]);
set(gca,'xticklabel',{'pre','post'});
figure(4);
subplot(3,4,(condition-1)*4+deltaCondition);
hist(ratio);
ylabel('Median response 0-200ms (z-units)');
title(['p=' num2str(signrank(ratio))]);
xlabel(['ratio, median =' num2str(median(ratio))]);
figure(4);
response = [mean(matrices{1}(:,InIntervals(qt,[0 0.2])),2) mean(matrices{2}(:,InIntervals(qt,[0 0.2])),2)];
d = diff(response,[],2)./sum(response,2);
subplot(3,4,(condition-1)*4+deltaCondition);
hist(d,100);
ylabel('Median response 0-200ms (z-units)');
title([name ': p=' num2str(signrank(ratio))]);
xlabel(['difference, median =' num2str(median(d))]);
d = diff(response,[],2);
subplot(3,4,(condition-1)*4+deltaCondition);
hist(d,100);
ylabel('Median response 0-200ms (z-units)');
title([name ': p=' num2str(signrank(ratio))]);
xlabel(['difference, median =' num2str(median(d))]);
hist(d,50);
ylabel('Median response 0-200ms (z-units)');
title([name ': p=' num2str(signrank(ratio))]);
xlabel(['difference, median =' num2str(median(d))]);
2
2
open dd
g = Group(dd{:});
open g
figure; anovabar(g(:,1),g(:,2));
kruskalbar(g(:,1),g(:,2));
gg = g; gg(:,2) = rem(gg(:,2)-1,4)+1;
open gg
gg = g; gg(:,2) = rem(g(:,2)-1,4)+1; gg(:,3) = floor(g(:,2)/4)+1;
figure; anovabar(gg(:,2),g(:,2));
anovabar(gg(:,2:3),g(:,2));
gg = g; gg(:,2) = floor(g(:,2)/3)+1; gg(:,4) = rem(g(:,2)-1,4)+1;
anovabar(gg(:,2:3),g(:,2));
gg = g; gg(:,2) = floor(g(:,2)/3)+1; gg(:,4) = rem(g(:,2)-1,3)+1;
anovabar(gg(:,2:3),g(:,2));
gg = g; gg(:,2) = floor(g(:,2)/3)+1; gg(:,3) = rem(g(:,2)-1,3)+1;
anovabar(gg(:,2:3),g(:,2));
gg = g; gg(:,2) = ceil(g(:,2)/3); gg(:,3) = rem(g(:,2)-1,3)+1;
anovabar(gg(:,2:3),g(:,2));
gg = g; gg(:,2) = ceil(g(:,2)/4); gg(:,3) = rem(g(:,2)-1,3)+1;
anovabar(gg(:,2:3),g(:,2));
gg = g; gg(:,2) = ceil(g(:,2)/4); gg(:,3) = rem(g(:,2)-1,4)+1;
anovabar(gg(:,2:3),g(:,2));
ok = gg(:,2)==1; anovabar(g(ok,1),gg(ok,3));
ok = gg(:,2)==1; kruskalbar(g(ok,1),gg(ok,3));
ok = gg(:,2)==2; kruskalbar(g(ok,1),gg(ok,3));
ok = gg(:,2)==3; kruskalbar(g(ok,1),gg(ok,3));
ok = gg(:,2)==3; anovabar(g(ok,1),gg(ok,3));
m = Accumulate(gg(:,2:3),g(:,1),'mode','mean');
bar(m);
bar(m');
for i=1:4,subplot(2,2,i); ok = gg(:,end)==i; anovabar(g(ok,1),gg(ok,2)); end
for i=1:4,subplot(2,2,i); ok = gg(:,end)==i; kruskalbar(g(ok,1),gg(ok,2)); end
deltaConditionNames = {'awake ripples','all sws ripples','ripples followed by deltas','sws ripples not followed by deltas'};
conditionNames = {'short','extended','opto'};
for i=1:4,subplot(2,2,i); ok = gg(:,end)==i; kruskalbar(g(ok,1),gg(ok,2));
title(deltaConditionNames{i});
set(gca,'xtick',conditionNames);
end
for i=1:4,subplot(2,2,i); ok = gg(:,end)==i; kruskalbar(g(ok,1),gg(ok,2));
title(deltaConditionNames{i});
set(gca,'XTickLabel'',conditionNames);
end
for i=1:4,subplot(2,2,i); ok = gg(:,end)==i; kruskalbar(g(ok,1),gg(ok,2));
title(deltaConditionNames{i});
set(gca,'xticklabel',conditionNames);
end
for i=1:4,subplot(2,2,i); ok = gg(:,end)==i;
anovabar(g(ok,1),gg(ok,2));
title(deltaConditionNames{i});
set(gca,'xticklabel',conditionNames);
end
for i=1:4,
subplot(2,2,i);
ok = gg(:,end)==i;
anovabar(g(ok,1),gg(ok,2));
title(deltaConditionNames{i});
ylabel('mean PFC response (0-200ms) to ripples (z-units)');
set(gca,'xticklabel',conditionNames);
end
SaveFig(fullfile('M:\home\raly\results\PFC\reactivation\PFC_spikes_delta_no_delta_differences_histograms'));
SaveFig(fullfile('M:\home\raly\results\PFC\reactivation\PFC_spikes_delta_no_delta_barplots'));
close all
open spikes
PETH(spikes(:,1),ripples(:,1));
PETH(pfc(:,1),ripples(:,1));
[h,ht] = PETH(pfc(:,1),ripples(:,1));
upstate = [deltas(1:end-1,2) deltas(2:end,2)];
hist(diff(upstate,[],2),100)
hist(Restrict(diff(upstate,[],2),[0 5]),100)
upstate = [deltas(1:end-1,2) deltas(2:end,2)]; upstate(diff(upstate,[],2)>5,:) = [];
rt = RelativeTime(ripples(:,1),upstate);
hist(rt,100)
PETH(ripples(:,1),deltas(:,2))
[pfc,hpc,deltas,ripples,sleep,session,MergePoints,clickers,sws,basepath] = X{i,:};
ripples = Restrict(ripples,sws);
[h,ht] = PETH(pfc(:,1),ripples(:,1));
rt = RelativeTime(ripples(:,1),upstate);
hist(rt,100)
hist(rt,0:0.1:1)
hist(rt,1000)
mmax(rt)
hist([rt; rt+1],100)
hist([rt; rt+1],1000)
PETH(ripples(:,1),deltas(:,2))
hist([rt; rt+1],1000)
[hd,ht] = PETH(ripples(:,1),deltas(:,2))
[hd,ht] = PETH(ripples(:,1),deltas(:,2));
PETH(hd)
PlotColorMap(Shrink(sortby(hd,d),1000,1),'x',ht)
PlotColorMap(Shrink(sortby(hd,r(:,1)),1000,1),'x',ht)
PlotColorMap(Shrink(sortby(hd,rt),1000,1),'x',ht)
open rt
upstate = [deltas(1:end-1,2) deltas(2:end,2)]; upstate(diff(upstate,[],2)>5,:) = [];
PETH(deltas(:,2),upstate(:,1))
PETH(ripples(:,1),upstate(:,1))
PETH(ripples(:,1),upstate(:,2))
v1
session = sessionTemplate(basepath,'showGUI',true);
hm = diff(deltas(:,2));
[q,qt] = PETH(ripples(:,1),upstate(:,2))
[q,qt] = PETH(ripples(:,1),upstate(:,2));
[q1,qt] = PETH(ripples(:,1),upstate(:,1));
PlotColorMap(q)
[q,qt] = PETH(ripples(:,1),upstate(:,2));
hm = diff(upstate,[],2);
PlotColorMap(Shrink(sortby(q,upstate(:,1)),1000,1),'x',qt)
PlotColorMap(Shrink(sortby(q,hm),1000,1),'x',qt)
figure; hist(hm,100);
PlotColorMap(Shrink(sortby(q,hm),1000,1),'x',qt)
upstate = [deltas(1:end-1,2) deltas(2:end,2)]; upstate(~InIntervals(diff(upstate,[],2),[0.5 5]),:) = [];
rt = RelativeTime(ripples(:,1),upstate);
hist([rt; rt+1],1000)
subplot(2,2,1); PlotColorMap(joint,'x',t2,'y',t1);
subplot(2,2,2); PlotColorMap(expected,'x',t2,'y',t1);
subplot(2,2,3); PlotColorMap(difference,'x',t2,'y',t1);
45
plot(m,1:max(ylim),'k');
open m
plot(ht(m),1:max(ylim),'k');
open m
345
open deltas
batchL = StartBatch(@BatchLoadHPCPFCData,'OJR_longTraining.batch');
Xl = get(batchL,'UserData'); % long
condition
basepath
X =Xl;
i=4;
[pfc,hpc,deltas,ripples,sleep,session,MergePoints,clickers,sws,basepath] = X{i,:};
ripples = Restrict(ripples,sws);
open Xl
BatchLoadHPCPFCStimData
Xl = get(batchL,'UserData'); % long
batchL = StartBatch(@BatchLoadHPCPFCData,'OJR_longTraining.batch');
batchL = StartBatch(@BatchLoadHPCPFCData,'OJR_longTraining.batch');
batchL = StartBatch(@BatchLoadHPCPFCData,'OJR_longTraining.batch');
Xl = get(batchL,'UserData'); % long
cd N:\OJRproject\OJR43\day10
load('day10.deltaWaves.events.mat')
X = Xo;
i=4;
[pfc,hpc,deltas,ripples,sleep,session,MergePoints,clickers,sws,basepath] = X{i,:};
ripples = Restrict(ripples,sws);
figure; semplot(ht,h);
X = Xl;
i=3;
i=5;
[pfc,hpc,deltas,ripples,sleep,session,MergePoints,clickers,sws,basepath] = X{i,:};
ripples = Restrict(ripples,sws);
i=4;
[pfc,hpc,deltas,ripples,sleep,session,MergePoints,clickers,sws,basepath] = X{i,:};
ripples = Restrict(ripples,sws);
i=3;
ripples = Restrict(ripples,sws);
open X
ripples = Restrict(ripples,sws);
[pfc,hpc,deltas,ripples,sleep,session,MergePoints,clickers,sws,basepath] = X{i,:};
ripples = Restrict(ripples,sws);
ripples = Restrict(ripples,sws);
i=4
condition=1;
X = XX{condition};
[pfc,hpc,deltas,ripples,sleep,session,MergePoints,clickers,sws,basepath] = X{i,:};
ripples = Restrict(ripples,sws);
upstate = [deltas(1:end-1,2) deltas(2:end,2)]; upstate(~InIntervals(diff(upstate,[],2),[1 4]),:) = [];
rt = RelativeTime(ripples(:,1),upstate);
[h,ht] = PETH(pfc(:,1),ripples(:,1),'durations',[-1 1]*0.5,'nBins',201);
[hd,ht] = PETH(deltas(:,2),ripples(:,1),'durations',[-1 1]*0.5,'nBins',201);
clim
clim(clim/2)
clim([0.1 0.4])
clim([0.15 0.4])
clim([0.2 0.4])
load('day3.deltaWaves.events.mat')
cd M:\home\raly\results\PFC\reactivation
ylabel('1')
Browse('all');
X = Xl;
open X
[pfc,hpc,deltas,ripples,sleep,session,MergePoints,clickers,sws,basepath] = X{i,:};
figure; PETH(pfc(:,1),ripples(:,1));
PETH(hpc(:,1),ripples(:,1));
PETH(hpc(:,1),ripples(:,1),'durations',[-1 1]*5);
PETH(pfc(:,1),ripples(:,1),'durations',[-1 1]*5);
PETH(deltas(:,2),ripples(:,1),'durations',[-1 1]*5);
PETH(deltas(:,2),ripples(:,1),'durations',[-1 1]*1);
open SWRpipeline.m
v1code
basepath
open try
v1
raly
5
figure; for i=1:4,
subplot(2,2,i);
ok = gg(:,end)==i;
anovabar(g(ok,1),gg(ok,2));
title(deltaConditionNames{i});
ylabel('mean PFC response (0-200ms) to ripples (z-units)');
set(gca,'xticklabel',conditionNames);
end
figure;    g = Group(dd{:}); gg = g; gg(:,2) = ceil(g(:,2)/4); gg(:,3) = rem(g(:,2)-1,4)+1;
for i=1:4,
subplot(2,2,i);
ok = gg(:,end)==i;
anovabar(g(ok,1),gg(ok,2));
title(deltaConditionNames{i});
ylabel('mean PFC response (0-200ms) to ripples (z-units)');
set(gca,'xticklabel',conditionNames);
end
SaveFig(fullfile('M:\home\raly\results\PFC\reactivation\PFC_spikes_delta_no_delta_differences_anovabars'));
set(gca,'xticklabel',conditionNames,'TickDir','out','box','off');
name
strrep(strrep(name,':',''),' ','_')
if structure==1
name = ['HPC'];
else
name = ['PFC'];
end
if kind==1
name = [name ' spikes'];
else
name = [name ' 100ms'];
if kind == 2
name = [name ' training'];
else
name = [name ' clicker object exploration'];
end
name = [name ' ICs'];
end
strrep(strrep(name,':',''),' ','_')
SaveFig(fullfile(['M:\home\raly\results\PFC\reactivation\PFC_responses_to_ripples_colormaps_ordered_by_proximity_to_delta']));
raly
%-- 3/11/2022 4:38 PM --%
v1
v1code
v1
v1code
script_launch_preprocessing
%-- 3/11/2022 6:05 PM --%
% try to fix ripples of OJR33
% either impose an additional HPC spiking threshold
% or even try to find bursts of activity instead of ripples
batchS = StartBatch(@BatchLoadHPCPFCData,'OJR_shortTraining.batch');
batchO = StartBatch(@BatchLoadHPCPFCData,'OJR_shortTraining_optoRipples_delayedPFC.batch');
batchL = StartBatch(@BatchLoadHPCPFCData,'OJR_longTraining.batch');
Xo = get(batchO,'UserData'); % opto
Xs = get(batchS,'UserData'); % short
Xl = get(batchL,'UserData'); % long
XX = {Xs,Xl,Xo};
write = false;
6
[pfc,hpc,deltas,ripples,sleep,session,MergePoints,clickers,sws,basepath] = X{i,:};
ripples = Restrict(ripples,sws);
upstate = [deltas(1:end-1,2) deltas(2:end,2)]; upstate(~InIntervals(diff(upstate,[],2),[1 4]),:) = [];
rt = RelativeTime(ripples(:,1),upstate);
clf
hist(rt,100)
response = mean(h(:,InIntervals(ht,[0 0.1])),2);
plot(rt,response,'.');
nancorr(rt,response)
b = ceil(rt/10);
anovabar(response,b);
meanControl = Accumulate(b,response,'mode','mean');
size(response)
size(b)
size(meanControl(b))
mmax(b)
sum(isnan(b))
56
semplotBin(b(ok),response(ok))
semplotBin(b(ok),response(ok),10)
semplotBin(rt(ok),response(ok),10);
clf
semplotBin(rt(ok),response(ok),10);
semplotBin(rt(ok),response(ok),100);
clf
semplotBin([rt(ok); rt(ok)+1],repmat(response(ok),2,1),100);
semplotBin([rt(ok); rt(ok)+1],repmat(response(ok),2,1),100,2);
semplotBin([rt(ok); rt(ok)+1],repmat(response(ok),2,1),100,'r',2);
clf
semplotBin([rt(ok); rt(ok)+1],repmat(response(ok),2,1),100,'r',2);
semplotBin([rtC; rtC+1],repmat(responseC,2,1),100,'k',2);
clf
semplotBin([rt(ok); rt(ok)+1],repmat(response(ok),2,1),100,'r',2);
semplotBin([rtC; rtC+1],repmat(responseC,2,1),100,'k',2);
34
legend(handle,'r')
legend([handle handle0],{'data','control'});
(condition-1)*4
i
Xs(cellfun(@isempty,Xs(:,1)),:) = []; Xl(cellfun(@isempty,Xl(:,1)),:) = []; Xo(cellfun(@isempty,Xo(:,1)),:) = [];
XX = {Xs,Xl,Xo};
batchS = StartBatch(@BatchLoadHPCPFCData,'OJR_shortTraining.batch');
Xs = get(batchS,'UserData'); % short
Xs(cellfun(@isempty,Xs(:,1)),:) = []; Xl(cellfun(@isempty,Xl(:,1)),:) = []; Xo(cellfun(@isempty,Xo(:,1)),:) = [];
XX = {Xs,Xl,Xo};
set(gca,'xtick',-1:0.2:1,'xticklabel',{'delta','','','','','delta','','','','','delta'});
PlotHVLines(0,'v','k--','linewidth',2);
legend([handle handle0],{'ripples','random control'}); legend('location','best');
set(gca,'xtick',-1:0.2:1,'xticklabel',{'delta','early upstate','','','late upstate','delta','','','','','delta'});
set(gca,'xtick',-1:1/3:1,'xticklabel',{'delta','early upstate','mid-upstate','late upstate','delta','early upstate','mid-upstate','late upstate','delta'});
set(gca,'xtick',-1:1/4:1,'xticklabel',{'delta','early upstate','mid-upstate','late upstate','delta','early upstate','mid-upstate','late upstate','delta'});
clim
clim([0 2
])
clim([0.2 1.5])
clim([0.2 3])
clim([0.5 3])
clim([0.7 3])
clim([1 3])
PlotHVLines(0,'v','k--','linewidth',2);
PlotHVLines(0.2,'v','k--','linewidth',2);
PlotHVLines(0.1,'v','k--','linewidth',2);
legend('location','southeast','box','off');
legend([handle handle0],{'100ms following ripples','following random control'}); legend('location','southwest','box','off');
SaveFig(fullfile(['M:\home\raly\results\PFC\reactivation\PFC_responses_to_ripples_colormaps_ordered_by_relative_upstate_time']));
X = Xl;
i=1;
[pfc,hpc,deltas,ripples,sleep,session,MergePoints,clickers,sws,basepath] = X{i,:};
basepath
figure; PETH(deltas(:,2),ripples(:,1));
PETH(pfc(:,1),ripples(:,1));
hold all
PETH(deltas(:,2),ripples(:,1));
clf
[hd,ht] = PETH(deltas(:,2),ripples(:,1),'durations',[-1 1]*0.5,'nBins',201);
[h,ht] = PETH(pfc(:,1),ripples(:,1),'durations',[-1 1]*0.5,'nBins',201);
semplot(ht,zscore(h,[],2),'r');
semplot(ht,zscore(hd,[],2),'b');
[hh,ht] = PETH(hpc(:,1),ripples(:,1),'durations',[-1 1]*0.5,'nBins',201);
semplot(ht,zscore(hh,[],2),'j');
semplot(ht,zscore(hh,[],2),'k');
X = Xs;
i=1;
[pfc,hpc,deltas,ripples,sleep,session,MergePoints,clickers,sws,basepath] = X{i,:};
[hh,ht] = PETH(hpc(:,1),ripples(:,1),'durations',[-1 1]*0.5,'nBins',201);
[h,ht] = PETH(pfc(:,1),ripples(:,1),'durations',[-1 1]*0.5,'nBins',201);
[hd,ht] = PETH(deltas(:,2),ripples(:,1),'durations',[-1 1]*0.5,'nBins',201);
figure; semplot(ht,zscore(h,[],2),'r');
semplot(ht,zscore(hd,[],2),'b');
semplot(ht,zscore(hh,[],2),'k');
clf
PETH(pfc(:,1),deltas(:,2))
subplot(2,2,1);
semplot(ht,zscore(h,[],2),'r');
xlabel('time around ripples (s)');
ylabel('HPC multiunit activity (z-scored)');
set(gca,'fontsize',15,'TickDir','out');
subplot(2,2,2);
[hp,ht] = PETH(pfc(:,1),deltas(:,2),'durations',[-1 1]*0.5,'nBins',201);
semplot(ht,zscore(hp,[],2),'b');
xlabel('time around delta waves (s)');
ylabel('PFC multiunit activity (z-scored)');
set(gca,'fontsize',15,'TickDir','out');
subplot(2,2,3);
semplot(ht,zscore(hd,[],2),'b');
semplot(ht,hd/mode(diff(ht)),'k');
cla
semplot(ht,hd/mode(diff(ht)),'k');
basepath
X = Xl;
i = ;1
i=1;
[pfc,hpc,deltas,ripples,sleep,session,MergePoints,clickers,sws,basepath] = X{i,:};
clf
PETH(pfc(:,1),deltas(:,2))
subplot(2,2,1);
semplot(ht,zscore(h,[],2),'r');
xlabel('time around ripples (s)');
ylabel('HPC multiunit activity (z-scored)');
set(gca,'fontsize',15,'TickDir','out');
subplot(2,2,2);
[hp,ht] = PETH(pfc(:,1),deltas(:,2),'durations',[-1 1]*0.5,'nBins',201);
semplot(ht,zscore(hp,[],2),'b');
xlabel('time around delta waves (s)');
ylabel('PFC multiunit activity (z-scored)');
set(gca,'fontsize',15,'TickDir','out');
subplot(2,2,3);
semplot(ht,hd/mode(diff(ht)),'k');
[hd,ht] = PETH(deltas(:,2),ripples(:,1),'durations',[-1 1]*0.5,'nBins',201);
[hp,ht] = PETH(pfc(:,1),deltas(:,2),'durations',[-1 1]*0.5,'nBins',201);
[hh,ht] = PETH(hpc(:,1),ripples(:,1),'durations',[-1 1]*0.5,'nBins',201);
clf
subplot(2,2,1);
semplot(ht,zscore(h,[],2),'r');
xlabel('time around ripples (s)');
ylabel('HPC multiunit activity (z-scored)');
set(gca,'fontsize',15,'TickDir','out');
subplot(2,2,2);
[hp,ht] = PETH(pfc(:,1),deltas(:,2),'durations',[-1 1]*0.5,'nBins',201);
semplot(ht,zscore(hp,[],2),'b');
xlabel('time around delta waves (s)');
ylabel('PFC multiunit activity (z-scored)');
set(gca,'fontsize',15,'TickDir','out');
subplot(2,2,3);
semplot(ht,hd/mode(diff(ht)),'k');
cla
semplot(ht,hd/mode(diff(ht)),'k',1);
ylabel('delta wave rate (Hz)'
ylabel('delta wave rate (Hz)')
xlabel('time around ripples (s)');
set(gca,'fontsize',15,'TickDir','out');
basepath
X = Xs;
i=size(X,1);
[pfc,hpc,deltas,ripples,sleep,session,MergePoints,clickers,sws,basepath] = X{i,:};
clf
[h,ht] = PETH(pfc(:,1),ripples(:,1),'durations',[-1 1]*0.5,'nBins',201);
PlotColorMap(h);
d = ripples(:,1);
PlotColorMap(Smooth(Shrink(sortby(h,d),2000,1),[0 2]),'x',ht);
PlotColorMap(Smooth(Shrink(sortby(h,d),200,1),[0 2]),'x',ht);
PlotColorMap(Smooth(Shrink(sortby(h,d),20,1),[0 2]),'x',ht);
ripples = Restrict(ripples,sws);
[h,ht] = PETH(pfc(:,1),ripples(:,1),'durations',[-1 1]*0.5,'nBins',201);
PlotColorMap(Smooth(Shrink(sortby(h,d),20,1),[0 2]),'x',ht);
d = ripples(:,1);
PlotColorMap(Smooth(Shrink(sortby(h,d),20,1),[0 2]),'x',ht);
d = ripples(:,3)-ripples(:,1);
PlotColorMap(Smooth(Shrink(sortby(h,d),20,1),[0 2]),'x',ht);
PlotColorMap(Smooth(Shrink(sortby(h,d),50,1),[0 2]),'x',ht);
PlotColorMap(Smooth(Shrink(sortby(h,d),100,1),[0 2]),'x',ht);
PlotHVLines(0,'v','k--','linewidth',2);
plot(Shrink(sortby(d,d),100,1),1:max(ylim),'w','linewidth',2);
clim
clim([0 0.8]);
clim([0.2 0.8]);
clim([0.2 0.7]);
clim([0.2 0.8]);
count = CountInIntervals(pfc(:,1),ripples(:,[1 3]))./diff(ripples(:,[1 3]));
count = CountInIntervals(pfc(:,1),ripples(:,[1 3]))./diff(ripples(:,[1 3]),[],2);
open count
figure; plot(d,count);
plot(d,count,'.');
nancorr(d,count)
bins = Bins(0,MergePoints.timestamps(end,end),0.1);
bins = Bins(0,MergePoints.timestamps(end,end),0.1,0.01);
smoothed = CountInIntervals(pfc(:,1),bins);
open smoothed
t = mean(bins,2);
PETH([t smoothed],ripples(:,1))
PETH(pfc(:,1),ripples(:,1))
PETH(pfc(:,1),ripples(:,1),'durations',[-1 1]*0.2)
PETH([t smoothed],ripples(:,1),'durations',[-1 1]*0.2)
bins = Bins(0,MergePoints.timestamps(end,end),0.01,0.01);
PETH([t smoothed],ripples(:,1),'durations',[-1 1]*0.2)
smoothed = CountInIntervals(pfc(:,1),bins);
t = mean(bins,2);
PETH([t smoothed],ripples(:,1),'durations',[-1 1]*0.2)
[in,w] = InIntervals(t,ripples(:,[1 3]));
fr = Accumulate(w(in),smoothed(in));
nancorr(d,fr)
plot(d,fr,'.');
fr = Accumulate(w(in),smoothed(in),'mode','mean');
nancorr(d,fr)
plot(d,fr,'.');
DensityMap(d,fr);
DensityMap(tiedrank(d,fr);)
DensityMap(tiedrank(d),fr);
DensityMap(tiedrank(d),tiedrank(fr));
DensityMap(d,fr);
DensityMap(d,fr,'smooth',0);
DensityMap(d,fr,'smooth',0,'nBins',50);
DensityMap(d,fr,'smooth',0,'nBins',100);
DensityMap(d,fr,'smooth',5,'nBins',100);
DensityMap(d,fr,'smooth',[0 5],'nBins',100);
DensityMap(d,fr,'smooth',[5 0],'nBins',100);
DensityMap(d,fr,'smooth',[5 1],'nBins',100);
DensityMap(d,fr,'smooth',[5 1],'nBins',500);
DensityMap(d,fr,'smooth',[5 1]*10,'nBins',500);
DensityMap(d,fr,'smooth',[5 1]*5,'nBins',500);
clim
clim([0 0.2])
DensityMap(d,tiedrank(fr),'smooth',[5 1]*5,'nBins',500);
DensityMap(d,tiedrank(fr),'smooth',[5 5]*5,'nBins',500);
DensityMap(d,tiedrank(fr),'smooth',[5 2]*5,'nBins',500);
max(d)
DensityMap(d/0.022,tiedrank(fr),'smooth',[5 2]*5,'nBins',500);
DensityMap(d/0.022,tiedrank(fr),'smooth',[5 2]*5,'nBins',50);
DensityMap(d/0.022,tiedrank(fr),'smooth',[5 2],'nBins',50);
fr = Accumulate(w(in),smoothed(in),'mode','max');
nancorr(d,fr)
figure; plot(d,fr);
plot(d,fr,'.');
smoothed = Smooth(smoothed,1);
fr = Accumulate(w(in),smoothed(in),'mode','max');
plot(d,fr,'.');
nancorr(d,fr)
fr = Accumulate(w(in),smoothed(in),'mode','mean');
nancorr(d,fr)
plot(d,fr,'.');
[in,w] = InIntervals(t,[ripples(:,2)-0.05 ripples(:,2)+0.05]);
fr = Accumulate(w(in),smoothed(in),'mode','mean');
plot(d,fr,'.');
nancorr(d,fr)
[c,p] = nancorr(d,fr)
[in,w] = InIntervals(t,[ripples(:,2)-0.05 ripples(:,2)+0.2]);
[in,w] = InIntervals(t,[ripples(:,1) ripples(:,1)+0.1]);
fr = Accumulate(w(in),smoothed(in),'mode','mean');
nancorr(d,fr)
plot(d,fr,'.');
[in,w] = InIntervals(t,[ripples(:,2)-0.02 ripples(:,2)+0.05]);
fr = Accumulate(w(in),smoothed(in),'mode','mean');
plot(d,fr,'.');
[c,p] = nancorr(d,fr)
cd M:\home\raly\results\PFC\reactivation
kkeyboard
v1
screenid = 2;
% AssertOpenGL; % Configure OpenGL Psychtoolbox:
% Screen('Preference', 'Verbosity', 0)
% Screen('Preference', 'VisualDebugLevel', 0); % don't show start screen\
% Screen('Preference', 'SkipSyncTests', 1); % skip some tests that may trigger code errors
PsychImaging('PrepareConfiguration');
PsychImaging('AddTask', 'AllViews', 'EnableCLUTMapping'); % enable screen gamma changes is necessary for the left cue to move
Screen('Preference', 'VisualDebugLevel', 0); % don't show start screen\
Screen('Preference', 'SkipSyncTests', 1);
win = PsychImaging('OpenWindow', screenid, 0); % Open a fullscreen onscreen window on that display
[width, height] = Screen('WindowSize', win);
helper_YmazeGreen_prepare_left
% Draw same image into backbuffer, so they're identical:
% Screen('DrawTexture', win, left_texture);
display(1);
pause(0.01)
abortScript = false;
close all
sca
t0 = GetSecs;
t0
diff(t0)
diff(t1)
t0
open t0
open t1
diff(t0(1:2))
diff([6.46937859995524,6.47353939997265])
diff([6.46937859995524,6.47353939997265])*1000
open t1
t1 - t0(end)
open t
diff(t)
t = diff([t0(:) t1(:)]);
profile on
temp
profile viewer
raly
v1
temp
profile on
temp
profile on
temp
profile viewer
profile on
temp
profile viewer
profile on
temp
profile on
temp
open origLUT
figure; PlotColorMap(origLUT)
win = PsychImaging('OpenWindow', screenid, 0);
Screen('FillRect',win, 0);
Screen('DrawTexture', win, left_texture(1));
PsychImaging('PrepareConfiguration');
PsychImaging('AddTask', 'AllViews', 'EnableCLUTMapping');
win = PsychImaging('OpenWindow', screenid, 0);
LoadIdentityClut(win);
if size(origLUT, 1) ~= 256
origLUT = origLUT(1:256, :);
end
newLUT = origLUT; newLUT(1,:)=0.5;
left_s = floor(min(width, height)/2)*10;
left_x = meshgrid(-left_s:left_s, -left_s:left_s);
left_fintex = ones(2*left_s+1,2*left_s+1);
left_fintex(:,:) = mod(left_x,255)+1;
left_texture = Screen('MakeTexture', win, left_fintex);
% Black background:
Screen('FillRect',win, 0);
% Single static gray-level ramp drawn as texture.
Screen('DrawTexture', win, left_texture(1));
Screen('Flip', win);
sca
figure; PlotColorMap(left_s)
open left_s
PlotColorMap(left_x)
PlotColorMap(left_fintex)
v1
open t
temp
PlotColorMap(left_fintex)
sca
PlotColorMap(newLOT);
PlotColorMap(newLUT);
PlotColorMap(newLUT); colormap(bone);
PlotColorMap(origLUT); colormap(bone);
sca
open t0
diff(t0)
PlotColorMap(left_fintex);
left_fintex = circshift(left_fintex',speed)';
PlotColorMap(left_fintex);
PlotColorMap(circshift(left_fintex',speed)');
PlotColorMap(circshift(left_fintex',1)');
PlotColorMap(circshift(left_fintex,1));
left_fintex = circshift(left_fintex,1);
left_fintex = circshift(left_fintex',1)';
PlotColorMap(circshift(left_fintex',1)');
PlotColorMap(circshift(left_fintex',100)');
PlotColorMap(circshift(left_fintex',1000)');
floor(min(width, height)/2)*10
-left_s
left_s
profile on
temptemp
profile viewer
mmax(left_findex)
mmax(left_fintex)
left_fintex = mod(left_fintex+speed-1,255)+1;
profile on
temptemp
profile viewer
v1
which Psychtoolbox
which PsychLinuxConfiguration
cd C:\Users\Cornell\Documents\MATLAB\Psychtoolbox\
DriftingMaskedGratingTutorial;
DriftingMaskedGratingTutorial(1,1,1,1);
MovingLineDemo
newLUT = origLUT; newLUT(1,:)=0.5;
left_s = max(width, height);
left_x = meshgrid(-left_s:left_s, -left_s:left_s);
left_fintex = ones(2*left_s+1,2*left_s+1);
left_fintex(:,:) = mod(left_x,255)+1;
find(left_findtex(1,:)==121,5)
find(left_fintex(1,:)==121,5)
current
profile on
temptemp
profile viewer
cd 'M:\home\raly\temp'
profile on
temptemp
profile viewer
profile viewer
profile on
temptemp
profile viewer
profile viewer
display(current);
display(num2str(current));
temp
sca
temptemp
profile on
temp
profile viewer
profile viewer
LoadIdentityClut
profile on
temp
profile viewer
profile viewer
open t
t = diff([t0(:) t1(:)]);
v1
which temp
profile on
temp
profile viewer
keyCode
find(keyCode)
abortScript
sc
sca
profile on
temp
profile viewer
which temp
edit temp
profile on
temp
profile viewer
clear all
clc
close all
load('day10.animal.behavior.mat', 'behavior')
figure; plot(behavior.position.x,behavior.position.y);
cd N
cd N:
PlotXY(behavior.positionTrials{1})
PlotXY(behavior.positionTrials{1},'.')
hold all
PlotXY(behavior.positionTrials{2},'r.')
hold all
plot(behavior.timestamps,behavior.speed,'k');
PlotIntervals(behavior.trials);
plot(behavior.timestamps,behavior.speed,'y');
hold all
PlotXY(behavior.positionTrialsRun{1},'k.')
run = behavior.timestamps(FindInterval(behavior.speed>0.02));
ylim([0 100])
ylim([0 10])
ylim([0 1])
PlotIntervals(run);
PlotIntervals(run,'color','g');
[w,in] = InIntervals(behavior.timestamps(:),run);
peak = Accumulate(w(in),behavior.speed(in));
peak = Accumulate(w(in),behavior.speed(in)','mode','max');
[in,w] = InIntervals(behavior.timestamps(:),run);
peak = Accumulate(w(in),behavior.speed(in)','mode','max');
opean peak
open peak
clf
hist(peak,100)
hist(peak,10000)
hist(Restrict(peak,[0 10]),10000)
plot(behavior.timestamps,behavior.speed,'y');
ylim([0 1])
run = behavior.timestamps(FindInterval(behavior.speed>0.02));
[in,w] = InIntervals(behavior.timestamps(:),run);
peak = Accumulate(w(in),behavior.speed(in)','mode','max');
run(peak<0.2,:) = [];
clf
PlotXY(behavior.positionTrialsRun{1},'k.')
PlotXY(behavior.positionTrials{1},'.')
hold all
PlotXY(behavior.positionTrials{2},'r.')
PlotIntervals(run,'color','k');
raly
hist(peak,100)
hist(Restrict(peak,[0 10]),10000)
hist(Restrict(peak,[0 10]),1000)
PlotXY(behavior.positionTrials{2},'r.')
findmin(peak(peak>0.2))
sum(peak>0.2)
PlotIntervals(run(70,:),'color','k');
run(70,:)
xlim(ans + [-1 1]*10)p;
xlim(ans + [-1 1]*10);
xlim(ans + [-1 1]*20);
xlim(ans + [-1 1]*200);
hold all
PlotXY(behavior.positionTrials{1});
run = behavior.timestamps(FindInterval(behavior.speed>0.02));
d = run(2:end,1) - run(1:end-1,2);
hist(d,100)
hist(Restrict(d,[0 10]),100)
hist(Restrict(d,[0 10]),1000)
figure; PlotXY(behavior.positionTrials{1},'.')
hold all
PlotXY(behavior.positionTrials{2},'r.')
PlotIntervals(run,'color','y');
hold all
ok = behavior.speed<2;
plot(behavior.timestamps(ok),behavior.speed(ok),'y');
plot(behavior.timestamps(ok),behavior.speed(ok),'r');
run = ConsolidateIntervals(run,'epsilon',1);
[in,w] = InIntervals(behavior.timestamps(:),run);
peak = Accumulate(w(in),behavior.speed(in)','mode','max');
hist(Restrict(peak,[0 10]),1000)
hist(Restrict(peak,[0 1]),1000)
PlotIntervals(run,'color','y');
PlotHVLines(0.2,'h','k--','linewidth',2);
PlotHVLines(0.2,'v','k--','linewidth',2);
PlotHVLines(0.2,'h','k--','linewidth',2);
PlotHVLines(0.1,'h','k--','linewidth',2);
k = kmeans(peaks,2);
open k
mmax(peaks(k==1))
mmax(peaks)
[group,threshold,em] = Otsu(peaks);
threshold
size(peaks)
[group,threshold,em] = Otsu(peak);
q = RemoveOutliers(peak);
mmax(q)
clf
sum(isnan(q))
peak(isnan(q))
mmax(peak)
min(peak(isnan(q))
min(peak(isnan(q)))
figure; PlotXY(behavior.positionTrials{1},'.')
hold all
PlotXY(behavior.positionTrials{2},'r.')
PlotIntervals(run(isnan(q),:),'color','y');
[group,threshold,em] = Otsu(peak);
peak = RemoveOutliers(peak);
[group,threshold,em] = Otsu(peak);
threshold
load(fullfile(basepath,[basename '.cell_metrics.cellinfo.mat']),'cell_metrics');
raly
load('day10.ripples.events.mat')
spikesCell = cell_metrics.spikes.times';
spikes = Group(spikesCell{:});
PETH(spikes(:,1),ripples.timestamps(:,1))
pos = Restrict(behavior.positionTrials{1},run);
PlotXY(pos);
PlotXY(pos,'.');
map = FiringMap(pos);
map = FiringCurve(pos,spikes(spikes(:,2)==1));
open map
plot(map.count);
plot(map.rate);
pos1 = Restrict(behavior.positionTrials{1},run);
pos2 = Restrict(behavior.positionTrials{2},run);
map2 = FiringCurve(pos2,spikes(spikes(:,2)==1));
hold all
plot(map.rate);
plot(map2.rate);
max(spikes(:,2))
for i=1:50
map = FiringCurve(pos1,spikes(spikes(:,2)==i));
curves1(i,:) = map.rate;
map = FiringCurve(pos2,spikes(spikes(:,2)==i));
curves2(i,:) = map.rate;
end
PlotColorMap([curves1 curves2]);
PlotColorMap(zscore([curves1 curves2],[],2))
PlotColorMap(sortby(zscore([curves1 curves2],[],2),out2(@max,curves1,[],2)))
clim
clim([-1 1]*2);
PlotColorMap(sortby(zscore([curves1 curves2],[],2),out2(@max,curves1,[],2)));  clim([-1 1]*2);
PlotColorMap(sortby(zscore([curves1 curves2],[],2),out2(@max,curves2,[],2)));  clim([-1 1]*2);
pos = sortrows([pos; pos2(:,1) 2-pos2(:,2)]);
figure; PlotXY(pos1);
hold all
PlotXY(pos2,'r');
PlotXY(pos,'k');
clf
PlotXY(pos1,'.');
hold all
PlotXY(pos2,'r.');
PlotXY(pos,'k.');
for i=1:50
map = FiringCurve(pos1,spikes(spikes(:,2)==i));
curves1(i,:) = map.rate;
map = FiringCurve(pos2,spikes(spikes(:,2)==i));
curves2(i,:) = map.rate;
map = FiringCurve(pos,spikes(spikes(:,2)==i));
curves(i,:) = map.rate;
end
PlotColorMap(sortby(zscore([curves1 curves2 curves],[],2),out2(@max,curves2,[],2)));  clim([-1 1]*2);
PlotColorMap(sortby(zscore([curves1 curves2 curves],[],2),out2(@max,curves,[],2)));  clim([-1 1]*2);
PlotColorMap(repmat(sortby(zscore([curves1 curves2 curves],[],2),out2(@max,curves,[],2)),2,1));  clim([-1 1]*2);
PlotColorMap(repmat(sortby(zscore([curves],[],2),out2(@max,curves,[],2)),2,1));  clim([-1 1]*2);
load('thetaCyclesTask.mat')
clear all
56
load('thetaCyclesTask.mat')
r = ripples.timestamps(:,1:2);
[windows] = SplitIntervals(cycles,'nPieces',6);
id = repmat((1:6)',length(cycles),1);
[rWindows,rID] = SplitIntervals(r,'pieceSize',0.02);
raly
regions = cell_metrics.brainRegion;
regionNames = unique(regions);
regionCell = zeros(length(regions),1);
for i=1:length(regionNames)
regionCell(strcmp(regions,regionNames{i})) = i;
end
regionNames
[estimations,actual,errors,average] = ReconstructPosition(pos,spikes,windows,'training',intervals,'id',id);
pos0 = sortrows([behavior.positionTrials{1}; behavior.positionTrials{2}(:,1) 2-behavior.positionTrials{2}(:,2)]);
pos = Restrict(pos,run);
pos0 = sortrows([behavior.positionTrials{1}; behavior.positionTrials{2}(:,1) 2-behavior.positionTrials{2}(:,2)]);
pos = Restrict(pos0,run);
45
load(fullfile(basepath,[basename '.ripple.events.mat']),'ripples');
load(fullfile(basepath,[basename '.ripples.events.mat']),'ripples');
45
PlotColorMap(average)
open pos
plot(pos(:,1))
sum(diff(pos(:,1))<0)
sum(diff(pos(:,1))==0)
sum(diff(pos(:,1))<=0)
[estimations,actual,errors,average] = ReconstructPosition(pos,spikes,windows,'training',intervals,'id',id);
PlotColorMap([average average]);
PlotHVLines(100,'h','w--','linewidth',2);
in = InIntervals(windows,stim);
[estimations,actual,errors,average] = ReconstructPosition(pos,spikes,windows(~in,:),'training',run,'id',id(~in));
PlotColorMap([average average]);
[estimationsS,actualS,errorsS,averageS] = ReconstructPosition(pos,spikes,windows(in,:),'training',run,'id',id(in));
PlotColorMap([average average]);
PlotColorMap([average averageS]);
PlotColorMap([averageS averageS]);
[estimations,actual,errors,average] = ReconstructPosition(pos,spikes,windows,'training',run,'id',id);
size(errors)
PlotColorMap(errors)
PlotColorMap(mean(reshape(errors,size(errors,1),[],6),3))
PlotColorMap(mean(reshape(errors,size(errors,1),[],6),2))
PlotColorMap(nanmean(reshape(errors,size(errors,1),[],6),2))
PlotColorMap(average)
PlotColorMap(nanmean(reshape(errors,size(errors,1),[],6),2))
PlotColorMap(nanmean(reshape(errors,size(errors,1),6,[]),3))
PlotColorMap(average)
PlotColorMap(nanmean(reshape(errors(:,in),size(errors,1),6,[]),3))
size(in)
sum(in)
in = repelem(InIntervals(cycles,stim),6,1);
size(in)
sum(in)
PlotColorMap(nanmean(reshape(errors(:,in),size(errors,1),6,[]),3))
PlotColorMap(nanmean(reshape(errors(:,~in),size(errors,1),6,[]),3))
PlotColorMap(repmat(nanmean(reshape(errors(:,~in),size(errors,1),6,[]),3),1,2));
clim
q = repmat(nanmean(reshape(errors(:,~in),size(errors,1),6,[]),3),1,2);
q = [q repmat(nanmean(reshape(errors(:,in),size(errors,1),6,[]),3),1,2)];
PlotColorMap(q);
hold off
PlotColorMap(estimations);
climc
clim
clim([0 0.2])
clim([0 0.1])
clim([0 0.05])
hold all
plot(actual);
open actual
mmax(actual)
clf
PlotColorMap(estimations);
hold all
plot(actual*100);
plot(actual*100,'w');
Portion(InIntervals(cycles,run))
45
PlotColorMap(estimations);
hold all
plot(actual*100,'w');
t = mean(windows,2);
tt = interp1(t,(1:length(t))',behavior.timestamps);
open tt
tt = interp1(t,(1:length(t))',behavior.timestamps(:));
hold all
plot(tt,behavior.speed*200,'r--');
t = behavior.timestamps(:);
speed = behavior.speed(:);
mmax(speed(InIntervals(t,run))
mmax(speed(InIntervals(t,run)))
threshold
run = behavior.timestamps(FindInterval(behavior.speed>0.02));
run = ConsolidateIntervals(run,'epsilon',0.01);
[in,w] = InIntervals(behavior.timestamps(:),run);
peak = Accumulate(w(in),behavior.speed(in)','mode','max');
% remove outliers (data in between sessions gives outlier speeds)
peak = RemoveOutliers(peak);
% remove run epochs that don't reach the speed threshold
[group,threshold,em] = Otsu(peak);
run(group==0,:) = [];
threshold
size(run)
dur(run)
mmax(speed(InIntervals(t,run)))
pos0 = sortrows([behavior.positionTrials{1}; behavior.positionTrials{2}(:,1) 2-behavior.positionTrials{2}(:,2)]);
pos = Restrict(pos0,run);
56
clf
PlotColorMap(estimations);
hold all
plot(tt,behavior.speed*200,'r--');
tt = interp1(t,(1:length(t))',behavior.timestamps(:));
clf
PlotColorMap(estimations);
hold all
plot(tt,behavior.speed*200,'r--');
mmax(speed(InIntervals(t,run)))
mmax(speed(InIntervals(t,cycles)))
bad = IntervalsIntersect(cycles, SubtractIntervals([0 Inf],run));
size(bad)
size(cycles)
Portion(bad)
load('thetaCyclesTask.mat','cycles');
bad = IntervalsIntersect(cycles, SubtractIntervals([0 Inf],run));
cycles(bad,:) = [];
Portion(bad)
tt = interp1(t,(1:length(t))',behavior.timestamps(:));
t = behavior.timestamps(:);
tt = interp1(t,(1:length(t))',behavior.timestamps(:));
clf
PlotColorMap(estimations);
hold all
plot(tt,behavior.speed*200,'r--');
mmax(speed(InIntervals(t,cycles)))
0.02*100
plot(t,speed,'r--');
ok = behavior.speed<2;
plot(t(ok),speed(ok),'r--');
ok = behavior.speed<1;
plot(t(ok),speed(ok),'r--');
PlotIntervals(cycles);
hold all
PlotIntervals(run,'color','y');
threshold
mmax(peak)
[group,threshold,em] = Otsu(peak);
mmax(threshold)
45
clf
PlotColorMap(estimations);
clf
q = repmat(nanmean(reshape(errors(:,~in),size(errors,1),6,[]),3),1,2);
in = repelem(InIntervals(cycles,stim),6,1);
q = repmat(nanmean(reshape(errors(:,~in),size(errors,1),6,[]),3),1,2);
q = [q repmat(nanmean(reshape(errors(:,in),size(errors,1),6,[]),3),1,2)];
PlotColorMap(q);
clf
PlotColorMap(estimations);
hold all
plot(actual*100,'w');
mmax(sum(estimations))
clim([0 0.05])
figure; PlotColorMap(q);
465
pos1 = Restrict(behavior.positionTrials{1},run);
pos2 = Restrict(behavior.positionTrials{2},run);
PlotColorMap(repmat(nanmean(reshape(errors(:,~in),size(errors,1),6,[]),3),1,2));
PlotColorMap(repmat(nanmean(reshape(errors(:,in),size(errors,1),6,[]),3),1,2));
q = repmat(nanmean(reshape(errors(:,~in),size(errors,1),6,[]),3),1,2);
q = [q repmat(nanmean(reshape(errors(:,in),size(errors,1),6,[]),3),1,2)];
PlotColorMap(repmat(nanmean(reshape(errors(:,in),size(errors,1),6,[]),3),1,2));
PlotColorMap(q);
pos2 = Restrict(behavior.positionTrials{2},run); pos2(:,2) = 1-pos2(:,2);
PlotColorMap(q);
clim
clim([0 0.05])
clim([0 0.02])
clim([0 0.01])
clim([0 0.02])
PlotColorMap(zscore(q));
clim
colormap(inferno)
clim
clim([-1 1]*2);
clim([0 1]*2);
clim([-1 2]);
clim([-1 3]);
clim([-1 4]);
pos2 = Restrict(behavior.positionTrials{2},run);
pos2 = Restrict(behavior.positionTrials{2},run); pos2(:,2) = 1-pos2(:,2);
sum(InIntervals(pos1,stim))
sum(InIntervals(pos2,stim))
sum(InIntervals(pos2,stim)) > sum(InIntervals(pos1,stim))
clims
[~,~,a,b] = FindReplayScore(off,'circular','off')
hold on; plot([1 6 [1 6]+6],[a b a b]);
set(gca,'ytick',100,'yticklabel','actual position');
[score,p,a,b,rShuffled,~,~,c,cShuffled] = FindReplayScore(average,'circular','off');
45
c
figure; hist(cShuffled,100);
[score,p,a,b,~,~,~,c,cShuffled] = FindReplayScore(average,'circular','off','wcorr','on');
c
figure; hist(cShuffled,100);
Portion(cShuffled<c)
p
score
PlotColorMap(average)
PlotColorMap(repmat(average,1,2))
PlotColorMap(repmat(circshift(average,1),1,2))
PlotColorMap(repmat(circshift(average',1)',1,2))
PlotColorMap(repmat(average,1,2))
PlotHVLines(100,'h','w--','linewidth',2);
intervals = run;
45
PlotColorMap(rEstimates);
PlotColorMap(rEstimations)
empty = (max(rEstimations)==min(rEstimations))';
[score,pValue,seqStartStop,stops] = deal(nan(size(r,1),1)); seqStartStop(:,2) = nan;
parfor i=1:length(r)
if sum(rID==i & ~empty)>=3
[score(i),pValue(i),seqStartStop(i),stops(i)] = FindReplayScore(rEstimations(:,rID==i),'circular','off');
end
if rem(i,100)==0, display([num2str(i) '/' num2str(size(r,1))]); end
end
seqStartStop(:,2) = stops;
slope = diff(seqStartStop,[],2)./Accumulate(rID,1,'size',[length(r) 1]); slope(seqStartStop(:,1)==0) = nan;
maze = mmax(run)
subplot(2,3,3);
pre = r(:,1)<maze(1); post = r(:,1)>maze(end);
i=length(r);
scatter(slope(1:i),score(1:i),10,post-pre,'filled');
PlotHVLines(0,'v','k--');
xlim([-1 1]*60);  ylim([0.2 1]);
subplot(2,2,3);
in = pre;
DensityMap(slope(in),score(in));
xlim([-1 1]*60); ylim([0 1]);
PlotHVLines(0,'v','k--');
ylabel('replay score'); xlabel('slope (bins/20ms)');
title(['pre sleep: ' num2str(100*sum(pValue(in)<0.05 & slope(in)>0)/sum(in)) '% forward and ' ...
num2str(100*sum(pValue(in)<0.05 & slope(in)<0)/sum(in)) '% backward']);
subplot(2,2,4);
in = post;
DensityMap(slope(in),score(in));
xlim([-1 1]*60); ylim([0 1]);
PlotHVLines(0,'v','k--');
ylabel('replay score'); xlabel('slope (bins/20ms)');
title(['post sleep: ' num2str(100*sum(pValue(in)<0.05 & slope(in)>0)/sum(in)) '% forward and ' ...
num2str(100*sum(pValue(in)<0.05 & slope(in)<0)/sum(in)) '% backward']);
drawnow
mmax(score)
cld
clf
plot(slope,score,'.');
hist(slope,200)
hist(slope,1000)
bad = slope==0;
plot(slope(~bad),score(~bad),'.');
hist(slope,1000)
kruskalbar(slope,p<0.05)
size(p)
sizse(pValue)
size(pValue)
size(slope)
kruskalbar(slope,pValue<0.05)
clf
kruskalbar(slope,post)
s = [slope slope]; s(post,1) = nan; s(pre,2) = nan;
kruskalbar(s,pValue<0.05))
kruskalbar(s,pValue<0.05)
anovabar(s,pValue<0.05)
sum(s(:)==0)
s(s==0) = nan;
anovabar(s,pValue<0.05)
kruskalbar(s,pValue<0.05)
plot(slope(~bad),score(~bad),'.');
plot(slope(~bad),-log(pValue(~bad)),'.');
bad = slope==0 | pValue<0.05;
plot(slope(~bad),-log(pValue(~bad)),'.');
plot(slope(~bad),score(~bad)),'.');
plot(slope(~bad),score(~bad),'.');
Portion(pValue<0.05)./Portion(~isnan(pValue)))
Portion(pValue<0.05)./Portion(~isnan(pValue))
anovabar(pValue<0.05,post)
size(sig)
sig = pValue; sig(~isnan(sig)) = pValue(~isnan(sig))<0.05;
size(sig)
open sig
anovabar(sig,post)
size(post)
sum(post)
sum(pre)
plot(-log(pValue))
plot(-log(pValue(~bad)))
bad = slope==0 | ~(pValue<0.05);
plot(-log(pValue(~bad)))
sum(isnan(pValue))
plot(-log(pValue(~bad)+0.0001))
plot(-log(pValue(~bad)+0.001))
plot(-log(pValue(~bad)+0.01))
plot(Smooth(-log(pValue(~bad)+0.01),2))
bad = slope<=0 | ~(pValue<0.05);
plot(Smooth(-log(pValue(~bad)+0.01),2))
bad = slope==0 | ~(pValue<0.05);
anovabar(slope,bad)
anovabar(s,bad)
anovabar(s)
anovabar(s(:,[1 2 2]))
anovabar(s(:,[1 2 2]),bad)
anovabar(s(:,[1 2 2]),post)
anovabar(s(:,[1 2 2]),bad)
open s
s0 = s; s0(s0<0) = -1; s0(s0>0) = 1;
open s0
anovabar(s(:,[1 2 2]),bad)
anovabar(s0(:,[1 2 2]),bad)
bad = slope==0 | ~(pValue<0.05) | pre;;
plot(slope(~bad),score(~bad)),'.');
size(bad)
size(score)
size(slope)
plot(slope(~bad),score(~bad)),'.');
plot(slope(~bad),score(~bad),'.');
FindClosest(slope,xlim)
FindClosest(slope,min(xlim))
min(xlim)
FindClosest(slope(~isnan(slope0),min(xlim))
FindClosest(slope(~isnan(slope),min(xlim)))
FindClosest(slope(~isnan(slope)),min(xlim))
FindClosest(slope,min(xlim))
sum(isnan(reference))
reference(isnan(reference)) = Inf;
% Make sure values to be matched are unique
[u,i] = unique(reference);
open u
% Find closest index
if strcmp(mode,'higher'),
indices = ceil(interp1(u,(1:length(u))',query));
elseif strcmp(mode,'lower'),
indices = floor(interp1(u,(1:length(u))',query));
else
indices = round(interp1(u,(1:length(u))',query));
end
nans = isnan(reference);
open nans
cumsum(~nans);
find(nans,1)
reference(reference==Inf) = nan;
find(nans,1)
nans = isnan(reference);
cumsum(~nans);
nonnans = find(~nans);
open nonnans
FindClosest(slope,min(xlim))
slope(1201)
xlim
score(1201)
FindClosest(score,min(ylim))
size(u)
FindClosest(score,min(ylim))
slope(6687)
score(6687)
r(6687)
post)6687)
post(6687)
mmax(rID)
open rID
PlotColorMap(rEstimations(:,rID==6687));
PlotColorMap(ReconstructPosition(pos,spikes,rWindows(rID==6687,:),'training',run);
PlotColorMap(ReconstructPosition(pos,spikes,rWindows(rID==6687,:),'training',run));
PlotColorMap(ReconstructPosition(pos1,spikes,rWindows(rID==6687,:),'training',run));
PlotColorMap(ReconstructPosition(pos2,spikes,rWindows(rID==6687,:),'training',run));
PlotColorMap(ReconstructPosition(pos,spikes,rWindows(rID==6687,:),'training',run));
PlotColorMap(ReconstructPosition(pos2,spikes,rWindows(rID==6687,:),'training',run));
a(6687)
seqStartStop(6687,:)
r(6687,:)
RasterPlot(spikes,r(6687,:));
RasterPlot(Restrict(spikes,r(6687,:)))
[~,m] = max(curves2,[],2);
for i=1:50
map = FiringCurve(pos1,spikes(spikes(:,2)==i));
curves1(i,:) = map.rate;
map = FiringCurve(pos2,spikes(spikes(:,2)==i));
curves2(i,:) = map.rate;
map = FiringCurve(pos,spikes(spikes(:,2)==i));
curves(i,:) = map.rate;
end
[~,m] = max(curves2,[],2);
open m
code(m) = 1:50;
plot(m)
[~,order] = sort(m);
open order
code(order) = 1:50;
code(1)
code(8)
ss = spikes; ss(:,2) = code(ss(:,2));
RasterPlot(Restrict(ss,r(6687,:)))
size(r)
PlotColorMap(ReconstructPosition(pos2,spikes,SplitIntervals(r(6687,:),0.01),'training',run));
PlotColorMap(ReconstructPosition(pos2,spikes,SplitIntervals(r(6687,:),'pieceSize',0.01)),'training',run));
PlotColorMap(ReconstructPosition(pos2,spikes,SplitIntervals(r(6687,:),'pieceSize',0.01),'training',run));
PlotColorMap(ReconstructPosition(pos2,spikes,SplitIntervals(r(6687,:),'pieceSize',0.005),'training',run));
PlotColorMap(ReconstructPosition(pos,spikes,SplitIntervals(r(6687,:),'pieceSize',0.005),'training',run));
PlotColorMap(ReconstructPosition(pos,spikes,SplitIntervals(r(6687,:),'pieceSize',0.05),'training',run));
PlotColorMap(ReconstructPosition(pos,spikes,SplitIntervals(r(6687,:),'pieceSize',0.02),'training',run));
PlotColorMap(ReconstructPosition(pos,spikes,SplitIntervals(r(6687,:),'pieceSize',0.015),'training',run));
PlotColorMap(ReconstructPosition(pos2,spikes,SplitIntervals(r(6687,:),'pieceSize',0.015),'training',run));
PlotColorMap(ReconstructPosition(pos2,spikes,SplitIntervals(r(6687,:),'pieceSize',0.02),'training',run));
SplitIntervals(r(6687,:),'pieceSize',0.02);
rWindows(rID==6687,:)
mmax(rID)
size(r)
open rID
open r
max(rID)
open r
diff(r(end,2))
diff(r(end,:))
plot(r(:,1))
hold all
plot(r(:,2))
open rWindwos
open rWindows
[rWindows,rID] = SplitIntervals(r,'pieceSize',0.02);
sum(pieces==0)
sum(pieces==1)
diff(pieces(end-5:end,:),[],2)
reset = pieces>0; reset(1)=1;
time2add = ones(size(pieces))*pieceSize;
time2add = CumSum(time2add,reset) - pieceSize; %rese
reset = pieces>0; reset(1)=1;
time2add = ones(size(pieces))*pieceSize;
time2add = CumSum(time2add,reset) - pieceSize; %reset at new bin starts and make sure baseline is 0
time2add = CumSum(time2add,reset) - pieceSize(:,1);
[rWindows,rID] = SplitIntervals(r,'pieceSize',0.02);
time2add = ones(size(pieces))*pieceSize;
time2add = CumSum(time2add,reset) - pieceSize;
if pieceSize<1, time2add = round(time2add/pieceSize)*pieceSize; end %fix underflow errers: time2add should be divisible by pieceSize
CumSum(pieces,reset) + time2add;
sum(reset)
min(d)
findmin(d)
min(piecesPerInterval)
sum(piecesPerInterval==0)
time2add = ones(size(pieces))*pieceSize;
time2add = CumSum(time2add,reset) - pieceSize; %reset at new bin starts and make sure baseline is 0
if pieceSize<1, time2add = round(time2add/pieceSize)*pieceSize; end %fix underflow errers: time2add should be divisible by pieceSize
pieces = CumSum(pieces,reset) + time2add; %repeat first pieces for the duration of the interval
pieces(:,2) = pieces(:,1)+pieceSize;
ids = CumSum(reset);
indicesNotEmpty = find(piecesPerInterval>0);
ids = indicesNotEmpty(ids);
mmax(ids)
[rWindows,rID] = SplitIntervals(r,'pieceSize',0.02);
[rWindows,rID0] = SplitIntervals(r,'pieceSize',0.02);
[rWindows,rID] = SplitIntervals(r,'pieceSize',0.02);
open rID0
open rID
mmax(rID(rID0==6687))
PlotColorMap(ReconstructPosition(pos2,spikes,SplitIntervals(r(6696,:),'pieceSize',0.02),'training',run));
PlotColorMap(ReconstructPosition(pos,spikes,SplitIntervals(r(6696,:),'pieceSize',0.02),'training',run));
PlotColorMap(ReconstructPosition(pos2,spikes,SplitIntervals(r(6696,:),'pieceSize',0.02),'training',run));
PlotColorMap(ReconstructPosition(pos,spikes,SplitIntervals(r(6696,:),'pieceSize',0.02),'training',run));
subplot(1,2,1);
PlotColorMap(ReconstructPosition(pos2,spikes,SplitIntervals(r(6696,:),'pieceSize',0.02),'training',run));
subplot(1,2,2);
PlotColorMap(ReconstructPosition(pos,spikes,SplitIntervals(r(6696,:),'pieceSize',0.02),'training',run));
subplot(1,3,1);
PlotColorMap(ReconstructPosition(pos1,spikes,SplitIntervals(r(6696,:),'pieceSize',0.02),'training',run));
subplot(1,3,2);
PlotColorMap(ReconstructPosition(pos2,spikes,SplitIntervals(r(6696,:),'pieceSize',0.02),'training',run));
subplot(1,3,3);
PlotColorMap(ReconstructPosition(pos,spikes,SplitIntervals(r(6696,:),'pieceSize',0.02),'training',run));
clim
clims
clims([0 0.2]);
plot(mean(intervals,2),count);
CountInIntervals(spikes(:,1),intervals(1:1000,:))
sum(CountInIntervals(spikes(:,1),intervals(1:1000,:)))
intervals(1000,:)
diff(intervals(1000,:))
sum(InIntervals(spikes(:,1),diff(intervals(1000,:))))
sum(InIntervals(spikes(:,1),intervals(1000,:)))
mmax(diff(intervals(:,1)))
hm  = bsxfun(@plus,mean(intervals,2),[-1 1]*0.005/2);
spikes = Group(spikesCell{:});
hm = Restrict(r,mmax(intervals));
open hm
int = interp1(intervals(:,1),(1:length(intervals(:,1)))',hm);
fr = Accumulate(spikes(:,2))./spikes(end,1)
open fr
figure; hist(fr,100)
bad = fr>2.2;
CellExplorer
pyr = cellfun(@(x) ~isempty(strfind(x,'Pyramidal')), cell_metrics.putativeCellType);
open pyr
pyr = cellfun(@(x) ~isempty(strfind(x,'Pyramidal')), cell_metrics.putativeCellType)';
sum(bad)
sum(bad(pyr))
sum(bad(~pyr))
sum(~pyr)
mmax(fr(pyr))
mmax(fr(~pyr))
pyr = cellfun(@(x) ~isempty(strfind(x,'Pyramidal')), cell_metrics.putativeCellType)';
spikes = Group(spikesCell{pyr});
mmax(spikes(:,2))
bursts = FindBursts(spikes);
open bursts
PlotIntervals(int,'color','r');
mmax(diff(bursts(:,[1 3]),[],2))
bursts = FindBursts2(spikes);
mmax(diff(bursts(:,[1 3]),[],2))
subplot(2,1,1);
RasterPlot(Restrict(ss,mmax(intervals)));
[~,m] = max(curves2,[],2);
clear curves curves1 curves2
for i=1:max(spikes(:,2))
map = FiringCurve(pos1,spikes(spikes(:,2)==i));
curves1(i,:) = map.rate;
map = FiringCurve(pos2,spikes(spikes(:,2)==i));
curves2(i,:) = map.rate;
map = FiringCurve(pos,spikes(spikes(:,2)==i));
curves(i,:) = map.rate;
end
[~,m] = max(curves2,[],2);
[~,order] = sort(m);
clead code
clear code
code(order) = 1:50;
code(order) = 1:max(spikes(:,2))
code(order) = 1:max(spikes(:,2));
ss = spikes; ss(:,2) = code(ss(:,2));
RasterPlot(Restrict(ss,mmax(intervals)));
PlotColorMap(ReconstructPosition(pos2,spikes,intervals,'training',SubtractIntervals(run,stim)));
PlotHVLines(int,'v','w--','linewidth',2);
45
PlotColorMap(zscore([curves1 curves2],[],2))
PlotColorMap(zscore([curves1 curves2 curve],[],2))
PlotColorMap(zscore([curves1 curves2 curves],[],2))
PlotColorMap(sortby(zscore([curves1 curves2 curves],[],2),m))
[~,m] = max(Smooth(curves1,[0 2],'type','cc),[],2);
[~,m] = max(Smooth(curves1,[0 2],'type','cc'),[],2);
PlotColorMap(sortby(zscore([curves1 curves2 curves],[],2),m))
PlotColorMap(sortby(nanzscore([curves1 curves2 curves],[],2),m))
order(3)
[~,order] = sort(m);
order(3)
sum(spikes(:,2)==31)
figure; hist(spikes(spikes(:,2)==31))
hist(spikes(spikes(:,2)==31),1000)
ok = InIntervaps(spikes(:,1),mmax(run));
ok = InIntervals(spikes(:,1),mmax(run));
nSpikes = Accumulate(spikes(ok,2));
open nSpikes
spikes = Group(spikesCell{pyr & nSpikes>200});
ok = InIntervals(spikes(:,1),mmax(run));
nSpikes = Accumulate(spikes(ok,2));
pyr = cellfun(@(x) ~isempty(strfind(x,'Pyramidal')), cell_metrics.putativeCellType)';
spikes = Group(spikesCell{pyr & nSpikes>200});
nSpikes = cellfun(@(x) sum(x(:,1)>run(1) & x(:,1)<run(end)), spikesCell);
% load spikes
spikesCell = cell_metrics.spikes.times';
nSpikes = cellfun(@(x) sum(x(:,1)>run(1) & x(:,1)<run(end)), spikesCell);
pyr = cellfun(@(x) ~isempty(strfind(x,'Pyramidal')), cell_metrics.putativeCellType)';
spikes = Group(spikesCell{pyr & nSpikes>200});
max(spikes(:,2))
sum(pyr)
clear curves curves1 curves2
for i=1:max(spikes(:,2))
map = FiringCurve(Restrict(pos1,SubtractIntervals(run,SubtractIntervals([0 Inf],stim)),'shift','on'),Restrict(spikes(spikes(:,2)==i),SubtractIntervals(run,SubtractIntervals([0 Inf],stim)),'shift','on'));
curves1(i,:) = map.rate;
map = FiringCurve(Restrict(pos2,SubtractIntervals(run,stim),'shift','on'),Restrict(spikes(spikes(:,2)==i),SubtractIntervals(run,stim),'shift','on'));
curves2(i,:) = map.rate;
map = FiringCurve(Restrict(pos,run,'shift','on'),Restrict(spikes(spikes(:,2)==i),run,'shift','on'));
curves(i,:) = map.rate;
end
[~,m] = max(Smooth(curves1,[0 2],'type','cc'),[],2);
[~,order] = sort(m);
PlotColorMap(sortby(nanzscore([curves1 curves2 curves],[],2),m))
PlotColorMap(sortby(zBaseline([curves1 curves2 curves],51:150),m))
PlotColorMap(sortby(zBaseline([curves1 curves2 curves],51:150,2),m))
clim
clim([-1 1]*2);
PlotColorMap(sortby(nanzscore([curves1 curves2 curves],[],2),m))
clim([-1 1]*2);
[~,m] = max(Smooth(curves2,[0 2],'type','cc'),[],2);
PlotColorMap(sortby(nanzscore([curves1 curves2 curves],[],2),m))
[~,m] = max(Smooth(curves,[0 2],'type','cc'),[],2);
PlotColorMap(sortby(nanzscore([curves1 curves2 curves],[],2),m))
curves1 curves2
PlotColorMap(sortby(nanzscore([curves],[],2),m))
PlotHVLines(25.5,'v');
bursts = FindBursts2(spikes,'thresholds',[-1 5]);
bursts = FindBursts2(spikes,'thresholds',[1 5]);
size(bursts)
bursts = FindBursts2(spikes,'thresholds',[0.5 5]);
size(bursts)
(0.05/3)
open hm
[19251.2225000000,19251.3075000000]
bursts = FindBursts(spikes,'thresholds',[0.5 5]);
figure; plot(t,z);
size(t)
size(z)
[binned,t] = binspikes(spikes(:,1),1/binSize); t = t(:);
figure; plot(t,z);
xlim([19251.2225000000,19251.3075000000]+[-1 1]);
PlotIntervals(Restrict(bursts(:,[1 3]),xlim))
if thresholds(1)>-Inf
% Find periods depassing the low threshold
bursts(:,[1 3]) = t(ToIntervals(z>thresholds(1)));
else
bursts(:,[1 3]) = t(ToIntervals(binned>0));
end
PlotIntervals(Restrict(bursts(:,[1 3]),xlim),'color','r')
bursts(:,[1 3]) = t(ToIntervals(binned>0));
clear bursts
bursts(:,[1 3]) = t(ToIntervals(binned>0));
PlotIntervals(Restrict(bursts(:,[1 3]),xlim),'color','y')
binned = Smooth(binned,smooth);
clear bursts
bursts(:,[1 3]) = t(ToIntervals(binned>0));
figure; plot(t,binned);
xlim([19251.2225000000,19251.3075000000]+[-1 1]);
PlotHVLines(eps,'h','k--')
sum(smoothed>0)
sum(binned<=0)
sum(binned<=eps)
bursts(:,[1 3]) = t(ToIntervals(binned>eps*2));
PlotIntervals(Restrict(bursts(:,[1 3]),xlim),'color','y')
clear bursts
bursts(:,[1 3]) = t(ToIntervals(binned>eps*10));
PlotIntervals(Restrict(bursts(:,[1 3]),xlim),'color','g')
thresholds
thresholds = [-0.5 5];
clear burst
bursts(:,[1 3]) = t(ToIntervals(z>thresholds(1)));
clear bursts
bursts(:,[1 3]) = t(ToIntervals(z>thresholds(1)));
figure; plot(t,z);
xlim([19251.2225000000,19251.3075000000]+[-1 1]);
PlotIntervals(Restrict(bursts(:,[1 3]),xlim),'color','g')
smooth
smooth = 20;
[binned,t] = binspikes(spikes(:,1),1/binSize); t = t(:);
binned = Smooth(binned,smooth);
z = zscore(binned);
plot(t,z);
smooth = 50;
[binned,t] = binspikes(spikes(:,1),1/binSize); t = t(:);
binned = Smooth(binned,smooth);
z = zscore(binned);
plot(t,z);
PlotHVLines(0,'h','k--')
xlim([19251.2225000000,19251.3075000000]+[-1 1]*10);
plot(t,z,'linewidth'm2);
plot(t,z,'linewidth',2);
plot(t,z,'k','linewidth',2);
PlotHVLines(0,'h','k--')
smooth = 100;
[binned,t] = binspikes(spikes(:,1),1/binSize); t = t(:);
binned = Smooth(binned,smooth);
z = zscore(binned);
plot(t,z,'y','linewidth',2);
bursts = FindBursts(spikes,'thresholds',[0 2],'smooth',0.05);
isdvector([0 2])
isdvector([0 2],'#2')
isdvector([0 2],'#2','>')
isdvector([0 2],'#2','<')
bursts = FindBursts(spikes,'thresholds',[0 2],'smooth',0.05);
dbcont
bursts = FindBursts(spikes,'thresholds',[0 2],'smooth',0.05);
binSize
smooth
smooth = 0.05;
binSize
smooth./binSize
smooth = smooth/binSize;
spikes = spikes(:,1);
[binned,t] = binspikes(spikes(:,1),1/binSize); t = t(:);
binned = Smooth(binned,smooth);
z = zscore(binned);
% Find periods depassing the low threshold
bursts(:,[1 3]) = t(ToIntervals(z>thresholds(1)));
thresholds
figure; plot(t,z);
xlim([19251.2225000000,19251.3075000000]+[-1 1]*10);
PlotIntervals(Restrict(bursts(:,[1 3]),xlim),'color','g')
% Find maximal z value within putative bursts
[in,id] = InIntervals(t,bursts(:,[1 3]));
[peak,idx] = Accumulate(id(in),z(in),'mode','max');
t = t(in);
bursts(:,2) = t(idx);
in = InIntervals(bursts(:,1),xlim);
mmax(peak(in))
thresholds
pass = peak>thresholds(2);
sum(pass)
bursts = FindBursts(spikes,'thresholds',[0 2],'smooth',0.05);
intervals = Bins(19251.2225000000,19251.3075000000,0.02);
open intervals
45
34
open intervals
q = ReconstructPosition(pos2,spikes,intervals,'training',SubtractIntervals(run,stim))
q = ReconstructPosition(pos2,spikes,intervals,'training',SubtractIntervals(run,stim));
open q
[score0,pValue0,seqStartStop0,stops0] = FindReplayScore(rEstimations(:,rID==i),'circular','off');
score0
pValue0
seqStartStop0
stops0
size(r)
r = bursts(:,[1 3]);
[rWindows,rID] = SplitIntervals(bursts(:,[1 3]),'pieceSize',0.02);
[rEstimations] = ReconstructPosition(pos2,spikes,rWindows,'training',intervals);
45.
clf
PlotColorMap(rEstimations)l;
PlotColorMap(rEstimations);
clim
clim([0 0.3])
clim([0 0.2])
PlotHVLines(cumsum(Accumulate(rID))+0.5,'w--');
rID(40263)
[binned,t] = binspikes(spikes(:,1),1/0.001); t = t(:);
binned = Smooth(binned,50);
z = zscore(binned);
figure; plot(t,z);
xlim([-1 1]+bursts(2667,[1 3]))
PlotIntervals(Restrict(bursts(:,[1 3]),xlim),'color','g')
diff(bursts(2667,[1 3]))
[binned,t] = binspikes(spikes(:,1),1/0.001); t = t(:);
binned = Smooth(binned,20);
binned = Smooth(binned,30);
plot(t,z,'y','linewidth',2);
[binned,t] = binspikes(spikes(:,1),1/0.001); t = t(:);
binned = Smooth(binned,30);
plot(t,z,'y','linewidth',2);
[binned,t] = binspikes(spikes(:,1),1/0.001); t = t(:);
binned = Smooth(binned,10);
plot(t,z,'y','linewidth',2);
plot(t,z,'r','linewidth',2);
[binned,t] = binspikes(spikes(:,1),1/0.001); t = t(:);
clf
[binned,t] = binspikes(spikes(:,1),1/0.001); t = t(:);
binned = Smooth(binned,30); z = zscore(binned);
plot(t,z,'y','linewidth',2);
xlim([-1 1]+bursts(2667,[1 3]))
[binned,t] = binspikes(spikes(:,1),1/0.001); t = t(:);
binned = Smooth(binned,10); z = zscore(binned);
plot(t,z,'r','linewidth',2);
xlim([-1 1]+bursts(2667,[1 3]))
bursts = FindBursts(spikes,'thresholds',[2 3],'smooth',0.03);
FindClosest(bursts(:,2),16153)
PlotIntervals(Restrict(bursts(:,[1 3]),xlim),'color','g')
bursts = FindBursts(spikes,'thresholds',[2 3],'smooth',0.01);
PlotIntervals(Restrict(bursts(:,[1 3]),xlim),'color','g')
FindClosest(bursts(:,2),16153)
burst(6747,2)-16153
bursts(6747,2)-16153
[rWindows,rID] = SplitIntervals(bursts(:,[1 3]),'pieceSize',0.02);
[rEstimations] = ReconstructPosition(pos2,spikes,rWindows,'training',intervals);
clf
PlotColorMap(rEstimations);
find(rID==6747)
PlotHVLines(cumsum(Accumulate(rID))+0.5,'w--');
figure; plot(t,z);
xlim([-1 1]+bursts(2667,[1 3]))
xlim([-1 1]+bursts(6747,[1 3]))
PlotIntervals(Restrict(bursts(:,[1 3]),xlim),'color','g')
PlotHVLines(bursts(6747,2),'v','k--');
xlim([-1 1]*5+bursts(6747,[1 3]))
PlotIntervals(Restrict(bursts(:,[1 3]),xlim),'color','g')
intervals = bursts(:,[1 3]); intervals(:,2) = intervals(:,2)+0.01;
sum(intervals(1:end,2)>intervals(2:end,1))
sum(intervals(1:end-1,2)>intervals(2:end,1))
intervals = bursts(:,[1 3]); intervals(1:end-1,2) = min([intervals(1:end-1,2)+0.01 mean([intervals(1:end-1,2) intervals(2:end,1)],2)]); intervals(end,2) = intervals(end,2) + 0.01;
intervals = bursts(:,[1 3]); intervals(1:end-1,2) = min([intervals(1:end-1,2)+0.01 mean([intervals(1:end-1,2) intervals(2:end,1)],2)],[],2); intervals(end,2) = intervals(end,2) + 0.01;
intervals(1:end-1,2) = min([intervals(1:end-1,2)+0.01 mean([intervals(1:end-1,2) intervals(2:end,1)],2)],[],2); intervals(end,2) = intervals(end,2) + 0.01;
intervals(2:end,1) = max([intervals(2:end,1)-0.01 mean([intervals(2:end,1) intervals(1:end-1,2)],2)],[],2); intervals(1,1) = intervals(1,1) - 0.01;
intervals = bursts(:,[1 3]);
intervals0 = intervals;
intervals(1:end-1,2) = min([intervals(1:end-1,2)+0.01 mean([intervals(1:end-1,2) intervals(2:end,1)],2)],[],2); intervals(end,2) = intervals(end,2) + 0.01;
intervals(2:end,1) = max([intervals(2:end,1)-0.01 mean([intervals(2:end,1) intervals(1:end-1,2)],2)],[],2); intervals(1,1) = intervals(1,1) - 0.01;
d = intervals - intervals0
open d
intervals = bursts(:,[1 3]); intervals0 = intervals;
intervals(1:end-1,2) = min([intervals0(1:end-1,2)+0.01 mean([intervals0(1:end-1,2) intervals0(2:end,1)],2)],[],2); intervals(end,2) = intervals0(end,2) + 0.01;
intervals(2:end,1) = max([intervals0(2:end,1)-0.01 mean([intervals0(2:end,1) intervals0(1:end-1,2)],2)],[],2); intervals(1,1) = intervals0(1,1) - 0.01;
d = intervals - intervals0;
56
clf
PlotColorMap(rEstimations);
PlotIntervals(Restrict(bursts(:,[1 3]),xlim),'color','g')
FindClosest(bursts(:,2),16153)
clf
PlotColorMap(rEstimations);
PlotHVLines(bursts(6747,2),'v','k--');
PlotHVLines(cumsum(Accumulate(rID))+0.5,'w--');
find(rID==2667)
find(rID==2667);
find(rID==2667,1)
clim([0 0.2])
[binned,t] = binspikes(spikes(:,1),1/0.001); t = t(:);
binned = Smooth(binned,10); z = zscore(binned);
figure; plot(t,z);
xlim([-1 1]*5+bursts(6747,[1 3]))
xlim([-1 1]*5+bursts(2667,[1 3]))
PlotIntervals(Restrict(bursts(:,[1 3]),xlim),'color','g')
PlotIntervals(Restrict(bursts(:,[1 3]),xlim),'color','g');
diff(Restrict(bursts(:,[1 3]),xlim),[],2)
rID(43139)
xlim([-1 1]*5+bursts(2681,[1 3]))
PlotIntervals(Restrict(bursts(:,[1 3]),xlim),'color','g');
bursts = FindBursts(spikes,'thresholds',[2 3],'smooth',0.01);
clf
PlotColorMap(rEstimations);
PlotHVLines(cumsum(Accumulate(rID))+0.5,'w--');
FindClosest(bursts(:,2),16153)
rID(6747)
find(rID==6747,1)
xlim(19156 + [-1 1]*20);
clim([0 0.2])
FindClosest(bursts(:,2),16204)
find(rID==6780,1)
xlim(19253 + [-1 1]*20);
PlotIntervals(Restrict(bursts(:,[1 3]),xlim),'color','g');
xlim([-1 1]*5+bursts(6780,[1 3]))
xlim(19253 + [-1 1]*20);
PlotHVLines(bursts(6780,[1 3],'k');
PlotHVLines(bursts(6780,[1 3],'k'));
PlotHVLines(bursts(6780,[1 3]),'k');
PlotHVLines(bursts(6781,[1 3]),'r');
PlotIntervals(Restrict(bursts(:,[1 3]),xlim),'color','g');
PlotIntervals(Restrict(bursts(:,[1 3]),xlim),'color','r');
intervals = bursts(:,[1 3]); intervals0 = intervals; leeway = 0.01;
intervals(1:end-1,2) = min([intervals0(1:end-1,2)+leeway mean([intervals0(1:end-1,2) intervals0(2:end,1)],2)],[],2); intervals(end,2) = intervals0(end,2) + leeway;
intervals(2:end,1) = max([intervals0(2:end,1)-leeway mean([intervals0(2:end,1) intervals0(1:end-1,2)],2)],[],2); intervals(1,1) = intervals0(1,1) - leeway;
PlotHVLines(Restrict(intervals,xlim),'color','r');
PlotHVLines(Restrict(intervals,xlim),'color','y');
PlotHVLines(Restrict(intervals,xlim),'color','k','linewidth','--');
PlotHVLines(Restrict(intervals,xlim),'color','k--','linewidth','2');
PlotHVLines(Restrict(intervals,xlim),'v','k--','linewidth','2');
PlotHVLines(Restrict(intervals,xlim),'v','k--','linewidth',2);
xlim
rIntervals = Bins(min(xlim),max(xlim),0.02);
intervals = SubtractIntervals(run,stim);
figure; PlotColorMap(ReconstructPosition(pos2,spikes,rWindows,'training',intervals),'x',mean(rIntervals,2));
clim([0 0.2])
PlotHVLines(Restrict(bursts(:,[1 3]),xlim),'v','w--','linewidth',2);
figure; PlotColorMap(ReconstructPosition(pos2,spikes,rIntervals,'training',intervals),'x',mean(rIntervals,2));
PlotHVLines(Restrict(bursts(:,[1 3]),xlim),'v','w--','linewidth',2);
x = xlim
xlim(x)
[~,m] = max(Smooth(curves2,[0 2],'type','cc'),[],2);
[~,order] = sort(m);
clear code; code(order) = max(spikes(:,2));
ss = spikes; ss(:,2) = code(ss(:,2));
figure; RasterPlot(Restrict(ss,x));
clear code; code(order) = 1:max(spikes(:,2));
ss = spikes; ss(:,2) = code(ss(:,2));
figure; RasterPlot(Restrict(ss,x));
PlotHVLines(Restrict(bursts(:,[1 3]),xlim),'v','k--','linewidth',2);
burstsCeline = FindBursts(spikes,'thresholds',[0 3],'smooth',0.01);
PlotHVLines(Restrict(burstsCeline(:,[1 3]),xlim),'v','r--','linewidth',2);
diff(Restrict(burstsCeline(:,[1 3]),xlim),[],2)
PlotHVLines(Restrict(burstsCeline(:,[1 3]),xlim),'v','r--','linewidth',2);
bursts = FindBursts(spikes,'thresholds',[0 3],'smooth',0.01);
bursts = FindBursts(spikes,'thresholds',[0 3],'smooth',0.01);
intervals = bursts(:,[1 3]); intervals0 = intervals; leeway = 0.01;
intervals(1:end-1,2) = min([intervals0(1:end-1,2)+leeway mean([intervals0(1:end-1,2) intervals0(2:end,1)],2)],[],2); intervals(end,2) = intervals0(end,2) + leeway;
intervals(2:end,1) = max([intervals0(2:end,1)-leeway mean([intervals0(2:end,1) intervals0(1:end-1,2)],2)],[],2); intervals(1,1) = intervals0(1,1) - leeway;
[rWindows,rID] = SplitIntervals(intervals,'pieceSize',0.02);
dur(rWindows)
[rEstimations] = ReconstructPosition(pos2,spikes,rWindows,'training',intervals);
[parentFolder,dayName] = fileparts(folder);
[~,projectName] = fileparts(parentFolder);
sessionID = [projectName '_' dayName];
folder = basepath
[parentFolder,dayName] = fileparts(folder);
[~,projectName] = fileparts(parentFolder);
sessionID = [projectName '_' dayName];
name = [sessionID '_OML_off.mat'];
sessionID
empty = (max(rEstimations)==min(rEstimations))';
[score,pValue,seqStartStop,stops] = deal(nan(size(r,1),1)); seqStartStop(:,2) = nan;
parfor i=1:length(r)
if sum(rID==i & ~empty)>=3
[score(i),pValue(i),seqStartStop(i),stops(i)] = FindReplayScore(rEstimations(:,rID==i),'circular','off');
end
if rem(i,100)==0, display([num2str(i) '/' num2str(size(r,1))]); end
end
seqStartStop(:,2) = stops;
save(['C:\Users\Cornell\Documents\code\new\tempFiles\' name],'estimations','actual','errors','average','rEstimations','score','pValue','seqStartStop');
r = bursts(:,[1 3]);
empty = (max(rEstimations)==min(rEstimations))';
[score,pValue,seqStartStop,stops] = deal(nan(size(r,1),1)); seqStartStop(:,2) = nan;
parfor i=1:length(r)
if sum(rID==i & ~empty)>=3
[score(i),pValue(i),seqStartStop(i),stops(i)] = FindReplayScore(rEstimations(:,rID==i),'circular','off');
end
if rem(i,100)==0, display([num2str(i) '/' num2str(size(r,1))]); end
end
seqStartStop(:,2) = stops;
save(['C:\Users\Cornell\Documents\code\new\tempFiles\' name],'estimations','actual','errors','average','rEstimations','score','pValue','seqStartStop');
% 3000 = 1900
raly
save(['M:\home\raly\Documents\code\new\tempFiles\' name],'estimations','actual','errors','average','rEstimations','score','pValue','seqStartStop');
save(['M:\home\raly\Documents\code\new\tempFiles\' name],'estimations','actual','errors','average','rEstimations','score','pValue','seqStartStop','bursts');
figure; plot(slope,score,'.');
slope = diff(seqStartStop,[],2)./Accumulate(rID,1,'size',[length(r) 1]); slope(seqStartStop(:,1)==0) = nan;
plot(t,z,'r','linewidth',2);
c;f
clf
plot(slope,score,'.');
anovabar(slope>0,post);
s = [slope slope]; s(post,1) = nan; s(pre,2) = nan;
anovabar(s,pValue<0)
anovabar(s,pValue<0.05)
s0 = s; s0(s0<0) = -1; s0(s0>0) = 1;
anovabar(s0,pValue<0.05)
anovabar(s0(:,[1 2 2]),pValue<0.05)
kruskalbar(s0(:,[1 2 2]),pValue<0.05)
sum(s==0)
sum(s>0 & pValue<0.05)
sum(s<0 & pValue<0.05)
zBinomialComparison(847,799+847,0.5)
zBinomialComparison(1109,920+1109,0.5)
zBinomialComparison(1109,920+1109,847,799+847)
sum(pValue<0.05)
sum(~isnan(s) & repmat(pValue<0.05,1,2))
zBinomialComparison(1109,2089,847,1687)
size(post)
post = InIntervals(bursts,postSleep);
sleep = ConsolidateIntervals(MergePoints.timestamps(cellfun(@(x) any(strfind(lower(x),'sleep')),MergePoints.foldernames),:))
load(fullfile(basepath,[basename '.MergePoints.events.mat']),'MergePoints');
sleep = ConsolidateIntervals(MergePoints.timestamps(cellfun(@(x) any(strfind(lower(x),'sleep')),MergePoints.foldernames),:));
load(fullfile(basepath,[basename '.SleepState.states.mat']));
sws = SleepState.ints.NREMstate;
open sleep
PlotXY(pos1)
PlotIntervals(stim);
PlotIntervals(stim,'color','r'));
PlotIntervals(stim,'color','r');
preSleep = SubtractIntervals(sleep(1,:), SubtractIntervals([0 Inf],sws));
postSleep = SubtractIntervals(sleep(2:end,:), SubtractIntervals([0 Inf],sws));
pre = InIntervals(bursts,preSleep);
post = InIntervals(bursts,postSleep);
s = [slope slope]; s(post,1) = nan; s(pre,2) = nan;
figure;
s0 = s; s0(s0<0) = -1; s0(s0>0) = 1;
anovabar(s0(:,[1 2 2]),pValue<0.05)
open s
s = [slope slope]; s(post,1) = nan; s(pre,2) = nan;
open pre
s = [slope slope]; s(~post,2) = nan; s(~pre,1) = nan;
sum(~isnan(s))
s0 = s; s0(s0<0) = -1; s0(s0>0) = 1;
anovabar(s0(:,[1 2 2]),pValue<0.05)
sum(~isnan(s) & repmat(pValue<0.05,1,2))
zBinomialComparison(1109,2089,847,1687)
sum(~isnan(s) & repmat(pValue<0.05,1,2))
sum(x>0 & repmat(pValue<0.05,1,2))
sum(~isnan(s) & repmat(pValue<0.05,1,2))
sum(x>0 & repmat(pValue<0.05,1,2))
sum(s>0 & repmat(pValue<0.05,1,2))
sum(~isnan(s) & repmat(pValue<0.05,1,2))
sum(s>0 & repmat(pValue<0.05,1,2))
sum(s<0 & repmat(pValue<0.05,1,2))
zBinomialComparison(515,991,0.5)
size(pValue)
zBinomialComparison(515,450+515,0.5)
zBinomialComparison(83,83+72,0.5)
(83/(83+72))
(515/(515+450))
open pValue
sum(empty)
size(empty)
length(r)
size(r)
size(bursts)
45
tic
for i=1:100
rEstimationsCell{i,1} = rEstimations(:,rID==i);
end
toc
45
outputs(i,:) = FindReplayScore(rEstimationsCell{i},'circular','off');
clear outputs
i=1
outputs{i,:} = FindReplayScore(rEstimationsCell{i},'circular','off');
edit temptemp
raly
edit temp
clear outputs
outputs(i,:) = temp(rEstimationsCell{i},'circular','off');
45
size(r,1)
2
shuffledScores = cell2mat(outputs(:,5));
shuffledSlopes = cell2mat(outputs(:,7))-cell2mat(outputs(:,6));
shuffleID = ceil((1:size(shuffledScores,1))'/500);
size(shuffledScores,1)/500
size(scores,1)
size(score,1)
scores = cellmat(outputs(:,1));
scores = cell2mat(outputs(:,1));
bad = cellfun(@isempty,outputs(:,1));
size(bad)
size(scores)
sum(bad)
914+11482
shuffleID = repmat((1:500)',size(scores,1));
shuffleID = repmat((1:500)',size(scores,1),1);
size(outputs{find(~bad,1),5})
size(outputs{find(~bad,1),5},1)
mmaxShuffledSlopes
mmax(shuffledSlopes)
mmax(scores)
mmax(shuffledScores)
size(rEstimations,1)
normalizedShuffledSlopes = (shuffledSlopes+size(rEstimations,1))/(size(rEstimations,1)*2);
mmax(normalizedShuffledSlopes)
size(slopes_)
size(slope)
slopes = cell2mat(outputs(:,4))-cell2mat(outputs(:,3));
normalizedSlopes = (slopes+size(rEstimations,1))/(size(rEstimations,1)*2);
map = DensityMap(normalizedSlopes,scores,'nBins',nBins);
nBins = [50 50];
map = DensityMap(normalizedSlopes,scores,'nBins',nBins);
PlotColorMap(map);
map = DensityMap(normalizedSlopes,scores,'nBins',nBins,'smooth',0);
546
zmap = (map-nanmean(mapsShuffled,3))./nanstd(mapsShuffled,[],3));
zmap = (map-nanmean(mapsShuffled,3))./nanstd(mapsShuffled,[],3);
PlotColorMap(zmap);
clim
clim([-1 1]*5);
PlotColorMap(zmap,'x',linspace(-size(rEstimations,1),size(rEstimations,1),nBins(1)),'y',linspace(0,1,nBins(2)));
PlotHVLines(0,'v','w--','linewidth',2);
45
zmap = (map-nanmean(mapsShuffled,3))./nanstd(mapsShuffled,[],3);
PlotColorMap(zmap,'x',linspace(-size(rEstimations,1),size(rEstimations,1),nBins(1)),'y',linspace(0,1,nBins(2)));
PlotHVLines(0,'v','w--','linewidth',2);
map = DensityMap(normalizedSlopes,scores,'nBins',nBins([2 1]),'smooth',0);
zmap = (map-nanmean(mapsShuffled,3))./nanstd(mapsShuffled,[],3);
PlotColorMap(zmap,'x',linspace(-size(rEstimations,1),size(rEstimations,1),nBins(1)),'y',linspace(0,1,nBins(2)));
PlotHVLines(0,'v','w--','linewidth',2);
size(map)
6
clims
slopes = cell2mat(outputs(ok,4))-cell2mat(outputs(ok,3));
mmax(slopes)
mmax(slopes./diff(r,[],2))
mmax(slopes./diff(r(ok,:),[],2))
duration = diff(r(ok,:),[],2);
open duration
mmax(duration)
mmax(slopes./duration./200*4)
mmax(shuffledSlopes./duration./200*4)
mmax(shuffledSlopes./repelem(duration,500)./200*4)
size(rEstimations,1)
slopes = (cell2mat(outputs(ok,4))-cell2mat(outputs(ok,3)))/size(rEstimations,1)*4./duration; % in m/s
mmax(slopes)
figure; hist(slopes,100);
slopes(slopes>30) = 30; slopes(slopes<-30)=-30;
hist(slopes,100);
shuffledSlopes = (cell2mat(outputs(ok,7))-cell2mat(outputs(ok,6)))/size(rEstimations,1)*4./repelem(duration,nShuffles);
mmax(shuffledSlopes)
5
clims
xlim([-1 1]*20);
name
save(['M:\home\raly\Documents\code\new\tempFiles\' name '-Cell'],'bursts','rEstimationsCell','outputs','intervals','pos2');
d M:\home\raly\Documents\code\new\tempFiles\
cd M:\home\raly\Documents\code\new\tempFiles\
name = [sessionID '_OML_off_Cell.mat'];
save(['M:\home\raly\Documents\code\new\tempFiles\' name],'bursts','rEstimationsCell','outputs','intervals','pos2');
outputsOff = outputs;
outputsOff = outputs; rEstimationsCellOff = rEstimationsCell;
size(outputs
size(outputs)
name = [sessionID '_OML_on_Cell.mat'];
intervals = SubtractIntervals(run,SubtractIntervals([0 Inf],stim)); % only stim run periods
[rEstimations] = ReconstructPosition(pos1,spikes,rWindows,'training',intervals);
tic
rEstimationsCell = cell(size(r,1),1);
for i=1:size(r,1)
rEstimationsCell{i,1} = rEstimations(:,rID==i);
end
toc
outputs = cell(size(r,1),13);
parfor i=1:5000
if sum(rID==i & ~empty)>=3
outputs(i,:) = temp(rEstimationsCell{i},'circular','off');
end
if rem(i,100)==0, display([num2str(i) '/' num2str(size(r,1))]); end
end
toc
parfor i=5001:size(r,1)
if sum(rID==i & ~empty)>=3
outputs(i,:) = temp(rEstimationsCell{i},'circular','off');
end
if rem(i,100)==0, display([num2str(i) '/' num2str(size(r,1))]); end
end
toc
tic
outputs = cell(size(r,1),13);
parfor i=1:5000
if sum(rID==i & ~empty)>=3
outputs(i,:) = temp(rEstimationsCell{i},'circular','off');
end
if rem(i,100)==0, display([num2str(i) '/' num2str(size(r,1))]); end
end
toc
open scores
open score
sum(pre)
correctedScore = (score - nanmedian(score(pre)))./(quantile(score(pre,0.75))-quantile(score(pre,0.25)));
quantile(score(pre,0.75))
quantile(score(pre,0.75)
quantile(score(pre,0.75))
correctedScore = (score - nanmedian(score(pre)))./(quantile(score(pre),0.75)-quantile(score(pre),0.25));
ope correctedScore
open correctedScore
figure; anovabar(correctedScore,-pre+post);
figure; anovabar(correctedScore(pre|post),post(pre|post))
kruskalbar(correctedScore(pre|post),post(pre|post))
signrank(correctedScore(post))
signrank(correctedScore(pre))
ttest(correctedScore(pre))
kruskalbar(correctedScore(pre|post),post(pre|post))
anovabar(correctedScore(pre|post),post(pre|post))
kruskalbar(correctedScore(pre|post),post(pre|post))
sig = pValue; sig(~isnan(sig)) = pValue(~isnan(sig))<0.05;
sum(sig)
nansum(sig==1)
nansum(sig==0)
sigp = sig; sigp(slope<0) = 0;
nansum(sigp==0)
nansum(sigp==1)
ok = (pre|post) & ~isnan(sigp);
Accumulate([post(ok)+1,sigp(ok)+2])
88/820
541/4886
(541/4886) ./ ( 88/820)
hist(score(pre),100);
correctedScore = (score - nanmean(score(pre)))./nanstd(score(pre));
anovabar(correctedScore(pre|post),post(pre|post))
ttest(correctedScore(pre))
signrank(correctedScore(pre))
signrank(correctedScore(post))
kruskalbar(correctedScore(pre|post),post(pre|post))
correctedScore = (score - nanmedian(score(pre)))./(quantile(score(pre),0.75)-quantile(score(pre),0.25));
kruskalbar(correctedScore(pre|post),post(pre|post))
name = [sessionID '_OML_on_Cell.mat'];
intervals = SubtractIntervals(run,SubtractIntervals([0 Inf],stim)); % only stim run periods
[rEstimations] = ReconstructPosition(pos1,spikes,rWindows,'training',intervals);
tic
rEstimationsCell = cell(size(r,1),1);
for i=1:size(r,1)
rEstimationsCell{i,1} = rEstimations(:,rID==i);
end
toc
outputs = cell(size(r,1),13);
parfor i=1:1000
if sum(rID==i & ~empty)>=3
outputs(i,:) = temp(rEstimationsCell{i},'circular','off');
end
if rem(i,100)==0, display([num2str(i) '/' num2str(size(r,1))]); end
end
toc
parfor i=1001:size(r,1)
if sum(rID==i & ~empty)>=3
outputs(i,:) = temp(rEstimationsCell{i},'circular','off');
end
if rem(i,100)==0, display([num2str(i) '/' num2str(size(r,1))]); end
end
toc
save(['M:\home\raly\Documents\code\new\tempFiles\' name],'bursts','rEstimationsCell','outputs','intervals','pos2');
outputsOn = outputs;
4
1/0.8
PlotColorMap(zmap,~isnan(zmap),'x',linspace(-100,100,nBins(1)),'y',linspace(0,1/0.8,nBins(2)));
PlotHVLines(0,'v','w--','linewidth',2);
xlim([-1 1]*20); ylim([0 1]);
clim
clims
clim
clims([-1 1]*2);
clims([-1 1]*5);
clims([-1 1]*10);
clims([-1 1]*20);
clims([-1 1]*40);
clims([0 1]*40);
clims([0 1]*50);
mmax(scores)
cellfun(@(x) x(1)>x(2), cat(2,outputsOn{:,1},outputsOff{:,1});
cellfun(@(x) x(1)>x(2), cat(2,outputsOn{:,1},outputsOff{:,1}));
this = cat(2,outputsOn(:,1),outputsOff(:,1));
open this
reOn = false(size(outputsOn,1),1);
reOff = false(size(outputsOn,1),1);
for i=1:size(outputsOn,1), reOn(i,1) = outputsOn{i,1}>outputsOff{i,1}; end
for i=1:size(outputsOn,1), try reOn(i,1) = outputsOn{i,1}>outputsOff{i,1}; end; end
for i=1:size(outputsOn,1), try reOff(i,1) = outputsOn{i,1}<outputsOff{i,1}; end; end
sum(reOn)
sum(reOff)
clims([0 10]);
clims([0 50]);
2
sum(sum(sum(mapsShuffled,3)==0))
sum(mapsShuffled(1,1,:))
sum(sum(map))
mShuffled = Smooth(nanmean(mapsShuffled,3),smooth); stdShuffled = Smooth(nanstd(mapsShuffled,[],3),smooth);
figure; hist(mShuffled(:),100)
hist(mShuffled(:),1000)
hist(stdShuffled(:),1000)
hist(Restrict(stdShuffled(:),[-1 0.01]),1000)
hist(Restrict(stdShuffled(:),[-1 0.0001]),1000)
hist(Restrict(stdShuffled(:),[-1 0.00001]),1000)
sum(stdShuffled(:)==0)
hist(Restrict(stdShuffled(:),[-1 0.0000001]),1000)
min(stdShuffled(:))
min(stdShuffled(:))>eps
eps
52
clims([0 10]);
script_Can_replay
clims([0 5]);
clims([0 2]);
figure; PlotColorMap(zmap);
clim
clim([0 10])
clim([0 5])
PlotColorMap(zmap,'x',linspace(-100,100,nBins(1)),'y',linspace(0,1/0.8,nBins(2)));
condition
clim([0 5])
figure; plot(map(:,251))
hold all
plot(mShuffled(:,251))
clf
plot(map(:,251)-mShuffled(:,251))
plot((map(:,251)-mShuffled(:,251)))
hold all
plot((stdShuffled(:,251)))
SideAxes(gca,'right',0.5);
plot((map(:,251)-mShuffled(:,251))./stdShuffled(:,251)))
plot((map(:,251)-mShuffled(:,251))./stdShuffled(:,251));
plot(0.15*(map(:,251)-mShuffled(:,251))./stdShuffled(:,251));
plot(0.15*(map(:,101)-mShuffled(:,101))./stdShuffled(:,101));
q = linspace(-100,100,nBins(1));
q(101)
FindClosest(q,-1.8p)
FindClosest(q,-1.890)
j=246; plot(0.15*(map(:,j)-mShuffled(:,j))./stdShuffled(:,j));
plot((stdShuffled(:,j)))
plot((map(:,j)-mShuffled(j)));
plot((map(:,j)-mShuffled(:,j)));
hold all
plot((stdShuffled(:,j)))
j=246; plot(0.15*(map(:,j)-mShuffled(:,j))./stdShuffled(:,j));
6
465
clims([0 10]);
clim
clims([0 2]);
clims([0 1]);
q = nansmooth(zmap,smooth);
open q
open smoothed
clims([0 1]);
clims([0 0.02]);
clims([0 0.01]);
clims([-1 0.01]);
clims([-5 0.01]);
PlotColorMap(Smooth(smoothed,0),Smooth(double(~isnan(zmap)),smooth),'x',linspace(-100,100,nBins(1)),'y',linspace(0,1/0.8,nBins(2)));
xlim([-1 1]*20); ylim([0 1]);
clim
clims([0 10]);
clims([0 5]);
clims([0 10]);
clims([0 20]);
sum(zmap(:)==Inf)
zmap = (map-mShuffled)./stdShuffled;
sum(zmap(:)==Inf)
max(zmap(zmap~=Inf))
condition = 2;
k=2;
if condition==1
outputs = outputsOn;
bad = cellfun(@isempty,outputs(:,1));
bad = bad | ~(reOn);
else
outputs = outputsOff;
bad = cellfun(@isempty,outputs(:,1));
bad = bad | ~(reOff);
end
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
pValues = cell2mat(outputs(ok,2));
figure; plot(pValue,scores,'.')_;
figure; plot(pValue,scores,'.');
plot(pValues,scores,'.');
corr(pValues,scores)
corr(pValues,slopes)
plot(pValues,slopes,'.');
plot(pValues,abs(slopes),'.');
mmax(slopes)
these{condition,k} = scores;
g = Group(these{:});
open these
Accumulate(g(:,2))
anovabar(g(:,1),g(:,2))
clf
anovabar(g(:,1),g(:,2))
anovabar(g(:,1)-0.5,g(:,2))
g = Group(these{:});
clf
anovabar(g(:,1)-0.5,g(:,end))
anovabar(g(:,2),g(:,end))
anovabar(g(:,3),g(:,end))
anovabar(g(:,3)>0,g(:,end))
anovabar(g(:,3)>0 & g(:,2)<0.05,g(:,end))
Accumulate(g(:,end),g(:,3)>0 & g(:,2)<0.05)
Accumulate(g(:,end),~(g(:,3)>0 & g(:,2)<0.05))
Accumulate(g(:,end),g(:,3)>0 & g(:,2)<0.05)./Accumulate(g(:,end),~(g(:,3)>0 & g(:,2)<0.05))
clf
Accumulate(g(:,end),g(:,3)>0 & g(:,2)<0.05)./Accumulate(g(:,end),~(g(:,3)>0 & g(:,2)<0.05))
g = Group(these{:});
Accumulate(g(:,end),g(:,3)>0 & g(:,2)<0.05)./Accumulate(g(:,end),~(g(:,3)>0 & g(:,2)<0.05))
Accumulate(g(:,end))
kruskalbar(g(:,1),g(:,end))
kruskalbar(g(:,1)-0.4,g(:,end))
g = Group(these{1,1},these{1,2});
56
anovabar(g(:,3)>0 & g(:,2)<0.05,g(:,end))
g = Group(these{2,1},these{2,2});
anovabar(g(:,3)>0 & g(:,2)<0.05,g(:,end))
anovabar(g(:,2)<0.05,g(:,end))
anovabar(g(:,3)>0,g(:,end))
open g
anovabar(g(:,2)<0.05,g(:,end))
kruskalbar(g(:,2)<0.05,g(:,end))
kruskalbar(g(:,1),g(:,end))
kruskalbar(g(:,1) - nanmedian(g(g(:,end)==1)),g(:,end))
anovabar(g(:,1) - nanmedian(g(g(:,end)==1)),g(:,end))
anovabar(zscore(g(:,1)),g(:,end))
kruskalbar(zscore(g(:,1)),g(:,end))
anovabar(zscore(g(:,1)),g(:,end))
Dist(100,g,'grouped');
Dist(1000,g,'grouped');
[h,ht] = Dist(1000,g,'grouped');
plot(ht,Smooth(h,[1 0]));
plot(ht,Smooth(h,[2 0]));
plot(ht,Smooth(h,[5 0]));
plot(ht,Smooth(h,[10 0]));
plot(ht,Smooth(h,[20 0]));
plot(ht,Smooth(h,[50 0]));
plot(ht,Smooth(h,[0 0]));
plot(ht,cumsum(Smooth(h,[0 0])));
g = Group(these{:});
[h,ht] = Dist(1000,g,'grouped');
plot(ht,cumsum(Smooth(h,[0 0])));
legend('on pre','off pre','on post','off post')
clf
anovabar(g(:,1) - nanmedian(g(g(:,end)==1)),g(:,end))
anovabar(g(:,1) - nanmedian(g(g(:,end)==2)),g(:,end))
kruskalbar(g(:,1) - nanmedian(g(g(:,end)==2)),g(:,end))
kruskalbar(g(:,1) - nanmedian(g(g(:,end)==1)),g(:,end))
gg = g; gg(ismember(gg(:,end),[1 3]),1) = (gg(ismember(gg(:,end),[1 3]),1) - nanmedian(gg(gg(:,end)==1,1)))./diff(quantile(gg(gg(:,end)==1,1),[0.25 0.75]));
gg(ismember(gg(:,end),[2 4]),1) = (gg(ismember(gg(:,end),[2 4]),1) - nanmedian(gg(gg(:,end)==2,1)))./diff(quantile(gg(gg(:,end)==2,1),[0.25 0.75]));
anovabar(gg(:,1),g(:,end))
kruskalbar(gg(:,1),g(:,end))
kruskalbar(gg(:,1)*100,g(:,end))
edit SaveFig
anovabar(gg(:,1)*100,g(:,end))
kruskalbar(gg(:,1)*100,g(:,end))
kruskalbar(abs(gg(:,3)),g(:,end))
anovabar(sign(g(:,3)),g(:,end))
ok = g(:,2)<0.05;
anovabar(sign(g(ok,3)),g(ok,end))
sum(reOn)
sum(reOff)
hm = cell2mat([outputsOn(:,1:2) outputsOff(:,1:2)]);
open hm
sum(hm(:,3)>hm(:,1))
sum(hm(:,2)>hm(:,4))
sum(hm(:,2)<hm(:,4))
sum(hm(:,[2 4])<0.05)
sum(sum(hm(:,2)<hm(:,4),2)==2)
sum(sum(hm(:,[2 4])<0.05,2)==2)
hmhm = hm(sum(hm(:,[2 4])<0.05,2)==2)
hmhm = hm(sum(hm(:,[2 4])<0.05,2)==2,:);
open hmhm
sum(~bad)
hm(:,5) = find(~bad);
hmhm = hm(sum(hm(:,[2 4])<0.05,2)==2,:);
intervals = r; intervals0 = intervals;
intervals(1:end-1,2) = min([intervals0(1:end-1,2)+leeway mean([intervals0(1:end-1,2) intervals0(2:end,1)],2)],[],2); intervals(end,2) = intervals0(end,2) + leeway;
intervals(2:end,1) = max([intervals0(2:end,1)-leeway mean([intervals0(2:end,1) intervals0(1:end-1,2)],2)],[],2); intervals(1,1) = intervals0(1,1) - leeway;
nonstimIntervals = SubtractIntervals(run,stim);
stimIntervals = SubtractIntervals(run,stim);
figure; PlotColorMap(ReconstructPosition(pos2,spikes,intervals,'training',nonstimIntervals));
figure; PlotColorMap(ReconstructPosition(pos2,spikes,rWindows(rID==44,:)));
hm = cell2mat([outputsOn(:,1:3) outputsOff(:,1:3)]);
hm = cell2mat([outputsOn(:,1:4) outputsOff(:,1:4)]);
hmhm = hm(sum(hm(:,[2 6])<0.05,2)==2,:);
hm(:,end+1) = find(~bad);
hmhm = hm(sum(hm(:,[2 6])<0.05,2)==2,:);
hold all
plot([1 41],[23 24],'w--','linewidth',2);
PlotColorMap(ReconstructPosition(pos2,spikes,rWindows(rID==44,:),'training',nonstimIntervals));
PlotColorMap(ReconstructPosition(pos1,spikes,rWindows(rID==44,:),'training',stimIntervals));
size(duration)
size(bursts0
size(bursts)
duration = diff(bursts(:,[1 3]),[],2);
anovabar(duration(hm(:,9)),hm(:,2)<0.05);
clf
anovabar(duration(hm(:,9)),hm(:,2)<0.05);
plot(duration(hm(:,9)),hm(:,2)<0.05);
plot(duration(hm(:,9)),hm(:,2)<0.05,'.');
duration(44)
duration(63)
clf
PlotColorMap(ReconstructPosition(pos1,spikes,rWindows(rID==63,:),'training',stimIntervals));
hold all
plot(xlim+[0.5 -0.5],[45 40],'w--','linewidth',2);
plot(xlim+[0.5 -0.5],[45 40]+15,'w--','linewidth',2);
plot(xlim+[0.5 -0.5],[45 40]-15,'w--','linewidth',2);
size(rWindows)\
size(rWindows)
max(rID)
max(hm(:,end))
PlotColorMap(ReconstructPosition(pos2,spikes,rWindows(rID==63,:),'training',nonstimIntervals));
hold all
plot(xlim+[0.5 -0.5],[36 136],'w--','linewidth',2);
PlotColorMap(ReconstructPosition(pos1,spikes,rWindows(rID==63,:),'training',stimIntervals));
plot(xlim+[0.5 -0.5],[45 40],'w--','linewidth',2);
this = ReconstructPosition(pos1,spikes,rWindows(rID==63,:),'training',stimIntervals);
open this
PlotColorMap(rEstimationsCell{63})
PlotColorMap(rEstimationsCellOff{63})
PlotColorMap(ReconstructPosition(pos2,spikes,rWindows(rID==63,:),'training',nonstimIntervals));
rEstimationsCellOn = rEstimationsCell;
nonstimIntervals = SubtractIntervals(run,stim);
stimIntervals = SubtractIntervals(run,SubtractIntervals([0 Inf],stim)); % only stim run periods
PlotColorMap(ReconstructPosition(pos2,spikes,rWindows(rID==63,:),'training',nonstimIntervals));
PlotColorMap(rEstimationsCellOff{63})
PlotColorMap(rEstimationsCellOn{63})
PlotColorMap(ReconstructPosition(pos1,spikes,rWindows(rID==63,:),'training',stimIntervals));
hold on
plot(xlim+[0.5 -0.5],[45 40],'w--','linewidth',2);
15/200*6
15/200*4
10/200*4
PlotColorMap(rEstimationsCellOff{63})
plot(xlim+[0.5 -0.5],[36 136],'w--','linewidth',2);
PlotColorMap(ReconstructPosition(pos2,spikes,rWindows(rID==63,:),'training',nonstimIntervals));
PlotColorMap(ReconstructPosition(pos2,spikes,rWindows(rID==63,:),'training',stimIntervals));
PlotColorMap(ReconstructPosition(pos2,spikes,rWindows(rID==63,:),'training',run));
PlotColorMap(rEstimationsCellOff{63})
PlotColorMap(ReconstructPosition(pos2,spikes,rWindows(rID==63,:),'training',run));
PlotColorMap(rEstimationsCellOff{63})
PlotColorMap(ReconstructPosition(pos2,spikes,rWindows(rID==63,:),'training',run));
PlotColorMap(ReconstructPosition(pos2,spikes,rWindows(rID==63,:),'training',rWindows));
PlotColorMap(rEstimationsCellOff{63})
%% First, off
name = [sessionID '_OML_off_Cell.mat'];
nonstimIntervals = SubtractIntervals(run,stim);
[rEstimations] = ReconstructPosition(pos2,spikes,rWindows,'training',nonstimIntervals);
rEstimationsCell = cell(size(r,1),1);
for i=1:size(r,1)
rEstimationsCell{i,1} = rEstimations(:,rID==i);
end
outputs = cell(size(r,1),13);
parfor i=1:size(r,1)
if sum(rID==i & ~empty)>=3
outputs(i,:) = temp(rEstimationsCell{i},'circular','off');
end
if rem(i,100)==0, display([num2str(i) '/' num2str(size(r,1))]); end
end
toc
save(['M:\home\raly\Documents\code\new\tempFiles\' name],'bursts','rEstimationsCell','outputs','intervals','pos2');
% check if ok, then replace OutputsOff
edit clim
45
PlotColorMap(rEstimationsCellOff{63})
PlotColorMap(rEstimationsCell{63})
PlotColorMap(ReconstructPosition(pos2,spikes,rWindows(rID==63,:),'training',rWindows));
PlotColorMap(ReconstructPosition(pos2,spikes,rWindows(rID==63,:),'training',nonstimIntervals));
PlotColorMap(rEstimationsCell{63})
PlotColorMap(rEstimationsCellOff{63})
outputsOff = outputs; rEstimationsCellOff = rEstimationsCell;
4
g = Group(these{:});
clf
anovabar(sign(g(ok,3)),g(ok,end))
ok = g(:,2)<0.05;
anovabar(sign(g(ok,3)),g(ok,end))
gg = g; gg(ismember(gg(:,end),[1 3]),1) = (gg(ismember(gg(:,end),[1 3]),1) - nanmedian(gg(gg(:,end)==1,1)))./diff(quantile(gg(gg(:,end)==1,1),[0.25 0.75]));
kruskalbar(abs(gg(:,3)),g(:,end))
kruskalbar(gg(:,1)*100,g(:,end))
gg = g; gg(ismember(gg(:,end),[1 3]),1) = (gg(ismember(gg(:,end),[1 3]),1) - nanmedian(gg(gg(:,end)==1,1)))./diff(quantile(gg(gg(:,end)==1,1),[0.25 0.75]));
gg(ismember(gg(:,end),[2 4]),1) = (gg(ismember(gg(:,end),[2 4]),1) - nanmedian(gg(gg(:,end)==2,1)))./diff(quantile(gg(gg(:,end)==2,1),[0.25 0.75]));
kruskalbar(gg(:,1)*100,g(:,end))
%-- 3/17/2022 10:23 AM --%
45
56
45
ylabel('Decoded position');
load(fullfile(basepath,[basename '.MergePoints.events.mat']),'MergePoints');
sleep = ConsolidateIntervals(MergePoints.timestamps(cellfun(@(x) any(strfind(lower(x),'sleep')),MergePoints.foldernames),:));
load(fullfile(basepath,[basename '.SleepState.states.mat']));
sws = SleepState.ints.NREMstate;
preSleep = SubtractIntervals(sleep(1,:), SubtractIntervals([0 Inf],sws));
postSleep = SubtractIntervals(sleep(2:end,:), SubtractIntervals([0 Inf],sws));
bursts = FindBursts(spikes,'thresholds',[0 3],'smooth',0.01);
bursts0 = bursts;
bursts = FindBursts(Restrict(spikes,sws,'shift','on'),'thresholds',[0 3],'smooth',0.01);
open bursts
bursts(:) = Unshift(bursts(:),sws);
PETH(bursts(:,2),bursts0(:,2))
open bursts0
g = Group(bursts0(:,2), bursts(:,2));
RasterPlot(g);
g = Group(bursts0(:,1), bursts(:,1));
RasterPlot(g,'g');
open g
clf
RasterPlot(Group(bursts0(:,2), bursts(:,2)));
RasterPlot(Group(bursts0(:,2), bursts(:,2)),1);
RasterPlot(Group(bursts0(:,1), bursts(:,1)),1,'g');
hold all
RasterPlot(Group(bursts0(:,2), bursts(:,2)),1);
RasterPlot(Group(bursts0(:,2), bursts(:,2)),1,'b');
RasterPlot(Group(bursts0(:,3), bursts(:,3)),1,'r');
close = abs(bursts(FindClosest(bursts(:,2),bursts0(:,2)),2) - bursts0(:,2));
open close
close = abs(bursts(FindClosest(bursts(:,2),bursts0(:,2)),2) - bursts0(:,2))<0.02;
Portion(close)
sum(close)
notclose = abs(bursts0(FindClosest(bursts0(:,2),bursts(:,2)),2) - bursts(:,2))>0.02;
sum(notclose)
find(notclose,5)
xlim([-5 5]+burst0(9,2))
xlim([-5 5]+bursts0(9,2))
bursts0(9,2)
size(bursts0)
size(notclose)
xlim([-5 5]+bursts(9,2))
bursts(9,2)
xlim([-1 1]*2+bursts(9,2))
sum(~notclose)
sum(close)
d = [diff(bursts(~notclose,[1 3]),[],2) diff(bursts0(close,[1 3]),[],2)];
open d
nanmean(d)
diff(nanmean(d),[],2)
diff(nanmedian(d),[],2)
figure; PlotXY(d);
PlotXY(d,'.');
PlotXY(Restrict(d,[0 5]),'.');
nanmean(Restrict(d,[0 5]))
sum(pre)
pre = InIntervals(bursts0,preSleep);
post = InIntervals(bursts0,postSleep);
sum(pre|post)
clear all
45
5
clear intervals
546
figure; PlotXY(pos1);
hold all
PlotXY(pos2);
close
figure; PlotXY(pos1,'.');
hold all
PlotXY(pos2,'r.');
PlotIntervals(stim,'color','g');
figure; PlotXY(positionTrials{1});
hold all
PlotXY(positionTrials{2});
figure; PlotXY(behavior.positionTrials{1});
hold all
PlotXY(behavior.positionTrials{2});
figure; PlotXY(behavior.positionTrials{1},'.');
hold all
PlotXY(behavior.positionTrials{2},'r.');
PlotIntervals(stim,'color','g');
PlotIntervals(run,'color','k');
stimIntervals = SubtractIntervals(run,SubtractIntervals([0 Inf],stim)); % only stim run periods
[estimations,actual,errors,average] = ReconstructPosition(pos1,spikes,windows,'training',stimIntervals,'id',id);
clf
PlotColorMap(estimations)
hold all
plot(actual*100,'w');
figure; hist(windows);
hist(windows(:,1));
sum(InIntervals(windows,run))
sum(~InIntervals(windows,run))
clf
PlotXY(Restrict(pos1,windows));
[in,w] = InIntervals(pos1,windows0;
[in,w] = InIntervals(pos1,windows);
ii = interp1(pos1(:,1),pos1(:,1),mean(windows,2));
45
clims
[~,m] = max(Smooth(curves,[0 2],'type','cc'),[],2);
PlotColorMap(sortby(nanzscore([curves],[],2),m))
pos = sortrows([pos1; pos2(:,1) 1+pos2(:,2)]);
hist(pos(:,2))
kkeybaord
kkeyboard
cd D:
465
45
length(tolerances)
45
2
open pf
open fp
k
j
tolerant = TolerantDefragmentSSA(detected,[],timescale,-Inf,1,0,tolerance);
mean(~ismember(comb(:,ok),triplets(:,ok),'rows'))
nComb = nchoosek(sum(ok),3)
sum(sum(triplets(:,ok),2)==3)
1 - sum(sum(triplets(:,ok),2)==3)/nComb
tic; for i=1:100,  comb = SSA_SplitAssemblies(combinations(i,:),3); mean(~ismember(comb(:,ok),triplets(:,ok),'rows'))<=tolerance; end; toc
i=1;
tic; for j=1:100,  comb = SSA_SplitAssemblies(combinations(i,:),3); mean(~ismember(comb(:,ok),triplets(:,ok),'rows'))<=tolerance; end; toc
tic; for j=1:1000,  comb = SSA_SplitAssemblies(combinations(i,:),3); mean(~ismember(comb(:,ok),triplets(:,ok),'rows'))<=tolerance; end; toc
tic; for j=1:10000,  comb = SSA_SplitAssemblies(combinations(i,:),3); mean(~ismember(comb(:,ok),triplets(:,ok),'rows'))<=tolerance; end; toc
tic; for j=1:10000,  nComb = nchoosek(sum(ok),3); 1 - sum(sum(triplets(:,ok),2)==3)/nComb<=tolerance end; toc
tic; for j=1:10000,  nComb = nchoosek(sum(ok),3); 1 - sum(sum(triplets(:,ok),2)==3)/nComb<=tolerance; end; toc
tic; for j=1:100000,  nComb = nchoosek(sum(ok),3); 1 - sum(sum(triplets(:,ok),2)==3)/nComb<=tolerance; end; toc
tic; for j=1:100000,  comb = SSA_SplitAssemblies(combinations(i,:),3); mean(~ismember(comb(:,ok),triplets(:,ok),'rows'))<=tolerance; end; toc
kkeyboard
45
splitg = SSA_SplitAssemblies(groundTruthAssemblies',3);
45
dbcont
45
load('final_assemblies.mat')
size(assemblies)
load('Harris.mat', 'spikes')
max(spikes(:,2))
clear all
load('Harris.mat')
load('final_assemblies.mat')
load('workspace_cycle_5.mat', 'threshold')
load('workspace_cycle_5.mat', 'nMin')
load('workspace_cycle_5.mat', 'windowSize')
verbose = 1;
sum(spikes(spikes(:,2)==0))
sum(spikes(:,2)==0)
spikes(spikes(:,2)==0,:) = [];
clear ans
save('assemblies2degrafment.mat');
save('data4assemblies2degrafment.mat');
45
threshold
nMin
windowSize
clear all
load('assemblies.mat','assemblies')
load('data4assemblies2degrafment', 'threshold','nMin','spikes','verbose','windowSize');
tolerance = 0.1; %??
sum(spikes(:,2)==0)
PlotColorMap(sortby(nanzscore([curves],[],2),m))
% compute firing curves
clear curves curves1 curves2
nonstimIntervals = SubtractIntervals(run,stim); stimIntervals = SubtractIntervals(run,SubtractIntervals([0 Inf],stim)); % only stim run periods
for i=1:max(spikes(:,2))
map = FiringCurve(Restrict(pos1,stimIntervals,'shift','on'),Restrict(spikes(spikes(:,2)==i),stimIntervals,'shift','on'));
curves1(i,:) = map.rate;
map = FiringCurve(Restrict(pos2,nonstimIntervals,'shift','on'),Restrict(spikes(spikes(:,2)==i),nonstimIntervals,'shift','on'));
curves2(i,:) = map.rate;
map = FiringCurve(Restrict(pos,run,'shift','on'),Restrict(spikes(spikes(:,2)==i),run,'shift','on'));
curves(i,:) = map.rate;
end
[~,m] = max(Smooth(curves,[0 2],'type','cc'),[],2);
PlotColorMap(sortby(nanzscore([curves],[],2),m))
clf
PlotXY(pos2);
clf
PlotColorMap(sortby(nanzscore([curves2],[],2),m))
PlotColorMap(sortby(nanzscore([curves1 curves2 curves],[],2),m))
pos = sortrows([Restrict(pos1,stimIntervals) ; Restrict(pos2(:,1) 1+pos2(:,2),nonstimIntervals)]);
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
[~,m] = max(Smooth(curves,[0 2],'type','cc'),[],2);
PlotColorMap(sortby(nanzscore([curves],[],2),m))
PlotColorMap(sortby(nanzscore([curves1 curves2 curves],[],2),m))
pos = sortrows([Restrict(pos1,stimIntervals) ; Restrict(pos2(:,1) 1+pos2(:,2),nonstimIntervals)]);
pos = sortrows([Restrict(pos1,stimIntervals) ; Restrict([pos2(:,1) 1+pos2(:,2)],nonstimIntervals)]);
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
[~,m] = max(Smooth(curves,[0 2],'type','cc'),[],2);
PlotColorMap(sortby(nanzscore([curves],[],2),m))
PlotColorMap(sortby(nanzscore([curves1 curves2 curves],[],2),m))
clim
clims([0 0.02]);
clf
PlotXY(pos);
PlotIntervals(nonstimIntervals,'color','k');
PlotIntervals(run,'color','k');
sum(IntervalsIntersect,run,sws)
sum(IntervalsIntersect(run,sws))
sum(IntervalsIntersect(run,swt))
load(fullfile(basepath,[basename '.MergePoints.events.mat']),'MergePoints');
sleep = ConsolidateIntervals(MergePoints.timestamps(cellfun(@(x) any(strfind(lower(x),'sleep')),MergePoints.foldernames),:));
load(fullfile(basepath,[basename '.SleepState.states.mat']));
sws = SleepState.ints.NREMstate;
preSleep = SubtractIntervals(sleep(1,:), SubtractIntervals([0 Inf],sws));
postSleep = SubtractIntervals(sleep(2:end,:), SubtractIntervals([0 Inf],sws));
sum(IntervalsIntersect(run,sleep))
run(IntervalsIntersect(run,sleep),:) = [];
pos0 = sortrows([behavior.positionTrials{1}; behavior.positionTrials{2}(:,1) 2-behavior.positionTrials{2}(:,2)]);
pos1 = Restrict(behavior.positionTrials{1},run);
pos2 = Restrict(behavior.positionTrials{2},run); pos2(:,2) = 1-pos2(:,2);
% the first direction is stim (for labels in later analyses)
if sum(InIntervals(pos2,stim)) > sum(InIntervals(pos1,stim)) % otherwise, swap them
pos1 = pos2; pos2 = Restrict(behavior.positionTrials{1},run);  pos2(:,2) = 1-pos2(:,2);
end
stimIntervals = SubtractIntervals(run,SubtractIntervals([0 Inf],stim)); % only stim run periods
nonstimIntervals = SubtractIntervals(run,stim);
pos = sortrows([Restrict(pos1,stimIntervals) ; Restrict([pos2(:,1) 1+pos2(:,2)],nonstimIntervals)]);
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
[~,m] = max(Smooth(curves,[0 2],'type','cc'),[],2);
PlotColorMap(sortby(nanzscore([curves],[],2),m))
PlotColorMap(sortby(nanzscore([curves1 curves2 curves],[],2),m))
clim
[parentFolder,dayName] = fileparts(basepath);
[~,projectName] = fileparts(parentFolder);
sessionID = [projectName '_' dayName];
cd(basepath)
edit preprocessSession.m
SleepScoreMaster(basepath,'noPrompts',true,'rejectchannels',1+[29 7 12 6 13 5 14 4 15 8 11 9 10 46 19 45 18 44 20 42]); % takes lfp in base 0
name = [sessionID '_OML_off_Cell.mat'];
[rEstimations] = ReconstructPosition(pos2,spikes,rWindows,'training',nonstimIntervals);
rEstimationsCell = cell(size(r,1),1);
for i=1:size(r,1)
rEstimationsCell{i,1} = rEstimations(:,rID==i);
end
outputs = cell(size(r,1),13);
parfor i=1:size(r,1)
if sum(rID==i & ~empty)>=3
outputs(i,:) = temp(rEstimationsCell{i},'circular','off');
end
if rem(i,100)==0, display([num2str(i) '/' num2str(size(r,1))]); end
end
toc
save(['M:\home\raly\Documents\code\new\tempFiles\' name],'bursts','rEstimationsCell','outputs','intervals','pos2');
outputsOff = outputs; rEstimationsCellOff = rEstimationsCell;
name = [sessionID '_OML_on_Cell.mat'];
[rEstimations] = ReconstructPosition(pos1,spikes,rWindows,'training',stimIntervals);
tic
rEstimationsCell = cell(size(r,1),1);
for i=1:size(r,1)
rEstimationsCell{i,1} = rEstimations(:,rID==i);
end
toc
outputs = cell(size(r,1),13);
parfor i=1001:size(r,1)
if sum(rID==i & ~empty)>=3
outputs(i,:) = temp(rEstimationsCell{i},'circular','off');
end
if rem(i,100)==0, display([num2str(i) '/' num2str(size(r,1))]); end
end
toc
save(['M:\home\raly\Documents\code\new\tempFiles\' name],'bursts','rEstimationsCell','outputs','intervals','pos2');
outputsOn = outputs;
% [r,p,st,sp,rShuffled,aShuffled,bShuffled,c,cShuffled,jump,jumpShuffled,maxJump,maxJumpShuffled]
nBins = [101 25];
bad = cellfun(@isempty,outputs(:,1));
outputs = outputsOff;
smooth = [0 1];
clf
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
title([condition k]);
drawnow
end
end
clims([0 20]);
ThIRASA && strcmp(SWweights,'PSS')
display('Calculating Theta Metric above PSS')
%Put the LFP in the right structure format
lfp.data = thLFP;
lfp.timestamps = t_LFP;
lfp.samplingRate = sf_LFP;
%Calculate PSS
[specslope,spec] = bz_PowerSpectrumSlope(lfp,window,window-noverlap,...
'nfreqs',200,'frange',f_all,'IRASA',ThIRASA);
t_thclu = specslope.timestamps;
specdt = 1./specslope.samplingRate;
thFFTspec = specslope.resid';
thFFTspec(thFFTspec<0)=0;
IRASAsmooth_th = spec.IRASAsmooth';
thFFTspec_raw = 10.^spec.amp';
% Remove transients before calculating SW histogram
zFFTspec = NormToInt(spec.amp,'modZ');
totz = NormToInt(abs(sum(zFFTspec,2)),'modZ');
badtimes_TH = find(totz>3);
thFFTfreqs = specslope.freqs';
thfreqs = (thFFTfreqs>=f_theta(1) & thFFTfreqs<=f_theta(2));
thratio = max((thFFTspec(thfreqs,:)),[],1);
thratio(badtimes_TH) = nan;     %Remove transients
thratio = smooth(thratio,smoothfact./specdt);
display('Calculating Theta Metric above PSS')
%Put the LFP in the right structure format
lfp.data = thLFP;
lfp.timestamps = t_LFP;
lfp.samplingRate = sf_LFP;
%Calculate PSS
[specslope,spec] = bz_PowerSpectrumSlope(lfp,window,window-noverlap,...
'nfreqs',200,'frange',f_all,'IRASA',ThIRASA);
45
t_thclu = specslope.timestamps;
specdt = 1./specslope.samplingRate;
thFFTspec = specslope.resid';
thFFTspec(thFFTspec<0)=0;
IRASAsmooth_th = spec.IRASAsmooth';
thFFTspec_raw = 10.^spec.amp';
% Remove transients before calculating SW histogram
zFFTspec = NormToInt(spec.amp,'modZ');
totz = NormToInt(abs(sum(zFFTspec,2)),'modZ');
badtimes_TH = find(totz>3);
thFFTfreqs = specslope.freqs';
thfreqs = (thFFTfreqs>=f_theta(1) & thFFTfreqs<=f_theta(2));
thratio = max((thFFTspec(thfreqs,:)),[],1);
thratio(badtimes_TH) = nan;     %Remove transients
thratio = smooth(thratio,smoothfact./specdt);
ignoretimeIDX = InIntervals(t_clus,ignoretime) | isnan(broadbandSlowWave) | isnan(thratio);
dbcont
5
45
load(fullfile(basepath,[basename '.MergePoints.events.mat']),'MergePoints');
sleep = ConsolidateIntervals(MergePoints.timestamps(cellfun(@(x) any(strfind(lower(x),'sleep')),MergePoints.foldernames),:));
load(fullfile(basepath,[basename '.SleepState.states.mat']));
sws = SleepState.ints.NREMstate;
preSleep = SubtractIntervals(sleep(1,:), SubtractIntervals([0 Inf],sws));
postSleep = SubtractIntervals(sleep(2:end,:), SubtractIntervals([0 Inf],sws));
45
rt6
i
parfor i=3442:size(r,1)
rEstimationsCell{i,1} = rEstimations(:,rID==i);
end
outputs = cell(size(r,1),13);
parfor i=1:size(r,1)
if sum(rID==i & ~empty)>=3
outputs(i,:) = FindReplayScoreCell(rEstimationsCell{i},'circular','off');
end
if rem(i,100)==0, display([num2str(i) '/' num2str(size(r,1))]); end
end
toc
save(['M:\home\raly\Documents\code\new\tempFiles\' name],'bursts','rEstimationsCell','outputs','intervals','pos2');
outputsOff = outputs; rEstimationsCellOff = rEstimationsCell;
% Then, on
45
56
bad = cellfun(@isempty,outputs(:,1));
sum(bad)
checkout = find(bad);
parfor j=1:length(checkout)
i=checkout(j);
if sum(rID==i & ~empty)>=3
outputs(i,:) = FindReplayScoreCell(rEstimationsCell{i},'circular','off');
end
if rem(i,100)==0, display([num2str(i) '/' num2str(size(r,1))]); end
end
toc
plot(checkout)
CLF
clf
plot(checkout)
plot(checkout,'.')
for i=1:4,
i
end
size(find(checkout))
parfor i=find(checkout)'
if sum(rID==i & ~empty)>=3
outputs(i,:) = FindReplayScoreCell(rEstimationsCell{i},'circular','off');
end
if rem(i,100)==0, display([num2str(i) '/' num2str(size(r,1))]); end
end
toc
toc
save(['M:\home\raly\Documents\code\new\tempFiles\' name],'bursts','rEstimationsCell','outputs','intervals','pos2');
outputsOn = outputs; rEstimationsCellOn = rEstimationsCell;
clear all
basepath = 'M:\Data\Can\OLM21\day9';
script_Can_replay
cd M:\home\raly\Documents\code\new\tempFiles\
load('OLM21_day9_OML_on_Cell.mat', 'outputs')
outputsOn = outputs; rEstimationsCellOn = rEstimationsCell;
size(rEstimations,1)*4./repelem(duration,nShuffles)
size(duration)
duration
duration = diff(r(ok,:),[],2);
plot(bad)
condition
sum(ok)
Portion(bad)
bad = cellfun(@isempty,outputs(:,1));
sum(ok)
Portion(bad)
condition
load('OLM21_day9_OML_on_Cell.mat')
outputsOn = outputs; rEstimationsCellOn = rEstimationsCell;
34
g = Group(these{:});
open g
gg = g; gg(ismember(gg(:,end),[1 3]),1) = (gg(ismember(gg(:,end),[1 3]),1) - nanmedian(gg(gg(:,end)==1,1)))./diff(quantile(gg(gg(:,end)==1,1),[0.25 0.75]));
gg(ismember(gg(:,end),[2 4]),1) = (gg(ismember(gg(:,end),[2 4]),1) - nanmedian(gg(gg(:,end)==2,1)))./diff(quantile(gg(gg(:,end)==2,1),[0.25 0.75]));
clf
anovabar(gg(:,1)*100,g(:,end))
kruskalbar(gg(:,1)*100,g(:,end))
nanmedian(gg(gg(:,end)==4,1))
[nanmedian(gg(gg(:,end)==3,1)) nanmedian(gg(gg(:,end)==4,1))]
Accumulate(g(:,end),g(:,2)<0.05)./Accumulate(g(:,end),1)
q = Accumulate(g(:,end),g(:,2)<0.05)./Accumulate(g(:,end),1); bar(q(3:4)./q(1:2))
q = Accumulate(g(:,end),g(:,2)<0.05)./Accumulate(g(:,end),1); bar(q(3:4)./q(1:2) - 1)
subplot(1,2,1);
ok = g(:,end)>2;
kruskalbar(gg(ok,1),g(ok,end)-1)
title([nanmedian(gg(gg(:,end)==3,1)) nanmedian(gg(gg(:,end)==4,1))])
title('all bursts');
set(gca,'tickdir','out','xtick',1:2,'xticklabel',{'on','off'},'box','off');
ylabel('score relative to baseline');
subplot(1,2,2);
q = Accumulate(g(:,end),g(:,2)<0.05)./Accumulate(g(:,end),1);
q = q(3:4)./q(1:2); bar(q-1);
title(q);
set(gca,'tickdir','out','xtick',1:2,'xticklabel',{'on','off'},'box','off');
set(gca,'yticklabel',get(gca,'ytick')+1);
EquateScales(1,3)
sum(reOn)
sum(reOff)
corr(reOn,reOff)
sum(reOn==reOff)
bad = cellfun(@isempty,outputs(:,1));
sum(reOn==reOff & ~bad)
find(reOn==reOff & ~bad,1)
size(reOn)
figure; plot(reOn==reOff & ~bad,1)
plot(reOn==reOff & ~bad)
outputsOn
open outputsOff
hold all
plot(reOn)
find(post,1)
plot(reOff)
Portion(pre|post)
clf
plot(cell2mat(outputsOff(:,1)));
hold all
plot(cell2mat(outputsOn(:,1)));
bad = cellfun(@isempty,outputsOn(:,1));
sum(bad)
sum(bad | cellfun(@eq,outputsOff(:,1),outputsOn(:,1)))
sum(bad | cellfun(@eq,outputsOff(:,1),outputsOn(:,1),'UniformOutput',false))
sum(bad | cell2mat(cellfun(@eq,outputsOff(:,1),outputsOn(:,1),'UniformOutput',false)))
cellfun(@eq,outputsOff(:,1),outputsOn(:,1),'UniformOutput',false);
sum(bad | cellfun(@(x) ~isempty(x) && x==1, cellfun(@eq,outputsOff(:,1),outputsOn(:,1),'UniformOutput',false)))
sum(bad)
checkout = find(bad);
open checkout
tic
for i=find(checkout)'
if sum(rID==i & ~empty)>=3
outputs(i,:) = FindReplayScoreCell(rEstimationsCell{i},'circular','off');
end
if rem(i,100)==0, display([num2str(i) '/' num2str(size(r,1))]); end
end
toc
i
for i=checkout'
if sum(rID==i & ~empty)>=3
outputs(i,:) = FindReplayScoreCell(rEstimationsCell{i},'circular','off');
end
if rem(i,100)==0, display([num2str(i) '/' num2str(size(r,1))]); end
end
toc
save(['M:\home\raly\Documents\code\new\tempFiles\' name],'bursts','rEstimationsCell','outputs','intervals','pos2');
outputsOn = outputs; rEstimationsCellOn = rEstimationsCell;
plot(cell2mat(outputsOn(:,1)));
plot(cell2mat(outputsOff(:,1)));
clf
plot(cell2mat(outputsOff(:,1)));
hold on
plot(cell2mat(outputsOn(:,1)));
clf
plot(cell2mat(outputsOn(:,1)) - cell2mat(outputsOff(:,1)));
4
45
i=200;
outputs(i,:) = FindReplayScoreCell(rEstimationsCell{i},'circular','off');
% Then, on
name = [sessionID '_OML_on_Cell.mat'];
[rEstimations] = ReconstructPosition(pos1,spikes,rWindows,'training',stimIntervals);
tic
rEstimationsCell = cell(size(r,1),1);
for i=1:size(r,1)
rEstimationsCell{i,1} = rEstimations(:,rID==i);
end
outputs(i,:) = FindReplayScoreCell(rEstimationsCell{i},'circular','off');
open outputs
bad = cellfun(@isempty,outputsOn(:,1));
sum(bad)
name
bad = (bad | cellfun(@(x) ~isempty(x) && x==1, cellfun(@eq,outputsOff(:,1),outputsOn(:,1),'UniformOutput',false)));
open bad
sum(bad)
plot(bad)
% Then, on
name = [sessionID '_OML_on_Cell.mat'];
[rEstimations] = ReconstructPosition(pos1,spikes,rWindows,'training',stimIntervals);
tic
rEstimationsCell = cell(size(r,1),1);
for i=1:size(r,1)
rEstimationsCell{i,1} = rEstimations(:,rID==i);
end
toc
outputs = cell(size(r,1),13);
parfor i=1:size(r,1)
if sum(rID==i & ~empty)>=3
outputs(i,:) = FindReplayScoreCell(rEstimationsCell{i},'circular','off');
end
if rem(i,100)==0, display([num2str(i) '/' num2str(size(r,1))]); end
end
toc
save(['M:\home\raly\Documents\code\new\tempFiles\' name],'bursts','rEstimationsCell','outputs','intervals','pos2');
outputsOn = outputs; rEstimationsCellOn = rEstimationsCell;
text(1,2,'5')
clf
text(1,2,'5')
help text
text(1,2,'HELLO')
text(0.5,0.5,'HELLO')
text(0.5,0.5,'HELLO','FontSize',12)
subplot(2,4,4);
ns = Accumulate(g(:,end),1);
q = Accumulate(g(:,end),g(:,2)<0.05 & g(:,3)>0);
z = zBinomialComparison(q(3),ns(3),q(1),ns(1));
z(2) = zBinomialComparison(q(4),ns(4),q(2),ns(2));
q = q./ns; q = q(3:4)./q(1:2); bar(q-1);
title(['all bursts: ' num2str(round(q(1)*1000)/1000) ', ' num2str(round(1000*q(2))/1000)]);
set(gca,'tickdir','out','xtick',1:2,'xticklabel',{['on,p=' num2str(z2p(z(1)))],['off,p=' num2str(z2p(z(2)))]},'box','off');
set(gca,'yticklabel',get(gca,'ytick')+1);
ylabel('proportion forward replay (post/baseline)');
text(1,q(1),num2str(round(q(1)*1000)/1000),'FontSize',12)
text(1,q(1)-1,num2str(round(q(1)*1000)/1000),'FontSize',12)
text(1,q(1)-1,num2str(round(q(1)*1000)/1000),'FontSize',12,'HorizontalAlignment','center')
subplot(2,4,4);
ns = Accumulate(g(:,end),1);
q = Accumulate(g(:,end),g(:,2)<0.05 & g(:,3)>0);
z = zBinomialComparison(q(3),ns(3),q(1),ns(1));
z(2) = zBinomialComparison(q(4),ns(4),q(2),ns(2));
q = q./ns; q = q(3:4)./q(1:2);
bar(q-1);
text(1,q(1)-1,num2str(round(q(1)*1000)/1000),'FontSize',12,'HorizontalAlignment','center')
text(1,q(1)-1,num2str(round(q(1)*1000)/1000),'FontSize',12,'HorizontalAlignment','center','FontAngle',90)
text(1,q(1)-1,num2str(round(q(1)*1000)/1000),'FontSize',12,'HorizontalAlignment','center','Rotation',90
text(1,q(1)-1,num2str(round(q(1)*1000)/1000),'FontSize',12,'HorizontalAlignment','center','Rotation',90)
text(1,(q(1)-1)/2,num2str(round(q(1)*1000)/1000),'FontSize',12,'HorizontalAlignment','center','Rotation',90)
text(1,(q(1)-1)/2,num2str(round(q(1)*1000)/1000),'FontSize',12,'HorizontalAlignment','center','Rotation',90,'Color',,'w')
text(1,(q(1)-1)/2,num2str(round(q(1)*1000)/1000),'FontSize',12,'HorizontalAlignment','center','Rotation',90,'Color','w')
text(2,(q(2)-1)/2,num2str(round(q(2)*1000)/1000),'FontSize',12,'HorizontalAlignment','center','Rotation',90,'Color','w')
subplot(2,4,4);
ns = Accumulate(g(:,end),1);
q = Accumulate(g(:,end),g(:,2)<0.05 & g(:,3)>0);
z = zBinomialComparison(q(3),ns(3),q(1),ns(1));
z(2) = zBinomialComparison(q(4),ns(4),q(2),ns(2));
q = q./ns; q = q(3:4)./q(1:2);
bar(q-1);
text(1,(q(1)-1)/2,num2str(round(q(1)*1000)/1000),'FontSize',12,'HorizontalAlignment','center','Rotation',90,'Color','w')
text(2,(q(2)-1)/2,num2str(round(q(2)*1000)/1000),'FontSize',12,'HorizontalAlignment','center','Rotation',90,'Color','w')
title(['all bursts: ' num2str(round(q(1)*1000)/1000) ', ' num2str(round(1000*q(2))/1000)]);
set(gca,'tickdir','out','xtick',1:2,'xticklabel',{['on,p=' num2str(z2p(z(1)))],['off,p=' num2str(z2p(z(2)))]},'box','off');
set(gca,'yticklabel',get(gca,'ytick')+1);
ylabel('proportion forward replay (post/baseline)');
z2p(zBinomialComparison(q(1),ns(3),q(2),ns(4))
zBinomialComparison(q(1),ns(3),q(2),ns(4))
this = g(ismember(g(:,3),[1 3]),2);
open this
this = g(ismember(g(:,end),[1 3]),2);
this = g(ismember(g(:,end),[1 3]),:);
this = g(ismember(g(:,end),[1 3]),[2 4]);
this = g(ismember(g(:,end),[1 3]),[2 4]); this(:,1) = this(:,1)<0.05;
b = bootstrp(100,@(x) mean(x(x(:,2)==2))./mean(x(x(:,2)==1)),this);
open b
this = g(ismember(g(:,end),[1 3]),[2 4]); this(:,1) = this(:,1)<0.05; this(:,2) = (this(:,2)-1)./2+1;
b = bootstrp(100,@(x) mean(x(x(:,2)==2))./mean(x(x(:,2)==1)),this);
this = g(ismember(g(:,end),[2 4]),[2 4]); this(:,1) = this(:,1)<0.05; this(:,2) = (this(:,2))./2;
b(:,2) = bootstrp(100,@(x) mean(x(x(:,2)==2))./mean(x(x(:,2)==1)),this);
figure; hist(b,100);
hist(b,10000);
hist(b,10);
nanmean(b)
z
q
q = Accumulate(g(:,end),g(:,2)<0.05 & g(:,3)>0);
q
[q ns q./ns]
[q ns]
[q ns q./ns];
open ans
q = q./ns; q = q(3:4)./q(1:2);
q
q = Accumulate(g(:,end),g(:,2)<0.05 & g(:,3)>0);
q0 = q;
q = q./ns; q = q(3:4)./q(1:2);
q0
q
[q0 ns q./ns];
[q0 ns q0./ns];
q(1)*(q0(1)/ns(1))*ns(3)
[q0 ns q0./ns];
that4 = q(1)*(q0(1)/ns(1))*ns(4);
that4 = q(1)*(q0(2)/ns(2))*ns(4);
that4 = q(1)*(q0(2)/ns(2));
zBinomialComparison(q(1),ns(3),q(2)/ns(4))
z2p(zBinomialComparison(q(1),ns(3),q(2)/ns(4)))
z2p(zBinomialComparison(q(2),ns(4),q(1)/ns(3)))
zBinomialComparison(q(2),ns(4),q(1)/ns(3))
q
subplot(2,4,4);
ns = Accumulate(g(:,end),1);
q = Accumulate(g(:,end),g(:,2)<0.05 & g(:,3)>0);
z = zBinomialComparison(q(3),ns(3),q(1),ns(1));
z(2) = zBinomialComparison(q(4),ns(4),q(2),ns(2));
p = z2p(zBinomialComparison(q(2),ns(4),q(1)/ns(3)));
p
q
q0
ns
prediction = (q(3)./ns(3))/(q(1)./ns(1))
prediction = (q(3)./ns(3))/(q(1)./ns(1))*(q(2)/ns(2))
q
zBinomialComparison(q(4),ns(4),(q(3)./ns(3))/(q(1)./ns(1))*(q(2)/ns(2)))
zBinomialComparison(q(4),ns(4),predition);
zBinomialComparison(q(4),ns(4),prediction);
thatAll = [224 146 519 365; 2501 4629]'
thatAll = [224 146 519 365; 2501 2501 4629 4629]'
q = Accumulate(g(:,end),g(:,2)<0.05 & g(:,3)>0);
ns = Accumulate(g(:,end),1);
thatNow = [q ns]
1259/2501
2198/4629
2431/4629
thatAll(:,1)./thatAll(:,2)
thatNow(:,1)./thatNow(:,2)
[thatAll(:,1)./thatAll(:,2) thatNow(:,1)./thatNow(:,2)
]
[thatAll(:,1)./thatAll(:,2) thatNow(:,1)./thatNow(:,2)]
[thatAll(:,1)./thatAll(:,2)./thatNow(:,1)./thatNow(:,2)]
[thatAll(:,1)./thatAll(:,2) thatNow(:,1)./thatNow(:,2)]
ans(:,1)./ans(:,2)
pp = cell2mat([outputsOn(:,2) outputsOff(:,2)]);
open pp
Accumulate((pp<0.05)+1)
size(pp)
size(pre)
size(pre(~bad))
bad = cellfun(@isempty,outputs(:,1));
size(pre(~bad))
Accumulate((pp(pre,:)<0.05)+1)
Accumulate((pp(pre(~bad),:)<0.05)+1)
Accumulate((pp(post(~bad),:)<0.05)+1)
Accumulate((pp(post(~bad),:)<0.05)+1)./Accumulate((pp(pre(~bad),:)<0.05)+1)
ss = cell2mat([outputsOn(:,1) outputsOff(:,1)]);
open ss
sum((ss(:,1)<ss(:,2)) & pp(:,1)<pp(:,2))
sum((ss(:,1)<ss(:,2)) & pp(:,1)>pp(:,2))
sum((ss(:,2)<ss(:,1)) & pp(:,2)>pp(:,1))
sum((ss(:,2)<ss(:,1)) & pp(:,2)<pp(:,1))
edit LinearTrackBehavior
speed = LinearVelocity([behavior.timestamps(:) behavior.position.x(:) behavior.position.y(:)],5);
figure; PlotXY(speed);
5
behavior
basepath
load(fullfile(basepath,[basename '.animal.behavior.mat']),'behavior');
basename
behavior
cd(basepath)
34
34
23
pos0 = sortrows([behavior.positionTrials{1}; behavior.positionTrials{2}(:,1) 2-behavior.positionTrials{2}(:,2)]);
pos1 = Restrict(behavior.positionTrials{1},run);
pos2 = Restrict(behavior.positionTrials{2},run); pos2(:,2) = 1-pos2(:,2);
sum(InIntervals(pos2,stim))
sum(InIntervals(pos1,stim))
pos0 = sortrows([behavior.positionTrials{1}; behavior.positionTrials{2}(:,1) 2-behavior.positionTrials{2}(:,2)]);
pos1 = Restrict(behavior.positionTrials{1},run);
pos2 = Restrict(behavior.positionTrials{2},run); pos2(:,2) = 1-pos2(:,2);
PlotXY(pos);
behavior
load('day10.animal.behavior.mat')
load('day9.animal.behavior.mat')
45
1
ylim
ylim([50.5 150.5]);
15/200*6
15/200*4
15/400*4
15/200*8
10/200*8
% Together (merged positions, on is from 0 to 0.5, off is from 0.5 to 1)
name = [sessionID '_OML_Cell.mat'];
[rEstimations] = ReconstructPosition(pos,spikes,rWindows,'training',stimIntervals);
tic
rEstimationsCell = cell(size(r,1),1);
for i=1:size(r,1)
rEstimationsCell{i,1} = rEstimations(:,rID==i);
end
toc
outputs = cell(size(r,1),13);
parfor i=1:size(r,1)
if sum(rID==i & ~empty)>=3
outputs(i,:) = FindReplayScoreCell(rEstimationsCell{i},'circular','off','threshold',10);
end
if rem(i,100)==0, display([num2str(i) '/' num2str(size(r,1))]); end
end
toc
save(['M:\home\raly\Documents\code\new\tempFiles\' name],'bursts','rEstimationsCell','outputs','intervals','pos2');
outputsBoth = outputs; rEstimationsCellBoth = rEstimationsCell;
56
4
clf
hist(cell2mat(outputsBoth(:,3)))
hist(cell2mat(outputsBoth(:,3)),100)
hist(cell2mat(outputsBoth(:,4)),100)
hist(pos(:,2))
hist(pos(:,2),100)
DensityMap(cell2mat(outputsBoth(:,3)),cell2mat(outputsBoth(:,4)))
matrix = [(cell2mat(outputsBoth(:,3)),cell2mat(outputsBoth(:,4)))]
matrix = [(cell2mat(outputsBoth(:,3)),cell2mat(outputsBoth(:,4)))];
matrix = [cell2mat(outputsBoth(:,3)),cell2mat(outputsBoth(:,4));
matrix = [cell2mat(outputsBoth(:,3)),cell2mat(outputsBoth(:,4))];
open matrix
matrix = repmat(matrix,1,5);
matrix(:,2:end-1) = nan;
matrix = matrix'; matrix(isnan(matrix(:))) = interp1(find(~isnan(matrix(:))),matrix(~isnan(matrix(:))),find(isnan(matrix(:))));
matrix = matrix'l
matrix = matrix';
PlotColorMap(matrix);
hist(matrix(:,1))
b = Bin(matrix);
b = Bin(matrix,100);
b = matrix; b(:) = Bin(matrix(:),100);
open b
b = hist(matrix,100);
open h
size(b)
PlotColorMap(b);
climc
clim
clim([0 1000])
clim([0 500])
PlotXY(pos)
plot(pos(:,1))
plot(pos(:,2))
figure; semplot(rEstimations');
% Together (merged positions, on is from 0 to 0.5, off is from 0.5 to 1)
name = [sessionID '_OML_Cell.mat'];
[rEstimations] = ReconstructPosition(pos,spikes,rWindows,'training',run);
clf
set(gca,'yticklabel',get(gca,'ytick')+1);
semplot(rEstimations');
% Together (merged positions, on is from 0 to 0.5, off is from 0.5 to 1)
name = [sessionID '_OML_Cell.mat'];
[rEstimations] = ReconstructPosition(pos,spikes,rWindows,'training',run);
tic
rEstimationsCell = cell(size(r,1),1);
for i=1:size(r,1)
rEstimationsCell{i,1} = rEstimations(:,rID==i);
end
toc
outputs = cell(size(r,1),13);
parfor i=1:size(r,1)
if sum(rID==i & ~empty)>=3
outputs(i,:) = FindReplayScoreCell(rEstimationsCell{i},'circular','off','threshold',10);
end
if rem(i,100)==0, display([num2str(i) '/' num2str(size(r,1))]); end
end
toc
save(['M:\home\raly\Documents\code\new\tempFiles\' name],'bursts','rEstimationsCell','outputs','intervals','pos2');
outputsBoth = outputs; rEstimationsCellBoth = rEstimationsCell;
45
open outputsBoth
duration = diff(r(ok,:),[],2);
max(duration)
this = cell2mat(outputsBoth);
this = cell2mat(outputsBoth(:,3:4));
figure; hist(this);
hist(this,100);
hist(this(pre,:),100);
hist(this(post,:),100);
pre = InIntervals(bursts0,preSleep);
pre = InIntervals(bursts,preSleep);
post = InIntervals(bursts0,postSleep);
post = InIntervals(bursts,postSleep);
hist(this(pre,:),100);
hist(this(post,:),100);
clear all
close all
basepath = 'M:\Data\Can\OLM21\day9';
script_Can_replay
cd M:\home\raly\Documents\
fullfile('M:\home\raly\results\OML',[sessionID '-theta-sequences-On-Off'])
title(strrep([sessionID 'stim OFF'],'_','-'));
conditionNames = {strrep([sessionID ' ON'],'_','-'),strrep([sessionID ' OFF'],'_','-')}
title([session ID ' all bursts: ' num2str(round(nanmedian(gg(gg(:,end)==3,1))*1000)/1000) ', ' num2str(round(1000*nanmedian(gg(gg(:,end)==4,1)))/1000)]);
outputs = outputsBoth;
size(rEstimations,1)
size(r)
size(ok)
sum(ok
sum(ok)
clf
subplot(3,2,5);
ok = g(:,end)>2;
kruskalbar(gg(ok,1),g(ok,end)-2)
title(['better fit (on/off): ' num2str(round(nanmedian(gg(gg(:,end)==3,1))*1000)/1000) ', ' num2str(round(1000*nanmedian(gg(gg(:,end)==4,1)))/1000)]);
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
size(pre)
Portion(reOn)
Portion(reOff)
ns = Accumulate(g(:,end),1);
q = Accumulate(g(:,end),g(:,2)<0.05);
[q ns]
ans(:,1)./ans(:,2)
sum(reOn)
sum(reOff)
sum(~reOn & ~reOff)
bad = cellfun(@isempty,outputs(:,1));
sum(bad)
SaveFig(fullfile('M:\home\raly\results\OML',[sessionID '-Replay-barplots']))
script_Can_replay
close all
%-- 3/22/2022 2:34 PM --%
v1code
45
%-- 3/22/2022 5:39 PM --%
v1code
4
cd N:\V1test\
anovabar(correct,right+1);
dayNumber
2
close all
Hide('none');
5
456
i
length(listOfDays)
%-- 3/23/2022 12:11 PM --%
cd N:\V1test\AO52\
45
576
nDays
I
i
45
basepath
[parentFolder,dayName] = fileparts(basepath)
[~,projectName] = fileparts(parentFolder);
projectName
parentXml = fullfile(parentFolder,[projectName '.xml']);
parentXml
currentXml = fullfile(basepath,[projectName '.xml']);
currentXml
exist(currentXml,'file')
clear all
i
basepath
folders
45
load('digitalIn.events.mat')
sync = [digitalIn.timestampsOn{5} digitalIn.timestampsOff{5}];
open sync
g = Group(sync(:,1),sync(:,2));
figure; RasterPlot(g);
mode(diff(sync,[],2))
load('day53.MergePoints.events.mat')
[in,ses] = InIntervals(sync(:,1),MergePoints.timestamps);
Accumulate(ses)
sync = syn(ses==2,:);
sync = sync(ses==2,:);
open sync
mod(diff(sync(:,1))
mod(diff(sync(:,1)))
mode(diff(sync(:,1)))
1/mode(diff(sync(:,1)))
%-- 3/24/2022 4:30 PM --%
cd N:\V1test\V1Jean\day53
load('day53.MergePoints.events.mat')
load('digitalIn.events.mat')
sync = [digitalIn.timestampsOn{5} digitalIn.timestampsOff{5}];
[in,ses] = InIntervals(sync(:,1),MergePoints.timestamps);
sync = sync(ses==2,:);
video_path = 'N:\V1test\V1Jean\day53\day53_220324_124050_Ymaze\Basler_acA800-510uc__22592284__20220324_124104166.avi')
video_path = 'N:\V1test\V1Jean\day53\day53_220324_124050_Ymaze\Basler_acA800-510uc__22592284__20220324_124104166.avi';
videoObj = VideoReader(video_path);   % get video
numFrames = get(videoObj, 'NumFrames')
32707/32642
32707-32642
65*40
g = Group(digitalIn.timestampsOn{:});
RasterPlot(g);
clf
hist(sync(:,1));
hist(sync(:,1),1000);
hist(sync(:,1),10000);
hist(sync(:,1),100000);
hist(sync(:,1),100);
hist(diff(sync(:,1))
hist(diff(sync(:,1)))
hist(diff(sync(:,1)),100)
sum(diff(sync(:,1))<0.03)
g = Group(sync(:,1),sync(:,2));
RasterPlot(g);
32707-32642
sum(diff(sync(:,1))<0.01)
min(diff(sync(:,1)))
figure; hist(diff(sync(:,1)),100);
length(unique(diff(sync(:,1))))
[u,~,uu] = unique(diff(sync(:,1)));
open u
open uuu
open uu
bar(u,Accumulate(uu))
bar(u,log(Accumulate(uu)))
bar(u,Accumulate(uu)); set(gca,'YScale','log');
sum(diff(sync(:,1))<0.02222)
digitalInp = getDigitalIn('all','fs',session.extracellular.sr);
digitalInp = getDigitalIn('all','fs',20000);
disp('Loading digital channels...');
m = memmapfile(filename,'Format','uint16','writable',false);
digital_word2 = double(m.Data);
clear m
Nchan = 16;
Nchan2 = 17;
if isempty(filename)
filename=dir('digitalIn.dat');
filename = filename.name;
elseif exist('filename','var')
disp(['Using input: ',filename])
else
disp('No digitalIn file indicated...');
end
try [amplifier_channels, notes, aux_input_channels, spike_triggers,...
board_dig_in_channels, supply_voltage_channels, frequency_parameters,board_adc_channels] =...
read_Intan_RHD2000_file_bz;
fs = frequency_parameters.board_dig_in_sample_rate;
catch
disp('File ''info.rhd'' not found. (Type ''help <a href="matlab:help loadAnalog">loadAnalog</a>'' for details) ');
end
disp('Loading digital channels...');
m = memmapfile(filename,'Format','uint16','writable',false);
digital_word2 = double(m.Data);
clear m
Nchan = 16;
Nchan2 = 17;
Nchan
for k = 1:Nchan
tester(:,Nchan2-k) = (digital_word2 - 2^(Nchan-k))>=0;
digital_word2 = digital_word2 - tester(:,Nchan2-k)*2^(Nchan-k);
test = tester(:,Nchan2-k) == 1;
test2 = diff(test);
pulses{Nchan2-k} = find(test2 == 1);
pulses2{Nchan2-k} = find(test2 == -1);
data(k,:) = test;
end
digital_on = pulses;
digital_off = pulses2;
figure; PlotColorMap(tester+1);
45
figure; PlotColorMap(tester(:1:2000,:)+1);
figure; PlotColorMap(tester(1:2000,:)+1);
plot(terster(:,5));
plot(tester(:,5));
clf
345
plot(tester(1:100000,5));
plot(tester(1:1000000,5));
plot(tester(1:10000000,5));
10000000/length(tester)
length(tester)/2
plot(tester((1:10000000)+160535616,5));
plot(tester((1:100000)+160535616,5));
plot(tester((1:10000)+160535616,5));
%-- 3/24/2022 5:25 PM --%
edit SWRpipeline.m
cd N:\V1test\V1Jean\day53
%-- 3/24/2022 5:28 PM --%
edit preprocessSession.m
cd N:\V1test\V1Jean\day53
f = dir('Kilosort*')
f0 = dir('Kilosort*')
isempty(d0)
isempty(f0)
isempty(f)
1
45
1
cd N:\V1test\V1Jean\
45
i=44
basepath = folders{i}
cd(basepath);
f = dir('Kilosort*');
if ~isempty(f),
done(i,1) = true;
error('done');
end
%% Pull meta data
% Get session names
if strcmp(basepath(end),filesep)
basepath = basepath(1:end-1);
end
[~,basename] = fileparts(basepath);
% Get xml file from the parent folder
[parentFolder,dayName] = fileparts(basepath);
[~,projectName] = fileparts(parentFolder);
parentXml = fullfile(parentFolder,[projectName '.xml']);
xmlFile = fullfile(basepath,[projectName '.xml']);
if ~exist(xmlFile,'file'), copyfile(parentXml,xmlFile); end
%         % Get xml file in order
%         xmlFile = checkFile('fileType','.xml','searchSubdirs',true);
%         xmlFile = xmlFile(1);
%         if ~(strcmp(xmlFile.folder,basepath)&&strcmp(xmlFile.name(1:end-4),basename))
%             copyfile([xmlFile.folder,filesep,xmlFile.name],[basepath,filesep,basename,'.xml'])
%         end
% Check info.rhd
% (assumes this will be the same across subsessions)
rhdFile = checkFile('fileType','.rhd','searchSubdirs',true);
rhdFile = rhdFile(1);
if ~(strcmp(rhdFile.folder,basepath)&&strcmp(rhdFile.name(1:end-4),basename))
copyfile([rhdFile.folder,filesep,rhdFile.name],[basepath,filesep,basename,'.rhd'])
end
%% Make SessionInfo
% Manually ID bad channels at this point. automating it would be good
session = sessionTemplate(pwd,'showGUI',false); % show GUI only after concatenating data and getting MergePoints
save([basename '.session.mat'],'session');
% Get xml file from the parent folder
[parentFolder,dayName] = fileparts(basepath);
[~,projectName] = fileparts(parentFolder);
parentXml = fullfile(parentFolder,[projectName '.xml']);
xmlFile = fullfile(basepath,[projectName '.xml'])
~exist(xmlFile,'file')
parentXml
if ~exist(xmlFile,'file'), copyfile(parentXml,xmlFile); end
xmlFile
projectName
xmlFile = fullfile(basepath,[dayName '.xml'])
45
clear all
load('rez.mat')
kilosortFolder = pwd
rezToPhy_KSW(rez);
CleanRez(rez,'savepath',kilosortFolder);;
file = memmapfile('day15.dat','Format','int16','Writable',false);
data = reshape(file.data,64,[]);
plot(data(:,1:20000)');
t = (1:size(data,2))/20000;
bad = t>19956.7;
Portion(bad)
v1code
data(:,bad) = 0;
45
%-- 3/28/2022 5:52 PM --%
cd N:\V1test\V1Jean\day15
basepath = pwd
session = getSession('basepath',basepath);
nChannels = session.extracellular.nChannels;
sf = session.extracellular.srLfp;
dat_path = [basepath,filesep, session.general.name, '.dat'];
lfp = getLFP('all','basepath',basepath);
% Identify Noise intervals
mean_lfp = mean(lfp.data,2);
t = lfp.timestamps;
interval = [-1 1]*10+19974;
restricted = Restrict([double(t) double(mean_lfp)],interval);
[w,wt,wf] = WaveletSpectrogram(restricted);
PlotColorMap(w,'x',wt,'y',wf);
PlotColorMap(zscore(w,[],2),'x',wt,'y',wf);
clim
clim([-1 20]);
clim([-1 10]);
clim([-1 1]);
interval = [-1 1]*20+19974;
interval
restricted = Restrict([double(t) double(mean_lfp)],interval);
[w,wt,wf] = WaveletSpectrogram(restricted);
PlotColorMap(zscore(w,[],2),'x',wt,'y',wf);
clim([-1 1]);
plot(wt,nanmean(w));
[w,wt,wf] = WaveletSpectrogramRaw(restricted);
plot(wf)
max(wf)
plot(zscore(w'))
PlotColorMap(w);
PlotColorMap(zscore(w,[],2));
clim
clim([-1 1]);
PlotColorMap(zscore(Shrink(w,1,1250),[],2));
PlotColorMap(zscore(Shrink(w,1,1250),[],2),'x',Shrink(t,1250));
PlotColorMap(zscore(Shrink(w,1,1250),[],2),'x',Shrink(t,1,1250));
PlotColorMap(zscore(Shrink(w,1,1250),[],2),'x',Shrink(wt,1250,1));
wf(150)
wf(151)
wf(165)
wf(168)
wf(170)
wf(175)
which FilterLFP
filtered = FilterLFP([t mean_lfp],'passband',[300 625]);
filtered = FilterLFP([t mean_lfp],'passband',[300 624]);
[~,ah] = Phase(filtered);
PlotXY(ah);
hi_pow = ah;
thresh = mean(hi_pow) + 9 * std(hi_pow)
hi_pow = ah(:,2);
thresh = mean(hi_pow) + 9 * std(hi_pow)
above_thresh = hi_pow>thresh;
noiseIntsSlo = findIntervals(above_thresh);
noiseIntsFast = noiseIntsSlo*8; %adjust for 8x sampling r
dur(noiseIntsFast)/20000
45
fclose('all');
cd D:\KiloSort
SSD_path = 'D:\KiloSort2'
kilosortFolder = KiloSortWrapper('SSD_path',SSD_path);
kkeyboard
dbcont
load('rez.mat')
CleanRez(rez,'savepath',kilosortFolder);;
Portion(noiseIndsFlat)
noiseIndsFlat(1,:)
[27568,27640]/20000
diff([27568,27640]/20000)
intWindow
hist(noiseIntsFast(:,1),1000);
plot(noiseIntsFast(:,1),1000);
plot(noiseIntsFast(:,1),diff(noiseIntsFast,[],2));
plot(noiseIntsFast(:,1),diff(noiseIntsFast,[],2)/20000);
t = lfp.timestamps;
bad = t>19956.7;
badbad = repelem(bad,64);
open bad
find(bad,1)
find(badbad,1)
find(badbad,1)/24945877
size(badbad)
size(badbad)-size(m.data,1)
size(bad)
badbad = repelem(bad,64*8);
size(badbad)-size(m.data,1)
size(badbad)-size(m.Data,1)
size(badbad,1)-size(m.Data,1)
size(m.Data,1)/length(t)
size(m.Data,1)/length(t)/8
20000/1250
16*64
badbad = repelem(bad,1024);
size(badbad,1)-size(m.Data,1)
noiseIndsFlat = badbad;
noiseIntsFast0 = noiseIndsFlat;
find(t>19956.7,1)
length(t)
noiseIntsSlo0 = noiseIntsSlo;
noiseIntsSlo = [24945877 25032856];
diff(noiseIntsSlo,[],2)
diff(noiseIntsSlo,[],2)/1250
noiseIntsFast = noiseIntsSlo*8; %adjust for 8x sampling rate in dat file [FIXTHIS]
20000/1250
noiseIntsSlo = findIntervals(above_thresh);
noiseIntsSlo = sortrows([noiseIntsSlo;24945877 25032856]);
noiseIntsSlo = ConsolidateIntervals(noiseIntsSlo);
noiseIntsFast = noiseIntsSlo*16; %adjust for 8x sampling rate in dat file [FIXTHIS]
45
fclose('all');
56
kilosortFolder = KiloSortWrapper('SSD_path',SSD_path);
load(fullfile(kilosortFolder,'rez.mat'),'rez');
CleanRez(rez,'savepath',kilosortFolder);
basepath
cd(basepath);
session = sessionTemplate(pwd,'showGUI',true,'remove_folder_date',true);
f = dir('Kilosort*');
spikes = loadSpikes('session',session,'clusteringpath',[f.folder filesep f.name]);
cell_metrics = ProcessCellMetrics('session',session,'manualAdjustMonoSyn',false);
channel_mapping
close all
clear all
session = sessionTemplate(pwd,'showGUI',true,'remove_folder_date',true);
f = dir('Kilosort*');
spikes = loadSpikes('session',session,'clusteringpath',[f.folder filesep f.name]);
cell_metrics = ProcessCellMetrics('session',session,'manualAdjustMonoSyn',false);
channel_mapping
close all
session = sessionTemplate(pwd,'showGUI',true,'remove_folder_date',true);
f = dir('Kilosort*');
spikes = loadSpikes('session',session,'clusteringpath',[f.folder filesep f.name]);
cell_metrics = ProcessCellMetrics('session',session,'manualAdjustMonoSyn',false);
channel_mapping
close all
clear all
session = sessionTemplate(pwd,'showGUI',true,'remove_folder_date',true);
f = dir('Kilosort*');
spikes = loadSpikes('session',session,'clusteringpath',[f.folder filesep f.name]);
cell_metrics = ProcessCellMetrics('session',session,'manualAdjustMonoSyn',false);
channel_mapping
close all
clear all
session = sessionTemplate(pwd,'showGUI',true,'remove_folder_date',true);
f = dir('Kilosort*');
spikes = loadSpikes('session',session,'clusteringpath',[f.folder filesep f.name]);
cell_metrics = ProcessCellMetrics('session',session,'manualAdjustMonoSyn',false);
channel_mapping
close all
concatenateDats(pwd,0,0);
session = sessionTemplate(pwd,'showGUI',true,'remove_folder_date',true);
load('day17.session.mat')
session.extracellular
concatenateDats(pwd,0,0);
%Convert from number of samples to recording time of start/ends
for ff = 1:length(datpaths.time)
f = fopen(datpaths.time{ff},'r');
% Determine total number of samples in file
fileStart = ftell(f);
%Read the first time point
firsttimepoint = fread(f,1,'int32');
status = fseek(f,-4,'eof'); %int32 = 4 bytes
lasttimepoint = fread(f,1,'int32');
fileStop = ftell(f);
firstlasttimepoints(ff,:) = [firsttimepoint lasttimepoint];
numsamples(ff) = fileStop./4;
if ff==1
transitiontimes_samp = firstlasttimepoints(ff,:);
else
transitiontimes_samp(ff,:) = firstlasttimepoints(ff,:)+transitiontimes_samp(ff-1,2)+1;
end
end
transitiontimes_samp
transitiontimes_sec = transitiontimes_samp./20000;
transitiontimes_sec(end)/3600
eventsfilename = fullfile(basepath,[basename,'.MergePoints.events.mat']);
MergePoints.timestamps = transitiontimes_sec;
MergePoints.timestamps_samples = transitiontimes_samp;
MergePoints.firstlasttimpoints_samples = firstlasttimepoints;
recordingnames = cellfun(@(x) strrep('N:\V1test\V1Jean\day17\',''), datpaths.time);
recordingnames = cellfun(@(x) strrep(x,'N:\V1test\V1Jean\day17\',''), datpaths.time);
recordingnames = cellfun(@(x) strrep(x,'N:\V1test\V1Jean\day17\',''), datpaths.time, 'UniformOutput', false);
recordingnames = cellfun(@(x) strrep(strrep(x,'N:\V1test\V1Jean\day17\',''),x,'\time.dat',''), datpaths.time, 'UniformOutput', false);
recordingnames = cellfun(@(x) strrep(strrep(x,'N:\V1test\V1Jean\day17\',''),'\time.dat',''), datpaths.time, 'UniformOutput', false);
MergePoints.foldernames = recordingnames;
datpaths.amplifier = cellfun(@(x) strrep(x,'time.dat','amplifier.dat'), datpaths.time, 'UniformOutput', false);
datpaths.amplifier = cellfun(@(x) strrep(x,'time.dat','amplifier.dat'), datpaths.time, 'UniformOutput', false);
datpaths.analogin = cellfun(@(x) strrep(x,'time.dat','analogin.dat'), datpaths.time, 'UniformOutput', false);
datpaths.auxiliary = cellfun(@(x) strrep(x,'time.dat','auxiliary.dat'), datpaths.time, 'UniformOutput', false);
datpaths.digitalin = cellfun(@(x) strrep(x,'time.dat','digitalin.dat'), datpaths.time, 'UniformOutput', false);
MergePoints.filesmerged = datpaths;
% Check that size of resultant .dat is equal to the sum of the components
t = dir(newpaths.(otherdattypes{odidx}));
if t.bytes ~= sum(datsizes.(otherdattypes{odidx}))
error(['New ' otherdattypes{odidx} '.dat size not right.  Exiting after .dats converted.  Not deleting'])
deleteoriginaldatsbool = 0;
sizecheck.(otherdattypes{odidx}) = false;
else
sizecheck.(otherdattypes{odidx}) = true;
disp([otherdattypes{odidx} ' concatenated and size checked'])
end
odidx = 1
t = dir(newpaths.(otherdattypes{odidx}));
if t.bytes ~= sum(datsizes.(otherdattypes{odidx}))
error(['New ' otherdattypes{odidx} '.dat size not right.  Exiting after .dats converted.  Not deleting'])
deleteoriginaldatsbool = 0;
sizecheck.(otherdattypes{odidx}) = false;
else
sizecheck.(otherdattypes{odidx}) = true;
disp([otherdattypes{odidx} ' concatenated and size checked'])
end
otherdattypes = {'analogin';'digitalin';'auxiliary';'time';'supply'};
% Check that size of resultant .dat is equal to the sum of the components
t = dir(newpaths.(otherdattypes{odidx}));
if t.bytes ~= sum(datsizes.(otherdattypes{odidx}))
error(['New ' otherdattypes{odidx} '.dat size not right.  Exiting after .dats converted.  Not deleting'])
deleteoriginaldatsbool = 0;
sizecheck.(otherdattypes{odidx}) = false;
else
sizecheck.(otherdattypes{odidx}) = true;
disp([otherdattypes{odidx} ' concatenated and size checked'])
end
for odidx = 1:length(otherdattypes)
%eval(['new' otherdattypes{odidx} 'path = fullfile(basepath,''' otherdattypes{odidx} '.dat'');'])
newpaths.(otherdattypes{odidx}) = fullfile(basepath,[otherdattypes{odidx}, '.dat']);
end
basepath = pwd
for odidx = 1:length(otherdattypes)
%eval(['new' otherdattypes{odidx} 'path = fullfile(basepath,''' otherdattypes{odidx} '.dat'');'])
newpaths.(otherdattypes{odidx}) = fullfile(basepath,[otherdattypes{odidx}, '.dat']);
end
for odidx = 1:length(otherdattypes)
% Check that size of resultant .dat is equal to the sum of the components
t = dir(newpaths.(otherdattypes{odidx}));
if t.bytes ~= sum(datsizes.(otherdattypes{odidx}))
error(['New ' otherdattypes{odidx} '.dat size not right.  Exiting after .dats converted.  Not deleting'])
deleteoriginaldatsbool = 0;
sizecheck.(otherdattypes{odidx}) = false;
else
sizecheck.(otherdattypes{odidx}) = true;
disp([otherdattypes{odidx} ' concatenated and size checked'])
end
end
odidx
t
d2 = dir(datpaths(otherdattypes{odidx}))
odidx=1
rcount = 1
d2 = dir(datpaths(otherdattypes{odidx}))
datpaths(otherdattypes{odidx})
otherdattypes
odidx
otherdattypes{odidx}
d2 = dir(datpaths.(otherdattypes{odidx}))
datpaths.(otherdattypes{odidx})
odidx
otherdattypes
(otherdattypes{odidx})
d2 = dir(datpaths.(otherdattypes{odidx})(rcount))
datpaths.(otherdattypes{odidx})(rcount)
datpaths.(otherdattypes{odidx}){rcount}
d2 = dir(datpaths.(otherdattypes{odidx}){rcount})
for odidx=1:length(otherdattypes)
for rcount = 1:8
d2 = dir(datpaths.(otherdattypes{odidx}){rcount});
datsizes.(otherdattypes{odidx})(rcount) = d2(1).bytes;
end
end
otherdattypes
MergePoints.timestamps = transitiontimes_sec;
MergePoints.timestamps_samples = transitiontimes_samp;
MergePoints.firstlasttimpoints_samples = firstlasttimepoints;
MergePoints.foldernames = recordingnames;
MergePoints.filesmerged = datpaths;
MergePoints.filesizes = datsizes;
MergePoints.sizecheck = sizecheck;
concatenateDats
for odidx = 1:length(otherdattypes)
t = dir(newpaths.(otherdattypes{odidx}));
if t.bytes ~= sum(datsizes.(otherdattypes{odidx}))
error(['New ' otherdattypes{odidx} '.dat size not right.  Exiting after .dats converted.  Not deleting'])
deleteoriginaldatsbool = 0;
sizecheck.(otherdattypes{odidx}) = false;
else
sizecheck.(otherdattypes{odidx}) = true;
disp([otherdattypes{odidx} ' concatenated and size checked'])
end
end
odidx
MergePoints.sizecheck = sizecheck;
sizecheck
MergePoints.detectorinfo.detectiondate = datestr(now,'yyyy-mm-dd');
MergePoints.detectorinfo.detectorname = 'Manual implementation of bz_ConcatenateDats as it crashed on me before producing the metadata';
eventsfilename = fullfile(basepath,[basename,'.MergePoints.events.mat']);
save(eventsfilename,'MergePoints');
eventsfilename = fullfile(basepath,[basename,'.MergePoints.events.mat']);
basename = basenameFromBasepath(basepath);
eventsfilename = fullfile(basepath,[basename,'.MergePoints.events.mat']);
save(eventsfilename,'MergePoints');
session = sessionTemplate(pwd,'showGUI',true,'remove_folder_date',true);
cd(basepath);
session = sessionTemplate(pwd,'showGUI',true,'remove_folder_date',true);
f = dir('Kilosort*');
spikes = loadSpikes('session',session,'clusteringpath',[f.folder filesep f.name]);
cell_metrics = ProcessCellMetrics('session',session,'manualAdjustMonoSyn',false);
channel_mapping
close all
45
raly
%-- 3/29/2022 6:00 PM --%
v1
edit preprocessV1
startup
addpath 'M:\home\raly\Documents\code\new'
startup
edit preprocessV1
session = sessionTemplate(pwd,'showGUI',true,'remove_folder_date',true);
f = dir('Kilosort*');
spikes = loadSpikes('session',session,'clusteringpath',[f.folder filesep f.name]);
cell_metrics = ProcessCellMetrics('session',session,'manualAdjustMonoSyn',false);
channel_mapping
close all
load(fullfile(kilosortFolder,'rez.mat'),'rez');
kilosortFolder = pwd
load(fullfile(kilosortFolder,'rez.mat'),'rez');
CleanRez(rez,'savepath',kilosortFolder);
figure
[13980 14100]
NoiseRemoval(pwd)
load('noiseIntervals.events.mat')
datSamplingRate = 20000
lfpSamplingRate = 1250
noiseIntervalIndices = round(noiseIntervals*datSamplingRate);
datFile = 'N:\V1test\AO52\day5\day5.dat';
m = memmapfile(datFile, 'Format','int16','Writable',true);
nSamples = round(length(m.Data)/nChannels);
nChannels=64
nSamples = round(length(m.Data)/nChannels);
noiseIntervalIndices(noiseIntervalIndices<2) = 2; noiseIntervalIndices(noiseIntervalIndices>nSamples-1) = nSamples-1;
i = 1
badTimeIndices = linspaceVector(noiseIntervalIndices(:,1),noiseIntervalIndices(:,2));
goodTimeIndices = sort([noiseIntervalIndices(:,1)-1; noiseIntervalIndices(:,2)+1]);
badIndices = sub2ind([nChannels,nSamples],i*ones(size(badTimeIndices)),badTimeIndices);
goodIndices = sub2ind([nChannels,nSamples],i*ones(size(goodTimeIndices)),goodTimeIndices);
goodValues = m.Data(goodIndices);
interpolated = interp1(goodTimeIndices,double(goodValues),badTimeIndices);
45
for i = 1:nChannels
badTimeIndices = linspaceVector(noiseIntervalIndices(:,1),noiseIntervalIndices(:,2));
goodTimeIndices = sort([noiseIntervalIndices(:,1)-1; noiseIntervalIndices(:,2)+1]);
badIndices = sub2ind([nChannels,nSamples],i*ones(size(badTimeIndices)),badTimeIndices);
goodIndices = sub2ind([nChannels,nSamples],i*ones(size(goodTimeIndices)),goodTimeIndices);
goodValues = m.Data(goodIndices);
interpolated = interp1(goodTimeIndices,double(goodValues),badTimeIndices);
m.Data(badIndices) = int16(interpolated);
i
datestr((datenum(clock)))
end
kilosortFolder = KiloSortWrapper('SSD_path',SSD_path);
load(fullfile(kilosortFolder,'rez.mat'),'rez');
CleanRez(rez,'savepath',kilosortFolder);
SSD_path = 'D:\KiloSort'
kilosortFolder = KiloSortWrapper('SSD_path',SSD_path);
load(fullfile(kilosortFolder,'rez.mat'),'rez');
CleanRez(rez,'savepath',kilosortFolder);
%-- 3/29/2022 8:16 PM --%
cd N:\V1test\V1Jean\day12
load('day12.EMGFromLFP.LFP.mat')
plot(EMGFromLFP.timestamps,EMGFromLFP.data);
load('day12.EMGFromLFP.LFP.mat')
getEMGFromLFP(basepath,'rejectchannels',1+[22 13 15 19 47 51 49 57 55 43 54 39 50 52 41 48 37],'overwrite',true); % takes lfp in base 0
getEMGFromLFP(pwd,'rejectchannels',1+[22 13 15 19 47 51 49 57 55 43 54 39 50 52 41 48 37],'overwrite',true); % takes lfp in base 0
[datasetfolder,recordingname] = fileparts(basepath);
matfilename = fullfile(basepath,[recordingname,'.EMGFromLFP.LFP.mat']);
p = inputParser;
addParameter(p,'restrict',[0 inf],@isnumeric)
addParameter(p,'specialChannels',[],@isnumeric)
addParameter(p,'rejectChannels',[],@isnumeric)
addParameter(p,'restrictChannels',[],@isnumeric)
addParameter(p,'saveMat',true,@islogical)
addParameter(p,'saveLocation','',@isstr)
addParameter(p,'overwrite',false,@islogical)
addParameter(p,'samplingFrequency',2,@isnumeric)
addParameter(p,'noPrompts',false,@islogical);
addParameter(p,'fromDat',false,@islogical);
addParameter(p,'chInfo',[],@isstruct);
parse(p,varargin{:})
restrict = p.Results.restrict;
specialChannels = p.Results.specialChannels;
rejectChannels = p.Results.rejectChannels;
restrictChannels = p.Results.restrictChannels;
saveMat = p.Results.saveMat;
overwrite = p.Results.overwrite;
samplingFrequency = p.Results.samplingFrequency;
noPrompts = p.Results.noPrompts;
fromDat = p.Results.fromDat;
chInfo = p.Results.chInfo;
if ~isempty(p.Results.saveLocation)
matfilename = fullfile(p.Results.saveLocation,[recordingname,'.EMGFromLFP.LFP.mat']);
end
exist(matfilename,'file') && ~overwrite
load(fullfile(basepath,[recordingname,'.session.mat']))
nChannels = session.extracellular.nChannels;
SpkGrps = session.extracellular.spikeGroups.channels;
Fs = session.extracellular.srLfp;
lfpFile = checkFile('basepath',basepath,'fileTypes',{'.lfp','.eeg'});
lfpFile = [basepath filesep lfpFile(1).name];
binScootS = 1 ./ samplingFrequency;
binScootSamps = round(Fs*binScootS); % must be integer, or error on line 190
corrChunkSz = 20; %for batch-processed correlations
isempty(restrictChannels)
xcorr_chs = restrictChannels;
~isempty(restrictChannels)
usablechannels = [];
spkgrpstouse = [];
for gidx = 1:length(SpkGrps)
usableshankchannels{gidx} = setdiff(SpkGrps{gidx},rejectChannels);
usablechannels = cat(2,usablechannels,usableshankchannels{gidx});
if ~isempty(usableshankchannels{gidx})
spkgrpstouse = cat(2,spkgrpstouse,gidx);
end
end
% check for good/bad shanks and update here
% spkgrpstouse = unique(cat(1,spkgrpstouse,specialshanks)); % this is redundant with taking all shanks.
% get list of channels (1 from each good spike group)
xcorr_chs = [];
for gidx=1:length(usableshankchannels)
%Remove rejectChannels
%     usableshankchannels = setdiff(SpkGrps(spkgrpstouse(i)).Channels,rejectChannels);
%grab random channel on each shank
if ~isempty(usableshankchannels{gidx})
randChfromShank = usableshankchannels{gidx}(randi(length(usableshankchannels{gidx})));
xcorr_chs = [xcorr_chs,randChfromShank];
end
end
xcorr_chs = unique([xcorr_chs,specialChannels]);
fromDat
lfp = LoadBinary(lfpFile ,'nChannels',nChannels,'channels',xcorr_chs,...
'start',restrict(1),'duration',diff(restrict),'frequency',Fs); %read and convert to mV
1
figure
maxfreqband = floor(max([625 Fs/2]));
% xcorr_freqband = [275 300 600 625]; % Hz
xcorr_freqband = [275 300 maxfreqband-25 maxfreqband]; % Hz
lfp = filtsig_in(lfp, Fs, xcorr_freqband);
sum(isnan(lfp))
xcorr_window_samps = round(binScootS*Fs);
xcorr_window_inds = -xcorr_window_samps:xcorr_window_samps;%+- that number of ms in samples
% new version... batches of correlation calculated at once
timestamps = (1+xcorr_window_inds(end)):binScootSamps:(size(lfp,1)-xcorr_window_inds(end));
numbins = length(timestamps);
EMGCorr = zeros(numbins, 1);
% tic
counter = 1;
xcorr_chs
xcorr_chs = [];
for gidx=1:length(usableshankchannels)
%Remove rejectChannels
%     usableshankchannels = setdiff(SpkGrps(spkgrpstouse(i)).Channels,rejectChannels);
%grab random channel on each shank
if ~isempty(usableshankchannels{gidx})
randChfromShank = usableshankchannels{gidx}(randi(length(usableshankchannels{gidx})));
xcorr_chs = [xcorr_chs,randChfromShank];
end
end
gidx
length(usableshankchannels)
for gidx=1:length(usableshankchannels{1})
%Remove rejectChannels
%     usableshankchannels = setdiff(SpkGrps(spkgrpstouse(i)).Channels,rejectChannels);
%grab random channel on each shank
if ~isempty(usableshankchannels{gidx})
randChfromShank = usableshankchannels{gidx}(randi(length(usableshankchannels{gidx})));
xcorr_chs = [xcorr_chs,randChfromShank];
end
end
xcorr_chs = unique([xcorr_chs,specialChannels]);
xcorr_chs = [usableshankchannels{1}(2:10:end)];
45
clear all
58
% First, remove points that occur in noisy lfp periods
tl = (1:length(lfp))'/1250; % timestamps for lfp
t = featureTs/1250; % timestamps for candidate ripple events
bad = false(size(t));
try % optionally, remove periods of noisy LFP (large deflections)
[clean,~,badIntervals] = CleanLFP([tl,double(lfp(:,2))],'thresholds',[6 Inf],'manual',true);
bad = bad | InIntervals(t,badIntervals);
end
clf; plot(t,d); hold all; plot(xlim,ones(1,2)*threshold2,'r--'); ylabel('lfp derivative (z-units'); legend('derivative','threshold2');
disp('Change "threshold2" manually.')
PlotHVLines([-1 1],'h','r--');
PlotHVLines([-1 1]*5,'h','r--');
PlotHVLines([-1 1]*3,'h','r--');
figure; plot(t,z);
x = xlim
xlim(x)
clf; plot(t,d); hold all; plot(xlim,ones(1,2)*threshold2,'r--'); ylabel('lfp derivative (z-units'); legend('derivative','threshold2');
thresholds = [6 3];
clf; plot(t,z); hold all; plot(xlim,ones(1,2)*threshold1,'r--'); ylabel('lfp derivative (z-units'); legend('derivative','threshold2');
% parameters
threshold1 = thresholds(1); % in sigmas deviating from the mean
aroundArtefact1 = aroundArtefact(1); % 2, Big and long artefacts
threshold2 = thresholds(2); % for derivative of z-scored signal
aroundArtefact2 = aroundArtefact(2); % 0.1 Very fast fluctuations (short time scale)
threshold2
threshold1
dbcont
Portion(bad)
idx1 = idx1 & ~bad; idx2 = idx2 & ~bad; % remove bad points from "idx1" and "idx2"
figure(1); % plot a clean figure without these points and check indivudual ripples
clf
plot(swDiffAll(idx1 | idx2),ripPowerAll(idx1 | idx2),'k.','markersize',1); hold on
matrix = [swDiffAll ripPowerAll];
z = bsxfun(@rdivide,bsxfun(@minus,matrix,mean(matrix(~bad,:))),std(matrix(~bad,:)));
scores = nan(size(ripPowerAll,1),1);
while true % click "Ctrl+C" to exit when you feel you're done
figure(1);
colors = Bright(1000);
[x,y,button] = ginput(1);
[~,i] = min(abs(x-swDiffAll.*(1-bad))+abs(y-ripPowerAll.*(1-bad)));
t = featureTs/1250;
tl = (1:length(lfp))'/1250;
figure(2);
clf; interval = t(i) + [-1 1]*0.5; in = tl>interval(1) & tl<interval(end);
subplot(3,1,1);
plot(tl(in) - t(i),lfp(in,:)); title(num2str(t(i)));
hold on; plot(bsxfun(@times,ones(1,2),(Restrict(t,interval) - t(i))),ylim,'k--');
legend('ripple channel','SW channel','noise channel');
subplot(3,1,2);
plot(tl(in) - t(i),double(swDiff(in)));
hold on; plot(bsxfun(@times,ones(1,2),(Restrict(t,interval) - t(i))),ylim,'k--');
subplot(3,1,3);
plot(tl(in) - t(i),double(ripPower0(in)));
hold on; plot(bsxfun(@times,ones(1,2),(Restrict(t,interval) - t(i))),ylim,'k--');
%     x = input( prompt )
commandwindow % select the command window
score = str2double(input('good(1)/bad(0)?','s'));
scores(i) = score; % save scores to optionally save
figure(1);
xlims = xlim; ylims = ylim; clf
% extrapolate the score of each ripple based on the score of the nearest scored ripple
scored = find(~isnan(scores));
distances = sqrt(bsxfun(@minus,z(:,1),z(scored,1)').^2 + bsxfun(@minus,z(:,2),z(~isnan(scores),2)').^2);
[~,nearestNeighbour] = min(distances,[],2);
estimated = scores(scored(nearestNeighbour));
selected = estimated>0.1;
idx1 = selected & ~bad; % final ripples
idx2 = ~selected & ~bad; % non-ripples (yet free from noise as well)
scatter(matrix(~bad,1),matrix(~bad,2),1,estimated(~bad),'filled'); colormap(Bright); set(gca,'CLim',[0 1])
xlim(xlims); ylim(ylims);
hold on;
scatter(swDiffAll(scored),ripPowerAll(scored),20,scores(scored)); scatter(swDiffAll(scored),ripPowerAll(scored),15,scores(scored));
set(get(colorbar,'YLabel'),'String','Estimated ripple score');
xlabel('Sharp wave depth'); ylabel('Ripple power');
title({['Manual scoring for ' basepath],['called with channels: ' num2str(Channels) ' (ripple, sharp wave, and (optionally) noise)']});
if rem(i,10)==0
display('saving...');
save(fullfile(basepath,'DetectSWR_manual_scoring.mat'),'t','swDiffAll','ripPowerAll','bad','scores');
end
end
0
[clean,~,badIntervals] = CleanLFP([tl,double(lfp(:,2))],'thresholds',[6 2]);
bad = bad | InIntervals(t,badIntervals);
Portion(bad)
idx1 = idx1 & ~bad; idx2 = idx2 & ~bad; % remove bad points from "idx1" and "idx2"
figure(1); % plot a clean figure without these points and check indivudual ripples
clf
plot(swDiffAll(idx1 | idx2),ripPowerAll(idx1 | idx2),'k.','markersize',1); hold on
matrix = [swDiffAll ripPowerAll];
z = bsxfun(@rdivide,bsxfun(@minus,matrix,mean(matrix(~bad,:))),std(matrix(~bad,:)));
scores = nan(size(ripPowerAll,1),1);
while true % click "Ctrl+C" to exit when you feel you're done
figure(1);
colors = Bright(1000);
[x,y,button] = ginput(1);
[~,i] = min(abs(x-swDiffAll.*(1-bad))+abs(y-ripPowerAll.*(1-bad)));
t = featureTs/1250;
tl = (1:length(lfp))'/1250;
figure(2);
clf; interval = t(i) + [-1 1]*0.5; in = tl>interval(1) & tl<interval(end);
subplot(3,1,1);
plot(tl(in) - t(i),lfp(in,:)); title(num2str(t(i)));
hold on; plot(bsxfun(@times,ones(1,2),(Restrict(t,interval) - t(i))),ylim,'k--');
legend('ripple channel','SW channel','noise channel');
subplot(3,1,2);
plot(tl(in) - t(i),double(swDiff(in)));
hold on; plot(bsxfun(@times,ones(1,2),(Restrict(t,interval) - t(i))),ylim,'k--');
subplot(3,1,3);
plot(tl(in) - t(i),double(ripPower0(in)));
hold on; plot(bsxfun(@times,ones(1,2),(Restrict(t,interval) - t(i))),ylim,'k--');
%     x = input( prompt )
commandwindow % select the command window
score = str2double(input('good(1)/bad(0)?','s'));
scores(i) = score; % save scores to optionally save
figure(1);
xlims = xlim; ylims = ylim; clf
% extrapolate the score of each ripple based on the score of the nearest scored ripple
scored = find(~isnan(scores));
distances = sqrt(bsxfun(@minus,z(:,1),z(scored,1)').^2 + bsxfun(@minus,z(:,2),z(~isnan(scores),2)').^2);
[~,nearestNeighbour] = min(distances,[],2);
estimated = scores(scored(nearestNeighbour));
selected = estimated>0.1;
idx1 = selected & ~bad; % final ripples
idx2 = ~selected & ~bad; % non-ripples (yet free from noise as well)
scatter(matrix(~bad,1),matrix(~bad,2),1,estimated(~bad),'filled'); colormap(Bright); set(gca,'CLim',[0 1])
xlim(xlims); ylim(ylims);
hold on;
scatter(swDiffAll(scored),ripPowerAll(scored),20,scores(scored)); scatter(swDiffAll(scored),ripPowerAll(scored),15,scores(scored));
set(get(colorbar,'YLabel'),'String','Estimated ripple score');
xlabel('Sharp wave depth'); ylabel('Ripple power');
title({['Manual scoring for ' basepath],['called with channels: ' num2str(Channels) ' (ripple, sharp wave, and (optionally) noise)']});
if rem(i,10)==0
display('saving...');
save(fullfile(basepath,'DetectSWR_manual_scoring.mat'),'t','swDiffAll','ripPowerAll','bad','scores');
end
end
bad = false(size(t));
[clean,~,badIntervals] = CleanLFP([tl,double(lfp(:,2))],'thresholds',[6 3],'manual',true,'aroundArtefact',[2 0.5]);
bad = bad | InIntervals(t,badIntervals);
dbcont
Portion(bad)
[clean,~,badIntervals] = CleanLFP([tl,double(lfp(:,2))],'thresholds',[6 3],'manual',true,'aroundArtefact',[2 0.5]);
dbcont
bad = bad | InIntervals(t,badIntervals);
Portion(bad)
[clean,~,badIntervals] = CleanLFP([tl,double(lfp(:,1))],'thresholds',[6 3]);
bad = bad | InIntervals(t,badIntervals);
Portion(bad)
idx1 = idx1 & ~bad; idx2 = idx2 & ~bad; % remove bad points from "idx1" and "idx2"
figure(1); % plot a clean figure without these points and check indivudual ripples
clf
plot(swDiffAll(idx1 | idx2),ripPowerAll(idx1 | idx2),'k.','markersize',1); hold on
matrix = [swDiffAll ripPowerAll];
z = bsxfun(@rdivide,bsxfun(@minus,matrix,mean(matrix(~bad,:))),std(matrix(~bad,:)));
scores = nan(size(ripPowerAll,1),1);
while true % click "Ctrl+C" to exit when you feel you're done
figure(1);
colors = Bright(1000);
[x,y,button] = ginput(1);
[~,i] = min(abs(x-swDiffAll.*(1-bad))+abs(y-ripPowerAll.*(1-bad)));
t = featureTs/1250;
tl = (1:length(lfp))'/1250;
figure(2);
clf; interval = t(i) + [-1 1]*0.5; in = tl>interval(1) & tl<interval(end);
subplot(3,1,1);
plot(tl(in) - t(i),lfp(in,:)); title(num2str(t(i)));
hold on; plot(bsxfun(@times,ones(1,2),(Restrict(t,interval) - t(i))),ylim,'k--');
legend('ripple channel','SW channel','noise channel');
subplot(3,1,2);
plot(tl(in) - t(i),double(swDiff(in)));
hold on; plot(bsxfun(@times,ones(1,2),(Restrict(t,interval) - t(i))),ylim,'k--');
subplot(3,1,3);
plot(tl(in) - t(i),double(ripPower0(in)));
hold on; plot(bsxfun(@times,ones(1,2),(Restrict(t,interval) - t(i))),ylim,'k--');
%     x = input( prompt )
commandwindow % select the command window
score = str2double(input('good(1)/bad(0)?','s'));
scores(i) = score; % save scores to optionally save
figure(1);
xlims = xlim; ylims = ylim; clf
% extrapolate the score of each ripple based on the score of the nearest scored ripple
scored = find(~isnan(scores));
distances = sqrt(bsxfun(@minus,z(:,1),z(scored,1)').^2 + bsxfun(@minus,z(:,2),z(~isnan(scores),2)').^2);
[~,nearestNeighbour] = min(distances,[],2);
estimated = scores(scored(nearestNeighbour));
selected = estimated>0.1;
idx1 = selected & ~bad; % final ripples
idx2 = ~selected & ~bad; % non-ripples (yet free from noise as well)
scatter(matrix(~bad,1),matrix(~bad,2),1,estimated(~bad),'filled'); colormap(Bright); set(gca,'CLim',[0 1])
xlim(xlims); ylim(ylims);
hold on;
scatter(swDiffAll(scored),ripPowerAll(scored),20,scores(scored)); scatter(swDiffAll(scored),ripPowerAll(scored),15,scores(scored));
set(get(colorbar,'YLabel'),'String','Estimated ripple score');
xlabel('Sharp wave depth'); ylabel('Ripple power');
title({['Manual scoring for ' basepath],['called with channels: ' num2str(Channels) ' (ripple, sharp wave, and (optionally) noise)']});
if rem(i,10)==0
display('saving...');
save(fullfile(basepath,'DetectSWR_manual_scoring.mat'),'t','swDiffAll','ripPowerAll','bad','scores');
end
end
[clean,~,badIntervals] = CleanLFP([tl,double(lfp(:,1))],'thresholds',[6 2]);
bad = bad | InIntervals(t,badIntervals);
idx1 = idx1 & ~bad; idx2 = idx2 & ~bad; % remove bad points from "idx1" and "idx2"
figure(1); % plot a clean figure without these points and check indivudual ripples
clf
plot(swDiffAll(idx1 | idx2),ripPowerAll(idx1 | idx2),'k.','markersize',1); hold on
matrix = [swDiffAll ripPowerAll];
z = bsxfun(@rdivide,bsxfun(@minus,matrix,mean(matrix(~bad,:))),std(matrix(~bad,:)));
scores = nan(size(ripPowerAll,1),1);
while true % click "Ctrl+C" to exit when you feel you're done
figure(1);
colors = Bright(1000);
[x,y,button] = ginput(1);
[~,i] = min(abs(x-swDiffAll.*(1-bad))+abs(y-ripPowerAll.*(1-bad)));
t = featureTs/1250;
tl = (1:length(lfp))'/1250;
figure(2);
clf; interval = t(i) + [-1 1]*0.5; in = tl>interval(1) & tl<interval(end);
subplot(3,1,1);
plot(tl(in) - t(i),lfp(in,:)); title(num2str(t(i)));
hold on; plot(bsxfun(@times,ones(1,2),(Restrict(t,interval) - t(i))),ylim,'k--');
legend('ripple channel','SW channel','noise channel');
subplot(3,1,2);
plot(tl(in) - t(i),double(swDiff(in)));
hold on; plot(bsxfun(@times,ones(1,2),(Restrict(t,interval) - t(i))),ylim,'k--');
subplot(3,1,3);
plot(tl(in) - t(i),double(ripPower0(in)));
hold on; plot(bsxfun(@times,ones(1,2),(Restrict(t,interval) - t(i))),ylim,'k--');
%     x = input( prompt )
commandwindow % select the command window
score = str2double(input('good(1)/bad(0)?','s'));
scores(i) = score; % save scores to optionally save
figure(1);
xlims = xlim; ylims = ylim; clf
% extrapolate the score of each ripple based on the score of the nearest scored ripple
scored = find(~isnan(scores));
distances = sqrt(bsxfun(@minus,z(:,1),z(scored,1)').^2 + bsxfun(@minus,z(:,2),z(~isnan(scores),2)').^2);
[~,nearestNeighbour] = min(distances,[],2);
estimated = scores(scored(nearestNeighbour));
selected = estimated>0.1;
idx1 = selected & ~bad; % final ripples
idx2 = ~selected & ~bad; % non-ripples (yet free from noise as well)
scatter(matrix(~bad,1),matrix(~bad,2),1,estimated(~bad),'filled'); colormap(Bright); set(gca,'CLim',[0 1])
xlim(xlims); ylim(ylims);
hold on;
scatter(swDiffAll(scored),ripPowerAll(scored),20,scores(scored)); scatter(swDiffAll(scored),ripPowerAll(scored),15,scores(scored));
set(get(colorbar,'YLabel'),'String','Estimated ripple score');
xlabel('Sharp wave depth'); ylabel('Ripple power');
title({['Manual scoring for ' basepath],['called with channels: ' num2str(Channels) ' (ripple, sharp wave, and (optionally) noise)']});
if rem(i,10)==0
display('saving...');
save(fullfile(basepath,'DetectSWR_manual_scoring.mat'),'t','swDiffAll','ripPowerAll','bad','scores');
end
end
plot(lf,zscore(lfp(:,1)));
plot(tl,zscore(lfp(:,1)));
plot(tl,zscore(double(lfp(:,1))));
figure; plot(tl,zscore(double(lfp(:,1))));
xlim([-1 1] + 16284.42);
plot(tl(2:end),diff(zscore(double(lfp(:,1)))));;
xlim([-1 1] + 16284.42);
figure; plot(tl,zscore(double(lfp(:,1))));
xlim([-1 1] + 16284.42);
[clean,~,badIntervals] = CleanLFP([tl,double(lfp(:,1))],'thresholds',[6 1]);
bad = bad | InIntervals(t,badIntervals);
Portion(bad)
xlim([-1 1] + 16284.42);
clf
plot(tl(2:end),diff(zscore(double(lfp(:,1)))));;
xlim([-1 1] + 16284.42);
PlotIntervals(Restrict(badIntervals,xlim),'color','r');
idx1 = idx1 & ~bad; idx2 = idx2 & ~bad; % remove bad points from "idx1" and "idx2"
figure(1); % plot a clean figure without these points and check indivudual ripples
clf
plot(swDiffAll(idx1 | idx2),ripPowerAll(idx1 | idx2),'k.','markersize',1); hold on
matrix = [swDiffAll ripPowerAll];
z = bsxfun(@rdivide,bsxfun(@minus,matrix,mean(matrix(~bad,:))),std(matrix(~bad,:)));
scores = nan(size(ripPowerAll,1),1);
while true % click "Ctrl+C" to exit when you feel you're done
figure(1);
colors = Bright(1000);
[x,y,button] = ginput(1);
[~,i] = min(abs(x-swDiffAll.*(1-bad))+abs(y-ripPowerAll.*(1-bad)));
t = featureTs/1250;
tl = (1:length(lfp))'/1250;
figure(2);
clf; interval = t(i) + [-1 1]*0.5; in = tl>interval(1) & tl<interval(end);
subplot(3,1,1);
plot(tl(in) - t(i),lfp(in,:)); title(num2str(t(i)));
hold on; plot(bsxfun(@times,ones(1,2),(Restrict(t,interval) - t(i))),ylim,'k--');
legend('ripple channel','SW channel','noise channel');
subplot(3,1,2);
plot(tl(in) - t(i),double(swDiff(in)));
hold on; plot(bsxfun(@times,ones(1,2),(Restrict(t,interval) - t(i))),ylim,'k--');
subplot(3,1,3);
plot(tl(in) - t(i),double(ripPower0(in)));
hold on; plot(bsxfun(@times,ones(1,2),(Restrict(t,interval) - t(i))),ylim,'k--');
%     x = input( prompt )
commandwindow % select the command window
score = str2double(input('good(1)/bad(0)?','s'));
scores(i) = score; % save scores to optionally save
figure(1);
xlims = xlim; ylims = ylim; clf
% extrapolate the score of each ripple based on the score of the nearest scored ripple
scored = find(~isnan(scores));
distances = sqrt(bsxfun(@minus,z(:,1),z(scored,1)').^2 + bsxfun(@minus,z(:,2),z(~isnan(scores),2)').^2);
[~,nearestNeighbour] = min(distances,[],2);
estimated = scores(scored(nearestNeighbour));
selected = estimated>0.1;
idx1 = selected & ~bad; % final ripples
idx2 = ~selected & ~bad; % non-ripples (yet free from noise as well)
scatter(matrix(~bad,1),matrix(~bad,2),1,estimated(~bad),'filled'); colormap(Bright); set(gca,'CLim',[0 1])
xlim(xlims); ylim(ylims);
hold on;
scatter(swDiffAll(scored),ripPowerAll(scored),20,scores(scored)); scatter(swDiffAll(scored),ripPowerAll(scored),15,scores(scored));
set(get(colorbar,'YLabel'),'String','Estimated ripple score');
xlabel('Sharp wave depth'); ylabel('Ripple power');
title({['Manual scoring for ' basepath],['called with channels: ' num2str(Channels) ' (ripple, sharp wave, and (optionally) noise)']});
if rem(i,10)==0
display('saving...');
save(fullfile(basepath,'DetectSWR_manual_scoring.mat'),'t','swDiffAll','ripPowerAll','bad','scores');
end
end
0
1
0
smoothed = Smooth(abs(swDiff),1250*10);
noisyPeriods = ConsolidateIntervals(tl(FindInterval(smoothed>500)),'epsilon',1);
Portion(bad)
bad0 = bad;
Portion(InIntervals(t,noisyPeriods);)
Portion(InIntervals(t,noisyPeriods))
bad = bad | InIntervals(t,noisyPeriods);
idx1 = idx1 & ~bad; idx2 = idx2 & ~bad; % remove bad points from "idx1" and "idx2"
figure(1); % plot a clean figure without these points and check indivudual ripples
clf
plot(swDiffAll(idx1 | idx2),ripPowerAll(idx1 | idx2),'k.','markersize',1); hold on
matrix = [swDiffAll ripPowerAll];
z = bsxfun(@rdivide,bsxfun(@minus,matrix,mean(matrix(~bad,:))),std(matrix(~bad,:)));
scores = nan(size(ripPowerAll,1),1);
while true % click "Ctrl+C" to exit when you feel you're done
figure(1);
colors = Bright(1000);
[x,y,button] = ginput(1);
[~,i] = min(abs(x-swDiffAll.*(1-bad))+abs(y-ripPowerAll.*(1-bad)));
t = featureTs/1250;
tl = (1:length(lfp))'/1250;
figure(2);
clf; interval = t(i) + [-1 1]*0.5; in = tl>interval(1) & tl<interval(end);
subplot(3,1,1);
plot(tl(in) - t(i),lfp(in,:)); title(num2str(t(i)));
hold on; plot(bsxfun(@times,ones(1,2),(Restrict(t,interval) - t(i))),ylim,'k--');
legend('ripple channel','SW channel','noise channel');
subplot(3,1,2);
plot(tl(in) - t(i),double(swDiff(in)));
hold on; plot(bsxfun(@times,ones(1,2),(Restrict(t,interval) - t(i))),ylim,'k--');
subplot(3,1,3);
plot(tl(in) - t(i),double(ripPower0(in)));
hold on; plot(bsxfun(@times,ones(1,2),(Restrict(t,interval) - t(i))),ylim,'k--');
%     x = input( prompt )
commandwindow % select the command window
score = str2double(input('good(1)/bad(0)?','s'));
scores(i) = score; % save scores to optionally save
figure(1);
xlims = xlim; ylims = ylim; clf
% extrapolate the score of each ripple based on the score of the nearest scored ripple
scored = find(~isnan(scores));
distances = sqrt(bsxfun(@minus,z(:,1),z(scored,1)').^2 + bsxfun(@minus,z(:,2),z(~isnan(scores),2)').^2);
[~,nearestNeighbour] = min(distances,[],2);
estimated = scores(scored(nearestNeighbour));
selected = estimated>0.1;
idx1 = selected & ~bad; % final ripples
idx2 = ~selected & ~bad; % non-ripples (yet free from noise as well)
scatter(matrix(~bad,1),matrix(~bad,2),1,estimated(~bad),'filled'); colormap(Bright); set(gca,'CLim',[0 1])
xlim(xlims); ylim(ylims);
hold on;
scatter(swDiffAll(scored),ripPowerAll(scored),20,scores(scored)); scatter(swDiffAll(scored),ripPowerAll(scored),15,scores(scored));
set(get(colorbar,'YLabel'),'String','Estimated ripple score');
xlabel('Sharp wave depth'); ylabel('Ripple power');
title({['Manual scoring for ' basepath],['called with channels: ' num2str(Channels) ' (ripple, sharp wave, and (optionally) noise)']});
if rem(i,10)==0
display('saving...');
save(fullfile(basepath,'DetectSWR_manual_scoring.mat'),'t','swDiffAll','ripPowerAll','bad','scores');
end
end
0
1
0
1
0
1
0
% Save your progress!
display('saving...');
save(fullfile(basepath,'DetectSWR_manual_scoring.mat'),'t','swDiffAll','ripPowerAll','idx1','idx2','bad','selected','scores');
saveas(figure(1),fullfile(basepath,'DetectSWR_manual_scoring.fig'));
dbcont
lfpFile = fullfile(basepath,[basename '.lfp']);
copyfile(lfpFile,fullfile(basepath,[basename '.lfp0']));
file = memmapfile(lfpFile,'Format','int16','Writable',true);
data = reshape(file.Data,64,[]);
m = int16(mean(data(1+[24 11 20 7 9 21 8 25 12 10 23 14 27 26 5 18 28 3 31 62 0 1 30 16 2 29 6 4 36 59 44 34 61 45 33 17 63 32 46 60 35 56 58 40 38 53 42],:)));
newData = bsxfun(@minus,data,m);
file.data = newData(:);
56
ripples = DetectSWR(1+[54 35],'saveMat',true,'check',true,'useSPW',false);
tl = (1:length(lfp))'/1250; % timestamps for lfp
t = featureTs/1250; % timestamps for candidate ripple events
bad = false(size(t));
try % optionally, remove periods of noisy LFP (large deflections)
[clean,~,badIntervals] = CleanLFP([tl,double(lfp(:,2))],'thresholds',[6 Inf],'manual',true);
bad = bad | InIntervals(t,badIntervals);
end
clf; plot(t,d); hold all; plot(xlim,ones(1,2)*threshold2,'r--'); ylabel('lfp derivative (z-units'); legend('derivative','threshold2');
disp('Change "threshold2" manually.')
PlotHVLines([-1 1]*3,'h','r--');
PlotHVLines([-1 1]*1,'h','r--');
x = xlim;
figure; plot(t,z); xlim(x);
threshold1
threshold2
threshold2 = 1;
dbcont
bad = bad | InIntervals(t,badIntervals);
Portion(bad)
[clean,~,badIntervals] = CleanLFP([tl,double(lfp(:,1))],'thresholds',[6 3],'manual',true);
threshold1 =2;
PlotHVLines([-1 1]*1,'h','r--');
PlotHVLines([-1 1]*2,'h','r--');
aroundArtefact1 = 0.5; % 2, Big and long artefacts
threshold1 = 2
clf; plot(t,d); hold all; plot(xlim,ones(1,2)*threshold2,'r--'); ylabel('lfp derivative (z-units'); legend('derivative','threshold2');
threshold2 = 1;
clf; plot(t,d); hold all; plot(xlim,ones(1,2)*threshold2,'r--'); ylabel('lfp derivative (z-units'); legend('derivative','threshold2');
threshold2 = 0.5;
clf; plot(t,d); hold all; plot(xlim,ones(1,2)*threshold2,'r--'); ylabel('lfp derivative (z-units'); legend('derivative','threshold2');
PlotHVLines([-1 1]*0.2,'h','k--');
threshold2 = 0.5;
dbcont
bad = bad | InIntervals(t,badIntervals);
Portion(bad)
smoothed = Smooth(abs(swDiff),1250*10);
noisyPeriods = ConsolidateIntervals(tl(FindInterval(smoothed>500)),'epsilon',1);
bad = bad | InIntervals(t,noisyPeriods);
Portion(bad)
idx1 = idx1 & ~bad; idx2 = idx2 & ~bad; % remove bad points from "idx1" and "idx2"
figure(1); % plot a clean figure without these points and check indivudual ripples
clf
plot(swDiffAll(idx1 | idx2),ripPowerAll(idx1 | idx2),'k.','markersize',1); hold on
matrix = [swDiffAll ripPowerAll];
z = bsxfun(@rdivide,bsxfun(@minus,matrix,mean(matrix(~bad,:))),std(matrix(~bad,:)));
scores = nan(size(ripPowerAll,1),1);
while true % click "Ctrl+C" to exit when you feel you're done
figure(1);
colors = Bright(1000);
[x,y,button] = ginput(1);
[~,i] = min(abs(x-swDiffAll.*(1-bad))+abs(y-ripPowerAll.*(1-bad)));
t = featureTs/1250;
tl = (1:length(lfp))'/1250;
figure(2);
clf; interval = t(i) + [-1 1]*0.5; in = tl>interval(1) & tl<interval(end);
subplot(3,1,1);
plot(tl(in) - t(i),lfp(in,:)); title(num2str(t(i)));
hold on; plot(bsxfun(@times,ones(1,2),(Restrict(t,interval) - t(i))),ylim,'k--');
legend('ripple channel','SW channel','noise channel');
subplot(3,1,2);
plot(tl(in) - t(i),double(swDiff(in)));
hold on; plot(bsxfun(@times,ones(1,2),(Restrict(t,interval) - t(i))),ylim,'k--');
subplot(3,1,3);
plot(tl(in) - t(i),double(ripPower0(in)));
hold on; plot(bsxfun(@times,ones(1,2),(Restrict(t,interval) - t(i))),ylim,'k--');
%     x = input( prompt )
commandwindow % select the command window
score = str2double(input('good(1)/bad(0)?','s'));
scores(i) = score; % save scores to optionally save
figure(1);
xlims = xlim; ylims = ylim; clf
% extrapolate the score of each ripple based on the score of the nearest scored ripple
scored = find(~isnan(scores));
distances = sqrt(bsxfun(@minus,z(:,1),z(scored,1)').^2 + bsxfun(@minus,z(:,2),z(~isnan(scores),2)').^2);
[~,nearestNeighbour] = min(distances,[],2);
estimated = scores(scored(nearestNeighbour));
selected = estimated>0.1;
idx1 = selected & ~bad; % final ripples
idx2 = ~selected & ~bad; % non-ripples (yet free from noise as well)
scatter(matrix(~bad,1),matrix(~bad,2),1,estimated(~bad),'filled'); colormap(Bright); set(gca,'CLim',[0 1])
xlim(xlims); ylim(ylims);
hold on;
scatter(swDiffAll(scored),ripPowerAll(scored),20,scores(scored)); scatter(swDiffAll(scored),ripPowerAll(scored),15,scores(scored));
set(get(colorbar,'YLabel'),'String','Estimated ripple score');
xlabel('Sharp wave depth'); ylabel('Ripple power');
title({['Manual scoring for ' basepath],['called with channels: ' num2str(Channels) ' (ripple, sharp wave, and (optionally) noise)']});
if rem(i,10)==0
display('saving...');
save(fullfile(basepath,'DetectSWR_manual_scoring.mat'),'t','swDiffAll','ripPowerAll','bad','scores');
end
end
0
1
0
1
0
1
% Save your progress!
display('saving...');
save(fullfile(basepath,'DetectSWR_manual_scoring.mat'),'t','swDiffAll','ripPowerAll','idx1','idx2','bad','selected','scores');
saveas(figure(1),fullfile(basepath,'DetectSWR_manual_scoring.fig'));
dbcont
clear all
load('0day13.ripples.events.mat')
load('day13.ripples.events.mat')
close all
SaveFig(['G:\My Drive\Doctorado\Estancia Cornell\Lab\Proyecto Mental Imagery\Labmeeting\figures\ripple_PETHs_' sessionID ])
write = true; % don't save the figures
34
cd M:\home\raly\results\V1\
write
write = true;
load('noiseIntervals.events.mat')
SaveCustomEvents('noise.noi.evt',noiseIntervals(:,1:2),{'noise start','noise stop'});
clear all
45
58.69/16
ripples = DetectSWR(1+[54 35],'saveMat',true,'check',true,'useSPW',false);
% First, remove points that occur in noisy lfp periods
tl = (1:length(lfp))'/1250; % timestamps for lfp
t = featureTs/1250; % timestamps for candidate ripple events
bad = false(size(t));
try % optionally, remove periods of noisy LFP (large deflections)
[clean,~,badIntervals] = CleanLFP([tl,double(lfp(:,2))],'thresholds',[6 Inf],'manual',true);
bad = bad | InIntervals(t,badIntervals);
end
clf; plot(t,d); hold all; plot(xlim,ones(1,2)*threshold2,'r--'); ylabel('lfp derivative (z-units'); legend('derivative','threshold2');
threshold2 = 1;
dbcont
Portion(bad)
[clean,~,badIntervals] = CleanLFP([tl,double(lfp(:,1))],'thresholds',[6 Inf],'manual',true);
threshold1 = 2;
clf; plot(t,d); hold all; plot(xlim,ones(1,2)*threshold2,'r--'); ylabel('lfp derivative (z-units'); legend('derivative','threshold2');
x = xlim;
figure; plot(t,z); xlim(x);
threshold2 = 1;
dbcont
bad = bad | InIntervals(t,badIntervals);
Portion(bad)
smoothed = Smooth(abs(swDiff),1250*10);
noisyPeriods = ConsolidateIntervals(tl(FindInterval(smoothed>500)),'epsilon',1);
bad = bad | InIntervals(t,noisyPeriods);
Portion(bad)
idx1 = idx1 & ~bad; idx2 = idx2 & ~bad; % remove bad points from "idx1" and "idx2"
figure(1); % plot a clean figure without these points and check indivudual ripples
clf
plot(swDiffAll(idx1 | idx2),ripPowerAll(idx1 | idx2),'k.','markersize',1); hold on
matrix = [swDiffAll ripPowerAll];
z = bsxfun(@rdivide,bsxfun(@minus,matrix,mean(matrix(~bad,:))),std(matrix(~bad,:)));
scores = nan(size(ripPowerAll,1),1);
while true % click "Ctrl+C" to exit when you feel you're done
figure(1);
colors = Bright(1000);
[x,y,button] = ginput(1);
[~,i] = min(abs(x-swDiffAll.*(1-bad))+abs(y-ripPowerAll.*(1-bad)));
t = featureTs/1250;
tl = (1:length(lfp))'/1250;
figure(2);
clf; interval = t(i) + [-1 1]*0.5; in = tl>interval(1) & tl<interval(end);
subplot(3,1,1);
plot(tl(in) - t(i),lfp(in,:)); title(num2str(t(i)));
hold on; plot(bsxfun(@times,ones(1,2),(Restrict(t,interval) - t(i))),ylim,'k--');
legend('ripple channel','SW channel','noise channel');
subplot(3,1,2);
plot(tl(in) - t(i),double(swDiff(in)));
hold on; plot(bsxfun(@times,ones(1,2),(Restrict(t,interval) - t(i))),ylim,'k--');
subplot(3,1,3);
plot(tl(in) - t(i),double(ripPower0(in)));
hold on; plot(bsxfun(@times,ones(1,2),(Restrict(t,interval) - t(i))),ylim,'k--');
%     x = input( prompt )
commandwindow % select the command window
score = str2double(input('good(1)/bad(0)?','s'));
scores(i) = score; % save scores to optionally save
figure(1);
xlims = xlim; ylims = ylim; clf
% extrapolate the score of each ripple based on the score of the nearest scored ripple
scored = find(~isnan(scores));
distances = sqrt(bsxfun(@minus,z(:,1),z(scored,1)').^2 + bsxfun(@minus,z(:,2),z(~isnan(scores),2)').^2);
[~,nearestNeighbour] = min(distances,[],2);
estimated = scores(scored(nearestNeighbour));
selected = estimated>0.1;
idx1 = selected & ~bad; % final ripples
idx2 = ~selected & ~bad; % non-ripples (yet free from noise as well)
scatter(matrix(~bad,1),matrix(~bad,2),1,estimated(~bad),'filled'); colormap(Bright); set(gca,'CLim',[0 1])
xlim(xlims); ylim(ylims);
hold on;
scatter(swDiffAll(scored),ripPowerAll(scored),20,scores(scored)); scatter(swDiffAll(scored),ripPowerAll(scored),15,scores(scored));
set(get(colorbar,'YLabel'),'String','Estimated ripple score');
xlabel('Sharp wave depth'); ylabel('Ripple power');
title({['Manual scoring for ' basepath],['called with channels: ' num2str(Channels) ' (ripple, sharp wave, and (optionally) noise)']});
if rem(i,10)==0
display('saving...');
save(fullfile(basepath,'DetectSWR_manual_scoring.mat'),'t','swDiffAll','ripPowerAll','bad','scores');
end
end
0
1
0
1
0
1
0
1
0
1
0
% Save your progress!
display('saving...');
save(fullfile(basepath,'DetectSWR_manual_scoring.mat'),'t','swDiffAll','ripPowerAll','idx1','idx2','bad','selected','scores');
saveas(figure(1),fullfile(basepath,'DetectSWR_manual_scoring.fig'));
clf
plot(swDiffAll(idx2),ripPowerAll(idx2),'k.','markersize',1); hold on
plot(swDiffAll(idx1),ripPowerAll(idx1),'r.','markersize',1); legend('non-ripples','ripples');
disp(['Please review the figure and change the idx1 (ripples) ' newline 'and idx2 (calm non-ripple periods) groups.' newline ...
'Note that noisy periods need not participate in either group.']);
idx2(swDiffAll>1532) = 0;
dbcont
fclose('all');
clear file
clear all
write = true;
clear all
% don't forget to rename lfp file to lfpN for normalized
1
raly
% First, remove points that occur in noisy lfp periods
tl = (1:length(lfp))'/1250; % timestamps for lfp
t = featureTs/1250; % timestamps for candidate ripple events
bad = false(size(t));
try % optionally, remove periods of noisy LFP (large deflections)
[clean,~,badIntervals] = CleanLFP([tl,double(lfp(:,2))],'thresholds',[6 Inf],'manual',true);
bad = bad | InIntervals(t,badIntervals);
end
clf
clf; plot(t,d); hold all; plot(xlim,ones(1,2)*threshold2,'r--'); ylabel('lfp derivative (z-units'); legend('derivative','threshold2');
PlotHVLines([-1 1]*1,'h','k--');
threshold2 = 2;
dbcont
[clean,~,badIntervals] = CleanLFP([tl,double(lfp(:,1))],'thresholds',[6 Inf],'manual',true);
bad = bad | InIntervals(t,badIntervals);
threshold1 = 2;
clf; plot(t,d); hold all; plot(xlim,ones(1,2)*threshold2,'r--'); ylabel('lfp derivative (z-units'); legend('derivative','threshold2');
threshold2 = 1;
clf; plot(t,d); hold all; plot(xlim,ones(1,2)*threshold2,'r--'); ylabel('lfp derivative (z-units'); legend('derivative','threshold2');
plot(xlim,-ones(1,2)*threshold2,'r--'); ylabel('lfp derivative (z-units'); legend('derivative','threshold2');
threshold2 = 0.3;
plot(xlim,-ones(1,2)*threshold2,'r--'); ylabel('lfp derivative (z-units'); legend('derivative','threshold2');
plot(xlim,ones(1,2)*threshold2,'r--'); ylabel('lfp derivative (z-units'); legend('derivative','threshold2');
threshold2
plot(xlim,-ones(1,2)*threshold2,'r--'); ylabel('lfp derivative (z-units'); legend('derivative','threshold2');
plot(xlim,ones(1,2)*threshold2,'r--'); ylabel('lfp derivative (z-units'); legend('derivative','threshold2');
threshold1
threshold2
dbcont
Portion(bad)
idx1 = idx1 & ~bad; idx2 = idx2 & ~bad; % remove bad points from "idx1" and "idx2"
figure(1); % plot a clean figure without these points and check indivudual ripples
clf
plot(swDiffAll(idx1 | idx2),ripPowerAll(idx1 | idx2),'k.','markersize',1); hold on
matrix = [swDiffAll ripPowerAll];
z = bsxfun(@rdivide,bsxfun(@minus,matrix,mean(matrix(~bad,:))),std(matrix(~bad,:)));
scores = nan(size(ripPowerAll,1),1);
while true % click "Ctrl+C" to exit when you feel you're done
figure(1);
colors = Bright(1000);
[x,y,button] = ginput(1);
[~,i] = min(abs(x-swDiffAll.*(1-bad))+abs(y-ripPowerAll.*(1-bad)));
t = featureTs/1250;
tl = (1:length(lfp))'/1250;
figure(2);
clf; interval = t(i) + [-1 1]*0.5; in = tl>interval(1) & tl<interval(end);
subplot(3,1,1);
plot(tl(in) - t(i),lfp(in,:)); title(num2str(t(i)));
hold on; plot(bsxfun(@times,ones(1,2),(Restrict(t,interval) - t(i))),ylim,'k--');
legend('ripple channel','SW channel','noise channel');
subplot(3,1,2);
plot(tl(in) - t(i),double(swDiff(in)));
hold on; plot(bsxfun(@times,ones(1,2),(Restrict(t,interval) - t(i))),ylim,'k--');
subplot(3,1,3);
plot(tl(in) - t(i),double(ripPower0(in)));
hold on; plot(bsxfun(@times,ones(1,2),(Restrict(t,interval) - t(i))),ylim,'k--');
%     x = input( prompt )
commandwindow % select the command window
score = str2double(input('good(1)/bad(0)?','s'));
scores(i) = score; % save scores to optionally save
figure(1);
xlims = xlim; ylims = ylim; clf
% extrapolate the score of each ripple based on the score of the nearest scored ripple
scored = find(~isnan(scores));
distances = sqrt(bsxfun(@minus,z(:,1),z(scored,1)').^2 + bsxfun(@minus,z(:,2),z(~isnan(scores),2)').^2);
[~,nearestNeighbour] = min(distances,[],2);
estimated = scores(scored(nearestNeighbour));
selected = estimated>0.1;
idx1 = selected & ~bad; % final ripples
idx2 = ~selected & ~bad; % non-ripples (yet free from noise as well)
scatter(matrix(~bad,1),matrix(~bad,2),1,estimated(~bad),'filled'); colormap(Bright); set(gca,'CLim',[0 1])
xlim(xlims); ylim(ylims);
hold on;
scatter(swDiffAll(scored),ripPowerAll(scored),20,scores(scored)); scatter(swDiffAll(scored),ripPowerAll(scored),15,scores(scored));
set(get(colorbar,'YLabel'),'String','Estimated ripple score');
xlabel('Sharp wave depth'); ylabel('Ripple power');
title({['Manual scoring for ' basepath],['called with channels: ' num2str(Channels) ' (ripple, sharp wave, and (optionally) noise)']});
if rem(i,10)==0
display('saving...');
save(fullfile(basepath,'DetectSWR_manual_scoring.mat'),'t','swDiffAll','ripPowerAll','bad','scores');
end
end
0
1
0
1
0
1
0
1
0
1
0
display('saving...');
save(fullfile(basepath,'DetectSWR_manual_scoring.mat'),'t','swDiffAll','ripPowerAll','idx1','idx2','bad','selected','scores');
saveas(figure(1),fullfile(basepath,'DetectSWR_manual_scoring.fig'));
idx2(swDiffAll>2000) = 0;
clf
plot(swDiffAll(idx2),ripPowerAll(idx2),'k.','markersize',1); hold on
plot(swDiffAll(idx1),ripPowerAll(idx1),'r.','markersize',1); legend('non-ripples','ripples');
dbcont
clear file
[spikes,regionID,regionNames] = GetAyaSpikes(basepath);
% Get a single string to include in figure labels
regionListString = regionNames; regionListString(2,:) = repmat({'->'},1,size(regionListString,2)); regionListString{end} = '';
regionListString = strcat(regionListString{:});
write = false; % don't save the figures
45
write = true; % don't save the figures
45
write
write = true; % don't save the figures
sleep
clf
bins = Bins(0,sleep(end),60*10,60*10); rippleRate = CountInIntervals(r(:,1),bins)/diff(bins(1,:));
bar(mean(bins,2),rippleRate);
handle = PlotIntervals(SubtractIntervals([0 Inf],sleep),'color','r');
legend(handle,'maze','Location','northwestsleep','box','off','fontsize',15);
set(gca,'box','off','tickdir','out','fontsize',12);
xlabel('time (s)')
ylabel('ripple rate (per second)')
max(ripples(:,2))
max(r(:,2))
sleep(end)
PlotIntervals(SubtractIntervals([0 Inf],sleep),'color','r');
SubtractIntervals([0 Inf],sleep)
clear all
nNeurons = max(spikes(:,2));
r = ripples.timestamps;
pre = InIntervals(r(:,1),sleep(1,:));
post = InIntervals(r(:,1),sleep(2,:));
clear all
try load(fullfile(basepath,[basename '.MergePoints.events.mat']),'MergePoints');
postID = find(cellfun(@(x) ~isempty(strfind(x,'post')), MergePoints.foldernames));
sleep = ConsolidateIntervals(MergePoints.timestamps(cellfun(@(x) any(strfind(lower(x),'sleep')),MergePoints.foldernames),:));
catch
MergePoints = [];
display('No merge points saved. Keeping all ripples (impossible to restrict to post-task sleep');
end
nNeurons = max(spikes(:,2));
r = ripples.timestamps;
pre = InIntervals(r(:,1),sleep(1,:));
post = InIntervals(r(:,1),sleep(2,:));
write
SaveFig(['M:\home\raly\results\V1\rippleRate_' sessionID])
clear all
456
% First, remove points that occur in noisy lfp periods
tl = (1:length(lfp))'/1250; % timestamps for lfp
t = featureTs/1250; % timestamps for candidate ripple events
bad = false(size(t));
try % optionally, remove periods of noisy LFP (large deflections)
[clean,~,badIntervals] = CleanLFP([tl,double(lfp(:,1))],'thresholds',[6 Inf],'manual',true);
bad = bad | InIntervals(t,badIntervals);
end
ripples = DetectSWR(1+[54 35],'saveMat',true,'check',true,'useSPW',false,'useEEG',true);
close all
% First, remove points that occur in noisy lfp periods
tl = (1:length(lfp))'/1250; % timestamps for lfp
t = featureTs/1250; % timestamps for candidate ripple events
bad = false(size(t));
try % optionally, remove periods of noisy LFP (large deflections)
[clean,~,badIntervals] = CleanLFP([tl,double(lfp(:,1))],'thresholds',[6 Inf],'manual',true);
bad = bad | InIntervals(t,badIntervals);
end
Channels
dbcont
bad = false(size(t));
Channels
bad = false(size(t));
try % optionally, remove periods of noisy LFP (large deflections)
[clean,~,badIntervals] = CleanLFP([tl,double(lfp(:,1))],'thresholds',[6 Inf],'manual',true);
bad = bad | InIntervals(t,badIntervals);
end
threshold1 = 3;
PlotHVLines([-1 1]*3,'h','k--');
PlotHVLines([-1 1]*2,'h','r--');
threshold1 = 2;
threshold2 = 1;
clf; plot(t,d); hold all; plot(xlim,ones(1,2)*threshold2,'r--'); ylabel('lfp derivative (z-units'); legend('derivative','threshold2');
plot(xlim,-ones(1,2)*threshold2,'r--'); ylabel('lfp derivative (z-units'); legend('derivative','threshold2');
threshold2 = 0.3;
plot(xlim,-ones(1,2)*threshold2,'r--'); ylabel('lfp derivative (z-units'); legend('derivative','threshold2');
plot(xlim,ones(1,2)*threshold2,'r--'); ylabel('lfp derivative (z-units'); legend('derivative','threshold2');
threshold2 = 0.4;
PlotHVLines([-1 1]*0.4,'h','y--');
PlotHVLines([-1 1]*0.3,'h','w--');
PlotHVLines([-1 1]*0.4,'h','r--');
PlotHVLines([-1 1]*1,'h','w--');
dbcont
Portion(bad)
[clean,~,badIntervals] = CleanLFP([tl,double(lfp(:,2))],'thresholds',[6 Inf],'manual',true);
bad = bad | InIntervals(t,badIntervals);
threshold2 = 1;
clf; plot(t,d); hold all; plot(xlim,ones(1,2)*threshold2,'r--'); ylabel('lfp derivative (z-units'); legend('derivative','threshold2');
threshold1
threshold2 = 2;
dbcont
Portion(bad)
idx1 = idx1 & ~bad; idx2 = idx2 & ~bad; % remove bad points from "idx1" and "idx2"
figure(1); % plot a clean figure without these points and check indivudual ripples
clf
plot(swDiffAll(idx1 | idx2),ripPowerAll(idx1 | idx2),'k.','markersize',1); hold on
matrix = [swDiffAll ripPowerAll];
z = bsxfun(@rdivide,bsxfun(@minus,matrix,mean(matrix(~bad,:))),std(matrix(~bad,:)));
scores = nan(size(ripPowerAll,1),1);
while true % click "Ctrl+C" to exit when you feel you're done
figure(1);
colors = Bright(1000);
[x,y,button] = ginput(1);
[~,i] = min(abs(x-swDiffAll.*(1-bad))+abs(y-ripPowerAll.*(1-bad)));
t = featureTs/1250;
tl = (1:length(lfp))'/1250;
figure(2);
clf; interval = t(i) + [-1 1]*0.5; in = tl>interval(1) & tl<interval(end);
subplot(3,1,1);
plot(tl(in) - t(i),lfp(in,:)); title(num2str(t(i)));
hold on; plot(bsxfun(@times,ones(1,2),(Restrict(t,interval) - t(i))),ylim,'k--');
legend('ripple channel','SW channel','noise channel');
subplot(3,1,2);
plot(tl(in) - t(i),double(swDiff(in)));
hold on; plot(bsxfun(@times,ones(1,2),(Restrict(t,interval) - t(i))),ylim,'k--');
subplot(3,1,3);
plot(tl(in) - t(i),double(ripPower0(in)));
hold on; plot(bsxfun(@times,ones(1,2),(Restrict(t,interval) - t(i))),ylim,'k--');
%     x = input( prompt )
commandwindow % select the command window
score = str2double(input('good(1)/bad(0)?','s'));
scores(i) = score; % save scores to optionally save
figure(1);
xlims = xlim; ylims = ylim; clf
% extrapolate the score of each ripple based on the score of the nearest scored ripple
scored = find(~isnan(scores));
distances = sqrt(bsxfun(@minus,z(:,1),z(scored,1)').^2 + bsxfun(@minus,z(:,2),z(~isnan(scores),2)').^2);
[~,nearestNeighbour] = min(distances,[],2);
estimated = scores(scored(nearestNeighbour));
selected = estimated>0.1;
idx1 = selected & ~bad; % final ripples
idx2 = ~selected & ~bad; % non-ripples (yet free from noise as well)
scatter(matrix(~bad,1),matrix(~bad,2),1,estimated(~bad),'filled'); colormap(Bright); set(gca,'CLim',[0 1])
xlim(xlims); ylim(ylims);
hold on;
scatter(swDiffAll(scored),ripPowerAll(scored),20,scores(scored)); scatter(swDiffAll(scored),ripPowerAll(scored),15,scores(scored));
set(get(colorbar,'YLabel'),'String','Estimated ripple score');
xlabel('Sharp wave depth'); ylabel('Ripple power');
title({['Manual scoring for ' basepath],['called with channels: ' num2str(Channels) ' (ripple, sharp wave, and (optionally) noise)']});
if rem(i,10)==0
display('saving...');
save(fullfile(basepath,'DetectSWR_manual_scoring.mat'),'t','swDiffAll','ripPowerAll','bad','scores');
end
end
0
1
0
1
0
1
0
1
0
1
idx2(swDiffAll>2000) = 0;
% Save your progress!
display('saving...');
save(fullfile(basepath,'DetectSWR_manual_scoring.mat'),'t','swDiffAll','ripPowerAll','idx1','idx2','bad','selected','scores');
saveas(figure(1),fullfile(basepath,'DetectSWR_manual_scoring.fig'));
clf
plot(swDiffAll(idx2),ripPowerAll(idx2),'k.','markersize',1); hold on
plot(swDiffAll(idx1),ripPowerAll(idx1),'r.','markersize',1); legend('non-ripples','ripples');
dbcont
basepath
write=true
45
clear all
46
% First, remove points that occur in noisy lfp periods
tl = (1:length(lfp))'/1250; % timestamps for lfp
t = featureTs/1250; % timestamps for candidate ripple events
bad = false(size(t));
try % optionally, remove periods of noisy LFP (large deflections)
[clean,~,badIntervals] = CleanLFP([tl,double(lfp(:,2))],'thresholds',[6 Inf],'manual',true);
bad = bad | InIntervals(t,badIntervals);
end
threshold2 = 2;
clf; plot(t,d); hold all; plot(xlim,ones(1,2)*threshold2,'r--'); ylabel('lfp derivative (z-units'); legend('derivative','threshold2');
dbcont
try % optionally, remove periods of noisy LFP (large deflections)
[clean,~,badIntervals] = CleanLFP([tl,double(lfp(:,1))],'thresholds',[6 Inf],'manual',true);
bad = bad | InIntervals(t,badIntervals);
end
threshold2 = 2;
clf; plot(t,d); hold all; plot(xlim,ones(1,2)*threshold2,'r--'); ylabel('lfp derivative (z-units'); legend('derivative','threshold2');
threshold2 = 0.4;
dbcont
Portion(bad)
idx1 = idx1 & ~bad; idx2 = idx2 & ~bad; % remove bad points from "idx1" and "idx2"
figure(1); % plot a clean figure without these points and check indivudual ripples
clf
plot(swDiffAll(idx1 | idx2),ripPowerAll(idx1 | idx2),'k.','markersize',1); hold on
matrix = [swDiffAll ripPowerAll];
z = bsxfun(@rdivide,bsxfun(@minus,matrix,mean(matrix(~bad,:))),std(matrix(~bad,:)));
scores = nan(size(ripPowerAll,1),1);
while true % click "Ctrl+C" to exit when you feel you're done
figure(1);
colors = Bright(1000);
[x,y,button] = ginput(1);
[~,i] = min(abs(x-swDiffAll.*(1-bad))+abs(y-ripPowerAll.*(1-bad)));
t = featureTs/1250;
tl = (1:length(lfp))'/1250;
figure(2);
clf; interval = t(i) + [-1 1]*0.5; in = tl>interval(1) & tl<interval(end);
subplot(3,1,1);
plot(tl(in) - t(i),lfp(in,:)); title(num2str(t(i)));
hold on; plot(bsxfun(@times,ones(1,2),(Restrict(t,interval) - t(i))),ylim,'k--');
legend('ripple channel','SW channel','noise channel');
subplot(3,1,2);
plot(tl(in) - t(i),double(swDiff(in)));
hold on; plot(bsxfun(@times,ones(1,2),(Restrict(t,interval) - t(i))),ylim,'k--');
subplot(3,1,3);
plot(tl(in) - t(i),double(ripPower0(in)));
hold on; plot(bsxfun(@times,ones(1,2),(Restrict(t,interval) - t(i))),ylim,'k--');
%     x = input( prompt )
commandwindow % select the command window
score = str2double(input('good(1)/bad(0)?','s'));
scores(i) = score; % save scores to optionally save
figure(1);
xlims = xlim; ylims = ylim; clf
% extrapolate the score of each ripple based on the score of the nearest scored ripple
scored = find(~isnan(scores));
distances = sqrt(bsxfun(@minus,z(:,1),z(scored,1)').^2 + bsxfun(@minus,z(:,2),z(~isnan(scores),2)').^2);
[~,nearestNeighbour] = min(distances,[],2);
estimated = scores(scored(nearestNeighbour));
selected = estimated>0.1;
idx1 = selected & ~bad; % final ripples
idx2 = ~selected & ~bad; % non-ripples (yet free from noise as well)
scatter(matrix(~bad,1),matrix(~bad,2),1,estimated(~bad),'filled'); colormap(Bright); set(gca,'CLim',[0 1])
xlim(xlims); ylim(ylims);
hold on;
scatter(swDiffAll(scored),ripPowerAll(scored),20,scores(scored)); scatter(swDiffAll(scored),ripPowerAll(scored),15,scores(scored));
set(get(colorbar,'YLabel'),'String','Estimated ripple score');
xlabel('Sharp wave depth'); ylabel('Ripple power');
title({['Manual scoring for ' basepath],['called with channels: ' num2str(Channels) ' (ripple, sharp wave, and (optionally) noise)']});
if rem(i,10)==0
display('saving...');
save(fullfile(basepath,'DetectSWR_manual_scoring.mat'),'t','swDiffAll','ripPowerAll','bad','scores');
end
end
0
1
0
1
0
1
0
1
0
1
0
scores(i) = 1;
idx1 = idx1 & ~bad; idx2 = idx2 & ~bad; % remove bad points from "idx1" and "idx2"
figure(1); % plot a clean figure without these points and check indivudual ripples
clf
plot(swDiffAll(idx1 | idx2),ripPowerAll(idx1 | idx2),'k.','markersize',1); hold on
matrix = [swDiffAll ripPowerAll];
z = bsxfun(@rdivide,bsxfun(@minus,matrix,mean(matrix(~bad,:))),std(matrix(~bad,:)));
scores = nan(size(ripPowerAll,1),1);
while true % click "Ctrl+C" to exit when you feel you're done
figure(1);
colors = Bright(1000);
[x,y,button] = ginput(1);
[~,i] = min(abs(x-swDiffAll.*(1-bad))+abs(y-ripPowerAll.*(1-bad)));
t = featureTs/1250;
tl = (1:length(lfp))'/1250;
figure(2);
clf; interval = t(i) + [-1 1]*0.5; in = tl>interval(1) & tl<interval(end);
subplot(3,1,1);
plot(tl(in) - t(i),lfp(in,:)); title(num2str(t(i)));
hold on; plot(bsxfun(@times,ones(1,2),(Restrict(t,interval) - t(i))),ylim,'k--');
legend('ripple channel','SW channel','noise channel');
subplot(3,1,2);
plot(tl(in) - t(i),double(swDiff(in)));
hold on; plot(bsxfun(@times,ones(1,2),(Restrict(t,interval) - t(i))),ylim,'k--');
subplot(3,1,3);
plot(tl(in) - t(i),double(ripPower0(in)));
hold on; plot(bsxfun(@times,ones(1,2),(Restrict(t,interval) - t(i))),ylim,'k--');
%     x = input( prompt )
commandwindow % select the command window
score = str2double(input('good(1)/bad(0)?','s'));
scores(i) = score; % save scores to optionally save
figure(1);
xlims = xlim; ylims = ylim; clf
% extrapolate the score of each ripple based on the score of the nearest scored ripple
scored = find(~isnan(scores));
distances = sqrt(bsxfun(@minus,z(:,1),z(scored,1)').^2 + bsxfun(@minus,z(:,2),z(~isnan(scores),2)').^2);
[~,nearestNeighbour] = min(distances,[],2);
estimated = scores(scored(nearestNeighbour));
selected = estimated>0.1;
idx1 = selected & ~bad; % final ripples
idx2 = ~selected & ~bad; % non-ripples (yet free from noise as well)
scatter(matrix(~bad,1),matrix(~bad,2),1,estimated(~bad),'filled'); colormap(Bright); set(gca,'CLim',[0 1])
xlim(xlims); ylim(ylims);
hold on;
scatter(swDiffAll(scored),ripPowerAll(scored),20,scores(scored)); scatter(swDiffAll(scored),ripPowerAll(scored),15,scores(scored));
set(get(colorbar,'YLabel'),'String','Estimated ripple score');
xlabel('Sharp wave depth'); ylabel('Ripple power');
title({['Manual scoring for ' basepath],['called with channels: ' num2str(Channels) ' (ripple, sharp wave, and (optionally) noise)']});
if rem(i,10)==0
display('saving...');
save(fullfile(basepath,'DetectSWR_manual_scoring.mat'),'t','swDiffAll','ripPowerAll','bad','scores');
end
end
close all
lfp= GetAyaLFP(54);
[clean,~,badIntervals] = CleanLFP(lfp,'thresholds',[2 2],'manual',true);
clf; plot(t,d); hold all; plot(xlim,ones(1,2)*threshold2,'r--'); ylabel('lfp derivative (z-units'); legend('derivative','threshold2');
threshold2 = 1;
dbcont
goodIntervals = SubtractIntervals([0 lfp(end,1)],badIntervals);
dur(badIntervals)
clf
PlotXY(lfp);
PlotIntervals(badIntervals);
basepath
raly
lfp0 = lfp;
lfp = GetAyaEEG(54);
hold all
PlotXY(lfp);
close
[clean,~,badIntervals] = CleanLFP(lfp,'thresholds',[2 2],'manual',true);
aroundArtefact1 = 1; % 2, Big and long artefacts
clf; plot(t,d); hold all; plot(xlim,ones(1,2)*threshold2,'r--'); ylabel('lfp derivative (z-units'); legend('derivative','threshold2');
threshold2 = 1;
dbcont
goodIntervals = SubtractIntervals([0 lfp(end,1)],badIntervals);
dur(badIntervals)
dur(badIntervals)/3600
goodIntervals = SubtractIntervals([0 lfp(end,1)],badIntervals);
ripples = FindRipples(basepath,1+54,'noise',1+7,'saveMat',true,'restrict',goodIntervals);
ripples = FindRipples(lfp(:,2),lfp(:,1),'noise',1+7,'saveMat',true,'restrict',goodIntervals);
[datasetfolder,recordingname] = fileparts(basepath);
matfilename = fullfile(basepath,[recordingname,'.EMGFromLFP.LFP.mat']);
%% xmlPameters
p = inputParser;
addParameter(p,'restrict',[0 inf],@isnumeric)
addParameter(p,'specialChannels',[],@isnumeric)
addParameter(p,'rejectChannels',[],@isnumeric)
addParameter(p,'restrictChannels',[],@isnumeric)
addParameter(p,'saveMat',true,@islogical)
addParameter(p,'saveLocation','',@isstr)
addParameter(p,'overwrite',false,@islogical)
addParameter(p,'samplingFrequency',2,@isnumeric)
addParameter(p,'noPrompts',false,@islogical);
addParameter(p,'fromDat',false,@islogical);
addParameter(p,'chInfo',[],@isstruct);
parse(p,varargin{:})
restrict = p.Results.restrict;
specialChannels = p.Results.specialChannels;
rejectChannels = p.Results.rejectChannels;
restrictChannels = p.Results.restrictChannels;
saveMat = p.Results.saveMat;
overwrite = p.Results.overwrite;
samplingFrequency = p.Results.samplingFrequency;
noPrompts = p.Results.noPrompts;
fromDat = p.Results.fromDat;
chInfo = p.Results.chInfo;
if ~isempty(p.Results.saveLocation)
matfilename = fullfile(p.Results.saveLocation,[recordingname,'.EMGFromLFP.LFP.mat']);
end
%% Check if EMGCorr has already been calculated for this recording
%If the EMGCorr file already exists, load and return with EMGCorr in hand
if exist(matfilename,'file') && ~overwrite
display('EMGFromLFP Correlation already calculated - loading from EMGFromLFP.LFP.mat')
load(matfilename)
if exist('EMGCorr','var')%for backcompatability
EMGFromLFP = EMGCorr;
end
if ~exist('EMGFromLFP','var')
display([matfilename,' does not contain a variable called EMGFromLFP'])
end
return
end
display('Calculating EMGFromLFP from High Frequency LFP Correlation')
load(fullfile(basepath,[recordingname,'.session.mat']))
nChannels = session.extracellular.nChannels;
SpkGrps = session.extracellular.spikeGroups.channels;
Fs = session.extracellular.srLfp;
lfpFile = checkFile('basepath',basepath,'fileTypes',{'.lfp','.eeg'});
lfpFile = [basepath filesep lfpFile(1).name];
%% get basics about.lfp/lfp file
% if ~isempty(chInfo)
%     nChannels = chInfo.nChannel;
%     SpkGrps = chInfo.one.AnatGrps;
%     Fs = chInfo.lfpSR;
%     lfpFile = checkFile('basepath',basepath,'fileType','.lfp');
%     lfpFile = [basepath filesep lfpFile(1).name];
% end
% sessionInfo = bz_getSessionInfo(basePath,'noPrompts',noPrompts); % now using the updated version
% switch fromDat
%     case false
%         if exist([basePath filesep sessionInfo.FileName '.lfp'])
%             lfpFile = [basePath filesep sessionInfo.FileName '.lfp'];
%         elseif exist([basePath filesep sessionInfo.FileName '.eeg'])
%             lfpFile = [basePath filesep sessionInfo.FileName '.eeg'];
%         else
%             error('could not find an LFP or EEG file...')
%         end
%
%         Fs = sessionInfo.lfpSampleRate; % Hz, LFP sampling rate
%
%
%     case true
%         if exist([basePath filesep sessionInfo.FileName '.dat'])
%             datFile = [basePath filesep sessionInfo.FileName '.dat'];
%         else
%             error('could not find a dat file...')
%         end
%
%         datFs = sessionInfo.rates.wideband;
%         Fs = sessionInfo.lfpSampleRate; % Hz, LFP sampling rate
% end
% nChannels = sessionInfo.nChannels;
%
% if isfield(sessionInfo,'SpkGrps')
%     SpkGrps = sessionInfo.SpkGrps;
% elseif isfield(sessionInfo,'AnatGrps')
%     SpkGrps = sessionInfo.AnatGrps;
%     display('No SpikeGroups, Using AnatomyGroups')
% else
%     error('No SpikeGroups...')
% end
binScootS = 1 ./ samplingFrequency;
binScootSamps = round(Fs*binScootS); % must be integer, or error on line 190
corrChunkSz = 20; %for batch-processed correlations
%% Pick channels and to analyze
% get spike groups,
% pick every other one... unless specialshanks, in which case pick non-adjacent
%This is potentially dangerous in combination with rejectChannels... i.e.
%what if you pick every other shank but then the ones you pick are all
%reject because noisy shank.
% xcorrs_chs is a list of channels that will be loaded
% spkgrpstouse is a list of spike groups to find channels from
if ~isempty(restrictChannels)    % If restrict channel case:
xcorr_chs = restrictChannels;
else
% get list of spike groups (aka shanks) that should be used
usablechannels = [];
spkgrpstouse = [];
for gidx = 1:length(SpkGrps)
usableshankchannels{gidx} = setdiff(SpkGrps{gidx},rejectChannels);
usablechannels = cat(2,usablechannels,usableshankchannels{gidx});
if ~isempty(usableshankchannels{gidx})
spkgrpstouse = cat(2,spkgrpstouse,gidx);
end
end
% check for good/bad shanks and update here
% spkgrpstouse = unique(cat(1,spkgrpstouse,specialshanks)); % this is redundant with taking all shanks.
% get list of channels (1 from each good spike group)
xcorr_chs = [];
for gidx=1:length(usableshankchannels)
%Remove rejectChannels
%     usableshankchannels = setdiff(SpkGrps(spkgrpstouse(i)).Channels,rejectChannels);
%grab random channel on each shank
if ~isempty(usableshankchannels{gidx})
randChfromShank = usableshankchannels{gidx}(randi(length(usableshankchannels{gidx})));
xcorr_chs = [xcorr_chs,randChfromShank];
end
end
xcorr_chs = unique([xcorr_chs,specialChannels]);
end
xcorr_chs
xcorr_chs = usableshankchannels{1}(1:5:end)
54
xcorr_chs = usableshankchannels{1}(1:10:end)
xcorr_chs = 1+[11 27 0 44 35 53];
45
%% xcorr 'strength' is the summed correlation coefficients between channel
% pairs for a sliding window of 25 ms
xcorr_window_samps = round(binScootS*Fs);
xcorr_window_inds = -xcorr_window_samps:xcorr_window_samps;%+- that number of ms in samples
% new version... batches of correlation calculated at once
timestamps = (1+xcorr_window_inds(end)):binScootSamps:(size(lfp,1)-xcorr_window_inds(end));
numbins = length(timestamps);
EMGCorr = zeros(numbins, 1);
% tic
counter = 1;
for j=1:(length(xcorr_chs))
for k=(j+1):length(xcorr_chs)
%disp([num2str(counter*2 ./ (length(xcorr_chs)*length(xcorr_chs)*length(timestamps)))])
bz_Counter(counter,(length(xcorr_chs)*(length(xcorr_chs)-1))./2,'Channel Pair')
c1 = [];
c2 = [];
binind = 0;
binindstart = 1;
for i = timestamps
binind = binind+1;
s1 =lfp(i + xcorr_window_inds, j);
s2 =lfp(i + xcorr_window_inds, k);
c1 = cat(2,c1,s1);
c2 = cat(2,c2,s2);
if size(c1,2) == corrChunkSz || i == timestamps(end)
binindend = binind;
tmp = corr(c1,c2);
tmp = diag(tmp);
EMGCorr(binindstart:binindend) = EMGCorr(binindstart:binindend) + tmp;
c1 = [];
c2 = [];
binindstart = binind+1;
end
end
counter = counter+1;
end
end
% toc
EMGCorr = EMGCorr/(length(xcorr_chs)*(length(xcorr_chs)-1)/2); % normalize
EMGFromLFP.timestamps = timestamps'./Fs;
EMGFromLFP.data = EMGCorr;
EMGFromLFP.channels = xcorr_chs;
EMGFromLFP.detectorName = 'getEMGFromLFP';
EMGFromLFP.samplingFrequency = samplingFrequency;
if saveMat
%Save in buzcodeformat
save(matfilename,'EMGFromLFP');
end
dbcont
[spikes,regionID,regionNames] = GetAyaSpikes(basepath);
% Get a single string to include in figure labels
regionListString = regionNames; regionListString(2,:) = repmat({'->'},1,size(regionListString,2)); regionListString{end} = '';
regionListString = strcat(regionListString{:});
try load(fullfile(basepath,[basename '.MergePoints.events.mat']),'MergePoints');
postID = find(cellfun(@(x) ~isempty(strfind(x,'post')), MergePoints.foldernames));
sleep = ConsolidateIntervals(MergePoints.timestamps(cellfun(@(x) any(strfind(lower(x),'sleep')),MergePoints.foldernames),:));
catch
MergePoints = [];
display('No merge points saved. Keeping all ripples (impossible to restrict to post-task sleep');
end
nNeurons = max(spikes(:,2));
r = ripples.timestamps;
pre = InIntervals(r(:,1),sleep(1,:));
post = InIntervals(r(:,1),sleep(2,:));
45
max(spikes(:,2))
write = true
raly
[fList,pList] = matlab.codetools.requiredFilesAndProducts('GLMgainMean.m');
edit Scramble
edit nanmean
which mean
cd C:\Program Files\MATLAB\R2021b\toolbox\matlab\datafun\
cd C:\Program Files\MATLAB\R2021b\toolbox\matlab\datafun
cd('C:\Program Files\MATLAB\R2021b\toolbox\matlab\datafun')
raly
which nanmean
cd C:
which SWRpipeline
cd C:\Users\Cornell\Documents\GitHub\neurocode\
edit ReconstructPosition
edit mean
mean([5 nan 6])
mean([5 nan 6],'omitnan')
edit median
raly
q = rand(1000,2);
kkeyboard
q = rand(1000,2);
[gain,shGain,predictions,er,sh,w] = GLMgain(q,sum(q,2)+rand(1000,1));
qq = sum(q,2)+rand(1000,1);
[gain,shGain,predictions,er,sh,w] = GLMgain(q,sum(q,2)+rand(1000,1));
figure; plot(prediction,qq,'.');
figure; plot(predictions,qq,'.');
[gain,shGain,predictions,er,sh,w] = GLMgain(q,sum(q,2)+rand(1000,1),'mode','median');
plot(tl(in) - t(i),double(ripPower0(in)));
plot(predictions,qq,'.');
[gain,shGain,predictions,er,sh,w] = GLMgain(q,q,'mode','median');
[gain,shGain,predictions,er,sh,w] = GLMgain(q,qq,'mode','median');
plot(predictions,qq,'.');
[gain,shGain,predictions,er,sh,w] = GLMgain(q,qq,'mode','mean');
plot(predictions,qq,'.');
rng(0); [gain,shGain,predictions,er,sh,w] = GLMgain(q,qq,'mode','mean');
plot(predictions,qq,'.');
rng(0); [gain,shGain,predictions,er,sh,w] = GLMgain(q,qq,'mode','median');
plot(predictions,qq,'.');
rng(0); [gain,shGain,predictions,er,sh,w] = GLMgain(q,qq,'mode','median');
rng(0); [gain,shGain,predictions,er,sh,w] = GLMgain(q,qq,'mode','mean');
clf
hist(shGain)
q = exp(rand(1000,2));
rng(0); [gain,shGain,predictions,er,sh,w] = GLMgain(q,qq,'mode','mean');
rng(0); [gain,shGain,predictions,er,sh,w] = GLMgain(q,qq,'mode','median');
q = exp(rand(1000,2)); qq = sum(q,2)+rand(1000,1);
rng(0); [gain,shGain,predictions,er,sh,w] = GLMgain(q,qq,'mode','median');
rng(0); [gain,shGain,predictions,er,sh,w] = GLMgain(q,qq,'mode','mean');
q = exp(rand(1000,2)); qq = exp(sum(q,2)+rand(1000,1));
rng(0); [gain,shGain,predictions,er,sh,w] = GLMgain(q,qq,'mode','mean');
rng(0); [gain,shGain,predictions,er,sh,w] = GLMgain(q,qq,'mode','median');
rng(0); [gain,shGain,predictions,er,sh,w] = GLMgain(q,qq,'mode','mean');
plot(predictions,qq,'.');
rng(0); [gain,shGain,predictions,er,sh,w] = GLMgain(q,qq,'mode','median');
plot(predictions,qq,'.');
edit ReconstructPosition.m
average
which average
rng(0); [gain,gainsShuffledData,predictions,errors,errorsShuffledData,weights] = GLMgain(q,qq,'mode','median');
p = mean(gainsShuffledData<gain)
p = mean(gainsShuffledData>=gain)
open erros
open errors
p = mean(errors>=errorsShuffledData)
edit glmfit
rng(0); [gain,shGain,predictions,er,sh,w] = GLMgain(q,qq,'mode','mean');
rng(0); [gain,shGain,predictions,er,sh,w] = GLMgain(q,qq,'mode','mean','link','identity');
open er
rng(0); [gain,shGain,predictions,er,sh,w] = GLMgain(q,qq,'mode','mean','link','identity','dist','normal');
v1code
cd(basepath);
which GLMgain
edit ReconstructPosition.m
edit ActivityTemplatesICA.m
edit ActivityTemplates
fmat
edit getLFP
isfolder
isfolder(pwd)
kkeyboard; p.Results.basepath = pwd; p.Results.rejectChannels = [1]
p.Results.basepath = pwd; p.Results.rejectChannels = [1]
field=fieldnames(p.Results)
fields=fieldnames(p.Results)
open p
cellfun(@(x,y) assignin('base', x, y), fields, p.Results)
struct2cell
struct2cell(p.Results)
cellfun(@(x,y) assignin('base', x, y), fields, struct2cell(p.Results))
struct2cell(p.Results)
[fields, struct2cell(p.Results)]
assignin('base','a',1)
which a
assignin([],'a',1)
assignin('caller','a',1)
which a
assignin('caller',a,1)
raly
SaveEEG
struct2cell(p.Results)
fields = fieldnames(p.Results);
cellfun(@(x,y) assignin('base', x, y), fields, struct2cell(p.Results))
clear basepath
SaveEEG
clear basepath
SaveEEG
which basenameFromBasepath
load(fullfile(basepath,[basename '.session.mat'])
load(fullfile(basepath,[basename '.session.mat']))
basename = basenameFromBasepath(basepath);
load(fullfile(basepath,[basename '.session.mat'])
load(fullfile(basepath,[basename '.session.mat'])
load(fullfile(basepath,[basename '.session.mat']))
cd N:\V1test\V1Jean\day17
basepath = pwd
basename = basenameFromBasepath(basepath);
load(fullfile(basepath,[basename '.session.mat'])
load(fullfile(basepath,[basename '.session.mat']))
nChannels = session.extracellular.nChannels
okChannels = find(~ismember((1:nChannels)',rejectChannels(:)));
SaveEEG
hist(spikes(:,1))
hist(r(:,1))
clf
hist(r(:,1))
hist(spikes(:,1))
[h,ht] = PETH(spikes(:,1),r(:,1));
PlotColorMap(Shrink(h,100,1),'x',ht);
write
clear all
cd M:\home\raly\results\V1
v1code
456
write = true
cd N:\V1test\V1Jean\day12
edit GLMgain
clear all
clc
[fList,pList] = matlab.codetools.requiredFilesAndProducts('script_Can_replay.m');
fList = fList';
open fList
which FindReplayScore
which WeightedCorrCirc
which Unshift
raly
fmat
which FindBusrst
which FindBusrsts
which FindBursts
which unshiftdata
which kruskalbar
which nansmooth
which repelem
fList(57) = [];
which clim
which PlotIntervals
raly
mkdir private
for i=49:71, [~,this] = fileparts(fList{i});
copyfile(fList{i},['M:\home\raly\Documents\code\new\private\' this '.m']);
end
edit script_Can_replay.m
45
which ReconstructPosition
fmat
2
PlotXY(pos1,'.');
hopld all
hold all
PlotXY(pos2,'r.');
PlotIntervals(run);
pos0 = sortrows([behavior.positionTrials{1}; behavior.positionTrials{2}(:,1) 2-behavior.positionTrials{2}(:,2)]);
pos1 = Restrict(behavior.positionTrials{1},run);
pos2 = Restrict(behavior.positionTrials{2},run); pos2(:,2) = 1-pos2(:,2);
% the first direction is stim (for labels in later analyses)
if sum(InIntervals(pos2,stim)) > sum(InIntervals(pos1,stim)) % otherwise, swap them
pos1 = pos2; pos2 = Restrict(behavior.positionTrials{1},run);  pos2(:,2) = 1-pos2(:,2);
end
PlotXY(pos1,'.');
clf
PlotXY(pos1,'.');
pos0 = sortrows([behavior.positionTrials{1}; behavior.positionTrials{2}(:,1) 2-behavior.positionTrials{2}(:,2)]);
PlotXY(pos0)
PlotXY(pos0,'.')
xlim([8310 8320]);
pos1 = Restrict(behavior.positionTrials{1},run);
pos2 = Restrict(behavior.positionTrials{2},run); pos2(:,2) = 1-pos2(:,2);
hold all
PlotXY(pos1,'k.');
PlotXY(pos2,'r.');
PlotIntervals(run);
PlotIntervals(sleep,'color','b');
xlim([8310 8320]+[-1 1]*50);
PlotXY(pos2,'r.');
pos0 = sortrows([behavior.positionTrials{1}; behavior.positionTrials{2}(:,1) 2-behavior.positionTrials{2}(:,2)]);
clf
PlotXY(pos0,'k.');
hold all
PlotXY(behavior.positionTrials{1},'r.');
PlotXY(behavior.positionTrials{2},'b.');
mimimum = min([min(behavior.positionTrials{1}(:,2)) min(behavior.positionTrials{2}(:,2)])
mimimum = min([min(behavior.positionTrials{1}(:,2)) min(behavior.positionTrials{2}(:,2))])
mimimum = min([min(behavior.positionTrials{1}(:,2)) min(behavior.positionTrials{2}(:,2))]);
maximum = max([max(behavior.positionTrials{1}(:,2)) max(behavior.positionTrials{2}(:,2))]);
behavior.positionTrials{1}(:,2) = (behavior.positionTrials{1}(:,2) - minimum)./(maximum-minimum);
behavior.positionTrials{2}(:,2) = (behavior.positionTrials{2}(:,2) - minimum)./(maximum-minimum);
minimum = min([min(behavior.positionTrials{1}(:,2)) min(behavior.positionTrials{2}(:,2))]);
maximum = max([max(behavior.positionTrials{1}(:,2)) max(behavior.positionTrials{2}(:,2))]);
behavior.positionTrials{1}(:,2) = (behavior.positionTrials{1}(:,2) - minimum)./(maximum-minimum);
behavior.positionTrials{2}(:,2) = (behavior.positionTrials{2}(:,2) - minimum)./(maximum-minimum);
PlotXY(behavior.positionTrials{1},'r.');
hold all
PlotXY(behavior.positionTrials{2},'b.');
pos0 = sortrows([behavior.positionTrials{1}; behavior.positionTrials{2}(:,1) 2-behavior.positionTrials{2}(:,2)]);
clf
PlotXY(pos0,'k.');
pos1 = Restrict(behavior.positionTrials{1},run);
pos2 = Restrict(behavior.positionTrials{2},run); pos2(:,2) = 1-pos2(:,2);
hold on
PlotXY(pos,'b.');
45
PlotXY(pos,'b.');
PlotXY(pos1,'b.');
hold all
PlotXY(pos2,'r.');
plot(behavior.position.x,behavior.position.y);
plot(behavior.position.timestamps,behavior.position.y);
plot(behavior.timestamps,behavior.position.y);
hold all
plot(behavior.timestamps,behavior.position.linearized);
PlotXY(behavior.positionTrials{1});
plot(behavior.positionTrials{1}(:,1),behavior.positionTrials{1}(:,2)*800);
plot(behavior.positionTrials{1}(:,1),behavior.positionTrials{1}(:,2)*800,'.');
plot(behavior.positionTrials{2}(:,1),behavior.positionTrials{2}(:,2)*800,'.');
plot(behavior.positionTrials{2}(:,1),behavior.positionTrials{2}(:,2)*800,'r.');
ok = ~isnan(behavior.position.x(:)) & ~isnan(behavior.position.y(:));
speed = LinearVelocity([behavior.timestamps(ok)' behavior.position.x(ok)' behavior.position.y(ok)'],5);
run = behavior.timestamps(FindInterval(speed(:,2)>100));
run = ConsolidateIntervals(run,'epsilon',0.01);
[in,w] = InIntervals(behavior.timestamps(:),run);
peak = Accumulate(w(in),behavior.speed(in)','mode','max');
x = xlim;
figure; PlotXY(speed); xlim(x);
SideAxes(gca,'bottom',0.5)
PlotXY(speed);
PlotXY(speed,'.');
plot(behavior.timestamps,behavior.position.y,'k.');
plot(behavior.timestamps,behavior.position.x,'y.');
xlim(x);
PlotIntervals(Restrict(run,xlim));
ok = ~isnan(behavior.position.x(:)) & ~isnan(behavior.position.y(:)); t = behavior.timestamps(ok);
speed = LinearVelocity([behavior.timestamps(ok)' behavior.position.x(ok)' behavior.position.y(ok)'],5);
run = t(FindInterval(speed(:,2)>100));
run = ConsolidateIntervals(run,'epsilon',0.01);
[in,w] = InIntervals(behavior.timestamps(:),run);
peak = Accumulate(w(in),behavior.speed(in)','mode','max');
PlotIntervals(Restrict(run,xlim),'color','r');
figure; PlotXY(pos2,'r.');
PlotIntervals(stim,'color','k');
PlotXY(pos1,'b.');
x = xlim;
figure; plot(behavior.timestamps,behavior.position.linearized,'.');
xlim(x);
PlotIntervals(stim,'color','y');
clf
PlotIntervals(stim,'color','y');
PlotXY(pos1,'b.');
PlotIntervals(stim,'color','y');
cd M:\home\raly\results\OML
2
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
SaveFig(fullfile('M:\home\raly\results\OML',[sessionID '-JointReplayScoreSlope-On-Off']))
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
kruskalbar(gg(ok,1),g(ok,end)-2)
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
kruskalbar(gg(ok,1),g(ok,end)-2)
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
kruskalbar(gg(ok,1),g(ok,end)-2)
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
SaveFig(fullfile('M:\home\raly\results\OML',[sessionID '-Replay-barplots']))
cd N:\V1test\V1Jean\day12
clear all
script_removeAO52_datNoise
basepath = pwd
basename = basenameFromBasepath(basepath);
datFile = [basepath,filesep, basename, '.dat'];
m = memmapfile(datFile, 'Format','int16','Writable',true);
data = reshape(m.data,64,[]);
nSamples = size(data,2);
okChannels = ~ismember((1:size(data,1))',rejectChannels);
sum(okChannels)
rejectChannels = 1+[22 13 15 24 47 19 51 49 57 55 43 54 39 50 52 41 48 37]; % for Jean\day54
okChannels = ~ismember((1:size(data,1))',rejectChannels);
sum(okChannels)
signal = mean(data(okChannels,:))';
datestr((datenum(clock)))
%-- 4/4/2022 4:34 PM --%
cd N:\OJRproject\OJR48\day7
raly
basepath = pwd
[parentFolder,basename] = fileparts(basepath);
[~,projectName] = fileparts(parentFolder);
sessionID = [projectName '_' basename];
cd(basepath);
load(fullfile(basepath,[basename '.session.mat']),'session');
channel = 48;
lfp = GetAyaLFP(channel);
display(['loaded! @' num2str(toc)]);
[clean,~,badIntervals] = CleanLFP(lfp,'thresholds',[6 10]);
PlotXY(lfp);
PlotIntervals(badIntervals);
figure; PlotXY(clean);
[clean,~,badIntervals] = CleanLFP(lfp,'thresholds',[6 10],'manual',true);
clf
plot(t,d);
threshold2 = 1.5;
dbcont
close all
PlotXY(clean);
deltas0 = FindDeltaPeaks(clean);
hist(deltas0(:,5)-deltas0(:,6),1000);
deltas = deltas0(deltas0(:,5)-deltas0(:,6)>3,:);
deltaWaves.timestamps = deltas(:,[1 3]); deltaWaves.peaks = deltas(:,2);  deltaWaves.peakNormedPower = deltas(:,5); deltaWaves.detectorName = ['channel ' num2str(channel) '(+1), CleanLFP, FindDeltaPeaks, peak-trough>3'];
save(fullfile(basepath,[basename '.deltaWaves.events.mat']),'deltaWaves');
ripples = DetectSWR(1+[20 17],'basepath',basepath,'saveMat',true,'forceDetect',true,'check',true,'useSPW',false);
% First, remove points that occur in noisy lfp periods
tl = (1:length(lfp))'/1250; % timestamps for lfp
t = featureTs/1250; % timestamps for candidate ripple events
bad = false(size(t));
[clean,~,badIntervals] = CleanLFP([tl,double(lfp(:,1))],'thresholds',[6 Inf],'manual',true);
threshold1 = 4;
clf
plot(t,d);
threshold2 = 2;
dbcont
bad = bad | InIntervals(t,badIntervals);
Portion(bad)
idx1 = idx1 & ~bad; idx2 = idx2 & ~bad; % remove bad points from "idx1" and "idx2"
figure(1); % plot a clean figure without these points and check indivudual ripples
clf
plot(swDiffAll(idx1 | idx2),ripPowerAll(idx1 | idx2),'k.','markersize',1); hold on
matrix = [swDiffAll ripPowerAll];
z = bsxfun(@rdivide,bsxfun(@minus,matrix,mean(matrix(~bad,:))),std(matrix(~bad,:)));
scores = nan(size(ripPowerAll,1),1);
while true % click "Ctrl+C" to exit when you feel you're done
figure(1);
colors = Bright(1000);
[x,y,button] = ginput(1);
[~,i] = min(abs(x-swDiffAll.*(1-bad))+abs(y-ripPowerAll.*(1-bad)));
t = featureTs/1250;
tl = (1:length(lfp))'/1250;
figure(2);
clf; interval = t(i) + [-1 1]*0.5; in = tl>interval(1) & tl<interval(end);
subplot(3,1,1);
plot(tl(in) - t(i),lfp(in,:)); title(num2str(t(i)));
hold on; plot(bsxfun(@times,ones(1,2),(Restrict(t,interval) - t(i))),ylim,'k--');
legend('ripple channel','SW channel','noise channel');
subplot(3,1,2);
plot(tl(in) - t(i),double(swDiff(in)));
hold on; plot(bsxfun(@times,ones(1,2),(Restrict(t,interval) - t(i))),ylim,'k--');
subplot(3,1,3);
plot(tl(in) - t(i),double(ripPower0(in)));
hold on; plot(bsxfun(@times,ones(1,2),(Restrict(t,interval) - t(i))),ylim,'k--');
%     x = input( prompt )
commandwindow % select the command window
score = str2double(input('good(1)/bad(0)?','s'));
scores(i) = score; % save scores to optionally save
figure(1);
xlims = xlim; ylims = ylim; clf
% extrapolate the score of each ripple based on the score of the nearest scored ripple
scored = find(~isnan(scores));
distances = sqrt(bsxfun(@minus,z(:,1),z(scored,1)').^2 + bsxfun(@minus,z(:,2),z(~isnan(scores),2)').^2);
[~,nearestNeighbour] = min(distances,[],2);
estimated = scores(scored(nearestNeighbour));
selected = estimated>0.1;
idx1 = selected & ~bad; % final ripples
idx2 = ~selected & ~bad; % non-ripples (yet free from noise as well)
scatter(matrix(~bad,1),matrix(~bad,2),1,estimated(~bad),'filled'); colormap(Bright); set(gca,'CLim',[0 1])
xlim(xlims); ylim(ylims);
hold on;
scatter(swDiffAll(scored),ripPowerAll(scored),20,scores(scored)); scatter(swDiffAll(scored),ripPowerAll(scored),15,scores(scored));
set(get(colorbar,'YLabel'),'String','Estimated ripple score');
xlabel('Sharp wave depth'); ylabel('Ripple power');
title({['Manual scoring for ' basepath],['called with channels: ' num2str(Channels) ' (ripple, sharp wave, and (optionally) noise)']});
if rem(i,10)==0
display('saving...');
save(fullfile(basepath,'DetectSWR_manual_scoring.mat'),'t','swDiffAll','ripPowerAll','bad','scores');
end
end
0
1
0
1
0
1
0
1
0
1
0
1
0
1
0
1
0
1
0
1
0
1
0
1
0
% Save your progress!
display('saving...');
save(fullfile(basepath,'DetectSWR_manual_scoring.mat'),'t','swDiffAll','ripPowerAll','idx1','idx2','bad','selected','scores');
saveas(figure(1),fullfile(basepath,'DetectSWR_manual_scoring.fig'));
scored = find(~isnan(scores));
distances = sqrt(bsxfun(@minus,z(:,1),z(scored,1)').^2 + bsxfun(@minus,z(:,2),z(~isnan(scores),2)').^2);
[~,nearestNeighbour] = min(distances,[],2);
estimated = scores(scored(nearestNeighbour));
selected = estimated>0.1;
idx1 = selected & ~bad; % final ripples
idx2 = ~selected & ~bad; % non-ri
clf
plot(swDiffAll(idx2),ripPowerAll(idx2),'k.','markersize',1); hold on
plot(swDiffAll(idx1),ripPowerAll(idx1),'r.','markersize',1); legend('non-ripples','ripples');
idx2(ripPowerAll>60) = 0;
clf
plot(swDiffAll(idx2),ripPowerAll(idx2),'k.','markersize',1); hold on
plot(swDiffAll(idx1),ripPowerAll(idx1),'r.','markersize',1); legend('non-ripples','ripples');
idx2(swDiffAll>400) = 0;
clf
plot(swDiffAll(idx2),ripPowerAll(idx2),'k.','markersize',1); hold on
plot(swDiffAll(idx1),ripPowerAll(idx1),'r.','markersize',1); legend('non-ripples','ripples');
save(fullfile(basepath,'DetectSWR_manual_scoring.mat'),'t','swDiffAll','ripPowerAll','idx1','idx2','bad','selected','scores');
dbcont
r = ripples.timestamps;
PETH(r(:,1),deltas(:,2));
PETH(r(:,1),deltas(:,2),'nBins',501);
PETH(r(:,1),deltas(:,1),'nBins',501);
PETH(r(:,1),deltas(:,3),'nBins',501);
PETH(r(:,1),deltas(:,1),'nBins',501);
PETH(r(:,1),deltas(:,2),'nBins',501);
PETH(r(:,1),deltas(:,3),'nBins',501);
PETH(r(:,1),deltas(:,2),'nBins',501);
PETH(mean(r,[],2),deltas(:,2),'nBins',501);
PETH(mean(r,2),deltas(:,2),'nBins',501);
[h,ht] = PETH(mean(r,2),deltas(:,2),'nBins',501);
semplot(ht,h,'k',1);
clf
h = h/diff(ht(1:2));
semplot(ht,h,'k',1);
ylabel('ripple rate (Hz)')
xlabel('time from delta wave peaks (s)');
set(gca,'tickdir','out','fontsize
set(gca,'tickdir','out','fontsize',15);
SaveCustomEvents('deltas.del.evt',deltas(:,1:3),{'deltas start','delta peak','deltas stop'});
% example at 1642.4
clear all
basepath = pwd
[parentFolder,basename] = fileparts(basepath);
[~,projectName] = fileparts(parentFolder);
sessionID = [projectName '-' basename];
lfp = GetAyaLFP(13);
[clean,~,badIntervals] = CleanLFP(lfp,'thresholds',[2 2],'manual',true);
clf
pllt(t,d);
clf; plot(t,d); hold all; plot(xlim,ones(1,2)*threshold2,'r--'); ylabel('lfp derivative (z-units'); legend('derivative','threshold2');
threshold2 = 1.5;
dbcont
goodIntervals = SubtractIntervals([0 lfp(end,1)],badIntervals);
ripples = FindRipples(basepath,1+13,'noise',1+7,'saveMat',true,'restrict',goodIntervals);
ripples = FindRipples(basepath,1+13,'noise',1+48,'saveMat',true,'restrict',goodIntervals);
basepath
ripples = FindRipples(basepath,1+13,'noise',1+48,'saveMat',true,'restrict',goodIntervals);
ripples = FindRipples('basepath',basepath,'channel',1+13,'noise',1+48,'saveMat',true,'restrict',goodIntervals);
p
if isstr(varargin{1})  % if first arg is basepath
%     addRequired(p, 'basepath',@isstr);
%     addRequired(p,'channel',@isnumeric)
parse(p,varargin{:})
basename = basenameFromBasepath(p.Results.basepath);
passband = p.Results.passband;
EMGThresh = p.Results.EMGThresh;
lfp = getLFP(p.Results.channel,'basepath',p.Results.basepath,'basename',basename);
signal = bz_Filter(lfp,'filter','butter','passband',passband,'order',3);
timestamps = lfp.timestamps;
basepath = p.Results.basepath;
elseif isnumeric(varargin{1}) % if first arg is filtered LFP
addRequired(p,'lfp',@isnumeric)
addRequired(p,'timestamps',@isnumeric)
parse(p,varargin{:})
passband = p.Results.passband;
EMGThresh = p.Results.EMGThresh;
basepath = p.Results.basepath;
timestamps = p.Results.timestamps;
% package into struct for bz_Filter, so we can return struct
samples.samplingRate = p.Results.SR;
samples.data = double(p.Results.lfp);
samples.timestamps = timestamps;
signal = bz_Filter(samples,'filter','butter','passband',passband,'order',3);
if ~exist('basepath','var')
basepath = pwd;
end
basename = basenameFromBasepath(basepath);
end
ripples = FindRipples('basepath',basepath,'channel',1+13,'noise',1+48,'saveMat',true,'restrict',goodIntervals);
dbcont
ripples = FindRipples('basepath',basepath,'channel',1+13,'noise',1+48,'saveMat',true,'restrict',goodIntervals);
ripples = DetectSWR(1+[9 14],'saveMat',true,'check',true,'useSPW',false,'overwrite',true);
ripples = DetectSWR(1+[9 14],'saveMat',true,'check',true,'useSPW',false,'forceDetect',true);
% First, remove points that occur in noisy lfp periods
tl = (1:length(lfp))'/1250; % timestamps for lfp
t = featureTs/1250; % timestamps for candidate ripple events
bad = false(size(t));
[clean,~,badIntervals] = CleanLFP([tl,double(lfp(:,1))],'thresholds',[2 1.5],'manual',true);
clf; plot(t,d); hold all; plot(xlim,ones(1,2)*threshold2,'r--'); ylabel('lfp derivative (z-units'); legend('derivative','threshold2');
threshold2 = 1;
dbcont
Portion(bad)
bad = bad | InIntervals(t,badIntervals);
Portion(bad)
[clean,~,badIntervals] = CleanLFP([tl,double(lfp(:,2))],'thresholds',[2 1.5],'manual',true);
threshold1 = 3;
clf; plot(t,d); hold all; plot(xlim,ones(1,2)*threshold2,'r--'); ylabel('lfp derivative (z-units'); legend('derivative','threshold2');
dbcont
bad = bad | InIntervals(t,badIntervals);
Portion(bad)
idx1 = idx1 & ~bad; idx2 = idx2 & ~bad; % remove bad points from "idx1" and "idx2"
figure(1); % plot a clean figure without these points and check indivudual ripples
clf
plot(swDiffAll(idx1 | idx2),ripPowerAll(idx1 | idx2),'k.','markersize',1); hold on
matrix = [swDiffAll ripPowerAll];
z = bsxfun(@rdivide,bsxfun(@minus,matrix,mean(matrix(~bad,:))),std(matrix(~bad,:)));
scores = nan(size(ripPowerAll,1),1);
while true % click "Ctrl+C" to exit when you feel you're done
figure(1);
colors = Bright(1000);
[x,y,button] = ginput(1);
[~,i] = min(abs(x-swDiffAll.*(1-bad))+abs(y-ripPowerAll.*(1-bad)));
t = featureTs/1250;
tl = (1:length(lfp))'/1250;
figure(2);
clf; interval = t(i) + [-1 1]*0.5; in = tl>interval(1) & tl<interval(end);
subplot(3,1,1);
plot(tl(in) - t(i),lfp(in,:)); title(num2str(t(i)));
hold on; plot(bsxfun(@times,ones(1,2),(Restrict(t,interval) - t(i))),ylim,'k--');
legend('ripple channel','SW channel','noise channel');
subplot(3,1,2);
plot(tl(in) - t(i),double(swDiff(in)));
hold on; plot(bsxfun(@times,ones(1,2),(Restrict(t,interval) - t(i))),ylim,'k--');
subplot(3,1,3);
plot(tl(in) - t(i),double(ripPower0(in)));
hold on; plot(bsxfun(@times,ones(1,2),(Restrict(t,interval) - t(i))),ylim,'k--');
%     x = input( prompt )
commandwindow % select the command window
score = str2double(input('good(1)/bad(0)?','s'));
scores(i) = score; % save scores to optionally save
figure(1);
xlims = xlim; ylims = ylim; clf
% extrapolate the score of each ripple based on the score of the nearest scored ripple
scored = find(~isnan(scores));
distances = sqrt(bsxfun(@minus,z(:,1),z(scored,1)').^2 + bsxfun(@minus,z(:,2),z(~isnan(scores),2)').^2);
[~,nearestNeighbour] = min(distances,[],2);
estimated = scores(scored(nearestNeighbour));
selected = estimated>0.1;
idx1 = selected & ~bad; % final ripples
idx2 = ~selected & ~bad; % non-ripples (yet free from noise as well)
scatter(matrix(~bad,1),matrix(~bad,2),1,estimated(~bad),'filled'); colormap(Bright); set(gca,'CLim',[0 1])
xlim(xlims); ylim(ylims);
hold on;
scatter(swDiffAll(scored),ripPowerAll(scored),20,scores(scored)); scatter(swDiffAll(scored),ripPowerAll(scored),15,scores(scored));
set(get(colorbar,'YLabel'),'String','Estimated ripple score');
xlabel('Sharp wave depth'); ylabel('Ripple power');
title({['Manual scoring for ' basepath],['called with channels: ' num2str(Channels) ' (ripple, sharp wave, and (optionally) noise)']});
if rem(i,10)==0
display('saving...');
save(fullfile(basepath,'DetectSWR_manual_scoring.mat'),'t','swDiffAll','ripPowerAll','bad','scores');
end
end
0
1
0
1
0
1
0
scores(i) = nan;
while true % click "Ctrl+C" to exit when you feel you're done
figure(1);
colors = Bright(1000);
[x,y,button] = ginput(1);
[~,i] = min(abs(x-swDiffAll.*(1-bad))+abs(y-ripPowerAll.*(1-bad)));
t = featureTs/1250;
tl = (1:length(lfp))'/1250;
figure(2);
clf; interval = t(i) + [-1 1]*0.5; in = tl>interval(1) & tl<interval(end);
subplot(3,1,1);
plot(tl(in) - t(i),lfp(in,:)); title(num2str(t(i)));
hold on; plot(bsxfun(@times,ones(1,2),(Restrict(t,interval) - t(i))),ylim,'k--');
legend('ripple channel','SW channel','noise channel');
subplot(3,1,2);
plot(tl(in) - t(i),double(swDiff(in)));
hold on; plot(bsxfun(@times,ones(1,2),(Restrict(t,interval) - t(i))),ylim,'k--');
subplot(3,1,3);
plot(tl(in) - t(i),double(ripPower0(in)));
hold on; plot(bsxfun(@times,ones(1,2),(Restrict(t,interval) - t(i))),ylim,'k--');
%     x = input( prompt )
commandwindow % select the command window
score = str2double(input('good(1)/bad(0)?','s'));
scores(i) = score; % save scores to optionally save
figure(1);
xlims = xlim; ylims = ylim; clf
% extrapolate the score of each ripple based on the score of the nearest scored ripple
scored = find(~isnan(scores));
distances = sqrt(bsxfun(@minus,z(:,1),z(scored,1)').^2 + bsxfun(@minus,z(:,2),z(~isnan(scores),2)').^2);
[~,nearestNeighbour] = min(distances,[],2);
estimated = scores(scored(nearestNeighbour));
selected = estimated>0.1;
idx1 = selected & ~bad; % final ripples
idx2 = ~selected & ~bad; % non-ripples (yet free from noise as well)
scatter(matrix(~bad,1),matrix(~bad,2),1,estimated(~bad),'filled'); colormap(Bright); set(gca,'CLim',[0 1])
xlim(xlims); ylim(ylims);
hold on;
scatter(swDiffAll(scored),ripPowerAll(scored),20,scores(scored)); scatter(swDiffAll(scored),ripPowerAll(scored),15,scores(scored));
set(get(colorbar,'YLabel'),'String','Estimated ripple score');
xlabel('Sharp wave depth'); ylabel('Ripple power');
title({['Manual scoring for ' basepath],['called with channels: ' num2str(Channels) ' (ripple, sharp wave, and (optionally) noise)']});
if rem(i,10)==0
display('saving...');
save(fullfile(basepath,'DetectSWR_manual_scoring.mat'),'t','swDiffAll','ripPowerAll','bad','scores');
end
end
1
0
1
0
1
0
1
0
idx2(swDiffAll>800) = 0;
idx2(ripPowerAll>60) = 0;
save(fullfile(basepath,'DetectSWR_manual_scoring.mat'),'t','swDiffAll','ripPowerAll','idx1','idx2','bad','selected','scores');
saveas(figure(1),fullfile(basepath,'DetectSWR_manual_scoring.fig'));
clf
plot(swDiffAll(idx2),ripPowerAll(idx2),'k.','markersize',1); hold on
plot(swDiffAll(idx1),ripPowerAll(idx1),'r.','markersize',1); legend('non-ripples','ripples');
dbcont
close all
eegFile = fullfile(basepath,[basename '.eeg']); % This is where I will store the normalized LFP file
lfpFile = fullfile(basepath,[basename '.lfp']);
copyfile(lfpFile,eegFile);
file = memmapfile(eegFile,'Format','int16','Writable',true);
data = reshape(file.Data,64,[]);
newData = bsxfun(@minus,data,m);
m = int16(mean(data(1+[9:20 41:44 47:52],:))); % OJR47
newData = bsxfun(@minus,data,m);
file.data = newData(:);
clear file
ripples = DetectSWR(1+[54 35],'saveMat',true,'check',true,'useSPW',false,'useEEG',true,'forceDetect',true);
dbquit
copyfile(lfpFile,eegFile);
file = memmapfile(eegFile,'Format','int16','Writable',true);
data = reshape(file.Data,64,[]);
newData = bsxfun(@minus,data,m);
okChannels = Unfind(1+[9:20 41:44 47:52],64);
m = int16(mean(data(okChannels,:))); % OJR47
mBad = int16(mean(data(~okChannels,:))); % OJR47
newData = data;
newData(okChannels,:) = bsxfun(@minus,data(okChannels,:),m);
newData(~okChannels,:) = bsxfun(@minus,data(~okChannels,:),m);
file.data = newData(:);
45
ripples = DetectSWR(1+[10 13],'saveMat',true,'check',true,'useSPW',false,'forceDetect',true);
close all
clf
plot(swDiffAll(idx2),ripPowerAll(idx2),'k.','markersize',1); hold on
plot(swDiffAll(idx1),ripPowerAll(idx1),'r.','markersize',1); legend('non-ripples','ripples');
tl = (1:length(lfp))'/1250; % timestamps for lfp
t = featureTs/1250; % timestamps for candidate ripple events
bad = false(size(t));
[clean,~,badIntervals] = CleanLFP([tl,double(lfp(:,2))],'thresholds',[2 1.5],'manual',true);
ripples = DetectSWR(1+[10 13],'saveMat',true,'check',true,'useSPW',false,'forceDetect',true,'useEEG');
ripples = DetectSWR(1+[10 13],'saveMat',true,'check',true,'useSPW',false,'forceDetect',true,'useEEG',true);
% First, remove points that occur in noisy lfp periods
tl = (1:length(lfp))'/1250; % timestamps for lfp
t = featureTs/1250; % timestamps for candidate ripple events
bad = false(size(t));
try % optionally, remove periods of noisy LFP (large deflections)
[clean,~,badIntervals] = CleanLFP([tl,double(lfp(:,2))],'thresholds',[2 1.5],'manual',true);
bad = bad | InIntervals(t,badIntervals);
end
threshold1 = 10;
clf; plot(t,d); hold all; plot(xlim,ones(1,2)*threshold2,'r--'); ylabel('lfp derivative (z-units'); legend('derivative','threshold2');
threshold2 = 5;
dbcont
Portion(bad)
[clean,~,badIntervals] = CleanLFP([tl,double(lfp(:,1))],'thresholds',[2 1.5],'manual',true);
bad = bad | InIntervals(t,badIntervals);
threshold1 = 10;
clf; plot(t,d); hold all; plot(xlim,ones(1,2)*threshold2,'r--'); ylabel('lfp derivative (z-units'); legend('derivative','threshold2');
dbcont
Portion(bad)
idx1 = idx1 & ~bad; idx2 = idx2 & ~bad; % remove bad points from "idx1" and "idx2"
figure(1); % plot a clean figure without these points and check indivudual ripples
clf
plot(swDiffAll(idx1 | idx2),ripPowerAll(idx1 | idx2),'k.','markersize',1); hold on
matrix = [swDiffAll ripPowerAll];
z = bsxfun(@rdivide,bsxfun(@minus,matrix,mean(matrix(~bad,:))),std(matrix(~bad,:)));
scores = nan(size(ripPowerAll,1),1);
while true % click "Ctrl+C" to exit when you feel you're done
figure(1);
colors = Bright(1000);
[x,y,button] = ginput(1);
[~,i] = min(abs(x-swDiffAll.*(1-bad))+abs(y-ripPowerAll.*(1-bad)));
t = featureTs/1250;
tl = (1:length(lfp))'/1250;
figure(2);
clf; interval = t(i) + [-1 1]*0.5; in = tl>interval(1) & tl<interval(end);
subplot(3,1,1);
plot(tl(in) - t(i),lfp(in,:)); title(num2str(t(i)));
hold on; plot(bsxfun(@times,ones(1,2),(Restrict(t,interval) - t(i))),ylim,'k--');
legend('ripple channel','SW channel','noise channel');
subplot(3,1,2);
plot(tl(in) - t(i),double(swDiff(in)));
hold on; plot(bsxfun(@times,ones(1,2),(Restrict(t,interval) - t(i))),ylim,'k--');
subplot(3,1,3);
plot(tl(in) - t(i),double(ripPower0(in)));
hold on; plot(bsxfun(@times,ones(1,2),(Restrict(t,interval) - t(i))),ylim,'k--');
%     x = input( prompt )
commandwindow % select the command window
score = str2double(input('good(1)/bad(0)?','s'));
scores(i) = score; % save scores to optionally save
figure(1);
xlims = xlim; ylims = ylim; clf
% extrapolate the score of each ripple based on the score of the nearest scored ripple
scored = find(~isnan(scores));
distances = sqrt(bsxfun(@minus,z(:,1),z(scored,1)').^2 + bsxfun(@minus,z(:,2),z(~isnan(scores),2)').^2);
[~,nearestNeighbour] = min(distances,[],2);
estimated = scores(scored(nearestNeighbour));
selected = estimated>0.1;
idx1 = selected & ~bad; % final ripples
idx2 = ~selected & ~bad; % non-ripples (yet free from noise as well)
scatter(matrix(~bad,1),matrix(~bad,2),1,estimated(~bad),'filled'); colormap(Bright); set(gca,'CLim',[0 1])
xlim(xlims); ylim(ylims);
hold on;
scatter(swDiffAll(scored),ripPowerAll(scored),20,scores(scored)); scatter(swDiffAll(scored),ripPowerAll(scored),15,scores(scored));
set(get(colorbar,'YLabel'),'String','Estimated ripple score');
xlabel('Sharp wave depth'); ylabel('Ripple power');
title({['Manual scoring for ' basepath],['called with channels: ' num2str(Channels) ' (ripple, sharp wave, and (optionally) noise)']});
if rem(i,10)==0
display('saving...');
save(fullfile(basepath,'DetectSWR_manual_scoring.mat'),'t','swDiffAll','ripPowerAll','bad','scores');
end
end
0
1
0
1
0
1
0
1
idx2(ripPowerAll>60) = 0;
idx2(swDiffAll>800) = 0;
display('saving...');
save(fullfile(basepath,'DetectSWR_manual_scoring.mat'),'t','swDiffAll','ripPowerAll','idx1','idx2','bad','selected','scores');
saveas(figure(1),fullfile(basepath,'DetectSWR_manual_scoring.fig'));
clf
plot(swDiffAll(idx2),ripPowerAll(idx2),'k.','markersize',1); hold on
plot(swDiffAll(idx1),ripPowerAll(idx1),'r.','markersize',1); legend('non-ripples','ripples');
dbcont
channel = 49;
lfp = GetAyaLFP(channel);
lfp = GetAyaEEG(channel);
display(['loaded! @' num2str(toc)]);
[clean,~,badIntervals] = CleanLFP(lfp,'thresholds',[6 10],'manual',true);
clf; plot(t,d); hold all; plot(xlim,ones(1,2)*threshold2,'r--'); ylabel('lfp derivative (z-units'); legend('derivative','threshold2');
threshold2 = 1.5;
dbcont
deltas0 = FindDeltaPeaks(clean);
[h,ht] = PETH(mean(r,2),deltas(:,2),'nBins',501);
r = ripples.timestamps;
[h,ht] = PETH(mean(r,2),deltas(:,2),'nBins',501);
[h,ht] = PETH(mean(r,2),deltas0(:,2),'nBins',501);
PlotColorMap(Shrink(sortby(h,deltas0(:,5)-deltas0(:,6)),100,1),'x',ht);
PlotColorMap(Shrink(sortby(h,deltas0(:,5)-deltas0(:,6)),500,1),'x',ht);
semplot(ht,h(deltas0(:,5)-deltas0(:,6) < 3,:));
semplot(ht,h(deltas0(:,5)-deltas0(:,6) > 3,:));
clf
semplot(ht,h(deltas0(:,5)-deltas0(:,6) > 3,:));
clf
PETH(mean(r,2),deltas0(:,2),'nBins',501);
PETH(mean(r,2),deltas0(:,1),'nBins',501);
PETH(mean(r,2),deltas0(:,2),'nBins',501);
deltas = deltas0(deltas0(:,5)-deltas0(:,6)>3,:); % these th
deltaWaves.timestamps = deltas(:,[1 3]); deltaWaves.peaks = deltas(:,2);  deltaWaves.peakNormedPower = deltas(:,5); deltaWaves.detectorName = ['channel ' num2str(channel) '(+1), CleanLFP, FindDeltaPeaks, peak-trough>3'];
save(fullfile(basepath,[basename '.deltaWaves.events.mat']),'deltaWaves');
SaveCustomEvents('deltas.del.evt',deltas(:,1:3),{'deltas start','delta peak','deltas stop'});
figure; PlotXY(Restrict(clean,[0 1]+926));
PlotHVLines(Restrict(deltas0(:,2),xlim),'k--');
figure; PlotXY(Restrict(clean,[-2 2]+926));
PlotHVLines(Restrict(deltas0(:,2),xlim),'k--');
deltas0 = FindDeltaPeaks(clean);
% Default values
highPeak = 3; % Threshold for filtered signal (number of SDs)
lowPeak = 1;
highTrough = 1.5;
lowTrough = 0;
minDuration = 150; % min time between successive zero crossings (in ms)
maxDuration = 650; % max time between successive zero crossings (in ms)
% Differentiate, filter and z-score signal
filtered = FilterLFP(lfp,'passband',[0 9],'order',8);
z = filtered;
z(:,2) = [diff(filtered(:,2));0];
z = FilterLFP(z,'passband',[0 6],'order',8);
z(:,2) = zscore(z(:,2));
figure; PlotXY(Restrict(z,[-2 2]+926));
hold all; PlotXY(Restrict(filtered,[-2 2]+926));
% Find positions (in # samples) of zero crossings
[up,down] = ZeroCrossings(z);
down = find(down);
up = find(up);
if down(1) < up(1), down(1) = []; end
% List positions (in # samples) of successive up,down,up crossings in an Nx3 matrix
n = length(up);
where = [up(1:n-1) down(1:n-1) up(2:n)];
% List positions but also z-scored amplitudes in an Nx6 matrix (positions then amplitudes)
z = filtered;
z(:,2) = zscore(filtered(:,2));
deltas = z(where,:);
deltas = reshape(deltas,size(where,1),6);
open deltas
PlotHVLines(Restrict(deltas(:,2),xlim),'y--');
% Discard waves that are too long or too short
duration = deltas(:,3) - deltas(:,1);
deltas(duration<minDuration/1000|duration>maxDuration/1000,:) = [];
PlotHVLines(Restrict(deltas(:,2),xlim),'r--');
% Threshold z-scored peak and trough amplitudes
base = deltas(:,4);
peak = deltas(:,5);
trough = deltas(:,6);
PlotHVLines(Restrict(deltas(:,1:3),xlim),'k-');
figure; PlotXY(Restrict(z,[-2 2]+926));
PlotHVLines(Restrict(deltas(:,1:3),xlim),'k-');
% Differentiate, filter and z-score signal
filtered = FilterLFP(lfp,'passband',[0 9],'order',8);
z = filtered;
z(:,2) = [diff(filtered(:,2));0];
z = FilterLFP(z,'passband',[0 6],'order',8);
z(:,2) = zscore(z(:,2));
figure; PlotXY(Restrict(z,[-2 2]+926));
xlim = x
x = xlim
xlim(x)
PlotHVLines(Restrict(deltas(:,1:3),xlim),'k-');
hold all
PlotXY(z2(Restrict(filtered,[-2 2]+926)));
PlotXY(z2(Restrict(filtered,[-2 2]+926)),'k');
z = filtered;
z(:,2) = [diff(filtered(:,2));0];
figure; PlotXY(Restrict(filtered,x));
hold all
PlotXY(Restrict(z,x));
clf
figure; PlotXY(z2(Restrict(filtered,x)));
hold all
PlotXY(z2(Restrict(z,x)));
PlotHVLines(0,'h','k--');
z1 = FilterLFP(z,'passband',[0 6],'order',8);
PlotXY(z2(Restrict(z1,x)));
z1 = FilterLFP(z,'passband',[0 9],'order',8);
PlotXY(z2(Restrict(z1,x)),'y');
% Differentiate, filter and z-score signal
filtered = FilterLFP(lfp,'passband',[0 9],'order',8);
z = filtered;
z(:,2) = [diff(filtered(:,2));0];
z = FilterLFP(z,'passband',[0 9],'order',8);
z(:,2) = zscore(z(:,2));
% Find positions (in # samples) of zero crossings
[up,down] = ZeroCrossings(z);
down = find(down);
up = find(up);
if down(1) < up(1), down(1) = []; end
% List positions (in # samples) of successive up,down,up crossings in an Nx3 matrix
n = length(up);
where = [up(1:n-1) down(1:n-1) up(2:n)];
% List positions but also z-scored amplitudes in an Nx6 matrix (positions then amplitudes)
z = filtered;
z(:,2) = zscore(filtered(:,2));
deltas = z(where,:);
deltas = reshape(deltas,size(where,1),6);
% Discard waves that are too long or too short
duration = deltas(:,3) - deltas(:,1);
deltas(duration<minDuration/1000|duration>maxDuration/1000,:) = [];
% Threshold z-scored peak and trough amplitudes
base = deltas(:,4);
peak = deltas(:,5);
trough = deltas(:,6);
clf
PlotXY(z2(Restrict(filtered,x)),'y');
PlotHVLines(Restrict(deltas(:,1:3),xlim),'k-');
highPeak
lowTrough
lowPeak
highTrough
close all
deltas00 = deltas0;
deltas0 = FindDeltaPeaks(clean);
dbcont
PETH(mean(r,2),deltas0(:,2),'nBins',501);
PETH(mean(r,2),deltas00(:,2),'nBins',501);
hold all
PETH(mean(r,2),deltas0(:,2),'nBins',501);
clf
PETH(mean(r,2),deltas0(:,1),'nBins',501);
PETH(mean(r,2),deltas0(:,3),'nBins',501);
hold all
PETH(mean(r,2),deltas0(:,1),'nBins',501);
PETH(mean(r,2),deltas0(:,2),'nBins',501);
clf
[h,ht] = PETH(mean(r,2),deltas0(:,2),'nBins',501);
size(deltas00)
size(deltas0)
[h,ht] = PETH(mean(r,2),deltas0(:,2),'nBins',501);
PlotColorMap(Shrink(sortby(h,deltas0(:,5)-deltas0(:,6)),500,1),'x',ht);
PlotColorMap(Shrink(sortby(h,deltas0(:,5)-deltas0(:,6)),100,1),'x',ht);
sum(deltas0(:,5)-deltas0(:,6)>3)
sum(deltas0(:,5)-deltas0(:,6)<3)
deltas = deltas0(deltas0(:,5)-deltas0(:,6)>3,:); % these
deltaWaves.timestamps = deltas(:,[1 3]); deltaWaves.peaks = deltas(:,2);  deltaWaves.peakNormedPower = deltas(:,5); deltaWaves.detectorName = ['channel ' num2str(channel) '(+1), CleanLFP, FindDeltaPeaks, peak-trough>3'];
save(fullfile(basepath,[basename '.deltaWaves.events.mat']),'deltaWaves');
SaveCustomEvents('deltas.del.evt',deltas(:,1:3),{'deltas start','delta peak','deltas stop'});
channel
lfp = GetAyaLFP(channel);
deltas0EEG = deltas0;
[clean,~,badIntervals] = CleanLFP(lfp,'thresholds',[6 10],'manual',true);
threshold1 = 3;
clf; plot(t,d); hold all; plot(xlim,ones(1,2)*threshold2,'r--'); ylabel('lfp derivative (z-units'); legend('derivative','threshold2');
threshold2 = 0.75;
dbcont
deltas0 = FindDeltaPeaks(clean);
PETH(mean(r,2),deltas0(:,2),'nBins',501);
size(deltas0)
PlotXY(clean);
deltas0 = FindDeltaPeaks(clean);
x = [-2 2]+2344;
PlotXY(Restrict(filtered,x));
deltas0 = FindDeltaPeaks(clean);
filtered = FilterLFP(lfp,'passband',[0 9]);
deltas0 = FindDeltaPeaks(clean);
x = [-2 2]+2344;
hold all; PlotXY(Restrict(filtered,x));
deltas0 = FindDeltaPeaks(clean);
x = [-2 2]+2344;
PlotXY(Restrict(filtered,x));
deltas0 = FindDeltaPeaks(clean);
x = [-2 2]+2344;
PlotXY(Restrict(filtered,x));
dbcont
size(deltas0)
clf
PETH(mean(r,2),deltas0(:,2),'nBins',501);
deltas0 = FindDeltaPeaks(clean);
x = [-2 2]+2344;
PlotXY(Restrict(filtered,x));
PlotHVLines(Restrict(deltas(:,1:3),xlim),'k-');
find(InIntervals(deltas(:,2),xlim))
deltas(7429,2)-2344
deltas(7430,2)-2344
deltas(7430,4:end)
x0 = x; x= xlim;
clf
PlotXY(Restrict(z,x));
highTrough
lowPeak
base = deltas(:,4);
peak = deltas(:,5);
trough = deltas(:,6);
case1 = peak > highPeak & trough <= -lowTrough;
case2 = peak >= lowPeak & trough < -highTrough;
peakabovebase = peak > 0.15 + base;
ok = (case1|case2) & peakabovebase;
ok(7430)
x = [-1 1]+669;
PlotXY(Restrict(z,x));
PlotHVLines(Restrict(deltas(:,1:3),xlim),'k-');
Restrict(deltas(:,1:3),xlim)
PlotHVLines(Restrict(deltas(:,2),xlim),'k-');
FindClosest(deltas(:,2),669.07)
deltas(1799,4:end)
ok(1799)
x = [-1 1]+818.7;
c;f
clf
PlotXY(Restrict(z,x));
PlotHVLines(Restrict(deltas(:,2),xlim),'k-');
FindClosest(deltas(:,2),818.7)
deltas(2298,4:end)
ok(2298)
FindClosest(deltas(:,2),819.4)
ok(2300)
PlotHVLines(Restrict(deltas(2300,2),xlim),'r-');
deltas(2300,4:end)
baseptaht
basepath
basepath = pwdf
basepath = pwd
basename = basenameFromBasepath(basepath);
45
2
45
data = reshape(file.Data,128,[]);
% okChannels = Unfind(1+[9:20 41:44 47:52],64);
% okChannels = true(64,1);
okChannels = ~Unfind(1+[15 5 6 7 64 71 96 114 113 117 118 119 103 120 121 122],128);
m = int16(mean(data(okChannels,:))); % OJR47
mBad = int16(mean(data(~okChannels,:))); % OJR47
newData = data;
newData(okChannels,:) = bsxfun(@minus,data(okChannels,:),m);
newData(~okChannels,:) = bsxfun(@minus,data(~okChannels,:),mBad);
file.data = newData(:);
clear file
2
ripples = DetectSWR(1+[3 0],'saveMat',true,'check',true,'useSPW',false,'useEEG',true,'forceDetect',true);
% First, remove points that occur in noisy lfp periods
tl = (1:length(lfp))'/1250; % timestamps for lfp
t = featureTs/1250; % timestamps for candidate ripple events
bad = false(size(t));
try % optionally, remove periods of noisy LFP (large deflections)
[clean,~,badIntervals] = CleanLFP([tl,double(lfp(:,1))],'thresholds',[2 1.5],'manual',true);
bad = bad | InIntervals(t,badIntervals);
end
threshold1= 10;
clf; plot(t,d); hold all; plot(xlim,ones(1,2)*threshold2,'r--'); ylabel('lfp derivative (z-units'); legend('derivative','threshold2');
dbcont
Portion(bad)
[clean,~,badIntervals] = CleanLFP([tl,double(lfp(:,2))],'thresholds',[2 1.5],'manual',true);
threshold1 = 6;
clf; plot(t,d); hold all; plot(xlim,ones(1,2)*threshold2,'r--'); ylabel('lfp derivative (z-units'); legend('derivative','threshold2');
threshold2 =2;
dbcont
bad = bad | InIntervals(t,badIntervals);
Portion(bad)
idx1 = idx1 & ~bad; idx2 = idx2 & ~bad; % remove bad points from "idx1" and "idx2"
figure(1); % plot a clean figure without these points and check indivudual ripples
clf
plot(swDiffAll(idx1 | idx2),ripPowerAll(idx1 | idx2),'k.','markersize',1); hold on
matrix = [swDiffAll ripPowerAll];
z = bsxfun(@rdivide,bsxfun(@minus,matrix,mean(matrix(~bad,:))),std(matrix(~bad,:)));
scores = nan(size(ripPowerAll,1),1);
while true % click "Ctrl+C" to exit when you feel you're done
figure(1);
colors = Bright(1000);
[x,y,button] = ginput(1);
[~,i] = min(abs(x-swDiffAll.*(1-bad))+abs(y-ripPowerAll.*(1-bad)));
t = featureTs/1250;
tl = (1:length(lfp))'/1250;
figure(2);
clf; interval = t(i) + [-1 1]*0.5; in = tl>interval(1) & tl<interval(end);
subplot(3,1,1);
plot(tl(in) - t(i),lfp(in,:)); title(num2str(t(i)));
hold on; plot(bsxfun(@times,ones(1,2),(Restrict(t,interval) - t(i))),ylim,'k--');
legend('ripple channel','SW channel','noise channel');
subplot(3,1,2);
plot(tl(in) - t(i),double(swDiff(in)));
hold on; plot(bsxfun(@times,ones(1,2),(Restrict(t,interval) - t(i))),ylim,'k--');
subplot(3,1,3);
plot(tl(in) - t(i),double(ripPower0(in)));
hold on; plot(bsxfun(@times,ones(1,2),(Restrict(t,interval) - t(i))),ylim,'k--');
%     x = input( prompt )
commandwindow % select the command window
score = str2double(input('good(1)/bad(0)?','s'));
scores(i) = score; % save scores to optionally save
figure(1);
xlims = xlim; ylims = ylim; clf
% extrapolate the score of each ripple based on the score of the nearest scored ripple
scored = find(~isnan(scores));
distances = sqrt(bsxfun(@minus,z(:,1),z(scored,1)').^2 + bsxfun(@minus,z(:,2),z(~isnan(scores),2)').^2);
[~,nearestNeighbour] = min(distances,[],2);
estimated = scores(scored(nearestNeighbour));
selected = estimated>0.1;
idx1 = selected & ~bad; % final ripples
idx2 = ~selected & ~bad; % non-ripples (yet free from noise as well)
scatter(matrix(~bad,1),matrix(~bad,2),1,estimated(~bad),'filled'); colormap(Bright); set(gca,'CLim',[0 1])
xlim(xlims); ylim(ylims);
hold on;
scatter(swDiffAll(scored),ripPowerAll(scored),20,scores(scored)); scatter(swDiffAll(scored),ripPowerAll(scored),15,scores(scored));
set(get(colorbar,'YLabel'),'String','Estimated ripple score');
xlabel('Sharp wave depth'); ylabel('Ripple power');
title({['Manual scoring for ' basepath],['called with channels: ' num2str(Channels) ' (ripple, sharp wave, and (optionally) noise)']});
if rem(i,10)==0
display('saving...');
save(fullfile(basepath,'DetectSWR_manual_scoring.mat'),'t','swDiffAll','ripPowerAll','bad','scores');
end
end
0
1
0
% Save your progress!
display('saving...');
save(fullfile(basepath,'DetectSWR_manual_scoring.mat'),'t','swDiffAll','ripPowerAll','idx1','idx2','bad','selected','scores');
saveas(figure(1),fullfile(basepath,'DetectSWR_manual_scoring.fig'));
dbcont
[parentFolder,basename] = fileparts(basepath);
[~,projectName] = fileparts(parentFolder);
sessionID = [projectName '_' basename];
cd(basepath);
load(fullfile(basepath,[basename '.session.mat']),'session');
load(fullfile(basepath,[basename '.cell_metrics.cellinfo.mat']),'cell_metrics');
nNeurons = length(cell_metrics.spikes.times);
regions = cell_metrics.brainRegion;
regionNames = unique(regions);
regionCell = zeros(length(regions),1);
for i=1:length(regionNames)
regionCell(strcmp(regions,regionNames{i})) = i;
end
spikesCell = cell_metrics.spikes.times';
PFCindex = [find(strcmp(regionNames,'ILA')) find(strcmp(regionNames,'PFC'))];
HPCindex = find(strcmp(regionNames,'CA1'));
cell_metrics = ProcessCellMetrics('session',session,'manualAdjustMonoSyn',false,'removeMetrics','deepSuperficial');
session = sessionTemplate(pwd,'showGUI',true,'remove_folder_date',true);
f = dir('Kilosort*');
spikes = loadSpikes('session',session,'clusteringpath',[f.folder filesep f.name]);
cell_metrics = ProcessCellMetrics('session',session,'manualAdjustMonoSyn',false,'removeMetrics','deepSuperficial');
spikes = Group(spikes.times{:});
clf
r = ripples.timestamps;
PETH(spikes(:,1),r(:,1));
lfp = GetEEG(123);
lfp = GetAyaEEG(123);
[clean,~,badIntervals] = CleanLFP(lfp,'thresholds',[6 10],'manual',true);
raly
clf; plot(t,d); hold all; plot(xlim,ones(1,2)*threshold2,'r--'); ylabel('lfp derivative (z-units'); legend('derivative','threshold2');
threshold2 =2;
clf; plot(t,d); hold all; plot(xlim,ones(1,2)*threshold2,'r--'); ylabel('lfp derivative (z-units'); legend('derivative','threshold2');
threshold2
threshold2 = 1;
clf; plot(t,d); hold all; plot(xlim,ones(1,2)*threshold2,'r--'); ylabel('lfp derivative (z-units'); legend('derivative','threshold2');
threshold2 = 1.5;
dbcont
PETH(r(:,1),deltas(:,2))
deltas = deltas0(deltas0(:,5)-deltas0(:,6)>3,:); % these
PETH(r(:,1),deltas(:,2))
deltas0 = FindDeltaPeaks(clean);
dbcont
open deltas0
PETH(r(:,1),deltas0(:,2))
[h,ht] = PETH(r(:,1),deltas0(:,2))
[h,ht] = PETH(r(:,1),deltas0(:,2));
PlotColorMap(h);
PlotColorMap(Shrink(h,100,1));
PlotColorMap(Shrink(sortby(h,deltas0(:,5)-deltas0(:,6))h,100,1));
PlotColorMap(Shrink(sortby(h,deltas0(:,5)-deltas0(:,6))*h,100,1));
PlotColorMap(Shrink(sortby(h,deltas0(:,5)-deltas0(:,6)),100,1));
cd(basepath);
SaveCustomEvents('deltas.del.evt',deltas0(:,1:3),{'deltas start','delta peak','deltas stop'});
deltas = deltas0(deltas0(:,5)-deltas0(:,6)>3,:); % these
SaveCustomEvents('deltas.del.evt',deltas(:,1:3),{'deltas start','delta peak','deltas stop'});
PETH(lfp,deltas0(:,2))
FindClosest(deltas0(:,2),2893.336)
PETH(lfp,deltas0(673,2))
PETH(lfp,deltas0(673:674,2))
[h,ht] = PETH(lfp,deltas0(673:674,2))
[h,ht] = PETH(lfp,deltas0(673:674,2));
plot(h(1,:))
PlotXY(Restrict(lfp,[0 1]+2893);
PlotXY(Restrict(lfp,[0 1]+2893));
PlotHVLines(deltas0(:,2),'v','k--');
clf
PlotHVLines(deltas0(:,2),'v','k--');
clf
PlotXY(Restrict(lfp,[0 1]+2893));
clf
PlotXY(lfp);
PlotHVLines(deltas0(:,2),'v','k--');
max(lfp(:,1))
size(lfp)
edit GetLFP
lfp2 = GetAyaEEG(123);
lfp2 = GetAyaLFP(123);
session.extracellular.nChannels
load('day7.session.mat')
session = sessionTemplate(pwd,'showGUI',true,'remove_folder_date',true);
cleaar all
clear all
lfp = GetAyaEEG(123);
[clean,~,badIntervals] = CleanLFP(lfp,'thresholds',[6 10],'manual',true);
threshold1 = 4;
clf; plot(t,d); hold all; plot(xlim,ones(1,2)*threshold2,'r--'); ylabel('lfp derivative (z-units'); legend('derivative','threshold2');
threshold2 = 1;
dbcont
clf
PlotXY(clean);
deltas0 = FindDeltaPeaks(clean);
dbcont
cd(basepath);
basepath = pwd
load('day3.ripples.events.mat')
load('day3.spikes.cellinfo.mat')
spikes = Group(spikes.times{:});
r = ripples.timestamps;
PETH(spikes(:,1),r(:,1))
PETH(spikes(:,1),r(:,1)); PlotHVLines(0,'v','r--');
PETH(spikes(:,1),deltas0(:,2)); PlotHVLines(0,'v','r--');
clf
PETH(r(:,1),deltas0(:,2))
[h,ht] = PETH(r(:,1),deltas0(:,2));
PlotColorMap(Shrink(sortby(h,deltas0(:,5)-deltas0(:,6)),100,1));
PlotColorMap(Shrink(sortby(h,deltas0(:,5)-deltas0(:,6)),500,1));
SaveCustomEvents('deltas.del.evt',deltas(:,1:3),{'deltas start','delta peak','deltas stop'});
SaveCustomEvents('deltas.del.evt',deltas0(:,1:3),{'deltas start','delta peak','deltas stop'});a
lfp = GetAyaLFP(123);
[clean,~,badIntervals] = CleanLFP(lfp,'thresholds',[6 10],'manual',true);
threshold1 = 3;
clf; plot(t,d); hold all; plot(xlim,ones(1,2)*threshold2,'r--'); ylabel('lfp derivative (z-units'); legend('derivative','threshold2');
threshold2 = 0.3;
clf; plot(t,d); hold all; plot(xlim,ones(1,2)*threshold2,'r--'); ylabel('lfp derivative (z-units'); legend('derivative','threshold2');
dbcont
clf
PlotXY(clean)
deltas = deltas0(deltas0(:,5)-deltas0(:,6)>3,:); % these thresholds should be manually refined for each session
deltas0 = FindDeltaPeaks(clean);
close all
PETH(r(:,1),deltas0(:,2))
[h,ht] = PETH(r(:,1),deltas0(:,2))
[h,ht] = PETH(r(:,1),deltas0(:,2));
PlotColorMap(Shrink(sortby(h,deltas0(:,5)-deltas0(:,6)),500,1));
PlotColorMap(Shrink(sortby(deltas0(:,5)-deltas0(:,6),deltas0(:,5)-deltas0(:,6)),500,1));
PlotColorMap(Shrink(sortby(h,deltas0(:,5)-deltas0(:,6)),500,1));
PlotColorMap(Shrink(sortby(h,deltas0(:,5)-deltas0(:,6)),1000,1));
semplot(ht,h);
deltas = deltas0(deltas0(:,5)-deltas0(:,6)>3,:); % these thresholds should be manually refined for each session
clf
[h,ht] = PETH(r(:,1),deltas(:,2));
semplot(ht,h);
raly
for i=1:max(spikes(:,2))
[h,ht] = PETH(spikes(spikes(:,2)==i),r(:,1));
m(i,:) = nanmean(h);
end
open m
PlotColorMap(zsocre(m,[],2));
PlotColorMap(zscore(m,[],2));
for i=1:max(spikes(:,2))
[h,ht] = PETH(spikes(spikes(:,2)==i),r(:,1),'durations',[-1 1]*0.4);
m(i,:) = nanmean(h);
end
PlotColorMap(zscore(m,[],2),'x',ht);
[~,mm] = max(Smooth(curves,[0 2],'type','cc'),[],2);
[~,mm] = max(Smooth(m,[0 2],'type','cc'),[],2);
PlotColorMap(sortby(zscore(m,[],2),mm),'x',ht);
load('day3.spikes.cellinfo.mat')
session = sessionTemplate(pwd,'showGUI',true,'remove_folder_date',true);
load('day10.session.mat')
session = sessionTemplate(pwd,'showGUI',true,'remove_folder_date',true);
f = dir('Kilosort*');
spikes = loadSpikes('session',session,'clusteringpath',[f.folder filesep f.name]);
cell_metrics = ProcessCellMetrics('session',session,'manualAdjustMonoSyn',false,'removeMetrics','deepSuperficial');
beta0 = [cell_metrics.general.chanCoords.x(bestChannels(1)),cell_metrics.general.chanCoords.y(bestChannels(1))]; % initial position
load('chanMap.mat')
load('day3.chanCoords.channelInfo.mat', 'chanCoords')
load('day3.chanCoords.channelInfo.mat')
cell_metrics = ProcessCellMetrics('session',session,'manualAdjustMonoSyn',false,'removeMetrics','deepSuperficial');
dbcont
channel_mapping
[pfc,hpc,deltas,r,sleep,session,MergePoints,clickers,sws,basepath] = BatchLoadHPCPFCData(basepath);
SaveCustomEvents('deltas.del.evt',deltas(:,1:3),{'deltas start','delta peak','deltas stop'});a
deltaWaves.timestamps = deltas(:,[1 3]); deltaWaves.peaks = deltas(:,2);  deltaWaves.peakNormedPower = deltas(:,5); deltaWaves.detectorName = ['channel ' num2str(channel) '(+1), CleanLFP, FindDeltaPeaks, peak-trough>3'];
save(fullfile(basepath,[basename '.deltaWaves.events.mat']),'deltaWaves');
X= {[pfc,hpc,deltas,r,sleep,session,MergePoints,clickers,sws,basepath]};
X= {pfc,hpc,deltas,r,sleep,session,MergePoints,clickers,sws,basepath};
basename = basenameFromBasepath(basepath);
save(fullfile(basepath,[basename '.deltaWaves.events.mat']),'deltaWaves');
deltaWaves.timestamps = deltas(:,[1 3]); deltaWaves.peaks = deltas(:,2);  deltaWaves.peakNormedPower = deltas(:,5); deltaWaves.detectorName = ['channel ' num2str(channel) '(+1), CleanLFP, FindDeltaPeaks, peak-trough>3'];
deltas = deltas0(deltas0(:,5)-deltas0(:,6)>3,:); % these thresholds should be manually refined for each session
deltaWaves.timestamps = deltas(:,[1 3]); deltaWaves.peaks = deltas(:,2);  deltaWaves.peakNormedPower = deltas(:,5); deltaWaves.detectorName = ['channel ' num2str(channel) '(+1), CleanLFP, FindDeltaPeaks, peak-trough>3'];
channel = 123;
deltaWaves.timestamps = deltas(:,[1 3]); deltaWaves.peaks = deltas(:,2);  deltaWaves.peakNormedPower = deltas(:,5); deltaWaves.detectorName = ['channel ' num2str(channel) '(+1), CleanLFP, FindDeltaPeaks, peak-trough>3'];
save(fullfile(basepath,[basename '.deltaWaves.events.mat']),'deltaWaves');
[pfc,hpc,deltas,r,sleep,session,MergePoints,clickers,sws,basepath] = BatchLoadHPCPFCData(basepath);
X= {pfc,hpc,deltas,r,sleep,session,MergePoints,clickers,sws,basepath};
XX{1} = X;
isempty(pfc)
channel_mapping
[pfc,hpc,deltas,r,sleep,session,MergePoints,clickers,sws,basepath] = BatchLoadHPCPFCData(basepath);
XX{1} = X;
X= {pfc,hpc,deltas,r,sleep,session,MergePoints,clickers,sws,basepath};
XX{1} = X;
open 'M:\home\raly\results\PFC\reactivation\PFC_responses_to_ripples_colormaps_ordered_by_proximity_to_delta
open M:\home\raly\results\PFC\reactivation\PFC_responses_to_ripples_colormaps_ordered_by_proximity_to_delta
cd M:\home\raly\results\PFC\reactivation\
open cell_metrics
max(pfc(:,2))
max(hpc(:,2))
ylabel('ripples (ordered by proximity to delta wave)');
xlabel('time from ripples (s)');
clabel('PFC firing rate (a.u.)');
set(gca,'tickdir','out','fontsize',15);
ColorMap(gca,[1 1 1],[1 0 0]);
ColorMap(gca,[1 1 1]*0,[1 0 0]);
magma;
raly
colormap(inferno)
colormap(summer)
colormap(autumn)
colormap(bone)
colormap(spring)
colormap(winter)
ColorMap(gca,[1 1 1]*0,[1 0 0]);
ColorMap(gcf,[1 1 1]*0,[1 0 0]);
ColorMap(gcf,[1 1 1]*0,,[1 1 1],[1 0 0]);
ColorMap(gcf,[1 1 1]*0,[1 1 1],[1 0 0]);
ColorMap(gcf,[1 1 1]*0,[1 1 0],[1 0 0]);
ColorMap(gcf,[1 1 1]*1,[1 1 0],[1 0 0]);
ColorMap(gcf,[0 0 1]*1,[1 1 0],[1 0 0]);
ColorMap(gcf,[0 0 1]*1,[0 1 1],[1 0 0]);
ColorMap(gcf,[0 0 1]*1,[0 1 1],[0 1 0]);
ColorMap(gcf,[0 0 1]*1,[0 1 0],[1 0 0]);
ColorMap(gcf,[0 0 1],[1 1 1],[1 0 0]);
ColorMap(gcf,[0 0 1],[1 1 1],[1 0.5 0]);
ColorMap(gcf,[0 0 1],[1 1 1],[1 0.2 0]);
ColorMap(gcf,[0 0 1],[1 1 1],[1 0 0]);
ColorMap(gcf,[0 0 1],[1 1 1],,[1 1 0],[1 0 0]);
ColorMap(gcf,[0 0 1],[1 1 1],[1 1 0],[1 0 0]);
ColorMap(gcf,[0 0 1],[1 1 1],[1 1 1],[1 1 0],[1 0 0]);
ColorMap(gcf,[0 0 1],[1 1 1],[1 0 0]);
ColorMap(gcf,[0 0.2 1],[1 1 1],[1 0 0]);
ColorMap(gcf,[0 0.5 1],[1 1 1],[1 0 0]);
ColorMap(gcf,[0 0.5 1],[1 1 1],[1 0.5 0]);
ColorMap(gcf,[0 0.5 1],[1 1 1],[0.65 0.5 0]);
ColorMap(gcf,[0 0.5 1],[1 1 1],[0.6 0.5 0]);
ColorMap(gcf,[0 0.5 1],[1 1 1],[0.6 0.6 0]);
ColorMap(gcf,[0 0.5 1],[1 1 1],[0.6 1 0]);
ColorMap(gcf,[0 0.5 1],[1 1 1],[1 1 0]);
ColorMap(gcf,[0 0.5 1],[1 1 1],[1 0.5 0]);
ColorMap(gcf,[0 0.5 1],[1 1 1],[1 0.2 0]);
ColorMap(gcf,[0 0.5 1],[1 1 1],[0.8 0.2 0]);
ColorMap(gcf,[0 0.5 1],[1 1 1],[0.8 0.1 0]);
ColorMap(gcf,[0 0.5 1],[1 1 1],[0.9 0.1 0]);
ColorMap(gcf,[0 0.4 1],[1 1 1],[0.9 0.1 0]);
colormap(jet0
colormap(jet)
colormap(hot)
colormap(flipdu(hot))
colormap(flipud(hot))
colormap(hot)
ColorMap(gcf,[0 0.5 1],[1 1 1],[0.9 0.1 0]);
ColorMap(gcf,[0 0 1],[1 1 1],[0.9 0.1 0]);
ColorMap(gcf,[0 0.5 1],[1 1 1],[0.9 0.1 0]);
basepath
close all
raly
clear all
basepath = 'N:\OJRproject\OJR49\day6';
cd C:
[parentFolder,basename] = fileparts(basepath);
[~,projectName] = fileparts(parentFolder);
sessionID = [projectName '_' basename];
cd(basepath);
load(fullfile(basepath,[basename '.session.mat']),'session');
try
load(fullfile(basepath,[basename '.cell_metrics.cellinfo.mat']),'cell_metrics');
nNeurons = length(cell_metrics.spikes.times);
regions = cell_metrics.brainRegion;
regionNames = unique(regions);
regionCell = zeros(length(regions),1);
for i=1:length(regionNames)
regionCell(strcmp(regions,regionNames{i})) = i;
end
spikesCell = cell_metrics.spikes.times';
PFCindex = [find(strcmp(regionNames,'ILA')) find(strcmp(regionNames,'PFC'))];
HPCindex = find(strcmp(regionNames,'CA1'));
pfc = []; if ~isempty(PFCindex), pfc = sortrows(Group(spikesCell{regionCell==PFCindex})); end
hpc = []; if ~isempty(HPCindex), hpc = sortrows(Group(spikesCell{regionCell==HPCindex})); end
neurons = true;
catch
neurons = false; pfc = zeros(0,2); hpc = zeros(0,2); spikesCell = {};
end
spikes = sortrows(Group(spikesCell{:}));
PFCindex
HPCindex
regionNames
channel_mapping
try
channel_mapping
load(fullfile(basepath,[basename '.cell_metrics.cellinfo.mat']),'cell_metrics');
nNeurons = length(cell_metrics.spikes.times);
regions = cell_metrics.brainRegion;
regionNames = unique(regions);
regionCell = zeros(length(regions),1);
for i=1:length(regionNames)
regionCell(strcmp(regions,regionNames{i})) = i;
end
spikesCell = cell_metrics.spikes.times';
PFCindex = [find(strcmp(regionNames,'ILA')) find(strcmp(regionNames,'PFC'))];
HPCindex = find(strcmp(regionNames,'CA1'));
pfc = []; if ~isempty(PFCindex), pfc = sortrows(Group(spikesCell{regionCell==PFCindex})); end
hpc = []; if ~isempty(HPCindex), hpc = sortrows(Group(spikesCell{regionCell==HPCindex})); end
neurons = true;
catch
neurons = false; pfc = zeros(0,2); hpc = zeros(0,2); spikesCell = {};
end
spikes = sortrows(Group(spikesCell{:}));
rippleChannels = 1+[18 30]; % the channel that is lower during the sharp-wave is second
% Detect ripples
if exist(fullfile(basepath,[basename '.ripples.events.mat']),'file')
load(fullfile(basepath,[basename '.ripples.events.mat']),'ripples');
else
load(fullfile(basepath,[basename '.session.mat']),'session');
rippleChannels = swrChannels('basepath',basepath); % or better yet, pick them on neuroscope (don't forget to add 1! rippleChannels.Ripple_Channel = neuroscope_Channel + 1)
ripples = DetectSWR(rippleChannels,'basepath',basepath,'saveMat',true,'forceDetect',true,'check',true,'useSPW',false);
end
% Detect ripples
if exist(fullfile(basepath,[basename '.ripples.events.mat']),'file')
load(fullfile(basepath,[basename '.ripples.events.mat']),'ripples');
else
load(fullfile(basepath,[basename '.session.mat']),'session');
%     rippleChannels = swrChannels('basepath',basepath); % or better yet, pick them on neuroscope (don't forget to add 1! rippleChannels.Ripple_Channel = neuroscope_Channel + 1)
ripples = DetectSWR(rippleChannels,'basepath',basepath,'saveMat',true,'forceDetect',true,'check',true,'useSPW',false);
end
% First, remove points that occur in noisy lfp periods
tl = (1:length(lfp))'/1250; % timestamps for lfp
t = featureTs/1250; % timestamps for candidate ripple events
bad = false(size(t));
try % optionally, remove periods of noisy LFP (large deflections)
[clean,~,badIntervals] = CleanLFP([tl,double(lfp(:,1))],'thresholds',[6 Inf],'manual',true);
bad = bad | InIntervals(t,badIntervals);
end
Channels
edit CleanLFP
clf; plot(t,d); hold all; plot(xlim,ones(1,2)*threshold2,'r--'); ylabel('lfp derivative (z-units'); legend('derivative','threshold2');
disp('Change "threshold2" manually.')
threshold2 = 2;
clf; plot(t,d); hold all; plot(xlim,ones(1,2)*threshold2,'r--'); ylabel('lfp derivative (z-units'); legend('derivative','threshold2');
disp('Change "threshold2" manually.')
threshold2 = 0.5;
clf; plot(t,d); hold all; plot(xlim,ones(1,2)*threshold2,'r--'); ylabel('lfp derivative (z-units'); legend('derivative','threshold2');
disp('Change "threshold2" manually.')
threshold2 = 1;
clf; plot(t,d); hold all; plot(xlim,ones(1,2)*threshold2,'r--'); ylabel('lfp derivative (z-units'); legend('derivative','threshold2');
dbcont
Portion(bad)
idx1 = idx1 & ~bad; idx2 = idx2 & ~bad; % remove bad points from "idx1" and "idx2"
figure(1); % plot a clean figure without these points and check indivudual ripples
clf
plot(swDiffAll(idx1 | idx2),ripPowerAll(idx1 | idx2),'k.','markersize',1); hold on
matrix = [swDiffAll ripPowerAll];
z = bsxfun(@rdivide,bsxfun(@minus,matrix,mean(matrix(~bad,:))),std(matrix(~bad,:)));
scores = nan(size(ripPowerAll,1),1);
while true % click "Ctrl+C" to exit when you feel you're done
figure(1);
colors = Bright(1000);
[x,y,button] = ginput(1);
[~,i] = min(abs(x-swDiffAll.*(1-bad))+abs(y-ripPowerAll.*(1-bad)));
t = featureTs/1250;
tl = (1:length(lfp))'/1250;
figure(2);
clf; interval = t(i) + [-1 1]*0.5; in = tl>interval(1) & tl<interval(end);
subplot(3,1,1);
plot(tl(in) - t(i),lfp(in,:)); title(num2str(t(i)));
hold on; plot(bsxfun(@times,ones(1,2),(Restrict(t,interval) - t(i))),ylim,'k--');
legend('ripple channel','SW channel','noise channel');
subplot(3,1,2);
plot(tl(in) - t(i),double(swDiff(in)));
hold on; plot(bsxfun(@times,ones(1,2),(Restrict(t,interval) - t(i))),ylim,'k--');
subplot(3,1,3);
plot(tl(in) - t(i),double(ripPower0(in)));
hold on; plot(bsxfun(@times,ones(1,2),(Restrict(t,interval) - t(i))),ylim,'k--');
%     x = input( prompt )
commandwindow % select the command window
score = str2double(input('good(1)/bad(0)?','s'));
scores(i) = score; % save scores to optionally save
figure(1);
xlims = xlim; ylims = ylim; clf
% extrapolate the score of each ripple based on the score of the nearest scored ripple
scored = find(~isnan(scores));
distances = sqrt(bsxfun(@minus,z(:,1),z(scored,1)').^2 + bsxfun(@minus,z(:,2),z(~isnan(scores),2)').^2);
[~,nearestNeighbour] = min(distances,[],2);
estimated = scores(scored(nearestNeighbour));
selected = estimated>0.1;
idx1 = selected & ~bad; % final ripples
idx2 = ~selected & ~bad; % non-ripples (yet free from noise as well)
scatter(matrix(~bad,1),matrix(~bad,2),1,estimated(~bad),'filled'); colormap(Bright); set(gca,'CLim',[0 1])
xlim(xlims); ylim(ylims);
hold on;
scatter(swDiffAll(scored),ripPowerAll(scored),20,scores(scored)); scatter(swDiffAll(scored),ripPowerAll(scored),15,scores(scored));
set(get(colorbar,'YLabel'),'String','Estimated ripple score');
xlabel('Sharp wave depth'); ylabel('Ripple power');
title({['Manual scoring for ' basepath],['called with channels: ' num2str(Channels) ' (ripple, sharp wave, and (optionally) noise)']});
if rem(i,10)==0
display('saving...');
save(fullfile(basepath,'DetectSWR_manual_scoring.mat'),'t','swDiffAll','ripPowerAll','bad','scores');
end
end
0
Channels
7/0.152
while true % click "Ctrl+C" to exit when you feel you're done
figure(1);
colors = Bright(1000);
[x,y,button] = ginput(1);
[~,i] = min(abs(x-swDiffAll.*(1-bad))+abs(y-ripPowerAll.*(1-bad)));
t = featureTs/1250;
tl = (1:length(lfp))'/1250;
figure(2);
clf; interval = t(i) + [-1 1]*0.5; in = tl>interval(1) & tl<interval(end);
subplot(3,1,1);
plot(tl(in) - t(i),lfp(in,:)); title(num2str(t(i)));
hold on; plot(bsxfun(@times,ones(1,2),(Restrict(t,interval) - t(i))),ylim,'k--');
legend('ripple channel','SW channel','noise channel');
subplot(3,1,2);
plot(tl(in) - t(i),double(swDiff(in)));
hold on; plot(bsxfun(@times,ones(1,2),(Restrict(t,interval) - t(i))),ylim,'k--');
subplot(3,1,3);
plot(tl(in) - t(i),double(ripPower0(in)));
hold on; plot(bsxfun(@times,ones(1,2),(Restrict(t,interval) - t(i))),ylim,'k--');
%     x = input( prompt )
commandwindow % select the command window
score = str2double(input('good(1)/bad(0)?','s'));
scores(i) = score; % save scores to optionally save
figure(1);
xlims = xlim; ylims = ylim; clf
% extrapolate the score of each ripple based on the score of the nearest scored ripple
scored = find(~isnan(scores));
distances = sqrt(bsxfun(@minus,z(:,1),z(scored,1)').^2 + bsxfun(@minus,z(:,2),z(~isnan(scores),2)').^2);
[~,nearestNeighbour] = min(distances,[],2);
estimated = scores(scored(nearestNeighbour));
selected = estimated>0.1;
idx1 = selected & ~bad; % final ripples
idx2 = ~selected & ~bad; % non-ripples (yet free from noise as well)
scatter(matrix(~bad,1),matrix(~bad,2),1,estimated(~bad),'filled'); colormap(Bright); set(gca,'CLim',[0 1])
xlim(xlims); ylim(ylims);
hold on;
scatter(swDiffAll(scored),ripPowerAll(scored),20,scores(scored)); scatter(swDiffAll(scored),ripPowerAll(scored),15,scores(scored));
set(get(colorbar,'YLabel'),'String','Estimated ripple score');
xlabel('Sharp wave depth'); ylabel('Ripple power');
title({['Manual scoring for ' basepath],['called with channels: ' num2str(Channels) ' (ripple, sharp wave, and (optionally) noise)']});
if rem(i,10)==0
display('saving...');
save(fullfile(basepath,'DetectSWR_manual_scoring.mat'),'t','swDiffAll','ripPowerAll','bad','scores');
end
end
0
while true % click "Ctrl+C" to exit when you feel you're done
figure(1);
colors = Bright(1000);
[x,y,button] = ginput(1);
[~,i] = min(abs(x-swDiffAll.*(1-bad))+abs(y-ripPowerAll.*(1-bad)));
t = featureTs/1250;
tl = (1:length(lfp))'/1250;
figure(2);
clf; interval = t(i) + [-1 1]*0.5; in = tl>interval(1) & tl<interval(end);
subplot(3,1,1);
plot(tl(in) - t(i),lfp(in,:)); title(num2str(t(i)));
hold on; plot(bsxfun(@times,ones(1,2),(Restrict(t,interval) - t(i))),ylim,'k--');
legend('ripple channel','SW channel','noise channel');
subplot(3,1,2);
plot(tl(in) - t(i),double(swDiff(in)));
hold on; plot(bsxfun(@times,ones(1,2),(Restrict(t,interval) - t(i))),ylim,'k--');
subplot(3,1,3);
plot(tl(in) - t(i),double(ripPower0(in)));
hold on; plot(bsxfun(@times,ones(1,2),(Restrict(t,interval) - t(i))),ylim,'k--');
%     x = input( prompt )
commandwindow % select the command window
score = str2double(input('good(1)/bad(0)?','s'));
scores(i) = score; % save scores to optionally save
figure(1);
xlims = xlim; ylims = ylim; clf
% extrapolate the score of each ripple based on the score of the nearest scored ripple
scored = find(~isnan(scores));
distances = sqrt(bsxfun(@minus,z(:,1),z(scored,1)').^2 + bsxfun(@minus,z(:,2),z(~isnan(scores),2)').^2);
[~,nearestNeighbour] = min(distances,[],2);
estimated = scores(scored(nearestNeighbour));
selected = estimated>0.1;
idx1 = selected & ~bad; % final ripples
idx2 = ~selected & ~bad; % non-ripples (yet free from noise as well)
scatter(matrix(~bad,1),matrix(~bad,2),1,estimated(~bad),'filled'); colormap(Bright); set(gca,'CLim',[0 1])
xlim(xlims); ylim(ylims);
hold on;
scatter(swDiffAll(scored),ripPowerAll(scored),20,scores(scored)); scatter(swDiffAll(scored),ripPowerAll(scored),15,scores(scored));
set(get(colorbar,'YLabel'),'String','Estimated ripple score');
xlabel('Sharp wave depth'); ylabel('Ripple power');
title({['Manual scoring for ' basepath],['called with channels: ' num2str(Channels) ' (ripple, sharp wave, and (optionally) noise)']});
if rem(i,10)==0
display('saving...');
save(fullfile(basepath,'DetectSWR_manual_scoring.mat'),'t','swDiffAll','ripPowerAll','bad','scores');
end
end
1
0
1
0
1
0
1
0
[clean,~,badIntervals] = CleanLFP([tl,double(lfp(:,2))],'thresholds',[6 Inf],'manual',true);
bad = bad | InIntervals(t,badIntervals);
clf; plot(t,d); hold all; plot(xlim,ones(1,2)*threshold2,'r--'); ylabel('lfp derivative (z-units'); legend('derivative','threshold2');
disp('Change "threshold2" manually.')
threshold2 = 1;
clf; plot(t,d); hold all; plot(xlim,ones(1,2)*threshold2,'r--'); ylabel('lfp derivative (z-units'); legend('derivative','threshold2');
disp('Change "threshold2" manually.')
clf; plot(t,d); hold all; plot([xlim;xlim],[ones(1,2)*threshold2; -ones(1,2)*threshold2],'r--'); ylabel('lfp derivative (z-units'); legend('derivative','threshold2');
clf; plot(t,d); hold all; plot(xlim,ones(1,2)*threshold2,'r--');plot(xlim,-ones(1,2)*threshold2,'r--'); ylabel('lfp derivative (z-units'); legend('derivative','threshold2');
threshold2 = 2;
clf; plot(t,d); hold all; plot(xlim,ones(1,2)*threshold2,'r--');plot(xlim,-ones(1,2)*threshold2,'r--'); ylabel('lfp derivative (z-units'); legend('derivative','threshold2');
threshold2 = 0.5;
clf; plot(t,d); hold all; plot(xlim,ones(1,2)*threshold2,'r--');plot(xlim,-ones(1,2)*threshold2,'r--'); ylabel('lfp derivative (z-units'); legend('derivative','threshold2');
threshold2 = 0.8;
dbcont
Portion(bad)
figure(1); % plot a clean figure without these points and check indivudual ripples
clf
plot(swDiffAll(idx1 | idx2),ripPowerAll(idx1 | idx2),'k.','markersize',1); hold on
matrix = [swDiffAll ripPowerAll];
z = bsxfun(@rdivide,bsxfun(@minus,matrix,mean(matrix(~bad,:))),std(matrix(~bad,:)));
while true % click "Ctrl+C" to exit when you feel you're done
figure(1);
colors = Bright(1000);
[x,y,button] = ginput(1);
[~,i] = min(abs(x-swDiffAll.*(1-bad))+abs(y-ripPowerAll.*(1-bad)));
t = featureTs/1250;
tl = (1:length(lfp))'/1250;
figure(2);
clf; interval = t(i) + [-1 1]*0.5; in = tl>interval(1) & tl<interval(end);
subplot(3,1,1);
plot(tl(in) - t(i),lfp(in,:)); title(num2str(t(i)));
hold on; plot(bsxfun(@times,ones(1,2),(Restrict(t,interval) - t(i))),ylim,'k--');
legend('ripple channel','SW channel','noise channel');
subplot(3,1,2);
plot(tl(in) - t(i),double(swDiff(in)));
hold on; plot(bsxfun(@times,ones(1,2),(Restrict(t,interval) - t(i))),ylim,'k--');
subplot(3,1,3);
plot(tl(in) - t(i),double(ripPower0(in)));
hold on; plot(bsxfun(@times,ones(1,2),(Restrict(t,interval) - t(i))),ylim,'k--');
%     x = input( prompt )
commandwindow % select the command window
score = str2double(input('good(1)/bad(0)?','s'));
scores(i) = score; % save scores to optionally save
figure(1);
xlims = xlim; ylims = ylim; clf
% extrapolate the score of each ripple based on the score of the nearest scored ripple
scored = find(~isnan(scores));
distances = sqrt(bsxfun(@minus,z(:,1),z(scored,1)').^2 + bsxfun(@minus,z(:,2),z(~isnan(scores),2)').^2);
[~,nearestNeighbour] = min(distances,[],2);
estimated = scores(scored(nearestNeighbour));
selected = estimated>0.1;
idx1 = selected & ~bad; % final ripples
idx2 = ~selected & ~bad; % non-ripples (yet free from noise as well)
scatter(matrix(~bad,1),matrix(~bad,2),1,estimated(~bad),'filled'); colormap(Bright); set(gca,'CLim',[0 1])
xlim(xlims); ylim(ylims);
hold on;
scatter(swDiffAll(scored),ripPowerAll(scored),20,scores(scored)); scatter(swDiffAll(scored),ripPowerAll(scored),15,scores(scored));
set(get(colorbar,'YLabel'),'String','Estimated ripple score');
xlabel('Sharp wave depth'); ylabel('Ripple power');
title({['Manual scoring for ' basepath],['called with channels: ' num2str(Channels) ' (ripple, sharp wave, and (optionally) noise)']});
if rem(i,10)==0
display('saving...');
save(fullfile(basepath,'DetectSWR_manual_scoring.mat'),'t','swDiffAll','ripPowerAll','bad','scores');
end
end
0
10/ 0.21735
filtered = FilterLFP([tl lfp(:,2)],'passband',[40 50]);
filtered = FilterLFP([tl double(lfp(:,2))],'passband',[40 50]);
[~,power45] = Phase(filtered);
x = [-1 1]+22502.17;
PlotXY(Restrict(power45,x));
figure; PlotXY(Restrict(power45,x));
clf
PlotXY(power45);
badIntervals = tl(FindInterval(power45(:,2)>500));
bad = bad | InIntervals(t,badIntervals);
Portion(bad)
clf
PlotXY([tl double(lfp(:,2))]);
PlotIntervals(badIntervals);
PlotIntervals(badIntervals,'color','r');
x = xlim;
figure; PlotXY(Restrict(power45,x));
SideAxes(gca,'bottom',0.5)
PlotXY(filtered,'y');
filtered = FilterLFP([tl double(lfp(:,2))],'passband',[30 60]);
PlotXY(filtered,'k');
power45_old = power45;
[~,power45] = Phase(filtered);
PlotXY(Restrict(power45,x));
PlotXY(power45);
PlotHVLines(500,'h');
PlotHVLines(1500,'h');
PlotHVLines(1000,'h','k--','linewidth',2);
badIntervals = tl(FindInterval(power45(:,2)>1000));
bad = bad | InIntervals(t,badIntervals);
Portion(bad)
x = [-1 1]+22502.17;
xlim(x)
PlotIntervals(badIntervals,'color','y');
badIntervals = tl(FindInterval(power45(:,2)>1000)); badIntervals = bsxfun(@plus,badIntervals,[-1 1]*0.15);
badIntervals = ConsolidateIntervals(badIntervals);
PlotIntervals(badIntervals,'color','g');
bad = bad | InIntervals(t,badIntervals);
Portion(bad)
mmax(t)
mmax(badIntervals)
clf
plot(swDiffAll(idx2),ripPowerAll(idx2),'k.','markersize',1); hold on
plot(swDiffAll(idx1),ripPowerAll(idx1),'r.','markersize',1); legend('non-ripples','ripples');
figure(1); % plot a clean figure without these points and check indivudual ripples
clf
plot(swDiffAll(idx1 | idx2),ripPowerAll(idx1 | idx2),'k.','markersize',1); hold on
matrix = [swDiffAll ripPowerAll];
z = bsxfun(@rdivide,bsxfun(@minus,matrix,mean(matrix(~bad,:))),std(matrix(~bad,:)));
while true % click "Ctrl+C" to exit when you feel you're done
figure(1);
colors = Bright(1000);
[x,y,button] = ginput(1);
[~,i] = min(abs(x-swDiffAll.*(1-bad))+abs(y-ripPowerAll.*(1-bad)));
t = featureTs/1250;
tl = (1:length(lfp))'/1250;
figure(2);
clf; interval = t(i) + [-1 1]*0.5; in = tl>interval(1) & tl<interval(end);
subplot(3,1,1);
plot(tl(in) - t(i),lfp(in,:)); title(num2str(t(i)));
hold on; plot(bsxfun(@times,ones(1,2),(Restrict(t,interval) - t(i))),ylim,'k--');
legend('ripple channel','SW channel','noise channel');
subplot(3,1,2);
plot(tl(in) - t(i),double(swDiff(in)));
hold on; plot(bsxfun(@times,ones(1,2),(Restrict(t,interval) - t(i))),ylim,'k--');
subplot(3,1,3);
plot(tl(in) - t(i),double(ripPower0(in)));
hold on; plot(bsxfun(@times,ones(1,2),(Restrict(t,interval) - t(i))),ylim,'k--');
%     x = input( prompt )
commandwindow % select the command window
score = str2double(input('good(1)/bad(0)?','s'));
scores(i) = score; % save scores to optionally save
figure(1);
xlims = xlim; ylims = ylim; clf
% extrapolate the score of each ripple based on the score of the nearest scored ripple
scored = find(~isnan(scores));
distances = sqrt(bsxfun(@minus,z(:,1),z(scored,1)').^2 + bsxfun(@minus,z(:,2),z(~isnan(scores),2)').^2);
[~,nearestNeighbour] = min(distances,[],2);
estimated = scores(scored(nearestNeighbour));
selected = estimated>0.1;
idx1 = selected & ~bad; % final ripples
idx2 = ~selected & ~bad; % non-ripples (yet free from noise as well)
scatter(matrix(~bad,1),matrix(~bad,2),1,estimated(~bad),'filled'); colormap(Bright); set(gca,'CLim',[0 1])
xlim(xlims); ylim(ylims);
hold on;
scatter(swDiffAll(scored),ripPowerAll(scored),20,scores(scored)); scatter(swDiffAll(scored),ripPowerAll(scored),15,scores(scored));
set(get(colorbar,'YLabel'),'String','Estimated ripple score');
xlabel('Sharp wave depth'); ylabel('Ripple power');
title({['Manual scoring for ' basepath],['called with channels: ' num2str(Channels) ' (ripple, sharp wave, and (optionally) noise)']});
if rem(i,10)==0
display('saving...');
save(fullfile(basepath,'DetectSWR_manual_scoring.mat'),'t','swDiffAll','ripPowerAll','bad','scores');
end
end
0
1
scores(i) = 0;
while true % click "Ctrl+C" to exit when you feel you're done
figure(1);
colors = Bright(1000);
[x,y,button] = ginput(1);
[~,i] = min(abs(x-swDiffAll.*(1-bad))+abs(y-ripPowerAll.*(1-bad)));
t = featureTs/1250;
tl = (1:length(lfp))'/1250;
figure(2);
clf; interval = t(i) + [-1 1]*0.5; in = tl>interval(1) & tl<interval(end);
subplot(3,1,1);
plot(tl(in) - t(i),lfp(in,:)); title(num2str(t(i)));
hold on; plot(bsxfun(@times,ones(1,2),(Restrict(t,interval) - t(i))),ylim,'k--');
legend('ripple channel','SW channel','noise channel');
subplot(3,1,2);
plot(tl(in) - t(i),double(swDiff(in)));
hold on; plot(bsxfun(@times,ones(1,2),(Restrict(t,interval) - t(i))),ylim,'k--');
subplot(3,1,3);
plot(tl(in) - t(i),double(ripPower0(in)));
hold on; plot(bsxfun(@times,ones(1,2),(Restrict(t,interval) - t(i))),ylim,'k--');
%     x = input( prompt )
commandwindow % select the command window
score = str2double(input('good(1)/bad(0)?','s'));
scores(i) = score; % save scores to optionally save
figure(1);
xlims = xlim; ylims = ylim; clf
% extrapolate the score of each ripple based on the score of the nearest scored ripple
scored = find(~isnan(scores));
distances = sqrt(bsxfun(@minus,z(:,1),z(scored,1)').^2 + bsxfun(@minus,z(:,2),z(~isnan(scores),2)').^2);
[~,nearestNeighbour] = min(distances,[],2);
estimated = scores(scored(nearestNeighbour));
selected = estimated>0.1;
idx1 = selected & ~bad; % final ripples
idx2 = ~selected & ~bad; % non-ripples (yet free from noise as well)
scatter(matrix(~bad,1),matrix(~bad,2),1,estimated(~bad),'filled'); colormap(Bright); set(gca,'CLim',[0 1])
xlim(xlims); ylim(ylims);
hold on;
scatter(swDiffAll(scored),ripPowerAll(scored),20,scores(scored)); scatter(swDiffAll(scored),ripPowerAll(scored),15,scores(scored));
set(get(colorbar,'YLabel'),'String','Estimated ripple score');
xlabel('Sharp wave depth'); ylabel('Ripple power');
title({['Manual scoring for ' basepath],['called with channels: ' num2str(Channels) ' (ripple, sharp wave, and (optionally) noise)']});
if rem(i,10)==0
display('saving...');
save(fullfile(basepath,'DetectSWR_manual_scoring.mat'),'t','swDiffAll','ripPowerAll','bad','scores');
end
end
0
1
0
1
0
1
0
1
0
1
0
1
0
1
0
1
0
1
0
10
scores(i) = 0;
while true % click "Ctrl+C" to exit when you feel you're done
figure(1);
colors = Bright(1000);
[x,y,button] = ginput(1);
[~,i] = min(abs(x-swDiffAll.*(1-bad))+abs(y-ripPowerAll.*(1-bad)));
t = featureTs/1250;
tl = (1:length(lfp))'/1250;
figure(2);
clf; interval = t(i) + [-1 1]*0.5; in = tl>interval(1) & tl<interval(end);
subplot(3,1,1);
plot(tl(in) - t(i),lfp(in,:)); title(num2str(t(i)));
hold on; plot(bsxfun(@times,ones(1,2),(Restrict(t,interval) - t(i))),ylim,'k--');
legend('ripple channel','SW channel','noise channel');
subplot(3,1,2);
plot(tl(in) - t(i),double(swDiff(in)));
hold on; plot(bsxfun(@times,ones(1,2),(Restrict(t,interval) - t(i))),ylim,'k--');
subplot(3,1,3);
plot(tl(in) - t(i),double(ripPower0(in)));
hold on; plot(bsxfun(@times,ones(1,2),(Restrict(t,interval) - t(i))),ylim,'k--');
%     x = input( prompt )
commandwindow % select the command window
score = str2double(input('good(1)/bad(0)?','s'));
scores(i) = score; % save scores to optionally save
figure(1);
xlims = xlim; ylims = ylim; clf
% extrapolate the score of each ripple based on the score of the nearest scored ripple
scored = find(~isnan(scores));
distances = sqrt(bsxfun(@minus,z(:,1),z(scored,1)').^2 + bsxfun(@minus,z(:,2),z(~isnan(scores),2)').^2);
[~,nearestNeighbour] = min(distances,[],2);
estimated = scores(scored(nearestNeighbour));
selected = estimated>0.1;
idx1 = selected & ~bad; % final ripples
idx2 = ~selected & ~bad; % non-ripples (yet free from noise as well)
scatter(matrix(~bad,1),matrix(~bad,2),1,estimated(~bad),'filled'); colormap(Bright); set(gca,'CLim',[0 1])
xlim(xlims); ylim(ylims);
hold on;
scatter(swDiffAll(scored),ripPowerAll(scored),20,scores(scored)); scatter(swDiffAll(scored),ripPowerAll(scored),15,scores(scored));
set(get(colorbar,'YLabel'),'String','Estimated ripple score');
xlabel('Sharp wave depth'); ylabel('Ripple power');
title({['Manual scoring for ' basepath],['called with channels: ' num2str(Channels) ' (ripple, sharp wave, and (optionally) noise)']});
if rem(i,10)==0
display('saving...');
save(fullfile(basepath,'DetectSWR_manual_scoring.mat'),'t','swDiffAll','ripPowerAll','bad','scores');
end
end
0
1
0
% Save your progress!
display('saving...');
save(fullfile(basepath,'DetectSWR_manual_scoring.mat'),'t','swDiffAll','ripPowerAll','idx1','idx2','bad','selected','scores');
saveas(figure(1),fullfile(basepath,'DetectSWR_manual_scoring.fig'));
clf
plot(swDiffAll(idx2),ripPowerAll(idx2),'k.','markersize',1); hold on
plot(swDiffAll(idx1),ripPowerAll(idx1),'r.','markersize',1); legend('non-ripples','ripples');
mean(swDiffAll(idx1))
sum(idx2)
idx2(swDiffAll>mean(swDiffAll(idx1))) = 0;
idx2(ripPowerAll>mean(ripPowerAll(idx1))) = 0;
sum(idx2)
dbcont
PETH(hpc(:,1),ripples.timestamps(:,1));
PETH(pfc(:,1),ripples.timestamps(:,1));
PlotHVLines(0,'v','k--','linewidth',2);
figure; PETH(hpc(:,1),ripples.timestamps(:,1));
PlotHVLines(0,'v','k--','linewidth',2);
channel = 1+85;
% Detect delta waves
if exist(fullfile(basepath,[basename '.deltaWaves.events.mat']),'file')
load(fullfile(basepath,[basename '.deltaWaves.events.mat']),'deltaWaves');
try deltas = [deltaWaves.timestamps(:,1) deltaWaves.peaks deltaWaves.timestamps(:,2) deltaWaves.peakNormedPower];
catch
deltas = repmat(deltaWaves.peaks,1,2);
end
else
lfp = GetAyaLFP(channel);
display(['loaded! @' num2str(toc)]);
[clean,~,badIntervals] = CleanLFP(lfp,'thresholds',[6 10],'manual',true);
deltas0 = FindDeltaPeaks(clean);
deltas = deltas0(deltas0(:,5)-deltas0(:,6)>3,:); % these thresholds should be manually refined for each session
deltaWaves.timestamps = deltas(:,[1 3]); deltaWaves.peaks = deltas(:,2);  deltaWaves.peakNormedPower = deltas(:,5); deltaWaves.detectorName = ['channel ' num2str(channel) '(+1), CleanLFP, FindDeltaPeaks, peak-trough>3'];
save(fullfile(basepath,[basename '.deltaWaves.events.mat']),'deltaWaves');
end
PlotHVLines([-1 1]*2.5,'v','k--','linewidth',2);
PlotHVLines([-1 1]*2.5,'h','k--','linewidth',2);
threshold1 = 3;
clf; plot(t,d); hold all; plot(xlim,ones(1,2)*threshold2,'r--');plot(xlim,-ones(1,2)*threshold2,'r--'); ylabel('lfp derivative (z-units'); legend('derivative','threshold2');
disp('Change "threshold2" manually.')
threshold2 = 1;
dbcont
clear all
basepath = 'N:\OJRproject\OJR49\day6';
channel = 1+85;
rippleChannels = 1+[18 30]; % the channel that is lower during the sharp-wave is second
[parentFolder,basename] = fileparts(basepath);
[~,projectName] = fileparts(parentFolder);
sessionID = [projectName '_' basename];
cd(basepath);
load(fullfile(basepath,[basename '.session.mat']),'session');
try
%     channel_mapping
load(fullfile(basepath,[basename '.cell_metrics.cellinfo.mat']),'cell_metrics');
nNeurons = length(cell_metrics.spikes.times);
regions = cell_metrics.brainRegion;
regionNames = unique(regions);
regionCell = zeros(length(regions),1);
for i=1:length(regionNames)
regionCell(strcmp(regions,regionNames{i})) = i;
end
spikesCell = cell_metrics.spikes.times';
PFCindex = [find(strcmp(regionNames,'ILA')) find(strcmp(regionNames,'PFC'))];
HPCindex = find(strcmp(regionNames,'CA1'));
pfc = []; if ~isempty(PFCindex), pfc = sortrows(Group(spikesCell{regionCell==PFCindex})); end
hpc = []; if ~isempty(HPCindex), hpc = sortrows(Group(spikesCell{regionCell==HPCindex})); end
neurons = true;
catch
neurons = false; pfc = zeros(0,2); hpc = zeros(0,2); spikesCell = {};
end
spikes = sortrows(Group(spikesCell{:}));
lfp = GetAyaLFP(channel);
display(['loaded! @' num2str(toc)]);
[clean,~,badIntervals] = CleanLFP(lfp,'thresholds',[3 1],'manual',false);
PlotXY(clean);
clf
PlotXY(clean);
deltas0 = FindDeltaPeaks(clean);
PETH(pfc(:,1),deltas0(:,2));
[h,ht] = PETH(pfc(:,1),deltas0(:,2));
clf
PlotColorMap(h,'x',ht);
PlotColorMap(sortby(h,deltas0(:,5)-deltas0(:,6)),'x',ht);
PlotColorMap(Shrink(sortby(h,deltas0(:,5)-deltas0(:,6)),100,1),'x',ht);
score = deltas0(:,5)-deltas0(:,6);
rank = tiedrank(score)/length(score);
open rank
PlotColorMap(Shrink(sortby(h,deltas0(:,5)-deltas0(:,6)),1,1),'x',ht);
PlotColorMap(Shrink(sortby(h,-(deltas0(:,5)-deltas0(:,6))),1,1),'x',ht);
rank = tiedrank(-score);
questionable = deltas0(rank<300,2);
mean(score(rank<300))
mean(score(rank>300))
clf
PETH(pfc(:,1),deltas0(rank>300,2));
hold all
PETH(pfc(:,1),deltas0(rank<300,2));
PlotXY(clean);
PlotHVLines(questionable,'v','r');
clf
PETH(deltas0(:,2),badIntervals,'duration',[-1 1]*30);
PETH(deltas0(:,2),badIntervals(:,1),'duration',[-1 1]*30);
hold all
PETH(questionable,badIntervals(:,1),'duration',[-1 1]*30);
PETH(deltas0(rank>300,2),badIntervals(:,1),'duration',[-1 1]*30);
badIntervals = bsxfun(@plus,badIntervals,[-1 1]*5);
badIntervals = ConsolidateIntervals(badIntervals);
clf
[clean,~,badIntervals] = CleanLFP(lfp,'thresholds',[3 1],'manual',false,'aroundArtefact',[5 1]);
clf
PlotXY(clean);
deltas00  = deltas0;
deltas0 = FindDeltaPeaks(clean);
[h,ht] = PETH(pfc(:,1),deltas0(:,2));
PlotColorMap(Shrink(sortby(h,deltas0(:,5)-deltas0(:,6)),100,1),'x',ht);
deltas0(end,1)
deltas0(end,1)/3600
SleepScoreMaster(basepath,'noPrompts',true,'rejectchannels',1+[22 23 31 71 96 109 114 113 117 118 104 119 103 120]); % takes lfp in base 0
%-- 4/12/2022 5:39 PM --%
cd C:
basepath = 'N:\OJRproject\OJR49\day6';
SleepScoreMaster(basepath,'noPrompts',true,'rejectchannels',1+[22 23 31 71 96 109 114 113 117 118 104 119 103 120]); % takes lfp in base 0
t_clus(end)
display('Calculating Theta Metric above PSS')
%Put the LFP in the right structure format
lfp.data = thLFP;
lfp.timestamps = t_LFP;
lfp.samplingRate = sf_LFP;
%Calculate PSS
[specslope,spec] = bz_PowerSpectrumSlope(lfp,window,window-noverlap,...
'nfreqs',200,'frange',f_all,'IRASA',ThIRASA);
t_thclu = specslope.timestamps;
specdt = 1./specslope.samplingRate;
thFFTspec = specslope.resid';
thFFTspec(thFFTspec<0)=0;
IRASAsmooth_th = spec.IRASAsmooth';
thFFTspec_raw = 10.^spec.amp';
% Remove transients before calculating SW histogram
zFFTspec = NormToInt(spec.amp,'modZ');
totz = NormToInt(abs(sum(zFFTspec,2)),'modZ');
badtimes_TH = find(totz>3);
thFFTfreqs = specslope.freqs';
thfreqs = (thFFTfreqs>=f_theta(1) & thFFTfreqs<=f_theta(2));
thratio = max((thFFTspec(thfreqs,:)),[],1);
thratio = thratio';
ignoretimeIDX = InIntervals(t_clus,ignoretime) | isnan(broadbandSlowWave) | isnan(thratio);
dbcont
cd(basepath);
load('day6.SleepScoreLFP.LFP.mat')
load('day6.SleepScoreLFP.LFP.mat', 'SleepScoreLFP')
load('day6.SleepStateEpisodes.states.mat')
channel = 1+85;
lfp = GetAyaLFP(channel);
display(['loaded! @' num2str(toc)]);
[clean,~,badIntervals] = CleanLFP(lfp,'thresholds',[3 1],'manual',false,'aroundArtefact',[5 1]);
sws = SleepStateEpisodes.ints.NREMepisode;
SleepStateEpisodes = getStruct(basepath,'SleepStateEpisodes');
sws = SleepStateEpisodes.ints.NREMepisode;
clean = Restrict(clean,sws);
PlotXY(clean);
close all
PlotXY(clean);
deltas0 = FindDeltaPeaks(clean);
try
%     channel_mapping
load(fullfile(basepath,[basename '.cell_metrics.cellinfo.mat']),'cell_metrics');
nNeurons = length(cell_metrics.spikes.times);
regions = cell_metrics.brainRegion;
regionNames = unique(regions);
regionCell = zeros(length(regions),1);
for i=1:length(regionNames)
regionCell(strcmp(regions,regionNames{i})) = i;
end
spikesCell = cell_metrics.spikes.times';
PFCindex = [find(strcmp(regionNames,'ILA')) find(strcmp(regionNames,'PFC'))];
HPCindex = find(strcmp(regionNames,'CA1'));
pfc = []; if ~isempty(PFCindex), pfc = sortrows(Group(spikesCell{regionCell==PFCindex})); end
hpc = []; if ~isempty(HPCindex), hpc = sortrows(Group(spikesCell{regionCell==HPCindex})); end
neurons = true;
catch
neurons = false; pfc = zeros(0,2); hpc = zeros(0,2); spikesCell = {};
end
[h,ht] = PETH(pfc(:,1),deltas0(:,2));
PlotColorMap(Shrink(sortby(h,deltas0(:,5)-deltas0(:,6)),100,1),'x',ht);
try
%     channel_mapping
load(fullfile(basepath,[basename '.cell_metrics.cellinfo.mat']),'cell_metrics');
nNeurons = length(cell_metrics.spikes.times);
regions = cell_metrics.brainRegion;
regionNames = unique(regions);
regionCell = zeros(length(regions),1);
for i=1:length(regionNames)
regionCell(strcmp(regions,regionNames{i})) = i;
end
spikesCell = cell_metrics.spikes.times';
PFCindex = [find(strcmp(regionNames,'ILA')) find(strcmp(regionNames,'PFC'))];
HPCindex = find(strcmp(regionNames,'CA1'));
pfc = []; if ~isempty(PFCindex), pfc = sortrows(Group(spikesCell{regionCell==PFCindex})); end
hpc = []; if ~isempty(HPCindex), hpc = sortrows(Group(spikesCell{regionCell==HPCindex})); end
neurons = true;
catch
neurons = false; pfc = zeros(0,2); hpc = zeros(0,2); spikesCell = {};
end
load(fullfile(basepath,[basename '.cell_metrics.cellinfo.mat']),'cell_metrics');
cell_metrics = getStruct(basepath,'cell_metrics']);
cell_metrics = getStruct(basepath,'cell_metrics');
try
%     channel_mapping
cell_metrics = getStruct(basepath,'cell_metrics');
nNeurons = length(cell_metrics.spikes.times);
regions = cell_metrics.brainRegion;
regionNames = unique(regions);
regionCell = zeros(length(regions),1);
for i=1:length(regionNames)
regionCell(strcmp(regions,regionNames{i})) = i;
end
spikesCell = cell_metrics.spikes.times';
PFCindex = [find(strcmp(regionNames,'ILA')) find(strcmp(regionNames,'PFC'))];
HPCindex = find(strcmp(regionNames,'CA1'));
pfc = []; if ~isempty(PFCindex), pfc = sortrows(Group(spikesCell{regionCell==PFCindex})); end
hpc = []; if ~isempty(HPCindex), hpc = sortrows(Group(spikesCell{regionCell==HPCindex})); end
neurons = true;
catch
neurons = false; pfc = zeros(0,2); hpc = zeros(0,2); spikesCell = {};
end
[h,ht] = PETH(pfc(:,1),deltas0(:,2));
PlotColorMap(Shrink(sortby(h,deltas0(:,5)-deltas0(:,6)),100,1),'x',ht);
PlotColorMap(Shrink(sortby(h,deltas0(:,5)-deltas0(:,6)),10,1),'x',ht);
PlotColorMap(Shrink(sortby(h,-(deltas0(:,5)-deltas0(:,6))),10,1),'x',ht);
questionable = deltas0(rank<100,2);
score = deltas0(:,5)-deltas0(:,6);
rank = tiedrank(-score);
PETH(pfc(:,1),questionable);
PETH(pfc(:,1),deltas0(rank<100,2));
hold all; PETH(pfc(:,1),deltas0(rank>100,2));
SaveCustomEvents('deltas.del.evt',deltas(rank<100,1:3),{'deltas start','delta peak','deltas stop'});a
SaveCustomEvents('deltas.del.evt',deltas0(rank<100,1:3),{'deltas start','delta peak','deltas stop'});a
SaveCustomEvents('deltas.del.evt',deltas0(rank<100,1:3),{'deltas start','delta peak','deltas stop'});
eegFile = fullfile(basepath,[basename '.eeg']); % This is where I will store the normalized LFP file
[parentFolder,basename] = fileparts(basepath);
[~,projectName] = fileparts(parentFolder);
sessionID = [projectName '_' basename];
cd(basepath);
eegFile = fullfile(basepath,[basename '.eeg']); % This is where I will store the normalized LFP file
lfpFile = fullfile(basepath,[basename '.lfp']);
copyfile(lfpFile,eegFile);
file = memmapfile(eegFile,'Format','int16','Writable',true);
data = reshape(file.Data,128,[]);
badChannels = 1+[22 23 31 71 96 109 114 113 117 118 104 119 103 120];
okChannels = true(128,1); okChannels(badChannels) = false;
okChannels = true(128,1); okChannels(badChannels) = false;
m = int16(mean(data(okChannels,:)));
newData = bsxfun(@minus,data,m);
file.data = newData(:);
clear file
eeg = GetAyaEEG(channel);
cleanLFP =clean;
clean = CleanLFP(eeg,'thresholds',[3 1],'manual',false,'aroundArtefact',[5 1],'manual',true);
threshold1
threshold1 = 5;
plot(t,z); hold all; plot(xlim,ones(1,2)*threshold1,'r--'); ylabel('lfp signal (z-units'); legend('signal','threshold1');
clf; plot(t,d); hold all; plot(xlim,ones(1,2)*threshold2,'r--');plot(xlim,-ones(1,2)*threshold2,'r--'); ylabel('lfp derivative (z-units'); legend('derivative','threshold2');
threshold2
threshold2 = 2;
clf; plot(t,d); hold all; plot(xlim,ones(1,2)*threshold2,'r--');plot(xlim,-ones(1,2)*threshold2,'r--'); ylabel('lfp derivative (z-units'); legend('derivative','threshold2');
dbcont
deltas0_1 = deltas0;
deltas0 = FindDeltaPeaks(clean);
clf
[h,ht] = PETH(pfc(:,1),deltas0(:,2));
PlotColorMap(Shrink(sortby(h,-(deltas0(:,5)-deltas0(:,6))),10,1),'x',ht);
deltas0 = deltas0_1;
[h,ht] = PETH(pfc(:,1),deltas0(:,2));
PlotColorMap(Shrink(sortby(h,-(deltas0(:,5)-deltas0(:,6))),10,1),'x',ht);
hist(deltas0(rank<100,2),1000);
clf
hist(deltas0(rank<100,2),1000);
PlotXY(clean);
xlim([-5 5]+405);
xlim([-1 1]*2+406.9);
channel
rippleLFP = GetAyaLFP(30);
clf
hist(deltas0(rank<100,2),1000);
badIntervals = ConsolidateIntervals([badIntervals; 390 450];
badIntervals = ConsolidateIntervals([badIntervals; 390 450]);
clean = Restrict(lfp,SubtractIntervals([0 Inf],badIntervals));
clf
PlotXY(clean);
hist(deltas0(rank<100,2),1000);
hist(deltas0(:,2),1000);
hist(deltas0(:,2),10000);
hist(deltas0(:,2),5000);
deltas0 = FindDeltaPeaks(clean);
[h,ht] = PETH(pfc(:,1),deltas0(:,2));
PlotColorMap(Shrink(sortby(h,-(deltas0(:,5)-deltas0(:,6))),10,1),'x',ht);
[h,ht] = PETH(pfc(:,1),deltas0(:,1));
PlotColorMap(Shrink(sortby(h,-(deltas0(:,5)-deltas0(:,6))),10,1),'x',ht);
[h,ht] = PETH(pfc(:,1),deltas0(:,1));
PlotColorMap(Shrink(sortby(h,-(deltas0(:,5)-deltas0(:,6))),100,1),'x',ht);
PlotColorMap(Shrink(sortby(h,-(deltas0(:,5)-deltas0(:,6))),70,1),'x',ht);
PlotColorMap(Shrink(sortby(h,-(deltas0(:,5)-deltas0(:,6))),72,1),'x',ht);
rank = tiedrank(-score)/length(score);
semplot(ht,h(rank>70));
semplot(ht,h(rank>70,:));
score = deltas0(:,5)-deltas0(:,6);
rank = tiedrank(-score)/length(score);
semplot(ht,h(rank>70,:));
[h,ht] = PETH(pfc(:,1),deltas0(:,1));
clf
semplot(ht,h(rank>70,:));
semplot(ht,h(rank>0.70,:));
hold all
figure; PlotColorMap(Shrink(sortby(h,-(deltas0(:,5)-deltas0(:,6))),72,1),'x',ht);
semplot(ht,h(rank<0.05,:),'r');
semplot(ht,h(rank<0.03,:),'b');
semplot(ht,h(rank>0.03 & rank<0.05,:),'k');
bad = rank<0.03 | rank>0.7;
clf
semplot(ht,h(~bad,:),'k');
deltas = deltas0(~bad,:); % these thresholds should be manually refined for each session
[h,ht] = PETH(ripples(:,2),deltas0(:,1));
r = ripples.timestamps;
% Detect ripples
if exist(fullfile(basepath,[basename '.ripples.events.mat']),'file')
load(fullfile(basepath,[basename '.ripples.events.mat']),'ripples');
else
load(fullfile(basepath,[basename '.session.mat']),'session');
%     rippleChannels = swrChannels('basepath',basepath); % or better yet, pick them on neuroscope (don't forget to add 1! rippleChannels.Ripple_Channel = neuroscope_Channel + 1)
ripples = DetectSWR(rippleChannels,'basepath',basepath,'saveMat',true,'forceDetect',true,'check',true,'useSPW',false);
end
r = ripples.timestamps;
clf
PETH(r(:,1),deltas0(:,2))
PETH(r(:,1),deltas0(~bad,2))
hold all
PETH(r(:,1),deltas0(bad,2))
deltaWaves.timestamps = deltas(:,[1 3]); deltaWaves.peaks = deltas(:,2);  deltaWaves.peakNormedPower = deltas(:,5); deltaWaves.detectorName = ['channel ' num2str(channel) '(+1), CleanLFP, FindDeltaPeaks, peak-trough>3'];
save(fullfile(basepath,[basename '.deltaWaves.events.mat']),'deltaWaves');
[h,ht] = PETH(pfc(:,1),deltas(:,1));
PlotColorMap(Shrink(sortby(h,-(deltas(:,5)-deltas(:,6))),72,1),'x',ht);
[h,ht] = PETH(r(:,1),deltas(:,1));
PlotColorMap(Shrink(sortby(h,-(deltas(:,5)-deltas(:,6))),72,1),'x',ht);
PlotColorMap(Shrink(sortby(h,-(deltas(:,5)-deltas(:,6))),200,1),'x',ht);
PlotColorMap(Shrink(sortby(h,-(deltas(:,5)-deltas(:,6))),500,1),'x',ht);
clf
semplot(ht,h);
clf
[h,ht] = PETH(r(:,1),deltas(:,2));
semplot(ht,h);
clf
[h,ht] = PETH(r(:,1),deltas(:,2),'durations',[-1 1]*5);
semplot(ht,h);
[h,ht] = PETH(r(:,1),deltas(:,2),'durations',[-1 1]*2);
clf
semplot(ht,h);
load('day6.MergePoints.events.mat')
sleep = ConsolidateIntervals(MergePoints.timestamps(cellfun(@(x) any(strfind(lower(x),'sleep')),MergePoints.foldernames),:));
sleep
clf
pre = deltas(:,2)<sleep(1,2);
post = deltas(:,2)>sleep(2,1);
semplot(ht,h(pre,:),'k');
semplot(ht,h(post,:),'r');
ripples = r;
ripples = Restrict(ripples,sws);
upstate = [deltas(1:end-1,2) deltas(2:end,2)]; upstate(~InIntervals(diff(upstate,[],2),[1 4]),:) = [];
rt = RelativeTime(ripples(:,1),upstate);
[h,ht] = PETH(pfc(:,1),ripples(:,1),'durations',[-1 1]*0.5,'nBins',201);
[hd,ht] = PETH(deltas(:,2),ripples(:,1),'durations',[-1 1]*0.5,'nBins',201);
subplot(3,4,(condition-1)*4 + i)
rt = ripples(:,1) - deltas(FindClosest(deltas(:,2),ripples(:,1)),2);
ok = InIntervals(rt,[-1 1]*0.5);
[~,m] = max(Shrink(sortby(hd(ok,:),rt(ok)),floor(sum(ok)/100),1),[],2);
PlotColorMap(Smooth(Shrink(sortby(h(ok,:),rt(ok)),floor(sum(ok)/100),1),2),'x',ht);
hold all
plot(ht(m),1:max(ylim),'w','linewidth',2);
PlotHVLines(0,'v','k--','linewidth',2);
title(basepath);
drawnow
rt = RelativeTime(ripples(:,1),upstate);
[h,ht] = PETH(pfc(:,1),ripples(:,1),'durations',[-1 1]*0.5,'nBins',201);
[hd,ht] = PETH(deltas(:,2),ripples(:,1),'durations',[-1 1]*0.5,'nBins',201);
rt = ripples(:,1) - deltas(FindClosest(deltas(:,2),ripples(:,1)),2);
ok = InIntervals(rt,[-1 1]*0.5);
[~,m] = max(Shrink(sortby(hd(ok,:),rt(ok)),floor(sum(ok)/100),1),[],2);
PlotColorMap(Smooth(Shrink(sortby(h(ok,:),rt(ok)),floor(sum(ok)/100),1),2),'x',ht);
hold all
plot(ht(m),1:max(ylim),'w','linewidth',2);
PlotHVLines(0,'v','k--','linewidth',2);
title(basepath);
clear all
56
raly
45
cd C:
SaveFig(fullfile(['M:\home\raly\results\PFC\reactivation\PFC_responses_to_ripples_colormaps_ordered_by_proximity_to_delta']));
raly
cd C:\Users\Cornell\Documents\GitHub\neurocode\projects\OJR
X = Xl;
i = 5;
[pfc,hpc,deltas,ripples,sleep,session,MergePoints,clickers,sws,basepath] = X{i,:};
[parentFolder,basename] = fileparts(basepath);
[~,projectName] = fileparts(parentFolder);
sessionID = [projectName '_' basename];
sessionIDs{i,1} = sessionID;
preSleep = SubtractIntervals(sleep(1,:), SubtractIntervals([0 Inf],sws));
postSleep = SubtractIntervals(sleep(2:end,:), SubtractIntervals([0 Inf],sws));
ripples = Restrict(ripples,sws);
deltas = Restrict(deltas,sws);
pre = InIntervals(deltas,preSleep);
post = InIntervals(deltas,postSleep);
% Restrict to 1h of sws only:
limit = Unshift(3600,preSleep); preSleep(preSleep(:,1)>limit,:) = []; if preSleep(end,2)>limit, preSleep(end,2) = limit; end
limit = Unshift(3600,postSleep); postSleep(postSleep(:,1)>limit,:) = []; if postSleep(end,2)>limit, postSleep(end,2) = limit; end
training = ConsolidateIntervals(MergePoints.timestamps(cellfun(@(x) ~isempty(strfind(x,'train')),MergePoints.foldernames),:)); % training epochs are named with "train" somewhere in the folder name
exploration = sortrows([clickers{2};clickers{3}]); exploration(any(isnan(exploration),2),:) = []; exploration = ConsolidateIntervals(exploration);
exploration = Restrict(exploration,training);
spikes = pfc;
% remove lost units
[h,ht] = Dist(0:600:spikes(end,1),spikes,'grouped'); % bin spikes in 10-minute bins
badUnits = sum(h==0,1)'>1; % any empty 10-minute bins are signs of cells that should be excluded. Allow a single exception bin
newIDs = cumsum(~badUnits); newIDs(badUnits) = 0;
spikes(:,2) = newIDs(spikes(:,2));
spikes(spikes(:,2)==0,:) = []; % remove bad units
intervals = training;
rng(0);
[templates,correlations,weights] = ActivityTemplatesICA(Restrict(spikes,intervals,'shift','on'),'binsize',0.1);
re = ReactivationStrength(spikes,templates,'step',0.01,'binSize',0.1);
[templates,correlations,weights] = ActivityTemplates(Restrict(spikes,intervals,'shift','on'),'binsize',0.1);
re = ReactivationStrength(spikes,templates,'step',0.01,'binSize',0.1);
for j=1:size(re,2)-1
[q,qt] = PETH(re(:,[1 1+j]),deltas(:,2),'durations',[-1 1]*2,'nBins',501);
kk = kk+1;
qq{kk,1} = nanmean(q(pre,:),1);
qq{kk,2} = nanmean(q(post,:),1);
qq{kk,3} = i;
%                             qqs{kk,1} = q;
end
kk=0;
rng(0);
[templates,correlations,weights] = ActivityTemplates(Restrict(spikes,intervals,'shift','on'),'binsize',0.1);
re = ReactivationStrength(spikes,templates,'step',0.01,'binSize',0.1);
for j=1:size(re,2)-1
[q,qt] = PETH(re(:,[1 1+j]),deltas(:,2),'durations',[-1 1]*2,'nBins',501);
kk = kk+1;
qq{kk,1} = nanmean(q(pre,:),1);
qq{kk,2} = nanmean(q(post,:),1);
qq{kk,3} = i;
%                             qqs{kk,1} = q;
end
q1 = cell2mat(qq(:,1));
q2 = cell2mat(qq(:,2));
PlotColorMap(q1);
PlotColorMap([q1 q2]);
%
nComponents = Accumulate(cell2mat(qq(:,end)));
addSlash = @(x) strrep(x,'_','\_');
sessions = cellfun(addSlash,sessionIDs(nComponents>0),'UniformOutput',0);
lines = cumsum(nComponents(nComponents>0))+0.5;
q = (cell2mat(qq(:,1:2)));
%         z = zBaseline(q,Unfind((length(qt)+1):(length(qt)*2),(length(qt)*2)),2);
z = nanzscore(q,[],2);
z1 = z(:,1:length(qt)); z2 = z(:,(length(qt)+1):end);
smooth = [0 2];
z1 = zscore(z1,[],2); z2 = zscore(z2,[],2);
matrices = {z1,z2,z2-z1};
j=2;
PlotColorMap(Smooth(matrices{j},smooth),~isnan(matrices{j}),'x',qt); PlotHVLines(0,'v','k--','linewidth',2);
set(gca,'ytick',mean([[0;lines(1:end-1)] lines],2),'yticklabel',sessions);
PlotHVLines(lines,'h','w--','linewidth',2);
title([name ', post'])
PlotHVLines(0,'v','k--','linewidth',2);
for i=1:3, subplot(2,3,i); PlotColorMap(matrices{i},'x',qt);
end
for i=1:3, subplot(2,3,i+3); semplot(qt,matrices{i}); end
ylim([-2 4]);
semplot(ht,matrices{2},'r');
semplot(qt,matrices{2},'r');
PlotHVLines(0,'h','k--','linewidth',2);
PlotHVLines(0,'v','k--','linewidth',2);
X = Xs;
open Xs
i=4;
figure; [pfc,hpc,deltas,ripples,sleep,session,MergePoints,clickers,sws,basepath] = X{i,:};
[parentFolder,basename] = fileparts(basepath);
[~,projectName] = fileparts(parentFolder);
sessionID = [projectName '_' basename];
sessionIDs{i,1} = sessionID;
preSleep = SubtractIntervals(sleep(1,:), SubtractIntervals([0 Inf],sws));
postSleep = SubtractIntervals(sleep(2:end,:), SubtractIntervals([0 Inf],sws));
ripples = Restrict(ripples,sws);
deltas = Restrict(deltas,sws);
pre = InIntervals(deltas,preSleep);
post = InIntervals(deltas,postSleep);
% Restrict to 1h of sws only:
limit = Unshift(3600,preSleep); preSleep(preSleep(:,1)>limit,:) = []; if preSleep(end,2)>limit, preSleep(end,2) = limit; end
limit = Unshift(3600,postSleep); postSleep(postSleep(:,1)>limit,:) = []; if postSleep(end,2)>limit, postSleep(end,2) = limit; end
training = ConsolidateIntervals(MergePoints.timestamps(cellfun(@(x) ~isempty(strfind(x,'train')),MergePoints.foldernames),:)); % training epochs are named with "train" somewhere in the folder name
exploration = sortrows([clickers{2};clickers{3}]); exploration(any(isnan(exploration),2),:) = []; exploration = ConsolidateIntervals(exploration);
exploration = Restrict(exploration,training);
spikes = pfc;
% remove lost units
[h,ht] = Dist(0:600:spikes(end,1),spikes,'grouped'); % bin spikes in 10-minute bins
badUnits = sum(h==0,1)'>1; % any empty 10-minute bins are signs of cells that should be excluded. Allow a single exception bin
newIDs = cumsum(~badUnits); newIDs(badUnits) = 0;
spikes(:,2) = newIDs(spikes(:,2));
spikes(spikes(:,2)==0,:) = []; % remove bad units
intervals = training;
rng(0);
[templates,correlations,weights] = ActivityTemplatesICA(Restrict(spikes,intervals,'shift','on'),'binsize',0.1);
re = ReactivationStrength(spikes,templates,'step',0.01,'binSize',0.1);
[templates,correlations,weights] = ActivityTemplates(Restrict(spikes,intervals,'shift','on'),'binsize',0.1);
re = ReactivationStrength(spikes,templates,'step',0.01,'binSize',0.1);
for j=1:size(re,2)-1
[q,qt] = PETH(re(:,[1 1+j]),deltas(:,2),'durations',[-1 1]*2,'nBins',501);
kk = kk+1;
qq{kk,1} = nanmean(q(pre,:),1);
qq{kk,2} = nanmean(q(post,:),1);
qq{kk,3} = i;
%                             qqs{kk,1} = q;
end
kk=0;
rng(0);
[templates,correlations,weights] = ActivityTemplates(Restrict(spikes,intervals,'shift','on'),'binsize',0.1);
re = ReactivationStrength(spikes,templates,'step',0.01,'binSize',0.1);
for j=1:size(re,2)-1
[q,qt] = PETH(re(:,[1 1+j]),deltas(:,2),'durations',[-1 1]*2,'nBins',501);
kk = kk+1;
qq{kk,1} = nanmean(q(pre,:),1);
qq{kk,2} = nanmean(q(post,:),1);
qq{kk,3} = i;
%                             qqs{kk,1} = q;
end
q1 = cell2mat(qq(:,1));
q2 = cell2mat(qq(:,2));
PlotColorMap(q1);
PlotColorMap([q1 q2]);
%
nComponents = Accumulate(cell2mat(qq(:,end)));
addSlash = @(x) strrep(x,'_','\_');
sessions = cellfun(addSlash,sessionIDs(nComponents>0),'UniformOutput',0);
lines = cumsum(nComponents(nComponents>0))+0.5;
q = (cell2mat(qq(:,1:2)));
%         z = zBaseline(q,Unfind((length(qt)+1):(length(qt)*2),(length(qt)*2)),2);
z = nanzscore(q,[],2);
z1 = z(:,1:length(qt)); z2 = z(:,(length(qt)+1):end);
smooth = [0 2];
z1 = zscore(z1,[],2); z2 = zscore(z2,[],2);
matrices = {z1,z2,z2-z1};
j=2;
PlotColorMap(Smooth(matrices{j},smooth),~isnan(matrices{j}),'x',qt); PlotHVLines(0,'v','k--','linewidth',2);
set(gca,'ytick',mean([[0;lines(1:end-1)] lines],2),'yticklabel',sessions);
PlotHVLines(lines,'h','w--','linewidth',2);
title([name ', post'])
PlotHVLines(0,'v','k--','linewidth',2);
for i=1:3, subplot(2,3,i); PlotColorMap(matrices{i},'x',qt);
end
for i=1:3, subplot(2,3,i+3); semplot(qt,matrices{i}); end
ylim([-2 4]);
semplot(ht,matrices{2},'r');
semplot(qt,matrices{2},'r');
PlotHVLines(0,'h','k--','linewidth',2);
PlotHVLines(0,'v','k--','linewidth',2);
[templates,correlations,weights] = ActivityTemplates(Restrict(spikes,intervals,'shift','on'),'binsize',0.1);
re = ReactivationStrength(spikes,templates,'step',0.01,'binSize',0.1);
for j=1:size(re,2)-1
[q,qt] = PETH(re(:,[1 1+j]),deltas(:,2),'durations',[-1 1]*2,'nBins',501);
kk = kk+1;
qq{kk,1} = nanmean(q(pre,:),1);
qq{kk,2} = nanmean(q(post,:),1);
qq{kk,3} = i;
%                             qqs{kk,1} = q;
end
q1 = cell2mat(qq(:,1));
q2 = cell2mat(qq(:,2));
PlotColorMap(q1);
PlotColorMap([q1 q2]);
%
nComponents = Accumulate(cell2mat(qq(:,end)));
addSlash = @(x) strrep(x,'_','\_');
sessions = cellfun(addSlash,sessionIDs(nComponents>0),'UniformOutput',0);
lines = cumsum(nComponents(nComponents>0))+0.5;
q = (cell2mat(qq(:,1:2)));
%         z = zBaseline(q,Unfind((length(qt)+1):(length(qt)*2),(length(qt)*2)),2);
z = nanzscore(q,[],2);
z1 = z(:,1:length(qt)); z2 = z(:,(length(qt)+1):end);
smooth = [0 2];
z1 = zscore(z1,[],2); z2 = zscore(z2,[],2);
matrices = {z1,z2,z2-z1};
j=2;
PlotColorMap(Smooth(matrices{j},smooth),~isnan(matrices{j}),'x',qt); PlotHVLines(0,'v','k--','linewidth',2);
set(gca,'ytick',mean([[0;lines(1:end-1)] lines],2),'yticklabel',sessions);
PlotHVLines(lines,'h','w--','linewidth',2);
title([name ', post'])
PlotHVLines(0,'v','k--','linewidth',2);
for i=1:3, subplot(2,3,i); PlotColorMap(matrices{i},'x',qt);
end
for i=1:3, subplot(2,3,i+3); semplot(qt,matrices{i}); end
ylim([-2 4]);
semplot(ht,matrices{2},'r');
semplot(qt,matrices{2},'r');
PlotHVLines(0,'h','k--','linewidth',2);
PlotHVLines(0,'v','k--','linewidth',2);
figure; [templates,correlations,weights] = ActivityTemplates(Restrict(spikes,intervals,'shift','on'),'binsize',0.1);
re = ReactivationStrength(spikes,templates,'step',0.01,'binSize',0.1);
for j=1:size(re,2)-1
[q,qt] = PETH(re(:,[1 1+j]),deltas(:,2),'durations',[-1 1]*2,'nBins',501);
kk = kk+1;
qq{kk,1} = nanmean(q(pre,:),1);
qq{kk,2} = nanmean(q(post,:),1);
qq{kk,3} = i;
%                             qqs{kk,1} = q;
end
q1 = cell2mat(qq(:,1));
q2 = cell2mat(qq(:,2));
PlotColorMap(q1);
PlotColorMap([q1 q2]);
%
nComponents = Accumulate(cell2mat(qq(:,end)));
addSlash = @(x) strrep(x,'_','\_');
sessions = cellfun(addSlash,sessionIDs(nComponents>0),'UniformOutput',0);
lines = cumsum(nComponents(nComponents>0))+0.5;
q = (cell2mat(qq(:,1:2)));
%         z = zBaseline(q,Unfind((length(qt)+1):(length(qt)*2),(length(qt)*2)),2);
z = nanzscore(q,[],2);
z1 = z(:,1:length(qt)); z2 = z(:,(length(qt)+1):end);
smooth = [0 2];
z1 = zscore(z1,[],2); z2 = zscore(z2,[],2);
matrices = {z1,z2,z2-z1};
j=2;
PlotColorMap(Smooth(matrices{j},smooth),~isnan(matrices{j}),'x',qt); PlotHVLines(0,'v','k--','linewidth',2);
set(gca,'ytick',mean([[0;lines(1:end-1)] lines],2),'yticklabel',sessions);
PlotHVLines(lines,'h','w--','linewidth',2);
title([name ', post'])
PlotHVLines(0,'v','k--','linewidth',2);
for i=1:3, subplot(2,3,i); PlotColorMap(matrices{i},'x',qt);
end
for i=1:3, subplot(2,3,i+3); semplot(qt,matrices{i}); end
ylim([-2 4]);
semplot(ht,matrices{2},'r');
semplot(qt,matrices{2},'r');
PlotHVLines(0,'h','k--','linewidth',2);
PlotHVLines(0,'v','k--','linewidth',2);
for i=1:3, subplot(2,3,i); PlotColorMap(matrices{i},'x',qt);
for i=1:3, subplot(2,3,i); PlotColorMap(matrices{i},'x',qt);
end
raly
v1code
fclose('all');
clear all
clc
SleepStateEpisodes = getStruct(basepath,'SleepStateEpisodes');
sws = SleepStateEpisodes.ints.NREMepisode;
cell_metrics = getStruct(basepath,'cell_metrics');
nNeurons = length(cell_metrics.spikes.times);
regions = cell_metrics.brainRegion;
regionNames = unique(regions);
regionCell = zeros(length(regions),1);
for i=1:length(regionNames)
regionCell(strcmp(regions,regionNames{i})) = i;
end
spikesCell = cell_metrics.spikes.times';
PFCindex = [find(strcmp(regionNames,'ILA')) find(strcmp(regionNames,'PFC'))];
HPCindex = find(strcmp(regionNames,'CA1'));
pfc = []; if ~isempty(PFCindex), pfc = sortrows(Group(spikesCell{regionCell==PFCindex})); end
hpc = []; if ~isempty(HPCindex), hpc = sortrows(Group(spikesCell{regionCell==HPCindex})); end
neurons = true;
channel_mapping
cell_metrics = getStruct(basepath,'cell_metrics');
nNeurons = length(cell_metrics.spikes.times);
regions = cell_metrics.brainRegion;
regionNames = unique(regions);
regionCell = zeros(length(regions),1);
for i=1:length(regionNames)
regionCell(strcmp(regions,regionNames{i})) = i;
end
spikesCell = cell_metrics.spikes.times';
PFCindex = [find(strcmp(regionNames,'ILA')) find(strcmp(regionNames,'PFC'))];
HPCindex = find(strcmp(regionNames,'CA1'));
pfc = []; if ~isempty(PFCindex), pfc = sortrows(Group(spikesCell{regionCell==PFCindex})); end
hpc = []; if ~isempty(HPCindex), hpc = sortrows(Group(spikesCell{regionCell==HPCindex})); end
neurons = true;
max(pfc(:,2))
max(hpc(:,2))
% Detect ripples
if exist(fullfile(basepath,[basename '.ripples.events.mat']),'file')
load(fullfile(basepath,[basename '.ripples.events.mat']),'ripples');
else
load(fullfile(basepath,[basename '.session.mat']),'session');
%     rippleChannels = swrChannels('basepath',basepath); % or better yet, pick them on neuroscope (don't forget to add 1! rippleChannels.Ripple_Channel = neuroscope_Channel + 1)
ripples = DetectSWR(rippleChannels,'basepath',basepath,'saveMat',true,'forceDetect',true,'check',true,'useSPW',false);
end
2
% First, remove points that occur in noisy lfp periods
tl = (1:length(lfp))'/1250; % timestamps for lfp
t = featureTs/1250; % timestamps for candidate ripple events
bad = false(size(t));
try % optionally, remove periods of noisy LFP (large deflections)
[clean,~,badIntervals] = CleanLFP([tl,double(lfp(:,2))],'thresholds',[6 Inf],'manual',true);
bad = bad | InIntervals(t,badIntervals);
end
threshold1 = 4;
PlotHVLines([-4 4],'h','k--','linewidth',2);
clf; plot(t,d); hold all; plot(xlim,ones(1,2)*threshold2,'r--');plot(xlim,-ones(1,2)*threshold2,'r--'); ylabel('lfp derivative (z-units'); legend('derivative','threshold2');
disp('Change "threshold2" manually.')
threshold2 = 1;
clf; plot(t,d); hold all; plot(xlim,ones(1,2)*threshold2,'r--');plot(xlim,-ones(1,2)*threshold2,'r--'); ylabel('lfp derivative (z-units'); legend('derivative','threshold2');
threshold2 = 2;
clf; plot(t,d); hold all; plot(xlim,ones(1,2)*threshold2,'r--');plot(xlim,-ones(1,2)*threshold2,'r--'); ylabel('lfp derivative (z-units'); legend('derivative','threshold2');
dbcont
Portion(bad)
[clean,~,badIntervals] = CleanLFP([tl,double(lfp(:,1))],'thresholds',[6 Inf],'manual',true);
bad = bad | InIntervals(t,badIntervals);
threshold2 = 2;
clf; plot(t,d); hold all; plot(xlim,ones(1,2)*threshold2,'r--');plot(xlim,-ones(1,2)*threshold2,'r--'); ylabel('lfp derivative (z-units'); legend('derivative','threshold2');
dbcont
Portion(bad)
idx1 = idx1 & ~bad; idx2 = idx2 & ~bad; % remove bad points from "idx1" and "idx2"
figure(1); % plot a clean figure without these points and check indivudual ripples
clf
plot(swDiffAll(idx1 | idx2),ripPowerAll(idx1 | idx2),'k.','markersize',1); hold on
matrix = [swDiffAll ripPowerAll];
z = bsxfun(@rdivide,bsxfun(@minus,matrix,mean(matrix(~bad,:))),std(matrix(~bad,:)));
scores = nan(size(ripPowerAll,1),1);
while true % click "Ctrl+C" to exit when you feel you're done
figure(1);
colors = Bright(1000);
[x,y,button] = ginput(1);
[~,i] = min(abs(x-swDiffAll.*(1-bad))+abs(y-ripPowerAll.*(1-bad)));
t = featureTs/1250;
tl = (1:length(lfp))'/1250;
figure(2);
clf; interval = t(i) + [-1 1]*0.5; in = tl>interval(1) & tl<interval(end);
subplot(3,1,1);
plot(tl(in) - t(i),lfp(in,:)); title(num2str(t(i)));
hold on; plot(bsxfun(@times,ones(1,2),(Restrict(t,interval) - t(i))),ylim,'k--');
legend('ripple channel','SW channel','noise channel');
subplot(3,1,2);
plot(tl(in) - t(i),double(swDiff(in)));
hold on; plot(bsxfun(@times,ones(1,2),(Restrict(t,interval) - t(i))),ylim,'k--');
subplot(3,1,3);
plot(tl(in) - t(i),double(ripPower0(in)));
hold on; plot(bsxfun(@times,ones(1,2),(Restrict(t,interval) - t(i))),ylim,'k--');
%     x = input( prompt )
commandwindow % select the command window
score = str2double(input('good(1)/bad(0)?','s'));
scores(i) = score; % save scores to optionally save
figure(1);
xlims = xlim; ylims = ylim; clf
% extrapolate the score of each ripple based on the score of the nearest scored ripple
scored = find(~isnan(scores));
distances = sqrt(bsxfun(@minus,z(:,1),z(scored,1)').^2 + bsxfun(@minus,z(:,2),z(~isnan(scores),2)').^2);
[~,nearestNeighbour] = min(distances,[],2);
estimated = scores(scored(nearestNeighbour));
selected = estimated>0.1;
idx1 = selected & ~bad; % final ripples
idx2 = ~selected & ~bad; % non-ripples (yet free from noise as well)
scatter(matrix(~bad,1),matrix(~bad,2),1,estimated(~bad),'filled'); colormap(Bright); set(gca,'CLim',[0 1])
xlim(xlims); ylim(ylims);
hold on;
scatter(swDiffAll(scored),ripPowerAll(scored),20,scores(scored)); scatter(swDiffAll(scored),ripPowerAll(scored),15,scores(scored));
set(get(colorbar,'YLabel'),'String','Estimated ripple score');
xlabel('Sharp wave depth'); ylabel('Ripple power');
title({['Manual scoring for ' basepath],['called with channels: ' num2str(Channels) ' (ripple, sharp wave, and (optionally) noise)']});
if rem(i,10)==0
display('saving...');
save(fullfile(basepath,'DetectSWR_manual_scoring.mat'),'t','swDiffAll','ripPowerAll','bad','scores');
end
end
0
1
0
1
0
1
0
1
0
1
0
% Save your progress!
display('saving...');
save(fullfile(basepath,'DetectSWR_manual_scoring.mat'),'t','swDiffAll','ripPowerAll','idx1','idx2','bad','selected','scores');
saveas(figure(1),fullfile(basepath,'DetectSWR_manual_scoring.fig'));
clf
plot(swDiffAll(idx2),ripPowerAll(idx2),'k.','markersize',1); hold on
plot(swDiffAll(idx1),ripPowerAll(idx1),'r.','markersize',1); legend('non-ripples','ripples');
dbcont
r = ripples.timestamps;
[h,ht] = PETH(hpc(:,1),r(:,1));
clf
PlotColorMap(Shrink(sortby(h,r(:,1)),500,1),'x',ht);
PlotColorMap(Shrink(sortby(h,ripples.peakNormedPower),500,1),'x',ht);
PlotColorMap(Shrink(sortby(h,ripples.RipMax),500,1),'x',ht);
PlotColorMap(Shrink(sortby(h,ripples.RipMax),100,1),'x',ht);
PlotColorMap(Shrink(sortby(h,diff(r,[],2)),100,1),'x',ht);
PlotColorMap(Shrink(sortby(h,diff(r,[],2)),10,1),'x',ht);
lfp = GetAyaLFP(channel);
display(['loaded! @' num2str(toc)]);
[clean,~,badIntervals] = CleanLFP(lfp,'thresholds',[6 10],'manual',true);
4/21.6969
threshold1 = 4;
PlotHVLines([-4 4],'h','k--','linewidth',2);
clf; plot(t,d); hold all; plot(xlim,ones(1,2)*threshold2,'r--');plot(xlim,-ones(1,2)*threshold2,'r--'); ylabel('lfp derivative (z-units'); legend('derivative','threshold2');
disp('Change "threshold2" manually.')
threshold2 = 2;
clf; plot(t,d); hold all; plot(xlim,ones(1,2)*threshold2,'r--');plot(xlim,-ones(1,2)*threshold2,'r--'); ylabel('lfp derivative (z-units'); legend('derivative','threshold2');
dbcont
deltas0 = FindDeltaPeaks(clean);
[h,ht] = PETH(pfc(:,1),deltas(:,2));
[h,ht] = PETH(pfc(:,1),deltas0(:,2));
PlotColorMap(Shrink(sortby(h,-(deltas0(:,5)-deltas0(:,6))),72,1),'x',ht);
semplot(ht,h);
clf
semplot(ht,h);
semplot(ht,h,'durations',[-1 1]*0.2);
PETH(pfc(:,1),deltas(:,2),'durations',[-1 1]*0.2);
PETH(pfc(:,1),deltas0(:,2),'durations',[-1 1]*0.2);
clf
PETH(pfc(:,1),deltas0(:,2),'durations',[-1 1]*0.2);
SaveCustomEvents('deltas.del.evt',deltas0(:,1:3),{'deltas start','delta peak','deltas stop'});a
kkeyboard
load('day6.deltaWaves.events.mat')
load('day6.cell_metrics.cellinfo.mat')
nNeurons = length(cell_metrics.spikes.times);
regions = cell_metrics.brainRegion;
regionNames = unique(regions);
regionCell = zeros(length(regions),1);
for i=1:length(regionNames)
regionCell(strcmp(regions,regionNames{i})) = i;
end
spikesCell = cell_metrics.spikes.times';
PFCindex = [find(strcmp(regionNames,'ILA')) find(strcmp(regionNames,'PFC'))];
HPCindex = find(strcmp(regionNames,'CA1'));
pfc = []; if ~isempty(PFCindex), pfc = sortrows(Group(spikesCell{regionCell==PFCindex})); end
hpc = []; if ~isempty(HPCindex), hpc = sortrows(Group(spikesCell{regionCell==HPCindex})); end
deltas=  deltaWaves.peaks;
deltas=  deltaWaves.peaks; deltas = [deltas deltas[;
deltas=  deltaWaves.peaks; deltas = [deltas deltas];
[h,ht] = PETH(pfc(:,1),deltas0(:,2));
[h,ht] = PETH(pfc(:,1),deltas(:,2));
figure; semplot(ht,h);
[h,ht] = PETH(pfc(:,1),deltas0(:,2));
figure; semplot(ht,h,'r');
semplot(ht,h,'r');
PlotHVLines(0,'v','k--','linewidth',2);
figure; PlotColorMap(Shrink(sortby(h,-(deltas0(:,5)-deltas0(:,6))),72,1),'x',ht);
median((deltas0(:,5)-deltas0(:,6)))
Portion((deltas0(:,5)-deltas0(:,6))>4)
deltas = deltas0(deltas0(:,5)-deltas0(:,6)>4,:); % these thresholds should be manually refined for each session
[h,ht] = PETH(pfc(:,1),deltas(:,2));
figure; PlotColorMap(Shrink(sortby(h,-(deltas(:,5)-deltas(:,6))),72,1),'x',ht);
PlotColorMap(Shrink(sortby(h,-(deltas(:,5)-deltas(:,6))),20,1),'x',ht);
deltas(findmax(deltas(:,5)-deltas(:,6)),2)
sec2min(deltas(findmax(deltas(:,5)-deltas(:,6)),2))
8/202
quantile(deltas(:,5)-deltas(:,6),0.05)
quantile(deltas(:,5)-deltas(:,6),0.95)
semplot(ht,h(deltas(:,5)-deltas(:,6) > 7,:),'r');
semplot(ht,h(deltas(:,5)-deltas(:,6) < 7,:),'b');
deltas(deltas(:,5)-deltas(:,6)>7,:) = [];
deltaWaves.timestamps = deltas(:,[1 3]); deltaWaves.peaks = deltas(:,2);  deltaWaves.peakNormedPower = deltas(:,5); deltaWaves.detectorName = ['channel ' num2str(channel) '(+1), CleanLFP, FindDeltaPeaks, peak-trough>3'];
save(fullfile(basepath,[basename '.deltaWaves.events.mat']),'deltaWaves');
SaveCustomEvents('deltas.del.evt',deltas0(:,1:3),{'deltas start','delta peak','deltas stop'});a
close all
[h,ht] = PETH(ripples.timestamps(:,1),deltas(:,2));
PlotColorMap(Shrink(sortby(h,-(deltas(:,5)-deltas(:,6))),20,1),'x',ht);
semplot(ht,h);
MergePoints = getStruct(basepath,'MergePoints');
basepath
sleep = ConsolidateIntervals(MergePoints.timestamps(cellfun(@(x) any(strfind(lower(x),'sleep')),MergePoints.foldernames),:));
MergePoints = getStruct(basepath,'MergePoints');
sleep = ConsolidateIntervals(MergePoints.timestamps(cellfun(@(x) any(strfind(lower(x),'sleep')),MergePoints.foldernames),:));
pre = InIntervals(deltas,preSleep);
post = InIntervals(deltas,postSleep);
[h,ht] = PETH(ripples.timestamps(:,1),deltas(:,2));
semplot(ht,h(pre,:),'b');
semplot(ht,h(post,:),'r');
MergePoints = getStruct(basepath,'MergePoints');
sleep = ConsolidateIntervals(MergePoints.timestamps(cellfun(@(x) any(strfind(lower(x),'sleep')),MergePoints.foldernames),:));
pre = InIntervals(deltas,sleep(1,:));
post = InIntervals(deltas,sleep(end,:));
[h,ht] = PETH(ripples.timestamps(:,1),deltas(:,2));
semplot(ht,h(pre,:),'b');
semplot(ht,h(post,:),'r');
raly
clear all
2
batchO = StartBatch(@BatchLoadHPCPFCData,'OJR_shortTraining_optoRipples_delayedPFC.batch');
Xo = get(batchO,'UserData'); % opto
XX = {Xs,Xl,Xo};
45
script_OJR_deltaCoupling
45
clear mm0 mm1
for k=1:length(XX) % for each of the three conditions
X = XX{k};
for i=1:size(X,1)
[pfc,hpc,deltas,r,sleep,session,MergePoints,clickers,sws,basepath] = X{i,:};
preSleep = SubtractIntervals(sleep(1,:), SubtractIntervals([0 Inf],sws));
postSleep = SubtractIntervals(sleep(2:end,:), SubtractIntervals([0 Inf],sws));
try
pred = InIntervals(deltas(:,2),preSleep);
postd = InIntervals(deltas(:,2),postSleep);
[hr,ht] = PETH(r(:,2),deltas(:,2),'nBins',501);
%             z = nanzscore([nanmean(hr(pred,:)) nanmean(hr(postd,:))]);
%             mm0{k}(i,:) = z(1:length(ht));
%             mm1{k}(i,:) = z((1:length(ht))+length(ht));
mm0{k}(i,:)  = zscore(nanmean(hr(pred,:)));
mm1{k}(i,:)  = zscore(nanmean(hr(postd,:)));
pre = InIntervals(r(:,2),preSleep);
post = InIntervals(r(:,2),postSleep);
nns{k}(i,1) = sum(InIntervals(deltas(:,2),[r(pre,2)+0.0 r(pre,2)+0.25]));%./dur(preSleep);
nns{k}(i,2) = sum(InIntervals(deltas(:,2),[r(post,2)+0.0 r(post,2)+0.25]));%./dur(postSleep);
end
end
end
open  Xs
open Xo
Unfind(3600,preSleep)
Unfind(3600,postSleep)
Unshift(3600,postSleep)
Unshift(3600,preSleep)
for k=1:length(XX) % for each of the three conditions
X = XX{k};
for i=1:size(X,1)
[pfc,hpc,deltas,r,sleep,session,MergePoints,clickers,sws,basepath] = X{i,:};
if sws(end)==Inf, thisisbad(i,k) = true;e nd
end
end
end
for k=1:length(XX) % for each of the three conditions
X = XX{k};
for i=1:size(X,1)
[pfc,hpc,deltas,r,sleep,session,MergePoints,clickers,sws,basepath] = X{i,:};
if sws(end)==Inf, thisisbad(i,k) = true;end
end
end
end
for k=1:length(XX) % for each of the three conditions
X = XX{k};
for i=1:size(X,1)
[pfc,hpc,deltas,r,sleep,session,MergePoints,clickers,sws,basepath] = X{i,:};
if sws(end)==Inf, thisisbad(i,k) = true;end
end
end
cd N:\OJRproject\OJR49\day3
kkeyboard
SleepScoreMaster(basepath,'noPrompts',true,'rejectchannels',1+[15 7 6 5 64 71 96 114 113 117 118 119 103 120 121 122]); % takes lfp in base 0
SleepScoreMaster(pwd,'noPrompts',true,'rejectchannels',1+[15 7 6 5 64 71 96 114 113 117 118 119 103 120 121 122]); % takes lfp in base 0
%-- 4/13/2022 5:46 PM --%
cd N:\OJRproject\OJR49\day3
SleepScoreMaster(pwd,'noPrompts',true,'rejectchannels',1+[15 7 6 5 64 71 96 114 113 117 118 119 103 120 121 122]); % takes lfp in base 0
t_clus(end)
dbcont
SleepScoreMaster(pwd,'noPrompts',true,'rejectchannels',1+[15 7 6 5 64 71 96 114 113 117 118 119 103 120 121 122]); % takes lfp in base 0
v1code
t_clus(end)
edit SleepScoreMaster.m
st = dbstack;
namestr = st.name
max(thratio)
size(thratio)
display('Calculating Theta Metric above PSS')
%Put the LFP in the right structure format
lfp.data = thLFP;
lfp.timestamps = t_LFP;
lfp.samplingRate = sf_LFP;
%Calculate PSS
[specslope,spec] = bz_PowerSpectrumSlope(lfp,window,window-noverlap,...
'nfreqs',200,'frange',f_all,'IRASA',ThIRASA);
t_thclu = specslope.timestamps;
specdt = 1./specslope.samplingRate;
thFFTspec = specslope.resid';
thFFTspec(thFFTspec<0)=0;
IRASAsmooth_th = spec.IRASAsmooth';
thFFTspec_raw = 10.^spec.amp';
% Remove transients before calculating SW histogram
zFFTspec = NormToInt(spec.amp,'modZ');
totz = NormToInt(abs(sum(zFFTspec,2)),'modZ');
badtimes_TH = find(totz>3);
thFFTfreqs = specslope.freqs';
thfreqs = (thFFTfreqs>=f_theta(1) & thFFTfreqs<=f_theta(2));
thratio = max((thFFTspec(thfreqs,:)),[],1);
size(thratio)
display('Calculating Theta Metric above PSS')
%Put the LFP in the right structure format
lfp.data = thLFP;
lfp.timestamps = t_LFP;
lfp.samplingRate = sf_LFP;
%Calculate PSS
[specslope,spec] = bz_PowerSpectrumSlope(lfp,window,window-noverlap,...
'nfreqs',200,'frange',f_all,'IRASA',ThIRASA);
t_thclu = specslope.timestamps;
specdt = 1./specslope.samplingRate;
thFFTspec = specslope.resid';
thFFTspec(thFFTspec<0)=0;
IRASAsmooth_th = spec.IRASAsmooth';
thFFTspec_raw = 10.^spec.amp';
% Remove transients before calculating SW histogram
zFFTspec = NormToInt(spec.amp,'modZ');
totz = NormToInt(abs(sum(zFFTspec,2)),'modZ');
badtimes_TH = find(totz>3);
thFFTfreqs = specslope.freqs';
thfreqs = (thFFTfreqs>=f_theta(1) & thFFTfreqs<=f_theta(2));
thratio = max((thFFTspec(thfreqs,:)),[],1);
size(thratio)
display('Calculating Theta Metric above PSS')
%Put the LFP in the right structure format
lfp.data = thLFP;
lfp.timestamps = t_LFP;
lfp.samplingRate = sf_LFP;
%Calculate PSS
[specslope,spec] = bz_PowerSpectrumSlope(lfp,window,window-noverlap,...
'nfreqs',200,'frange',f_all,'IRASA',ThIRASA);
t_thclu = specslope.timestamps;
specdt = 1./specslope.samplingRate;
thFFTspec = specslope.resid';
thFFTspec(thFFTspec<0)=0;
IRASAsmooth_th = spec.IRASAsmooth';
thFFTspec_raw = 10.^spec.amp';
% Remove transients before calculating SW histogram
zFFTspec = NormToInt(spec.amp,'modZ');
totz = NormToInt(abs(sum(zFFTspec,2)),'modZ');
badtimes_TH = find(totz>3);
thFFTfreqs = specslope.freqs';
thfreqs = (thFFTfreqs>=f_theta(1) & thFFTfreqs<=f_theta(2));
thratio = max((thFFTspec(thfreqs,:)),[],1);
size(thratio)
display('Calculating Theta Metric above PSS')
%Put the LFP in the right structure format
lfp.data = thLFP;
lfp.timestamps = t_LFP;
lfp.samplingRate = sf_LFP;
%Calculate PSS
[specslope,spec] = bz_PowerSpectrumSlope(lfp,window,window-noverlap,...
'nfreqs',200,'frange',f_all,'IRASA',ThIRASA);
t_thclu = specslope.timestamps;
specdt = 1./specslope.samplingRate;
thFFTspec = specslope.resid';
thFFTspec(thFFTspec<0)=0;
IRASAsmooth_th = spec.IRASAsmooth';
thFFTspec_raw = 10.^spec.amp';
% Remove transients before calculating SW histogram
zFFTspec = NormToInt(spec.amp,'modZ');
totz = NormToInt(abs(sum(zFFTspec,2)),'modZ');
badtimes_TH = find(totz>3);
thFFTfreqs = specslope.freqs';
thfreqs = (thFFTfreqs>=f_theta(1) & thFFTfreqs<=f_theta(2));
thratio = max((thFFTspec(thfreqs,:)),[],1);
specdt = 1./specslope.samplingRate;
size(thratio)
display('Calculating Theta Metric above PSS')
%Put the LFP in the right structure format
lfp.data = thLFP;
lfp.timestamps = t_LFP;
lfp.samplingRate = sf_LFP;
%Calculate PSS
[specslope,spec] = bz_PowerSpectrumSlope(lfp,window,window-noverlap,...
'nfreqs',200,'frange',f_all,'IRASA',ThIRASA);
t_thclu = specslope.timestamps;
specdt = 1./specslope.samplingRate;
thFFTspec = specslope.resid';
thFFTspec(thFFTspec<0)=0;
IRASAsmooth_th = spec.IRASAsmooth';
thFFTspec_raw = 10.^spec.amp';
% Remove transients before calculating SW histogram
zFFTspec = NormToInt(spec.amp,'modZ');
totz = NormToInt(abs(sum(zFFTspec,2)),'modZ');
badtimes_TH = find(totz>3);
thFFTfreqs = specslope.freqs';
thfreqs = (thFFTfreqs>=f_theta(1) & thFFTfreqs<=f_theta(2));
thratio = max((thFFTspec(thfreqs,:)),[],1);
size(thratio)
display('Calculating Theta Metric above PSS')
%Put the LFP in the right structure format
lfp.data = thLFP;
lfp.timestamps = t_LFP;
lfp.samplingRate = sf_LFP;
%Calculate PSS
[specslope,spec] = bz_PowerSpectrumSlope(lfp,window,window-noverlap,...
'nfreqs',200,'frange',f_all,'IRASA',ThIRASA);
t_thclu = specslope.timestamps;
specdt = 1./specslope.samplingRate;
thFFTspec = specslope.resid';
thFFTspec(thFFTspec<0)=0;
IRASAsmooth_th = spec.IRASAsmooth';
thFFTspec_raw = 10.^spec.amp';
% Remove transients before calculating SW histogram
zFFTspec = NormToInt(spec.amp,'modZ');
totz = NormToInt(abs(sum(zFFTspec,2)),'modZ');
badtimes_TH = find(totz>3);
thFFTfreqs = specslope.freqs';
thfreqs = (thFFTfreqs>=f_theta(1) & thFFTfreqs<=f_theta(2));
thratio = max((thFFTspec(thfreqs,:)),[],1);
size(thratio)
display('Calculating Theta Metric above PSS')
%Put the LFP in the right structure format
lfp.data = thLFP;
lfp.timestamps = t_LFP;
lfp.samplingRate = sf_LFP;
%Calculate PSS
[specslope,spec] = bz_PowerSpectrumSlope(lfp,window,window-noverlap,...
'nfreqs',200,'frange',f_all,'IRASA',ThIRASA);
t_thclu = specslope.timestamps;
specdt = 1./specslope.samplingRate;
thFFTspec = specslope.resid';
thFFTspec(thFFTspec<0)=0;
IRASAsmooth_th = spec.IRASAsmooth';
thFFTspec_raw = 10.^spec.amp';
% Remove transients before calculating SW histogram
zFFTspec = NormToInt(spec.amp,'modZ');
totz = NormToInt(abs(sum(zFFTspec,2)),'modZ');
badtimes_TH = find(totz>3);
thFFTfreqs = specslope.freqs';
thfreqs = (thFFTfreqs>=f_theta(1) & thFFTfreqs<=f_theta(2));
thratio = max((thFFTspec(thfreqs,:)),[],1);
size(thratio)
display('Calculating Theta Metric above PSS')
%Put the LFP in the right structure format
lfp.data = thLFP;
lfp.timestamps = t_LFP;
lfp.samplingRate = sf_LFP;
%Calculate PSS
[specslope,spec] = bz_PowerSpectrumSlope(lfp,window,window-noverlap,...
'nfreqs',200,'frange',f_all,'IRASA',ThIRASA);
t_thclu = specslope.timestamps;
specdt = 1./specslope.samplingRate;
thFFTspec = specslope.resid';
thFFTspec(thFFTspec<0)=0;
IRASAsmooth_th = spec.IRASAsmooth';
thFFTspec_raw = 10.^spec.amp';
% Remove transients before calculating SW histogram
zFFTspec = NormToInt(spec.amp,'modZ');
totz = NormToInt(abs(sum(zFFTspec,2)),'modZ');
badtimes_TH = find(totz>3);
thFFTfreqs = specslope.freqs';
thfreqs = (thFFTfreqs>=f_theta(1) & thFFTfreqs<=f_theta(2));
thratio = max((thFFTspec(thfreqs,:)),[],1);
size(thratio)
display('Calculating Theta Metric above PSS')
%Put the LFP in the right structure format
lfp.data = thLFP;
lfp.timestamps = t_LFP;
lfp.samplingRate = sf_LFP;
%Calculate PSS
[specslope,spec] = bz_PowerSpectrumSlope(lfp,window,window-noverlap,...
'nfreqs',200,'frange',f_all,'IRASA',ThIRASA);
t_thclu = specslope.timestamps;
specdt = 1./specslope.samplingRate;
thFFTspec = specslope.resid';
thFFTspec(thFFTspec<0)=0;
IRASAsmooth_th = spec.IRASAsmooth';
thFFTspec_raw = 10.^spec.amp';
% Remove transients before calculating SW histogram
zFFTspec = NormToInt(spec.amp,'modZ');
totz = NormToInt(abs(sum(zFFTspec,2)),'modZ');
badtimes_TH = find(totz>3);
thFFTfreqs = specslope.freqs';
thfreqs = (thFFTfreqs>=f_theta(1) & thFFTfreqs<=f_theta(2));
thratio = max((thFFTspec(thfreqs,:)),[],1);
size(thratio)
[specslope,spec] = bz_PowerSpectrumSlope(lfp,window,window-noverlap,...
'nfreqs',200,'frange',f_all,'IRASA',ThIRASA);
lfp.samplingRate
lfp.timestamps(end)
[specslope,spec] = bz_PowerSpectrumSlope(lfp,window,window-noverlap,...
'nfreqs',200,'frange',f_all,'IRASA',ThIRASA);
[specslope,spec] = bz_PowerSpectrumSlope(lfp,window,window-noverlap,...
'nfreqs',200,'frange',f_all,'IRASA',ThIRASA); % this is the function that sometimes cuts off the
[fList,pList] = matlab.codetools.requiredFilesAndProducts('bz_PowerSpectrumSlope.m');
size(t_clus)
ignoretimeIDX = InIntervals(t_clus,ignoretime) | isnan(broadbandSlowWave) | isnan(thratio);
open fList
fList = fList';
dbcont
SleepScoreMaster(pwd,'noPrompts',true,'rejectchannels',1+[15 7 6 5 64 71 96 114 113 117 118 119 103 120 121 122]); % takes lfp in base 0
v1code
2
for k=1:length(XX) % for each of the three conditions
X = XX{k};
for i=1:size(X,1)
[pfc,hpc,deltas,r,sleep,session,MergePoints,clickers,sws,basepath] = X{i,:};
if sws(end)==Inf, thisisbad(i,k) = true;end
end
end
cd 'N:\OJRproject\OJR49\day7'
SleepScoreMaster(pwd,'noPrompts',true,'rejectchannels',1+[22 23 31 71 96 109 114 113 117 118 119 103 120 121 122]); % takes lfp in base 0
t_clus(end)
size(t_clus)
dbquit
SleepScoreMaster(pwd,'noPrompts',true,'rejectchannels',1+[22 23 31 71 96 109 114 113 117 118 119 103 120 121 122]); % takes lfp in base 0
%-- 4/13/2022 6:05 PM --%
cd 'N:\OJRproject\OJR49\day7'
SleepScoreMaster(pwd,'noPrompts',true,'rejectchannels',1+[22 23 31 71 96 109 114 113 117 118 119 103 120 121 122]); % takes lfp in base 0
t_clus(end)
size(t_clus_
size(t_clus)
SleepScoreMaster(pwd,'noPrompts',true,'rejectchannels',1+[22 23 31 71 96 109 114 113 117 118 119 103 120 121 122]); % takes lfp in base 0
v1code
edit startup
5
2
32
close all
Hide('none');
2
k=1
X = XX{k};
i=1
[pfc,hpc,deltas,r,sleep,session,MergePoints,clickers,sws,basepath] = X{i,:};
preSleep = SubtractIntervals(sleep(1,:), SubtractIntervals([0 Inf],sws));
postSleep = SubtractIntervals(sleep(2:end,:), SubtractIntervals([0 Inf],sws));
try preSleep = SubtractIntervals(preSleep,[Unshift(3600,preSleep) Inf]); end
try postSleep = SubtractIntervals(postSleep,[Unshift(3600,preSleep) Inf]); end
try
pred = InIntervals(deltas(:,2),preSleep);
postd = InIntervals(deltas(:,2),postSleep);
[hr,ht] = PETH(r(:,2),deltas(:,2),'nBins',501);
%             z = nanzscore([nanmean(hr(pred,:)) nanmean(hr(postd,:))]);
%             mm0{k}(i,:) = z(1:length(ht));
%             mm1{k}(i,:) = z((1:length(ht))+length(ht));
mm0{k}(i,:)  = zscore(nanmean(hr(pred,:)));
mm1{k}(i,:)  = zscore(nanmean(hr(postd,:)));
pre = InIntervals(r(:,2),preSleep);
post = InIntervals(r(:,2),postSleep);
ns{k}(i,1) = sum(InIntervals(deltas(:,2),[r(pre,2)+0.0 r(pre,2)+0.25]))./dur(preSleep);
ns{k}(i,2) = sum(InIntervals(deltas(:,2),[r(post,2)+0.0 r(post,2)+0.25]))./dur(postSleep);
end
sleep
2
figure; anovabar(q(:,1),q(:,2));
anovabar(q(:,1)-1,q(:,2));
kruskalbar(q(:,1)-1,q(:,2));