folder = pwd;
[~,basename] = fileparts(folder);

% load('LED_timestamps.mat','LED');
load([basename '.cell_metrics.cellinfo.mat'])
load([basename '.chanCoords.channelInfo.mat'])
load([basename '.MergePoints.events.mat'])

regions = cell_metrics.brainRegion;
regionNames = unique(regions)
regionCell = zeros(length(regions),1);
for i=1:length(regionNames)
    regionCell(strcmp(regions,regionNames{i})) = i;
end


depth = chanCoords.y(cell_metrics.maxWaveformCh1)';
matrix = dlmread('stimuli_synced.times');
% Windows is imprecise when saving the "last modified" date of the file:
% it says the file was modified between 1.72s to 0.72 before it was actually finished!
% The real time is therefore AT LEAST 0.72s before the saved time
matrix(:,1) = matrix(:,1); % I'm subtracting 1, although the real number (between 1.72 and 0.72) is unknown
matrix(isnan(matrix(:,1)),:) = [];

spikesCell = cell_metrics.spikes.times';

% for i=1:length(spikesCell)
%     spikesCell{i} = spikesCell{i}*1;
% end
% matrix(:,1) = matrix(:,1)*1;

nUnits = length(spikesCell);
% ok = LED(:,1)>MergePoints.timestamps(end)-1440;
stim = matrix(diff([0;matrix(:,2)])~=0 & matrix(:,2)>-1,:);
% stim = matrix(matrix(:,2)==-2);

% stim = stim(10:20:end);
clf
for i=1:nUnits
    [h,ht]  = PETH(spikesCell{i}*1,stim(:,1)*1,'durations',[-3 3]*2,'nBins',501);
    % PlotColorMap(Smooth(Shrink(h,1,1),0),'x',ht)
    m(i,:) = mean(h,1);
end

z = zscore(m,[],2);
score = nanmean(z(:,InIntervals(ht,[0 0.2])),2) - nanmean(z(:,InIntervals(ht,[-1 0])),2);

subplot(1,2,1);
PlotColorMap(sortby(Smooth(zscore(m,[],2),[0 2]),(depth)),'x',ht);
PlotHVLines([-1:2],'v','k--','linewidth',2);
PlotHVLines(sum(depth<-900)+0.5,'h','w--','linewidth',2);
% xlim([-1 1]*5);
set(gca,'ytick',1:max(ylim),'yticklabel',sort(depth));
ylabel('unit depth (a.u. as screw turns were likely not taken into account)');
xlabel('time from LED pulse (s)');
clabel('Firing rate (z units)');
clim([-1 1]*1.5);
set(gca,'fontsize',12);


subplot(2,2,2);
semplot(ht,(zscore(m(depth<-750,:),[],2)),'k', 2);
semplot(ht,(zscore(m(depth>-750,:),[],2)),'r', 2);
y = ylim;
PlotHVLines(0,'v','k--','linewidth',2);
ylim(y);
% xlim([-1 1]*1);
ylabel('Mean firing rate (z units)');
xlabel('time from LED pulse (s)');
handle = get(gca,'children');
legend(handle([2 4]),'putative V1','putative HPC');
legend('box','off');
set(gca,'fontsize',12);

subplot(2,2,4);
cla
semplot(ht,abs(zscore(m(depth<-750,:),[],2)),'k', 2);
semplot(ht,abs(zscore(m(depth>-750,:),[],2)),'r', 2);
y = ylim;
PlotHVLines(0,'v','k--','linewidth',2);
ylim(y);
% xlim([-1 1]*1);
ylabel('Mean (absolute) change from baseline');
xlabel('time from LED pulse (s)');
handle = get(gca,'children');
legend(handle([2 4]),'putative V1','putative HPC');
legend('box','off');
set(gca,'fontsize',12);

% SaveFig('Effect');

%%
u = unique(matrix(:,2)); u(u==-1) = [];
stims = matrix(diff([0;matrix(:,2)])~=0,:);
hs = {};
for i=1:nUnits
    for j=1:length(u)
        angle = u(j);
        [h,ht]  = PETH(spikesCell{i}*1,stims(stims(:,2)==angle)*1,'durations',[-1 1]*5,'nBins',1001);
        hs{i,1}(j,:) = nanmean(h);
    end
end

stimIntervals = bsxfun(@plus,stims(stims(:,2)~=-1,1),[0 1]);
intAngles = stims(stims(:,2)~=-1,2);
lAngles = zeros(size(intAngles,1),length(u));
for i=1:size(lAngles,1)
    lAngles(i,:) = 1 - abs(wrap((intAngles(i)-u)/90*pi))/(pi); % true angle is 1, farthest angle (90deg.) is 0 (180deg from correct angle is identical, which is why we divide by 90 and not 180deg)
    lAngles(i,:) = intAngles(i)==u; % true angle is 1, farthest angle (90deg.) is 0 (180deg from correct angle is identical, which is why we divide by 90 and not 180deg)
end
lAngles(intAngles==-2,2:end) = 0; lAngles(intAngles~=-2,1) = 0;
counts = nan(size(stimIntervals,1),nUnits);
for i=1:nUnits
    counts(:,i) = CountInIntervals(spikesCell{i}*1,stimIntervals*2);
end

%% Train it one one half of the data and test it in the other half

rng(0); [~,order] = sort(rand(size(intAngles)));
trainingHalf = order(1:end/2); order(1:end/2) = []; testingHalf = order;
clear weights
for j=1:length(u)
    weights(:,j) = glmfit(counts(trainingHalf,:),lAngles(trainingHalf,j));
end
guess = bsxfun(@plus,counts(testingHalf,:)*weights(2:end,:),weights(1,:));

figure(3);
clf
q = zscore(guess);
for i=1:length(u)
    subplot(3,3,i);
    ok = intAngles(testingHalf)==u(i);
    bar(nanmean(q(ok,:)))
    hold all
    bar(i,nanmean(q(ok,i)),'facecolor',[0.9 0.7 0.2]);
    set(gca,'xtick',1:max(xlim),'xticklabel',u);
    ylabel('Prediction (z-units)');
    if i==1
        title(['For blank screen trials']);
    else
        title(['For ' num2str(u(i)) ' deg trials']);
    end
end
EquateScales

% A different representation
figure(2)

truth = FindClosest(u,intAngles(testingHalf));
[~,prediction] = max(guess,[],2);

subplot(1,2,1);
accuracy = Accumulate([truth prediction]);
PlotColorMap((accuracy)./sum(accuracy,2)*100)
ylabel('If the trial was...');
set(gca,'ytick',1:max(ylim),'yticklabel',u);
xlabel('... then the GLM would classify it as a [X degrees] trial')
ylabel('If the trial was a [Y degrees] trial, ...');
clabel('Probabity of guess (%)');
title('Classification matrix');
axis square
set(gca,'box','off','fontsize',12)
accuracyScore = [sum(diag(accuracy))/sum(accuracy(:)) sum(diag(accuracy(2:end,2:end)))/sum(sum(accuracy(2:end,2:end)))];

subplot(1,2,2);
ok = prediction>1;
accuracy = Accumulate([truth(ok) prediction(ok)]);
PlotColorMap((accuracy)./sum(accuracy,2)*100)
ylabel('If the trial was...');
set(gca,'ytick',1:max(ylim),'yticklabel',u);
xlabel('... then the GLM would classify it as a [X degrees] trial')
ylabel('If the trial was a [Y degrees] trial, ...');
clabel('Probabity of guess (%)');
axis square
set(gca,'box','off','fontsize',12)
title('Same as left but ignoring trials classified as blank');
set(gca,'xtick',1:max(xlim),'xticklabel',u);

%% Same but for 0 vs 90 classification only

lAngles = (intAngles==-2)+1-1;
lAngles(:,2) = ismember(intAngles,[157.5 0 22.5]);
lAngles(:,3) = ismember(intAngles,[67.5 90 112.5]);

rng(0); [~,order] = sort(rand(size(intAngles)));
trainingHalf = order(1:end/2); order(1:end/2) = []; testingHalf = order;
clear weights
for j=1:size(lAngles,2)
    weights(:,j) = glmfit(counts(trainingHalf,:),lAngles(trainingHalf,j));
end
guess = bsxfun(@plus,counts(testingHalf,:)*weights(2:end,:),weights(1,:));

clf
q = zscore(guess);
for i=1:size(lAngles,2)
    subplot(3,3,i);
    ok = lAngles(testingHalf,i) == 1;  
    bar(nanmean(q(ok,:)))
    hold all
    bar(i,nanmean(q(ok,i)),'facecolor',[0.9 0.7 0.2]);
    set(gca,'xtick',1:max(xlim)-1,'xticklabel',u([1 2 6]));
    ylabel('Prediction (z-units)');
    if i==1
        title(['For blank screen trials']);
    elseif i==2;
        title('For trials around 0 degrees');
    else
        title('For trials around 90 degrees');
%         title(['For ' num2str(u(i)) ' deg trials']);
    end
end
EquateScales

% A different representation
% figure(2)

truth = FindClosest(u,intAngles(testingHalf));
[~,prediction] = max(guess,[],2);

subplot(2,2,3);
accuracy = Accumulate([truth prediction]);
PlotColorMap((accuracy)./sum(accuracy,2)*100)
ylabel('If the trial was...');
set(gca,'ytick',1:max(ylim),'yticklabel',u);
set(gca,'xtick',1:max(xlim),'xticklabel',u([1 2 6]));
xlabel('... then the GLM would classify it as a [X degrees] trial')
ylabel('If the trial was a [Y degrees] trial, ...');
clabel('Probabity of guess (%)');
title('Classification matrix');
axis square
set(gca,'box','off','fontsize',12)

subplot(2,2,4);
ok = prediction>1;
accuracy = Accumulate([truth(ok) prediction(ok)]);
PlotColorMap((accuracy)./sum(accuracy,2)*100)
ylabel('If the trial was...');
set(gca,'ytick',1:max(ylim),'yticklabel',u);
xlabel('... then the GLM would classify it as a [X degrees] trial')
ylabel('If the trial was a [Y degrees] trial, ...');
clabel('Probabity of guess (%)');
axis square
set(gca,'box','off','fontsize',12)
title('Same as left but ignoring trials classified as blank');
set(gca,'xtick',1:max(xlim),'xticklabel',u([1 2 6]));


%%
figure(1); clf
spikes = Group(spikesCell{:});
shuffled = spikes; shuffled(:,2) = Scramble(shuffled(:,2));
for i=1:nUnits,shuffledCell{i} = shuffled(shuffled(:,2)==i); end

bins = Bins(0,MergePoints.timestamps(end,end),1,0.01);
for i=1:nUnits
    bCounts(:,i) = CountInIntervals(spikesCell{i},bins);
end
decoded = [ones(size(bCounts(:,1))) bCounts]*weights;
for i=1:size(decoded,2)
    [h,ht] = PETH([t decoded(:,i)],r(:,1),'durations',[-1 1]*10);
    hm(i,:) = nanmean(h);
end
for i=1:size(decoded,2)
    [h,ht] = PETH([t decoded(:,i)],r(:,1),'durations',[-1 1]*10);
    dh{i,1} = h;
end
hh = cat(3,dh{:});
[~,m] = max(decoded,[],2); % most likely position
mTranspose = m'; 
for i=1:9
    for j=1:9
        transitions(i,j) = numel(strfind(mTranspose(:)',[i j]));
        this = nan(100,1);
        for k=1:100
            rng(k); scrambled = Scramble(mTranspose);
        this(k,1) = numel(strfind(scrambled(:)',[i j]));
        end
        ctrlTransitions(i,j) = nanmean(this);
    end
end
d = @(x,y) (x-y)./(x+y);
k = [1:9 2:9]; PlotColorMap(d(transitions(k,k),ctrlTransitions(k,k)))
clims([-1 1]);
axis square
title(datestr(clock))

for i=1:9
    for j=1:9
        tTransitions{i,j} = t(strfind(mTranspose(:)',[i j]));
    end
end


%% Save figures
u = unique(matrix(:,2)); u(u==-2) = [];
ok = diff([0;matrix(:,2)])~=0 & matrix(:,end)>0;
stims = matrix(ok,:);
d = @(x,y) (x-y)./(x+y);
for j=1:nUnits
% for j=2
    i=j;
% for i=4
    clf
    rng(100*i);
%     i=1;
    for j=1:length(u)
        angle = u(j);
        [h,ht]  = PETH(spikesCell{i},stims(stims(:,2)==angle)*1+0.05,'durations',[-1 1]*5,'nBins',2001);
        if j==1, hs{i,1} = []; hs1{i,1} = []; end
        hs{i,1}(j,:) = nanmean(h);
        [~,order] = sort(rand(size(h,1),1));
        hs1{i,1}(j,:) = nanmean(h(order(1:ceil(end/2)),:)); order(1:ceil(end/2)) = [];
        hs1{i,1}(j+length(u),:) = nanmean(h(order,:));
    end
    smoothed = Smooth(hs1{i},[0 2]);
    tuningCurve = zscore(max(smoothed(:,InIntervals(ht,[00.05 0.15])),[],2));
    tuningCurves{i,1} = tuningCurve;
    score(i,2) = sum(tuningCurve(2:9).*tuningCurve((2:9)+9));
    subplot(2,1,1);
    PlotColorMap(Smooth(hs1{i}([1 10 2 11 3 12 4 13 5 14 6 15 7 16 8 17 9 18],:),[0 2]),'x',ht);
    set(gca,'ytick',1:max(ylim),'yticklabel',repelem(u,2));
    
    PlotColorMap(Smooth(hs1{i},[0 2]),'x',ht);
    set(gca,'ytick',1:max(ylim),'yticklabel',repmat(u,2,1));
    PlotHVLines([0:2],'v','k--','linewidth',2);
    xlim([-1 3]);
    xlabel('time from stimulus (stimulus=0s, flip=1s, off=2s)');
    ylabel({'angle (-2=OFF)','data split in half (analysed independently)'});
    title(['neuron ' num2str(i) ' recorded at relative depth ' num2str(depth(i))]);
    set(gca,'fontsize',12,'box','off');
    
    
    subplot(2,1,2);
    semplot(ht,hs{i});
    PlotHVLines([0:2],'v','k--','linewidth',2);
    xlim([-1 3]);
    title(score(i,:));
    drawnow
%     pause(2);
%     SaveFig(['y:\V1test\figures\' [projectName '_' dayName] '_unit_' num2str(i)]);
end

%%

load([basename '.ripples.events.mat']);
r = ripples.timestamps(:,1:2);
isHPC = find(depth<-900);

spikes = Group(spikesCell{:});
hpc = spikes(ismember(spikes(:,2),isHPC),:);
s = spikes(~ismember(spikes(:,2),isHPC),1);
bins = Bins(0,MergePoints.timestamps(end,end),1,0.01);
t = mean(bins,2);
[h,~] = hist(s,t);
silence = t(FindInterval(Smooth(h(:),1)<0.1));
silence(diff(silence,[],2)<0.05,:) = [];

f = InIntervals(r(:,1),silence(:,1)-[0.2 0]);
ff = InIntervals(silence(:,1),[r(f) r(f)+0.4]);

clf
PETH(silence(:,1),r(:,1),'durations',[-1 1]*1)


%%
dt = 0.1;
bins = Bins(0,spikes(end,1),0.2,dt);
counts = CountInIntervals(s,bins);
counts(:,2) = CountInIntervals(hpc,bins);

clf
inSWS = InIntervals(bins(:,1),SleepState.ints.NREMstate); inSWS(:) = true;
in = bins(:,1)<stim(1) & inSWS; [x,xt] = xcorr(counts(in,1),counts(in,2),'coeff');
plot(xt*dt,x);
hold all
in = InIntervals(bins(:,1),[stim(1) stim(end,1)]); [x,xt] = xcorr(counts(in,1),counts(in,2),'coeff');
plot(xt*dt,x,'linewidth',2);
in = bins(:,1)>stim(end,1) & inSWS; [x,xt] = xcorr(counts(in,1),counts(in,2),'coeff'); 
plot(xt*dt,x);
xlim([-1 1]*0.5)


%%

clf

for i=1:length(spikesCell)
    subplot(4,ceil(length(spikesCell)/4),i);
    response = CountInIntervals(spikesCell{i},bsxfun(@plus,stims(:,1),[0.05 0.3]));
    
    mm = nan(size(response,1),length(u)); mm(sub2ind(size(mm),(1:size(response,1))',FindClosest(u,stims(:,2)))) = response;
    mm(stims(:,end)==3,:) = [];
    mm = sort(mm); mm(sum(~isnan(mm),2)==0,:) = [];
    preferred(i,1) = u(findmax(nanmean(mm)));
    dd = u; dd(2:end) = wrap((dd(2:end)-preferred(i,1))/180*pi);
    
    
    colors = {'','k','b','r'};
    for k=2:4
        mm = nan(size(response,1),length(u)); mm(sub2ind(size(mm),(1:size(response,1))',FindClosest(u,stims(:,2)))) = response;
        mm(stims(:,end)~=k,:) = [];
        mm = sort(mm); mm(sum(~isnan(mm),2)==0,:) = [];
        
        x = dd(sum(~isnan(mm))>0); x = [x; x(2:end)+4; x(2:end)+8];
        x = u(sum(~isnan(mm))>0); x = [x; x(2:end)+180; x(2:end)+360];
        y = mm(:,sum(~isnan(mm))>0); y(:) = nanzscore(y(:)); y = [y(:,1) repmat(y(:,2:end),1,3)];
        semplot(x,y,colors{k},sum(sum(~isnan(mm))>0)/9);
    end
    
    PlotHVLines([0 180 360]+preferred(i,1),'v','k--');
    xlim([0 360]);
    title(i);
end


%%

templates = ActivityTemplates(Restrict(spikes,MergePoints.timestamps(3,:)));
re = ReactivationStrength(spikes,templates);

%%

[q,qt] = PETH(re(:,[1 1+1]),stims(:,1));
for i=1:9
    qq(i,:) = nanmean(q(idx==i,:));
end
PlotColorMap(qq,'x',qt);

















