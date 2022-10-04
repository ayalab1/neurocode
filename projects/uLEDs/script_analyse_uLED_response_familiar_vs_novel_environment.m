folder = pwd;

[parentFolder,dayName] = fileparts(folder);
[~,projectName] = fileparts(parentFolder);
sessionID = [projectName '_' dayName];
cd(folder);
load(fullfile(folder,[dayName '.session.mat']),'session');
load(fullfile(folder,[dayName '.ripples.events.mat']),'ripples');
load(fullfile(folder,[dayName '.MergePoints.events.mat']),'MergePoints');
load(fullfile(folder,[dayName '.cell_metrics.cellinfo.mat']),'cell_metrics');
r = ripples.timestamps(:,1:2);
nNeurons = length(cell_metrics.spikes.times);
regions = cell_metrics.brainRegion;
regionNames = unique(regions)
regionCell = zeros(length(regions),1);
for i=1:length(regionNames)
    regionCell(strcmp(regions,regionNames{i})) = i;
end
spikesCell = cell_metrics.spikes.times';

if ~exist(fullfile(folder,[dayName '.uLEDPulses.events.mat']),'file')
    load(fullfile(folder,[dayName '.pulses.events.mat']),'pulses');
    load(fullfile(folder,'digitalIn.events.mat'),'digitalIn');
    uLEDPulses = combineULEDPulses('analogPulses',pulses,'digitalPulses',digitalIn);
    save(fullfile(folder,[dayName '.uLEDPulses.events.mat']),'uLEDPulses');
else
    load(fullfile(folder,[dayName '.uLEDPulses.events.mat']),'uLEDPulses');
end
LED = uLEDPulses.timestamps;
[~,~,id] = unique([uLEDPulses.shank uLEDPulses.LED],'rows');

% for each cell, find optimal LED id
intervals = bsxfun(@plus,LED,[1 -1]/1000);
ctrlIntervals = bsxfun(@minus,intervals,diff(LED,[],2));

for i=1:length(spikesCell)
    for j=1:max(id)
        response = CountInIntervals(spikesCell{i},intervals(id==j,:));
        ctrlResponse = CountInIntervals(spikesCell{i},ctrlIntervals(id==j,:));
        mr(i,j) = mean(response);
        cr(i,j) = mean(ctrlResponse);
        p(i,j) = signrank(response,ctrlResponse,'tail','right');
    end
end
gain = (mr-cr)./(mr+cr);
[~,LEDidx] = max(gain,[],2);
nonResponsive = p(sub2ind(size(p),(1:length(spikesCell))',LEDidx))>0.05/max(id);

%%
novel = MergePoints.timestamps(find(~cellfun(@isempty,strfind(MergePoints.foldernames,'novel'))),:);
familiar = MergePoints.timestamps(find(~cellfun(@isempty,strfind(MergePoints.foldernames,'familiar'))),:);

for i=1:length(spikesCell)
    events = LED(id==LEDidx(i),1);
    response = CountInIntervals(spikesCell{i},intervals(id==LEDidx(i),:));
    ctrlResponse = CountInIntervals(spikesCell{i},ctrlIntervals(id==LEDidx(i),:));
    pre = events<novel(1);
    post = events>novel(1);
    mr2(i,1:2) = [mean(response(pre)) mean(response(post))];
    cr2(i,1:2) = [mean(ctrlResponse(pre)) mean(ctrlResponse(post))];
end

gain2 = (mr2-cr2)./(mr2+cr2);

% Get firing rates in familiar and novel environments:
spikes = sortrows(Group(spikesCell{:}));
fr = [Accumulate(spikes(InIntervals(spikes(:,1),familiar),2))/sum(diff(familiar,[],2)),...
    Accumulate(spikes(InIntervals(spikes(:,1),novel),2))/sum(diff(novel,[],2))];
novelPreference = (fr(:,2)-fr(:,1))./sum(fr,2);

PlotColorMap(sortby(gain2,novelPreference));

%%
i=2;
clf
subplot(2,2,1);
plot(fr);
PlotHVLines(i,'v','k--','linewidth',2);
title(novelPreference(i));

subplot(2,2,3);
plot(gain2);
PlotHVLines(i,'v','k--','linewidth',2);
title(gain2(i,:));

subplot(2,2,2);
events = LED(id==LEDidx(i),1);
[h,ht] = PETH(spikesCell{i},events,'nBins',501,'durations',[-1 1]*0.050);
pre = events<novel(1);
post = events>novel(1);
semplot(ht,h(pre,:),'k',1);
semplot(ht,h(post,:),'r',1);
title(i);

subplot(2,2,4);
plot(novelPreference(~nonResponsive),diff(cr2(~nonResponsive,:),[],2),'.','markersize',20);
title(p(i,LEDidx(i)));
hold on
plot(novelPreference(i),diff(gain2(i,:),[],2),'r.','markersize',20);
PlotHVLines(0,'h','k--','linewidth',2);
PlotHVLines(0,'v','k--','linewidth',2);
ylabel('post gain - pre gain');
xlabel('novelty preference');


%%
smooth = 2;
for i=1:length(spikesCell)
    events = LED(id==LEDidx(i),1);
    [h,ht] = PETH(spikesCell{i},events,'nBins',501,'durations',[-1 1]*0.050);
    h = h/diff(ht(1:2));
    pre = events<novel(1);
    post = events>novel(1);
    hPre(i,:) = nanmean(h(pre,:));
    hPost(i,:) = nanmean(h(post,:));
end
hd = hPost-hPre;
clf
subplot(1,2,1);
PlotColorMap(Smooth(sortby(hd(:,:),novelPreference(:)),[0 smooth]),'x',ht);
ylabel({'cell ID (ordered from firing rate difference','(familiar preferring -> white line -> novel preferring)'});
PlotHVLines([0 0.02],'v','k--','linewidth',2);
PlotHVLines(sum(novelPreference<0)+0.5,'h','w--','linewidth',2);
clabel('firing rate Post-Pre difference (Hz)');
xlabel('time from LED pulse (s)');
set(gca,'fontsize',12,'box','off');
title([projectName '_' dayName]);

subplot(1,2,2);
semplot(ht,hd(novelPreference<0,:),'b',smooth);
semplot(ht,hd(novelPreference>0,:),'r',smooth);
PlotHVLines([0 0.02],'v','k--','linewidth',2);
ylabel('firing rate Post-Pre difference (Hz)');
xlabel('time from LED pulse (s)');
handle = get(gca,'children');
legend(handle([5 3]),'familiar preferring','novel preferring');
set(gca,'fontsize',12,'box','off');
title([projectName '_' dayName]);

multiSaveFig(['uLED_responses_for_fam_vs_novel_preferring_cells-' projectName '_' dayName])

%% Same but in ripples

smooth = 5;

for i=1:length(spikesCell)
    events = LED(id==LEDidx(i),1);
    in = IntervalsIntersect(intervals(id==LEDidx(i),:),r);
    events = events(in);
    [h,ht] = PETH(spikesCell{i},events,'nBins',501,'durations',[-1 1]*0.05);
    h = h/diff(ht(1:2));
    pre = events<novel(1);
    post = events>novel(1);
    hPre(i,:) = nanmean(h(pre,:));
    hPost(i,:) = nanmean(h(post,:));
end
hd = hPost-hPre;
clf
subplot(1,2,1);
PlotColorMap(Smooth(sortby(hd(:,:),novelPreference(:)),[0 smooth]),'x',ht);
ylabel({'cell ID (ordered from firing rate difference','(familiar preferring -> white line -> novel preferring)'});
PlotHVLines([0 0.02],'v','k--','linewidth',2);
PlotHVLines(sum(novelPreference<0)+0.5,'h','w--','linewidth',2);
clabel('firing rate Post-Pre difference (Hz)');
xlabel('time from LED pulse (s)');
set(gca,'fontsize',12,'box','off');
title([projectName '_' dayName ' uLED restricted to ripples']);
clim([-1 1]*30);

subplot(1,2,2);

semplot(ht,hd(novelPreference<0,:),'b',smooth);
semplot(ht,hd(novelPreference>0,:),'r',smooth);
PlotHVLines([0 0.02],'v','k--','linewidth',2);
ylabel('firing rate Post-Pre difference (Hz)');
xlabel('time from LED pulse (s)');
handle = get(gca,'children');
legend(handle([5 3]),'familiar preferring','novel preferring');
set(gca,'fontsize',12,'box','off');
title([projectName '_' dayName ' uLED restricted to ripples']);

% multiSaveFig(['uLED_responses_inRipples_for_fam_vs_novel_preferring_cells-' projectName '_' dayName])

%%












