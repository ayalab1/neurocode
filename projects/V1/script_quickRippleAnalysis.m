basepath = pwd;
[parentFolder,basename] = fileparts(basepath);
[~,projectName] = fileparts(parentFolder);
sessionID = [projectName '-' basename];
ripples = DetectSWR(1+[24 17],'saveMat',true,'check',true,'useSPW',false);
% figure; plot(tl,(zscore(double(lfp(:,2))))); Browse('all');
load(fullfile(basepath,[basename '.cell_metrics.cellinfo.mat']),'cell_metrics');

nNeurons = length(cell_metrics.spikes.times);
regions = cell_metrics.brainRegion;
regionNames = unique(regions);
regionCell = zeros(length(regions),1);
for i=1:length(regionNames)
    regionCell(strcmp(regions,regionNames{i})) = i;
end
spikesCell = cell_metrics.spikes.times';
CA2index = [find(strcmp(regionNames,'CA2'))];
CA1index = find(strcmp(regionNames,'CA1'));
ca2 = []; if ~isempty(CA2index), ca2 = sortrows(Group(spikesCell{regionCell==CA2index})); end
ca1 = []; if ~isempty(CA1index), ca1 = sortrows(Group(spikesCell{regionCell==CA1index})); end
hpc = [ca1; ca2(:,1) ca2(:,2)+max(ca1(:,2))];

%%
try load(fullfile(basepath,[basename '.MergePoints.events.mat']),'MergePoints');
    postID = find(cellfun(@(x) ~isempty(strfind(x,'post')), MergePoints.foldernames));
    sleep = ConsolidateIntervals(MergePoints.timestamps(cellfun(@(x) any(strfind(lower(x),'sleep')),MergePoints.foldernames),:));
catch
    MergePoints = [];
    display('No merge points saved. Keeping all ripples (impossible to restrict to post-task sleep');
end

r = ripples.timestamps;
pre = InIntervals(r(:,1),sleep(1,:));
post = InIntervals(r(:,1),sleep(2,:));
nBins = 101;
m = nan(nNeurons,nBins); mPre = nan(nNeurons,nBins); mPost = nan(nNeurons,nBins);

for i=1:nNeurons
    [h,ht] = PETH(hpc(hpc(:,2)==i),ripples.timestamps(:,1),'nBins',nBins);
    m(i,:) = mean(h);
    mPre(i,:) = mean(h(pre,:));
    mPost(i,:) = mean(h(post,:));
end

z = [mPre mPost];
zPre = (mPre - mean(z,2))./std(z,[],2);
zPost = (mPost - mean(z,2))./std(z,[],2);
zDiff = zPost - zPre;
indexCA1 = 1:max(ca1(:,2));
indexCA2 = max(ca1(:,2))+1:max(hpc(:,2));

clf
subplot(2,3,1);
PlotColorMap(zPre,'x',ht); 
set(gca,'TickDir','out');
ylabel('Neuron ID (CA1 -> CA2)'); xlabel('time from ripple start (s)');
PlotHVLines(max(ca1(:,2))+0.5,'h','w--','linewidth',2);
clabel('firing rate (z units)');
title('Pre sleep ripples')

subplot(2,3,2);
PlotColorMap(zPost,'x',ht); 
set(gca,'TickDir','out');
ylabel('Neuron ID (CA1 -> CA2)'); xlabel('time from ripple start (s)');
PlotHVLines(max(ca1(:,2))+0.5,'h','w--','linewidth',2);
clabel('firing rate (z units)')
title('Post sleep ripples')
clims

subplot(2,3,3);
PlotColorMap(zDiff,'x',ht); 
set(gca,'TickDir','out');
ylabel('Neuron ID (CA1 -> CA2)'); xlabel('time from ripple start (s)');
PlotHVLines(max(ca1(:,2))+0.5,'h','w--','linewidth',2);
clabel('firing rate difference (z units)')
ColorMap(gca,[0 0 1],[1 1 1],[1 0 0]);
title('Post-Pre')
clim(max(abs(clim))*[-1 1]);

subplot(2,3,4);
semplot(ht,zPre(indexCA1,:),'k')
semplot(ht,zPre(indexCA2,:),'b')
y = ylim;

subplot(2,3,5);
semplot(ht,zPost(indexCA1,:),'k')
semplot(ht,zPost(indexCA2,:),'b')
PlotHVLines(0,'h','k--','linewidth',2); PlotHVLines(0,'v','k--','linewidth',2)
y = [min(y(1),min(ylim)),max(y(2),max(ylim))];
ylim(y); PlotHVLines(0,'h','k--','linewidth',2); PlotHVLines(0,'v','k--','linewidth',2)
subplot(2,3,4); ylim(y); PlotHVLines(0,'h','k--','linewidth',2); PlotHVLines(0,'v','k--','linewidth',2)


subplot(2,3,6);
semplot(ht,zDiff(indexCA1,:),'k')
semplot(ht,zDiff(indexCA2,:),'b')
ylim(y); PlotHVLines(0,'h','k--','linewidth',2); PlotHVLines(0,'v','k--','linewidth',2)

SaveFig(['G:\My Drive\Doctorado\Estancia Cornell\Lab\Proyecto Mental Imagery\Labmeeting\figures\ripple_PETHs_' sessionID ])

%%

clf
bins = Bins(0,sleep(end),60*10,60*10); rippleRate = CountInIntervals(r(:,1),bins)/diff(bins(1,:));
bar(mean(bins,2),rippleRate);
handle = PlotIntervals(SubtractIntervals([0 Inf],sleep),'color','r');
legend(handle,'maze','Location','northwest','box','off','fontsize',15);
set(gca,'box','off','tickdir','out','fontsize',12);
xlabel('time (s)')
ylabel('ripple rate (per second)')
SaveFig(['G:\My Drive\Doctorado\Estancia Cornell\Lab\Proyecto Mental Imagery\Labmeeting\figures\rippleRate_' sessionID ])




