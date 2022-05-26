basepath = pwd;
[parentFolder,basename] = fileparts(basepath);
[~,projectName] = fileparts(parentFolder);
sessionID = [projectName '-' basename];

% eegFile = fullfile(basepath,[basename '.eeg']); % This is where I will store the normalized LFP file
% if ~exist(eegFile,'file')
%     lfpFile = fullfile(basepath,[basename '.lfp']);
%     copyfile(lfpFile,eegFile);
%     file = memmapfile(eegFile,'Format','int16','Writable',true);
%     data = reshape(file.Data,64,[]);
%     m = int16(mean(data(1+[24 11 20 7 9 21 8 25 12 10 23 14 27 26 5 18 28 3 31 62 0 1 30 16 2 29 6 4 36 59 44 34 61 45 33 17 63 32 46 60 35 56 58 40 38 53 42],:)));
%     newData = bsxfun(@minus,data,m);
%     file.data = newData(:);
%     clear file
% end

ripples = DetectSWR(1+[54 35],'saveMat',true,'check',true,'useSPW',false,'useEEG',true);
% ripples = FindRipples0(1+[54 35],'saveMat',true,'check',true,'useSPW',false,'useEEG',true);

% lfp = GetAyaEEG(54);
% [clean,~,badIntervals] = CleanLFP(lfp,'thresholds',[2 2],'manual',true);
% goodIntervals = SubtractIntervals([0 lfp(end,1)],badIntervals);
% ripples = FindRipples(basepath,1+54,'noise',1+7,'saveMat',true,'restrict',goodIntervals);
% figure; plot(tl,(zscore(double(lfp(:,2))))); Browse('all');

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

write = false; % don't save the figures
%%

nBins = 101;
m = nan(nNeurons,nBins); mPre = nan(nNeurons,nBins); mPost = nan(nNeurons,nBins);

for i=1:nNeurons
    [h,ht] = PETH(spikes(spikes(:,2)==i),ripples.timestamps(:,1),'nBins',nBins);
    m(i,:) = mean(h);
    mPre(i,:) = mean(h(pre,:));
    mPost(i,:) = mean(h(post,:));
end

for kind=1:2,
    if kind==1 % zscore together; this also captures changes in the global firing rate
        z = [mPre mPost];
        zPre = (mPre - mean(z,2))./std(z,[],2);
        zPost = (mPost - mean(z,2))./std(z,[],2);
    else
        zPre = zscore(mPre,[],2);
        zPost = zscore(mPost,[],2);
    end
    zDiff = zPost - zPre;

    clf
    subplot(2,3,1);
    PlotColorMap(zPre,'x',ht);
    set(gca,'TickDir','out');
    ylabel(['Neuron ID: ' regionListString]); xlabel('time from ripple start (s)');
    PlotHVLines(0,'v','k--','linewidth',2);
    PlotHVLines(cumsum(Accumulate(regionID))+0.5,'h','w--','linewidth',2);
    clabel('firing rate (z units)');
    title('Pre sleep ripples')

    subplot(2,3,2);
    PlotColorMap(zPost,'x',ht);
    set(gca,'TickDir','out');
    ylabel(['Neuron ID: ' regionListString]); xlabel('time from ripple start (s)');
    PlotHVLines(cumsum(Accumulate(regionID))+0.5,'h','w--','linewidth',2);
    PlotHVLines(0,'v','k--','linewidth',2);
    clabel('firing rate (z units)')
    title('Post sleep ripples')
    clims

    subplot(2,3,3);
    PlotColorMap(zDiff,'x',ht);
    set(gca,'TickDir','out');
    ylabel(['Neuron ID: ' regionListString]); xlabel('time from ripple start (s)');
    PlotHVLines(cumsum(Accumulate(regionID))+0.5,'h','w--','linewidth',2);
    PlotHVLines(0,'v','k--','linewidth',2);
    clabel('firing rate difference (z units)')
    ColorMap(gca,[0 0 1],[1 1 1],[1 0 0]);
    title('Post-Pre')
    clim(max(abs(clim))*[-1 1]);

    subplot(2,3,4);
    colors = get(gca,'colororder'); clear handles
    for i=1:max(regionID),
        if sum(regionID==i)>1
            handles(i) = semplot(ht,zPre(regionID==i,:),colors(i,:));
        else
            handles(i) = plot(ht,zPre(regionID==i,:),'color',colors(i,:),'linewidth',2);
        end
        if i==1, hold on; end
    end
    handles4 = handles;
    y = ylim;

    subplot(2,3,5);
    colors = get(gca,'colororder'); clear handles
    for i=1:max(regionID),
        if sum(regionID==i)>1
            handles(i) = semplot(ht,zPost(regionID==i,:),colors(i,:));
        else
            handles(i) = plot(ht,zPost(regionID==i,:),'color',colors(i,:),'linewidth',2);
        end
        if i==1, hold on; end
    end
    y = ylim;

    y = [min(y(1),min(ylim)),max(y(2),max(ylim))];
    ylim(y); PlotHVLines(0,'h','k--','linewidth',2); PlotHVLines(0,'v','k--','linewidth',2); legend(handles,regionNames{:});  legend('box','off','location','southwest');
    subplot(2,3,4); ylim(y); PlotHVLines(0,'h','k--','linewidth',2); PlotHVLines(0,'v','k--','linewidth',2); legend(handles4,regionNames{:}); legend('box','off','location','southwest');

    subplot(2,3,6);
    colors = get(gca,'colororder'); clear handles
    for i=1:max(regionID),
        if sum(regionID==i)>1
            handles(i) = semplot(ht,zDiff(regionID==i,:),colors(i,:));
        else
            handles(i) = plot(ht,zDiff(regionID==i,:),'color',colors(i,:),'linewidth',2);
        end
        if i==1, hold on; end
    end
    ylim(y); PlotHVLines(0,'h','k--','linewidth',2); PlotHVLines(0,'v','k--','linewidth',2);
    legend(handles,regionNames{:});
    legend('box','off','location','southwest');

    drawnow
    if write
        if kind==1,
            SaveFig(['M:\home\raly\results\V1\ripple_PETHs_' sessionID])
        else
            SaveFig(['M:\home\raly\results\V1\ripple_PETHs_zscored_individually_' sessionID])
        end
    end
end
% SaveFig(['G:\My Drive\Doctorado\Estancia Cornell\Lab\Proyecto Mental Imagery\Labmeeting\figures\ripple_PETHs_' sessionID ])
%%

clf
bins = Bins(0,MergePoints.timestamps(end),60*10,60*10); rippleRate = CountInIntervals(r(:,1),bins)/diff(bins(1,:));
bar(mean(bins,2),rippleRate);
handle = PlotIntervals(SubtractIntervals([0 MergePoints.timestamps(end)],sleep),'color','r');
legend(handle,'maze','Location','northwest','box','off','fontsize',15);
set(gca,'box','off','tickdir','out','fontsize',12);
xlabel('time (s)')
ylabel('ripple rate (per second)')
% SaveFig(['G:\My Drive\Doctorado\Estancia Cornell\Lab\Proyecto Mental Imagery\Labmeeting\figures\rippleRate_' sessionID ])

SaveFig(['M:\home\raly\results\V1\rippleRate_' sessionID])
%%
for i=1:6
    subplot(2,3,i);
    [h,ht] = PETH(spikes(spikes(:,2)==i,1),r(:,1),'durations',[-1 1]*0.5);
    PlotColorMap(Shrink(sortby(zscore(h,[],2),ripples.RipMax),1000,1),'x',ht);
    PlotHVLines(0,'v','k','linewidth',2);
end

