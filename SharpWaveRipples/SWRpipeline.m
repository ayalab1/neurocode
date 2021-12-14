%% Ripple analysis pipeline
% This is a wrapper to concatenate all functions to perform basic ripple analysis in a single session
basepath = pwd;
basename = basenameFromBasepath(basepath);

%% 0- Pre processing
%   - Initial preprocessing and detection of sleep states is done in preprocessSession
%   - Spikes and cell_metrics are extraction in preprocessSpikes (after manual clustering)

%% 1 - Automatic estimation of best channels for swr detection
   swrCh = swrChannels('basepath',basepath);
   % this code often doesn't work well if you have dentate gyrus channels 
   
%% 2 - SWR detection

    % preferred method
     ripples = DetectSWR([swrCh.ripple swrCh.sharpwave],'saveMat',true);
     % only use this is you don't have sharp-wave
     ripples = FindRipples(basepath,swrCh.ripple,'noise',swrCh.noise,'saveMat',true);
     
% optional steps     
     
    % refine ripple detection using spiking level
    spikes = importSpikes('cellType', "Pyramidal Cell", 'brainRegion', "CA1");
    ripples = eventSpikingTreshold(ripples,'spikes',spikes,'spikingThreshold',0.5); 
    save([basename '.ripples.events.mat'],'ripples');saveas(gcf,'SWRmua.png');
    
    % remove very large amplitude events (likely artifacts)
    ripples = removeArtifactsFromEvents(ripples,'stdThreshold',10);
    
    % Add ripples to cell metrics. ProcessCellMetrics will locate newly created .events file
    load(fullfile(basepath,[basename,'.session.mat']))
    cell_metrics = ProcessCellMetrics('session',session,'manualAdjustMonoSyn',false);

    % restrict ripples to certain intervals (e.g. NREM sleep)
    load([basename '.SleepState.states.mat']);
    ripples = eventIntervals(ripples,SleepState.ints.NREMstate,1);
    
    % save a .evt file to inspect in Neuroscope 
    createEVT(ripples);

    % plot wavelet to check detection quality 
    mkdir('Ripple_Profile');
    lfpRip = getLFP(swrCh.ripple);
    [wavAvg,lfpAvg] = eventWavelet(lfpRip,ripples.peaks(1:500),'twin',[0.1 0.1]);
    saveas(gcf,['Ripple_Profile\swrWaveletSample.png']);
    
%% 3- Separate ripples by task epochs     
    % creat pre/task/post structure. Need to be improved
    [behavEpochs.int_samples,behavEpochs.int] = make_pre_task_post(pwd);
    save([basename '.behavEpochs.mat'],'behavEpochs')
 
    % create swr structures for each epoch
    inputLabels={'behavEpochs.int.pre','behavEpochs.int.task','behavEpochs.int.post'};
    outputLabels={'SWRepochs.pre','SWRepochs.task','SWRepochs.post'};

    SWRepochs=[];
for epochs= 1:3 
            eval(['stateInterval=' inputLabels{epochs} ]);
            if ~isempty(stateInterval)
            clear var tmp1 St;
            tmp1 = find(ripples.peaks>=stateInterval(1) & ripples.peaks<=stateInterval(2));
            St.timestamps = ripples.timestamps(tmp1,:);
            St.peaks = ripples.peaks(tmp1,:);
            St.stdev = ripples.stdev;
            try
            St.SvMax = ripples.SwMax(tmp1,:);
            St.RipMax = ripples.RipMax(tmp1,:);
            catch; end
            St.detectorinfo = ripples.detectorinfo; 
            for i = 1:length(St.peaks)
               St.duration(i,1) = St.timestamps(i,2)-St.timestamps(i,1);
            end
            eval([outputLabels{epochs} '=St;' ])
            else
            eval([outputLabels{epochs} '=[];' ])   
            end
end
    save([basename '.SWRepochs.mat'],'SWRepochs');   
 
    %optional. plot wavelet spectrogram for each epoch
    figure; names={'PRE','TASK','POST'}; clear wavAvg lfpAvg;
for epochs= 1:3
    if ~isempty(eval([outputLabels{epochs}]))
    peaks = eval([outputLabels{epochs} '.peaks;' ]);
    if numel(peaks) > 100
        peaks = peaks(1:100);
    end
    [wavT,lfpT]= eventWavelet(lfpRip,peaks,'twin',[0.1 0.1],'plotWave',false,'plotLFP',false);
    wavAvg{epochs} = wavT; lfpAvg{epochs} = lfpT; clear lfpT wavT;
    
    subplot(1,3,epochs);
    contourf(wavAvg{epochs}.timestamps*1000,wavAvg{epochs}.freqs,wavAvg{epochs}.data',30,'LineColor','none');hold on;
    set(gca,'YScale','log');
    ylim([wavAvg{epochs}.freqs(1) wavAvg{epochs}.freqs(end)]);
    colormap jet;
    xlabel('time (ms)'); ylabel('frequency (Hz)');
   
    fmin =  wavAvg{epochs}.freqs(round(length(wavAvg{epochs}.freqs)/4)); 
    fmax =  fmin*2;
    lfpSc =(lfpAvg{epochs}.data-min(lfpAvg{epochs}.data))*(fmax-fmin)/(max(lfpAvg{epochs}.data)-min(lfpAvg{epochs}.data))+fmin;
    plot(wavAvg{epochs}.timestamps*1000,lfpSc,'w','LineWidth',2);hold on;
    title(names{epochs});
    end
end   
    saveas(gcf,['swrWaveletSample.png']);

%% 4-Compare SWR properties (duration, freq, power, rate) in each task block 
   colorR = {'b','k','r'};
figure;
for fig = 1
subplot(2,3,1);
for epochs= 1:3
    if ~isempty(eval([outputLabels{epochs}]))
    tt = eval([outputLabels{epochs} '.timestamps(:,1);' ]);
    durtt=tt(end)-tt(1);
    rate(epochs)=numel(tt)/durtt;
    end
end
bar(rate);set(gca,'XTickLabel',{'PRE','TASK','POST'});ylabel('swr rate (Hz)');title('SWR rate');
subplot(2,3,2);
for epochs= 1:3
    if ~isempty(eval([outputLabels{epochs}]))
    dur = eval([outputLabels{epochs} '.duration;' ])*1000;
    h=histc(dur,20:5:200);
    plot([20:5:200],h/sum(h),'color',colorR{epochs},'LineWidth',2); hold on;
    set(gca,'XScale','log');xlabel('duration (ms)');ylabel('frac. of swr');
    end
end
title('duration');
subplot(2,3,3);
for epochs= 1:3
    if ~isempty(eval([outputLabels{epochs}]))
    dur = eval([outputLabels{epochs} '.duration;' ])*1000;
    long(epochs)=sum(dur>100)/numel(dur);
    end
end
bar(long);set(gca,'XTickLabel',{'PRE','TASK','POST'});ylabel('frac. swr > 100ms');
subplot(2,3,4);
for epochs= 1:3
    if ~isempty(eval([outputLabels{epochs}]))
    SvMax = eval([outputLabels{epochs} '.SvMax;' ]);
    h=histc(SvMax,2:0.3:10);
    plot([2:0.3:10],h/sum(h),'color',colorR{epochs},'LineWidth',2); hold on;
    set(gca,'XScale','log');xlabel('SW pow');ylabel('frac. of swr');
    end
end
title('SW pow');
subplot(2,3,5);
for epochs= 1:3
    if ~isempty(eval([outputLabels{epochs}]))
    RipMax = eval([outputLabels{epochs} '.RipMax;' ]);
    h=histc(RipMax,2:0.4:15);
    plot([2:0.4:15],h/sum(h),'color',colorR{epochs},'LineWidth',2); hold on;
    set(gca,'XScale','log');xlabel('rip pow');ylabel('frac. of swr');
    end
end
title('rip pow');
subplot(2,3,6);
plot(0,'color',colorR{1});hold on;plot(0,'color',colorR{2});hold on;plot(0,'color',colorR{3});hold on;
legend({'PRE','TASK','POST'});
end
saveas(gcf,['swrMetrics.png']);

%% 5- SWR psth pre/post
win = [-0.2 0.2];

        st1 = SWRepochs.pre.peaks;
        spikeResponse1 = [];
        for jj = 1:size(spikes.UID,2)
            [stccg, t] = CCG({spikes.times{jj} st1},[],'binSize',0.005,'duration',1);
            spikeResponse1 = [spikeResponse1; zscore(squeeze(stccg(:,end,1:end-1)))'];
        end
        st2 = SWRepochs.post.peaks;
        spikeResponse2 = [];
        for jj = 1:size(spikes.UID,2)
            [stccg, t] = CCG({spikes.times{jj} st2},[],'binSize',0.005,'duration',1);
            spikeResponse2 = [spikeResponse2; zscore(squeeze(stccg(:,end,1:end-1)))'];
        end      
        for s = 1:size(spikeResponse1,1)
            spikeResponse3(s,:) = smooth(spikeResponse2(s,:),3,'sgolay')-smooth(spikeResponse1(s,:),3,'sgolay');
            maxDif(s,1) = max(spikeResponse3(s,80:120));
        end
        [b,ind] = sort(maxDif,'descend');
        
        figure;
        subplot(1,4,1);
        imagesc([t(1) t(end)],[1 size(spikeResponse1,1)], spikeResponse1(ind,:)); caxis([-3 3]); colormap(jet);
        xlim([-.2 .2]); set(gca,'TickDir','out'); xlabel('Time'); ylabel('Cells');title('PRE');
        subplot(1,4,2)
        imagesc([t(1) t(end)],[1 size(spikeResponse2,1)], spikeResponse2(ind,:)); caxis([-3 3]); colormap(jet);
        xlim([-.2 .2]); set(gca,'TickDir','out'); xlabel('Time'); ylabel('Cells');title('POST');      
        subplot(1,4,3);
        imagesc([t(1) t(end)],[1 size(spikeResponse3,1)], spikeResponse3(ind,:)); caxis([-3 3]); colormap(jet);
        xlim([-.2 .2]); set(gca,'TickDir','out'); xlabel('Time'); ylabel('Cells');title('DIF');
        subplot(1,4,4);
        plot(smooth(mean(spikeResponse1),5,'sgolay'),'k','LineWidth',2);hold on;
        plot(smooth(mean(spikeResponse2),5,'sgolay'),'b','LineWidth',2);hold on;
        plot(smooth(mean(spikeResponse3),7,'sgolay'),'r','LineWidth',2);hold on;
        xlim([50 150]); 

%% 6 - Calculate SWR metrics for each cell
inputLabels={'SWRepochs.pre','SWRepochs.task','SWRepochs.post'};
outputLabels={'SWRspikes.pre','SWRspikes.task','SWRspikes.post'};
load([basename '.spikes.cellinfo.mat']);
SWRspikes=[];

    % get spikes within swr
    for states = 1:3
        clear var stateMatrix M
        eval( ['stateMatrix=' inputLabels{states} ';']);
       if ~isempty(stateMatrix)
          M = getRipSpikes('spikes',spikes,'events',stateMatrix  ,'saveMat',false); 
          eval([outputLabels{states} '=M;']);
       else
          eval([outputLabels{states} '=[];']); 
       end
    end
    save([basename '.SWRspikes.mat'],'SWRspikes');  

    % Calculate Cell ripple firing properties per task epoch
    inputLabels={'SWRspikes.pre','SWRspikes.task','SWRspikes.post'};
    outputLabels={'SWRunitMetrics.pre','SWRunitMetrics.task','SWRunitMetrics.post'};
    SWRunitMetrics=[];

    for states=1:3
        clear var stateMatrix M
        if ~isempty(eval([inputLabels{states}]))
            eval(['stateMatrix=' inputLabels{states} ';']);
            M = unitSWRmetrics(stateMatrix);
            eval([outputLabels{states} '=M;']);
        else
            eval([outputLabels{states} '=[];']);    
        end
    end
    save([basename '.SWRunitMetrics.mat'],'SWRunitMetrics');  

    % add some plots here...
    



