
%% Ripple analysis pipeline
% This is a wrapper to concatenate all functions to perform basic ripple analysis in a single session
% run it from session folder
basename = basenameFromBasepath(pwd);

%% 0- Pre processing
%   - Detection of sleep states and theta epochs is done in bz_PreprocessSession
%   - Spikes and cell_metrics are calculated in bz_PreprocessSpikes 

%% 1 - Atomatic best channel detection
   swrCh = bz_swrChannels('basepath',basepath,'Manual',true);
   % need to fix 0-1 channel indx,remove promts, save figure and varible
   % by default
   % Farnaz:  promts are remove ,  figure and varible  are fixed by by default 
   % and channel indexing is 0 now
%% 2 - SWR detection

if isnan(swrCh.sharpwave)==0
     ripples = bz_DetectSWR([swrCh.ripple swrCh.sharpwave],'saveMat',true);
elseif isnan(Noise)==0
    [ripples] = bz_FindRipples(basepath,swrCh.ripple,'noise',swrCh.noise,'saveMat',true);
else
    [ripples] = bz_FindRipples(basepath,swrCh.ripple,'saveMat',true);
end

%%
    % refine ripple detection using spiking level
    ripples = eventSpikingTreshold(ripples,'spikes',spikes,'spikingThreshold',0.5); 

    % plot wavelet to check detection quality 
    mkdir('Ripple_Profile');
    lfpRip = bz_GetLFP(swrCh.ripple);
    [wavAvg,lfpAvg]=bz_eventWavelet(lfpRip,ripples.peaks(1:500),'twin',[0.1 0.1]);
    saveas(gcf,['Ripple_Profile\ripWaveletSample.png']);
    
%% 3- Separate task epochs     
% creat pre/task/post 
    % needs to take into account different session structure
        % maybe with an option to manually provide the numbers of preSleep.
        % task and postSleep sessions to make_pre_task_post
    [behavEpochs.int_samples,behavEpochs.int] = make_pre_task_post(pwd);
    save([basename '.behavEpochs.mat'],'behavEpochs')
    
% plot behaviour (optional)
    clearvars;
    basename = basenameFromBasepath(pwd);
    load([basename '.behavEpochs.mat']);
    thetaEpochs(pwd);
    
    Plot_recording_States(pwd);  
    saveas(gcf,[basename '.states.png']);
    
% Add theta run states
    load([basename '.SleepState.states.mat'])
    Int_run=behavEpochs.int.task;  
    Thetastate = SleepState.ints.THETA;  
    ND=find(Thetastate(:,1)>=Int_run(1) & Thetastate(:,2)<=Int_run(2));
    Thetastate_Run=Thetastate(ND,:);
    SleepState.ints.THETAtask=Thetastate_Run;
    save([basename '.SleepState.states.mat'],'-append','SleepState')
    close all;
    
%% 4- Separate ripples by task epochs
clearvars;
basename = basenameFromBasepath(pwd);
load([basename '.behavEpochs.mat']);
load([basename '.ripples.events.mat']);
load([basename '.swrCh.mat']);
lfpRip = getLFP(swrCh.ripple);

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
            St.SvMax = ripples.SwMax(tmp1,:);
            St.RipMax = ripples.RipMax(tmp1,:);
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
    saveas(gcf,['Ripple_Profile\ripWaveletSample.png']);

% 5 - SWR properties per task epoch
%clearvars;
%basename = bz_BasenameFromBasepath(pwd);

   % Compared SWR properties (duration, freq, power, rate) in each task block 
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
saveas(gcf,['Ripple_Profile\ripMetrics.png']);

% 6 - Get spikes in ripples
inputLabels={'SWRepochs.pre','SWRepochs.task','SWRepochs.post'};
outputLabels={'SWRspikes.pre','SWRspikes.task','SWRspikes.post'};
load([basename '.spikes.cellinfo.mat']);
SWRspikes=[];

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

% 7 - Cell ripple firing properties per task epoch
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

%% 8 - Rank order 
% it will be good to save the rank as another matrix inside each SWRunitMetrics struct
inputLabels={'SWRspikes.pre','SWRspikes.task','SWRspikes.post'};
outputLabels={'rank_m.pre','rank_m.task','rank_m.post'};
rank_m=[];

for states=1:3
    clear var stateMatrix M
    eval( ['stateMatrix=' inputLabels{states} ';']);
    
    if isempty(stateMatrix)==1
       M=[];
    else 
    [M] = bz_RankOrder('basepath',pwd,'spkEventTimes',stateMatrix,'timeSpike','mean','minUnits',5, 'doPlot',false,'saveMat',false,'numRep',0);
    end
    eval([outputLabels{states} '=M;']);
end
    
%% 9 - Calculate place fields
   % 1-Calculate firing mals L-R and R-L
    load('posTrials.mat'); load('day20.spikes.cellinfo.mat')
    [firingMaps] = bz_firingMapAvg(posTrials,spikes);
    save([basename '.firingMapsAvg.cellinfo.mat'],'firingMaps');
   % 2- Detect (average) place fields
    [placeFieldStats] = bz_findPlaceFields1D('firingMaps',firingMaps,'minPeak',1.5,'sepEdge',0.04,'doPlot',false,'saveMat',false);
    save([basename '.placeFields.cellinfo.mat'],'placeFieldStats');
   % 3- Calculate template of place fields in the maze
    [placeFieldTemplate] = bz_findPlaceFieldsTemplate('placeFieldStats',placeFieldStats,'firingMaps',firingMaps,'saveMat',false);
    save([basename '.placeFieldTemplate.mat'],'placeFieldTemplate');
     
%% 10 - Replay (using rank order correlation)
   % 1 - Select only place cells - this is not workign well 
    %    (not sure this is necessary to do at all)
%      PCuid = unique([placeFieldTemplate.Peak{1}(:,3);placeFieldTemplate.Peak{2}(:,3)]);
%      for i = 1:length(placeFieldStats.UID)
%          if ~isempty(find(PCuid == placeFieldStats.UID(i)))
%              goodPCs(i,1) = 1;
%          else
%              goodPCs(i,1) = 0;
%          end
%      end
%      
%      % generate new SWRspikes only with PCs
%      inputLabels={'SWRepochs.pre','SWRepochs.task','SWRepochs.post'};
%      outputLabels={'SWRspikesPCs.pre','SWRspikesPCs.task','SWRspikesPCs.post'};
%      [spikesPC] = bz_ImportSpikes('UID',PCuid);
%      SWRspikesPCs=[];
%      for states = 1:3
%         clear var stateMatrix M
%         eval( ['stateMatrix=' inputLabels{states} ';']);
%         M = bz_getRipSpikes('spikes',spikesPC,'events',stateMatrix  ,'saveMat',false); 
%         eval([outputLabels{states} '=M;']);
%      end
%      save([basename '.SWRspikesPCs.mat'],'SWRspikesPCs');     
%     
%      % generate new placeFieldTemplate to match SWRspikesPC
%      [TfiringMaps] = bz_firingMapAvg(posTrials,spikesPC,'saveMat',false);
%      [TplaceFieldStats] = bz_findPlaceFields1D('firingMaps',TfiringMaps,'minPeak',2,'sepEdge',0.04,'saveMat',false,'doPlot',false);
%      [placeFieldTemplate] = bz_findPlaceFieldsTemplate('placeFieldStats',TplaceFieldStats,'firingMaps',TfiringMaps);                
        

   % 2- Rank order correlation using maze template
       % method from Diba & Buzsaki,NatNeur, 2006
    inputLabels={'SWRspikes.pre','SWRspikes.task','SWRspikes.post'};
    outputLabels={'rankStats.pre','rankStats.task','rankStats.post'};rankStats=[];
    for states=1:3
        clear var stateMatrix M
        eval( ['stateMatrix=' inputLabels{states} ';']);
        if isempty(stateMatrix)==1
           M=[];
        else 
           [M] = bz_RankOrder('spkEventTimes',stateMatrix,'templateType','Peak',...
                      'timeSpike','first','minUnits',5,'numRep',500);    end
        eval([outputLabels{states} '=M;']);
    end   

   % 3- Identify significant forward and reverse replay events
      % don't understand why no reverse replay events
        for m = 1:2
           fowReplayInd.pre{m} = find(rankStats.pre.pvalEvents(m,:)<0.05 & rankStats.pre.corrEvents(m,:)>0);
           revReplayInd.pre{m} = find(rankStats.pre.pvalEvents(m,:)<0.05 & rankStats.pre.corrEvents(m,:)<0);
           fowReplayInd.task{m} = find(rankStats.task.pvalEvents(m,:)<0.05 & rankStats.task.corrEvents(m,:)>0);
           revReplayInd.task{m} = find(rankStats.task.pvalEvents(m,:)<0.05 & rankStats.task.corrEvents(m,:)<0);
           fowReplayInd.post{m} = find(rankStats.post.pvalEvents(m,:)<0.05 & rankStats.post.corrEvents(m,:)>0);
           revReplayInd.post{m} = find(rankStats.post.pvalEvents(m,:)<0.05 & rankStats.post.corrEvents(m,:)<0);
        end
        
        % distribution of crr values for significant events
        figure; 
        subplot(1,2,1);
        h1=histc(cat(2,rankStats.pre.corrEvents(1,fowReplayInd.pre{1}),rankStats.pre.corrEvents(2,fowReplayInd.pre{2})),[0:0.1:1]);
        plot([0:0.1:1],h1/length(rankStats.pre.corrEvents)*100,'b','LineWidth',2);hold on;
        h2=histc(cat(2,rankStats.task.corrEvents(1,fowReplayInd.task{1}),rankStats.task.corrEvents(2,fowReplayInd.task{2})),[0:0.1:1]);
        plot([0:0.1:1],h2/length(rankStats.task.corrEvents)*100,'k','LineWidth',2);hold on;
        h3=histc(cat(2,rankStats.post.corrEvents(1,fowReplayInd.post{1}),rankStats.post.corrEvents(2,fowReplayInd.post{2})),[0:0.1:1]);
        plot([0:0.1:1],h3/length(rankStats.post.corrEvents)*100,'r','LineWidth',2);hold on;    
        subplot(1,2,2);
        bar((numel(fowReplayInd.task{1})+numel(fowReplayInd.task{2}))/length(rankStats.task.corrEvents));        
        
        % plot individual events
        for c = 1:2
            for i = 1:length(fowReplayInd.task{c})
                unitEvent{c}{i} = SWRspikes.task.EventRel{fowReplayInd.task{c}(i)}(2,:);
                for u = 1:length(unitEvent{c}{i})
                    indT=find(placeFieldTemplate.Peak{c}(:,2)==unitEvent{c}{i}(u));
                    if ~isempty(indT)
                        unitEvent{c}{i}(2,u) = placeFieldTemplate.Peak{c}(indT,1);
                    else
                        unitEvent{c}{i}(2,u) = NaN;
                    end
                end
            end
        end
        
        
%%   

        fowReplayInd{m}=[]; revReplayInd{m}=[];
        temp = find(rankStats.pvalEvents(m,:)<0.05);
        for i = 1:length(temp)
            if rankStats.corrEvents(m,temp(i))>0
               fowReplayInd{m} = cat(1,fowReplayInd{m},temp(i));
            elseif rankStats.corrEvents(m,temp(i))<0
               revReplayInd{m} = cat(1,revReplayInd{m},temp(i));
            end
        end






