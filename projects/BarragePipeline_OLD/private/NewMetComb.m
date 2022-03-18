%% Cumulative Metrics for Figures 

function NewMetComb(logicIn,combine)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NewMetComb
%
% This script is for calculating population metrics across many animals and
% many sessions.
%
% INPUTS
% logicIn -> array containing binary values denoting:
%   (1) NEWRIPS --> whether or not new ripples have been detected since 
%   the last run (default 1)
%   (2) NEWEVTS --> whether or not new barrages have been detected since 
%   the last run (default 1)
%   (3) skipPSTH --> whether or not we should calculate a population PSTH
%   (default 0, which DOES plot the PSTH)
% combine -> string denoting if we should be plotting mouse data ("mouse"),
% rat data ("rat"), or both ("both")
%
% Lindsay K, 3/22
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if nargin < 2
    combine = "both";
elseif nargin < 1
    NEWRIPS = 1;
    NEWEVTS = 1;
    skipPSTH = 0;
    combine = "both";
else
    NEWRIPS = logicIn(1);
    NEWEVTS = logicIn(2);
    skipPSTH = logicIn(3);
end
if ~((combine == "mouse")||(combine == "rat")||(combine == "both"))||(nargin < 1)
    warning('Must choose mouse, rat, or both. Defaulting to both');
    combine = "both";
end
savePath = ('Z:\home\Lindsay\Barrage\');
if combine == "both"
    useSess = readtable('Z:\home\Lindsay\Barrage\rat_sessions.csv');
    useSess_temp2 = readtable('Z:\home\Lindsay\Barrage\mouse_sessions.csv');
    useSess(size(useSess,1)+1:size(useSess,1)+size(useSess_temp2,1),:) = useSess_temp2; clear useSess_temp2
elseif combine == "mouse"
    useSess = readtable('Z:\home\Lindsay\Barrage\mouse_sessions.csv');
else
    useSess = readtable('Z:\home\Lindsay\Barrage\rat_sessions.csv');
end

%% Initialize cells and arrays
evtDur = []; %event durations
evtDur_tot = [];
evtDurRip_tot = [];
rateB = []; rateR = [];
absFRpx = []; %absolute firing rate (during bar) data pyramidal cells
absFRpg = []; %absolute firing rate (during bar) group pyramidal cells
absFRix = []; %absolute firing rate (during bar) data interneurons
absFRig = []; %absolute firing rate (during bar) group interneurons
absFR2x = []; %CA2 modulation
absFR2g = []; %CA2 modulation
catPSTH = cell(7,1);
catPSTH_rip = cell(7,1); %no unknown
catPSTHsc = cell(7,1);
catPSTHsc_rip = cell(7,1); %no unknown
gainFRpb = cell(8,1);
gainFRpbTYPE = cell(8,1);
gainFRps = cell(8,1);
gainMetp = cell(8,1);
gainMetpg = cell(8,1);
gainMetRip = cell(8,1);
gainFRib = cell(8,1);
gainFRis = cell(8,1);
% gainMeti = cell(8,1);
burstcp = cell(8,1);

%% Iterate through saved sessions
for p = 1:size(useSess,1) 
% for p = 1:2
%     disp(paths_save(p)); %just need some form of progress bar
%     fprintf('This message is sent at time %s\n', datestr(now,'HH:MM:SS'));
    ses = useSess.basepath{p};
    basepath = useSess.basepath{p};
    basename = useSess.basename{p};
    animName = useSess.animal{p};
    cd(ses);
    
    if contains(basepath, 'Z:\')
        label = "rat";
    elseif contains(basepath, 'Y:\')
        label = "mouse";
    else
        error('Check location labels');
    end
    if (label==combine)||(combine=="both")
%         load(strcat(basepath,'\Barrage_Files\Population\',basename,'.useSpk.props.mat'));
        load(strcat(basepath,'\',basename,'.cell_metrics.cellinfo.mat'));
        load(strcat(basepath,'\Barrage_Files\Cell\',basename,'.props.mat'),'regID','burstIndex','regKey');
        load(strcat(basepath,'\Barrage_Files\Population\',basename,'.useSpk.spkEventTimes.mat')); barSpkEventTimes = spkEventTimes; clear spkEventTimes
        load(strcat(basepath,'\Barrage_Files\Population\',basename,'.useSpk.unitBar.mat'));
        load(strcat(basepath,'\',basename,'.session.mat'));
        load([basepath '\' basename '.SleepState.states.mat']);
        load(strcat(basepath,'\',basename,'.spikes.cellinfo.mat')); spikesALL = spikes; clear spikes
        load(strcat(basepath,'\Barrage_Files\',basename,'.allpyr.cellinfo.mat')); spikesPYR = spikes; clear spikes
        load(strcat(basepath,'\Barrage_Files\',basename,'.allint.cellinfo.mat')); spikesINT = spikes; clear spikes
        load(strcat(basepath,'\',basename,'.ripples.events.mat'));
%         load(strcat(basepath,'\',basename,'.spkEventTimes.mat')); ripSpkEventTimes = spkEventTimes; clear spkEventTimes      
        intSave = strcat(basepath, '\Barrage_Files\Population\',basename,'.useSpk');
        %% Rerun mets for interneurons
        if (~exist(strcat(intSave,'.barSpkEventTimesINT.mat')))||(NEWEVTS)
            load([basepath '\Barrage_Files\' basename '.HSE.mat']);
            if isfield(HSE, 'keep')
                HSE.timestamps = HSE.timestamps(HSE.keep,:);
                HSE.peaks = HSE.peaks(HSE.keep);
                HSE.amplitudes = HSE.amplitudes(HSE.keep);
                HSE.duration = HSE.duration(HSE.keep);
            end
            HSEnrem = eventIntervals(HSE,SleepState.ints.NREMstate,1);
            [barSpkEventTimesINT] = getRipSpikes('basepath',pwd,'events',HSEnrem.timestamps,'spikes',spikesINT,'padding',0,'saveMat',false);
            save(strcat(intSave,'.barSpkEventTimesINT.mat'),'barSpkEventTimesINT', '-v7.3');
            [unitBarINT] = unitSWRmetrics(barSpkEventTimesINT,spikesINT);
            save(strcat(intSave,'.unitBarINT.mat'),'unitBarINT', '-v7.3');
        else
            load(strcat(intSave,'.barSpkEventTimesINT.mat'));
            load(strcat(intSave,'.unitBarINT.mat'));
        end
        
        
        %% Rerun mets for ripples?
        if ~exist(strcat(basepath,'\Barrage_Files\',basename,'.ripSpkEventTimes.mat'))||(NEWRIPS)
            ripSpkEventTimes = getRipSpikes('basepath', pwd, 'events', ripples.timestamps,'spikes',spikesPYR,'padding',0,'saveMat',false);
            save(strcat(basepath,'\Barrage_Files\',basename,'.ripSpkEventTimes.mat'),'ripSpkEventTimes','-v7.3');
        else
            load(strcat(basepath,'\Barrage_Files\',basename,'.ripSpkEventTimes.mat'));
        end
        if ~exist(strcat(basepath,'\Barrage_Files\',basename,'.unitRip.mat'))||(NEWRIPS)
            unitRip = unitSWRmetrics(ripSpkEventTimes,spikesPYR);
            save(strcat(basepath,'\Barrage_Files\',basename,'.unitRip.mat'),'unitRip','-v7.3');
        else
            load(strcat(basepath,'\Barrage_Files\',basename,'.unitRip.mat'));
        end

        %% Population Level Event Duration Histogram
        binsEvtDur = [0:0.01:2];
        evtDur(p,1:length(binsEvtDur)) = histc(barSpkEventTimes.EventDuration, binsEvtDur);
        evtDur(p,1:length(binsEvtDur)) = evtDur(p,1:length(binsEvtDur))/sum(evtDur(p,1:length(binsEvtDur)));
        evtDur_tot = [evtDur_tot; barSpkEventTimes.EventDuration];
        try
            sesDur = (session.epochs{1,end}.stopTime-session.epochs{1,1}.startTime);
        catch
            sesDur = session.general.duration;
        end
        
        evtDurRip(p,1:length(binsEvtDur)) = histc(ripSpkEventTimes.EventDuration, binsEvtDur);
        evtDurRip(p,1:length(binsEvtDur)) = evtDurRip(p,1:length(binsEvtDur))/sum(evtDurRip(p,1:length(binsEvtDur)));
        evtDurRip_tot = [evtDurRip_tot; ripSpkEventTimes.EventDuration];
        
        
        %% Population Level evt rates
        dur=0;
        for i = 1:size(SleepState.ints.NREMstate,1)
            dur=dur+(SleepState.ints.NREMstate(i,2)-SleepState.ints.NREMstate(i,1));
        end
        rateB = cat(1,rateB,numel(barSpkEventTimes.EventDuration)/dur); 

        % ripple duration and rate
        if isfield(SleepState.ints,'nonTHETA')
            for i = 1:size(SleepState.ints.nonTHETA,1)
                dur=dur+(SleepState.ints.nonTHETA(i,2)-SleepState.ints.nonTHETA(i,1));
            end
        end
        rateR = cat(1,rateR,numel(ripSpkEventTimes.EventDuration)/dur); clear dur;

        %% PYR RUNS
        % Absolute FR during barr (per region, box)
        meanPYR = NaN(1,size(spikesPYR.UID,2));
        stdPYR = NaN(1,size(spikesPYR.UID,2));
        spkhist = cell(1,size(spikesPYR.UID,2));
        useUn = cell(1,size(spikesPYR.UID,2));
        absFRpxt = cell(1,size(spikesPYR.UID,2));
        absFRpgt = cell(1,size(spikesPYR.UID,2));
        absFR2xt = cell(1,size(spikesPYR.UID,2));
        absFR2gt = cell(1,size(spikesPYR.UID,2));
        nSpkEach = unitBar.nSpkEach; EventDuration = barSpkEventTimes.EventDuration;
        nSpkEachRip = unitRip.nSpkEach; EventDurationRip = ripSpkEventTimes.EventDuration;
        try
            Narr = cell_metrics.tags.N; 
        catch
            Narr = []; %warning(['No N cells in session: ' basepath]);
        end
        try
            Parr = cell_metrics.tags.P;
        catch
            Parr = []; %warning(['No N cells in session: ' basepath]);
        end
        parfor u = 1:size(spikesPYR.UID,2)
            % Absolute FR during barr (per region, box)
            tSmooth = 0.02;
            binSz = 0.005;
            spkhist{u} = spkRtHist(spikesPYR.times{u}, 'tSmooth',tSmooth,'binSz',binSz,'ifz', false);
            meanPYR(u) = mean(spkhist{u}/binSz); %divide by binsz
            stdPYR(u) = std(spkhist{u}/binSz);
            useUn{u} = find(nSpkEach(u,:) ~= 0);
            if ~isempty(useUn{u})
                absFRpxt{u} = ((nSpkEach(u,useUn{u})'./EventDuration(useUn{u}))-meanPYR(u))/stdPYR(u);
                absFRpgt{u}= ones(length(useUn{u}),1)*regID(spikesPYR.UID(u));
                if regID(spikesPYR.UID(u))==3
                    absFR2xt{u} = absFRpxt{u};
                    absFR2gt{u} = sortMod(u, length(absFRpxt{u}), Narr, Parr);
                end
            end
        end
        absFRpx = [absFRpx; cat(1,absFRpxt{:})]; absFRpg = [absFRpg; cat(1,absFRpgt{:})]; 
        absFR2x = [absFR2x; cat(1,absFR2xt{:})]; absFR2g = [absFR2g; cat(1,absFR2gt{:})]; 
        clear absFRpxt absFRpgt absFR2xt absFR2gt useUn spkhist
        for u = 1:size(spikesPYR.UID,2)
            %barrage abs FR
            gainFRpb{regID(spikesPYR.UID(u))} = addArr(gainFRpb{regID(spikesPYR.UID(u))}, mean(nSpkEach(u,:)'./EventDuration(:)));
            
            %cellType for FR
            if ~isempty(find(cell_metrics.tags.P==spikesPYR.UID(u),1)) %unknown=1;P=2;N=3
                gainFRpbTYPE{regID(spikesPYR.UID(u))} = addArr(gainFRpbTYPE{regID(spikesPYR.UID(u))},2);
            elseif ~isempty(find(cell_metrics.tags.N==spikesPYR.UID(u),1))
                gainFRpbTYPE{regID(spikesPYR.UID(u))} = addArr(gainFRpbTYPE{regID(spikesPYR.UID(u))},3);
            else
                gainFRpbTYPE{regID(spikesPYR.UID(u))} = addArr(gainFRpbTYPE{regID(spikesPYR.UID(u))},1);
            end
            
            %session FR
            gainFRps{regID(spikesPYR.UID(u))} = addArr(gainFRps{regID(spikesPYR.UID(u))}, length(spikesPYR.times{u})/sesDur);
            
            %single gain metric
            gainMet = ((mean(nSpkEach(u,:)'./EventDuration(:)))-(length(spikesPYR.times{u})/sesDur))/((mean(nSpkEach(u,:)'./EventDuration(:)))+(length(spikesPYR.times{u})/sesDur));
            gainMetp{regID(spikesPYR.UID(u))} = addArr(gainMetp{regID(spikesPYR.UID(u))}, gainMet);
            
            gainMetpg{regID(spikesPYR.UID(u))} = addArr(gainMetpg{regID(spikesPYR.UID(u))}, regID(spikesPYR.UID(u)));

            %gain for rips
            gainMet = [];
            gainMet = ((mean(nSpkEachRip(u,:)'./EventDurationRip(:)))-(length(spikesPYR.times{u})/sesDur))/((mean(nSpkEachRip(u,:)'./EventDurationRip(:)))+(length(spikesPYR.times{u})/sesDur));
            gainMetRip{regID(spikesPYR.UID(u))} = addArr(gainMetRip{regID(spikesPYR.UID(u))}, gainMet);

            %burst props
            burstcp{regID(spikesPYR.UID(u))} = addArr(burstcp{regID(spikesPYR.UID(u))}, burstIndex(u));
        end 
        
        
        
        %% INT RUNS
%         fprintf('This message is sent at time %s\n', datestr(now,'HH:MM:SS'));
        % Absolute FR during barr (per region, box)
        meanINT = NaN(1,size(spikesINT.UID,2));
        stdINT = NaN(1,size(spikesINT.UID,2));
        spkhist = cell(1,size(spikesINT.UID,2));
        useUn = cell(1,size(spikesINT.UID,2));
        absFRixt = cell(1,size(spikesINT.UID,2));
        absFRigt = cell(1,size(spikesINT.UID,2));
        nSpkEach = unitBarINT.nSpkEach; EventDuration = barSpkEventTimesINT.EventDuration;
        parfor u = 1:size(spikesINT.UID,2)
            % Absolute FR during barr (per region, box)
            tSmooth = 0.02;
            binSz = 0.005;
            spkhist{u} = spkRtHist(spikesINT.times{u}, 'tSmooth',tSmooth,'binSz',binSz,'ifz', false);
            meanINT(u) = mean(spkhist{u}/binSz); %divide by binsz
            stdINT(u) = std(spkhist{u}/binSz);
            % session FR
            useUn{u} = find(nSpkEach(u,:) ~= 0);
            if ~isempty(useUn{u})
                absFRixt{u} = ((nSpkEach(u,useUn{u})'./EventDuration(useUn{u}))-meanINT(u))/stdINT(u);
                absFRigt{u}= ones(length(useUn{u}),1)*regID(spikesINT.UID(u));
            end
        end
        absFRix = [absFRix; cat(1,absFRixt{:})]; absFRig = [absFRig; cat(1,absFRigt{:})]; clear absFRixt absFRigt useUn spkhist
        for u = 1:size(spikesINT.UID,2)
            %barrage abs FR
            gainFRib{regID(spikesINT.UID(u))} = addArr(gainFRib{regID(spikesINT.UID(u))}, mean(nSpkEach(u,:)'./EventDuration(:)));
            
            %session FR
            gainFRis{regID(spikesINT.UID(u))} = addArr(gainFRis{regID(spikesINT.UID(u))}, length(spikesINT.times{u})/sesDur);
        end
    end 
    
    %% PSTH runs
    if ~skipPSTH
    if exist(strcat(basepath, '\Barrage_Files\PSTHmet\',basename,'.CA1PSTHmets.mat'))
        load(strcat(basepath, '\Barrage_Files\PSTHmet\',basename,'.CA1PSTHmets.mat'));
        tempArr = []; tempArr = catPSTH{1}; 
        startInd = size(tempArr,2)+1; stopInd = startInd + size(PSTHmets.PSTH_bar.responsecurve,2)-1;
        tempArr(1:200, startInd:stopInd) = PSTHmets.PSTH_bar.responsecurve;
        catPSTH{1} = tempArr;
        
        tempArr = []; tempArr = catPSTHsc{1};
        tempArr(1:200, startInd:stopInd) = zscore(PSTHmets.PSTH_bar.responsecurve);
        catPSTHsc{1} = tempArr;
        
        tempArr = []; tempArr = catPSTH_rip{1}; 
        startInd = size(tempArr,2)+1; stopInd = startInd + size(PSTHmets.PSTH_ripples.responsecurve,2)-1;
        tempArr(1:200, startInd:stopInd) = PSTHmets.PSTH_ripples.responsecurve;
        catPSTH_rip{1} = tempArr;
        
        tempArr = []; tempArr = catPSTHsc_rip{1};
        tempArr(1:200, startInd:stopInd) = zscore(PSTHmets.PSTH_ripples.responsecurve);
        catPSTHsc_rip{1} = tempArr;
    end
    if exist(strcat(basepath, '\Barrage_Files\PSTHmet\',basename,'.CA2PSTHmets.mat'))
        load(strcat(basepath, '\Barrage_Files\PSTHmet\',basename,'.CA2PSTHmets.mat'));
        tempArr = []; tempArr = catPSTH{2}; 
        startInd = size(tempArr,2)+1; stopInd = startInd + size(PSTHmets.PSTH_bar.responsecurve,2)-1;
        tempArr(1:200, startInd:stopInd) = PSTHmets.PSTH_bar.responsecurve;
        catPSTH{2} = tempArr;
        
        tempArr = []; tempArr = catPSTHsc{2};
        tempArr(1:200, startInd:stopInd) = zscore(PSTHmets.PSTH_bar.responsecurve);
        catPSTHsc{2} = tempArr;
        
        tempArr = []; tempArr = catPSTH_rip{2}; 
        startInd = size(tempArr,2)+1; stopInd = startInd + size(PSTHmets.PSTH_ripples.responsecurve,2)-1;
        tempArr(1:200, startInd:stopInd) = PSTHmets.PSTH_ripples.responsecurve;
        catPSTH_rip{2} = tempArr;
        
        tempArr = []; tempArr = catPSTHsc_rip{2};
        tempArr(1:200, startInd:stopInd) = zscore(PSTHmets.PSTH_ripples.responsecurve);
        catPSTHsc_rip{2} = tempArr;
    end
    if exist(strcat(basepath, '\Barrage_Files\PSTHmet\',basename,'.CA3PSTHmets.mat'))
        load(strcat(basepath, '\Barrage_Files\PSTHmet\',basename,'.CA3PSTHmets.mat'));
        tempArr = []; tempArr = catPSTH{3}; 
        startInd = size(tempArr,2)+1; stopInd = startInd + size(PSTHmets.PSTH_bar.responsecurve,2)-1;
        tempArr(1:200, startInd:stopInd) = PSTHmets.PSTH_bar.responsecurve;
        catPSTH{3} = tempArr;
        
        tempArr = []; tempArr = catPSTHsc{3};
        tempArr(1:200, startInd:stopInd) = zscore(PSTHmets.PSTH_bar.responsecurve);
        catPSTHsc{3} = tempArr;
        
        tempArr = []; tempArr = catPSTH_rip{3}; 
        startInd = size(tempArr,2)+1; stopInd = startInd + size(PSTHmets.PSTH_ripples.responsecurve,2)-1;
        tempArr(1:200, startInd:stopInd) = PSTHmets.PSTH_ripples.responsecurve;
        catPSTH_rip{3} = tempArr;
        
        tempArr = []; tempArr = catPSTHsc_rip{3};
        tempArr(1:200, startInd:stopInd) = zscore(PSTHmets.PSTH_ripples.responsecurve);
        catPSTHsc_rip{3} = tempArr;
    end
%     if exist(strcat(basepath, '\Barrage_Files\PSTHmet\',basename,'.CTXPSTHmets.mat'))
%         load(strcat(basepath, '\Barrage_Files\PSTHmet\',basename,'.CTXPSTHmets.mat'));
%         tempArr = []; tempArr = catPSTH{4}; 
%         startInd = size(tempArr,2)+1; stopInd = startInd + size(PSTHmets.PSTH_bar.responsecurve,2)-1;
%         tempArr(1:200, startInd:stopInd) = PSTHmets.PSTH_bar.responsecurve;
%         catPSTH{4} = tempArr;
%         
%         tempArr = []; tempArr = catPSTHsc{4};
%         tempArr(1:200, startInd:stopInd) = zscore(PSTHmets.PSTH_bar.responsecurve);
%         catPSTHsc{4} = tempArr;
%         
%         tempArr = []; tempArr = catPSTH_rip{4}; 
%         startInd = size(tempArr,2)+1; stopInd = startInd + size(PSTHmets.PSTH_ripples.responsecurve,2)-1;
%         tempArr(1:200, startInd:stopInd) = PSTHmets.PSTH_ripples.responsecurve;
%         catPSTH_rip{4} = tempArr;
%         
%         tempArr = []; tempArr = catPSTHsc_rip{4};
%         tempArr(1:200, startInd:stopInd) = zscore(PSTHmets.PSTH_ripples.responsecurve);
%         catPSTHsc_rip{4} = tempArr;
%     end
%     if exist(strcat(basepath, '\Barrage_Files\PSTHmet\',basename,'.DGPSTHmets.mat'))
%         load(strcat(basepath, '\Barrage_Files\PSTHmet\',basename,'.DGPSTHmets.mat'));
%         tempArr = []; tempArr = catPSTH{5}; 
%         startInd = size(tempArr,2)+1; stopInd = startInd + size(PSTHmets.PSTH_bar.responsecurve,2)-1;
%         tempArr(1:200, startInd:stopInd) = PSTHmets.PSTH_bar.responsecurve;
%         catPSTH{5} = tempArr;
%         
%         tempArr = []; tempArr = catPSTHsc{5};
%         tempArr(1:200, startInd:stopInd) = zscore(PSTHmets.PSTH_bar.responsecurve);
%         catPSTHsc{5} = tempArr;
%         
%         tempArr = []; tempArr = catPSTH_rip{5}; 
%         startInd = size(tempArr,2)+1; stopInd = startInd + size(PSTHmets.PSTH_ripples.responsecurve,2)-1;
%         tempArr(1:200, startInd:stopInd) = PSTHmets.PSTH_ripples.responsecurve;
%         catPSTH_rip{5} = tempArr;
%         
%         tempArr = []; tempArr = catPSTHsc_rip{5};
%         tempArr(1:200, startInd:stopInd) = zscore(PSTHmets.PSTH_ripples.responsecurve);
%         catPSTHsc_rip{5} = tempArr;
%     end
%     if exist(strcat(basepath, '\Barrage_Files\PSTHmet\',basename,'.MECPSTHmets.mat'))
%         load(strcat(basepath, '\Barrage_Files\PSTHmet\',basename,'.MECPSTHmets.mat'));
%         tempArr = []; tempArr = catPSTH{6}; 
%         startInd = size(tempArr,2)+1; stopInd = startInd + size(PSTHmets.PSTH_bar.responsecurve,2)-1;
%         tempArr(1:200, startInd:stopInd) = PSTHmets.PSTH_bar.responsecurve;
%         catPSTH{6} = tempArr;
%         
%         tempArr = []; tempArr = catPSTHsc{6};
%         tempArr(1:200, startInd:stopInd) = zscore(PSTHmets.PSTH_bar.responsecurve);
%         catPSTHsc{6} = tempArr;
%         
%         tempArr = []; tempArr = catPSTH_rip{6}; 
%         startInd = size(tempArr,2)+1; stopInd = startInd + size(PSTHmets.PSTH_ripples.responsecurve,2)-1;
%         tempArr(1:200, startInd:stopInd) = PSTHmets.PSTH_ripples.responsecurve;
%         catPSTH_rip{6} = tempArr;
%         
%         tempArr = []; tempArr = catPSTHsc_rip{6};
%         tempArr(1:200, startInd:stopInd) = zscore(PSTHmets.PSTH_ripples.responsecurve);
%         catPSTHsc_rip{6} = tempArr;
%     end
%     if exist(strcat(basepath, '\Barrage_Files\PSTHmet\',basename,'.LECPSTHmets.mat'))
%         load(strcat(basepath, '\Barrage_Files\PSTHmet\',basename,'.LECPSTHmets.mat'));
%         tempArr = []; tempArr = catPSTH{7}; 
%         startInd = size(tempArr,2)+1; stopInd = startInd + size(PSTHmets.PSTH_bar.responsecurve,2)-1;
%         tempArr(1:200, startInd:stopInd) = PSTHmets.PSTH_bar.responsecurve;
%         catPSTH{7} = tempArr;
%         
%         tempArr = []; tempArr = catPSTHsc{7};
%         tempArr(1:200, startInd:stopInd) = zscore(PSTHmets.PSTH_bar.responsecurve);
%         catPSTHsc{7} = tempArr;
%         
%         tempArr = []; tempArr = catPSTH_rip{7}; 
%         startInd = size(tempArr,2)+1; stopInd = startInd + size(PSTHmets.PSTH_ripples.responsecurve,2)-1;
%         tempArr(1:200, startInd:stopInd) = PSTHmets.PSTH_ripples.responsecurve;
%         catPSTH_rip{7} = tempArr;
%         
%         tempArr = []; tempArr = catPSTHsc_rip{7};
%         tempArr(1:200, startInd:stopInd) = zscore(PSTHmets.PSTH_ripples.responsecurve);
%         catPSTHsc_rip{7} = tempArr;
%     end
    end
end

regKeyUse = [];
for i = 1:size(regKey,2)
    regKeyUse(i) = str2double(regKey(2,i));
end

%% %%% Plotting %%% %%

%% Event duration histogram
figure('Position', get(0, 'Screensize'));
binsEvtDur = [0:0.01:2];
evt_mean = mean(evtDur,1);
evt_low = evt_mean - (std(evtDur,0,1)/sqrt(size(evtDur,1)));
evt_high = evt_mean + (std(evtDur,0,1)/sqrt(size(evtDur,1)));

hold on
plot(binsEvtDur, evt_mean, 'k');
plot(binsEvtDur, evt_low, 'Color', [0.25 0.25 0.25], 'LineStyle', ':');
plot(binsEvtDur, evt_high, 'Color', [0.25 0.25 0.25], 'LineStyle', ':');
xline(median(evtDur_tot),'r');
xlim([0 2]);
set(gca,'xscale','log');
title('Barrage Event Duration');
ylabel('Normalized Barrage Count');
xlabel('Event duration (s)');
saveas(gcf, ['Z:\home\Lindsay\Barrage\FigPlots\' convertStringsToChars(combine) '.evtDurHist.png']);

%% Ripple Event Duration Histogram
figure('Position', get(0, 'Screensize'));
binsEvtDur = [0:0.01:2];
evtRip_mean = mean(evtDurRip,1);
evtRip_low = evtRip_mean - (std(evtDurRip,0,1)/sqrt(size(evtDurRip,1)));
evtRip_high = evtRip_mean + (std(evtDurRip,0,1)/sqrt(size(evtDurRip,1)));

hold on
plot(binsEvtDur, evtRip_mean, 'k');
plot(binsEvtDur, evtRip_low, 'Color', [0.25 0.25 0.25], 'LineStyle', ':');
plot(binsEvtDur, evtRip_high, 'Color', [0.25 0.25 0.25], 'LineStyle', ':');
xline(median(evtDurRip_tot),'r');
xlim([0 2]);
set(gca,'xscale','log');
title('Ripple Event Duration');
ylabel('Normalized Ripple Count');
xlabel('Event duration (s)');
saveas(gcf, ['Z:\home\Lindsay\Barrage\FigPlots\' convertStringsToChars(combine) '.evtDurHistRip.png']);

%% Barr/Rip Event Duration Overlay
figure('Position', get(0, 'Screensize'));
hold on
plot(binsEvtDur, evt_mean, 'k');
plot(binsEvtDur, evt_low, 'Color', [0.25 0.25 0.25], 'LineStyle', ':');
plot(binsEvtDur, evt_high, 'Color', [0.25 0.25 0.25], 'LineStyle', ':');
xline(median(evtDur_tot),'k');
plot(binsEvtDur, evtRip_mean, 'b');
plot(binsEvtDur, evtRip_low, 'Color', [0.25 0.25 1], 'LineStyle', ':');
plot(binsEvtDur, evtRip_high, 'Color', [0.25 0.25 1], 'LineStyle', ':');
xline(median(evtDurRip_tot),'b');
xlim([0 2]);
set(gca,'xscale','log');
title('Barrage and Ripple Event Durations');
ylabel('Normaized Event Count');
xlabel('Event duration (s)');
saveas(gcf, ['Z:\home\Lindsay\Barrage\FigPlots\' convertStringsToChars(combine) '.evtDurHistComb.png']);

%% Event Rate
boxRatex=cat(1,rateB,rateR);
boxRateg=cat(1,ones(length(rateB),1), 2*ones(length(rateR),1));
figure('Position', get(0, 'Screensize'));
boxplot(boxRatex, boxRateg); title('Event Rate'); ylabel('Event Rate (1/s)');
set(gca,'XTickLabel',{'Barrage','SWR'});
saveas(gcf, ['Z:\home\Lindsay\Barrage\FigPlots\' convertStringsToChars(combine) '.evtRates.png']);

%% PYR Absolute FR during barrages (boxplot)
figure('Position', get(0, 'Screensize'));
boxplot(absFRpx, absFRpg);
[~,useLabels] = intersect(regKeyUse,absFRpg); 
xticklabels(regKey(1,useLabels));
zline = refline(0,0);
zline.LineStyle = '--'; zline.Color = 'k';
title('Pyramidal Cell absolute firing rate during barrages, per region');
ylabel('Absolute Firing Rate (# spikes/second)');
xlabel('Region');
saveas(gcf, ['Z:\home\Lindsay\Barrage\FigPlots\' convertStringsToChars(combine) '.PYRabsFR.png']);

%% INT Absolute FR during barrages (boxplot)
figure('Position', get(0, 'Screensize'));
boxplot(absFRix, absFRig);
[~,useLabels] = intersect(regKeyUse,absFRig); 
xticklabels(regKey(1,useLabels));
zline = refline(0,0);
zline.LineStyle = '--'; zline.Color = 'k';
title('Interneuron absolute firing rate during barrages, per region');
ylabel('Absolute Firing Rate (# spikes/second)');
xlabel('Region');
saveas(gcf, ['Z:\home\Lindsay\Barrage\FigPlots\' convertStringsToChars(combine) '.INTabsFR.png']);

%% CA2Mod Absolute FR during barrages (boxplot)
figure('Position', get(0, 'Screensize'));
boxplot(absFR2x, absFR2g);
modLabels = ["Unknown"; "P"; "N"];
xticklabels(modLabels(unique(absFR2g)));
zline = refline(0,0);
zline.LineStyle = '--'; zline.Color = 'k';
title('Pyramidal CA2 absolute firing rate during barrages, per modulation type');
ylabel('Absolute Firing Rate (# spikes/second)');
xlabel('CA2 Modulation Type');
saveas(gcf, ['Z:\home\Lindsay\Barrage\FigPlots\' convertStringsToChars(combine) '.CA2absFR.png']);

%% PYR FR Gain
figure('Position', get(0, 'Screensize'));
legendUse = [];
hold on
for i = 1:length(regKeyUse)
    if ~isempty(gainFRpb{regKeyUse(i)})
        scatter(gainFRpb{regKeyUse(i)}, gainFRps{regKeyUse(i)});
%         plot(10E-3:10E-3:50, log10(10E-3:10E-3:50), 'Color','k');
        set(gca,'xscale','log', 'yscale','log');
%         set(gca,'xscale','log');
        ylim([0 20]);
        xlim([0 100]);
        legendUse = [legendUse regKey(1,i)];
    end
end
% unity=refline(1,0);
% unity.Color = 'k';
title('Firing Rate Gain PYR');
xlabel('Absolute FR during Barrage');
ylabel('Average FR across session');
legend(legendUse);
plot(1E-3:10E-3:100, (1E-3:10E-3:100), 'Color','k');
hold off
saveas(gcf, ['Z:\home\Lindsay\Barrage\FigPlots\' convertStringsToChars(combine) '.PYRgainFR.png']);

%subplots
figure('Position', get(0, 'Screensize'));
for i = 1:length(regKeyUse)
    subplot(4,2,i);
    hold on
    title([regKey(1,i)]);
    if ~isempty(gainFRpb{regKeyUse(i)})
        c = zeros(length(gainFRpbTYPE{regKeyUse(i)}),3);
        for j = 1:length(gainFRpbTYPE{regKeyUse(i)})
            if gainFRpbTYPE{regKeyUse(i)}(j)==1 %unknown=1;P=2;N=3
                c(j,:) = [0 0 0];
            elseif gainFRpbTYPE{regKeyUse(i)}(j)==2
                c(j,:) = [0 1 0];
            else
                c(j,:) = [1 0 0];
            end
        end
        scatter(gainFRpb{regKeyUse(i)}, gainFRps{regKeyUse(i)}, [], c);
        set(gca,'xscale','log', 'yscale','log');
%         unity=refline(1,0);
%         unity.Color = 'k';
        xlabel('Absolute FR during Barrage');
        ylabel('Average FR across session');
    end
    plot(1E-3:10E-3:100, (1E-3:10E-3:100), 'Color','k');
    ylim([10E-4 20]);
    xlim([10E-4 50]);
    hold off
end
sgtitle('Firing Rate Gain PYR');

saveas(gcf, ['Z:\home\Lindsay\Barrage\FigPlots\' convertStringsToChars(combine) '.PYRgainFRsub.png']);

% single metric per region box
gainMetBar = gainMetp;
gainMetp = cat(1,gainMetp{:}); gainMetpg = cat(1,gainMetpg{:});
figure('Position', get(0, 'Screensize'));
boxplot(gainMetp, gainMetpg);
yline(0, '--');
[~,useLabels] = intersect(regKeyUse,gainMetpg); 
xticklabels(regKey(1,useLabels));
title('Gain Metric per region');
ylabel('(BarrFR-SessFR)/(BarrFR+SessFR)');
xlabel('Region');
saveas(gcf, ['Z:\home\Lindsay\Barrage\FigPlots\' convertStringsToChars(combine) '.PYRgainMetBox.png']);

%% Ripple and Barrage Gain
% gainMetRip = cat(1,gainMetRip{:});
% groupBoxGain = cell(1,length(useLabels));
% j=1;
% for i = 1:length(gainMetBar)
%     if ~isempty(gainMetBar{i})
%         tempArr = [];
%         tempArr(:,1) = gainMetBar{i};
%         tempArr(:,2) = gainMetRip{i};
%         groupBoxGain{j} = tempArr;
%         j= j+1;
%     end
% end
% figure('Position', get(0, 'Screensize'));
% boxplotGroup(groupBoxGain, 'PrimaryLabels',regKey(1,useLabels), 'SecondaryLabels',{'Barrage', 'Ripples'});
% saveas(gcf, ['Z:\home\Lindsay\Barrage\FigPlots\' convertStringsToChars(combine) '.PYRgainMetBoxComb.png']);

%% Try 2
groupBoxx = [];
groupBoxg = [];
xtickUse = [];
xtickLabUse = {};
k = 1;
groupCt = 1;
for i = 1:length(gainMetBar)
    if ~isempty(gainMetBar{i})
        for j = 1:length(gainMetBar{i})
           groupBoxx(k) = gainMetBar{i}(j);
           groupBoxg(k) = groupCt;
           k = k+1;
        end
        xtickUse = [xtickUse; groupCt];
        xtickLabUse{groupCt} = strcat(regKey(1,i),'.Bar');
        groupCt = groupCt+1;
    end
    if ~isempty(gainMetRip{i})
        for j = 1:length(gainMetRip{i})
           groupBoxx(k) = gainMetRip{i}(j);
           groupBoxg(k) = groupCt;
           k = k+1;
        end
        xtickUse = [xtickUse; groupCt];
        xtickLabUse{groupCt} = strcat(regKey(1,i),'.Rip');
%         xtickLabUse{groupCt+1} = '';
        groupCt = groupCt+1;
    end
end
figure('Position', get(0, 'Screensize'));
boxplot(groupBoxx, groupBoxg);
xticks([1:length(xtickLabUse)]); xticklabels(xtickLabUse);
for i = 2:length(xtickLabUse)
    if ~mod(i,2)
        xline(i+0.5,'k--');
    end
end
yline(0,'k--');
xlim([0.5 length(xtickLabUse)+0.5]); ylim([-1.25 1.25]); 
title('Gain Metric for Barrages vs Ripples');
saveas(gcf, ['Z:\home\Lindsay\Barrage\FigPlots\' convertStringsToChars(combine) '.PYRgainMetBoxComb.png']);
%% Interneuron 
figure('Position', get(0, 'Screensize'));
legendUse = [];
hold on
for i = 1:length(regKeyUse)
    if ~isempty(gainFRib{regKeyUse(i)})
        scatter(gainFRib{regKeyUse(i)}, gainFRis{regKeyUse(i)});
        set(gca,'xscale','log', 'yscale','log');
        legendUse = [legendUse regKey(1,i)];
    end
end
% unity=refline(1,0);
% unity.Color = 'k';
title('Firing Rate Gain INT');
xlabel('Absolute FR during Barrage');
ylabel('Average FR across session');
legend(legendUse);
plot(1E-3:10E-3:100, (1E-3:10E-3:100), 'Color','k');
% ylim([0 15]);
% xlim([0 100]);
hold off
saveas(gcf, ['Z:\home\Lindsay\Barrage\FigPlots\' convertStringsToChars(combine) '.INTgainFR.png']);

%subplots
figure('Position', get(0, 'Screensize'));
for i = 1:length(regKeyUse)
    subplot(4,2,i);
    hold on
    title([regKey(1,i)]);
    if ~isempty(gainFRib{regKeyUse(i)})
        scatter(gainFRib{regKeyUse(i)}, gainFRis{regKeyUse(i)});
        set(gca,'xscale','log', 'yscale','log');
%         unity=refline(1,0);
%         unity.Color = 'k';
        xlabel('Absolute FR during Barrage');
        ylabel('Average FR across session');
    end
    plot(1E-3:10E-3:100, (1E-3:10E-3:100), 'Color','k');
%     ylim([0 15]);
%     xlim([0 100]);
    hold off
end

sgtitle('Firing Rate Gain INT');
saveas(gcf, ['Z:\home\Lindsay\Barrage\FigPlots\' convertStringsToChars(combine) '.INTgainFRsub.png']);

%% Burst Properties PYR ONLY
figure('Position', get(0, 'Screensize'));
for i = 1:length(regKeyUse)
    subplot(2,4,i);
    hold on
    title([regKey(1,i)]);
    if ~isempty(burstcp{regKeyUse(i)})
        scatter(gainFRpb{regKeyUse(i)}, burstcp{regKeyUse(i)});
        [m, b, r2] = simpleLin(length(gainFRpb{regKeyUse(i)}), gainFRpb{regKeyUse(i)}, burstcp{regKeyUse(i)});
        plot([10E-4:10E-3:100], (m*[10E-4:10E-3:100]) + b, '--k');
        text(15E-4, 8, ['m = ' num2str(m) '  r2 = ' num2str(r2)]);
        set(gca,'xscale','log', 'yscale','log','PlotBoxAspectRatio',[1,1,1]);
        xlabel('Absolute FR during Barrage');
        ylabel('Burst Index');
        xlim([10E-4 100]); 
        ylim([10E-4 100]);
    else
        set(gca,'PlotBoxAspectRatio',[1,1,1]);
    end
    hold off
end
sgtitle('Absolute FR during Barrage vs Burst Index of Cells');
saveas(gcf, ['Z:\home\Lindsay\Barrage\FigPlots\' convertStringsToChars(combine) '.PYRcellBurstsub.png']);


%% Big combined PSTH
if ~skipPSTH
load(['Z:\home\Lindsay\Barrage\FigPlots\bigPSTH.mat']);
load(['Z:\home\Lindsay\Barrage\FigPlots\bigPSTH_rip.mat']);
xlab = 0:10:size(bigPSTH.z_score_CA1_2plot);
actTime = [-0.5 -0.3 -0.1 0.1 0.3 0.5];
xticklab=[];
for xCt = 1:6
    xticklab{xCt} = num2str(actTime(xCt));
end
figure('Position', get(0, 'Screensize'));
subplot(3,3,1);
plot(bigPSTH.z_score_CA1_2plot,'k');hold on;plot(bigPSTH.z_score_CA1_2plot+bigPSTH.z_score_CA1_2plot_err,'--k');hold on;
plot(bigPSTH.z_score_CA1_2plot-bigPSTH.z_score_CA1_2plot_err,'--k');hold on;
plot(bigPSTH_rip.z_score_CA1_rip_2plot,'b');hold on;plot(bigPSTH_rip.z_score_CA1_rip_2plot+bigPSTH_rip.z_score_CA1_rip_2plot_err,'--b');hold on;
plot(bigPSTH_rip.z_score_CA1_rip_2plot-bigPSTH_rip.z_score_CA1_rip_2plot_err,'--b');hold on;
xlim([1 49]);ylim([-0.1 0.3]);
xticks(xlab); xticklabels(xticklab); title('CA1');hold on;
line([round(length(bigPSTH.bin_time)/2)-1 round(length(bigPSTH.bin_time)/2)-1],[-0.1 0.3],'Color','red','LineStyle','--');

subplot(3,3,2);
plot(bigPSTH.z_score_CA2_2plot,'k');hold on;plot(bigPSTH.z_score_CA2_2plot+bigPSTH.z_score_CA2_2plot_err,'--k');hold on;
plot(bigPSTH.z_score_CA2_2plot-bigPSTH.z_score_CA2_2plot_err,'--k');hold on;
plot(bigPSTH_rip.z_score_CA2_rip_2plot,'b');hold on;plot(bigPSTH_rip.z_score_CA2_rip_2plot+bigPSTH_rip.z_score_CA2_rip_2plot_err,'--b');hold on;
plot(bigPSTH_rip.z_score_CA2_rip_2plot-bigPSTH_rip.z_score_CA2_rip_2plot_err,'--b');hold on;
xlim([1 49]);ylim([-0.1 0.3]);
xticks(xlab); xticklabels(xticklab);title('CA2');hold on;
line([round(length(bigPSTH.bin_time)/2)-1 round(length(bigPSTH.bin_time)/2)-1],[-0.1 0.3],'Color','red','LineStyle','--');

subplot(3,3,3);
plot(bigPSTH.z_score_CA3_2plot,'k');hold on;plot(bigPSTH.z_score_CA3_2plot+bigPSTH.z_score_CA3_2plot_err,'--k');hold on;
plot(bigPSTH.z_score_CA3_2plot-bigPSTH.z_score_CA3_2plot_err,'--k');hold on;
plot(bigPSTH_rip.z_score_CA3_rip_2plot,'b');hold on;plot(bigPSTH_rip.z_score_CA3_rip_2plot+bigPSTH_rip.z_score_CA3_rip_2plot_err,'--b');hold on;
plot(bigPSTH_rip.z_score_CA3_rip_2plot-bigPSTH_rip.z_score_CA3_rip_2plot_err,'--b');hold on;
xlim([1 49]);ylim([-0.1 0.3]);
xticks(xlab); xticklabels(xticklab);title('CA3');hold on;
line([round(length(bigPSTH.bin_time)/2)-1 round(length(bigPSTH.bin_time)/2)-1],[-0.1 0.3],'Color','red','LineStyle','--');

cmin = mean([min(bigPSTH.z_score_CA1_2plot) min(bigPSTH.z_score_CA2_2plot) min(bigPSTH.z_score_CA3_2plot) min(bigPSTH_rip.z_score_CA1_rip_2plot) min(bigPSTH_rip.z_score_CA2_rip_2plot) min(bigPSTH_rip.z_score_CA3_rip_2plot)]);
cmax = mean([max(bigPSTH_rip.z_score_CA1_rip_2plot) max(bigPSTH.z_score_CA2_2plot) max(bigPSTH_rip.z_score_CA3_rip_2plot)]);

subplot(3,3,4);imagesc(bigPSTH.z_score_CA1); xlim([1 49]); xticks(xlab); xticklabels(xticklab); caxis([cmin cmax]); xline(25, 'r--');
subplot(3,3,5);imagesc(bigPSTH.z_score_CA2); xlim([1 49]); xticks(xlab); xticklabels(xticklab); caxis([cmin cmax]); xline(25, 'r--'); title('Barrage PSTH');
subplot(3,3,6);imagesc(bigPSTH.z_score_CA3); xlim([1 49]); xticks(xlab); xticklabels(xticklab); caxis([cmin cmax]); xline(25, 'r--');
subplot(3,3,7);imagesc(bigPSTH_rip.z_score_CA1_rip); xlim([1 49]); xticks(xlab); xticklabels(xticklab); caxis([cmin cmax]); xline(25, 'r--');
subplot(3,3,8);imagesc(bigPSTH_rip.z_score_CA2_rip); xlim([1 49]); xticks(xlab); xticklabels(xticklab); caxis([cmin cmax]); xline(25, 'r--'); title('Ripple PSTH');
subplot(3,3,9);imagesc(bigPSTH_rip.z_score_CA3_rip); xlim([1 49]); xticks(xlab); xticklabels(xticklab); caxis([cmin cmax]); xline(25, 'r--'); colormap(parula); colormap(viridis);
saveas(gcf, ['Z:\home\Lindsay\Barrage\FigPlots\' convertStringsToChars(combine) '.bigPSTHcomb.png']);
end
%% computePSTH Big combined PSTH
if ~skipPSTH
xlab = 0:25:200;
actTime = [-2:0.5:2];
xticklab=[];
for xCt = 1:5
    xticklab{xCt} = num2str(actTime(xCt));
end
figure('Position', get(0, 'Screensize'));
subplot(3,3,1);
plot(PSTHmets.PSTH_bar.time, mean(catPSTH{1}'),'k');hold on;plot(PSTHmets.PSTH_bar.time, mean(catPSTH{1}')+std(catPSTH{1}'),'--k');hold on;
plot(PSTHmets.PSTH_bar.time, mean(catPSTH{1}')-std(catPSTH{1}'),'--k');hold on;
plot(PSTHmets.PSTH_bar.time, mean(catPSTH_rip{1}'),'b');hold on;plot(PSTHmets.PSTH_bar.time, mean(catPSTH_rip{1}')+std(catPSTH_rip{1}'),'--b');hold on;
plot(PSTHmets.PSTH_bar.time, mean(catPSTH_rip{1}')-std(catPSTH_rip{1}'),'--b');hold on;
xline(0,'Color','red','LineStyle','--');
xlim([-.5 .5]);ylim([-5 15]);
title('CA1');hold on;

subplot(3,3,2);
plot(PSTHmets.PSTH_bar.time, mean(catPSTH{2}'),'k');hold on;plot(PSTHmets.PSTH_bar.time, mean(catPSTH{2}')+std(catPSTH{2}'),'--k');hold on;
plot(PSTHmets.PSTH_bar.time, mean(catPSTH{2}')-std(catPSTH{2}'),'--k');hold on;
plot(PSTHmets.PSTH_bar.time, mean(catPSTH_rip{2}'),'b');hold on;plot(PSTHmets.PSTH_bar.time, mean(catPSTH_rip{2}')+std(catPSTH_rip{2}'),'--b');hold on;
plot(PSTHmets.PSTH_bar.time, mean(catPSTH_rip{2}')-std(catPSTH_rip{2}'),'--b');hold on;
xline(0,'Color','red','LineStyle','--');
xlim([-.5 .5]);ylim([-5 15]);
title('CA2');hold on;

subplot(3,3,3);
plot(PSTHmets.PSTH_bar.time, mean(catPSTH{3}'),'k');hold on;plot(PSTHmets.PSTH_bar.time, mean(catPSTH{3}')+std(catPSTH{3}'),'--k');hold on;
plot(PSTHmets.PSTH_bar.time, mean(catPSTH{3}')-std(catPSTH{3}'),'--k');hold on;
plot(PSTHmets.PSTH_bar.time, mean(catPSTH_rip{3}'),'b');hold on;plot(PSTHmets.PSTH_bar.time, mean(catPSTH_rip{3}')+std(catPSTH_rip{3}'),'--b');hold on;
plot(PSTHmets.PSTH_bar.time, mean(catPSTH_rip{3}')-std(catPSTH_rip{3}'),'--b');hold on;
xline(0,'Color','red','LineStyle','--');
xlim([-.5 .5]);ylim([-5 15]);
title('CA3');hold on;

% cmin = mean([min(mean(catPSTH{1}')) min(mean(catPSTH{2}')) min(mean(catPSTH{3}')) min(mean(catPSTH_rip{1}')) min(mean(catPSTH_rip{2}')) min(mean(catPSTH_rip{3}'))]);
cmax = mean([max(mean(catPSTH_rip{1}')) max(mean(catPSTH{2}')) max(mean(catPSTH_rip{3}'))]);

subplot(3,3,4);imagesc(catPSTHsc{1}'); caxis([-3 cmax]); xlim([75 125]); xticks(xlab); xticklabels(xticklab); xline(100, 'r--');
subplot(3,3,5);imagesc(catPSTHsc{2}'); caxis([-3 cmax]); xlim([75 125]); xticks(xlab); xticklabels(xticklab); xline(100, 'r--'); title('Barrage PSTH');
subplot(3,3,6);imagesc(catPSTHsc{3}'); caxis([-3 cmax]); xlim([75 125]); xticks(xlab); xticklabels(xticklab); xline(100, 'r--');
subplot(3,3,7);imagesc(catPSTHsc_rip{1}'); caxis([-3 cmax]); xlim([75 125]); xticks(xlab); xticklabels(xticklab); xline(100, 'r--');
subplot(3,3,8);imagesc(catPSTHsc_rip{2}'); caxis([-3 cmax]); xlim([75 125]); xticks(xlab); xticklabels(xticklab); xline(100, 'r--'); title('Ripple PSTH');
subplot(3,3,9);imagesc(catPSTHsc_rip{3}'); caxis([-3 cmax]); xlim([75 125]); xticks(xlab); xticklabels(xticklab); xline(100, 'r--'); colormap(parula); colormap(viridis);

saveas(gcf, ['Z:\home\Lindsay\Barrage\FigPlots\' convertStringsToChars(combine) '.bigCompPSTHcomb.png']);
end

close all

end

function updatedArr = addArr(original, addition)
updatedArr = []; 
updatedArr = original;
updatedArr = [updatedArr; addition];
end
function group = sortMod(u, lenGroup, Narr, Parr)
    if ~isempty(find(Narr==u,1))
        group = ones(lenGroup,1)*3;
    elseif ~isempty(find(Parr==u,1))
        group = ones(lenGroup,1)*2;
    else
        group = ones(lenGroup,1)*1;
    end
end

function [m,b,r2] = simpleLin(n, x, y)
sum_xy = 0; sum_x = 0; sum_y = 0; sum_x2 = 0;
for i = 1:length(x)
    sum_xy = sum_xy + (x(i)*y(i));
    sum_x = sum_x + x(i);
    sum_y = sum_y + y(i);
    sum_x2 = sum_x2 + (x(i)*x(i));
end
m = ((n*sum_xy)-(sum_x*sum_y))/((n*sum_x2)-(sum_x*sum_x));
b = (sum_y - (m*sum_x))/n;

y_resid = 0;
for i = 1:length(x)
    y_resid = y(i) - ((m*x(i)) + b);
end
SS_resid = sum(y_resid.^2);
SS_total = (length(y)-1)*var(y);
r2 = 1-(SS_resid/SS_total);
end