%% Clear prior saves for new cumulative metrics
% load('Z:\home\Lindsay\Barrage\combinedPaths.mat');
% paths_save = [];
% save('Z:\home\Lindsay\Barrage\combinedPaths.mat', 'paths_save');
% load('Z:\home\Lindsay\Barrage\ratPaths.mat');
% paths_save = [];
% save('Z:\home\Lindsay\Barrage\ratPaths.mat', 'paths_save');
% load('Z:\home\Lindsay\Barrage\mousePaths.mat');
% paths_save = [];
% save('Z:\home\Lindsay\Barrage\mousePaths.mat', 'paths_save');

%% Set paths to be run
useSess = readtable('Z:\home\Lindsay\Barrage\rat_sessions.csv');
% readtable('Z:\home\Lindsay\Barrage\rat_sessions.csv');

%% Set main parameters
ifHSE = 0; %Run detection
ifUseMet = 1; %Load in previous metrics
ifDetUn = 1; %Choose only high firing rate units
ifPickUn = 0; %Choose only units that fire many spikes during events
ifPare = 1; %Keep only the events with so many spikes
ifPSTH = 0; %Get analytical plots for each run
ifAnalysis = 0; %Run analytical plots
ifBigPSTH = 0;
ifCum = 1; %Run cumulative metrics analysis
    ifNewRip = 0;
    ifNewBar = 0;
    ifSkipPSTH = 0;
bigCCG = 0; %Get population CCG
bigDur = 0; %Get population event duration distribution (box plot)

if ~ifUseMet
    useMet.nSigma = 4;
    useMet.tSmooth = 0.02;
    useMet.binSz = 0.005;
    useMet.tSepMax = 0.01; %was playing with 0.005, 0.1
    useMet.mindur = 0.05;
    useMet.maxdur = 10;
    useMet.lastmin = 0.05;
    useMet.sstd = -1*(useMet.nSigma-0.1); %sets start to nSig+sstd %-1*(useMet.nSigma-2);
    useMet.estd = (useMet.nSigma-0.1); %sets end to nSig-sstd %(useMet.nSigma-2);
    useMet.EMGThresh = 0.8;
    
    useMet.Hz = 45; useMet.ft = 0.1; useMet.numEvt = 2;
    useMet.bound = 2.5; useMet.spkThresh = 0.3;
    
    useMet.unMin = 2; useMet.spkNum = 5; useMet.spkHz = 70; useMet.unMax = 0;
    useMet.DetUn = ifDetUn; useMet.PickUn = ifPickUn; useMet.ifPare = ifPare;
end

%% Iterate through our sessions
for p = 1:size(useSess,1)
    % Set naming conventions for session p
    cd(useSess.basepath{p});
    curPath = [useSess.animal{p} '\' useSess.basename{p}];
    basepath = useSess.basepath{p};
    basename = useSess.basename{p};
    animName = useSess.animal{p};
    if ~exist(strcat(basepath,'\','Barrage_Files'))
        mkdir('Barrage_Files');
    end
    savePath = convertStringsToChars(strcat(basepath, '\Barrage_Files\', basename, '.'));
    if ifHSE
        mkSpks([basepath '\Barrage_Files']); %this is still bulky, add check to see if already made
    end
    if exist(strcat(basepath,'\Barrage_Files\',basename,'.CA2pyr.cellinfo.mat'))
        load(strcat(basepath,'\Barrage_Files\',basename,'.CA2pyr.cellinfo.mat'));
        if ~isempty(spikes.UID)
            CA2pyrRun = 1;
        else
            CA2pyrRun = 0;
        end
    else
        CA2pyrRun = 0;
    end
    if ifHSE
        if CA2pyrRun %this only loads if it exists, is it empty
            if ifUseMet
                load([savePath 'useMetNew.mat']);
                if ~isfield(useMet,'DetUn') %catch if this previously wasn't saved
                    useMet.DetUn = 1;
                    useMet.PickUn = 0;
                    useMet.ifPare = 1;
                    save([savePath 'useMetNew.mat'],'useMet');
                end
            end
            if ~exist(strcat(basepath,'\',basename,'.SleepState.states.mat'))
                SleepState = SleepScoreMaster(basepath);
            else
                load(strcat(basepath,'\',basename,'.SleepState.states.mat'));
            end
            
            % Detect high firing units
            if ifDetUn||(ifUseMet&&useMet.DetUn)
                unitsForDetection(useMet.Hz,useMet.ft,useMet.numEvt);
                loadPath = strcat(basepath,'\Barrage_Files\',basename,'.useSpk.cellinfo.mat');
                
                %Only keep units that fire enough spikes in detected events
                if ifPickUn||(ifUseMet&&useMet.PickUn)
                    bound = useMet.bound; useSpkThresh = useMet.spkThresh;
                    [~,useSpkThresh,~,bound] = pickUns(bound,useSpkThresh,0.05,useMet,loadPath, useMet.Hz, useMet.ft, useMet.numEvt);
                    [~,normSpkThresh,~,bound] = pickUns(bound,useSpkThresh,-0.01,useMet,loadPath, useMet.Hz, useMet.ft, useMet.numEvt);
                    checking = burstCellDetection(normSpkThresh,loadPath);
                    useMet.bound = bound;
                    useMet.spkThresh = normSpkThresh;
                    note_all = "PickUns";
                else
                    note_all = "DetectUn";
                end
            else
                note_all = "No pre-vetting";
            end
            
            if exist([savePath 'HSEfutEVT.mat'])
                load([savePath 'HSEfutEVT.mat']);
                runLast = size(evtSave,1);
            else
                runLast = 0;
            end
            runNum = runLast+1;
            load([basepath '\' basename '.cell_metrics.cellinfo.mat']);
            if ifDetUn
                load([savePath 'useSpk.UIDkeep.mat']);
                load([savePath 'useSpk.cellinfo.mat']); %load in the spikes that we've picked and run to save
            else
                load([savePath 'CA2pyr.cellinfo.mat']);
                UIDkeep = spikes.UID;
                save([savePath 'useSpk.UIDkeep.mat'],'UIDkeep');
            end
            if ifDetUn&&isempty(spikes.UID)
                load([savePath 'CA2pyr.cellinfo.mat']);
                UIDkeep = spikes.UID;
                save([savePath 'useSpk.UIDkeep.mat'],'UIDkeep');
                warning('unitsForDetection did not return any units, defaulting to full set');
                note_all = "Failed DetectUn";
            end
            HSE = find_HSE_b('spikes',spikes,...
                'nSigma',useMet.nSigma,'binSz',useMet.binSz,...
                'tSmooth',useMet.tSmooth,'tSepMax',useMet.tSepMax,...
                'mindur',useMet.mindur,'maxdur',useMet.maxdur,...
                'lastmin',useMet.lastmin,'EMGThresh',useMet.EMGThresh,...
                'Notes',note_all,'sstd',useMet.sstd,'estd',useMet.estd,...
                'recordMetrics',true,'neuro2',true,'runNum',runNum);
            if ifPare||(ifUseMet&&useMet.ifPare)
                HSE = pareEvt(HSE, spikes, useMet.unMin, useMet.spkNum, useMet.spkHz, useMet.unMax);
            else
                HSE = pareEvt(HSE, spikes, 0, 0, 0, 0);
            end
            % Save
            save([savePath 'useMetNew.mat'],'useMet');
        else
            warning(strcat('No CA2pyr for ', paths(p), ' :('));
        end
    end
    
    if ~ifHSE||(CA2pyrRun&&(~isempty(SleepState.ints.NREMstate))&&(~isempty(HSE.NREM)))
        if ifAnalysis
            BarAnalysis(useSess.basepath{p});
        end
        if ifPSTH
            regionPSTH();
            close all
        end
        if bigDur
            load(strcat(basepath,'\Barrage_Files\',basename,'.HSE.mat'));
            if p==1
                useGroup = [];
                ugc = 1;
                totDurX = HSE.timestamps(:,2)-HSE.timestamps(:,1);
                useGroup{ugc} = animName;
                totDurG = ugc*ones(size(totDurX,1),1);
                totDurLab = [curPath];
            else
                tempName = [];
                tempName = animName;
                if find(contains(useGroup, tempName))
                    ugc = find(contains(useGroup, tempName));
                else
                    ugc = length(useGroup) + 1;
                    useGroup{ugc} = animName;
                end
                tempDur = []; tempDur = HSE.timestamps(:,2)-HSE.timestamps(:,1);
                totDurX = [totDurX; tempDur];
                totDurG = [totDurG; ugc*ones(size(tempDur,1),1)];
                totDurLab = [totDurLab curPath];
            end
        end
        if bigCCG
            load([savePath 'CCG_dat.mat']);
            if p == 1
                CCGsum = zeros(length(CCG_dat.time),size(useSess,1));
            end
            CCGsum(:,p) = zscore(CCG_dat.y);
        end
    end
end

if ifBigPSTH
    
end

if ifCum
    if contains(useSess.basepath{1},'Y:\')
        NewMetComb([ifNewRip ifNewBar ifSkipPSTH], "mouse");
    else
        NewMetComb([ifNewRip ifNewBar ifSkipPSTH], "rat");
    end
end
if bigDur
    figure(2);
    boxplot(totDurX,totDurG);title('Barrage Duration Across Animals');ylabel('Duration (s)');
    xticklabels(useGroup);
    if contains(useSess.basepath{1}, 'Z:\')
        saveas(gcf,'Z:\home\Lindsay\Barrage\ratEvtBox.png');
    else
        saveas(gcf,'Z:\home\Lindsay\Barrage\mouseEvtBox.png');
    end
end
if bigCCG
    CCG_mean = mean(CCGsum,2);
    CCG_low = CCG_mean - (std(CCGsum,0,2)/sqrt(size(CCGsum,2)));
    CCG_high = CCG_mean + (std(CCGsum,0,2)/sqrt(size(CCGsum,2)));
    figure(3);
    hold on
    plot(CCG_dat.time, CCG_mean, 'k');
    plot(CCG_dat.time,CCG_low, 'Color', [0.25 0.25 0.25], 'LineStyle', ':');
    plot(CCG_dat.time,CCG_high, 'Color', [0.25 0.25 0.25], 'LineStyle', ':');
    if contains(useSess.basepath{1}, 'Z:\')
        saveas(gcf,'Z:\home\Lindsay\Barrage\ratCCG.png');
    else
        saveas(gcf,'Z:\home\Lindsay\Barrage\mouseCCG.png');
    end
end
