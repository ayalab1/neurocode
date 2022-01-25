%% MassSessRun
% Most up-to-date metrics for running analysis through find_HSE_b

% Just run rats or mice together for ease to start

% If we need to clear our path_saves for a completely new cumulative run
% load('Z:\home\Lindsay\Barrage\combinedPaths.mat');
% paths_save = [];
% save('Z:\home\Lindsay\Barrage\combinedPaths.mat', 'paths_save');
% load('Z:\home\Lindsay\Barrage\ratPaths.mat');
% paths_save = [];
% save('Z:\home\Lindsay\Barrage\ratPaths.mat', 'paths_save');
% load('Z:\home\Lindsay\Barrage\mousePaths.mat');
% paths_save = [];
% save('Z:\home\Lindsay\Barrage\mousePaths.mat', 'paths_save');


% No classifications:
% ["AYA9\day12"; "AYA10\day31";
% No CA2:
% ["AB1\day1"];

% DONE:
main = "Z:\Data\AYAold\";
paths = ["AYA7\day19"];

% paths = ["AB4\day03"; "AB4\day07"; "AB4\day08"; "AB4\day09"; "AB4\day11";...
%     "AYA6\day17"; "AYA6\day19"; "AYA6\day20"; "AYA7\day19";...
%     "AYA7\day20"; "AYA7\day22"; "AYA7\day24"; "AYA7\day25";...
%     "AYA7\day27"; "AYA7\day30";  %broke on AYA7/day27

%     paths = ["AYA9\day15"; "AYA9\day16"; "AYA9\day17";...
%     "AYA9\day20"; "AYA10\day25"; "AYA10\day27"; "AYA10\day32"; "AYA10\day34"];


% main = "Y:\SMproject\";
% paths = ["AZ1\day13"; "AO50\day20"; "AO50\day21"; "AO50\day22"; "AO50\day23";...
%     "AO51\day18"; "AO51\day19"; "AO51\day20"; "AO51\day21"];
% paths = ["AZ1\day13"];

load('C:\Users\Cornell\Documents\GitHub\neurocode\projects\BarragePipeline\curRatMet.mat');
bigSave = 'Z:\home\Lindsay\Barrage\mousePaths.mat'; %change to mouse, potentially - change below as well
comSave = 'Z:\home\Lindsay\Barrage\combinedPaths.mat';

ifHSE = 1;
ifAnalysis = 1;
ifCum = 0;
ifPSTH = 1;
ifNeuro1 = 0;
manual = [];
for p = 1:length(paths)
    cd(strcat(main,paths(p)));
    curPath = convertStringsToChars(paths(p));
    basepath = convertStringsToChars(strcat(main,curPath));
    basename = convertStringsToChars(basenameFromBasepath(basepath));

    if ~exist(strcat(basepath,'\','Barrage_Files'))
        mkdir('Barrage_Files');
    end

    savePath = convertStringsToChars(strcat(basepath, '\Barrage_Files\', basename, '.'));
    if ifHSE 
        %% Set metrics to use
        useMet.nSigma = 3;
        useMet.tSmooth = 0.02;
        useMet.binSz = 0.005;
        useMet.tSepMax = 0.005; %was playing with 0.005, 0.1
        useMet.mindur = 0.2;
        useMet.maxdur = 10;
        useMet.lastmin = 0.2;
        useMet.sstd = -1*(useMet.nSigma-2);
        useMet.estd = (useMet.nSigma-2);
        useMet.EMGThresh = 0.8;
        %% New set of spikes and pick units to use for detection
        if ~isempty(manual)
            load([savePath 'CA2pyr.cellinfo.mat']);
            keep = [];
            ki = 1;
            for i = 1:length(manual)
                if ismember(spikes.UID, manual(i))
                    kInd = find(spikes.UID==manual(i),1);
                    UIDkeep(ki) = spikes.UID(kInd);
                    keep{ki} = spikes.times{find(spikes.UID==manual(i),1)};
                    ki = ki+1;
                else
                    warning('Selected unit is not a CA2 pyr');
                end
            end
            clear spikes
            spikes = keep;
            save([savePath 'brstDt.cellinfo.mat'], 'spikes');
            save([savePath 'brstDt.UIDkeep.mat'],'UIDkeep');
        else
            mkSpks([basepath '\Barrage_Files']);
            Hz = 60; ft = 0.2;
            unitsForDetection(Hz, ft); %threshold with FR
            loadPath = strcat(basepath,'\Barrage_Files\',basename,'.brstDt.cellinfo.mat');

            %% Pare detection units
            bound = 3;
            useSpkThresh = 0.3;
            [~,useSpkThresh,~,bound] = pickUns(bound,useSpkThresh,0.05,useMet,loadPath, Hz, ft);
            % if it's passed, let's keep pushing the boundary to see if we can make it better
            [~,normSpkThresh,~,bound] = pickUns(bound,useSpkThresh,-0.01,useMet,loadPath, Hz, ft);
            checking = burstCellDetection(normSpkThresh,loadPath); %output for my own sanity
        end
        
        useMet.Hz = Hz;
        useMet.ft = ft; 
        useMet.bound = bound;
        useMet.spkThresh = normSpkThresh;
        save([savePath 'useMet.mat'],'useMet');
        
        %% Real detection
        if exist([savePath 'HSEfutEVT.mat'])
            load([savePath 'HSEfutEVT.mat']);
            runLast = size(evtSave,1);
        else
            runLast = 0;
        end
        runNum = runLast+1;
        load([basepath '\' basename '.cell_metrics.cellinfo.mat']);
        note_all = "Burst detected";
        load([savePath 'brstDt.UIDkeep.mat']);
%         for n = 1:length(UIDkeep)
            load([savePath 'brstDt.cellinfo.mat']); %load in the spikes that we've picked and run to save
%             keep.times = spikes.times{n};
%             keep.UID = spikes.UID(n);
%             spikes = []; spikes = keep;
            
            nSigma = useMet.nSigma;
            tSmooth = useMet.tSmooth;
            binSz = useMet.binSz;
            tSepMax = useMet.tSepMax; %was playing with 0.005, 0.1
            mindur = useMet.mindur;
            maxdur = useMet.maxdur;
            lastmin = useMet.lastmin;
            sstd = useMet.sstd;
            estd = useMet.estd;
            EMGThresh = useMet.EMGThresh;
            recordMetrics = true;
            neuro2 = true;
            notes = note_all;

            HSE = find_HSE_b('spikes',spikes,...
                        'nSigma',nSigma,'binSz',binSz,'tSmooth',tSmooth,'tSepMax',tSepMax,...
                        'mindur',mindur,'maxdur',maxdur,'lastmin',lastmin,...
                        'EMGThresh',EMGThresh,'Notes',notes,'sstd',sstd,'estd',estd,...
                        'recordMetrics',recordMetrics,'neuro2',neuro2,'runNum',runNum);
                    
            if ifNeuro1
                load([savePath 'HSEfutEVT.mat']);
                trial = size(evtSave,1);
                current_time = evtSave{trial,1};
                current_peak = evtSave{trial,2};
                createEVT(current_time(:,1), current_peak, current_time(:,2), 'saveName', 'H', 'savePath', strcat(pwd,'\Barrage_Files'));
            end
%             keepTimes{n} = HSE.timestamps;
%             keepPeaks{n} = HSE.peaks;
%         end
%         j=1;
%         goodTimes = []; goodPeaks = [];
%         for i = 1:length(keepTimes)
%             if ~isempty(keepTimes{i})
%                 goodTimes{j} = keepTimes{i};
%                 goodPeaks{j} = keepPeaks{i};
%                 j = j+1;
%             end
%         end
%         keepTimes = []; keepPeaks = []; keepTimes = goodTimes; keepPeaks = goodPeaks;
%         tSepMax = 0.005;
%         for n = 1:size(keepTimes,2)-1
%             if n==1
%                 evtstart = cat(1, keepTimes{n}(:,1), keepTimes{n+1}(:,1));
%                 evtstop = cat(1, keepTimes{n}(:,2), keepTimes{n+1}(:,2));
%                 evtpeak = cat(1, keepPeaks{n}, keepPeaks{n+1});
%                 [evtstart, ind] = sort(evtstart);
%                 tempstop = []; temppeak = [];
%                 for e = 1:length(ind)
%                     tempstop(e) = evtstop(ind(e));
%                     temppeak(e) = evtpeak(ind(e));
%                 end
%                 evtstop = tempstop; evtpeak = temppeak;
%                 evtdur = evtstop-evtstart;
%                 evtamp = zeros(length(evtstart),1);
%                 diffConc = [];
%                 for e = 1:length(evtstart)-1
%                     diffConc(e) = evtstart(e+1)-evtstop(e); %start of next-stop of last
%                 end
%                 flagConc = (diffConc <= tSepMax); %flag indices where event distance is good
%                 [start,stop,~,peak,numCat] = CatCon(evtstart,evtstop,evtpeak,evtamp,flagConc);
%                 start = start'; stop = stop'; peak = peak';
%             else
%                 evtstart = cat(1, start, keepTimes{n+1}(:,1));
%                 evtstop = cat(1, stop, keepTimes{n+1}(:,2));
%                 evtpeak = cat(1, peak, keepPeaks{n+1});
%                 [evtstart, ind] = sort(evtstart);
%                 tempstop = []; temppeak = [];
%                 for e = 1:length(ind)
%                     tempstop(e) = evtstop(ind(e));
%                     temppeak(e) = evtpeak(ind(e));
%                 end
%                 evtstop = tempstop; evtpeak = temppeak;
%                 evtdur = evtstop-evtstart;
%                 evtamp = zeros(length(evtstart),1);
%                 diffConc = [];
%                 for e = 1:length(evtstart)-1
%                     diffConc(e) = evtstart(e+1)-evtstop(e); %start of next-stop of last
%                 end
%                 flagConc = (diffConc <= tSepMax); %flag indices where event distance is good
%                 [start,stop,~,peak,numCat] = CatCon(evtstart,evtstop,evtpeak,evtamp,flagConc);
%                 start = start'; stop = stop'; peak = peak';
%             end
%         end
%         finThresh = 0.4;
%         finDur = stop-start;
%         start = start(finDur >= finThresh);
%         stop = stop(finDur >= finThresh);
%         peak = peak(finDur >= finThresh);
%         createEVT(start, peak, stop, 'saveName', 'H', 'savePath', strcat(pwd,'\Barrage_Files'));
        load(bigSave);
        paths_save = [paths_save; strcat(main,paths(p))];
        save(bigSave, 'paths_save');
        paths_save = [];
        load(comSave);
        paths_save = [paths_save; strcat(main,paths(p))];
        save(comSave, 'paths_save');
        
    end
    if ifAnalysis
        BarAnalysis(strcat(main,paths(p)));
    end
    if ifPSTH
        regionPSTH();
    end
end

if ifCum
    CumMetRun("mouse");
    CumMetRun("rat");
    CumMetRun("both");
end