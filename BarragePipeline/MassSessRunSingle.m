%% MassSessRun
% Most up-to-date metrics for running analysis through find_HSE_b
% Pick blocks based on what you need

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
% paths = ["AYA7\day24"];

paths = ["AB4\day03"; "AB4\day07"; "AB4\day08"; "AB4\day09"; "AB4\day11";...
    "AYA6\day17"; "AYA6\day19"; "AYA6\day20"; "AYA7\day19";...
    "AYA7\day20"; "AYA7\day22"; "AYA7\day24"; "AYA7\day25";...
    "AYA7\day27"; "AYA7\day30"; "AYA9\day15"; "AYA9\day16"; "AYA9\day17";...
    "AYA9\day20"; "AYA10\day25"; "AYA10\day27"; "AYA10\day32"; "AYA10\day34"];
% main = "Y:\SMproject\";
% paths = ["AZ1\day13"];

load('C:\Users\Cornell\Documents\GitHub\neurocode\BarragePipeline\curRatMet.mat');
bigSave = 'Z:\home\Lindsay\Barrage\ratPaths.mat'; %change to mouse, potentially - change below as well
comSave = 'Z:\home\Lindsay\Barrage\combinedPaths.mat';

ifHSE = 1;
ifAnalysis = 1;
ifCum = 0;
ifPSTH = 0;
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
        %% New set of spikes
        unitsForDetection(); %threshold with FR
        loadPath = strcat(basepath,'\Barrage_Files\',basename,'.brstDt.cellinfo.mat');
        bound = 3;
        useSpkThresh = 0.3;
        [~,useSpkThresh,~] = pickUns(bound,useSpkThresh,0.05,loadPath);
        % if it's passed, let's keep pushing the boundary to see if we can make it better
        [~,normSpkThresh,~] = pickUns(bound,useSpkThresh,-0.01,loadPath);
        checking = burstCellDetection(normSpkThresh,loadPath); %output for my own sanity
        
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
        for n = 1:length(UIDkeep)
            load([savePath 'brstDt.cellinfo.mat']);
            keep.times = spikes.times{n};
            keep.UID = spikes.UID(n);
            spikes = []; spikes = keep;
            
            nSigma = 5;
            tSmooth = 0.02;
            binSz = 0.005;
            tSepMax = 0.005; %was playing with 0.005, 0.1
            mindur = 0.3;
            maxdur = 10;
            lastmin = 0.3;
            sstd = -1*(nSigma-1);
            estd = (nSigma-1);
            EMGThresh = 0.8;
            save_evts = false;
            neuro2 = false;
            notes = note_all;

            HSE = find_HSE_b('spikes',spikes,...
                        'nSigma',nSigma,'binSz',binSz,'tSmooth',tSmooth,'tSepMax',tSepMax,...
                        'mindur',mindur,'maxdur',maxdur,'lastmin',lastmin,...
                        'EMGThresh',EMGThresh,'Notes',notes,'sstd',sstd,'estd',estd,...
                        'save_evts',save_evts,'neuro2',neuro2,'runNum',runNum);
            keepTimes{n} = HSE.timestamps;
            keepPeaks{n} = HSE.peaks;
        end
        j=1;
        goodTimes = []; goodPeaks = [];
        for i = 1:length(keepTimes)
            if ~isempty(keepTimes{i})
                goodTimes{j} = keepTimes{i};
                goodPeaks{j} = keepPeaks{i};
                j = j+1;
            end
        end
        keepTimes = []; keepPeaks = []; keepTimes = goodTimes; keepPeaks = goodPeaks;
        tSepMax = 0.005;
        for n = 1:size(keepTimes,2)-1
            if n==1
                evtstart = cat(1, keepTimes{n}(:,1), keepTimes{n+1}(:,1));
                evtstop = cat(1, keepTimes{n}(:,2), keepTimes{n+1}(:,2));
                evtpeak = cat(1, keepPeaks{n}, keepPeaks{n+1});
                [evtstart, ind] = sort(evtstart);
                tempstop = []; temppeak = [];
                for e = 1:length(ind)
                    tempstop(e) = evtstop(ind(e));
                    temppeak(e) = evtpeak(ind(e));
                end
                evtstop = tempstop; evtpeak = temppeak;
                evtdur = evtstop-evtstart;
                evtamp = zeros(length(evtstart),1);
                diffConc = [];
                for e = 1:length(evtstart)-1
                    diffConc(e) = evtstart(e+1)-evtstop(e); %start of next-stop of last
                end
                flagConc = (diffConc <= tSepMax); %flag indices where event distance is good
                [start,stop,~,peak,numCat] = CatCon(evtstart,evtstop,evtpeak,evtamp,flagConc);
                start = start'; stop = stop'; peak = peak';
            else
                evtstart = cat(1, start, keepTimes{n+1}(:,1));
                evtstop = cat(1, stop, keepTimes{n+1}(:,2));
                evtpeak = cat(1, peak, keepPeaks{n+1});
                [evtstart, ind] = sort(evtstart);
                tempstop = []; temppeak = [];
                for e = 1:length(ind)
                    tempstop(e) = evtstop(ind(e));
                    temppeak(e) = evtpeak(ind(e));
                end
                evtstop = tempstop; evtpeak = temppeak;
                evtdur = evtstop-evtstart;
                evtamp = zeros(length(evtstart),1);
                diffConc = [];
                for e = 1:length(evtstart)-1
                    diffConc(e) = evtstart(e+1)-evtstop(e); %start of next-stop of last
                end
                flagConc = (diffConc <= tSepMax); %flag indices where event distance is good
                [start,stop,~,peak,numCat] = CatCon(evtstart,evtstop,evtpeak,evtamp,flagConc);
                start = start'; stop = stop'; peak = peak';
            end
        end
        finThresh = 0.4;
        finDur = stop-start;
        start = start(finDur >= finThresh);
        stop = stop(finDur >= finThresh);
        peak = peak(finDur >= finThresh);
        createEVT(start, peak, stop, 'saveName', 'H', 'savePath', strcat(pwd,'\Barrage_Files'));
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
