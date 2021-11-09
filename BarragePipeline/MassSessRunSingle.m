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
paths = "AB4\day03";
% paths = ["AB4\day03"; "AB4\day07"; "AB4\day08"; "AB4\day09"; "AB4\day11";...
%     "AYA6\day17"; "AYA6\day19"; "AYA6\day20"; "AYA7\day19";...
%     "AYA7\day20"; "AYA7\day22"; "AYA7\day24"; "AYA7\day25";...
%     "AYA7\day27"; "AYA7\day30"; "AYA9\day15"; "AYA9\day16"; "AYA9\day17";...
%     "AYA9\day20"; "AYA10\day25"; "AYA10\day27"; "AYA10\day32"; "AYA10\day34"];
% main = "Y:\SMproject\";
% paths = ["AZ1\day13"];

load('C:\Users\Cornell\Documents\GitHub\neurocode\BarragePipeline\curRatMet.mat');
bigSave = 'Z:\home\Lindsay\Barrage\ratPaths.mat'; %change to mouse, potentially - change below as well
comSave = 'Z:\home\Lindsay\Barrage\combinedPaths.mat';

ifHSE = 1;
ifAnalysis = 0;
ifCum = 0;
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
        singleCellDetection();
        if exist([savePath 'HSEfutEVT.mat'])
            load([savePath 'HSEfutEVT.mat']);
            runLast = size(evtSave,1);
        else
            runLast = 0;
        end
            runNum = runLast+1;
            load([basepath '\' basename '.cell_metrics.cellinfo.mat']);
            load([savePath 'brstDt.cellinfo.mat']);
            note_all = "Burst detected";
            
            nSigma = 5;
            tSmooth = 0.03;
            binSz = 0.005;
            tSepMax = 1;
            mindur = 0.2;
            maxdur = 10;
            lastmin = 0.25;
            sstd = 3;
            estd = 1;
            EMGThresh = ratMet.EMGThresh;
            save_evts = true;
            neuro2 = false;
            notes = note_all;

            find_HSE_b('spikes',spikes,...
                                'nSigma',nSigma,'binSz',binSz,'tSmooth',tSmooth,'tSepMax',tSepMax,'mindur',mindur,...
                                'maxdur',maxdur,'lastmin',lastmin,'EMGThresh',EMGThresh,...
                                'Notes',notes,'sstd',sstd,'estd',estd,...
                                'save_evts',save_evts,'neuro2',neuro2,'runNum',runNum);
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
end

if ifCum
    CumMetRun("mouse");
    CumMetRun("rat");
    CumMetRun("both");
end