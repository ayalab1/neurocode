%% MassSessRun
% Most up-to-date metrics for running analysis through find_HSE_b
% Pick blocks based on what you need

% Just run rats or mice together for ease to start
paths = ["AB4\day03";"AB4\day07"; "AB4\day08"; "AB4\day09"; "AB4\day11";...
    "AYA6\day17"; "AYA6\day19"; "AYA6\day20";...
    "AYA7\day19"; "AYA7\day20"; "AYA7\day22"; "AYA7\day24"; "AYA7\day25";...
    "AYA7\day27"; "AYA7\day30"]; %not all but a lot to start
    
%done: "AB4\day03";"AB4\day07"; "AB4\day08"; "AB4\day09"; "AB4\day11";...
%     "AYA6\day17"; "AYA6\day19"; "AYA6\day20"; 
%"AYA6\day21"; 'AYA9\day12', 'AYA9\day15', 'AYA9\day16', 'AYA9\day17',
%'AYA9\day20' "AYA6\day24";...
main = 'Z:\Data\AYAold\';

load('C:\Users\Cornell\Documents\GitHub\neurocode\BarragePipeline\curRatMet.mat');
bigSave = 'Z:\home\Lindsay\Barrage\ratPaths.mat'; %change to mouse, potentially - change below as well
cumSave = 'Z:\home\Lindsay\Barrage\combinedPaths.mat';

ifHSE = 0;
ifAnalysis = 0;
ifCum = 0;
if ifHSE
    for p = 1:length(paths)
        cd(strcat(main,paths(p)));
        curPath = convertStringsToChars(paths(p));
        basepath = strcat(main,curPath);
        basename = basenameFromBasepath(basepath);

        if ~exist(strcat(basepath,'\','Barrage_Files'))
            mkdir('Barrage_Files');
        end

        savePath = strcat(basepath, '\Barrage_Files\', basename, '.');

        %% New set of spikes
        note_all = [];
        nac = 1;
        file_n =[];
        load([basepath '\' basename '.cell_metrics.cellinfo.mat']);
        regions = unique(cell_metrics.brainRegion);
        check = ["CA1" "CA2" "CA3"];
        for i = 1:length(check)
            for j=1:length(regions)
                regCheck = regions{j};
                if contains(regCheck, check(i))
                    br = convertStringsToChars(check(i));
                    spikes = importSpikes('brainRegion', br, 'cellType', 'Pyramidal Cell');
                    save([savePath br 'pyr.cellinfo.mat'], 'spikes');
                    note_all{nac} = check(i);
                    file_n(nac,1:6) = strcat(br,'pyr');
                    nac = nac+1;
                end
            end
        end
        spikes = importSpikes('cellType', 'Pyramidal Cell');
        save([savePath 'allpyr.cellinfo.mat'], 'spikes');
        note_all{nac} = "All pyr";
        file_n(nac,1:6) = 'allpyr';

        if exist([savePath 'HSEfutEVT.mat'])
            load([savePath 'HSEfutEVT.mat']);
            runLast = size(evtSave,1);
        else
            runLast = 0;
        end
        for i = 1:length(note_all)
            runNum = runLast+i;
            spikes = load([savePath file_n(i,:) '.cellinfo.mat']);
            spikes = spikes.spikes;
            nSigma = ratMet.nSigma;
            tSepMax = ratMet.tSepMax;
            mindur = ratMet.mindur;
            maxdur = ratMet.maxdur;
            lastmin = ratMet.lastmin;
            sstd = ratMet.sstd;
            estd = ratMet.estd;
            EMGThresh = ratMet.EMGThresh;
            save_evts = false;
            neuro2 = false;
            notes = note_all{i};

            find_HSE_b('spikes',spikes,...
                                'nSigma',nSigma,'tSepMax',tSepMax,'mindur',mindur,...
                                'maxdur',maxdur,'lastmin',lastmin,'EMGThresh',EMGThresh,...
                                'Notes',notes,'sstd',sstd,'estd',estd,...
                                'save_evts',save_evts,'neuro2',neuro2,'runNum',runNum);
        end
        load(bigSave);
        ratPaths(end+1) = strcat(main,paths(p));
        save(bigSave, 'ratPaths');
        load(cumSave);
        combinedPaths(end+1) = strcat(main,paths(p)); 
        save(cumSave, 'combinedPaths');
        
        if ifAnalysis
            BarAnalysis;
        end
    end
end

if ifCum
    CumMetRun("mouse");
    CumMetRun("rat");
    CumMetRun("both");
end