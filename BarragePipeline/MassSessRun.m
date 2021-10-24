%% MassSessRun
% Most up-to-date metrics for running analysis through find_HSE_b
% Pick blocks based on what you need

paths = ["AYA7\day19"; "AYA7\day20"; "AYA7\day22"; "AYA7\day24"; "AYA7\day25";...
    "AYA7\day27"; "AYA7\day30"]; %not all but a lot to start
    
%done: "AB4\day03";"AB4\day07"; "AB4\day08"; "AB4\day09"; "AB4\day11";...
%     "AYA6\day17"; "AYA6\day19"; "AYA6\day20"; 
%"AYA6\day21"; 'AYA9\day12', 'AYA9\day15', 'AYA9\day16', 'AYA9\day17',
%'AYA9\day20' "AYA6\day24";...
main = 'Z:\Data\AYAold\';

load('C:\Users\Cornell\Documents\GitHub\neurocode\BarragePipeline\curRatMet.mat');

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
%                 disp(check(i));
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
    
    %% Run our tests
    % note_all = ["CA1" "CA2" "All pyr"]; %change these based on the parameters to run
    % file_n = ['CA1pyr'; 'CA2pyr'; 'allpyr'];
%     note_all = ["CA2"];
%     file_n = ['CA2pyr'];
    % load(strcat(basepath,'.',basename,'cell_metrics.cellinfo.mat'));
    % note_all = ["CA2 and bursty CA3"];
    % file_n = ['CA2CA3'];

    % for i = 1:length(note_all)
    %     load([savePath file_n(i,:) '.cellinfo.mat']);
    %     barCh(i) = findBarCh(cell_metrics.maxWaveformCh1, spikes.times, spikes.UID);
    % end
    if exist([savePath 'HSEfutEVT.mat'])
        load([savePath 'HSEfutEVT.mat']);
        runLast = size(evtSave,1);
    else
        runLast = 0;
    end
    % WaitMessage = parfor_wait(length(note_all),'Waitbar',true);
    for i = 1:length(note_all)
        runNum = runLast+i;
        spikes = load([savePath file_n(i,:) '.cellinfo.mat']);
        spikes = spikes.spikes;
        nSigma = ratMet.nSigma;
    %     nSigma = 1.75; %+-1.5 => [-0.5 1.5]
        tSepMax = ratMet.tSepMax;
    %     tSepMax = 1;
        mindur = ratMet.mindur;
    %     mindur = 0.3;
        maxdur = ratMet.maxdur;
    %     maxdur = 15;
        lastmin = ratMet.lastmin;
    %     lastmin = 0.4;
        sstd = ratMet.sstd;
        estd = ratMet.estd;
        EMGThresh = ratMet.EMGThresh;
    %     EMGThresh = 1.1;
    % DO NOT use loadspkhist if running multiple different spike files
        save_evts = false;
        neuro2 = false;
        notes = note_all{i};

        find_HSE_b('spikes',spikes,...
                            'nSigma',nSigma,'tSepMax',tSepMax,'mindur',mindur,...
                            'maxdur',maxdur,'lastmin',lastmin,'EMGThresh',EMGThresh,...
                            'Notes',notes,'sstd',sstd,'estd',estd,...
                            'save_evts',save_evts,'neuro2',neuro2,'runNum',runNum);
    %     load([savePath 'HSEfutEVT.mat']);
    %     current = evtSave{end,2};
    % %     barCh = findBarCh(cell_metrics.maxWaveformCh1, spikes.times, spikes.UID);
    % %     lfpBar = getLFP(barCh(i));
    %     if length(current) >= 500
    %         [wavAvg,lfpAvg] = eventWavelet(lfpBar,current(1:500),'twin',[1 1], 'frange', [10 250]);
    %     else
    %         [wavAvg,lfpAvg] = eventWavelet(lfpBar,current(1:length(current)),'twin',[1 1], 'frange', [10 250]);
    %     end
    %     title(notes);
    %     yticks([10:30:250]);
    %     if ~exist(strcat(basepath,'\','Barrage_Profile'))
    %         mkdir('Barrage_Profile');
    %     end
    %     saveas(gcf,['Barrage_Profile\barWaveletSample' num2str(size(evtSave,1)) '.png']);
    % WaitMessage.Send;
    end
end