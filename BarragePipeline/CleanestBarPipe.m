%% Cleanest Barrage Pipeline
% Most up-to-date metrics for running analysis through find_HSE_b
% Pick blocks based on what you need

basepath = pwd;
basename = basenameFromBasepath(basepath);
savePath = strcat(basepath, '\Barrage_Files\', basename, '.');

%% New set of spikes if needed
% spikes = importSpikes('brainRegion', 'CA3', 'cellType', 'Pyramidal Cell');
spikes = importSpikes('cellType', 'Pyramidal Cell');
% basepath = pwd;
% basename = basenameFromBasepath(basepath);
save([savePath 'allpyr.cellinfo.mat'], 'spikes');
    
%% Run with presaved metrics
load('C:\Users\Cornell\Documents\GitHub\neurocode\BarragePipeline\curRatMet.mat');
% note_all = ["CA1" "CA2" "CA3" "All pyr"];
% file_n = ['ca1pyr'; 'ca2pyr'; 'ca3pyr'; 'allpyr'];

note_all = ["CA2 and bursty CA3"];
file_n = ['CA2CA3'];

% for i = 1:length(note_all)
%     load([savePath file_n(i,:) '.cellinfo.mat']);
%     barCh(i) = findBarCh(cell_metrics.maxWaveformCh1, spikes.times, spikes.UID);
% end

% note_all = ["CA2"];
% file_n = ['CA2pyr'];
for i = 1:length(note_all)
    load([savePath file_n(i,:) '.cellinfo.mat']);

%     nSigma = ratMet.nSigma;
    nSigma = 3.5;
%     tSepMax = ratMet.tSepMax;
    tSepMax = 1.4;
%     mindur = ratMet.mindur;
    mindur = 0.25;
%     maxdur = ratMet.maxdur;
    maxdur = 5;
%     lastmin = ratMet.lastmin;
    lastmin = 0.3;
    EMGThresh = ratMet.EMGThresh;
    % EMGThresh = 0.8;
    loadspkhist = false;
    save_evts = false;
    neuro2 = false;
    notes = note_all(i);

    find_HSE_b('spikes',spikes,...
                        'nSigma',nSigma,'tSepMax',tSepMax,'mindur',mindur,'maxdur',maxdur,'lastmin',lastmin,...
                        'EMGThresh',EMGThresh,'Notes',notes,'loadspkhist',loadspkhist,...
                        'save_evts',save_evts,'neuro2',neuro2);
    load([savePath 'HSEfutEVT.mat']);
    current = evtSave{end,2};
%     barCh = findBarCh(cell_metrics.maxWaveformCh1, spikes.times, spikes.UID);
%     lfpBar = getLFP(barCh);
    if length(current) >= 500
        [wavAvg,lfpAvg] = eventWavelet(lfpBar,current(1:500),'twin',[1 1], 'frange', [10 250]);
    else
        [wavAvg,lfpAvg] = eventWavelet(lfpBar,current(1:length(current)),'twin',[1 1], 'frange', [10 250]);
    end
    title(notes);
    yticks([10:30:250]);
    saveas(gcf,['Barrage_Profile\barWaveletSample' num2str(size(evtSave,1)) '.png']);
end
load([savePath 'HSEmetrics.mat']);
%NeuroScope2

%% Create Neuroscope2 file for past events
trial = 29;
HSEn2.timestamps = evtSave{trial,1};
HSEn2.peaktimes = evtSave{trial,2};
save([basename '.HSE.events.mat'], 'HSEn2');
NeuroScope2

%% Run with previous metrics (or mostly previous metrics)
note_all = ["CA1" "CA2" "CA3" "All pyr"];
file_n = ['ca1pyr'; 'ca2pyr'; 'ca3pyr'; 'allpyr'];
for i = 1:length(note_all)
    load([basename '.' file_n(i,:) '.cellinfo.mat']);
    load([basename '.HSEmetrics.mat']);
    trial = 1;

%     nSigma = HSEmetrics.nSigma(trial);
    nSigma = 6;
%     tSepMax = HSEmetrics.tSepMax(trial);
    tSepMax = 1.4;
    % mindur = HSEmetrics.mindur(trial);
    mindur = 0.2;
%     maxdur = HSEmetrics.maxdur(trial);
    maxdur = 5;
%     lastmin = HSEmetrics.lastmin(trial);
    lastmin = 0.4;
    EMGThresh = HSEmetrics.EMGThresh(trial);
    % EMGThresh = 0.8;
    loadspkhist = false;
    save_evts = false;
    neuro2 = false;
    notes = note_all(i);

    find_HSE_b('spikes',spikes, ...
                        'nSigma',nSigma,'tSepMax',tSepMax,'mindur',mindur,'maxdur',maxdur,'lastmin',lastmin,...
                        'EMGThresh',EMGThresh,'Notes',notes,'loadspkhist',loadspkhist,...
                        'save_evts',save_evts,'neuro2',neuro2);
    load([basename '.HSEfutEVT.mat']);
    current = evtSave{end,2};
%     barCh = findBarCh(cell_metrics.maxWaveformCh1, spikes.times, spikes.UID);
%     lfpBar = getLFP(172);
    if length(current) >= 500
        [wavAvg,lfpAvg] = eventWavelet(lfpBar,current(1:500),'twin',[1 1], 'frange', [10 250]);
    else
        [wavAvg,lfpAvg] = eventWavelet(lfpBar,current(1:length(current)),'twin',[1 1], 'frange', [10 250]);
    end
    title(notes);
    yticks([10:30:250]);
    saveas(gcf,['Barrage_Profile\barWaveletSample' num2str(size(evtSave,1)) '.png']);
end
load([basename '.HSEmetrics.mat']);
%NeuroScope2

%% Create EVT for past events
trial = 118;
current_time = evtSave{trial,1};
current_peak = evtSave{trial,2};
createEVT(current_time(:,1), current_peak, current_time(:,2), 'saveName', 'T');

%% Create spikes struct with multiple brain regions
load([savePath 'CA2pyr.cellinfo.mat']);
spikes_1 = spikes;
spikes_1_UID = spikes.UID;
spikes_1_times = spikes.times;

clear spikes
load([savePath 'barCA3.cellinfo.mat']);
spikes_2 = spikes;
spikes_2_UID = spikes.UID;
spikes_2_times = spikes.times;

all_UID = [spikes_1_UID spikes_2_UID];
all_times = [spikes_1_times spikes_2_times];
[all_UID idx] = sort([spikes_1_UID spikes_2_UID]);
for i = 1:length(idx)
    sort_times(i) = all_times(idx(i));
end
all_times = sort_times;
clear sort_times;
clear spikes;

spikes.UID = all_UID;
spikes.times = all_times;
save([savePath 'CA2CA3.cellinfo.mat'], 'spikes');

%% Another iteration section
it = 0;
for i = 4:0.5:6
    for j = 0.7:0.1:0.9
        for k = 0.6:0.2:1.4
            nSigma = i;
            tSepMax = k;
            mindur = 0.2;
            maxdur = 3.5;
            lastmin = 0.5;
            EMGThresh = j;
            loadspkhist = true;
            save_evts = false;
            neuro2 = false;

            find_HSE_b('spikes',spikes, ...
                                'nSigma',nSigma,'tSepMax',tSepMax,'mindur',mindur,'maxdur',maxdur,'lastmin',lastmin,...
                                'EMGThresh',EMGThresh,'loadspkhist',loadspkhist,...
                                'save_evts',save_evts,'neuro2',neuro2);
            load('day13.HSEfutEVT.mat');
            
            % lfpBar = getLFP(13);
            for i = 1:10
                runnum = 50;
                current = evtSave{runnum,2};
                evtnum = randi([1 length(current)]);
    %             if length(current) >= 500
    %                 [wavAvg,lfpAvg] = eventWavelet(lfpBar,current(1:500),'twin',[1 1], 'frange', [10 250]);
    %             else
                    [wavAvg,lfpAvg] = eventWavelet(lfpBar,current(randi([1 length(current)])),'twin',[1 1], 'frange', [10 250]);
    %             end
                yticks([10:30:250]);
                saveas(gcf,['Barrage_Profile\SingleEvent\barWaveletSample' num2str(runnum) '_' num2str(evtnum) '.png']);
                close all
            end
            it = it+1;
            disp([num2str(i) ' sigma, ' num2str(j) ' EMGthresh, ' num2str(k) ' tSepMax']);
            disp([num2str(it) '/75 iterations completed']);
        end
    end
end
