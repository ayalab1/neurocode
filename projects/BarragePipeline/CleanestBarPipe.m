%% Cleanest Barrage Pipeline
% Most up-to-date metrics for running analysis through find_HSE_b
% Pick blocks based on what you need

basepath = pwd;
basename = basenameFromBasepath(basepath);

if ~exist(strcat(basepath,'\','Barrage_Files'))
    mkdir('Barrage_Files');
end

savePath = strcat(basepath, '\Barrage_Files\', basename, '.');

%% New set of spikes if needed
% spikes = importSpikes('brainRegion', 'CA2', 'cellType', 'Pyramidal Cell');
spikes = [];
spikes = importSpikes('cellType', 'Pyramidal Cell');
save([savePath 'allpyr.cellinfo.mat'], 'spikes');
    
%% Run with presaved metrics
disp('starting');
load('C:\Users\Cornell\Documents\GitHub\neurocode\BarragePipeline\curRatMet.mat');
% note_all = ["CA1" "CA2" "All pyr"]; %change these based on the parameters to run
% file_n = ['CA1pyr'; 'CA2pyr'; 'allpyr'];
note_all = ["CA1" "CA2" "All pyr"];
file_n = ['CA1pyr'; 'CA2pyr'; 'allpyr'];
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
%     tSepMax = ratMet.tSepMax;
    tSepMax = 0.75;
%     mindur = ratMet.mindur;
    mindur = 0.3;
    maxdur = ratMet.maxdur;
%     maxdur = 15;
    lastmin = ratMet.lastmin;
%     lastmin = 0.01;
%     sstd = ratMet.sstd;
    sstd = 3.5;
    estd = ratMet.estd;
    EMGThresh = ratMet.EMGThresh;
%     EMGThresh = 0;
% DO NOT use loadspkhist if running multiple different spike files
    save_evts = false;
    neuro2 = false;
    notes = note_all(i);

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
load([savePath 'HSEmetrics.mat']);
% barSumFigs(HSEmetrics);
% WaitMessage.Destroy
% NeuroScope2

%% Create EVT for past events
load([savePath 'HSEfutEVT.mat']);
trial = 76;
current_time = evtSave{trial,1};
current_peak = evtSave{trial,2};
createEVT(current_time(:,1), current_peak, current_time(:,2), 'saveName', 'H', 'savePath', strcat(pwd,'\Barrage_Files'));

%% Create Neuroscope2 file for past events
savePath = strcat(basepath, '\Barrage_Files\', basename, '.');
load([savePath 'HSEfutEVT.mat']);
load([savePath 'HSE.mat']);
trial = size(evtSave,1);
HSEn2.timestamps = evtSave{trial,1};
HSEn2.timestamps = HSEn2.timestamps(HSE.keep,:);
HSEn2.peaktimes = evtSave{trial,2};
HSEn2.peaktimes = HSEn2.peaktimes(HSE.keep,:);
save([basename '.HSE.events.mat'], 'HSEnrem');
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
    mindur = 0.4;
%     maxdur = HSEmetrics.maxdur(trial);
    maxdur = 5;
%     lastmin = HSEmetrics.lastmin(trial);
    lastmin = 0.3;
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

%% Load what we think are the barrage-y CA3 units
% we'll need to check in neuroscope first to identify
neurokeep = [167:172]; %pulled ALL channels that seemed to have barrage activity
neurokeep = [neurokeep 207:213]; %USE BASE 0 - DIRECTLY FROM NEUROSCOPE

check = [];
for i = 1:length(cell_metrics.brainRegion)
    if contains(cell_metrics.brainRegion{i}, 'CA3')
        check(end+1) = i;
    end
end
chan = NaN(length(check),1);
chan = cell_metrics.maxWaveformCh(check);

keep = [];
for i = 1:length(neurokeep)
    keep = [keep find(chan==neurokeep(i))];
end
load(strcat(basepath,'\',basename,'.spikes.cellinfo.mat'));
spikes = importSpikes('spikes', spikes, 'UID', keep, 'brainRegion', 'CA3', 'cellType', 'Pyramidal Cell');
save([savePath 'barCA3.cellinfo.mat'], 'spikes');

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
clear sort_times spikes

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

%% Check an interval for nSigma value

%5452.148
load('Z:\Data\AYAold\AYA7\day20\Barrage_Files\day20.CA2pyr.cellinfo.mat');
timepoints = 2276;
estd = 2.5;
sstd = -2.5;
% timepoints = [1040.819];
nSigma = 3;
for i = 1:length(timepoints)
    timepoint = timepoints(i);
    binsz = 0.001;
    tSmooth = 0.15;


    allspk = cat(1,spikes.times{:});
    allspk = sort(allspk);
    [spkhist,spkmean,spkstd] = spkRtHist(allspk);
    spkTimes = ((1:length(spkhist))*binsz)-binsz;

    difspk = spkTimes-timepoint;
    ind = find(difspk==min(abs(difspk)));
    tsur = 20; %in seconds
    figure(1);
    hold on
    plot((spkTimes((ind-(tsur/binsz)):(ind+(tsur/binsz)))), spkhist((ind-(tsur/binsz)):(ind+(tsur/binsz))));
    xlabel('Time (s)');
    ylabel('nSigma');
    ylim([-1 6]);
    xline(ind*binsz, '--');
    yline(nSigma,'k');
    yline((spkmean+(sstd*spkstd)), 'g');
    yline((spkmean-(estd*spkstd)), 'r');
%     yline(spkmean+sstd, 'g');
%     yline(spkmean+estd, 'r');

end
hold off

%% I'm losing my mind
maxdur = 15;
mindur = 0.3;
evtidx = spkhist(30000:32000)>nSigma; %flag events outside nSigma stds

figure(4);
    hold on
    plot(spkTimes(30000:32000), spkhist(30000:32000));
    xlabel('Time (s)');
    ylabel('nSigma');
    ylim([-1 6]);
    xline(ind*binsz, '--');
    yline(1.75,'k');
    yline((spkmean+(sstd*spkstd)), 'g');
    yline((spkmean-(estd*spkstd)), 'r');
    plot(spkTimes(30567),spkhist(30567),'.');
    plot(spkTimes(30682),spkhist(30682),'.');
    plot(spkTimes(30693),spkhist(30693),'.');
    plot(spkTimes(30814),spkhist(30814),'.');
hold off
%disp([' >>> Number of events initially flagged: ' num2str(length(evtidx))]);
evtidx = find(diff(evtidx)==1)+1; %find where we switch from above/below our acceptance threshold, mark start
% belowmstart = spkhist<(spkmean+(sstd*spkstd)); % Logical to run faster, what threshold to start an event
% belowmstop = spkhist<(spkmean-(estd*spkstd)); %what threshold to end an event
belowmstart = spkhist(ind-(tsur/binsz)):(ind+(tsur/binsz))<(spkmean+(sstd*spkstd)); % Logical to run faster, what threshold to start an event
belowmstop = spkhist(ind-(tsur/binsz)):(ind+(tsur/binsz))<(spkmean-(estd*spkstd)); %what threshold to end an event
[startID, stopID, evtstart, evtstop, evtdur, evtamp, evtpeak] = deal(zeros(1,length(evtidx)));
%startID = 1; % Initialize to 1 for the first comparison to work

for e = 1:length(evtidx)
    startID(e) = max([1 find(belowmstart(1:evtidx(e)),1,'last')]); %set startID(e) = 1 or the last index below the spike mean
    if startID(e)>max(stopID) %need to find the end
        stopID(e) = min([length(belowmstop) evtidx(e)+find(belowmstop(evtidx(e):end),1,'first')]); %find the earliest point below mean
        evtstart(e) = startID(e)*binsz - binsz; %get time of start (subtract binsz to set to 0)
        evtstop(e) = stopID(e)*binsz - binsz; %get end time
        evtdur(e) = (stopID(e) - startID(e))*binsz; %length of event
        
        % Get amplitude and peak
        [amp, peakID] = max(spkhist(startID(e):stopID(e))); %max fr and location
        evtamp(e) = amp;
        peakID = peakID + startID(e); %update location in reference to whole array
        evtpeak(e) = peakID*binsz - binsz; %get peak time
    end
end
    
% Not sure what this is below, but keeping in case it's useful later
%for e = length(evtidx):-1:1%length(evtidx)-1000
%    %if e==length(evtidx) || ~InIntervals(evtidx(e),[evtstart(e+1) evtstop(e+1)])
%   %singular = evtidx(e) > startID; % Compare to previous start ID to exclude double detection of evts.
%   startID(e) = max([1 find(belowm(1:evtidx(e)),1,'last')]);
%   stopID(e) = min([length(evtidx) evtidx(e)+find(belowm(evtidx(e):end),1,'first')]);
%   
%   if ~isempty(startID) && ~isempty(stopID) %&& singular
%       if e==length(evtidx) || stopID(e) < max(startID)
%       evtstart(e) = startID*binsz - binsz;
%       evtstop(e) = stopID*binsz - binsz;
%       evtdur(e) = (stopID - startID)*binsz;
%       
%       % Get amplitude and peak
%       [evtamp(e), peakID] = max(spkhist(startID:stopID));
%       peakID = peakID + startID;
%       evtpeak(e) = peakID*binsz - binsz;
%   end
% end

% tempavg = mean(evtdur)
% tempstd = std(evtdur)

%disp([' >>> Number of events after initial pull: ' num2str(length(evtstart))]);
goodHSE = find((evtdur<maxdur)&(evtdur>mindur)); %keep events within our bounds
evtstart = evtstart(goodHSE);
evtstop = evtstop(goodHSE);
evtdur = evtdur(goodHSE);
evtpeak = evtpeak(goodHSE);
evtamp = evtamp(goodHSE);


%%
figure(2);
    hold on
    plot((spkTimes((ind-(tsur/binsz)):(ind+(tsur/binsz)))), spkhist((ind-(tsur/binsz)):(ind+(tsur/binsz))));
    xlabel('Time (s)');
    ylabel('nSigma');
    ylim([-1 6]);
    xline(ind*binsz, '--');
    yline(nSigma,'k');
    yline((spkmean+(sstd*spkstd)), 'g');
    yline((spkmean-(estd*spkstd)), 'r');

%% Summary check figures

close all
CA1ses = find(HSEmetrics.Notes == "CA1");
CA2ses = find(HSEmetrics.Notes == "CA2");
CA3ses = find(HSEmetrics.Notes == "CA3");
PYRses = find(HSEmetrics.Notes == "All pyr");

% CA1ses = CA1ses(start:stop);
% CA2ses = CA2ses(start:stop);
% CA3ses = CA3ses(start:stop);
% PYRses = PYRses(start:stop);

avgISI_tot(1) = mean(HSEmetrics.avgISI(CA1ses))*1000;
avgISI_tot(2) = mean(HSEmetrics.avgISI(CA2ses))*1000;
avgISI_tot(3) = mean(HSEmetrics.avgISI(CA3ses))*1000;
avgISI_tot(4) = mean(HSEmetrics.avgISI(PYRses))*1000;
avgFR_tot(1) = mean(HSEmetrics.avgFR(CA1ses));
avgFR_tot(2) = mean(HSEmetrics.avgFR(CA2ses));
avgFR_tot(3) = mean(HSEmetrics.avgFR(CA3ses));
avgFR_tot(4) = mean(HSEmetrics.avgFR(PYRses));
stdISI_tot(1) = std(HSEmetrics.avgISI(CA1ses))*1000;
stdISI_tot(2) = std(HSEmetrics.avgISI(CA2ses))*1000;
stdISI_tot(3) = std(HSEmetrics.avgISI(CA3ses))*1000;
stdISI_tot(4) = std(HSEmetrics.avgISI(PYRses))*1000;
stdFR_tot(1) = std(HSEmetrics.avgFR(CA1ses));
stdFR_tot(2) = std(HSEmetrics.avgFR(CA2ses));
stdFR_tot(3) = std(HSEmetrics.avgFR(CA3ses));
stdFR_tot(4) = std(HSEmetrics.avgFR(PYRses));

figure(1);
hold on
isi = plot(1:4,avgISI_tot, '.');
set(gca,'xtick',[1:4],'xticklabel',{'CA1';'CA2';'CA3';'Pyr'},'xlim',[0 5])
ylabel('ISI (ms)');
title('Average ISI');
er = errorbar(1:4,avgISI_tot,-1*stdISI_tot,stdISI_tot);
er.Color = [0 0 0];                            
er.LineStyle = 'none'; 
hold off

figure(2);
hold on
plot(1:4,avgFR_tot, '.');
set(gca,'xtick',[1:4],'xticklabel',{'CA1';'CA2';'CA3';'Pyr'},'xlim',[0 5])
ylabel('FR (#/s)');
title('Average FR');
er = errorbar(1:4,avgFR_tot,-1*stdFR_tot,stdFR_tot);
er.Color = [0 0 0];                            
er.LineStyle = 'none'; 
hold off
