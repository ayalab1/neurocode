
basepath =  'N:\OJRproject\OJR49\day11'; 
channel = 86; % consider 96
rippleChannels = 1+[18 30]; % the channel that is lower during the sharp-wave is second 

[parentFolder,basename] = fileparts(basepath);
[~,projectName] = fileparts(parentFolder);
sessionID = [projectName '_' basename];
cd(basepath);
load(fullfile(basepath,[basename '.session.mat']),'session');

% Optionally, make an EEG file
% 
% eegFile = fullfile(basepath,[basename '.eeg']); % This is where I will store the normalized LFP file
% if ~exist(eegFile,'file')
%     lfpFile = fullfile(basepath,[basename '.lfp']);
%     copyfile(lfpFile,eegFile);
%     file = memmapfile(eegFile,'Format','int16','Writable',true);
%     data = reshape(file.Data,128,[]);
%     badChannels = 1+[22 23 31 71 96 109 114 113 117 118 104 119 103 120];
%     okChannels = true(128,1); okChannels(badChannels) = false;
%     m = int16(mean(data(okChannels,:)));
%     newData = bsxfun(@minus,data,m);
%     file.data = newData(:);
%     clear file
% end

%%
try
    SleepStateEpisodes = getStruct(basepath,'SleepStateEpisodes');
    sws = SleepStateEpisodes.ints.NREMepisode;        
end

try
%     channel_mapping
    cell_metrics = getStruct(basepath,'cell_metrics');
    nNeurons = length(cell_metrics.spikes.times);
    regions = cell_metrics.brainRegion;
    regionNames = unique(regions);
    regionCell = zeros(length(regions),1);
    for i=1:length(regionNames)
        regionCell(strcmp(regions,regionNames{i})) = i;
    end
    spikesCell = cell_metrics.spikes.times';
    PFCindex = [find(strcmp(regionNames,'ILA')) find(strcmp(regionNames,'PFC'))];
    HPCindex = find(strcmp(regionNames,'CA1'));
    pfc = []; if ~isempty(PFCindex), pfc = sortrows(Group(spikesCell{regionCell==PFCindex})); end
    hpc = []; if ~isempty(HPCindex), hpc = sortrows(Group(spikesCell{regionCell==HPCindex})); end
    neurons = true;
catch
    neurons = false; pfc = zeros(0,2); hpc = zeros(0,2); spikesCell = {};
end

spikes = sortrows(Group(spikesCell{:}));

% Detect delta waves
if exist(fullfile(basepath,[basename '.deltaWaves.events.mat']),'file')
    load(fullfile(basepath,[basename '.deltaWaves.events.mat']),'deltaWaves');
    try deltas = [deltaWaves.timestamps(:,1) deltaWaves.peaks deltaWaves.timestamps(:,2) deltaWaves.peakNormedPower];
    catch
        deltas = repmat(deltaWaves.peaks,1,2);
    end
else
    tic;
    lfp = GetAyaLFP(channel);
    try display(['loaded! @' num2str(toc)]); end
    [clean,~,badIntervals] = CleanLFP(lfp,'thresholds',[6 10],'manual',true);
    deltas0 = FindDeltaPeaks(clean);
    try EMG = getStruct(basepath,'EMG');
        immobility = EMG.timestamps(FindInterval(EMG.data<0.6)); immobility(diff(immobility,[],2)<1,:) = [];
        deltas0 = Restrict(deltas0,immobility);
        immobilityCheck = true;
    catch
        disp('No EMG detected. No immobility restriction'); immobilityCheck = false;
    end

    % Optionally, view the PFC firing around the detected events
    %     [h,ht] = PETH(pfc(:,1),deltas0(:,2));
    %     PlotColorMap(Shrink(sortby(h,-(deltas0(:,5)-deltas0(:,6))),72,1),'x',ht);
    threshold = 3;
    deltas = deltas0(deltas0(:,5)-deltas0(:,6)>threshold,:); % these thresholds should be manually refined for each session
    deltaWaves.timestamps = deltas(:,[1 3]); deltaWaves.peaks = deltas(:,2);  deltaWaves.peakNormedPower = deltas(:,5);
    if immobilityCheck
        deltaWaves.detectorName = ['channel ' num2str(channel) '(+1), CleanLFP, FindDeltaPeaks, EMG.data<0.6, peak-trough>' num2str(threshold)];
    else
        deltaWaves.detectorName = ['channel ' num2str(channel) '(+1), CleanLFP, FindDeltaPeaks, peak-trough>' num2str(threshold)];
    end
    save(fullfile(basepath,[basenameFromBasepath(basepath) '.deltaWaves.events.mat']),'deltaWaves');
    SaveCustomEvents(fullfile(basepath,'deltas.del.evt'),deltas(:,1:3),{'deltas start','delta peak','deltas stop'});
end


% Detect ripples
if exist(fullfile(basepath,[basename '.ripples.events.mat']),'file')
    load(fullfile(basepath,[basename '.ripples.events.mat']),'ripples');
else
    load(fullfile(basepath,[basename '.session.mat']),'session');
%     rippleChannels = swrChannels('basepath',basepath); % or better yet, pick them on neuroscope (don't forget to add 1! rippleChannels.Ripple_Channel = neuroscope_Channel + 1)
    ripples = DetectSWR(rippleChannels,'basepath',basepath,'saveMat',true,'forceDetect',true,'check',true,'useSPW',false);
end

MergePoints = getStruct(basepath,'MergePoints');
sleep = ConsolidateIntervals(MergePoints.timestamps(cellfun(@(x) any(strfind(lower(x),'sleep')),MergePoints.foldernames),:));
pre = InIntervals(deltas,sleep(1,:));
post = InIntervals(deltas,sleep(end,:));
[h,ht] = PETH(ripples.timestamps(:,1),deltas(:,2));
semplot(ht,h(pre,:),'b');
semplot(ht,h(post,:),'r');



