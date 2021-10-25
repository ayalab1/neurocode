function HSE = find_HSE_b(varargin)
% Find high-sychrony events among spikes of a set of units (e.g. for
% decoding).
% For each session the combined spiking of all recorded 
% CA1 pyramidal cells were binned in 1ms bins and convolved with a 15 ms Gaussian 
% kernel (5). For each session a trigger rate was defined as being 3 standard deviations 
% above the mean of all 1 ms bins within NREM epochs of both PRE and POST epochs 
% combined (Grosmark 2016)
%
% INPUTS
%   'spikes'    buzcode compatible 'spikes.cellinfo' struct
%
%   (optional)
%	'algorithm'	Currently supported: 'bayes', 'PVcorr', default 'bayes'
%
% OUTPUT
%   Outvar:   	Description
%
% EXAMPLE CALLS
% [] = DefaultTemplateBuzcode(pwd)
%
% Thomas Hainmueller, 2020, Buzsakilab
% Edited by Lindsay Karaba, 2021, AYA Lab

%%%% This is still VERY messy - I'll work on cleaning this up and splitting 
%%%% into functions throughout the next week or so in order to make it more
%%%% readable.

%% Input handling
p = inputParser;
%loading info
addParameter(p,'spikes',[],@isstruct); 
addParameter(p,'runNum',0,@isnumeric); %number of the current run (so we can use a parfor)
%addParameter(p,'UIDs',[],@isnumeric); %generally not needed
%relevant variables
addParameter(p,'nSigma',3,@isnumeric); %originally 3
addParameter(p,'tSmooth',.015,@isnumeric); % in s
addParameter(p,'binsz',.001,@isnumeric); % in s, originally 0.001
addParameter(p,'tSepMax',1,@isnumeric); %max separation between events
addParameter(p,'mindur',0,@isnumeric); %originally 0.05
addParameter(p,'maxdur',3.5,@isnumeric); %originally 0.5
addParameter(p,'lastmin',0.2,@isnumeric); %last pass for minimum duration
addParameter(p,'EMGThresh',0.9,@isnumeric); %threshold for removing EMG events
addParameter(p,'Notes',[],@isstring); %any relevant run notes
addParameter(p,'sstd',0,@isnumeric); %multiplier for standard deviation to taper detection start time
addParameter(p,'estd',0,@isnumeric); %multiplier for standard deviation to taper detection end time
% logical options
addParameter(p,'remRip',false,@islogical); %are we going to remove ripple times from our spike data
addParameter(p,'loadspkhist',false,@islogical); %option to load zscored firing rate and skip this step
addParameter(p,'saveSpkHist',true,@islogical); %do we want to save our spike histogram array?
addParameter(p,'ifCat',true,@islogical); %should we try to concatenate events?
addParameter(p,'save_evts',false,@islogical);
addParameter(p,'neuro2',false,@islogical); %create a file for use in neuroscope2
addParameter(p,'recordMetrics',true,@islogical); %keep track of metrics for measuring fit
addParameter(p,'futEVT',true,@islogical); %do we want to save timestamp info for future use?
%generally not needed file directions/naming conventions
addParameter(p,'basename',[],@ischar);
addParameter(p,'basepath',pwd,@ischar);
addParameter(p,'name',[],@ischar);

parse(p,varargin{:})

spikes = p.Results.spikes;
runNum = p.Results.runNum;
%UIDs = p.Results.UIDs;
nSigma = p.Results.nSigma;
tSmooth = p.Results.tSmooth;
binsz = p.Results.binsz;
tSepMax = p.Results.tSepMax;
mindur = p.Results.mindur;
maxdur = p.Results.maxdur;
lastmin = p.Results.lastmin;
EMGThresh = p.Results.EMGThresh;
Notes = p.Results.Notes;
sstd = p.Results.sstd;
estd = p.Results.estd;
remRip = p.Results.remRip;
loadspkhist = p.Results.loadspkhist;
saveSpkHist = p.Results.saveSpkHist;
ifCat = p.Results.ifCat;
save_evts = p.Results.save_evts;
neuro2 = p.Results.neuro2;
recordMetrics = p.Results.recordMetrics;
futEVT = p.Results.futEVT;
basename = p.Results.basename;
basepath = p.Results.basepath;
name = p.Results.name;

%% Set defaults
if isempty(name)
    name = 'HSE';
end

if isempty(basename)
    basename = basenameFromBasepath(basepath);
end

if isempty(spikes)&&(~loadspkhist) %if we don't plan to load spikes but they weren't inputted
    load(strcat(basepath, '\', basename, '.spikes.cellinfo.mat'));
end

if ~exist(strcat(basepath,'\','Barrage_Files'))
    mkdir('Barrage_Files');
end
savePath = strcat(basepath, '\Barrage_Files\', basename, '.');

fprintf('\n');
if (recordMetrics&loadspkhist)
    warning('It is not recommended to load prior spike histogram data when recording trial metrics as the readout may not be correctly updated');
end

if ~runNum %if runNum isn't set, find out what it should be
    load([savePath 'HSEfutEVT.mat']);
    runNum = size(evtSave,1) + 1;
    clear evtSave
end

%% Get spike rate over time and save
if ~loadspkhist
    if remRip
        allspk = remSWR(basepath, basename, spikes);
        [spkhist,spkmean,spkstd] = spkRtHist(allspk, tSmooth, binsz);
    else
        allspk = cat(1,spikes.times{:});
        allspk = sort(allspk);
        [spkhist,spkmean,spkstd] = spkRtHist(allspk, tSmooth, binsz);
    end
else
    if remRip
        warning('remRip and loadspkhist has not yet been implemented');
    end
    %disp('Loading in spike rate information...');
    load([savePath 'HSEspkHist.mat']);
    spkhist = HSEspkHist.spkhist;
    spkmean = HSEspkHist.spkmean;
    spkstd = HSEspkHist.spkstd;
    allspk = cat(1,spikes.times{:});
    allspk = sort(allspk);
end

if saveSpkHist
    HSEspkHist.spkhist = spkhist;
    HSEspkHist.spkmean = spkmean;
    HSEspkHist.spkstd = spkstd;
    save([savePath 'HSEspkHist.mat'], 'HSEspkHist');
end

%% Flag first set of events
[evtstart,evtstop,evtdur,evtpeak,evtamp] = fpEvt(spkhist,nSigma,spkmean,spkstd,sstd,estd,binsz,maxdur,mindur);

%% Remove EMG noise
% Pulled from FindRipples.m
if EMGThresh
    %disp('Removing EMG info...');
    sessionInfo = getSession('basepath',basepath); %NEED TO CHANGE 
    EMGfilename = fullfile(basepath,[sessionInfo.general.name '.EMGFromLFP.LFP.mat']);
    if exist(EMGfilename)
        EMGFromLFPvar = load(EMGfilename);   %should use a bz_load script here
        EMGFromLFPvar = EMGFromLFPvar.EMGFromLFP;
    else
        [EMGFromLFPvar] = EMGFromLFP(basepath,'samplingFrequency',10,'savemat',false,'noPrompts',true);
    end
    excluded = logical(zeros(size(evtstart,2),1));
    for i = 1:size(evtstart,2)
       [~, ts] = min(abs(evtstart(i)-EMGFromLFPvar.timestamps)); % get closest sample
       if EMGFromLFPvar.data(ts) > EMGThresh
           excluded(i) = 1;           
       end
    end
    bad = sortrows(evtstart(excluded));
    evtstart = evtstart(~excluded);
    evtstop = evtstop(~excluded);
    evtdur = evtdur(~excluded);
    evtpeak = evtpeak(~excluded);
    evtamp = evtamp(~excluded);
%     %disp(['After EMG noise removal: ' num2str(size(ripples,1)) ' events.']);
    %disp([' >>> Number of events after EMG removal: ' num2str(length(evtstart))]);
end


%% TODO: Remove events that are too short evtdur < mindur
% MIN FOR DECODING: 100 ms; split in 20 ms, non overlapping bins.

%% TODO(LK): Concatenate events that do not satisfy certain delay
% Need to concatenate events that are within a certain threshold to call

%disp('Concatenating close enough events...');
if ifCat
    diffConc = NaN(1, length(evtstop)-1);
    for i = 1:length(evtstart)-1
        diffConc(i) = evtstart(i+1)-evtstop(i); %start of next-stop of last
    end
    flagConc = (diffConc <= tSepMax); %flag indices where event distance is good

    [evtstart, evtstop, evtdur, evtpeak] = CatCon(evtstart,evtstop,evtpeak,evtamp,flagConc);
    %disp([' >>> Number of events after concatenation: ' num2str(length(evtstart))]);
end

%% Remove overlapping evts.!
% Need to put this into a function since it's essentially the above code
%disp('Concatenating overlapping events...');
for i = 1:length(evtstart)-1
    flagConc(i) = (evtstart(i+1) <= evtstop(i)); %flag indices where events overlap
end

[evtstart, evtstop, evtdur, evtpeak] = CatCon(evtstart,evtstop,evtpeak,evtamp,flagConc);
%disp([' >>> Number of events after overlap concatenation: ' num2str(length(evtstart))]);

%% Final pass to remove any short events that did not get concatenated
%disp(['Removing events shorter than ' num2str(lastmin*1000) ' milliseconds...']); %should make this factor an input later
shortPass = find(evtdur>lastmin); %keep events within our bounds
evtstart = evtstart(shortPass);
evtstop = evtstop(shortPass);
evtdur = evtdur(shortPass);
evtpeak = evtpeak(shortPass);
evtamp = evtamp(shortPass);
%disp([' >>> Number of events after min thresholding: ' num2str(length(evtstart))]);

%% Create buzcode event structure and save it
HSE.timestamps = cat(2,evtstart',evtstop');
HSE.peaks = evtpeak';
HSE.amplitudes = evtamp;
HSE.amplitudeUnits = 'spikes';
HSE.eventID = ones(size(evtpeak));
HSE.eventIDlabels = repmat({name},length(evtpeak),1);
HSE.eventIDbinary = false(length(evtpeak),1);
HSE.duration = evtdur;
HSE.center = evtstart + evtdur/2;

HSE.detectorinfo.detectorname = 'find_HSE';
HSE.detectorinfo.detectionparms = [];
HSE.detectorinfo.detectionintervals = [0 Inf];
HSE.detectorinfo.detectiondate = datetime('today');

HSE.detectorinfo.spkhist = spkhist;
HSE.detectorinfo.nSigma = nSigma;
HSE.detectorinfo.tSmooth = tSmooth;
HSE.detectorinfo.binsz = binsz;
HSE.detectorinfo.tSepMax = tSepMax;
HSE.detectorinfo.mindur = mindur;
HSE.detectorinfo.maxdur = maxdur;
HSE.detectorinfo.lastmin = lastmin;
HSE.detectorinfo.ifCat = ifCat;

if save_evts
    %disp('Saving HSE struct...');
    save([savePath name '.mat'],'HSE');
end

%% Create FMA .evt structure and save it
% .evt (FMA standard)
if save_evts
    %disp('Saving .evt...');
    n = length(evtstart);
    d1 = cat(1,evtstart,evtpeak,evtstop);%DS1triad(:,1:3)';
    events1.time = d1(:);
    for i = 1:3:3*n
        events1.description{i,1} = [name ' start'];
        events1.description{i+1,1} = [name ' peak'];
        events1.description{i+2,1} = [name ' stop'];
    end
    
%     SaveEvents([basename '_' name '.HSE.evt'],events1);
    createEVT(HSE.timestamps(:,1), HSE.peaks, HSE.timestamps(:,2), 'saveName', 'H', 'savePath', strcat(pwd,'\Barrage_Files'));
end

if neuro2
    %disp('Saving neuroscope2 file...');
    HSEn2.timestamps = HSE.timestamps;
    HSEn2.peaktimes = HSE.peaks;
    save([basepath '\' basename '.' name '.events.mat'], 'HSEn2'); %save in main path for neuroscope2
end

if futEVT
    %disp('Saving data for creating .evt files in the future...');
    evtFiles = dir([savePath 'HSEfutEVT.mat']);
    if isempty(evtFiles) %check if this exists yet, if not, create and fill
        evtSave = cell(1,2);
        evtSave{1,1} = HSE.timestamps;
        evtSave{1,2} = HSE.peaks;
    else
        load([savePath 'HSEfutEVT.mat']);
        evtSave{runNum,1} = HSE.timestamps;
        evtSave{runNum,2} = HSE.peaks;
    end
    save([savePath 'HSEfutEVT.mat'], 'evtSave');
end

% DS2peakIDCs = round(DS2triad(:,2)*1250);
% HSE.amplitudes = (res_lfp_h(DS2peakIDCs,2)-res_lfp_m(DS2peakIDCs,2))*.00038; % .00038 = conversion to mV
% % amplitude: amplitude of each event (Px1).
% % amplitudeUnits: specify the units of the amplitude vector.
% HSE.amplitudeUnits = 'mV';
% % eventID: numeric ID for classifying various event types (Px1).
% HSE.eventID = ones(size(DS2triad,1),1);
% % eventIDlabels: cell array with labels for classifying various event types defined in stimID (cell array, Px1).
% HSE.eventIDlabels = repmat({'DS2'},size(DS2triad,1),1);
% % eventIDbinary: boolean specifying if eventID should be read as binary values (default: false).
% HSE.eventIDbinary = false(size(DS2triad,1),1);
% % center: center time-point of event (in seconds; calculated from timestamps; Px1).
% % duration: duration of event (in seconds; calculated from timestamps; Px1).
% HSE.duration = DS2triad(:,2)-DS2triad(:,1);
% HSE.center = DS2triad(:,1)+HSE.duration;

% figure; hold on;
% ts = 1:length(tshist);
% plot(ts,tshist,'k');
% plot(ts(evts),tshist(evts),'r');
% xlim([0 10000]);

%% Save our metrics for understanding what works vs not
nEvents = length(HSE.peaks);
avgLen = mean(HSE.duration);
sep = NaN(1,length(evtstart)-1);
for i = 1:length(evtstart)-1
    sep(i) = evtstart(i+1)-evtstop(i);
end
avgSep = mean(sep);
FR = NaN(length(evtstart),1);
barspk = cell(length(evtstart),1);
for i = 1:length(evtstart)
    temp_ct = [];
    temp_ct = find((allspk<=evtstop(i))&(allspk>=evtstart(i)));
    barspk{i} = allspk(temp_ct);
    FR(i) = length(temp_ct)/(evtstop(i)-evtstart(i));
end
avgFR = mean(FR);
temp_ct = [];
for i = 1:length(barspk)
    barT = barspk{i};
    for j = 1:(length(barT)-1)
        temp_ct(end+1) = barT(j+1)-barT(j);
    end
end
avgISI = mean(temp_ct);

if recordMetrics
    %disp('Saving metrics for later comparison...');
    HSEmetFiles = dir([savePath 'HSEmetrics.mat']);
    if isempty(HSEmetFiles) %check if this exists yet, if not, create and fill
        HSEmetrics = table(datetime('now'),Notes,nEvents,avgISI,avgFR,...
            nSigma,tSepMax,mindur,maxdur,lastmin,sstd,estd,...
            EMGThresh,tSmooth,binsz,avgLen,avgSep,ifCat,...
            'VariableNames',{'date','Notes','nEvents','avgISI','avgFR'...
            'nSigma','tSepMax','mindur','maxdur','lastmin','sstd','estd',...
            'EMGThresh','tSmooth','binsz','avgLen','avgSep','ifCat'});
    else
        load([savePath 'HSEmetrics.mat']);
        newMetrics = {datetime('now') Notes nEvents avgISI avgFR...
        nSigma tSepMax mindur maxdur lastmin sstd estd...
        EMGThresh tSmooth binsz avgLen avgSep ifCat};
        HSEmetrics(runNum,:) = newMetrics;
    end
    save([savePath 'HSEmetrics.mat'], 'HSEmetrics');
end

%disp('Done!');
fprintf('\n');


%% Functions
function [allspk] = remSWR(basepath, basename, spikes)
%% Remove sharp wave ripples
% Used for removing sharp wave ripple events from an event set in order to
% improve identification of non-SWR events
%
% Should be improved in the future to take in the SWR structure,
% considering this structure may change names based on how it is
% calculated.
%
% L Karaba 10/21
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%disp('Removing Ripples...');
load(strcat(basepath, '\', basename, '.SWR.events.mat'));
spiket = cat(1,spikes.times{:}); 
spiket = sort(spiket);
SWRt = SWR.timestamps;
excludeInd = [];
for i = 1:size(SWRt, 1);
    tempFind = [];
    tempFind = find((spiket > SWRt(i, 1))&(spiket <= SWRt(i,2)));
    excludeInd = [excludeInd; tempFind];
end

allspk = [];
for i = 1:length(spiket)
    if isempty(find(excludeInd == i))
        allspk = [allspk; spiket(i)];
    end
end
end

function [evtstart,evtstop,evtdur,evtpeak,evtamp] = fpEvt(spkhist,nSigma,spkmean,spkstd,sstd,estd,binsz,maxdur,mindur)
%% First pass for flagging events
%
% This function is for taking a first pass at flagging events based on the
% set parameters.
%
% Note that sstd and estd are not fully implemented yet, but could
% potentially be used to taper the boundary thresholds
%
% LKaraba 10/21
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

evtidx = spkhist>nSigma; %flag events outside nSigma stds
%disp([' >>> Number of events initially flagged: ' num2str(length(evtidx))]);
evtidx = find(diff(evtidx)==1)+1; %find where we switch from above/below our acceptance threshold, keep start index
% belowmstart = spkhist<(spkmean+(sstd*spkstd)); % Logical to run faster, what threshold to start an event
% belowmstop = spkhist<(spkmean-(estd*spkstd)); %what threshold to end an event
belowmstart = spkhist<(spkmean+(sstd*spkstd)); % Logical to run faster, what threshold to start an event
belowmstop = spkhist<(spkmean-(estd*spkstd)); %what threshold to end an event
[startID, stopID, evtstart, evtstop, evtdur, evtamp, evtpeak] = deal(zeros(1,length(evtidx)));
%startID = 1; % Initialize to 1 for the first comparison to work

for e = 1:length(evtidx)
    %choose last index that is below our threshold, or the first timepoint
    startID(e) = max([1 find(belowmstart(1:evtidx(e)),1,'last')]); %set startID(e) = 1 or the last index below the spike mean
    %if our starting index is longer than our last stop index, we need to
    %find an endpoint for this event
    if startID(e)>max(stopID) %need to find the end
        %find earliest point when we go below our threshold, or pick end
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

%if an event is filled with 0s, the event was invalid
%disp([' >>> Number of events after initial pull: ' num2str(length(evtstart))]);
goodHSE = find((evtdur<maxdur)&(evtdur>mindur)); %keep events within our bounds
evtstart = evtstart(goodHSE);
evtstop = evtstop(goodHSE);
evtdur = evtdur(goodHSE);
evtpeak = evtpeak(goodHSE);
evtamp = evtamp(goodHSE);
%disp([' >>> Number of events after first length thresholding: ' num2str(length(evtstart))]);
end

end