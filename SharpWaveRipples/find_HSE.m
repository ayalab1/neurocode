function HSE = find_HSE(spikes,varargin)
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

%% Input handling
p = inputParser;
addParameter(p,'nSigma',3,@isnumeric);
addParameter(p,'tSmooth',.015,@isnumeric); % in s
addParameter(p,'binsz',.001,@isnumeric); % in s
addParameter(p,'UIDs',[],@isnumeric);
addParameter(p,'name',[],@ischar);
addParameter(p,'basename',[],@ischar);
addParameter(p,'save_evts',true,@islogical);
addParameter(p,'mindur',.05,@isnumeric);
addParameter(p,'maxdur',.5,@isnumeric);

parse(p,varargin{:})

nSigma = p.Results.nSigma;
tSmooth = p.Results.tSmooth;
binsz = p.Results.binsz;
UIDs = p.Results.UIDs;
name = p.Results.name;
basename = p.Results.basename;
save_evts = p.Results.save_evts;
mindur = p.Results.mindur;
maxdur = p.Results.maxdur;

%% Set defaults
if isempty(UIDs)
    UIDs = spikes.UID;
end

if isempty(name)
    name = 'HSE';
end

if isempty(basename)
    basename = basenameFromBasepath(pwd);
end

%% Get spike rate over time
allspk = cat(1,spikes.times{UIDs});
allspk = sort(allspk);
ts = 0:binsz:allspk(end);
spkhist = hist(allspk,ts);
spkmean = mean(spkhist);

fsize = tSmooth/binsz;
gfilt = fspecial('gauss',[10*fsize 1],fsize);

spkhist = conv(spkhist,gfilt,'same');
spkhist = zscore(spkhist);

evtidx = spkhist>nSigma;
evtidx = find(diff(evtidx)==1)+1;
belowm = spkhist<spkmean; % Logical to run faster
[startID, stopID, evtstart, evtstop, evtdur, evtamp, evtpeak] = deal(zeros(1,length(evtidx)));
%startID = 1; % Initialize to 1 for the first comparison to work

for e = 1:length(evtidx)
    startID(e) = max([1 find(belowm(1:evtidx(e)),1,'last')]);
    if startID(e)>max(stopID)
        stopID(e) = min([length(belowm) evtidx(e)+find(belowm(evtidx(e):end),1,'first')]);
        evtstart(e) = startID(e)*binsz - binsz;
        evtstop(e) = stopID(e)*binsz - binsz;
        evtdur(e) = (stopID(e) - startID(e))*binsz;
        
        % Get amplitude and peak
        [amp, peakID] = max(spkhist(startID(e):stopID(e)));
        evtamp(e) = amp;
        peakID = peakID + startID(e);
        evtpeak(e) = peakID*binsz - binsz;
    end
end
    

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

goodHSE = find(evtdur>mindur & evtdur<maxdur); % Add mindur here, if desired
evtstart = evtstart(goodHSE);
evtstop = evtstop(goodHSE);
evtdur = evtdur(goodHSE);
evtpeak = evtpeak(goodHSE);

%% TODO: Remove events that are too short evtdur < mindur
% MIN FOR DECODING: 100 ms; split in 20 ms, non overlapping bins.

%% Remove overlapping evts.!

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

if save_evts
    save([basename '.' name '.mat'],'HSE');
end

%% Create FMA .evt structure and save it
% .evt (FMA standard)
if save_evts
    n = length(evtstart);
    d1 = cat(1,evtstart,evtpeak,evtstop);%DS1triad(:,1:3)';
    events1.time = d1(:);
    for i = 1:3:3*n
        events1.description{i,1} = [name ' start'];
        events1.description{i+1,1} = [name ' peak'];
        events1.description{i+2,1} = [name ' stop'];
    end
    
%     SaveEvents([basename '_' name '.HSE.evt'],events1);
    createEVT(HSE.timestamps(:,1), HSE.peaks, HSE.timestamps(:,2), 'saveName', 'H');
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

end