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
addParameter(p,'name','HSE',@ischar);
addParameter(p,'basepath',pwd,@isfolder);
addParameter(p,'save_evts',true,@islogical);
addParameter(p,'mindur',.05,@isnumeric);
addParameter(p,'maxdur',.5,@isnumeric);

parse(p,varargin{:})

nSigma = p.Results.nSigma;
tSmooth = p.Results.tSmooth;
binsz = p.Results.binsz;
UIDs = p.Results.UIDs;
name = p.Results.name;
basepath = p.Results.basepath;
save_evts = p.Results.save_evts;
mindur = p.Results.mindur;
maxdur = p.Results.maxdur;

basename = basenameFromBasepath(basepath);

% pull out UIDs
if isempty(UIDs)
    UIDs = spikes.UID;
    % UIDs is used as an index later on so check and correct
    if ~all(diff(UIDs)==1)
        UIDs = 1:length(UIDs);
    end
end

%% Get spike rate over time
% concat all spikes together and sort
allspk = cat(1,spikes.times{UIDs});
allspk = sort(allspk);

% bin over time
ts = 0:binsz:allspk(end);
spkhist = histcounts(allspk, ts);

% smooth binned spikes
fsize = tSmooth/binsz;
gfilt = fspecial('gauss', [10*fsize 1], fsize);
spkhist = conv(spkhist, gfilt, 'same');

% zscore for later thresholds
spkhist = zscore(spkhist);

% locate above thres
evtidx = spkhist > nSigma;
evtidx = find(diff(evtidx) == 1) + 1;
% locate below thres
belowm = spkhist < 0;

[startID, stopID, evtstart, evtstop, evtdur, evtamp, evtpeak] = deal(zeros(1, length(evtidx)));

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

goodHSE = find(evtdur>mindur & evtdur<maxdur); % Add mindur here, if desired
evtstart = evtstart(goodHSE);
evtstop = evtstop(goodHSE);
evtdur = evtdur(goodHSE);
evtpeak = evtpeak(goodHSE);

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
    save(fullfile(basepath,[basename,'.',name,'.events.mat']),'HSE');
end

%% Create FMA .evt structure and save it
% .evt (FMA standard)
if save_evts
    n = length(evtstart);
    d1 = cat(1,evtstart,evtpeak,evtstop);
    events1.time = d1(:);
    for i = 1:3:3*n
        events1.description{i,1} = [name ' start'];
        events1.description{i+1,1} = [name ' peak'];
        events1.description{i+2,1} = [name ' stop'];
    end
    createEVT(HSE.timestamps(:,1), HSE.peaks, HSE.timestamps(:,2), 'saveName', 'H','basepath',basepath);
end

end