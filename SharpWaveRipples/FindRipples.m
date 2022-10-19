function [ripples] = FindRipples(varargin)
%FindRipples - Find hippocampal ripples when no SW channel is present
%
% USAGE
%    [ripples] = FindRipples(lfp(:,2),lfp(:,1),<options>)
%    OR
%    [ripples] = FindRipples(basepath,channel,<options>)
%
%    Ripples are detected using the normalized squared signal (NSS) by
%    thresholding the baseline, merging neighboring events, thresholding
%    the peaks, and discarding events with excessive duration.
%    Thresholds are computed as multiples of the standard deviation of
%    the NSS. Alternatively, one can use explicit values, typically obtained
%    from a previous call.  The estimated EMG can be used as an additional
%    exclusion criteria
%
% INPUTS - note these are also name-value pairs.
%    lfp            unfiltered LFP (one channel) to use
%	 timestamps	    timestamps to match filtered variable
%    <options>      optional list of property-value pairs (see tables below)
%
%    OR
%
%    basepath       path to a single session to run findRipples on
%    channel      	Ripple channel to use for detection (1-indexed)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'thresholds'  thresholds for ripple beginning/end and peak, in multiples
%                   of the stdev (default = [2 5]); must be integer values
%     'durations'   min inter-ripple interval and max ripple duration, in ms
%                   (default = [50 500]). 
%     'minDuration' min ripple duration. Keeping this input nomenclature for backwards
%                   compatibility
%     'restrict'    interval used to compute normalization (default = all)
%     'SR'          sampling rate (in Hz) (default = 1250Hz)
%     'stdev'       reuse previously computed stdev
%     'show'        plot results (default = 'off')
%     'noise'       noisy unfiltered channel used to exclude ripple-
%                   like noise (events also present on this channel are
%                   discarded)
%     'passband'    N x 2 matrix of frequencies to filter for ripple detection 
%                   (default = [80 250])
%     'EMGThresh'   0-1 threshold of EMG to exclude noise
%     'saveMat'     logical (default=false) to save in buzcode format
%     'EVENTFILE'   logical (default=true) to save a .evt file for viewing
%                   in Neuroscope
%     'plotType'    1=original version (several plots); 2=only raw lfp
%    =========================================================================
%
% OUTPUT
%
%    ripples        buzcode format .event. struct with the following fields
%                   .timestamps        Nx2 matrix of start/stop times for
%                                      each ripple
%                   .detectorName      string ID for detector function used
%                   .peaks             Nx1 matrix of peak power timestamps 
%                   .stdev             standard dev used as threshold
%                   .noise             candidate ripples that were
%                                      identified as noise and removed
%                   .peakNormedPower   Nx1 matrix of peak power values
%                   .detectorParams    struct with input parameters given
%                                      to the detector
% SEE ALSO
%
%    See also bz_Filter, bz_RippleStats, bz_SaveRippleEvents, bz_PlotRippleStats.

% Copyright (C) 2004-2011 by Michaël Zugaro, initial algorithm by Hajime Hirase
% edited by David Tingley, 2017
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.


% Default values
p = inputParser;
addParameter(p,'basepath',pwd,@isstr)
addParameter(p,'channel',1,@isnumeric)
addParameter(p,'thresholds',[0.5 2.5],@isnumeric)
addParameter(p,'durations',[50 500],@isnumeric)
addParameter(p,'restrict',[],@isnumeric)
addParameter(p,'SR',1250,@isnumeric)
addParameter(p,'stdev',[],@isnumeric)
addParameter(p,'show','off',@isstr)
addParameter(p,'noise',[],@ismatrix)
addParameter(p,'passband',[80 250],@isnumeric)
addParameter(p,'EMGThresh',.9,@isnumeric);
addParameter(p,'saveMat',true,@islogical);
addParameter(p,'EVENTFILE',true,@islogical);
addParameter(p,'minDuration',25,@isnumeric)
addParameter(p,'plotType',2,@isnumeric)

if isstr(varargin{1})  % if first arg is basepath
%     addRequired(p, 'basepath',@isstr);
%     addRequired(p,'channel',@isnumeric) 
    parse(p,varargin{:})
    basename = basenameFromBasepath(p.Results.basepath);
    passband = p.Results.passband;
    EMGThresh = p.Results.EMGThresh;
    lfp = getLFP(p.Results.channel,'basepath',p.Results.basepath);
    signal = FilterLFP(lfp,'passband',passband);
    signal(:,1) = []; % remove timestamps
    timestamps = lfp(:,1);
    basepath = p.Results.basepath;
    channel = p.Results.channel;
elseif isnumeric(varargin{1}) % if first arg is filtered LFP
    addRequired(p,'lfp',@isnumeric)
    addRequired(p,'timestamps',@isnumeric)
    parse(p,varargin{:})
    passband = p.Results.passband;
    EMGThresh = p.Results.EMGThresh;
    basepath = p.Results.basepath;
    timestamps = p.Results.timestamps;
    
    % package into struct for bz_Filter, so we can return struct
    samples.samplingRate = p.Results.SR;
    samples.data = double(p.Results.lfp);
    samples.timestamps = timestamps;
    signal = bz_Filter(samples,'filter','butter','passband',passband,'order',3);
    
    if ~exist('basepath','var')
        basepath = pwd;
    end
    basename = basenameFromBasepath(basepath);
end

% assign parameters (either defaults or given)
SR = p.Results.SR;
show = p.Results.show;
restrict = p.Results.restrict;
sd = p.Results.stdev;
noise = p.Results.noise;
lowThresholdFactor = p.Results.thresholds(1);
highThresholdFactor = p.Results.thresholds(2);
minInterRippleInterval = p.Results.durations(1);
maxRippleDuration = p.Results.durations(2);
minRippleDuration = p.Results.minDuration;
plotType = p.Results.plotType;
saveMat = p.Results.saveMat;
EVENTFILE = p.Results.EVENTFILE;

%% filter and calculate noise


% Parameters
windowLength = SR/SR*11;

% Square and normalize signal
squaredSignal = signal.^2;
% squaredSignal = abs(opsignal);
window = ones(windowLength,1)/windowLength;
keep = [];
if ~isempty(restrict)
    for i=1:size(restrict,1)
        keep = InIntervals(timestamps,restrict);
    end
end
keep = logical(keep); 

[normalizedSquaredSignal,sd] = unity(Filter0(window,sum(squaredSignal,2)),sd,keep);

% Detect ripple periods by thresholding normalized squared signal
thresholded = normalizedSquaredSignal > lowThresholdFactor;
start = find(diff(thresholded)>0);
stop = find(diff(thresholded)<0);
% Exclude last ripple if it is incomplete
if length(stop) == length(start)-1
	start = start(1:end-1);
end
% Exclude first ripple if it is incomplete
if length(stop)-1 == length(start)
    stop = stop(2:end);
end
% Correct special case when both first and last ripples are incomplete
if start(1) > stop(1)
	stop(1) = [];
	start(end) = [];
end
firstPass = [start,stop];
if isempty(firstPass)
	disp('Detection by thresholding failed');
	return
else
	disp(['After detection by thresholding: ' num2str(length(firstPass)) ' events.']);
end

% Merge ripples if inter-ripple period is too short
minInterRippleSamples = minInterRippleInterval/1000*SR;
secondPass = [];
ripple = firstPass(1,:);
for i = 2:size(firstPass,1)
	if firstPass(i,1) - ripple(2) < minInterRippleSamples
		% Merge
		ripple = [ripple(1) firstPass(i,2)];
	else
		secondPass = [secondPass ; ripple];
		ripple = firstPass(i,:);
	end
end
secondPass = [secondPass ; ripple];
if isempty(secondPass)
	disp('Ripple merge failed');
	return
else
	disp(['After ripple merge: ' num2str(length(secondPass)) ' events.']);
end

% Discard ripples with a peak power < highThresholdFactor
thirdPass = [];
peakNormalizedPower = [];
for i = 1:size(secondPass,1)
	[maxValue,maxIndex] = max(normalizedSquaredSignal([secondPass(i,1):secondPass(i,2)]));
	if maxValue > highThresholdFactor
		thirdPass = [thirdPass ; secondPass(i,:)];
		peakNormalizedPower = [peakNormalizedPower ; maxValue];
	end
end
if isempty(thirdPass),
	disp('Peak thresholding failed.');
	return
else
	disp(['After peak thresholding: ' num2str(length(thirdPass)) ' events.']);
end

% Detect negative peak position for each ripple
peakPosition = zeros(size(thirdPass,1),1);
for i=1:size(thirdPass,1)
	[minValue,minIndex] = min(signal(thirdPass(i,1):thirdPass(i,2)));
	peakPosition(i) = minIndex + thirdPass(i,1) - 1;
end

% Discard ripples that are way too long
ripples = [timestamps(thirdPass(:,1)) timestamps(peakPosition) ...
           timestamps(thirdPass(:,2)) peakNormalizedPower];
duration = ripples(:,3)-ripples(:,1);
ripples(duration>maxRippleDuration/1000,:) = NaN;
%disp(['After duration test: ' num2str(size(ripples,1)) ' events.']);

% Discard ripples that are too short
ripples(duration<minRippleDuration/1000,:) = NaN;
ripples = ripples((all((~isnan(ripples)),2)),:);

disp(['After duration test: ' num2str(size(ripples,1)) ' events.']);

% If a noise channel was provided, find ripple-like events and exclude them
bad = [];
if ~isempty(noise)
    if length(noise) == 1 % you gave a channel number
       noiselfp = getLFP(p.Results.noise,'basepath',p.Results.basepath);%currently cannot take path inputs
       squaredNoise = FilterLFP(noiselfp,'passband',passband).^2;
    else
            
	% Filter, square, and pseudo-normalize (divide by signal stdev) noise
	squaredNoise = bz_Filter(double(noise),'filter','butter','passband',passband,'order', 3).^2;
    end
    
	window = ones(windowLength,1)/windowLength;
	normalizedSquaredNoise = unity(Filter0(window,sum(squaredNoise,2)),sd,[]);
	excluded = logical(zeros(size(ripples,1),1));
	% Exclude ripples when concomittent noise crosses high detection threshold
	previous = 1;
	for i = 1:size(ripples,1)
		j = FindInInterval([timestamps],[ripples(i,1),ripples(i,3)],previous);
		previous = j(2);
		if any(normalizedSquaredNoise(j(1):j(2))>highThresholdFactor)
			excluded(i) = 1;
        end
	end
	bad = ripples(excluded,:);
	ripples = ripples(~excluded,:);
	disp(['After ripple-band noise removal: ' num2str(size(ripples,1)) ' events.']);
end
    %% lets try to also remove EMG artifact?
if EMGThresh
    basepath = p.Results.basepath
    sessionInfo = getSession('basepath',basepath); %NEED TO CHANGE 
    EMGfilename = fullfile(basepath,[sessionInfo.general.name '.EMGFromLFP.LFP.mat']);
    if exist(EMGfilename)
        load(EMGfilename)   %should use a bz_load script here
    else
        [EMGFromLFP] = getEMGFromLFP(basepath,'samplingFrequency',10,'savemat',false,'noPrompts',true);
    end
    excluded = logical(zeros(size(ripples,1),1));
    for i = 1:size(ripples,1)
       [a ts] = min(abs(ripples(i,1)-EMGFromLFP.timestamps)); % get closest sample
       if EMGFromLFP.data(ts) > EMGThresh
           excluded(i) = 1;           
       end
    end
    bad = sortrows([bad; ripples(excluded,:)]);
    ripples = ripples(~excluded,:);
    disp(['After EMG noise removal: ' num2str(size(ripples,1)) ' events.']);
end


% Optionally, plot results
if strcmp(show,'on')
  if plotType == 1
	figure;
	if ~isempty(noise)
		MultiPlotXY([timestamps, signal],[timestamps squaredSignal],...
            [timestamps normalizedSquaredSignal],[timestamps squaredNoise],...
            [timestamps squaredNoise],[timestamps normalizedSquaredNoise]);
		nPlots = 6;
		subplot(nPlots,1,3);
 		ylim([0 highThresholdFactor*1.1]);
		subplot(nPlots,1,6);
  		ylim([0 highThresholdFactor*1.1]);
	else
		MultiPlotXY([timestamps, signal],[timestamps squaredSignal],[timestamps normalizedSquaredSignal]);
%  		MultiPlotXY(time,signal,time,squaredSignal,time,normalizedSquaredSignal);
		nPlots = 3;
		subplot(nPlots,1,3);
  		ylim([0 highThresholdFactor*1.1]);
	end
	for i = 1:nPlots
		subplot(nPlots,1,i);
		hold on;
  		yLim = ylim;
		for j=1:size(ripples,1)
			plot([ripples(j,1) ripples(j,1)],yLim,'g-');
			plot([ripples(j,2) ripples(j,2)],yLim,'k-');
			plot([ripples(j,3) ripples(j,3)],yLim,'r-');
			if i == 3
				plot([ripples(j,1) ripples(j,3)],[ripples(j,4) ripples(j,4)],'k-');
			end
		end
		for j=1:size(bad,1)
			plot([bad(j,1) bad(j,1)],yLim,'k-');
			plot([bad(j,2) bad(j,2)],yLim,'k-');
			plot([bad(j,3) bad(j,3)],yLim,'k-');
			if i == 3
				plot([bad(j,1) bad(j,3)],[bad(j,4) bad(j,4)],'k-');
			end
		end
		if mod(i,3) == 0
			plot(xlim,[lowThresholdFactor lowThresholdFactor],'k','linestyle','--');
			plot(xlim,[highThresholdFactor highThresholdFactor],'k-');
		end
    end
  elseif plotType == 2
      lfpPlot = FilterLFP(lfp,'passband',[50 300]);
     
		plot(timestamps,lfpPlot);hold on;
  		yLim = ylim;
		for j=1:size(ripples,1)
			plot([ripples(j,1) ripples(j,1)],yLim,'g-');
			plot([ripples(j,2) ripples(j,2)],yLim,'k-');
			plot([ripples(j,3) ripples(j,3)],yLim,'r-');
			if i == 3
				plot([ripples(j,1) ripples(j,3)],[ripples(j,4) ripples(j,4)],'k-');
			end
		end
		for j=1:size(bad,1)
			plot([bad(j,1) bad(j,1)],yLim,'k-');
			plot([bad(j,2) bad(j,2)],yLim,'k-');
			plot([bad(j,3) bad(j,3)],yLim,'k-');
			if i == 3
				plot([bad(j,1) bad(j,3)],[bad(j,4) bad(j,4)],'k-');
			end
		end
		if mod(i,3) == 0
			plot(xlim,[lowThresholdFactor lowThresholdFactor],'k','linestyle','--');
			plot(xlim,[highThresholdFactor highThresholdFactor],'k-');
		end  
  end
end


%% BUZCODE Struct Output
rips = ripples; clear ripples

ripples.timestamps = rips(:,[1 3]);
ripples.peaks = rips(:,2);            %peaktimes? could also do these as timestamps and then ripples.ints for start/stops?
ripples.peakNormedPower = rips(:,4);  %amplitudes?
ripples.stdev = sd;
if ~isempty(bad)
    ripples.noise.times = bad(:,[1 3]);
    ripples.noise.peaks = bad(:,[2]);
    ripples.noise.peakNormedPower = bad(:,[4]);
else
    ripples.noise.times = [];
    ripples.noise.peaks = [];
    ripples.noise.peakNormedPower = [];
end

% Compute instantaneous frequency
hilb = hilbert(signal);
unwrapped = unwrap(angle(hilb));
frequency = bz_Diff(medfilt1(unwrapped,12),timestamps,'smooth',0);
frequency = frequency/(2*pi);

ripples.amplitude = interp1(timestamps,abs(hilb),ripples.peaks,'linear');
ripples.frequency = interp1(timestamps,frequency,ripples.peaks,'linear');
ripples.duration = ripples.timestamps(:,2) - ripples.timestamps(:,1);


%The detectorinto substructure
detectorinfo.detectorname = 'FindRipples';
detectorinfo.detectiondate = datestr(floor(datenum(clock)));
detectorinfo.detectionintervals = restrict;
detectorinfo.detectionparms = p.Results;
detectorinfo.detectionparms = rmfield(detectorinfo.detectionparms,'noise');
if isfield(detectorinfo.detectionparms,'timestamps')  
    detectorinfo.detectionparms = rmfield(detectorinfo.detectionparms,'timestamps');
end
if isfield(p.Results,'channel')
    detectorinfo.detectionchannel = p.Results.channel;
end
if ~isempty(noise)
    detectorinfo.noisechannel = noise;
else
    detectorinfo.noisechannel = nan;
end


%Put it into the ripples structure
ripples.detectorinfo = detectorinfo;

%Save
if p.Results.saveMat
    save(fullfile(basepath, [basename '.ripples.events.mat']),'ripples')
end

%Create .evt file
if EVENTFILE
    
    Filebase = basenameFromBasepath(basepath);
    Filebase = fullfile(basepath,Filebase);
    [pathname, filename, extname] = fileparts(Filebase);
    if isempty(pathname)
        pathname = pwd;
    end
    
    rippleFiles = dir('*.R*.evt');
    if isempty(rippleFiles)
        fileN = 1;
    else
        %set file index to next available value\
        pat = '.R[0-9].';
        fileN = 0;
        for ii = 1:length(rippleFiles)
            token  = regexp(rippleFiles(ii).name,pat);
            val    = str2double(rippleFiles(ii).name(token+2:token+4));
            fileN  = max([fileN val]);
        end
        fileN = fileN + 1;
    end
    fid = fopen(sprintf('%s%s%s.R%02d.evt',pathname,filesep,filename,fileN),'w');
    % convert detections to milliseconds
    SWR_start   = ripples.timestamps(:,1).*(1000);
    SWR_peak    = ripples.peaks.*(1000);
    SWR_end     = ripples.timestamps(:,2).*(1000);
    fprintf(1,'Writing event file ...\n');
    for ii = 1:size(SWR_peak)
        fprintf(fid,'%9.1f\tstart\n',SWR_start(ii));
        fprintf(fid,'%9.1f\tpeak\n',SWR_peak(ii));
        fprintf(fid,'%9.1f\tstop\n',SWR_end(ii));
    end
    fclose(fid);
end

function y = Filter0(b,x)

if size(x,1) == 1
	x = x(:);
end

if mod(length(b),2)~=1
	error('filter order should be odd');
end

shift = (length(b)-1)/2;

[y0 z] = filter(b,1,x);

y = [y0(shift+1:end,:) ; z(1:shift,:)];

function [U,stdA] = unity(A,sd,restrict)

if ~isempty(restrict),
	meanA = mean(A(restrict));
	stdA = std(A(restrict));
else
	meanA = mean(A);
	stdA = std(A);
end
if ~isempty(sd),
	stdA = sd;
end

U = (A - meanA)/stdA;

