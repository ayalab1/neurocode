function bursts = FindBursts(spikes,varargin)

% FindBursts - Find bursts of action potentials in a spike train
% This function implements the method used by Wierzynski 2009. I.e. smooths
% the mua spike train and identified periods where it surpasses a low-threshold,
% provided that a high threshold is surpassed somewhere within the burst
%
%  USAGE
%
%    bursts = FindBursts(spikes,<options>)
%
%    Find spike bursts by smoothing and thresholding the input single unit or
%    MUA spike train. Two thresholds are used: the upper threshold defines a
%    burst, while the lower threshold identifies the beginning and end of the
%    burst. See e.g. Wierzynski et al. (2009
%
% INPUTS
%
% spikes            spike train (either single unit or MUA)
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties         Values
%    -------------------------------------------------------------------------
%     'thresholds'      thresholds for burst beginning and end, and for burst peak,
%                       in multiples of the standard deviation (default = [0 3])
%     'binSize'         bin size (in s, default = 0.001)
%     'smooth'          smoothing kernel width (in ms, default = 0.0167)
%     'intervals'       [start stop] intervals of interests (e.g. slow wave sleep)
%                       Bursts will be confined to these windows and the function
%                       would ignore all activity outside those windows.
%                       (default = [0 Inf])
%     'durations'       [min max] durations. Bursts outside of this range will
%                       be discarded (default = [0 Inf]).
%    =========================================================================
%
% Copyright (C) 2016-2023 by Ralitsa Todorova, Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.
%
%------------------------------------------------------------------------

% Default values

binSize = 0.001; % in s ; 50 ms in Wierzynski
thresholds = [0 3]; % [2 2] in Wierzynski
smooth = (0.05/3); % 3σ = 50 ms in Wierzynski
intervals = [0 Inf];
durations = [0 Inf];
basepath = [];
save_mat = false;

for i = 1:2:length(varargin),
    if ~ischar(varargin{i}),
        error(['Parameter ' num2str(i+2) ' is not a property (type ''help <a href="matlab:help FindBursts">FindBursts</a>'' for details).']);
    end
    switch(lower(varargin{i})),
        case 'thresholds',
            thresholds = varargin{i+1};
            if ~isdvector(thresholds,'#2','<'),
                error('Incorrect value for property ''thresholds'' (type ''help <a href="matlab:help FindBursts">FindBursts</a>'' for details).');
            end
        case 'binsize'
            binSize = varargin{i+1};
            if ~isdvector(binSize,'#1'),
                error('Incorrect value for property ''binSize'' (type ''help <a href="matlab:help FindBursts">FindBursts</a>'' for details).');
            end
        case 'smooth',
            smooth = varargin{i+1};
            if ~isdvector(smooth,'#1','>0'),
                error('Incorrect value for property ''smooth'' (type ''help <a href="matlab:help FindBursts">FindBursts</a>'' for details).');
            end
        case 'intervals',
            intervals = varargin{i+1};
            if ~isdmatrix(intervals) || size(intervals,2) ~= 2,
                error('Incorrect value for property ''intervals'' (type ''help <a href="matlab:help FindBursts">FindBursts</a>'' for details).');
            end
        case 'durations',
            durations = varargin{i+1};
            if ~isdmatrix(durations) || size(durations,2) ~= 2,
                error('Incorrect value for property ''durations'' (type ''help <a href="matlab:help FindBursts">FindBursts</a>'' for details).');
            end
        case 'basepath'
            basepath = varargin{i+1};
            if ~exist(basepath,'dir')
                error('Incorrect value for property ''basepath'' (type ''help <a href="matlab:help FindBursts">FindBursts</a>'' for details).');
            end
        case 'save_mat'
            save_mat = varargin{i+1};
            if ~islogical(save_mat)
                error('Incorrect value for property ''save_mat'' (type ''help <a href="matlab:help FindBursts">FindBursts</a>'' for details).');
            end
        otherwise,
            error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help FindBursts">FindBursts</a>'' for details).']);
    end
end

smooth = smooth/binSize;
% remove spikes outside of the intervals of interest:
spikes = Restrict(spikes(:,1),intervals,'shift','on');
[binned,t] = binspikes(spikes(:,1),1/binSize); t = t(:);
binned = Smooth(binned,smooth);

z = zscore(binned);

% Find periods depassing the low threshold
bursts(:,[1 3]) = t(ToIntervals(z>thresholds(1)));

% Find maximal z value within putative bursts
[in,id] = InIntervals(t,bursts(:,[1 3]));
[peak,idx] = Accumulate(id(in),z(in),'mode','max');
t = t(in);
bursts(:,2) = t(idx);
bursts(:,4) = peak;

% Apply high threshold
pass = peak>thresholds(2);
bursts = bursts(pass,:);

% Move the start of the burst to the first spike:
bursts(:,1) = spikes(FindClosest(spikes,bursts(:,1),'higher'),1);
% Move the end of the burst to the last spike
bursts(:,3) = spikes(FindClosest(spikes,bursts(:,3),'lower'),1);
% Merge overlapping bursts:
[~,target] = ConsolidateIntervals(bursts(:,[1 3]));
overlapping = find(ismember(target,find(Accumulate(target)>1)));
for i=1:2:length(overlapping)
    i1 = overlapping(i);
    i2 = overlapping(i+1);
    bursts(i1,3) = bursts(i2,3);
    bursts(i1,4) = max([bursts(i1,4) bursts(i2,4)]);
end
bursts(overlapping(2:2:end),:) = [];

% shift timestamps back to original time
bursts0 = bursts; bursts(:) = Unshift(bursts0(:),intervals); bursts(:,4) = bursts0(:,4);

% duration threshold
duration = bursts(:,3)-bursts(:,1);
tooShort = duration<durations(1);
tooLong = duration>durations(2);
bursts(tooShort | tooLong,:) = [];

% save event file
if save_mat && ~isempty(basepath)
    save_mat_file(bursts,basepath,intervals,binSize,durations,thresholds,smooth)
end

end

%% HELPER FUNCTIONS

function save_mat_file(bursts,basepath,intervals,binSize,durations,thresholds,smooth)
% package cellexplorer format

bursts_.timestamps = bursts(:,[1, 3]);
bursts_.peaks = bursts(:,2);
bursts_.amplitude = bursts(:,4);
bursts_.amplitudeUnits = 'zscore';
bursts_.eventID = [];
bursts_.eventIDlabels = [];
bursts_.eventIDbinary = [];
bursts_.duration = bursts_.timestamps(:,2) - bursts_.timestamps(:,1);
bursts_.center = bursts_.timestamps(:,1) + bursts_.duration/2;

detectorinfo.detectorname = 'FindBurst.m';
detectorinfo.detectiondate = datetime("today");
detectorinfo.detectionintervals = intervals;
detectorinfo.detectionparms.binSize = binSize;
detectorinfo.detectionparms.durations = durations;
detectorinfo.detectionparms.thresholds = thresholds;
detectorinfo.detectionparms.smooth = smooth;
bursts_.detectorinfo = detectorinfo;

bursts = bursts_;
basename = basenameFromBasepath(basepath);
save(fullfile(basepath,[basename,'.bursts.events.mat']),'bursts')
end

% Chronux' binspikes
function [dN,t]=binspikes(data,Fs,t)
% bin spikes at a specified frequency sampling i.e. sampling rate 1/sampling
% eg: 1ms accuracy use sampling = 1000
% Usage: [dN,t]=binspikes(data,Fs,t)
% Inputs:
% data   (data as a structure array of spike times; or as a single
%        vector of spike times)
% Fs     (binning frequency)
% t      (the minimum and maximum times to be used to form the bins - [mint maxt]
%            - optional. Default use the spike times themselves to
%              determine the location of the bins.
% Note: the times in data can be in any units. However, it is important
% that all units are chosen consistently. So, if spike times are in secs,
% Fs and t (if present) have to be in Hz and secs respectively. If spike
% times are in number of samples, Fs has to be 1, and t has to be in number
% of samples.
% Outputs:
% dN     (output binned spike counts as a matrix defined on bins starting with the
%         earliest spike across all channels and ending with the latest spike)
% t      (lower limit of each bin)
if nargin < 2; error('Need at least two input arguments'); end;
binSize=1/Fs;
binSizemp='';
if isstruct(data);
    C=length(data);
    fnames=fieldnames(data);
    if nargin <3 || isempty(t);
        mintime=zeros(1,C);
        maxtime=zeros(1,C);
        for ch=1:C
            eval(['binSizemp=data(ch).' fnames{1} ';'])
            mintime(ch)=min(binSizemp);
            maxtime(ch)=max(binSizemp);
        end
        mintime=min(mintime);
        maxtime=max(maxtime);
    else
        %        maxtimech=zeros(1,C);
        %        for ch=1:C
        %          eval(['binSizemp=data(ch).' fnames{1} ';'])
        % %          mintimech(ch)=min(binSizemp);
        %          maxtimech(ch)=max(binSizemp);
        %        end
        mintime=t(1);
        maxtime=t(end);
        %        mintimech=min(mintimech);
        %        maxtimech=max(maxtimech);
        %        if maxtimech > max(t); t=[t maxtimech+binSize]; end;
    end
    t=linspace(mintime,maxtime,1+(maxtime-mintime)/binSize);
    for ch=1:C;
        eval(['binSizemp=data(ch).' fnames{1} ';'])
        x=histc(binSizemp,t);
        dN(:,ch)=x(:);
    end
else
    binSizemp=data;
    if nargin < 3;
        mintime=min(binSizemp);
        maxtime=max(binSizemp);
    else
        mintime=t(1);
        maxtime=t(end);
    end
    t=linspace(mintime,maxtime,1+(maxtime-mintime)/binSize);
    if max(binSizemp)>max(t); t=[t maxtime+binSize]; end;
    x=histc(binSizemp,t);
    dN=x(:);
end
end

function timestamps = Unshift(timestamps, intervals),

% The opposite of the option 'shift' in Restrict
%
% Copyright (C) 2018 Ralitsa Todorova
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

shiftedIntervals(:,2) = CumSum(diff(intervals,[],2));
shiftedIntervals(:,1) = [0; shiftedIntervals(1:end-1,2)];
toAdd = intervals(:,1) - shiftedIntervals(:,1);
[ok,w] = InIntervals(timestamps,shiftedIntervals);
timestamps(~ok) = nan;
timestamps(ok) = timestamps(ok) + toAdd(w(ok));

end
