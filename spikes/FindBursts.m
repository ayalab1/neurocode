function bursts = FindBursts(spikes,varargin)

%FindBursts - Find bursts of action potentials in a spike train.
% This function implements the method used by Wierzynski 2009. I.e. smooths
% the mua spike train and identified periods where it surpasses a low-threshold,
% provided that a high threshold is surpassed somewhere within the burst.
%
%  USAGE
%
%    bursts = FindBursts(spikes,<options>)
%
%    Find spike bursts by smoothing and thresholding the input single unit or
%    MUA spike train. Two thresholds are used: the upper threshold defines a
%    burst, while the lower threshold identifies the beginning and end of the
%    burst. See e.g. Wierzynski et al. (2009).
%
%    spikes         spike train (either single unit or MUA)
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'thresholds'  thresholds for burst beginning and end, in multiples of
%                   the stdev (default = [2 5])
%     'binSize'     bin size (in ms, default = 1)
%     'smooth'      smoothing kernel width (in ms, default = 16.667)
%    =========================================================================
%
% Copyright (C) 2016 by Ralitsa Todorova, Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Default values

binSize = 0.001; % in s ; 50 ms in Wierzynski
thresholds = [2 5]; % [2 2] in Wierzynski
smooth = (0.05/3); % 3σ = 50 ms in Wierzynski

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
        case 'window'
            binSize = varargin{i+1};
            if ~isdvector(binSize,'#1'),
                error('Incorrect value for property ''window'' (type ''help <a href="matlab:help FindBursts">FindBursts</a>'' for details).');
            end
        case 'smooth',
            smooth = varargin{i+1};
            if ~isdvector(smooth,'#1','>0'),
                error('Incorrect value for property ''smooth'' (type ''help <a href="matlab:help FindBursts">FindBursts</a>'' for details).');
            end
        otherwise,
            error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help FindBursts">FindBursts</a>'' for details).']);
    end
end

smooth = smooth/binSize;
spikes = spikes(:,1);
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

end


%% HELPER FUNCTIONS
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