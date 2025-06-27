function [spkhist,spkmean,spkstd,ts] = spkRtHist(allspk, varargin)
p = inputParser;
addParameter(p,'tSmooth',0.015,@isnumeric);
addParameter(p,'binsz',0.001,@isnumeric);
addParameter(p,'ifz',true,@islogical);
addParameter(p,'tmax',[],@isnumeric);
parse(p,varargin{:})
tSmooth = p.Results.tSmooth;
binsz = p.Results.binsz;
ifz = p.Results.ifz;
tmax = p.Results.tmax;
%% Get spike rate over time (primarily used for HSE_b)
if isempty(tmax)
    ts = 0:binsz:allspk(end); %bins for sorting spikes
else
    ts = 0:binsz:tmax;
    allspk = Restrict(allspk,[0 tmax]);
end

% hist_plot = histogram(allspk,ts);
% spkhist = hist_plot.BinCounts; close all
spkhist = hist(allspk,ts); %sort spike times into those bins we created
spkmean = mean(spkhist/binsz); %convert to FR but get mean
spkstd = std(spkhist/binsz); %convert to FR but get std

if ~(tSmooth==0)
    fsize = tSmooth/binsz; 
    gfilt = fspecial('gaussian',[10*fsize 1],fsize); %create our filter
    spkhist = conv(spkhist,gfilt,'same'); %filter our spikes
end
if ifz
    spkhist = zscore(spkhist); %normalize spike times
end
end