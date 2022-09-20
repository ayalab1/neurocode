function [spkhist,spkmean,spkstd,ts] = spkRtHist(allspk, varargin)
p = inputParser;
addParameter(p,'tSmooth',0.015,@isnumeric);
addParameter(p,'binsz',0.001,@isnumeric);
addParameter(p,'ifz',true,@islogical);
parse(p,varargin{:})
tSmooth = p.Results.tSmooth;
binsz = p.Results.binsz;
ifz = p.Results.ifz;
%% Get spike rate over time (primarily used for HSE_b)
ts = 0:binsz:allspk(end); %bins for sorting spikes
spkhist = hist(allspk,ts); %sort spike times into those bins we created
spkmean = mean(spkhist);
spkstd = std(spkhist);

fsize = tSmooth/binsz; 
gfilt = fspecial('gaussian',[10*fsize 1],fsize); %create our filter
spkhist = conv(spkhist,gfilt,'same'); %filter our spikes
if ifz
    spkhist = zscore(spkhist); %normalize spike times
end
end