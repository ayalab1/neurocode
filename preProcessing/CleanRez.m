function rez = CleanRez(rez,varargin)

% Remove noisy clusters from the rez file. This is an optional preprocessing
% step the user may choose to perform before exporting the Kilosort results to
% phy. Alternatively, one may provide a path (of the Kilosort folder), where
% the function would export the "noise" labels for the noisy clusters in
% phy format. If prior labels exist, "noise" labels will be appended to the %
% existing cluster_group.tsv file.
%
%  USAGE
%
%    [cleanRez] = CleanRez(rez, <options>)
%
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'savepath'        folder where phy files have already been exported
%                       (by default, undefined as you may run CleanRez before
%                       exporting the clean rez to phy). If provided, cluster
%                       labels will be saved in a cluster_groups.tsv file for phy.
%     'minNumberOfSpikes'   clusters with nSpikes<minNumberOfSpikes will be 
%                       removed (default = 20).
%    =========================================================================
%
%  OUTPUT
%
%    rez                cleaned rez file, without the noisy clusters
%
%
% Copyright (C) 2022 by Ralitsa Todorova and Ryan Harvey
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

minNumberOfSpikes = 20;

% Parse parameter list
for i = 1:2:length(varargin)
    if ~ischar(varargin{i})
        error(['Parameter ' num2str(i+2) ' is not a property (type ''help <a href="matlab:help CleanRez">CleanRez</a>'' for details).']);
    end
    switch(lower(varargin{i}))
        case 'savepath'
            savepath = varargin{i+1};
            if ~isfolder(savepath)
                error('Incorrect value for property ''savepath'' (type ''help <a href="matlab:help CleanRez">CleanRez</a>'' for details).');
            end
            
        otherwise
            error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help CleanRez">CleanRez</a>'' for details).']);
    end
end

nClusters = rez.ops.Nfilt;
nChannels = size(rez.U,1);
% Ignore temporal bins which are always empty (seems to be a Kilosort bug, I haven't investigated)
empty = ~any(any(rez.W,3),2);
nBins = sum(~empty);
% Get the average waveform for each electrode group. We will be comparing these.
meanWaveform = nan(nChannels, nBins, nClusters);
for i = 1:nClusters, meanWaveform(:,:,i) = squeeze(rez.U(:,i,:)) * squeeze(rez.W(~empty,i,:))';end

% Define the waveform range (difference between peak value and trough value).
% The electrode closest to the spike will have the widest range.
mWaveformRange = squeeze(max(meanWaveform,[],2) - min(meanWaveform,[],2));

[signal,noise] = deal(nan(nClusters,1));
for i=1:nClusters
    thismax = mWaveformRange(:,i)==max(mWaveformRange(:,i)); % logical
    % Get a matrix of differences between the ranges of the signal in each group.
    % If the spike is a real spike, the range of the real electrode group will be A LOT higher than the rest.
    % This is what we will use to detect noise clusters, where the ranges of electrode groups are more similar.
    d = abs(bsxfun(@minus,mWaveformRange(:,i),mWaveformRange(:,i)'));
    d(eye(size(d))==1) = nan;
    % Difference between range of the electrode closest to the spike with the ranges of all the other electrodes
    % When the spike is clearly visible in one elecrode group, this should be very large.
    signal(i,1) = nanmean(reshape(d(thismax,~thismax),[],1));
    % Difference between ranges of all the other electrodes (where no real spike should be visible)
    % This is our estimate of differences to be expected due to noise, rather than a real spike
    noise(i,1) = nanmean(reshape(d(~thismax,~thismax),[],1));
end

% signal to noise ratio (snr) below 4 is noisy (imperically defined)
snr = signal./noise;
detectedOnAllElectrodes = snr<3.5;

% Make sure the trough/peak is not a single deviation, but part of at least a couple of points deviating together
% (peak/trough on first or last bin not allowed)
d = nan(nClusters,1);
for i=1:nClusters
    % Define the (closest) channel with the widest waveform range
    thismax = mWaveformRange(:,i)==max(mWaveformRange(:,i)); % logical
    groupID = find(thismax,1); % index
    if isempty(groupID), continue; end
    waveform = (meanWaveform(groupID,:,i))'; waveform = (waveform-median(waveform))./diff(quantile(waveform,[0.25 0.75]));
    
    [maxValue,idx] = max(abs(waveform));
    si = sign(waveform(idx));
    if idx==1 || idx==length(waveform), continue; end
    % make sure the nearest bins are also the same sign (not a single bin fluctuation)
    if abs(waveform(idx-1)-waveform(idx))<abs(waveform(idx+1)-waveform(idx))
        idx1 = idx-1; else; idx1 = idx+1;
    end
    d(i,1) = waveform(idx1)/waveform(idx);
end
singleBin = ~ (d>0);

% Some waveforms look like a noisy oscillation: no single peak/trough:
nPeaksAndTroughs = nan(nClusters,1);
for i=1:nClusters
    thismax = mWaveformRange(:,i)==max(mWaveformRange(:,i)); % logical
    groupID = find(thismax,1); % index
    waveform = nanzscore(meanWaveform(groupID,:,i));
    if any(~isnan(waveform))
        troughs = strfind([nan;diff(Smooth(waveform(:),1))>0]',[0 1])';
        maxTrough = min(waveform(troughs));
        troughs(waveform(troughs)>maxTrough/2) = []; % ignore local minima that are negligible fluctuations compared to the real trough
        nTroughs = length(troughs);
        peaks = strfind([nan;diff(Smooth(waveform(:),1))>0]',[1 0])';
        maxPeak = max(waveform(peaks));
        peaks(waveform(peaks)<maxPeak/2) = []; % ignore local minima that are negligible fluctuations compared to the real trough
        nPeaks = length(peaks);
        nPeaksAndTroughs(i,1) = min([nPeaks nTroughs]);
    end
end
multiplePeaksAndTroughs = nPeaksAndTroughs>1;

isiViolation = false(nClusters,1);
for i=1:nClusters
    spikes = rez.st3(rez.st3(:,2)==i)/rez.ops.fs;
    [acg,acg_t] = CCG(spikes,ones(size(spikes)),'binSize',0.003,'duration',0.031);
    [~,idxMax] = max(acg);
    [~,idx0] = min(abs(acg_t)); % bin t=0
    if idxMax==idx0
        isiViolation(i) = true;
    end
end

% We can also delete clusters with less than a certain number of spikes, e.g. 5
nSpikes = accumarray(rez.st3(:,2),1);
tooFew = nSpikes<minNumberOfSpikes;

detectedOnAllElectrodes(end+1:nClusters) = false;
singleBin(end+1:nClusters) = false;
multiplePeaksAndTroughs(end+1:nClusters) = false;
tooFew(end+1:nClusters) = false;

% Clusters selected by any of the above criteria are marked as "noisy"
noisy = detectedOnAllElectrodes | singleBin | multiplePeaksAndTroughs | isiViolation | tooFew;

% delete noisy clusters to clean the rez file
rez.cProj(ismember(rez.st3(:,2),find(noisy)),:) = [];
rez.st3(ismember(rez.st3(:,2),find(noisy)),:) = [];

% Optionally, export labels for phy
if exist('savepath','var') && exist(savepath,'dir') % if the user has already exported to phy and wants added labels to noise clusters
    display(['Exporting results to ' fullfile(savepath,'cluster_group.tsv')]);
    if ~exist(fullfile(savepath,'cluster_group.tsv'))
        fid = fopen(fullfile(savepath,'cluster_group.tsv'),'w');
        fwrite(fid, sprintf('cluster_id\t%s\r\n', 'group'));
    else
        fid = fopen(fullfile(savepath,'cluster_group.tsv'),'a'); % the file already exists. Add labels to the end
    end
    indices = find(noisy)-1; % phy notation starts from 0
    for i=1:length(indices)
        fwrite(fid, sprintf('%d\t%s\r\n', indices(i), 'noise'));
    end
    fclose(fid);
end