function Phy2Neurosuite(clustering_path,savePath,datFile)

% Generate Neurosuite files (.clu, .res., .spk, .fet) for
% each shank given the phy files directory and the .dat file.
%
% Hint: this works much faster if the "savePath" is local, transfer your files later
%
% Copyright (C) 2020-2023 by Ralitsa Todorova
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.


if ~exist('clustering_path','var')
    clustering_path = pwd;
end

if ~exist('savePath','var')
    savePath = fileparts(clustering_path);
end

if ~exist('datFile','var')
    [basepath,~] = fileparts(clustering_path);
    [~,basename] = fileparts(basepath);
    datFile = fullfile(basepath,[basename '.dat']);
end

[basepath,basename] = fileparts(datFile);

xmlFile = fullfile(basepath,[basename '.xml']);
if ~exist(xmlFile,'file')
    error(['Please make sure file ' xmlFile ' exists.']);
end

%% Initialize variables
session = import_xml2session(xmlFile,struct);
nChannels = session.extracellular.nChannels;

tsv_check = fullfile(clustering_path,'cluster_info.tsv');
if ~exist(tsv_check, 'file') 
    error('Cluster_info.tsv does not exist! Open up phy and save to create it. Did you open phy earlier and local "Recluster Local PCAs all?"'); 
end
cluster_info_table = importdata(fullfile(clustering_path,'cluster_info.tsv'));
shankID = cluster_info_table.data(:,end);
clusterIDlist = cellfun(@str2double,cluster_info_table.textdata(2:end,1));
if isempty(clusterIDlist)
    clusterIDlist = cluster_info_table.data(:,1);
end

cluster_info_table = readtable(fullfile(clustering_path,'cluster_info.tsv'), 'FileType','text','Delimiter', '\t');
shankID = table2array(cluster_info_table(:,end));
clusterIDlist = table2array(cluster_info_table(:,1));

noisyIDs = clusterIDlist(cellfun(@(x) contains(x,'noise'), table2cell(cluster_info_table(:,6))));
muaIDs = clusterIDlist(cellfun(@(x) contains(x,'mua'), table2cell(cluster_info_table(:,6))));
allSpiketimes = double(readNPY(fullfile(clustering_path, 'spike_times.npy')));
channels = double(readNPY(fullfile(clustering_path,'channel_map.npy')));
channelShanks = double(readNPY(fullfile(clustering_path, 'channel_shanks.npy')))';
spikeIDs = double(readNPY(fullfile(clustering_path, 'spike_clusters.npy')));

%% load dat file into memory map and shift spike times

d = dir(datFile);
nSamplesInDat = d.bytes/2/nChannels;
mmf = memmapfile(datFile, 'Format', {'int16', [nChannels, nSamplesInDat], 'data'});
nSamples = 32;
nFeatures = 4;
borders = [true(5,1); false(nSamples,1); true(5,1)];

% Correct for shifts because Kilosort sometimes fails to center the spikes.
% Compute the mean waveform for each cluster, find the center, and shift the spike times accordingly.
nSelection = 200;
mean_waveforms = nan(length(clusterIDlist),nSamples);
for i=1:length(clusterIDlist),
    ok = find(spikeIDs==clusterIDlist(i));
    if length(ok)<nSelection, continue; end
    selected = ok(randi(length(ok),[nSelection 1]));
    theseChannels = channels(channelShanks==shankID(i)) + 1;

    temp = zeros(length(selected),nSamples);
    for k=1:length(selected)
        spk = mmf.Data.data(theseChannels,allSpiketimes(selected(k))+(-(round(nSamples/2)+4):(round(nSamples/2)+5)));

        % Filter
        filtered = spk;
        for j=5+1:size(spk,2)-5
            filtered(:,j) = bsxfun(@minus,filtered(:,j), median(spk(:,j-5:j+5),2));
        end
        if size(filtered,2) < length(borders)
            filtered(:,length(borders)) = 0;
        end
        filtered(:,borders) = [];
        temp(k,:) = mean(filtered);
    end
    mean_waveforms(i,:) = median(temp);
end
[~,peak] = max(abs(mean_waveforms),[],2); % peak deviation is the center of the spike waveform
peak(peak==1) = nan; % absent waveforms
shift = peak - nSamples/2;
originalSpiketimes = allSpiketimes;
for i=find(abs(shift) >= 1)' % deviations of at least 1 bin
    ok = (spikeIDs==clusterIDlist(i));
    allSpiketimes(ok) = allSpiketimes(ok) + round(shift(i));
end
disp([datestr(clock) ': Shifted spiketimes of  ' num2str(sum(abs(shift)>=1)) ' clusters.']);

%% Loop all shanks and save neurosuite files

try
    for i = unique(shankID(:))'
        currentIDs = clusterIDlist(shankID==i); % clusters on current octrode
        ok = ismember(spikeIDs,currentIDs); % spikes on current octrode

        % Save .clu file (cluster id)
        clusterID = spikeIDs(ok);
        if any(ismember(currentIDs,noisyIDs))
            noisy = ismember(clusterID,noisyIDs);
        else
            noisy = false(size(clusterID));
        end

        if any(ismember(currentIDs,muaIDs))
            mua = ismember(clusterID,muaIDs);
        else
            mua = false(size(clusterID));
        end

        % Change cluster numbers to start with 2 (in a clu-file, cluster 1 is for noise)
        clusterID(noisy) = 0;
        clusterID(mua) = 1;
        nonnoisy = ~noisy & ~mua;

        [~,~,clusterID(nonnoisy)] = unique(clusterID(nonnoisy));
        clusterID(nonnoisy) = clusterID(nonnoisy) + 1;

        dlmwrite(fullfile(savePath,[basename '.clu.' num2str(i)]),[max(clusterID)+1;clusterID]);

        % Save .res file (spikes times)
        spiketimes = originalSpiketimes(ok);
        dlmwrite(fullfile(savePath,[basename '.res.' num2str(i)]),spiketimes,'precision','%.0f');

        spkFilename = fullfile(savePath,[basename '.spk.' num2str(i)]);
        theseChannels = channels(channelShanks==i) + 1;
        pointer = fopen(spkFilename,'w');

        disp([datestr(clock) ': Starting .spk. ' num2str(i) '...']);

        for k=1:length(spiketimes)

            if rem(k,100000)==0
                disp([datestr(clock), ': Done with ', num2str(k),...
                    ' out of ', num2str(length(spiketimes)),...
                    ' spikes for spk.', num2str(i),...
                    ' (' num2str(100*k/length(spiketimes)), '%)...']);
            end

            % load current spike with padding
            indices = spiketimes(k) + (-round(nSamples/2+4):round(nSamples/2+5));
            if indices(1)<1 || indices(end)>nSamplesInDat,
                okIndices = InIntervals(indices,[1 nSamplesInDat]);
                spk = zeros(length(theseChannels),length(indices));
                spk(:,okIndices) = mmf.Data.data(theseChannels,indices(okIndices));
            else
                spk = mmf.Data.data(theseChannels,indices);
            end

            % Filter
            filtered = spk;
            for j=5+1:size(spk,2)-5
                filtered(:,j) = bsxfun(@minus,filtered(:,j), median(spk(:,j-5:j+5),2));
            end
            if size(filtered,2) < length(borders)
                filtered(:,length(borders)) = 0;
            end
            filtered(:,borders) = [];

            % save to file
            fwrite(pointer,filtered(:),'int16');
        end
        fclose(pointer);

        m = memmapfile(spkFilename, 'Format', 'int16','Writable',false);
        data = reshape(m.Data,length(theseChannels),nSamples,length(spiketimes));
        fet = [];
        for ch = 1:size(data,1)
            c = corr(squeeze(double(data(ch,:,:)))');
            [v,~] = eig(c);
            v = v(:,1:nFeatures);
            proj = (v' *double(squeeze(data(ch,:,:))))';
            fet = [fet, proj];
        end
        fet = round((fet)/max(abs(fet(:))) * (2^15 - 1));
        fet(:,end+1) = spiketimes;

        fid = fopen(fullfile(savePath,[basename '.fet.' num2str(i)]),'w');
        fprintf(fid, '%d\n', size(fet,2));
        fprintf(fid,['%d' repmat('\t%d',1,size(fet,2)-1) '\n'],fet');
        fclose(fid);
    end
catch
    keyboard
end