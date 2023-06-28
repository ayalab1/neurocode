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


if ~exist('clustering_path','var'), clustering_path = pwd; end
if ~exist('savePath','var'), savePath = clustering_path; end
if ~exist('datFile','var'), [basepath,~] = fileparts(clustering_path); [~,basename] = fileparts(basepath); datFile = fullfile(basepath,[basename '.dat']); end
[basepath,basename] = fileparts(datFile);

xmlFile = fullfile(basepath,[basename '.xml']);
if ~exist(xmlFile,'file'),error(['Please make sure file ' xmlFile ' exists.']); end

%% Initialize variables
session = import_xml2session(xmlFile,struct);
nChannels = session.extracellular.nChannels;

cluster_info_table = importdata(fullfile(clustering_path,'cluster_info.tsv')); 
shankID = cluster_info_table.data(:,end); 
clusterIDlist = cellfun(@str2double,cluster_info_table.textdata(2:end,1)); 
noisyIDs = clusterIDlist(cellfun(@(x) contains(x,'noise'),cluster_info_table.textdata(2:end,6)));
muaIDs = clusterIDlist(cellfun(@(x) contains(x,'mua'),cluster_info_table.textdata(2:end,6)));
allSpiketimes = double(readNPY(fullfile(clustering_path, 'spike_times.npy')));
channels = double(readNPY('channel_map.npy'));
channelShanks = double(readNPY(fullfile(clustering_path, 'channel_shanks.npy'))';
spikeIDs = double(readNPY(fullfile(clustering_path, 'spike_clusters.npy')));
nShanks = max(shankID);

%% Loop all shanks
nSamples = 32;
nFeatures = 4;
borders = [true(5,1); false(nSamples,1); true(5,1)];
for i = 1:nShanks
    currentIDs = clusterIDlist(shankID==i); % clusters on current octrode
    ok = ismember(spikeIDs,currentIDs); % spikes on current octrode
    
    % Save .clu file (cluster id)
    clusterID = spikeIDs(ok);
    if any(ismember(currentIDs,noisyIDs)), noisy = ismember(clusterID,noisyIDs); else, noisy = false(size(clusterID)); end
    if any(ismember(currentIDs,muaIDs)), mua = ismember(clusterID,muaIDs); else, mua = false(size(clusterID)); end
    % Change cluster numbers to start with 2 (in a clu-file, cluster 1 is for noise)
    clusterID(noisy) = 0; clusterID(mua) = 1; nonnoisy = ~noisy & ~mua; 
    [~,~,clusterID(nonnoisy)] = unique(clusterID(nonnoisy));
    clusterID(nonnoisy) = clusterID(nonnoisy)+1;
    dlmwrite(fullfile(savePath,[basename '.clu.' num2str(i)]),[max(clusterID)+1;clusterID]);
%     
%     % Save .res file (spikes times)
    spiketimes = allSpiketimes(ok);
    dlmwrite(fullfile(savePath,[basename '.res.' num2str(i)]),spiketimes,'precision','%.0f');

    spkFilename = fullfile(savePath,[basename '.spk.' num2str(i)]);
    theseChannels = channels(channelShanks==i)+1;
    pointer = fopen(spkFilename,'w');
    disp([datestr(clock) ': Starting .spk. ' num2str(i) '...']);
    for k=1:length(spiketimes)
        if rem(k,100000)==0,
            disp([datestr(clock) ': Done with ' num2str(k) ' out of ' num2str(length(spiketimes)) ' spikes for spk.' num2str(i) ' (' num2str(100*k/length(spiketimes)) '%)...']);
        end
        spk = loadBinary(datFile,'start',spiketimes(k)/20000-(16+5)/20000,'duration',(nSamples+10)/20000,'nChannels',nChannels,'channels',theseChannels)';
        % consider loading a bigger chunk and correcting for small shifts (also noting them, so you can correct the .res file)
        % Filter
        filtered = spk;
        for j=5+1:size(spk,2)-5
            filtered(:,j) = bsxfun(@minus,filtered(:,j),median(spk(:,j-5:j+5),2));
        end
        filtered(:,borders) = []; spk(:,borders) = [];
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
        fet = [fet proj];
    end
    fet = round((fet)/max(abs(fet(:))) * (2^15 - 1));
    fet(:,end+1) = spiketimes;

    fid = fopen(fullfile(savePath,[basename '.fet.' num2str(i)]),'w');
    fprintf(fid, '%d\n', size(fet,2));
    fprintf(fid,['%d' repmat('\t%d',1,size(fet,2)-1) '\n'],fet');
    fclose(fid);
end







