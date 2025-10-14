function UpdatePhyFromNeurosuite(clustering_path,neurosuite_path)

% After running the AutomatedCuration, use this to export the results
% back to a phy format.
% The old phy files will be moved.
%
% Copyright (C) 2023 by Ralitsa Todorova
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

if ~exist('clustering_path','var'), clustering_path = pwd; end
try rmdir(fullfile(clustering_path, '.phy'), 's'); end

%% Load automatic curator data

channels = double(readNPY('channel_map.npy'));
channelShanks = double(readNPY('channel_shanks.npy'))';
possibleShanks = unique(channelShanks(:))';
nShanks = max(channelShanks);
xml = dir([neurosuite_path filesep '*.xml']);
[~,basename] = fileparts(fullfile(xml.folder,xml.name));

%%
pickleCell = cell(nShanks,1);
for shank = possibleShanks
    try
        fid=py.open(fullfile(neurosuite_path,[basename '.' num2str(shank) '.pkl']),'rb');
        data=py.pickle.load(fid);
        for j=1:length(data)
            pickleCell{shank}{j} = double(int64(data(j).pop)');
        end
    end
end

% Decide whether to merge clusters after applying CCG criterion.
% Mark noise clusters as noise
% Concatenate all cluster IDs
add = 1; clus = []; noiseIndices = [];
for shank = possibleShanks
    clu = dlmread(fullfile(neurosuite_path,[basename '.clu.' num2str(shank)])); clu(1) = [];
    res = dlmread(fullfile(neurosuite_path,[basename '.res.' num2str(shank)]));
    [ccg,t] = CCG(res(clu>0)/20000,clu(clu>0),'duration',0.03*2,'binSize',0.001);
    nSpikes = Accumulate(clu(clu>0));
    for group = 2:length(pickleCell{shank})
        proposedGroup = pickleCell{shank}{group}; proposedGroup(proposedGroup==0) = [];
        proposedGroup = sortby(proposedGroup,-nSpikes(proposedGroup));
        proposedGroup1 = proposedGroup;
        i = proposedGroup(findmax(nSpikes(proposedGroup)));
        isDifferent = false(numel(proposedGroup));
        for ii=1:length(proposedGroup)
            for jj=1:length(proposedGroup)
                i = proposedGroup(ii); j = proposedGroup(jj);
                try
                    isDifferent(ii,jj) = kstest2(repelem(t,ccg(:,i,j)),repelem(t,ccg(:,i,i)),'alpha',0.001)...
                        || kstest2(repelem(t,ccg(:,i,j)),repelem(t,ccg(:,j,j)),'alpha',0.001) ||...
                        kstest2(repelem(t,ccg(:,j,j)),repelem(t,ccg(:,i,i)),'alpha',0.001);
                end
            end
        end

        if ~any(isDifferent(:)), % No cluster is different, merge the whole group
            clu(ismember(clu,proposedGroup)) = i; proposedGroup1(ismember(proposedGroup1,proposedGroup)) = i;
        else
            similarity = nan(size(isDifferent));
            for ii=1:length(proposedGroup)
                for jj=1:length(proposedGroup)
                    i = proposedGroup(ii); j = proposedGroup(jj);
                    similarity(ii,jj) = min([nancorr(ccg(:,j,j),ccg(:,i,i)),...
                        nancorr(ccg(:,i,j),ccg(:,i,i)),nancorr(ccg(:,i,j),ccg(:,j,j))]);
                end
            end
            similarity(isnan(similarity)) = 0;
            similarity = min(cat(3,similarity,similarity'),[],3);
            d = isDifferent(~triu(ones(size(similarity)))) - similarity(~triu(ones(size(similarity))))*0.5;
            idx = cluster(linkage(d','complete'),'criterion','distance','cutoff',eps);
            for ii=1:max(idx)
                if any(isDifferent(idx==ii,idx==ii))
                    error('These weren''t supposed to be grouped together...');
                end
                clu(ismember(clu,proposedGroup(idx==ii))) = proposedGroup(find(idx==ii,1));
                proposedGroup1(ismember(proposedGroup1,proposedGroup(idx==ii))) = proposedGroup(find(idx==ii,1));
            end
        end
        yaa{shank,group} = proposedGroup1;
    end
    clu = clu+add; 
    if ~isempty(pickleCell{shank})
        noiseIndices = [noiseIndices; pickleCell{shank}{1}(:)+add];
    end
    add = max(clu)+1;
    clu(:,2) = shank;
    clu(:,3) = res;
    clus = [clus; clu];
end

%% Rename original phy files
rawPath = fullfile(clustering_path,['clustering_copy_' datestr(clock,'yyyy-mm-dd_HHMMSS')]);
try mkdir(rawPath); end
try movefile(fullfile(clustering_path, 'cluster_info.tsv'),fullfile(rawPath, 'cluster_info.tsv')); end
try movefile(fullfile(clustering_path, 'spike_templates.npy'),fullfile(rawPath, 'spike_templates.npy')); end
try movefile(fullfile(clustering_path, 'spike_clusters.npy'),fullfile(rawPath, 'spike_clusters.npy')); end
try movefile(fullfile(clustering_path, 'cluster_group.tsv'),fullfile(rawPath, 'cluster_group.tsv')); end
try movefile(fullfile(clustering_path, 'templates.npy'),fullfile(rawPath, 'templates.npy')); end
try movefile(fullfile(clustering_path, 'template_feature_ind.npy'),fullfile(rawPath, 'template_feature_ind.npy')); end
try movefile(fullfile(clustering_path, 'pc_features.npy'),fullfile(rawPath, 'pc_features.npy')); end
try movefile(fullfile(clustering_path, 'pc_feature_ind.npy'),fullfile(rawPath, 'pc_feature_ind.npy')); end
try movefile(fullfile(clustering_path, 'similar_templates.npy'),fullfile(rawPath, 'similar_templates.npy')); end

%% Write new phy files

clus = sortrows(clus,[3 2 1]);
shankID = clus(:,2);
nUnits = max(clus(:,1));
writeNPY(uint32(clus(:,1))-1, fullfile(clustering_path, 'spike_templates.npy')); % zero indexing
writeNPY(uint32(clus(:,1))-1, fullfile(clustering_path, 'spike_clusters.npy')); % zero indexing

% Get mean spike waveform instead of a "template"
meanWaveforms = single(zeros(nUnits,44,length(channels)));
for shank = possibleShanks
    clu = clus(shankID==shank);
    m = memmapfile(fullfile(neurosuite_path,[basename '.spk.' num2str(shank)]),'Format', 'int16','Writable',false);
    spk = reshape(m.Data,[],32,size(clu,1));
    for j=unique(clu)',
        mm = mean(spk(:,:,clu==j),3);
        meanWaveforms(j,7:end-6,channelShanks==shank) = mm'/1000;
    end
end
writeNPY(meanWaveforms, fullfile(clustering_path, 'templates.npy')); % [nUnits,44,nChannels] single array
writeNPY(meshgrid(0:length(channels)-1,1:nUnits), fullfile(clustering_path, 'templates_ind.npy'));

% first, compute the similiarity between individual templates. Take the top 16 nearest
% neighbours for each template, making up templateFeatureInds.
% e.g. templateFeatureInds(3,1) = 258 would mean that the the 3rd nearest neighbour of
% "template" 1 is "template" 258
flat = double(reshape(meanWaveforms,nUnits,[]))';
c = nancorr(flat); c0 = c; c(eye(size(c))==1) = nan;
nearestNeighbours = uint32(zeros(16,nUnits));
for j=1:nUnits
    [~,top] = sort(-c(:,j));
    nearestNeighbours(:,j) = top(1:16);
end
writeNPY(nearestNeighbours'-1, fullfile(clustering_path, 'template_feature_ind.npy'));% -1 for zero indexing

c0(isnan(c0)) = 0; c0(c0<0) = 0;
writeNPY(single(c0), fullfile(clustering_path, 'similar_templates.npy'));

fets = single(zeros(size(clus,1),3,16)); tops = uint32(zeros(16,nUnits));
for shank = possibleShanks
    fet = dlmread(fullfile(neurosuite_path,[basename '.fet.' num2str(shank)]));
    fet(1,:) = [];
    these_channels = find(channelShanks==shank);
    fet = fet(:,1:length(these_channels)*3);
    fet = permute(reshape(fet,size(fet,1),length(these_channels),3),[1 3 2]);
    fetComplete = fet; if size(fet,3)>16,fet(:,:,17:end) = []; elseif size(fet,3)<16, fet(:,:,end+1:16) = 0; end
    clu = clus(shankID==shank);
    for j=unique(clu)',
        ok = clu==j;
        m = max(squeeze(abs(mean(fetComplete(ok,:,:),1))));
        [~,top] = sort(-m,2);
        if size(top,2)>16, top(:,17:end) = []; elseif size(top,2)<16, top(:,end+1:16) = 1; end
        tops(:,j) = these_channels(top);
        fet(ok,:,:) = fetComplete(ok,:,top);
    end
    fets(shankID==shank,:,:) = single(fet);
end

writeNPY(fets, fullfile(clustering_path, 'pc_features.npy'));
writeNPY(tops'-1, fullfile(clustering_path, 'pc_feature_ind.npy'));% -1 for zero indexing

% Tell phy which clusters were marked as noisy
fid = fopen(fullfile(clustering_path,'cluster_group.tsv'),'w');

fwrite(fid, sprintf('cluster_id\t%s\r\n', 'group'));
indices = noiseIndices-1; % phy notation starts from 0
for i=1:length(indices)
    fwrite(fid, sprintf('%d\t%s\r\n', indices(i), 'noise'));
end
fclose(fid);
end